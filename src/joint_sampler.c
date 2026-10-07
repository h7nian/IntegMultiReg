#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <math.h>
#include <string.h>

/* One chain of the joint pMOM/MRF posterior. All scratch allocations belong to
 * R; an unwind handler returns RNG ownership on interrupts and numerical errors. */
typedef struct {
    int n, p, forced, outcome, predictor_valid;
    const double *x, *observed, *scale;
    const int *status;
    double *beta, *latent, *column_squares, *residual, *predictor;
    double variance, residual_shape, residual_rate;
    double *saved_beta, *saved_latent;
} imr_group;

typedef struct {
    int groups, features, states, pairs;
    const int *group, *column;
    unsigned char *pattern;
    double nu, log_normalizer;
    double *theta, *energy, *weights, *cross, *precision, *log_bf;
    double *saved_theta;
} imr_platform;

typedef struct {
    int groups, platforms, draws, burnin, thin, sharing, keep_latent, verbose;
    int rng_active;
    double theta_step, swap_rate, interaction_shape, interaction_rate;
    double swaps_proposed, swaps_accepted, theta_proposed, theta_accepted;
    imr_group *group;
    imr_platform *platform;
    SEXP result;
    double *saved_variance, *saved_density;
} imr_joint;

static SEXP element(SEXP x, const char *name)
{
    SEXP names = getAttrib(x, R_NamesSymbol);
    if (!isNewList(x) || TYPEOF(names) != STRSXP || XLENGTH(names) != XLENGTH(x))
        Rf_error("Invalid native sampler input list");
    for (R_xlen_t i = 0; i < XLENGTH(x); ++i)
        if (!strcmp(CHAR(STRING_ELT(names, i)), name)) return VECTOR_ELT(x, i);
    Rf_error("Missing native sampler input '%s'", name);
    return R_NilValue;
}

static int integer_scalar(SEXP x, const char *name, int minimum)
{
    if (TYPEOF(x) != INTSXP || XLENGTH(x) != 1 || INTEGER(x)[0] < minimum ||
        INTEGER(x)[0] == NA_INTEGER) Rf_error("Invalid native integer '%s'", name);
    return INTEGER(x)[0];
}

static double positive_scalar(SEXP x, const char *name)
{
    if (!isReal(x) || XLENGTH(x) != 1 || !R_FINITE(REAL(x)[0]) || REAL(x)[0] <= 0)
        Rf_error("Invalid native positive scalar '%s'", name);
    return REAL(x)[0];
}

static void finite_vector(SEXP x, R_xlen_t n, const char *name)
{
    if (!isReal(x) || XLENGTH(x) != n) Rf_error("Invalid native vector '%s'", name);
    for (R_xlen_t i = 0; i < n; ++i)
        if (!R_FINITE(REAL(x)[i])) Rf_error("Non-finite native vector '%s'", name);
}

static double pmom_normal(double mean, double sd)
{
    if (!R_FINITE(mean) || !R_FINITE(sd) || sd <= 0)
        Rf_error("Non-finite pMOM conditional; check predictor and prior scales");
    double scale = fmax(fabs(mean), sd), a = mean / scale, b = sd / scale;
    double mixture = a * a / (a * a + b * b);
    for (int trial = 0; ; ++trial) {
        if (trial % 1024 == 0) R_CheckUserInterrupt();
        double z;
        if (runif(0, 1) < mixture) z = rnorm(0, 1);
        else {
            double sign = runif(0, 1) < .5 ? -1 : 1;
            z = sign * sqrt(rchisq(3));
        }
        double numerator = a + b * z;
        double acceptance = numerator * numerator / (2 * (a * a + b * b * z * z));
        if (runif(0, 1) < acceptance) {
            double value = mean + sd * z;
            if (!R_FINITE(value)) Rf_error("Coefficient draw exceeds the numerical range; rescale the model");
            if (value != 0) return value;
        }
    }
}

static double lower_normal(double mean, double sd, double lower)
{
    double log_tail = pnorm5((lower - mean) / sd, 0, 1, 0, 1);
    double value = mean + sd * qnorm5(log(runif(0, 1)) + log_tail, 0, 1, 0, 1);
    if (!R_FINITE(value)) Rf_error("Non-finite augmented response; check response and predictor scales");
    /* Cancellation can put a remote-tail draw one ulp below its support. */
    return fmax(value, lower);
}

static double log_bayes_factor(double precision, double cross, double variance, double scale)
{
    double conditional_variance = variance / precision;
    double mean = cross / precision;
    double z = mean / sqrt(conditional_variance);
    double log_prior_variance = log(scale) + log(variance);
    double value = .5 * (log(conditional_variance) - log_prior_variance) +
        .5 * z * z + 2 * log(hypot(mean, sqrt(conditional_variance))) - log_prior_variance;
    if (!R_FINITE(value)) Rf_error("Non-finite conditional Bayes factor; rescale predictors or priors");
    return value;
}

static double normalizer(const imr_platform *p, int first, int second, double change)
{
    double largest = R_NegInf, total = 0;
    for (int k = 0; k < p->states; ++k) {
        double value = p->energy[k];
        if (change != 0) value += change * p->pattern[k * p->groups + first] *
            p->pattern[k * p->groups + second];
        if (value > largest) largest = value;
    }
    for (int k = 0; k < p->states; ++k) {
        double value = p->energy[k];
        if (change != 0) value += change * p->pattern[k * p->groups + first] *
            p->pattern[k * p->groups + second];
        total += exp(value - largest);
    }
    return largest + log(total);
}

static void refresh_platform_energy(imr_platform *p)
{
    for (int k = 0; k < p->states; ++k) {
        double value = 0;
        for (int i = 0; i < p->groups; ++i) if (p->pattern[k * p->groups + i]) {
            value += p->nu;
            for (int j = 0; j < i; ++j) if (p->pattern[k * p->groups + j])
                value += 2 * p->theta[i + p->groups * j];
        }
        if (!R_FINITE(value)) Rf_error("MRF energy exceeds the numerical range");
        p->energy[k] = value;
    }
    p->log_normalizer = normalizer(p, 0, 0, 0);
}

static double residual_crossproduct(const imr_group *g, int column)
{
    double value = 0;
    for (int i = 0; i < g->n; ++i) value += g->x[i + (R_xlen_t)g->n * column] * g->residual[i];
    return value;
}

static void set_coefficient(imr_group *g, int column, double value)
{
    double delta = value - g->beta[column];
    if (delta == 0) return;
    for (int i = 0; i < g->n; ++i)
        g->residual[i] -= g->x[i + (R_xlen_t)g->n * column] * delta;
    g->beta[column] = value;
    /* A retained density evaluation can cache X beta for the next response update. */
    g->predictor_valid = 0;
}

static void update_response_and_forced(imr_group *g)
{
    for (int i = 0; i < g->n; ++i) {
        double mean = 0;
        if (g->predictor_valid) mean = g->predictor[i];
        else for (int j = 0; j < g->p; ++j) mean += g->x[i + (R_xlen_t)g->n * j] * g->beta[j];
        if (g->outcome == 1) {
            double sign = g->observed[i] == 1 ? 1 : -1;
            g->latent[i] = sign * lower_normal(sign * mean, sqrt(g->variance), 0);
        } else if (g->outcome == 2 && !g->status[i]) {
            g->latent[i] = lower_normal(mean, sqrt(g->variance), g->observed[i]);
        }
        g->residual[i] = g->latent[i] - mean;
    }
    for (int j = 0; j < g->forced; ++j) {
        double diagonal = g->column_squares[j];
        double precision = diagonal + 1 / g->scale[j];
        double cross = residual_crossproduct(g, j) + diagonal * g->beta[j];
        set_coefficient(g, j, pmom_normal(cross / precision, sqrt(g->variance / precision)));
    }
}

static void update_feature(imr_joint *state, imr_platform *p, int feature)
{
    for (int s = 0; s < p->groups; ++s) {
        imr_group *g = state->group + p->group[s];
        int j = p->column[s + p->groups * feature];
        double diagonal = g->column_squares[j];
        p->precision[s] = diagonal + 1 / g->scale[j];
        p->cross[s] = residual_crossproduct(g, j) + diagonal * g->beta[j];
        p->log_bf[s] = log_bayes_factor(p->precision[s], p->cross[s], g->variance, g->scale[j]);
    }
    int selected = 0;
    if (state->sharing) {
        double largest = R_NegInf, sum = 0;
        for (int k = 0; k < p->states; ++k) {
            double value = p->energy[k];
            for (int s = 0; s < p->groups; ++s)
                value += p->pattern[k * p->groups + s] * p->log_bf[s];
            if (!R_FINITE(value)) Rf_error("Selection log weight exceeds the numerical range");
            p->weights[k] = value;
            if (value > largest) largest = value;
        }
        for (int k = 0; k < p->states; ++k) sum += (p->weights[k] = exp(p->weights[k] - largest));
        double target = runif(0, 1) * sum;
        for (selected = 0; selected < p->states - 1; ++selected) {
            target -= p->weights[selected];
            if (target <= 0) break;
        }
    } else {
        /* With zero interactions the full-column conditional factorizes. */
        for (int s = 0; s < p->groups; ++s)
            if (runif(0, 1) < plogis(p->nu + p->log_bf[s], 0, 1, 1, 0)) selected |= 1 << s;
    }
    for (int s = 0; s < p->groups; ++s) {
        imr_group *g = state->group + p->group[s];
        int j = p->column[s + p->groups * feature];
        double value = (selected & (1 << s)) ?
            pmom_normal(p->cross[s] / p->precision[s], sqrt(g->variance / p->precision[s])) : 0;
        set_coefficient(g, j, value);
    }
}

static void swap_features(imr_joint *state, imr_platform *p)
{
    if (p->features < 2) return;
    int attempts = (int)ceil(state->swap_rate * p->features);
    for (int attempt = 0; attempt < attempts; ++attempt) {
        int first = (int)R_unif_index(p->features);
        int second = (int)R_unif_index(p->features - 1);
        if (second >= first) ++second;
        double ratio = 0;
        int affected = 0;
        for (int s = 0; s < p->groups; ++s) {
            imr_group *g = state->group + p->group[s];
            int a = p->column[s + p->groups * first], b = p->column[s + p->groups * second];
            if ((g->beta[a] != 0) == (g->beta[b] != 0)) continue;
            if (g->beta[b] != 0) { int tmp = a; a = b; b = tmp; }
            double cross_columns = 0;
            for (int i = 0; i < g->n; ++i)
                cross_columns += g->x[i + (R_xlen_t)g->n * a] * g->x[i + (R_xlen_t)g->n * b];
            double old_cross = residual_crossproduct(g, a) + g->column_squares[a] * g->beta[a];
            double new_cross = residual_crossproduct(g, b) + cross_columns * g->beta[a];
            double old_precision = g->column_squares[a] + 1 / g->scale[a];
            double new_precision = g->column_squares[b] + 1 / g->scale[b];
            ratio += log_bayes_factor(new_precision, new_cross, g->variance, g->scale[b]) -
                log_bayes_factor(old_precision, old_cross, g->variance, g->scale[a]);
            p->cross[s] = new_cross; p->precision[s] = new_precision;
            ++affected;
        }
        if (!affected) continue;
        if (!R_FINITE(ratio)) Rf_error("Feature-swap log ratio exceeds the numerical range");
        ++state->swaps_proposed;
        if (log(runif(0, 1)) >= fmin(0, ratio)) continue;
        for (int s = 0; s < p->groups; ++s) {
            imr_group *g = state->group + p->group[s];
            int a = p->column[s + p->groups * first], b = p->column[s + p->groups * second];
            if ((g->beta[a] != 0) == (g->beta[b] != 0)) continue;
            if (g->beta[b] != 0) { int tmp = a; a = b; b = tmp; }
            double value = pmom_normal(p->cross[s] / p->precision[s], sqrt(g->variance / p->precision[s]));
            set_coefficient(g, a, 0); set_coefficient(g, b, value);
        }
        ++state->swaps_accepted;
    }
}

static void update_variance(imr_group *g)
{
    double squares = 0;
    int active = 0;
    for (int i = 0; i < g->n; ++i) {
        double residual = g->latent[i];
        for (int j = 0; j < g->p; ++j) residual -= g->x[i + (R_xlen_t)g->n * j] * g->beta[j];
        squares += residual * residual;
    }
    for (int j = 0; j < g->p; ++j) if (g->beta[j] != 0) {
        squares += g->beta[j] * g->beta[j] / g->scale[j];
        ++active;
    }
    double shape = g->residual_shape + g->n / 2.0 + 1.5 * active;
    double rate = g->residual_rate + .5 * squares;
    g->variance = 1 / rgamma(shape, 1 / rate);
    if (!R_FINITE(g->variance) || g->variance <= 0)
        Rf_error("Non-finite residual variance; check the data and prior scales");
}

static void update_interactions(imr_joint *state, imr_platform *p)
{
    for (int j = 1; j < p->groups; ++j) for (int i = 0; i < j; ++i) {
        double old = p->theta[i + p->groups * j];
        double proposal = exp(log(old) + rnorm(0, state->theta_step));
        ++state->theta_proposed;
        if (!R_FINITE(proposal) || proposal <= 0) continue;
        double change = 2 * (proposal - old);
        if (!R_FINITE(change)) continue;
        double new_normalizer = normalizer(p, i, j, change);
        if (!R_FINITE(new_normalizer)) continue;
        int shared = 0;
        for (int f = 0; f < p->features; ++f)
            shared += state->group[p->group[i]].beta[p->column[i + p->groups * f]] != 0 &&
                      state->group[p->group[j]].beta[p->column[j + p->groups * f]] != 0;
        double ratio = change * shared - p->features * (new_normalizer - p->log_normalizer) +
            dgamma(proposal, state->interaction_shape, 1 / state->interaction_rate, 1) -
            dgamma(old, state->interaction_shape, 1 / state->interaction_rate, 1) + log(proposal) - log(old);
        if (ISNAN(ratio)) Rf_error("Undefined interaction acceptance ratio");
        if (log(runif(0, 1)) < fmin(0, ratio)) {
            p->theta[i + p->groups * j] = p->theta[j + p->groups * i] = proposal;
            for (int k = 0; k < p->states; ++k)
                p->energy[k] += change * p->pattern[k * p->groups + i] * p->pattern[k * p->groups + j];
            p->log_normalizer = new_normalizer;
            ++state->theta_accepted;
        }
    }
    /* Rebuild from current theta to prevent cumulative roundoff in cached energies. */
    refresh_platform_energy(p);
}

static double joint_log_density(const imr_joint *state)
{
    double value = 0;
    for (int s = 0; s < state->groups; ++s) {
        imr_group *g = state->group + s;
        for (int i = 0; i < g->n; ++i) {
            double mean = 0;
            for (int j = 0; j < g->p; ++j) mean += g->x[i + (R_xlen_t)g->n * j] * g->beta[j];
            g->predictor[i] = mean;
            value += dnorm4(g->latent[i], mean, sqrt(g->variance), 1);
        }
        /* Reuse only while coefficients stay unchanged; preserve the original
         * arithmetic loop above rather than changing its floating-point contraction. */
        g->predictor_valid = 1;
        for (int j = 0; j < g->p; ++j) if (g->beta[j] != 0)
            value += 2 * log(fabs(g->beta[j])) - log(g->scale[j]) - log(g->variance) +
                dnorm4(g->beta[j], 0, sqrt(g->scale[j] * g->variance), 1);
        value += g->residual_shape * log(g->residual_rate) - lgammafn(g->residual_shape) -
            (g->residual_shape + 1) * log(g->variance) - g->residual_rate / g->variance;
    }
    for (int l = 0; l < state->platforms; ++l) {
        const imr_platform *p = state->platform + l;
        for (int f = 0; f < p->features; ++f) {
            int pattern = 0;
            for (int s = 0; s < p->groups; ++s)
                if (state->group[p->group[s]].beta[p->column[s + p->groups * f]] != 0) pattern |= 1 << s;
            value += p->energy[pattern] - p->log_normalizer;
        }
        if (state->sharing) for (int j = 1; j < p->groups; ++j) for (int i = 0; i < j; ++i)
            value += dgamma(p->theta[i + p->groups * j], state->interaction_shape,
                            1 / state->interaction_rate, 1);
    }
    if (!R_FINITE(value)) Rf_error("Non-finite joint log density; check data and prior scales");
    return value;
}

static SEXP run_chain(void *data)
{
    imr_joint *state = data;
    GetRNGstate(); state->rng_active = 1;
    int total = state->burnin + state->draws * state->thin;
    int report_every = total / 10; if (report_every < 1) report_every = 1;
    for (int iteration = 0; iteration < total; ++iteration) {
        R_CheckUserInterrupt();
        for (int s = 0; s < state->groups; ++s) {
            imr_group *g = state->group + s;
            update_response_and_forced(g);
        }
        for (int l = 0; l < state->platforms; ++l) {
            imr_platform *p = state->platform + l;
            if (!p->groups) continue;
            for (int f = 0; f < p->features; ++f) {
                if (f % 64 == 0) R_CheckUserInterrupt();
                update_feature(state, p, f);
            }
        }
        for (int l = 0; l < state->platforms; ++l) swap_features(state, state->platform + l);
        for (int s = 0; s < state->groups; ++s) update_variance(state->group + s);
        if (state->sharing) for (int l = 0; l < state->platforms; ++l)
            update_interactions(state, state->platform + l);
        if (iteration >= state->burnin && (iteration - state->burnin + 1) % state->thin == 0) {
            int row = (iteration - state->burnin + 1) / state->thin - 1;
            for (int s = 0; s < state->groups; ++s) {
                imr_group *g = state->group + s;
                for (int j = 0; j < g->p; ++j) g->saved_beta[row + (R_xlen_t)state->draws * j] = g->beta[j];
                state->saved_variance[row + (R_xlen_t)state->draws * s] = g->variance;
                if (g->saved_latent) for (int i = 0; i < g->n; ++i)
                    g->saved_latent[row + (R_xlen_t)state->draws * i] = g->latent[i];
            }
            for (int l = 0; l < state->platforms; ++l) {
                imr_platform *p = state->platform + l;
                int column = 0;
                for (int j = 1; j < p->groups; ++j) for (int i = 0; i < j; ++i)
                    if (state->sharing) p->saved_theta[row + (R_xlen_t)state->draws * column++] = p->theta[i + p->groups * j];
            }
            state->saved_density[row] = joint_log_density(state);
        }
        if (state->verbose && ((iteration + 1) % report_every == 0 || iteration + 1 == total))
            Rprintf("Iteration %d of %d\n", iteration + 1, total);
    }
    SEXP counts = VECTOR_ELT(state->result, 5);
    REAL(counts)[0] = state->swaps_proposed; REAL(counts)[1] = state->swaps_accepted;
    REAL(counts)[2] = state->theta_proposed; REAL(counts)[3] = state->theta_accepted;
    return state->result;
}

static void release_rng(void *data, Rboolean jump)
{
    imr_joint *state = data;
    (void)jump;
    if (state->rng_active) { PutRNGstate(); state->rng_active = 0; }
}

SEXP imr_joint_sample(SEXP groups, SEXP platforms, SEXP settings, SEXP initial)
{
    if (!isNewList(groups) || !isNewList(platforms) || !XLENGTH(groups) || !XLENGTH(platforms) ||
        XLENGTH(groups) > INT_MAX || XLENGTH(platforms) > INT_MAX)
        Rf_error("Invalid joint sampler groups/platforms");
    imr_joint state = {0};
    state.groups = (int)XLENGTH(groups); state.platforms = (int)XLENGTH(platforms);
    state.draws = integer_scalar(element(settings, "draws"), "draws", 1);
    state.burnin = integer_scalar(element(settings, "burnin"), "burnin", 0);
    state.thin = integer_scalar(element(settings, "thin"), "thin", 1);
    if ((double)state.draws * state.thin + state.burnin > INT_MAX) Rf_error("Too many MCMC iterations");
    state.sharing = integer_scalar(element(settings, "sharing"), "sharing", 0);
    state.keep_latent = integer_scalar(element(settings, "keep_latent"), "keep_latent", 0);
    state.verbose = integer_scalar(element(settings, "verbose"), "verbose", 0);
    if (state.sharing > 1 || state.keep_latent > 1 || state.verbose > 1) Rf_error("Invalid sampler flag");
    state.theta_step = positive_scalar(element(settings, "theta_step"), "theta_step");
    SEXP swap = element(settings, "swap_rate"); finite_vector(swap, 1, "swap_rate");
    state.swap_rate = REAL(swap)[0];
    if (state.swap_rate < 0 || state.swap_rate > 10) Rf_error("Invalid swap proposal rate");
    SEXP interaction_prior = element(settings, "interaction_prior"); finite_vector(interaction_prior, 2, "interaction_prior");
    state.interaction_shape = REAL(interaction_prior)[0]; state.interaction_rate = REAL(interaction_prior)[1];
    if (state.interaction_shape <= 0 || state.interaction_rate <= 0) Rf_error("Invalid interaction prior");
    SEXP beta = element(initial, "coefficients"), variance = element(initial, "variance"), theta = element(initial, "interaction");
    if (!isNewList(beta) || XLENGTH(beta) != state.groups || !isNewList(theta) || XLENGTH(theta) != state.platforms)
        Rf_error("Invalid initial state dimensions");
    finite_vector(variance, state.groups, "initial variance");
    state.group = (imr_group *)R_alloc(state.groups, sizeof(imr_group));
    memset(state.group, 0, (size_t)state.groups * sizeof(imr_group));
    state.platform = (imr_platform *)R_alloc(state.platforms, sizeof(imr_platform));
    memset(state.platform, 0, (size_t)state.platforms * sizeof(imr_platform));
    int **covered = (int **)R_alloc(state.groups, sizeof(int *));
    state.result = PROTECT(allocVector(VECSXP, 6));
    SEXP names = PROTECT(allocVector(STRSXP, 6));
    const char *labels[] = {"coefficients", "variance", "interaction", "latent", "log_density", "acceptance"};
    for (int i = 0; i < 6; ++i) SET_STRING_ELT(names, i, mkChar(labels[i]));
    setAttrib(state.result, R_NamesSymbol, names);
    SET_VECTOR_ELT(state.result, 0, allocVector(VECSXP, state.groups));
    SET_VECTOR_ELT(state.result, 1, allocMatrix(REALSXP, state.draws, state.groups));
    SET_VECTOR_ELT(state.result, 2, allocVector(VECSXP, state.platforms));
    SET_VECTOR_ELT(state.result, 3, state.keep_latent ? allocVector(VECSXP, state.groups) : R_NilValue);
    SET_VECTOR_ELT(state.result, 4, allocVector(REALSXP, state.draws));
    SET_VECTOR_ELT(state.result, 5, allocVector(REALSXP, 4));
    state.saved_variance = REAL(VECTOR_ELT(state.result, 1)); state.saved_density = REAL(VECTOR_ELT(state.result, 4));
    for (int s = 0; s < state.groups; ++s) {
        imr_group *g = state.group + s; SEXP spec = VECTOR_ELT(groups, s);
        SEXP x = element(spec, "design"), dim = getAttrib(x, R_DimSymbol);
        if (!isReal(x) || !isInteger(dim) || LENGTH(dim) != 2) Rf_error("Invalid design matrix");
        g->n = INTEGER(dim)[0]; g->p = INTEGER(dim)[1];
        if (g->n < 1 || g->p < 1) Rf_error("Empty subgroup design");
        finite_vector(x, (R_xlen_t)g->n * g->p, "design"); g->x = REAL(x);
        SEXP y = element(spec, "response"); finite_vector(y, g->n, "response"); g->observed = REAL(y);
        SEXP scale = element(spec, "prior_scale"); finite_vector(scale, g->p, "prior_scale"); g->scale = REAL(scale);
        g->forced = integer_scalar(element(spec, "forced"), "forced", 1);
        g->outcome = integer_scalar(element(spec, "outcome"), "outcome", 0);
        if (g->forced > g->p || g->outcome > 2) Rf_error("Invalid subgroup specification");
        SEXP residual = element(spec, "residual_prior"); finite_vector(residual, 2, "residual_prior");
        g->residual_shape = REAL(residual)[0]; g->residual_rate = REAL(residual)[1];
        if (g->residual_shape <= 0 || g->residual_rate <= 0) Rf_error("Invalid residual prior");
        SEXP status = element(spec, "status");
        if (g->outcome == 2) {
            if (!isInteger(status) || XLENGTH(status) != g->n) Rf_error("Invalid survival status");
            g->status = INTEGER(status);
        }
        g->variance = REAL(variance)[s]; if (g->variance <= 0) Rf_error("Invalid initial variance");
        SEXP start = VECTOR_ELT(beta, s); finite_vector(start, g->p, "initial coefficients");
        g->beta = (double *)R_alloc(g->p, sizeof(double)); memcpy(g->beta, REAL(start), (size_t)g->p * sizeof(double));
        g->latent = (double *)R_alloc(g->n, sizeof(double));
        g->residual = (double *)R_alloc(g->n, sizeof(double));
        g->predictor = (double *)R_alloc(g->n, sizeof(double));
        g->column_squares = (double *)R_alloc(g->p, sizeof(double));
        covered[s] = (int *)R_alloc(g->p, sizeof(int)); memset(covered[s], 0, (size_t)g->p * sizeof(int));
        for (int j = 0; j < g->p; ++j) {
            if (g->scale[j] <= 0 || (j < g->forced && g->beta[j] == 0)) Rf_error("Invalid prior scale or forced coefficient start");
            double value = 0;
            for (int i = 0; i < g->n; ++i) {
                double x = g->x[i + (R_xlen_t)g->n * j];
                value += x * x;
            }
            if (!R_FINITE(value)) Rf_error("Non-finite design sum of squares; rescale predictors");
            g->column_squares[j] = value;
        }
        for (int i = 0; i < g->n; ++i) {
            if (g->outcome == 1 && g->observed[i] != 0 && g->observed[i] != 1) Rf_error("Invalid binary response");
            if (g->outcome == 2 && g->status[i] != 0 && g->status[i] != 1) Rf_error("Invalid event status");
            g->latent[i] = g->outcome == 1 ? (g->observed[i] == 1 ? .5 : -.5) : g->observed[i];
        }
        SEXP saved = PROTECT(allocMatrix(REALSXP, state.draws, g->p));
        SET_VECTOR_ELT(VECTOR_ELT(state.result, 0), s, saved); g->saved_beta = REAL(saved); UNPROTECT(1);
        if (state.keep_latent && g->outcome != 0) {
            saved = PROTECT(allocMatrix(REALSXP, state.draws, g->n));
            SET_VECTOR_ELT(VECTOR_ELT(state.result, 3), s, saved); g->saved_latent = REAL(saved); UNPROTECT(1);
        }
    }
    for (int l = 0; l < state.platforms; ++l) {
        imr_platform *p = state.platform + l; SEXP spec = VECTOR_ELT(platforms, l);
        SEXP group = element(spec, "groups"), column = element(spec, "columns"), dim = getAttrib(column, R_DimSymbol);
        if (!isInteger(group) || XLENGTH(group) > 16 || !isInteger(column) || !isInteger(dim) || LENGTH(dim) != 2)
            Rf_error("Invalid platform mapping or MRF size (maximum 16 subgroups)");
        p->groups = (int)XLENGTH(group); p->group = INTEGER(group); p->column = INTEGER(column);
        p->features = INTEGER(dim)[1];
        if (INTEGER(dim)[0] != p->groups || p->features < 1 || state.swap_rate * p->features > INT_MAX)
            Rf_error("Invalid platform columns");
        for (int s = 0; s < p->groups; ++s) {
            int index = p->group[s]; if (index < 0 || index >= state.groups) Rf_error("Invalid subgroup index");
            for (int t = 0; t < s; ++t) if (p->group[t] == index) Rf_error("Duplicated subgroup index");
            imr_group *g = state.group + index;
            for (int f = 0; f < p->features; ++f) {
                int j = p->column[s + p->groups * f];
                if (j < g->forced || j >= g->p || covered[index][j]++) Rf_error("Invalid or duplicated feature mapping");
            }
        }
        SEXP nu = element(spec, "nu"); finite_vector(nu, 1, "nu"); p->nu = REAL(nu)[0];
        SEXP interaction = VECTOR_ELT(theta, l); finite_vector(interaction, (R_xlen_t)p->groups * p->groups, "initial interaction");
        p->theta = (double *)R_alloc((size_t)p->groups * p->groups, sizeof(double));
        if (p->groups) memcpy(p->theta, REAL(interaction), (size_t)p->groups * p->groups * sizeof(double));
        for (int i = 0; i < p->groups; ++i) for (int j = 0; j < p->groups; ++j) {
            double value = p->theta[i + p->groups * j];
            if ((i == j && value != 0) || value != p->theta[j + p->groups * i] ||
                (i != j && (state.sharing ? value <= 0 : value != 0))) Rf_error("Invalid initial interaction matrix");
        }
        p->states = 1 << p->groups; p->pairs = state.sharing ? p->groups * (p->groups - 1) / 2 : 0;
        p->pattern = (unsigned char *)R_alloc((size_t)p->states * p->groups, sizeof(unsigned char));
        p->energy = (double *)R_alloc(p->states, sizeof(double));
        p->weights = (double *)R_alloc(p->states, sizeof(double));
        p->cross = (double *)R_alloc(p->groups, sizeof(double));
        p->precision = (double *)R_alloc(p->groups, sizeof(double));
        p->log_bf = (double *)R_alloc(p->groups, sizeof(double));
        for (int k = 0; k < p->states; ++k)
            for (int i = 0; i < p->groups; ++i) p->pattern[k * p->groups + i] = (k >> i) & 1;
        refresh_platform_energy(p);
        SEXP saved = PROTECT(allocMatrix(REALSXP, state.draws, p->pairs));
        SET_VECTOR_ELT(VECTOR_ELT(state.result, 2), l, saved); p->saved_theta = REAL(saved); UNPROTECT(1);
    }
    for (int s = 0; s < state.groups; ++s) for (int j = state.group[s].forced; j < state.group[s].p; ++j)
        if (!covered[s][j]) Rf_error("A molecular column has no platform mapping");
    SEXP answer = R_UnwindProtect(run_chain, &state, release_rng, &state, NULL);
    UNPROTECT(2);
    return answer;
}
