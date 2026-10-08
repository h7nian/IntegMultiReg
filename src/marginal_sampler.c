#include "marginal_sampler.h"

SEXP imr_marginal_element(SEXP list, const char *name)
{
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    for (int i = 0; i < Rf_length(names); ++i)
        if (!strcmp(CHAR(STRING_ELT(names, i)), name)) return VECTOR_ELT(list, i);
    Rf_error("Missing marginal sampler field: %s", name);
    return R_NilValue;
}

/* Keep R's matrix-product policy (including options(matprod)) and BLAS.
   Scalar updates and traversal live in C; these calls perform dense algebra. */
SEXP imr_marginal_product(imr_marginal_state *state, SEXP x, SEXP y, int cross)
{
    SEXP call = PROTECT(Rf_lang3(cross ? state->crossprod : state->multiply, x, y));
    SEXP value = Rf_eval(call, R_BaseEnv);
    UNPROTECT(1);
    return value;
}

double imr_marginal_sum(const double *values, int n)
{
    long double sum = 0;
    for (int i = 0; i < n; ++i) sum += values[i];
    return (double)sum;
}

double imr_marginal_log_sum_exp(double *values, int n)
{
    double maximum = R_NegInf;
    for (int i = 0; i < n; ++i) if (values[i] > maximum) maximum = values[i];
    long double sum = 0;
    for (int i = 0; i < n; ++i) sum += exp(values[i] - maximum);
    return maximum + log((double)sum);
}

void imr_marginal_energy(imr_marginal_state *state, imr_marginal_platform *p, double *out)
{
    SEXP product = PROTECT(imr_marginal_product(state, p->patterns, p->theta, 0));
    for (int k = 0; k < p->patterns_count; ++k) {
        long double interaction = 0;
        int included = 0;
        for (int i = 0; i < p->members; ++i) {
            int index = k + p->patterns_count * i;
            included += REAL(p->patterns)[index];
            double term = REAL(product)[index] * REAL(p->patterns)[index];
            interaction += term;
        }
        out[k] = included * p->nu + (double)interaction;
    }
    UNPROTECT(1);
}

int imr_marginal_select(imr_marginal_state *state, imr_marginal_platform *p, double *log_bf)
{
    SEXP factors = PROTECT(Rf_allocVector(REALSXP, p->members));
    memcpy(REAL(factors), log_bf, p->members * sizeof(double));
    SEXP product = PROTECT(imr_marginal_product(state, p->patterns, factors, 0));
    imr_marginal_energy(state, p, p->energy);
    double maximum = R_NegInf;
    for (int k = 0; k < p->patterns_count; ++k) {
        p->weights[k] = p->energy[k] + REAL(product)[k];
        if (p->weights[k] > maximum) maximum = p->weights[k];
    }
    /* Match sample.int(size = 1, prob = ...), including its sorting of ties
       and double-precision probability normalization. */
    double total = 0;
    for (int k = 0; k < p->patterns_count; ++k) {
        p->weights[k] = exp(p->weights[k] - maximum);
        total += p->weights[k];
        p->permutation[k] = k;
    }
    if (!R_FINITE(total) || total <= 0) Rf_error("Invalid marginal selection weights");
    for (int k = 0; k < p->patterns_count; ++k) p->weights[k] /= total;
    revsort(p->weights, p->permutation, p->patterns_count);
    double uniform = unif_rand(), cumulative = 0;
    int selected = p->patterns_count - 1;
    for (int k = 0; k < p->patterns_count - 1; ++k) {
        cumulative += p->weights[k];
        if (uniform <= cumulative) { selected = k; break; }
    }
    UNPROTECT(2);
    return p->permutation[selected];
}

void imr_marginal_pair(imr_marginal_platform *p, int *pair)
{
    for (int i = 0; i < p->features; ++i) p->pair_pool[i] = i;
    int remaining = p->features;
    for (int i = 0; i < 2; ++i) {
        int index = (int)R_unif_index(remaining);
        pair[i] = p->pair_pool[index];
        p->pair_pool[index] = p->pair_pool[--remaining];
    }
}

int imr_marginal_accept(double log_ratio)
{
    if (ISNAN(log_ratio)) Rf_error("Undefined marginal Metropolis ratio");
    return log(unif_rand()) < fmin(0, log_ratio);
}

static void refresh_selection(imr_marginal_state *state, imr_marginal_platform *p)
{
    for (int f = 0; f < p->features; ++f) for (int i = 0; i < p->members; ++i) {
        imr_marginal_group *g = state->groups + p->groups[i];
        int j = p->columns[i + p->members * f];
        REAL(p->selected_matrix)[i + p->members * f] = g->selected[j];
    }
}

void imr_marginal_interaction(imr_marginal_state *state, imr_marginal_platform *p)
{
    if (!state->sharing) return;
    refresh_selection(state, p);
    double *theta = REAL(p->theta);
    for (int j = 1; j < p->members; ++j) for (int i = 0; i < j; ++i) {
        double old = theta[i + p->members * j];
        double proposal = exp(log(old) + rnorm(0, state->theta_step));
        state->acceptance[2]++;
        if (!R_FINITE(proposal) || proposal <= 0) continue;
        int shared = 0;
        for (int f = 0; f < p->features; ++f)
            shared += REAL(p->selected_matrix)[i + p->members * f] *
                REAL(p->selected_matrix)[j + p->members * f];
        theta[i + p->members * j] = theta[j + p->members * i] = proposal;
        imr_marginal_energy(state, p, p->energy);
        double proposed_normalizer = imr_marginal_log_sum_exp(p->energy, p->patterns_count);
        theta[i + p->members * j] = theta[j + p->members * i] = old;
        imr_marginal_energy(state, p, p->energy);
        double old_normalizer = imr_marginal_log_sum_exp(p->energy, p->patterns_count);
        double ratio = 2 * (proposal - old) * shared - p->features *
            (proposed_normalizer - old_normalizer) +
            dgamma(proposal, state->interaction_shape, 1 / state->interaction_rate, 1) -
            dgamma(old, state->interaction_shape, 1 / state->interaction_rate, 1) + log(proposal) - log(old);
        if (imr_marginal_accept(ratio)) {
            theta[i + p->members * j] = theta[j + p->members * i] = proposal;
            state->acceptance[3]++;
        }
    }
}

void imr_marginal_store(imr_marginal_state *state, int row, double score)
{
    for (int s = 0; s < state->groups_count; ++s) {
        imr_marginal_group *g = state->groups + s;
        for (int j = 0; j < g->p; ++j) g->coefficient_draws[row + (R_xlen_t)state->draws * j] = g->beta[j];
        g->variance_draws[row] = g->variance;
        if (g->latent_draws) for (int i = 0; i < g->n; ++i)
            g->latent_draws[row + (R_xlen_t)state->draws * i] = g->response[i];
    }
    for (int l = 0; l < state->platforms_count; ++l) {
        imr_marginal_platform *p = state->platforms + l;
        refresh_selection(state, p);
        long double baseline = 0;
        for (int f = 0; f < p->features; ++f) {
            int count = 0;
            for (int i = 0; i < p->members; ++i) count += REAL(p->selected_matrix)[i + p->members * f];
            baseline += count * p->nu;
        }
        SEXP product = PROTECT(imr_marginal_product(state, p->theta, p->selected_matrix, 0));
        long double interaction = 0;
        for (int i = 0; i < p->members * p->features; ++i) {
            double term = REAL(p->selected_matrix)[i] * REAL(product)[i];
            interaction += term;
        }
        imr_marginal_energy(state, p, p->energy);
        score = score + (double)baseline + (double)interaction -
            p->features * imr_marginal_log_sum_exp(p->energy, p->patterns_count);
        if (state->sharing) {
            long double prior = 0;
            int column = 0;
            for (int j = 1; j < p->members; ++j) for (int i = 0; i < j; ++i) {
                double value = REAL(p->theta)[i + p->members * j];
                p->theta_draws[row + (R_xlen_t)state->draws * column++] = value;
                prior += dgamma(value, state->interaction_shape, 1 / state->interaction_rate, 1);
            }
            score += (double)prior;
        }
        UNPROTECT(1);
    }
    state->log_density[row] = score;
}

void imr_marginal_progress(imr_marginal_state *state, int iteration)
{
    int total = state->burnin + state->draws * state->thin;
    int every = total / 10;
    if (every < 1) every = 1;
    if (state->verbose && (iteration % every == 0 || iteration == total))
        Rprintf("Iteration %d of %d\n", iteration, total);
}

static double scalar(SEXP x, const char *name)
{
    SEXP value = imr_marginal_element(x, name);
    if (Rf_length(value) != 1 || !Rf_isNumeric(value) || !R_FINITE(Rf_asReal(value)))
        Rf_error("Invalid marginal setting: %s", name);
    return Rf_asReal(value);
}

static int integer_setting(SEXP x, const char *name, int minimum, int maximum)
{
    double value = scalar(x, name);
    if (value < minimum || value > maximum || value != floor(value))
        Rf_error("Invalid marginal setting: %s", name);
    return (int)value;
}

static double *numeric_data(SEXP value, R_xlen_t size, const char *name)
{
    if (TYPEOF(value) != REALSXP || XLENGTH(value) != size) Rf_error("Invalid marginal %s", name);
    double *data = REAL(value);
    for (R_xlen_t i = 0; i < size; ++i) if (!R_FINITE(data[i])) Rf_error("Non-finite marginal %s", name);
    return data;
}

static SEXP run(void *data)
{
    imr_marginal_state *state = data;
    if (state->coefficients_marginalized) imr_coefficients_run(state); else imr_variance_run(state);
    return state->result;
}

static void finish(void *data) { (void)data; PutRNGstate(); }

SEXP imr_marginal_sample(SEXP groups, SEXP platforms, SEXP settings, SEXP initial, SEXP algebra)
{
    imr_marginal_state state = {0};
    state.groups_count = Rf_length(groups);
    state.platforms_count = Rf_length(platforms);
    if (!Rf_isNewList(groups) || !Rf_isNewList(platforms) || !state.groups_count || !state.platforms_count)
        Rf_error("Invalid marginal sampler specification");
    state.draws = integer_setting(settings, "draws", 1, INT_MAX);
    state.burnin = integer_setting(settings, "burnin", 0, INT_MAX);
    state.thin = integer_setting(settings, "thin", 1, INT_MAX);
    state.sharing = integer_setting(settings, "sharing", 0, 1);
    state.verbose = integer_setting(settings, "verbose", 0, 1);
    state.coefficients_marginalized = integer_setting(settings, "coefficients_marginalized", 0, 1);
    if (state.draws < 1 || state.burnin < 0 || state.thin < 1 ||
        (double)state.burnin + (double)state.draws * state.thin > INT_MAX)
        Rf_error("Invalid marginal iteration budget");
    state.theta_step = scalar(settings, "theta_step"); state.swap_rate = scalar(settings, "swap_rate");
    state.interaction_shape = scalar(settings, "interaction_shape");
    state.interaction_rate = scalar(settings, "interaction_rate");
    if (state.theta_step <= 0 || state.swap_rate < 0 || state.swap_rate > 10 ||
        state.interaction_shape <= 0 || state.interaction_rate <= 0)
        Rf_error("Invalid marginal proposal or prior settings");
    int keep_latent = integer_setting(settings, "keep_latent", 0, 1);
    double *variance_steps = numeric_data(imr_marginal_element(settings, "variance_step"), state.groups_count, "variance steps");
    double *nu = numeric_data(imr_marginal_element(settings, "nu"), state.platforms_count, "inclusion prior");
    state.multiply = imr_marginal_element(algebra, "multiply");
    state.crossprod = imr_marginal_element(algebra, "crossprod");
    state.largest_eigenvalue = imr_marginal_element(algebra, "largest_eigenvalue");
    state.prepare_model = imr_marginal_element(algebra, "prepare_model");
    state.roots = PROTECT(Rf_allocVector(VECSXP, 4));
    SEXP starts = PROTECT(Rf_duplicate(imr_marginal_element(initial, "coefficients")));
    SEXP interactions = PROTECT(Rf_duplicate(imr_marginal_element(initial, "interaction")));
    SEXP variance = imr_marginal_element(initial, "variance");
    numeric_data(variance, state.groups_count, "initial variance");
    if (Rf_length(starts) != state.groups_count || Rf_length(interactions) != state.platforms_count)
        Rf_error("Invalid marginal initial dimensions");
    SET_VECTOR_ELT(state.roots, 0, starts); SET_VECTOR_ELT(state.roots, 1, interactions);
    SET_VECTOR_ELT(state.roots, 2, Rf_allocVector(VECSXP, state.platforms_count));
    SET_VECTOR_ELT(state.roots, 3, Rf_allocVector(VECSXP, state.groups_count));
    state.result = PROTECT(Rf_allocVector(VECSXP, 6));
    const char *names[] = {"coefficients", "variance", "interaction", "latent", "log_density", "acceptance"};
    SEXP result_names = PROTECT(Rf_allocVector(STRSXP, 6));
    for (int i = 0; i < 6; ++i) SET_STRING_ELT(result_names, i, Rf_mkChar(names[i]));
    Rf_setAttrib(state.result, R_NamesSymbol, result_names);
    SET_VECTOR_ELT(state.result, 0, Rf_allocVector(VECSXP, state.groups_count));
    Rf_setAttrib(VECTOR_ELT(state.result, 0), R_NamesSymbol, Rf_getAttrib(groups, R_NamesSymbol));
    SET_VECTOR_ELT(state.result, 1, Rf_allocMatrix(REALSXP, state.draws, state.groups_count));
    SET_VECTOR_ELT(state.result, 2, Rf_allocVector(VECSXP, state.platforms_count));
    if (keep_latent) {
        SET_VECTOR_ELT(state.result, 3, Rf_allocVector(VECSXP, state.groups_count));
        Rf_setAttrib(VECTOR_ELT(state.result, 3), R_NamesSymbol, Rf_getAttrib(groups, R_NamesSymbol));
    }
    SET_VECTOR_ELT(state.result, 4, Rf_allocVector(REALSXP, state.draws));
    SET_VECTOR_ELT(state.result, 5, Rf_allocVector(REALSXP, state.coefficients_marginalized ? 8 : 4));
    state.log_density = REAL(VECTOR_ELT(state.result, 4));
    state.acceptance = REAL(VECTOR_ELT(state.result, 5));
    memset(state.acceptance, 0, (state.coefficients_marginalized ? 8 : 4) * sizeof(double));
    const char *acceptance_names[] = {"swap_proposals", "swap_accepts", "interaction_proposals", "interaction_accepts",
        "variance_proposals", "variance_accepts", "latent_proposals", "latent_accepts"};
    SEXP labels = PROTECT(Rf_allocVector(STRSXP, state.coefficients_marginalized ? 8 : 4));
    for (int i = 0; i < Rf_length(labels); ++i) SET_STRING_ELT(labels, i, Rf_mkChar(acceptance_names[i]));
    Rf_setAttrib(VECTOR_ELT(state.result, 5), R_NamesSymbol, labels);
    state.groups = (imr_marginal_group *)R_alloc(state.groups_count, sizeof(imr_marginal_group));
    memset(state.groups, 0, state.groups_count * sizeof(imr_marginal_group));
    for (int s = 0; s < state.groups_count; ++s) {
        imr_marginal_group *g = state.groups + s;
        SEXP source = VECTOR_ELT(groups, s);
        g->design = imr_marginal_element(source, "X");
        SEXP dim = Rf_getAttrib(g->design, R_DimSymbol);
        if (Rf_length(dim) != 2) Rf_error("Invalid marginal design dimensions");
        g->n = INTEGER(dim)[0]; g->p = INTEGER(dim)[1];
        numeric_data(g->design, (R_xlen_t)g->n * g->p, "design");
        g->forced = integer_setting(source, "n_forced", 1, g->p);
        if (g->n < 1 || g->forced < 1 || g->forced > g->p) Rf_error("Invalid marginal group dimensions");
        const char *outcome = CHAR(Rf_asChar(imr_marginal_element(source, "outcome_type")));
        if (strcmp(outcome, "continuous") && strcmp(outcome, "binary") && strcmp(outcome, "right.censored"))
            Rf_error("Invalid marginal outcome type");
        g->outcome = !strcmp(outcome, "continuous") ? 0 : !strcmp(outcome, "binary") ? 1 : 2;
        g->observed = numeric_data(imr_marginal_element(source, "y"), g->n, "response");
        g->tau = numeric_data(imr_marginal_element(source, "prior_scale"), g->p, "prior scale");
        double *prior = numeric_data(imr_marginal_element(source, "residual_prior"), 2, "variance prior");
        g->shape = prior[0]; g->rate = prior[1];
        if (g->shape <= 0 || g->rate <= 0) Rf_error("Invalid marginal variance prior");
        g->beta = numeric_data(VECTOR_ELT(starts, s), g->p, "initial coefficients");
        g->variance = REAL(variance)[s];
        if (g->variance <= 0) Rf_error("Invalid marginal initial variance");
        g->variance_step = variance_steps[s];
        if (g->variance_step <= 0) Rf_error("Invalid marginal variance step");
        g->response = (double *)R_alloc(g->n, sizeof(double));
        g->mean = (double *)R_alloc(g->n, sizeof(double));
        g->other = (double *)R_alloc(g->p, sizeof(double));
        g->scratch = (double *)R_alloc(g->n + g->p, sizeof(double));
        g->selected = (int *)R_alloc(g->p, sizeof(int));
        g->status = (int *)R_alloc(g->n, sizeof(int));
        SEXP status = imr_marginal_element(source, "status");
        if (g->outcome == 2 && (TYPEOF(status) != INTSXP || Rf_length(status) != g->n))
            Rf_error("Invalid marginal event indicators");
        for (int j = 0; j < g->p; ++j) {
            if (g->tau[j] <= 0 || (j < g->forced && g->beta[j] == 0)) Rf_error("Invalid marginal coefficient state");
            g->selected[j] = g->beta[j] != 0;
        }
        for (int i = 0; i < g->n; ++i) {
            g->response[i] = g->outcome == 1 ? (g->observed[i] == 1 ? .5 : -.5) : g->observed[i];
            g->status[i] = g->outcome == 2 ? INTEGER(status)[i] : 1;
        }
        SEXP saved = PROTECT(Rf_allocMatrix(REALSXP, state.draws, g->p));
        SET_VECTOR_ELT(VECTOR_ELT(state.result, 0), s, saved); g->coefficient_draws = REAL(saved);
        g->variance_draws = REAL(VECTOR_ELT(state.result, 1)) + (R_xlen_t)state.draws * s;
        if (keep_latent && g->outcome) {
            SEXP latent = PROTECT(Rf_allocMatrix(REALSXP, state.draws, g->n));
            SET_VECTOR_ELT(VECTOR_ELT(state.result, 3), s, latent); g->latent_draws = REAL(latent);
            UNPROTECT(1);
        }
        UNPROTECT(1);
    }
    state.platforms = (imr_marginal_platform *)R_alloc(state.platforms_count, sizeof(imr_marginal_platform));
    for (int l = 0; l < state.platforms_count; ++l) {
        imr_marginal_platform *p = state.platforms + l;
        SEXP source = VECTOR_ELT(platforms, l);
        SEXP members = imr_marginal_element(source, "members"), columns = imr_marginal_element(source, "columns");
        p->members = Rf_length(members); p->features = Rf_length(imr_marginal_element(source, "features"));
        if (p->members > 7 || TYPEOF(members) != INTSXP || TYPEOF(columns) != INTSXP ||
            Rf_length(columns) != p->members * p->features) Rf_error("Invalid marginal platform mapping");
        p->groups = (int *)R_alloc(p->members, sizeof(int));
        p->columns = (int *)R_alloc(p->members * p->features, sizeof(int));
        for (int i = 0; i < p->members; ++i) {
            p->groups[i] = INTEGER(members)[i] - 1;
            if (p->groups[i] < 0 || p->groups[i] >= state.groups_count) Rf_error("Invalid marginal subgroup index");
            imr_marginal_group *g = state.groups + p->groups[i];
            for (int f = 0; f < p->features; ++f) {
                int index = i + p->members * f;
                p->columns[index] = INTEGER(columns)[index] - 1;
                if (p->columns[index] < g->forced || p->columns[index] >= g->p) Rf_error("Invalid marginal feature index");
            }
        }
        p->nu = nu[l];
        p->patterns_count = 1 << p->members;
        p->patterns = imr_marginal_element(source, "patterns");
        numeric_data(p->patterns, p->patterns_count * p->members, "selection patterns");
        p->theta = VECTOR_ELT(interactions, l);
        numeric_data(p->theta, p->members * p->members, "initial interactions");
        p->selected_matrix = Rf_allocMatrix(REALSXP, p->members, p->features);
        SET_VECTOR_ELT(VECTOR_ELT(state.roots, 2), l, p->selected_matrix);
        p->energy = (double *)R_alloc(p->patterns_count, sizeof(double));
        p->weights = (double *)R_alloc(p->patterns_count, sizeof(double));
        p->permutation = (int *)R_alloc(p->patterns_count, sizeof(int));
        p->pair_pool = (int *)R_alloc(p->features, sizeof(int));
        SEXP saved = PROTECT(Rf_allocMatrix(REALSXP, state.draws, state.sharing ? p->members * (p->members - 1) / 2 : 0));
        SET_VECTOR_ELT(VECTOR_ELT(state.result, 2), l, saved); p->theta_draws = REAL(saved);
        UNPROTECT(1);
    }
    GetRNGstate();
    SEXP result = R_ExecWithCleanup(run, &state, finish, NULL);
    UNPROTECT(6);
    return result;
}
