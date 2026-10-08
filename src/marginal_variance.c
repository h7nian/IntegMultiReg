#include "marginal_sampler.h"

typedef struct { double location, scale, df, log_bf; } imr_variance_conditional;

static void fitted_mean(imr_marginal_state *state, imr_marginal_group *g, const double *beta)
{
    SEXP coefficients = PROTECT(Rf_allocVector(REALSXP, g->p));
    memcpy(REAL(coefficients), beta, g->p * sizeof(double));
    SEXP mean = PROTECT(imr_marginal_product(state, g->design, coefficients, 0));
    memcpy(g->mean, REAL(mean), g->n * sizeof(double));
    UNPROTECT(2);
}

static double prior_squares(imr_marginal_group *g, const double *beta)
{
    for (int j = 0; j < g->p; ++j) g->scratch[j] = beta[j] * beta[j] / g->tau[j];
    return imr_marginal_sum(g->scratch, g->p);
}

static int active_count(imr_marginal_group *g, const double *beta)
{
    int count = 0;
    for (int j = 0; j < g->p; ++j) count += beta[j] != 0;
    return count;
}

static imr_variance_conditional conditional(imr_marginal_state *state, imr_marginal_group *g,
                                         int column, int omit_first, int omit_second)
{
    memcpy(g->other, g->beta, g->p * sizeof(double));
    g->other[omit_first] = 0;
    if (omit_second >= 0) g->other[omit_second] = 0;
    fitted_mean(state, g, g->other);
    const double *x = REAL(g->design) + (R_xlen_t)g->n * column;
    double tau = g->tau[column];
    for (int i = 0; i < g->n; ++i) g->scratch[i] = x[i] * x[i];
    double precision = imr_marginal_sum(g->scratch, g->n) + 1 / tau;
    for (int i = 0; i < g->n; ++i) g->scratch[i] = x[i] * (g->response[i] - g->mean[i]);
    double location = imr_marginal_sum(g->scratch, g->n) / precision;
    double shape = g->shape + g->n / 2.0 + 1.5 * active_count(g, g->other);
    for (int i = 0; i < g->n; ++i) {
        double residual = (g->response[i] - g->mean[i]) - x[i] * location;
        g->scratch[i] = residual * residual;
    }
    double residual_squares = imr_marginal_sum(g->scratch, g->n);
    double rate = g->rate + .5 * (residual_squares + prior_squares(g, g->other) + location * location / tau);
    double gain = .5 * precision * (location * location);
    imr_variance_conditional value;
    value.location = location;
    value.scale = sqrt(rate / (precision * (shape + 1)));
    value.df = 2 * shape + 2;
    value.log_bf = -1.5 * log(tau) - .5 * log(precision) + shape * log1p(gain / rate) +
        log(shape * (location * location) / rate + 1 / precision);
    if (!R_FINITE(value.log_bf) || !R_FINITE(value.scale) || value.scale <= 0)
        Rf_error("Non-finite variance-marginal conditional");
    return value;
}

static double draw_coefficient(imr_variance_conditional value)
{
    double scale = fmax(fabs(value.location), value.scale * sqrt(value.df / (value.df - 2)));
    double a = value.location / scale, b = value.scale / scale;
    double mixture = a * a / (a * a + b * b * value.df / (value.df - 2));
    for (int attempt = 0; ; ++attempt) {
        if (attempt % 1024 == 0) R_CheckUserInterrupt();
        double proposal;
        if (unif_rand() < mixture) {
            proposal = rt(value.df);
        } else {
            int sign = unif_rand() < .5 ? -1 : 1;
            double numerator = rchisq(3);
            double denominator = rchisq(value.df - 2);
            proposal = sign * sqrt(value.df * numerator / denominator);
        }
        double shifted = a + b * proposal;
        double acceptance = shifted * shifted / (2 * (a * a + b * b * (proposal * proposal)));
        if (!R_FINITE(acceptance)) Rf_error("Non-finite weighted-t proposal");
        if (unif_rand() < acceptance) {
            double result = value.location + value.scale * proposal;
            if (!R_FINITE(result)) Rf_error("Non-finite coefficient draw");
            if (result != 0) return result;
        }
    }
}

static void update_latent(imr_marginal_state *state, imr_marginal_group *g)
{
    fitted_mean(state, g, g->beta);
    double shape = g->shape + g->n / 2.0 + 1.5 * active_count(g, g->beta);
    double prior = prior_squares(g, g->beta);
    for (int i = 0; i < g->n; ++i) {
        if (!g->outcome || (g->outcome == 2 && g->status[i])) continue;
        long double squares = 0;
        for (int k = 0; k < g->n; ++k) if (k != i) {
            double residual = g->response[k] - g->mean[k];
            double square = residual * residual;
            squares += square;
        }
        double rate = g->rate + .5 * ((double)squares + prior);
        double scale = sqrt(rate / (shape - .5)), df = 2 * shape - 1;
        int sign = g->outcome == 1 && g->observed[i] == 0 ? -1 : 1;
        double lower = g->outcome == 1 ? 0 : g->observed[i];
        double location = sign * g->mean[i];
        double tail = pt((lower - location) / scale, df, 0, 1);
        double value = location + scale * qt(log(unif_rand()) + tail, df, 0, 1);
        if (!R_FINITE(value)) Rf_error("Non-finite latent draw");
        g->response[i] = sign * fmax(lower, value);
    }
}

static double variance_rate(imr_marginal_state *state, imr_marginal_group *g)
{
    fitted_mean(state, g, g->beta);
    for (int i = 0; i < g->n; ++i) {
        double residual = g->response[i] - g->mean[i];
        g->scratch[i] = residual * residual;
    }
    double squares = imr_marginal_sum(g->scratch, g->n);
    return g->rate + .5 * (squares + prior_squares(g, g->beta));
}

void imr_variance_run(imr_marginal_state *state)
{
    imr_variance_conditional parameters[IMR_MARGINAL_MAX_SUBGROUPS], proposed[IMR_MARGINAL_MAX_SUBGROUPS];
    int old_columns[IMR_MARGINAL_MAX_SUBGROUPS], new_columns[IMR_MARGINAL_MAX_SUBGROUPS];
    double log_bf[IMR_MARGINAL_MAX_SUBGROUPS];
    int total = state->burnin + state->draws * state->thin;
    for (int sweep = 0; sweep < total; ++sweep) {
        int iteration = sweep + 1;
        R_CheckUserInterrupt();
        for (int s = 0; s < state->groups_count; ++s) {
            imr_marginal_group *g = state->groups + s;
            update_latent(state, g);
            for (int j = 0; j < g->forced; ++j) g->beta[j] = draw_coefficient(conditional(state, g, j, j, -1));
        }
        for (int l = 0; l < state->platforms_count; ++l) {
            imr_marginal_platform *p = state->platforms + l;
            for (int f = 0; f < p->features; ++f) {
                for (int i = 0; i < p->members; ++i) {
                    int j = p->columns[i + p->members * f];
                    parameters[i] = conditional(state, state->groups + p->groups[i], j, j, -1);
                    log_bf[i] = parameters[i].log_bf;
                }
                int pattern = imr_marginal_select(state, p, log_bf);
                for (int i = 0; i < p->members; ++i) {
                    imr_marginal_group *g = state->groups + p->groups[i];
                    int j = p->columns[i + p->members * f];
                    g->selected[j] = REAL(p->patterns)[pattern + p->patterns_count * i];
                    g->beta[j] = g->selected[j] ? draw_coefficient(parameters[i]) : 0;
                }
            }
            int attempts = p->features >= 2 ? (int)ceil(state->swap_rate * p->features) : 0;
            for (int attempt = 0; attempt < attempts; ++attempt) {
                int pair[2], moves = 0;
                imr_marginal_pair(p, pair);
                double ratio = 0;
                for (int i = 0; i < p->members; ++i) {
                    imr_marginal_group *g = state->groups + p->groups[i];
                    int first = p->columns[i + p->members * pair[0]], second = p->columns[i + p->members * pair[1]];
                    old_columns[i] = -1;
                    if (g->selected[first] == g->selected[second]) continue;
                    int old = g->selected[first] ? first : second, next = g->selected[first] ? second : first;
                    imr_variance_conditional from = conditional(state, g, old, first, second);
                    proposed[i] = conditional(state, g, next, first, second);
                    old_columns[i] = old; new_columns[i] = next; moves++;
                    ratio = ratio + proposed[i].log_bf - from.log_bf;
                }
                if (moves) {
                    state->acceptance[IMR_SWAP_PROPOSALS]++;
                    if (imr_marginal_accept(ratio)) {
                        for (int i = 0; i < p->members; ++i) if (old_columns[i] >= 0) {
                            imr_marginal_group *g = state->groups + p->groups[i];
                            g->selected[old_columns[i]] = 0; g->beta[old_columns[i]] = 0;
                            g->selected[new_columns[i]] = 1; g->beta[new_columns[i]] = draw_coefficient(proposed[i]);
                        }
                        state->acceptance[IMR_SWAP_ACCEPTS]++;
                    }
                }
            }
            imr_marginal_interaction(state, p);
        }
        /* Recovery consumes the same draws on every sweep, including warmup. */
        for (int s = 0; s < state->groups_count; ++s) {
            imr_marginal_group *g = state->groups + s;
            double shape = g->shape + g->n / 2.0 + 1.5 * active_count(g, g->beta);
            g->variance = 1 / rgamma(shape, 1 / variance_rate(state, g));
        }
        if (iteration > state->burnin && (iteration - state->burnin) % state->thin == 0) {
            double score = 0;
            for (int s = 0; s < state->groups_count; ++s) {
                imr_marginal_group *g = state->groups + s;
                int active = active_count(g, g->beta);
                double shape = g->shape + g->n / 2.0 + 1.5 * active;
                double rate = variance_rate(state, g);
                long double scales = 0, coefficients = 0;
                for (int j = 0; j < g->p; ++j) if (g->beta[j] != 0) {
                    scales += log(g->tau[j]); coefficients += log(fabs(g->beta[j]));
                }
                score = score + g->shape * log(g->rate) - lgammafn(g->shape) -
                    (g->n + active) / 2.0 * log(2 * M_PI) - 1.5 * (double)scales +
                    2 * (double)coefficients + lgammafn(shape) - shape * log(rate);
            }
            imr_marginal_store(state, (iteration - state->burnin) / state->thin - 1, score);
        }
        imr_marginal_progress(state, iteration);
    }
}
