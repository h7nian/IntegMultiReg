#include "marginal_sampler.h"

typedef struct imr_coefficient_model {
    struct imr_coefficient_model *next;
    int d, nodes_count, *active, *selection;
    SEXP design, covariance, normal_factor;
    double *tau, *nodes, *log_weights, *latent_precision, *moment_terms;
    double H, log_D;
} imr_coefficient_model;

typedef struct {
    imr_coefficient_model *model;
    double *mean, B;
} imr_coefficient_density;

typedef struct {
    imr_marginal_state *state;
    imr_coefficient_model **models;
    double **first_mean, **second_mean, **trial_response;
    int **trial_selection;
} imr_coefficient_workspace;

static imr_coefficient_model *model_for(imr_coefficient_workspace *work, int s, const int *selection)
{
    imr_marginal_state *state = work->state;
    imr_marginal_group *g = state->groups + s;
    for (imr_coefficient_model *model = work->models[s]; model; model = model->next)
        if (!memcmp(model->selection, selection, g->p * sizeof(int))) return model;
    SEXP group_index = PROTECT(Rf_ScalarInteger(s + 1));
    SEXP selected = PROTECT(Rf_allocVector(INTSXP, g->p - g->forced));
    for (int j = g->forced; j < g->p; ++j) INTEGER(selected)[j - g->forced] = selection[j];
    SEXP call = PROTECT(Rf_lang3(state->prepare_model, group_index, selected));
    SEXP prepared = PROTECT(Rf_eval(call, R_BaseEnv));
    SEXP roots = VECTOR_ELT(state->roots, 3);
    SEXP saved = PROTECT(Rf_cons(prepared, VECTOR_ELT(roots, s)));
    SET_VECTOR_ELT(roots, s, saved);
    imr_coefficient_model *model = (imr_coefficient_model *)R_alloc(1, sizeof(imr_coefficient_model));
    model->next = work->models[s]; work->models[s] = model;
    model->selection = (int *)R_alloc(g->p, sizeof(int));
    memcpy(model->selection, selection, g->p * sizeof(int));
    SEXP active = imr_marginal_element(prepared, "active");
    model->d = Rf_length(active);
    model->active = (int *)R_alloc(model->d, sizeof(int));
    for (int j = 0; j < model->d; ++j) model->active[j] = INTEGER(active)[j] - 1;
    model->design = imr_marginal_element(prepared, "X");
    model->covariance = imr_marginal_element(prepared, "covariance");
    model->normal_factor = imr_marginal_element(prepared, "normal_factor");
    model->tau = REAL(imr_marginal_element(prepared, "tau"));
    SEXP nodes = imr_marginal_element(prepared, "nodes");
    model->nodes_count = INTEGER(Rf_getAttrib(nodes, R_DimSymbol))[0];
    model->nodes = REAL(nodes);
    model->log_weights = REAL(imr_marginal_element(prepared, "log_weight"));
    model->moment_terms = (double *)R_alloc(model->nodes_count, sizeof(double));
    model->latent_precision = g->outcome ? REAL(imr_marginal_element(prepared, "latent_precision")) : NULL;
    model->H = Rf_asReal(imr_marginal_element(prepared, "H"));
    model->log_D = Rf_asReal(imr_marginal_element(prepared, "log_D"));
    UNPROTECT(5);
    return model;
}

static imr_coefficient_density statistics(imr_coefficient_workspace *work, int s, const int *selected,
                                        const double *response, double *mean)
{
    imr_marginal_state *state = work->state;
    imr_marginal_group *g = state->groups + s;
    imr_coefficient_model *model = model_for(work, s, selected);
    SEXP y = PROTECT(Rf_allocVector(REALSXP, g->n));
    memcpy(REAL(y), response, g->n * sizeof(double));
    SEXP xy = PROTECT(imr_marginal_product(state, model->design, y, 1));
    SEXP coefficient_mean = PROTECT(imr_marginal_product(state, model->covariance, xy, 0));
    memcpy(mean, REAL(coefficient_mean), model->d * sizeof(double));
    SEXP fitted = PROTECT(imr_marginal_product(state, model->design, coefficient_mean, 0));
    long double residual_squares = 0, prior_squares = 0;
    for (int i = 0; i < g->n; ++i) {
        double residual = response[i] - REAL(fitted)[i];
        double square = residual * residual;
        residual_squares += square;
    }
    for (int j = 0; j < model->d; ++j) {
        double term = mean[j] * mean[j] / model->tau[j];
        prior_squares += term;
    }
    imr_coefficient_density value = {model, mean,
        g->rate + .5 * ((double)residual_squares + (double)prior_squares)};
    UNPROTECT(4);
    return value;
}

static double log_kernel(imr_coefficient_density density, double variance)
{
    imr_coefficient_model *model = density.model;
    double scale = sqrt(variance);
    for (int i = 0; i < model->nodes_count; ++i) {
        if (i % 16384 == 0) R_CheckUserInterrupt();
        long double sum = 0;
        for (int j = 0; j < model->d; ++j) {
            double value = scale * model->nodes[i + (R_xlen_t)model->nodes_count * j] + density.mean[j];
            sum += log(fabs(value));
        }
        model->moment_terms[i] = model->log_weights[i] + 2 * (double)sum;
    }
    double moment = imr_marginal_log_sum_exp(model->moment_terms, model->nodes_count);
    double value = model->log_D - (model->H + 1) * log(variance) - density.B / variance + moment;
    if (!R_FINITE(value)) Rf_error("Non-finite exact coefficient-integrated density");
    return value;
}

static void restore_coefficients(imr_marginal_state *state, imr_marginal_group *g, imr_coefficient_density density)
{
    imr_coefficient_model *model = density.model;
    int d = model->d;
    double *scales = g->other;
    long double sum = 0;
    for (int j = 0; j < d; ++j) {
        scales[j] = sqrt(density.mean[j] * density.mean[j] + g->variance * REAL(model->covariance)[j + d * j]);
        double ratio = density.mean[j] / scales[j];
        double square = ratio * ratio;
        sum += square;
    }
    double mean_norm = sqrt((double)sum);
    SEXP scaled_covariance = PROTECT(Rf_allocMatrix(REALSXP, d, d));
    for (int j = 0; j < d; ++j) for (int i = 0; i < d; ++i)
        REAL(scaled_covariance)[i + d * j] = (g->variance * REAL(model->covariance)[i + d * j]) / (scales[i] * scales[j]);
    SEXP call = PROTECT(Rf_lang2(state->largest_eigenvalue, scaled_covariance));
    SEXP eigenvalue = PROTECT(Rf_eval(call, R_BaseEnv));
    double largest = Rf_asReal(eigenvalue) * (1 + 64 * DBL_EPSILON * d);
    double a = .25, b = a * mean_norm / sqrt(largest);
    double radius = 2 * d / (b + sqrt(b * b + 4 * a * d));
    double bound = 2 * d * log(mean_norm + sqrt(largest) * radius) -
        d * log((double)d) - a * (radius * radius) + 64 * DBL_EPSILON * d;
    SEXP normal = PROTECT(Rf_allocMatrix(REALSXP, 16, d));
    SEXP proposed = PROTECT(Rf_allocMatrix(REALSXP, 16, d));
    for (;;) {
        R_CheckUserInterrupt();
        /* Draw the entire proposal batch, then all 16 uniforms, even if an
           early proposal succeeds. This is the reference's RNG contract. */
        for (int i = 0; i < 16 * d; ++i) REAL(normal)[i] = rnorm(0, 1);
        SEXP product = PROTECT(imr_marginal_product(state, normal, model->normal_factor, 0));
        double scale = sqrt(2 * g->variance);
        for (int j = 0; j < d; ++j) for (int i = 0; i < 16; ++i)
            REAL(proposed)[i + 16 * j] = scale * REAL(product)[i + 16 * j] + density.mean[j];
        double log_acceptance[16];
        for (int i = 0; i < 16; ++i) {
            long double logs = 0, squares = 0;
            for (int j = 0; j < d; ++j) {
                logs += log(fabs(REAL(proposed)[i + 16 * j] / scales[j]));
                double value = REAL(normal)[i + 16 * j];
                double square = value * value;
                squares += square;
            }
            log_acceptance[i] = 2 * (double)logs - .5 * (double)squares - bound;
            if (ISNAN(log_acceptance[i]) || log_acceptance[i] > 1e-10)
                Rf_error("Invalid coefficient rejection envelope");
        }
        int accepted = -1;
        for (int i = 0; i < 16; ++i) {
            double uniform = unif_rand();
            if (log(uniform) < fmin(0, log_acceptance[i]) && accepted < 0) accepted = i;
        }
        UNPROTECT(1);
        if (accepted >= 0) {
            memset(g->beta, 0, g->p * sizeof(double));
            for (int j = 0; j < d; ++j) g->beta[model->active[j]] = REAL(proposed)[accepted + 16 * j];
            break;
        }
    }
    UNPROTECT(5);
}

static void update_latent(imr_coefficient_workspace *work, int s)
{
    imr_marginal_state *state = work->state;
    imr_marginal_group *g = state->groups + s;
    for (int i = 0; i < g->n; ++i) {
        if (!g->outcome || (g->outcome == 2 && g->status[i])) continue;
        imr_coefficient_density current = statistics(work, s, g->selected, g->response, work->first_mean[s]);
        double *precision = current.model->latent_precision;
        long double sum = 0;
        for (int j = 0; j < g->n; ++j) if (j != i) {
            double term = precision[i + g->n * j] * g->response[j];
            sum += term;
        }
        double location = -(double)sum / precision[i + g->n * i];
        double scale = sqrt(g->variance / precision[i + g->n * i]);
        int sign = g->outcome == 1 && g->observed[i] == 0 ? -1 : 1;
        double lower = g->outcome == 1 ? 0 : g->observed[i];
        double tail = pnorm5((lower - sign * location) / scale, 0, 1, 0, 1);
        double proposed = sign * location + scale * qnorm5(log(unif_rand()) + tail, 0, 1, 0, 1);
        if (!R_FINITE(proposed)) Rf_error("Non-finite latent proposal");
        double value = sign * fmax(lower, proposed);
        memcpy(work->trial_response[s], g->response, g->n * sizeof(double));
        work->trial_response[s][i] = value;
        imr_coefficient_density candidate = statistics(work, s, g->selected, work->trial_response[s], work->second_mean[s]);
        double proposed_score = log_kernel(candidate, g->variance);
        double current_score = log_kernel(current, g->variance);
        double ratio = proposed_score - current_score + dnorm4(g->response[i], location, scale, 1) - dnorm4(value, location, scale, 1);
        state->acceptance[IMR_LATENT_PROPOSALS]++;
        if (imr_marginal_accept(ratio)) { g->response[i] = value; state->acceptance[IMR_LATENT_ACCEPTS]++; }
    }
}

void imr_coefficients_run(imr_marginal_state *state)
{
    imr_coefficient_workspace work = {0};
    work.state = state;
    work.models = (imr_coefficient_model **)R_alloc(state->groups_count, sizeof(imr_coefficient_model *));
    work.first_mean = (double **)R_alloc(state->groups_count, sizeof(double *));
    work.second_mean = (double **)R_alloc(state->groups_count, sizeof(double *));
    work.trial_response = (double **)R_alloc(state->groups_count, sizeof(double *));
    work.trial_selection = (int **)R_alloc(state->groups_count, sizeof(int *));
    for (int s = 0; s < state->groups_count; ++s) {
        imr_marginal_group *g = state->groups + s;
        work.models[s] = NULL;
        work.first_mean[s] = (double *)R_alloc(g->p, sizeof(double));
        work.second_mean[s] = (double *)R_alloc(g->p, sizeof(double));
        work.trial_response[s] = (double *)R_alloc(g->n, sizeof(double));
        work.trial_selection[s] = (int *)R_alloc(g->p, sizeof(int));
    }
    int total = state->burnin + state->draws * state->thin;
    for (int sweep = 0; sweep < total; ++sweep) {
        int iteration = sweep + 1;
        R_CheckUserInterrupt();
        for (int s = 0; s < state->groups_count; ++s) update_latent(&work, s);
        for (int l = 0; l < state->platforms_count; ++l) {
            imr_marginal_platform *p = state->platforms + l;
            for (int f = 0; f < p->features; ++f) {
                double log_bf[IMR_MARGINAL_MAX_SUBGROUPS];
                for (int i = 0; i < p->members; ++i) {
                    int s = p->groups[i], j = p->columns[i + p->members * f];
                    imr_marginal_group *g = state->groups + s;
                    int *trial = work.trial_selection[s];
                    memcpy(trial, g->selected, g->p * sizeof(int)); trial[j] = 1;
                    imr_coefficient_density included = statistics(&work, s, trial, g->response, work.first_mean[s]);
                    double included_score = log_kernel(included, g->variance);
                    trial[j] = 0;
                    imr_coefficient_density excluded = statistics(&work, s, trial, g->response, work.second_mean[s]);
                    log_bf[i] = included_score - log_kernel(excluded, g->variance);
                }
                int pattern = imr_marginal_select(state, p, log_bf);
                for (int i = 0; i < p->members; ++i)
                    state->groups[p->groups[i]].selected[p->columns[i + p->members * f]] =
                        REAL(p->patterns)[pattern + p->patterns_count * i];
            }
            int attempts = p->features >= 2 ? (int)ceil(state->swap_rate * p->features) : 0;
            for (int attempt = 0; attempt < attempts; ++attempt) {
                int pair[2], moves[IMR_MARGINAL_MAX_SUBGROUPS] = {0}, count = 0;
                imr_marginal_pair(p, pair);
                double ratio = 0;
                for (int i = 0; i < p->members; ++i) {
                    int s = p->groups[i];
                    imr_marginal_group *g = state->groups + s;
                    int first = p->columns[i + p->members * pair[0]], second = p->columns[i + p->members * pair[1]];
                    if (g->selected[first] == g->selected[second]) continue;
                    int *trial = work.trial_selection[s];
                    memcpy(trial, g->selected, g->p * sizeof(int));
                    trial[first] = g->selected[second]; trial[second] = g->selected[first];
                    imr_coefficient_density candidate = statistics(&work, s, trial, g->response, work.first_mean[s]);
                    double candidate_score = log_kernel(candidate, g->variance);
                    imr_coefficient_density current = statistics(&work, s, g->selected, g->response, work.second_mean[s]);
                    ratio = ratio + candidate_score - log_kernel(current, g->variance);
                    moves[i] = 1; count++;
                }
                if (count) {
                    state->acceptance[IMR_SWAP_PROPOSALS]++;
                    if (imr_marginal_accept(ratio)) {
                        for (int i = 0; i < p->members; ++i) if (moves[i]) {
                            int s = p->groups[i];
                            memcpy(state->groups[s].selected, work.trial_selection[s], state->groups[s].p * sizeof(int));
                        }
                        state->acceptance[IMR_SWAP_ACCEPTS]++;
                    }
                }
            }
            imr_marginal_interaction(state, p);
        }
        for (int s = 0; s < state->groups_count; ++s) {
            imr_marginal_group *g = state->groups + s;
            imr_coefficient_density current = statistics(&work, s, g->selected, g->response, work.first_mean[s]);
            double proposal = exp(log(g->variance) + rnorm(0, g->variance_step));
            state->acceptance[IMR_VARIANCE_PROPOSALS]++;
            if (!R_FINITE(proposal) || proposal <= 0) continue;
            double proposed_score = log_kernel(current, proposal), old_score = log_kernel(current, g->variance);
            double ratio = proposed_score - old_score + log(proposal) - log(g->variance);
            if (imr_marginal_accept(ratio)) { g->variance = proposal; state->acceptance[IMR_VARIANCE_ACCEPTS]++; }
        }
        for (int s = 0; s < state->groups_count; ++s) {
            imr_marginal_group *g = state->groups + s;
            imr_coefficient_density current = statistics(&work, s, g->selected, g->response, work.first_mean[s]);
            restore_coefficients(state, g, current);
        }
        if (iteration > state->burnin && (iteration - state->burnin) % state->thin == 0) {
            double score = 0;
            for (int s = 0; s < state->groups_count; ++s) {
                imr_marginal_group *g = state->groups + s;
                imr_coefficient_density current = statistics(&work, s, g->selected, g->response, work.first_mean[s]);
                score += log_kernel(current, g->variance);
            }
            imr_marginal_store(state, (iteration - state->burnin) / state->thin - 1, score);
        }
        imr_marginal_progress(state, iteration);
    }
}
