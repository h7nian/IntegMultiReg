/* Independent enumeration of proposals and the symmetric MRF target. */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "sampler_math.h"

static int count(int state, int n) {
    int result = 0;
    for (int j = 0; j < n; ++j) result += (state >> j) & 1;
    return result;
}

static double proposal(int from, int to, int n) {
    int selected = count(from, n);
    int changed = count(from ^ to, n);
    if (changed == 1) return (selected == 0 || selected == n ? 1.0 : .5) / n;
    if (changed == 2 && count(to, n) == selected)
        return .5 / (selected * (n - selected));
    return 0;
}

int main(void) {
    for (int n = 1; n <= 3; ++n) {
        int size = 1 << n;
        double target[8], transition[8][8] = {{0}}, normalizer = 0;
        double theta[2] = {0, .4};
        _Bool current[3], neighbor[3] = {1, 0, 1};
        _Bool *selection[2] = {current, neighbor};
        for (int state = 0; state < size; ++state) {
            double log_target = .13 * state - .03 * state * state;
            for (int j = 0; j < n; ++j)
                log_target += ((state >> j) & 1) * (-.7 + 2 * theta[1] * neighbor[j]);
            target[state] = exp(log_target);
            normalizer += target[state];
        }
        for (int state = 0; state < size; ++state) target[state] /= normalizer;
        for (int from = 0; from < size; ++from) {
            double leaving = 0;
            for (int j = 0; j < n; ++j) current[j] = (from >> j) & 1;
            for (int to = 0; to < size; ++to) {
                double q = proposal(from, to, n);
                if (!q) continue;
                double ratio = .13 * (to - from) - .03 * (to * to - from * from);
                for (int j = 0; j < n; ++j) {
                    double odds = imr_gamma_log_odds(2, 0, j, theta, selection, -.7, IMR_SAMPLER_PAPER);
                    assert(fabs(odds - (-.7 + .8 * neighbor[j])) < 1e-14);
                    assert(fabs(imr_gamma_log_odds(2, 0, j, theta, selection, -.7, IMR_SAMPLER_LEGACY) -
                                (-.7 + .4 * neighbor[j])) < 1e-14);
                    ratio += (((to >> j) & 1) - current[j]) * odds;
                }
                ratio += imr_gamma_log_hastings(n, count(from, n), count(to, n), .5);
                transition[from][to] = q * fmin(1, exp(ratio));
                leaving += transition[from][to];
            }
            assert(leaving <= 1 + 1e-14);
            transition[from][from] = 1 - leaving;
        }
        for (int to = 0; to < size; ++to) {
            double stationary = 0;
            for (int from = 0; from < size; ++from) {
                assert(fabs(target[from] * transition[from][to] -
                            target[to] * transition[to][from]) < 1e-13);
                stationary += target[from] * transition[from][to];
            }
            assert(fabs(stationary - target[to]) < 1e-13);
        }
    }
    puts("MRF increments and detailed balance passed for 1, 2 and 3 features.");
    return 0;
}
