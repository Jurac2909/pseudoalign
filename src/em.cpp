#include "em.hpp"
#include <cmath>
#include <numeric>

// Expected number of positions where a randomly-sampled fragment can start
// inside a transcript of length L, under fragment length distribution
// N(mean_fl, sd_fl) truncated to [1, L].
//
// eff_len = Σ_{l=1}^{L} f(l)·(L-l+1) / Σ_{l=1}^{L} f(l)
//
// Normalising only over [1,L] (not [1,∞)) gives the conditional expectation
// E[L-l+1 | fragment fits in transcript], matching kallisto's formula.
// For L >> mean_fl this reduces to L - mean_fl + 1.
// For L < mean_fl the denominator is the small tail mass below L, which
// makes eff_len small (e.g. ~7 for a 71 bp transcript with mean_fl=200),
// strongly down-weighting transcripts that are shorter than a typical fragment.
static double compute_eff_len(int L, int /*k*/, double mean_fl, double sd_fl)
{
    double num = 0.0, denom = 0.0;
    // For long transcripts the Gaussian weight is negligible outside
    // [mean-4σ, mean+4σ], so cap the upper bound to avoid iterating
    // millions of zero-weight steps.  For short transcripts (L < mean)
    // we must go all the way to L to get the correct conditional mass.
    int lmax = std::min(L, static_cast<int>(mean_fl + 4.0 * sd_fl) + 1);
    for (int l = 1; l <= lmax; ++l) {
        double w = std::exp(-0.5 * std::pow((l - mean_fl) / sd_fl, 2));
        num   += w * (L - l + 1);
        denom += w;
    }
    return (denom < 1e-15) ? 1.0 : std::max(1.0, num / denom);
}

EMResult run_em(const ECCounts& ec_counts, const DBG& g,
                double mean_fl, double sd_fl,
                int max_iter, double tol)
{
    const int T = static_cast<int>(g.transcript_names.size());

    std::vector<double> eff_len(T);
    for (int t = 0; t < T; ++t)
        eff_len[t] = compute_eff_len(static_cast<int>(g.transcript_lengths[t]),
                                     g.k, mean_fl, sd_fl);

    uint64_t total = 0;
    for (const auto& [tids, cnt] : ec_counts) total += cnt;

    // Initialise alpha proportionally to effective read share:
    // for each EC, split its read count equally among member transcripts
    // (then weight by 1/eff_len so short transcripts don't dominate).
    // This gives every transcript that appears in any EC a non-zero start,
    // while still suppressing transcripts whose k-mers only appear in
    // highly-specific ECs where other transcripts have much larger weight.
    std::vector<double> alpha(T, 0.0);
    for (const auto& [tids, cnt] : ec_counts) {
        double share = static_cast<double>(cnt) / tids.size();
        for (uint32_t t : tids) alpha[t] += share;
    }
    double alpha_sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
    if (alpha_sum == 0.0)
        std::fill(alpha.begin(), alpha.end(), 1.0 / T);
    else
        for (double& a : alpha) a /= alpha_sum;

    std::vector<double> alpha_new(T);

    for (int iter = 0; iter < max_iter; ++iter) {
        std::fill(alpha_new.begin(), alpha_new.end(), 0.0);

        for (const auto& [tids, cnt] : ec_counts) {
            double denom = 0.0;
            for (uint32_t t : tids) denom += alpha[t] / eff_len[t];
            if (denom == 0.0) continue;
            for (uint32_t t : tids)
                alpha_new[t] += static_cast<double>(cnt) * (alpha[t] / eff_len[t]) / denom;
        }

        double sum = std::accumulate(alpha_new.begin(), alpha_new.end(), 0.0);
        if (sum == 0.0) break;
        for (double& a : alpha_new) a /= sum;

        double delta = 0.0;
        for (int t = 0; t < T; ++t) delta += std::fabs(alpha_new[t] - alpha[t]);
        alpha = alpha_new;
        if (delta < tol) break;
    }

    std::vector<double> counts(T);
    for (int t = 0; t < T; ++t) counts[t] = alpha[t] * static_cast<double>(total);

    std::vector<double> rpk(T);
    for (int t = 0; t < T; ++t) rpk[t] = alpha[t] / eff_len[t];
    double scale = std::accumulate(rpk.begin(), rpk.end(), 0.0) / 1e6;

    std::vector<double> tpm(T, 0.0);
    if (scale > 0.0)
        for (int t = 0; t < T; ++t) tpm[t] = rpk[t] / scale;

    return {counts, tpm};
}
