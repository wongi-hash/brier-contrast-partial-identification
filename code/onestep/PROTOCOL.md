# One-step endpoint estimator — prespecified protocol

The protocol was fixed before any simulation or application run of the one-step
estimator. The one-step analysis is secondary. It was added after the primary
plug-in analyses had been finalized, and it does not change any primary grid,
seed, or decision.

## Estimand (unchanged)

For an endpoint E ∈ {L, U} and a sensitivity value γ:

- E_γ = E_T[f_E(Z, m_S(Z))], where f_E(z, m) = a(z) − 2 b(z) expit{logit m + s_E(z) γ}.
- s_L(z) = +1 if b(z) ≥ 0, and −1 otherwise. s_U(z) = −s_L(z).
- ∂_m f_E(z, m) = −2 b(z) g′(logit m + s_E γ) / {m(1 − m)}, with g′ = expit (1 − expit).

## Estimator

Ê_os = (1/n_T) Σ_T f_E(Z_j, m̂(Z_j)) + (1/n_S) Σ_S r̂(Z_i) ∂_m f_E(Z_i, m̂(Z_i)) ω̂_i {Y_i − m̂(Z_i)}

The plug-in estimator (the primary analysis) is the first term only.

**Outcome and censoring nuisances (m̂, ω̂).** These are identical to the primary implementation:

- conditional Cox–Breslow censoring model on H = Z, with censoring probabilities floored at 0.02;
- IPCW-weighted logistic regression on an intercept and standardized (q₀, q₁);
- m̂ clipped to [10⁻⁴, 1 − 10⁻⁴].

**Density ratio (r̂).** r̂ is estimated by an unpenalized logistic domain classifier
(target versus source) on an intercept and standardized (logit q₀, logit q₁), pooled
over the source-calibration and target Z. It is fitted by Newton's method with a
pseudo-inverse. Then:

1. r̂(z) = exp{η̂(z)} · (n_S / n_T);
2. r̂ is normalized to mean 1 over the source sample;
3. r̂ is capped at 20, and each activation of the cap is recorded.

## Inference

- **Bootstrap:** joint pairs bootstrap, resampling source and target independently and refitting G, m_S, and the classifier in every draw.
- **Envelopes:** method-specific studentized sup-t critical value and outer envelope [Ê_L − ĉ ŝ_L, Ê_U + ĉ ŝ_U]. The decision rule is unchanged.
- **Simulations:** the critical value is taken over 𝒢 = {0, 0.1, 0.2, 0.3, 0.5}.
- **Applications:** the bootstrap draws and the critical-value grid {0, 0.01, …, 1.00} are those of the reported plug-in analysis.

## Simulation design (ADEMP)

**Data-generating mechanism.** The primary-implementation mechanism:

- X_S ~ N(0, 1) and X_T ~ N(0.5, 1);
- q₀ = expit(−0.15 + 0.45x) and q₁ = expit(0.10 + 0.80x);
- Weibull event times with shape 1.5 and τ = 1;
- exponential censoring with multiplier exp{z_std · (1, −0.75)}, with the baseline calibrated to the target censoring rate.

**Source laws.**

| Law | Source logit | Role |
|---|---|---|
| A | −0.35 + 1.25 q₀ − 0.70 q₁ | correctly specified |
| B | law A + 2.5 (q₁ − 0.55)² | misspecified working model |

The omitted term in law B is the form used in the outcome-model misspecification
scenario of the application-scale comparator; it was not tuned.

**Scales.**

| Scale | n_S | n_T | Censoring |
|---|---|---|---|
| NSCLC | 114 | 80 | 65% |
| Rotterdam–GBSG | 516 | 686 | 38% |

**Replication.**

- Cells: 2 laws × 2 scales = 4.
- Each cell: 1,000 replications, B = 199.
- Replication seed: `stable_seed("onestep-sim-20261028", law, scale_key, rep)`.
- Truth: endpoints under the true m_S on X_T ~ N(0.5, 1), computed from 2,000,000 Monte Carlo draws (seed 990011).

**Metrics.** Both estimators are computed on the same data and the same bootstrap draws:

- bias, empirical SD, and RMSE of L̂_γ and Û_γ at each γ;
- simultaneous outer coverage with a Wilson 95% CI;
- envelope width;
- fatal fits, density-ratio cap activation, and censoring-floor activation.

## Applications (secondary)

- **NSCLC:** primary PET-radiomic candidate, seeds 42/52/62/72, Stanford and VA, τ = 730 days. The institutional inputs are restricted.
- **Rotterdam–GBSG:** supporting clinical candidate, seeds 42/52/62/72, τ = 1,826.25 days. The inputs are produced by the public pipeline.
- The one-step identified sets, envelopes, and decisions are reported alongside the plug-in values computed from the same draws.
- Target outcomes are not read.
