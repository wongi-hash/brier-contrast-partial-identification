# Rotterdam–GBSG supporting survival analysis specification

This supporting analysis evaluates whether the outcome-free procedure can
certify a sufficiently separated, literature-defined predictor contrast in a
public censored-survival dataset. It is not presented as prospectively blinded
clinical validation: GBSG outcomes had been accessed in antecedent exploratory
work. The predictor definitions below were fixed before the finalized
outcome-free computation and all prespecified decisions are reported.

## Data and endpoint

- Source: Rotterdam node-positive subset, n=1,546.
- Target: GBSG, n=686.
- Endpoint: recurrence-free survival, defined by recurrence or death.
- Prediction horizon: 5 years (1,826.25 days).
- Administrative censoring: 7 years.
- Source split: 2:1 training/calibration, event-stratified.
- Seeds: 42, 52, 62, and 72.

## Predictor contrast

Both Cox predictors use transformations and standardization fitted in the
source training sample and a training-sample Breslow baseline.

The reference predictor contains age fractional-polynomial terms, menopausal
status, tumor-size categories, and hormonal treatment. The candidate adds
positive nodes transformed as nodes^(-1/2) and estrogen receptor. These added
variables are part of the established Royston–Altman model and were not chosen
by inspecting the finalized contrast.

## Identification and inference

- Deployment variables: Z=(q0,q1).
- Source outcome regression: IPCW logistic regression in the source calibration
  sample.
- Conditional censoring model: Cox proportional hazards given Z, with Breslow
  baseline and nuisance refitting in every bootstrap draw.
- Sensitivity class: absolute source-to-target logit shift bounded by gamma.
- Reported sensitivity parameters: 0, 0.1, 0.2, 0.3, and 0.5.
- Simultaneous grid: 0 to 1 by 0.01.
- Joint pairs bootstrap: 999 draws; predictors remain frozen.
- Decision rule: ADOPT CANDIDATE if the upper simultaneous limit is below zero;
  KEEP REFERENCE if the lower simultaneous limit is above zero; otherwise
  DEFER.

The public pipeline preserves the separation between outcome-free computation
and retrospective outcome comparison. Dataset, predictor definitions,
sensitivity values, horizon, and seeds are not modified after inspecting the
reported decisions.
