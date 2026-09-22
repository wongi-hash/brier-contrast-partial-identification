# Rotterdam–GBSG supporting survival analysis

This directory reconstructs the source predictors and runs the outcome-free
analysis using the public `rotterdam` and `gbsg` datasets in the R package
`survival`.

```bash
python run_public_pipeline.py --bootstrap 999
```

The default pipeline exports the public data, prepares source-trained locked
predictors, runs the outcome-free analysis, and compares direct identification
with separate-risk subtraction. Add `--open-outcomes` only for the explicitly
retrospective outcome comparison.

GBSG outcomes are stored separately and are rejected by the outcome-free target
loader. The analysis is a supporting external demonstration, not prospectively
blinded clinical validation; target outcomes had been accessed in antecedent
exploratory work. Predictor definitions and decision rules used here are fully
documented in `supporting_analysis_specification.md`.

The data export was verified with R package `survival` 3.8-3. Raw exported CSV
files are generated locally and are not committed.
