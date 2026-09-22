# NSCLC outcome-free evaluation input schema

Institutional patient-level inputs are restricted and are not distributed.
The public runner expects one calibration file and one outcome-free target file
for each prespecified predictor, site, and seed.

```text
<input-root>/pet_radiomic_candidate/seed42/calibration_Stanford.csv
<input-root>/pet_radiomic_candidate/seed42/target_Stanford.csv
```

The corresponding folders for the SRDO sensitivity analysis use
`pet_radiomic_candidate_srdo`. Sites are `Stanford` and `VA`; seeds are 42, 52,
62, and 72.

Calibration files contain exactly:

```text
patient_id,time,event,q0_surv2y,q1_surv2y
```

Outcome-free target files contain exactly:

```text
patient_id,q0_surv2y,q1_surv2y
```

The target loader rejects any outcome-like column. Prediction probabilities
must lie strictly between zero and one. The synthetic fixture follows the same
schema but contains no real patient information and is intended only to test
the computational path.
