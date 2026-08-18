# Resume bullets

- Tested measurement invariance of the nine PHQ-9 items across sex, age and race/ethnicity in
  NHANES 2017-2018 using ordinal CFA in lavaan, n = 5,068, finding the scale fails scalar
  invariance across sex and that 11% of the observed female-male gap in sum scores is item
  non-invariance rather than symptom burden.
- Validated the pipeline against 30 random group splits (0 of 60 false positives) and
  simulated data with a planted item shift, which caught a lavaan identification default that
  had made threshold invariance untestable and had reversed the primary result.
- Extended the analysis to item response theory, fitting a graded response model in `mirt` and
  testing differential item functioning by sex against a purified anchor, then measured the
  DIF procedure's own error rate rather than reporting the flags alone: 30 random splits put
  the false positive rate at 0.4 to 2.2% per item test against a nominal 5%, and a swept
  magnitude study put 80% detection at 0.054 item-score points, which the largest real effect
  in the data (0.064) only just clears.
- Falsified my own first explanation of the result. The default all-other-item anchor flagged
  6 of 9 items and I attributed it to anchor contamination; the null controls could not
  support that, because a random split contains no real DIF to contaminate an anchor with.
  Built the condition that does test it, planting bias on a growing number of items, and
  showed the default procedure's false positive rate on unbiased items rises from 0.006 to
  0.870 while a purified anchor holds at 0.267.
