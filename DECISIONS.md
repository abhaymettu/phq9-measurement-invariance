# Judgment calls

## WLSMV primary, MLR reported alongside

The items are four-category and skewed, sparsest cell 0.45%, so `ordered=` with WLSMV is the
defensible estimator and MLR treating them as continuous is a shortcut. I reported both
rather than picking one, because the whole point of the project is that a modelling choice
nobody reports can change the published conclusion, and here it does: for race/ethnicity
WLSMV says invariance holds and MLR says it fails. Rejected alternative: report only WLSMV
and mention MLR in a limitation. That would have hidden the most interesting result.

**What would change my mind:** if the two estimators had agreed everywhere, reporting both
would be padding and I would have kept only WLSMV.

## Latent response intercepts fixed at zero for the ordinal models

Not a style choice, a correctness fix. lavaan's multigroup default for ordered indicators
frees the item intercepts in every group but the first the moment thresholds are constrained,
and a free intercept absorbs a threshold shift on the same item exactly, so threshold
invariance becomes untestable. Location has to be carried by the thresholds alone. The
simulated-non-invariance control is what exposed this: a planted +0.6 shift was producing a
chi-square rise of 9 on 9 df, which is nothing. The fix moved the WLSMV verdict for sex from
"passes" to "fails", so the first version of this repo would have published the wrong answer.

**What would change my mind:** nothing about the direction, but `semTools::measEq.syntax`
implements the full Wu and Estabrook sequence properly and if it were installed I would use
it instead of hand-specifying, and would add residual-variance invariance as a fourth rung.

## Delta CFI and delta RMSEA as the decision rule, not chi-square

At n = 5,068 the scaled chi-square difference test rejects every rung, including on random
splits where nothing is different. I print it once to make the point concrete (144.4 on 8 df,
p < 2.2e-16 for sex) and then do not use it. Conventional cutoffs of 0.010 for CFI and 0.015
for RMSEA are the rule.

**What would change my mind:** these cutoffs were derived under continuous-indicator ML
simulations and are known to behave differently for WLSMV with skewed ordinal items. If I
had to defend the exact 0.0130 vs 0.010 margin for sex, the honest answer is a specification
search bootstrap, not a citation. The random-split control is my substitute: it shows what
delta CFI looks like in this data when there is genuinely nothing to find, and -0.0130 is
outside that range while the age result is inside it.

## Delta CFI is the wrong criterion for the simulated control, and I said so

Data generated from a true one-factor model fits at CFI = 1.000 at every rung, so delta CFI
cannot move and reporting "pass" there would be meaningless. For the planted-non-invariance
check the criterion is the chi-square rise and the score-test ranking. Using different
criteria in different places needs to be stated out loud, so it is stated in the script
output, the README, and here.

## Two age bands, six race categories

18-49 vs 50+ rather than the five bands used in the earlier `nhanes-depression` work.
Category 3 of the suicidal-ideation item is 0.45% of the sample; five age bands means cells
in the single digits and a polychoric correlation that is not estimable. Coarse bands are
the price of testing an ordinal model on a rare-endorsement item. Race/ethnicity is kept at
the six `RIDRETH3` categories because collapsing them would defeat the purpose.

**What would change my mind:** if the goal were the age gradient rather than the measurement
question, I would collapse the response categories to 0 vs 1+ instead of collapsing the age
bands.

## Unweighted CFA, weighted descriptives

Measurement invariance is a property of the instrument's parameters. Weighting changes the
target of estimation to a population covariance matrix, which is a different and defensible
question, but `lavaan.survey` is not installed and hand-rolling replicate-weight SEs for 130
model fits was not worth the runtime. The population quantity the bias actually contaminates,
the female-male gap, is reported with the correct MEC weights, strata and PSUs. The README
does not claim representativeness for the CFA.

**What would change my mind:** if the invariance verdict were marginal, the weighted fit
could plausibly flip it and the shortcut would not be defensible. Sex fails by 0.003 to
0.006 past the cutoff under both estimators, which is a margin thin enough that I flagged
it as a limitation rather than dismissing it.

## Decomposition done in the MLR parameterisation

The consequence question is "how many sum-score points of the observed gap are artifact",
and the sum score is linear in the parameters only in the continuous parameterisation. In
the ordinal model the expected sum is a nonlinear function of thresholds and would need
numerical integration for the same answer. I chose the exact linear decomposition in the
labelled-as-secondary parameterisation over an approximate one in the primary. The item
released (item 4) is in both estimators' offender set, so the two are not telling different
stories.

## Complete-case, not imputation

465 of 5,533 adults dropped. 440 of them answered zero items, which is unit nonresponse on
the module and not something imputation should invent nine correlated item responses for.
Only 25 answered some but not all. Imputing 25 people cannot move an invariance test on
5,068. Reported at each step instead of imputed.

## Reference-group robustness

With `std.lv` the reference group's factor variance is fixed at 1, so the split between the
latent and measurement parts of the gap is not guaranteed to be reference invariant. I
refit with the rows reversed, which flips which group lavaan calls group 1, and the split is
unchanged to three decimals on the real data. In the test suite the same check on simulated
data moves by about 0.02 points, which is why the test tolerance is 0.05 and not 1e-6. A
test that asserted exact invariance there would have been asserting something false.

## Group order is resolved from the fit, never from factor levels

lavaan numbers groups by order of first appearance in the data frame, not by factor level
order. In this data the first respondent is female, so `group == 1` is Female even though
the factor levels are Male then Female. Reading group 1 as the first level silently flipped
the sign of the entire decomposition, which is how the bug was found: the model said females
were 0.19 SD lower on latent depression while the observed sum score said they were 0.93
points higher. `decompose_gap` resolves identity through `lavInspect(fit, "group.label")`
and a test shuffles the rows to prove the answer does not move.

## The IRT extension: why a second method at all

The CFA already answers "is this instrument invariant across sex". The IRT stage answers a
narrower question the CFA cannot: holding the trait fixed, how many raw points of a 0 to 3
item do the two groups differ by. That effect size, not another p-value, is what the CFA side
was missing, and at n = 5,068 it is the only part of the output that carries information. The
two methods are kept side by side in 07 rather than one replacing the other.

**What would change my mind:** if the two rankings had disagreed sharply, the interesting
report would be the disagreement and I would have led with it. They agree on item 4, so the
comparison is a robustness check rather than the finding.

## Anchor choice is reported as a choice, not defended as correct

DIF is defined relative to the items assumed clean, so every DIF result is conditional on an
anchor and most papers never say which. The all-other-item default flags 6 of 9 items here.
The purified rank-based anchor of 4 items (Woods 2009) flags 5. Sweeping the anchor size from
1 to 6 moves the flagged set from six items down to three, and only items 4, 5 and 9 survive
every anchor size. Those three are what 05 reports as stable; the rest are reported as
anchor-dependent, which is what they are.

**What would change my mind:** an anchor set justified by content rather than by statistics.
If there were a substantive reason to believe specific items are clean across sex, that would
beat any purification rule, and the sweep would become a sensitivity analysis rather than the
result.

## The contamination claim had to be tested, not asserted

The first version of 05 explained 6-of-9 flagging as anchor contamination and left it there.
The null runs in 06 then showed both procedures are *conservative* when there is no DIF: 0.004
and 0.013 false positives per item test against a nominal 0.05, because Benjamini-Hochberg
across nine tests is strict. So the random-split control cannot support the contamination
story at all. It contains no real DIF, and there is nothing to contaminate an anchor with.

Section 4 of 06 exists because of that gap. Planting DIF on a growing number of items and
watching the false-alarm rate on the innocent ones is the condition that actually tests the
claim, and it is the only place the two anchor strategies separate. This is the same failure
mode the covary work cost eight rounds to learn: a check that passes for a reason unrelated to
the property it is named after.

## Detection floor reported in item-score points, not logits

"The procedure has 80% power at 0.19 logits" is unreadable. Every threshold is converted to
the same expected-score effect size 05 reports in, using a large simulated sample so the
conversion reflects the model rather than sampling noise. That conversion is what lets the
central comparison be made at all: the largest real sex effect in this data sits just above
the floor, which is a much weaker statement than "five items show DIF" and is the honest one.

## Unweighted here too, for the same reason as the CFA

Same argument as the unweighted CFA above, and the same limitation. Item parameters are
properties of the instrument; the population quantity that the bias contaminates is reported
with MEC weights in 02.

## Every DIF flag is reported with the anchor it depends on

The single most repeated failure in this literature is a flagged-item list with no anchor
named. Item 7 flags here and ranks last of nine on the CFA, with an effect of +0.010
item-score points. It is in the output because n is large, and it is labelled that way rather
than dropped, because dropping it would be a specification search and reporting it plainly
would be a false finding. Only items 4, 5 and 9 survive every anchor size and only those are
described as the result.
