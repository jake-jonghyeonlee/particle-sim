# arch: uncertainty sampler

`NormalSpec` is converted by the inverse normal CDF; `UniformSpec` is converted linearly. A Sobol point is
uniform before transformation, so a NormalSpec produces a finite quasi-random approximation to the
specified normal distribution.

Scenario schema 2.0 accepts independent x/y overlay distributions per layer. Deterministic dimension order
uses sorted layer keys, then `.x` and `.y`, followed by CD, diameter, center offset, embedding, planarization
and contact margin variables. `build_scenarios` uses scrambled Sobol `random_base2(m)` and stores seed,
seed-local sample index and global scenario ID. Degenerate distributions produce one fixed scenario per
seed request path.

`simulate-probability --seed N --additional-seed M ...` combines all requested seeds while convergence is
reported separately for each seed. The run manifest records every seed and the scenario evidence.
Current correlations are explicitly independent; covariance and empirical distributions are deferred.
