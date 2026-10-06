# Reading the simulation diagnostics

Each saved dataset is used for SDR and logistic TMLE at both outcome times. Failures and pending replications are reported separately. Coverage among successful fits and coverage counting a failed fit as producing no interval should both be reported; unfinished datasets cannot be counted as completed.

The columns `mean_sequential_eif` and `mean_bridge_adjoint_eif` are sample averages of the two contributions, with the true parameter subtracted in the sequential contribution. Their sum is `estimate - truth`. Individual influence values centered at the estimated parameter instead have total sample mean zero.

`population_bias` evaluates the already fitted contribution on every point of the generating distribution. `sequential_remainder` compares the sequential contribution with the parameter formed using the fitted bridge; `bridge_remainder` accounts for the bridge error and its correction. They add to `population_bias`. These are different from the sample mean EIF contributions. TMLE targeting uses validation observations, so the fixed-function population evaluation is not a literal conditional expectation given only the training observations. No true function is supplied to a fitted learner.

`bridge_equation_max_error` and `adjoint_equation_max_error` are the largest absolute conditional-equation errors over reachable predictor combinations. A large maximum in a rare combination need not have a large contribution to estimation error. The full conditional equations and weighted remainder must be examined together.

`lambda_mse_to_selected_solution` compares with one selected adjoint solution. It is not by itself an error criterion when the solution is nonunique. The defining equations are the relevant check in those cases.

`beta_upper_fraction` is a historical diagnostic: the probability-weighted fraction of measured observations with fitted bridge at least 6 (up to numerical tolerance). It does **not** count clipping at the current application limit 100. `beta_lower_fraction` counts fitted values within numerical tolerance of 1. Raw cmbridge fits are not clipped after fitting.
