using Test
using Random
using DataFrames
using Distributions

include(joinpath(@__DIR__, "..", "src", "BayesREProbit.jl"))
using .BayesREProbit

rng = MersenneTwister(7301)
n = 120
J = 8
id = repeat(1:J, inner=n ÷ J)
x2 = randn(rng, n)
x3 = rand(rng, Bernoulli(0.45), n)
x4 = randn(rng, n)
a = rand(rng, Normal(0, 0.7), J)
eta = 0.2 .- 0.9 .* x2 .+ 0.5 .* x3 .+ a[id]
y = Int.(randn(rng, n) .< eta)

dat = DataFrame(y=y, id=id, x2=x2, x3=Float64.(x3), x4=x4)
allowmissing!(dat, [:x2, :x3, :x4])
dat.x2[randperm(rng, n)[1:24]] .= missing
dat.x3[randperm(rng, n)[1:18]] .= missing
dat.x4[randperm(rng, n)[1:20]] .= missing

reference_fit = BayesREProbit.run_medoid_reference_imputation(
    dat, ["x2", "x3"], ["x2", "x3", "x4"];
    rep=1, mechanism="MCAR",
    n_iter=160, burnin=40, thin=4, seed=91,
    impute_mode="cartbb_every", cart_every=5,
    cart_max_depth=3, cart_min_samples_leaf=5,
    verbose_every=0
)

medoid = BayesREProbit.select_medoid_imputation(reference_fit)
completed = BayesREProbit.completed_dataframe_from_design(
    dat, ["x2", "x3", "x4"], medoid.X
)

@test !any(ismissing, completed.x2)
@test !any(ismissing, completed.x3)
@test !any(ismissing, completed.x4)
@test reference_fit["analysis_terms"] == ["(Intercept)", "x2", "x3"]
@test 1 <= medoid.draw_id <= size(reference_fit["imputation_draws"], 1)

fit = BayesREProbit.run_gibbs_re_probit(
    completed, ["x2", "x3"];
    model="M_test", rep=1, mechanism="MCAR", data_version="medoid",
    n_iter=220, burnin=60, thin=4, seed=92, impute_mode="none",
    chib_reduced_iter=120, chib_reduced_burnin=30,
    chib_reduced_thin=2, medoid_draw_id=medoid.draw_id,
    verbose_every=0
)

evidence = BayesREProbit.marginal_likelihood_dataframe(fit)
@test nrow(evidence) == 3
@test all(isfinite, evidence.log_marginal_likelihood)
@test all(evidence.evidence_target .== "conditional_on_fixed_X")

@test_throws ErrorException BayesREProbit.chib_marginal_likelihood(
    reference_fit; fixed_point="mean"
)

println(evidence[:, [:fixed_point, :log_marginal_likelihood]])
