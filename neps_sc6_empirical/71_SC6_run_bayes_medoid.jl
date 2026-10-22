#!/usr/bin/env julia

# 71_SC6_run_bayes_medoid.jl
#
# NEPS SC6 Wave-16 Bayesian CART-BB workflow:
#   1. read numeric R export
#   2. select one common CART-BB medoid completion using a reference chain
#   3. expand unemployment categories into treatment-coded dummies after medoid selection
#   4. fit M1--M5 random-effects probit models on the same fixed medoid
#   5. compute conditional Chib evidence p(y | X_med, M_k)
#
# This script intentionally keeps unemployment duration as one categorical code
# during imputation and creates dummies only after the common medoid is fixed.

using CSV
using DataFrames
using Random
using Statistics
using Printf

# CLI parsing independent of the module, so --module can be used before include.
function parse_cli_args_local(args::Vector{String})
    out = Dict{String,String}()
    i = 1
    while i <= length(args)
        key = replace(args[i], "--" => "")
        if i == length(args) || startswith(args[i + 1], "--")
            out[key] = "true"
            i += 1
        else
            out[key] = args[i + 1]
            i += 2
        end
    end
    return out
end

get_arg(args::Dict{String,String}, key::String, default::String) = get(args, key, default)

args = parse_cli_args_local(ARGS)

module_path = get_arg(args, "module", joinpath(@__DIR__, "..", "julia", "src", "BayesREProbit.jl"))
isfile(module_path) || error("BayesREProbit module not found: $module_path")
include(module_path)
using .BayesREProbit

function parse_terms_arg(s::AbstractString)
    isempty(strip(s)) && return String[]
    return String.(strip.(split(s, ",")))
end

function coerce_numeric_missing!(dat::DataFrame)
    for nm in names(dat)
        dat[!, nm] = Union{Missing,Float64}[
            ismissing(x) ? missing : Float64(x) for x in dat[!, nm]
        ]
    end
    return dat
end

function rounded_code(x; lo::Int=1, hi::Int=5)
    if ismissing(x)
        return missing
    end
    v = Int(round(Float64(x)))
    return min(max(v, lo), hi)
end

function complete_from_design_rounded(dat::DataFrame,
                                      terms::Vector{String},
                                      X_completed::Matrix{Float64})
    out = copy(dat)
    for (j, term) in enumerate(terms)
        col = BayesREProbit.find_column_name(out, term)
        miss_rows = ismissing.(out[!, col])
        if any(miss_rows)
            vals = X_completed[miss_rows, j + 1]
            if term == "unemp_code"
                out[miss_rows, col] = Float64.(round.(clamp.(vals, 1.0, 5.0)))
            else
                out[miss_rows, col] = vals
            end
        end
    end
    return out
end

function add_unemployment_dummies!(dat::DataFrame)
    code = [rounded_code(x) for x in dat.unemp_code]
    if any(ismissing.(code))
        error("unemp_code still contains missing values after medoid completion.")
    end
    dat.unemp_lt1  = Float64.(code .== 2)
    dat.unemp_1to2 = Float64.(code .== 3)
    dat.unemp_2to5 = Float64.(code .== 4)
    dat.unemp_ge5  = Float64.(code .== 5)
    return dat
end

function chib_row_dataframe(fit::Dict, fixed_point::String)
    cm = BayesREProbit.chib_marginal_likelihood(fit; fixed_point=fixed_point)
    return DataFrame(
        method = [BayesREProbit.method_label(fit)],
        model = [fit["model"]],
        rep = [fit["rep"]],
        mechanism = [fit["mechanism"]],
        fixed_point = [cm.fixed_point],
        log_likelihood_at_theta = [cm.log_likelihood_at_theta],
        log_prior_at_theta = [cm.log_prior_at_theta],
        log_posterior_ordinate = [cm.log_posterior_ordinate],
        log_ordinate_beta = [cm.log_ordinate_beta],
        log_ordinate_a_given_beta = [cm.log_ordinate_a_given_beta],
        log_ordinate_sigma2_given_beta_a = [cm.log_ordinate_sigma2_given_beta_a],
        reduced_draws = [cm.reduced_draws],
        log_marginal_likelihood = [cm.log_marginal_likelihood],
        beta_intercept = [cm.beta_intercept],
        sigma2_a = [cm.sigma2_a],
        impute_mode = [fit["impute_mode"]],
        data_version = [fit["data_version"]],
        evidence_target = [fit["evidence_target"]],
        medoid_draw_id = [fit["medoid_draw_id"]],
        tau2_beta = [fit["tau2_beta"]],
        nobs = [fit["nobs"]],
        n_groups = [fit["n_groups"]],
        n_iter = [fit["n_iter"]],
        burnin = [fit["burnin"]],
        thin = [fit["thin"]]
    )
end

function pairwise_logbf(logml::DataFrame)
    models = String.(logml.model)
    vals = Float64.(logml.log_marginal_likelihood)
    mat = Matrix{Float64}(undef, length(models), length(models))
    for i in eachindex(models), j in eachindex(models)
        mat[i, j] = vals[i] - vals[j]
    end
    wide = DataFrame(model_i=models)
    for (j, m) in enumerate(models)
        wide[!, m] = mat[:, j]
    end

    long = DataFrame(
        model_i = String[],
        model_j = String[],
        logBF_i_vs_j = Float64[]
    )
    for i in eachindex(models), j in eachindex(models)
        push!(long, (models[i], models[j], mat[i, j]))
    end
    return wide, long
end

data_file = get_arg(args, "data", "")
isempty(data_file) && error("Missing --data path.")
isfile(data_file) || error("Input data not found: $data_file")

outdir = get_arg(args, "outdir", joinpath(pwd(), "results", "sc6_bayes"))
mkpath(outdir)

medoid_reference_model = get_arg(args, "medoid_reference_model", "M5")
fixed_points = parse_terms_arg(get_arg(args, "fixed_points", "mean"))

medoid_n_iter = parse(Int, get_arg(args, "medoid_n_iter", "20000"))
medoid_burnin = parse(Int, get_arg(args, "medoid_burnin", "5000"))
medoid_thin = parse(Int, get_arg(args, "medoid_thin", "10"))
n_iter = parse(Int, get_arg(args, "n_iter", "20000"))
burnin = parse(Int, get_arg(args, "burnin", "4000"))
thin = parse(Int, get_arg(args, "thin", "2"))
chib_reduced_iter = parse(Int, get_arg(args, "chib_reduced_iter", "4000"))
chib_reduced_burnin = parse(Int, get_arg(args, "chib_reduced_burnin", "1000"))
chib_reduced_thin = parse(Int, get_arg(args, "chib_reduced_thin", "2"))

tau2_beta = parse(Float64, get_arg(args, "tau2_beta", "1.0"))
c0 = parse(Float64, get_arg(args, "c0", "2.0"))
d0 = parse(Float64, get_arg(args, "d0", "1.0"))
seed = parse(Int, get_arg(args, "seed", "20260628"))

cart_every = parse(Int, get_arg(args, "cart_every", "10"))
cart_max_depth = parse(Int, get_arg(args, "cart_max_depth", "5"))
cart_min_samples_leaf = parse(Int, get_arg(args, "cart_min_samples_leaf", "10"))
cart_max_donors = parse(Int, get_arg(args, "cart_max_donors", "50"))
verbose_every = parse(Int, get_arg(args, "verbose_every", "1000"))

dat = CSV.read(data_file, DataFrame; missingstring="")
coerce_numeric_missing!(dat)

# Convert outcome and group back to numeric values without missing.
dat.y = Float64.(dat.y)
dat.id = Float64.(dat.id)

imputation_terms = [
    "age_z", "age_z2", "female",
    "years_education", "ict_wle", "ict_meta", "math_wle", "reading_wle",
    "unemp_code",
    "migration_background", "father_unemployed", "mother_unemployed"
]

model_terms = Dict(
    "M1" => ["age_z", "age_z2", "female"],
    "M2" => ["age_z", "age_z2", "female",
             "years_education", "ict_wle", "ict_meta", "math_wle", "reading_wle"],
    "M3" => ["age_z", "age_z2", "female",
             "unemp_lt1", "unemp_1to2", "unemp_2to5", "unemp_ge5"],
    "M4" => ["age_z", "age_z2", "female",
             "migration_background", "father_unemployed", "mother_unemployed"],
    "M5" => ["age_z", "age_z2", "female",
             "years_education", "ict_wle", "ict_meta", "math_wle", "reading_wle",
             "unemp_lt1", "unemp_1to2", "unemp_2to5", "unemp_ge5",
             "migration_background", "father_unemployed", "mother_unemployed"]
)

# Reference chain uses the categorical unemployment code as one imputation target.
reference_terms_by_model = Dict(
    "M1" => ["age_z", "age_z2", "female"],
    "M2" => ["age_z", "age_z2", "female",
             "years_education", "ict_wle", "ict_meta", "math_wle", "reading_wle"],
    "M3" => ["age_z", "age_z2", "female", "unemp_code"],
    "M4" => ["age_z", "age_z2", "female",
             "migration_background", "father_unemployed", "mother_unemployed"],
    "M5" => imputation_terms
)

haskey(reference_terms_by_model, medoid_reference_model) ||
    error("Unknown medoid_reference_model=$medoid_reference_model")

println("SC6 Bayesian CART-BB | data=$data_file")
println("SC6 Bayesian CART-BB | outdir=$outdir")
println("SC6 Bayesian CART-BB | reference model=$medoid_reference_model")

reference_seed = seed + BayesREProbit.deterministic_offset("SC6_wave16_medoid_" * medoid_reference_model)

reference_fit = BayesREProbit.run_medoid_reference_imputation(
    dat,
    reference_terms_by_model[medoid_reference_model],
    imputation_terms;
    rep = 1,
    mechanism = "SC6_wave16",
    n_iter = medoid_n_iter,
    burnin = medoid_burnin,
    thin = medoid_thin,
    tau2_beta = tau2_beta,
    c0 = c0,
    d0 = d0,
    seed = reference_seed,
    impute_mode = "cartbb_every",
    cart_every = cart_every,
    cart_max_depth = cart_max_depth,
    cart_min_samples_leaf = cart_min_samples_leaf,
    cart_max_donors = cart_max_donors,
    verbose_every = verbose_every
)

medoid = BayesREProbit.select_medoid_imputation(reference_fit)
completed = complete_from_design_rounded(dat, imputation_terms, medoid.X)
add_unemployment_dummies!(completed)

CSV.write(joinpath(outdir, "sc6_medoid_completed.csv"), completed)

medoid_meta = DataFrame(
    rep = [1],
    mechanism = ["SC6_wave16"],
    reference_model = [medoid_reference_model],
    reference_analysis_terms = [join(reference_terms_by_model[medoid_reference_model], ",")],
    imputation_terms = [join(imputation_terms, ",")],
    medoid_draw_id = [medoid.draw_id],
    total_scaled_distance = [medoid.total_distance],
    n_imputation_draws = [size(reference_fit["imputation_draws"], 1)],
    n_missing_cells = [size(reference_fit["imputation_draws"], 2)],
    medoid_n_iter = [medoid_n_iter],
    medoid_burnin = [medoid_burnin],
    medoid_thin = [medoid_thin],
    cart_every = [cart_every],
    cart_max_depth = [cart_max_depth],
    cart_min_samples_leaf = [cart_min_samples_leaf],
    cart_max_donors = [cart_max_donors],
    seed = [reference_seed]
)
CSV.write(joinpath(outdir, "bayes_medoid_metadata.csv"), medoid_meta)

all_model_results = DataFrame()
all_estimates = DataFrame()
all_marglik = DataFrame()
all_predictive = DataFrame()

for model in ["M1", "M2", "M3", "M4", "M5"]
    println("SC6 Bayesian CART-BB | fitting $model")
    fit_seed = seed + BayesREProbit.deterministic_offset("SC6_wave16_fit_" * model)
    fit = BayesREProbit.run_gibbs_re_probit(
        completed,
        model_terms[model];
        model = model,
        rep = 1,
        mechanism = "SC6_wave16",
        data_version = "medoid",
        yvar = "y",
        idvar = "id",
        n_iter = n_iter,
        burnin = burnin,
        thin = thin,
        tau2_beta = tau2_beta,
        c0 = c0,
        d0 = d0,
        seed = fit_seed,
        impute_mode = "none",
        chib_reduced_iter = chib_reduced_iter,
        chib_reduced_burnin = chib_reduced_burnin,
        chib_reduced_thin = chib_reduced_thin,
        evidence_target = "conditional_on_common_SC6_wave16_medoid",
        medoid_draw_id = medoid.draw_id,
        verbose_every = verbose_every
    )

    marglik_model = DataFrame()
    for fp in fixed_points
        println("SC6 Bayesian CART-BB | Chib evidence $model | fixed point=$fp")
        append!(marglik_model, chib_row_dataframe(fit, fp))
    end
    append!(all_marglik, marglik_model)

    mean_ml = marglik_model[marglik_model.fixed_point .== "mean", :]
    if nrow(mean_ml) == 0
        mean_ml = marglik_model[1:1, :]
    end

    modelres, est = BayesREProbit.summarize_fit(fit; marglik=mean_ml, compute_evidence=true)
    pred = BayesREProbit.predictive_criteria_dataframe(fit)

    append!(all_model_results, modelres)
    append!(all_estimates, est)
    append!(all_predictive, pred)

    CSV.write(joinpath(outdir, "fit_$(model)_estimates.csv"), est)
    CSV.write(joinpath(outdir, "fit_$(model)_marginal_likelihood.csv"), marglik_model)
end

CSV.write(joinpath(outdir, "bayes_model_results.csv"), all_model_results)
CSV.write(joinpath(outdir, "bayes_estimates.csv"), all_estimates)
CSV.write(joinpath(outdir, "bayes_marginal_likelihood.csv"), all_marglik)
CSV.write(joinpath(outdir, "bayes_predictive_criteria.csv"), all_predictive)

logml_mean = all_marglik[all_marglik.fixed_point .== "mean", :]
if nrow(logml_mean) == 0
    error("No mean fixed-point marginal likelihood rows were produced.")
end
sort!(logml_mean, :model)

bf_wide, bf_long = pairwise_logbf(logml_mean)
CSV.write(joinpath(outdir, "bayes_conditional_logbf_matrix.csv"), bf_wide)
CSV.write(joinpath(outdir, "bayes_conditional_logbf_pairwise.csv"), bf_long)

ranking = select(logml_mean, :model, :fixed_point, :log_marginal_likelihood,
                 :log_likelihood_at_theta, :log_prior_at_theta,
                 :log_posterior_ordinate, :sigma2_a, :nobs, :n_groups,
                 :n_iter, :burnin, :thin)
sort!(ranking, order(:log_marginal_likelihood, rev=true))
ranking.rank_conditional_logml = collect(1:nrow(ranking))
best = ranking.log_marginal_likelihood[1]
ranking.delta_logml_from_best = ranking.log_marginal_likelihood .- best
ranking.logBF_vs_best = ranking.log_marginal_likelihood .- best
CSV.write(joinpath(outdir, "bayes_conditional_logml_ranking.csv"), ranking)

println("SC6 Bayesian CART-BB finished. Outputs written to: $outdir")
