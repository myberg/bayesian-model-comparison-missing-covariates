module BayesREProbit

using CSV
using DataFrames
using Random
using LinearAlgebra
using Statistics
using Distributions
using DecisionTree
using Printf

include(joinpath(@__DIR__, "diagnostic_draw_export.jl"))


export parse_cli_args, run_from_cli

# -----------------------------------------------------------------------------
# CLI helpers
# -----------------------------------------------------------------------------
function parse_cli_args(args::Vector{String})
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

function parse_terms(s)
    if ismissing(s) || String(s) == ""
        return String[]
    end
    return clean_token.(String.(split(String(s), ",")))
end

function parse_only_models(s::String)
    isempty(s) && return nothing
    return Set(clean_token.(String.(split(s, ","))))
end

function parse_only_reps(s::String)
    isempty(s) && return nothing
    vals = Int[]
    for x in split(s, ",")
        xx = strip(String(x))
        isempty(xx) && continue
        push!(vals, parse(Int, xx))
    end
    return Set(vals)
end

function parse_only_mechanisms(s::String)
    isempty(s) && return nothing
    return Set(clean_token.(String.(split(s, ","))))
end

function deterministic_offset(s::AbstractString)
    # Julia's built-in hash can be session-dependent.  Use a stable offset for seeds.
    return sum(Int(c) for c in codeunits(String(s)))
end

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------
function read_csv_dataframe(path::AbstractString)
    return CSV.read(path, DataFrame)
end

function numeric_vector_from_col(col)
    return Float64[ismissing(x) ? NaN : Float64(x) for x in col]
end

# DataFrames.jl stores column names as strings internally.  Older versions and
# CSV imports may also carry whitespace or quotes from model_specs.csv.  These
# helpers avoid Symbol-vs-String lookup problems.
clean_token(x) = String(strip(replace(String(x), "\"" => "", "'" => "", "`" => "")))

function find_column_name(dat::DataFrame, requested)
    target = clean_token(requested)
    cols = String.(names(dat))
    idx = findfirst(c -> clean_token(c) == target, cols)
    if idx === nothing
        error("Variable '$target' not found in dataset. Available columns: $(join(cols, ", "))")
    end
    return cols[idx]
end

function make_design_matrix(dat::DataFrame, terms::AbstractVector{<:AbstractString}; yvar="y", idvar="id")
    n = nrow(dat)
    terms = clean_token.(terms)
    p = length(terms) + 1
    X = Matrix{Float64}(undef, n, p)
    X[:, 1] .= 1.0
    colnames = String["(Intercept)"]

    for (j, term) in enumerate(terms)
        col = find_column_name(dat, term)
        X[:, j + 1] .= numeric_vector_from_col(dat[!, col])
        push!(colnames, clean_token(term))
    end

    miss = isnan.(X)
    miss[:, 1] .= false
    return X, miss, colnames
end

function get_y_id(dat::DataFrame; yvar="y", idvar="id")
    ycol = find_column_name(dat, yvar)
    idcol = find_column_name(dat, idvar)

    y = Int.(round.(Float64.(dat[!, ycol])))
    raw_id = Int.(round.(Float64.(dat[!, idcol])))

    levels = sort(unique(raw_id))
    map_id = Dict(v => i for (i, v) in enumerate(levels))
    group = [map_id[v] for v in raw_id]

    return y, group, length(levels)
end

function is_binary_observed(x::AbstractVector{Float64}, miss::AbstractVector{Bool})
    obs = x[.!miss]
    isempty(obs) && return false
    vals = unique(round.(obs; digits=10))
    return all(v -> v == 0.0 || v == 1.0, vals)
end

function initialize_missing!(X::Matrix{Float64}, miss::BitMatrix; binary_cols=nothing)
    n, p = size(X)
    for j in 1:p
        any(miss[:, j]) || continue
        obs = X[.!miss[:, j], j]
        fill_value = isempty(obs) ? 0.0 : mean(obs)
        X[miss[:, j], j] .= fill_value
        if binary_cols !== nothing && binary_cols[j]
            X[miss[:, j], j] .= round.(clamp.(X[miss[:, j], j], 0.0, 1.0))
        end
    end
    return X
end


# -----------------------------------------------------------------------------
# CART-BB Alternative 2: terminal-node donor pools + Bayesian bootstrap
# -----------------------------------------------------------------------------
# This block combines the generic R-export workflow with the user's Julia donor-pool
# idea.  It deliberately removes the z-likelihood weighting component.  The latent
# z and y may still be used as CART predictors to keep the imputation closer to the
# substantive probit model, but donor selection within a terminal node is purely
# Bayesian-bootstrap based.

function weighted_index(rng::AbstractRNG, w::AbstractVector{<:Real})
    s = sum(w)
    if !(isfinite(s)) || s <= 0
        return rand(rng, 1:length(w))
    end
    u = rand(rng) * s
    cs = 0.0
    @inbounds for i in eachindex(w)
        cs += w[i]
        if u <= cs
            return i
        end
    end
    return length(w)
end

mutable struct SimpleTree
    leaf::Bool
    feature::Int
    threshold::Float64
    left::Any
    right::Any
    leafid::Int
end

function impurity(y::AbstractVector{<:Real}, isclass::Bool)
    n = length(y)
    n == 0 && return 0.0
    if isclass
        # Multiclass Gini impurity. This is used both for binary dummies and
        # for joint dummy-pair states such as A/B/C encoded as 0/1/2.
        counts = Dict{Int,Int}()
        @inbounds for v in y
            key = Int(round(v))
            counts[key] = get(counts, key, 0) + 1
        end
        g = 1.0
        for c in values(counts)
            pr = c / n
            g -= pr^2
        end
        return n * g
    else
        m = mean(y)
        return sum((y .- m).^2)
    end
end

function best_split(F::Matrix{Float64}, y::AbstractVector{<:Real};
                    isclass::Bool=false, minleaf::Int=20)
    n, p = size(F)
    base = impurity(y, isclass)
    best_gain = 0.0
    best_j = 0
    best_thr = 0.0

    for j in 1:p
        col = F[:, j]
        qs = unique(quantile(col, [0.2, 0.4, 0.6, 0.8]))
        for thr in qs
            left = col .<= thr
            nl = count(left)
            nr = n - nl
            if nl < minleaf || nr < minleaf
                continue
            end
            gain = base - impurity(y[left], isclass) - impurity(y[.!left], isclass)
            if gain > best_gain
                best_gain = gain
                best_j = j
                best_thr = Float64(thr)
            end
        end
    end
    return best_j, best_thr, best_gain
end

function build_simple_tree(F::Matrix{Float64}, y::AbstractVector{<:Real};
                           isclass::Bool=false, minleaf::Int=20,
                           maxdepth::Int=4, depth::Int=0, leafcounter=Ref(0))
    if depth >= maxdepth || size(F, 1) < 2 * minleaf
        leafcounter[] += 1
        return SimpleTree(true, 0, 0.0, nothing, nothing, leafcounter[])
    end

    j, thr, gain = best_split(F, y; isclass=isclass, minleaf=minleaf)
    if j == 0 || gain <= 1e-12
        leafcounter[] += 1
        return SimpleTree(true, 0, 0.0, nothing, nothing, leafcounter[])
    end

    leftidx = F[:, j] .<= thr
    left = build_simple_tree(F[leftidx, :], y[leftidx]; isclass=isclass,
                             minleaf=minleaf, maxdepth=maxdepth,
                             depth=depth+1, leafcounter=leafcounter)
    right = build_simple_tree(F[.!leftidx, :], y[.!leftidx]; isclass=isclass,
                              minleaf=minleaf, maxdepth=maxdepth,
                              depth=depth+1, leafcounter=leafcounter)
    return SimpleTree(false, j, thr, left, right, 0)
end

function predict_leaf(t::SimpleTree, row::AbstractVector{Float64})
    if t.leaf
        return t.leafid
    elseif row[t.feature] <= t.threshold
        return predict_leaf(t.left, row)
    else
        return predict_leaf(t.right, row)
    end
end

function cartbb_feature_matrix(X::Matrix{Float64}, y::Vector{Int}, z::Vector{Float64}, target::Int)
    return cartbb_feature_matrix_excluding(X, y, z, Set([target]))
end

function cartbb_feature_matrix_excluding(X::Matrix{Float64},
                                         y::Vector{Int},
                                         z::Vector{Float64},
                                         excluded::Set{Int})
    n, p = size(X)
    # Features: y, current latent z, and all current X columns except excluded
    # targets. The intercept column is not informative for splitting and is omitted.
    pred_cols = [j for j in 2:p if !(j in excluded)]
    F = Matrix{Float64}(undef, n, 2 + length(pred_cols))
    F[:, 1] .= Float64.(y)
    F[:, 2] .= z
    for (k, j) in enumerate(pred_cols)
        F[:, 2 + k] .= X[:, j]
    end
    return F
end


# Compatibility overloads: if an older/experimental call site still passes
# group, a and eta, ignore them and use the baseline CART-BB feature set.
# This keeps CART-BB based on y, z and X_except_target only.
function cartbb_feature_matrix(X::Matrix{Float64}, y::Vector{Int}, z::Vector{Float64},
                               group::Vector{Int}, a::Vector{Float64}, eta::Vector{Float64},
                               target::Int)
    return cartbb_feature_matrix(X, y, z, target)
end

function cartbb_feature_matrix_excluding(X::Matrix{Float64}, y::Vector{Int}, z::Vector{Float64},
                                         group::Vector{Int}, a::Vector{Float64}, eta::Vector{Float64},
                                         excluded::Set{Int})
    return cartbb_feature_matrix_excluding(X, y, z, excluded)
end

function known_dummy_pairs(colnames::Vector{String})
    requested = [("x5B", "x5C"), ("x8B", "x8C")]
    out = Tuple{Int,Int}[]
    for (b, c) in requested
        ib = findfirst(==(b), colnames)
        ic = findfirst(==(c), colnames)
        if !(ib === nothing || ic === nothing)
            push!(out, (ib, ic))
        end
    end
    return out
end

dummy_state_from_values(b::Real, c::Real) = b >= 0.5 ? 1 : (c >= 0.5 ? 2 : 0)

function dummy_values_from_state(s::Integer)
    if s == 1
        return 1.0, 0.0
    elseif s == 2
        return 0.0, 1.0
    else
        return 0.0, 0.0
    end
end

function dummy_state_consistent(s::Integer,
                                X::Matrix{Float64},
                                miss::BitMatrix,
                                i::Int, ib::Int, ic::Int)
    bval, cval = dummy_values_from_state(s)
    if !miss[i, ib] && round(X[i, ib]) != bval
        return false
    end
    if !miss[i, ic] && round(X[i, ic]) != cval
        return false
    end
    return true
end

function fallback_dummy_state(X::Matrix{Float64}, miss::BitMatrix, i::Int, ib::Int, ic::Int)
    # Use observed parts of the dummy pair whenever they identify the category.
    if !miss[i, ib] && round(X[i, ib]) == 1.0
        return 1
    elseif !miss[i, ic] && round(X[i, ic]) == 1.0
        return 2
    elseif !miss[i, ib] && !miss[i, ic]
        return dummy_state_from_values(X[i, ib], X[i, ic])
    else
        return 0
    end
end

function enforce_known_dummy_pairs!(X::Matrix{Float64},
                                    colnames::Vector{String};
                                    rng::AbstractRNG=Random.GLOBAL_RNG)
    # The simulated categorical predictors are represented by dummy pairs.
    # This guard avoids impossible states such as x5B=x5C=1.
    for (ib, ic) in known_dummy_pairs(colnames)
        @inbounds for i in 1:size(X, 1)
            X[i, ib] = X[i, ib] >= 0.5 ? 1.0 : 0.0
            X[i, ic] = X[i, ic] >= 0.5 ? 1.0 : 0.0
            if X[i, ib] == 1.0 && X[i, ic] == 1.0
                if rand(rng) < 0.5
                    X[i, ic] = 0.0
                else
                    X[i, ib] = 0.0
                end
            end
        end
    end
    return X
end

function cartbb_update_dummy_pair!(rng::AbstractRNG,
                                   X::Matrix{Float64},
                                   miss::BitMatrix,
                                   y::Vector{Int},
                                   z::Vector{Float64},
                                   colnames::Vector{String},
                                   ib::Int,
                                   ic::Int;
                                   max_depth::Int=4,
                                   min_samples_leaf::Int=20,
                                   max_donors::Int=50)
    # Joint CART-BB imputation for dummy pairs representing one categorical
    # variable, e.g. x5B/x5C = A/B/C. This avoids treating the two dummies as
    # independent continuous variables and respects partially observed pairs.
    target_idx = findall(miss[:, ib] .| miss[:, ic])
    isempty(target_idx) && return X

    train_idx = findall((.!miss[:, ib]) .& (.!miss[:, ic]))
    isempty(train_idx) && return X

    states = Vector{Float64}(undef, size(X, 1))
    @inbounds for i in 1:size(X, 1)
        states[i] = Float64(dummy_state_from_values(X[i, ib], X[i, ic]))
    end

    F = cartbb_feature_matrix_excluding(X, y, z, Set([ib, ic]))
    tree = build_simple_tree(F[train_idx, :], states[train_idx];
                             isclass=true,
                             minleaf=min_samples_leaf,
                             maxdepth=max_depth)

    train_leaf = [predict_leaf(tree, view(F, i, :)) for i in train_idx]
    target_leaf = [predict_leaf(tree, view(F, i, :)) for i in target_idx]

    for (jj, i) in enumerate(target_idx)
        leaf = target_leaf[jj]
        donor_pool_local = train_idx[train_leaf .== leaf]
        if isempty(donor_pool_local)
            donor_pool_local = train_idx
        end

        # Respect observed parts of the dummy pair. If x5B is observed as 1,
        # only state B is allowed; if x5C is observed as 1, only state C is allowed.
        donor_states_all = Int.(round.(states[donor_pool_local]))
        keep = [dummy_state_consistent(s, X, miss, i, ib, ic) for s in donor_states_all]
        donor_states = donor_states_all[keep]

        if isempty(donor_states)
            # Try all training states before using a deterministic fallback.
            all_states = Int.(round.(states[train_idx]))
            keep_all = [dummy_state_consistent(s, X, miss, i, ib, ic) for s in all_states]
            donor_states = all_states[keep_all]
        end

        if isempty(donor_states)
            snew = fallback_dummy_state(X, miss, i, ib, ic)
        else
            if length(donor_states) > max_donors
                perm = randperm(rng, length(donor_states))[1:max_donors]
                donor_states = donor_states[perm]
            end
            wbb = randexp(rng, length(donor_states))
            snew = donor_states[weighted_index(rng, wbb)]
        end

        bval, cval = dummy_values_from_state(snew)
        X[i, ib] = bval
        X[i, ic] = cval
    end

    return X
end

function cartbb_update!(rng::AbstractRNG,
                        X::Matrix{Float64},
                        miss::BitMatrix,
                        y::Vector{Int},
                        z::Vector{Float64},
                        colnames::Vector{String};
                        binary_cols=nothing,
                        max_depth::Int=4,
                        min_samples_leaf::Int=20,
                        max_donors::Int=50)
    n, p = size(X)

    pairs = known_dummy_pairs(colnames)
    dummy_pair_cols = Set{Int}()
    for (ib, ic) in pairs
        push!(dummy_pair_cols, ib)
        push!(dummy_pair_cols, ic)
    end

    for target in 2:p  # skip intercept
        # Dummy pairs such as x5B/x5C are updated jointly below.
        target in dummy_pair_cols && continue
        any(miss[:, target]) || continue

        train_idx = findall(.!miss[:, target])
        miss_idx = findall(miss[:, target])
        isempty(train_idx) && continue
        isempty(miss_idx) && continue

        F = cartbb_feature_matrix(X, y, z, target)
        ytarget = X[:, target]
        isclass = binary_cols !== nothing && binary_cols[target]

        tree = build_simple_tree(F[train_idx, :], ytarget[train_idx];
                                 isclass=isclass,
                                 minleaf=min_samples_leaf,
                                 maxdepth=max_depth)

        train_leaf = [predict_leaf(tree, view(F, i, :)) for i in train_idx]
        miss_leaf = [predict_leaf(tree, view(F, i, :)) for i in miss_idx]

        for (jj, i) in enumerate(miss_idx)
            leaf = miss_leaf[jj]
            donor_pool_local = train_idx[train_leaf .== leaf]
            if isempty(donor_pool_local)
                donor_pool_local = train_idx
            end

            if length(donor_pool_local) > max_donors
                perm = randperm(rng, length(donor_pool_local))[1:max_donors]
                donor_pool_local = donor_pool_local[perm]
            end

            donor_vals = ytarget[donor_pool_local]
            # Bayesian bootstrap donor draw: Exp(1) weights, no likelihood weighting.
            wbb = randexp(rng, length(donor_vals))
            val = donor_vals[weighted_index(rng, wbb)]

            if isclass
                val = round(clamp(val, 0.0, 1.0))
            end
            X[i, target] = Float64(val)
        end
    end

    for (ib, ic) in pairs
        if any(miss[:, ib]) || any(miss[:, ic])
            cartbb_update_dummy_pair!(rng, X, miss, y, z, colnames, ib, ic;
                                      max_depth=max_depth,
                                      min_samples_leaf=min_samples_leaf,
                                      max_donors=max_donors)
        end
    end

    enforce_known_dummy_pairs!(X, colnames; rng=rng)
    return X
end

# Fast approximate CART-assisted imputation.
# It follows the spirit of the user's old Julia code: initial deterministic fill,
# then tree predictions plus stochastic residual draws for originally missing cells.
function cart_update!(rng::AbstractRNG,
                      X::Matrix{Float64},
                      miss::BitMatrix;
                      binary_cols=nothing,
                      max_depth::Int=3,
                      min_samples_leaf::Int=20)
    n, p = size(X)

    for target in 2:p  # skip intercept
        any(miss[:, target]) || continue

        obs_idx = findall(.!miss[:, target])
        mis_idx = findall(miss[:, target])
        isempty(obs_idx) && continue
        isempty(mis_idx) && continue

        pred_cols = [j for j in 1:p if j != target]
        X_train = X[obs_idx, pred_cols]
        y_train = X[obs_idx, target]
        X_test = X[mis_idx, pred_cols]

        tree = DecisionTreeRegressor(max_depth=max_depth, min_samples_leaf=min_samples_leaf)
        DecisionTree.fit!(tree, X_train, y_train)

        pred_test = DecisionTree.predict(tree, X_test)
        pred_train = DecisionTree.predict(tree, X_train)
        resid = y_train .- pred_train

        for k in eachindex(mis_idx)
            draw = pred_test[k]
            if !isempty(resid)
                draw += rand(rng, resid)
            end
            if binary_cols !== nothing && binary_cols[target]
                draw = round(clamp(draw, 0.0, 1.0))
            end
            X[mis_idx[k], target] = draw
        end
    end

    return X
end

# -----------------------------------------------------------------------------
# Gibbs sampler for random-effects probit
# -----------------------------------------------------------------------------
function pointwise_loglik_probit(y::Vector{Int}, eta::Vector{Float64})
    d0 = Normal()
    lli = Vector{Float64}(undef, length(y))
    @inbounds for i in eachindex(y)
        if y[i] == 1
            lli[i] = logcdf(d0, eta[i])
        else
            lli[i] = logccdf(d0, eta[i])
        end
    end
    return lli
end

function loglik_probit(y::Vector{Int}, eta::Vector{Float64})
    return sum(pointwise_loglik_probit(y, eta))
end

function sample_truncated_z!(rng::AbstractRNG,
                             z::Vector{Float64},
                             y::Vector{Int},
                             eta::Vector{Float64})
    @inbounds for i in eachindex(y)
        if y[i] == 1
            z[i] = rand(rng, truncated(Normal(eta[i], 1.0), 0.0, Inf))
        else
            z[i] = rand(rng, truncated(Normal(eta[i], 1.0), -Inf, 0.0))
        end
    end
    return z
end

function update_random_effects!(rng::AbstractRNG,
                                a::Vector{Float64},
                                z::Vector{Float64},
                                Xβ::Vector{Float64},
                                group_indices::Vector{Vector{Int}},
                                sigma2_a::Float64)
    resid = z .- Xβ
    @inbounds for j in eachindex(group_indices)
        idx = group_indices[j]
        nj = length(idx)
        vj = 1.0 / (nj + 1.0 / sigma2_a)
        mj = vj * sum(resid[idx])
        a[j] = rand(rng, Normal(mj, sqrt(vj)))
    end
    return a
end

function run_gibbs_re_probit(dat::DataFrame,
                             terms::AbstractVector{<:AbstractString};
                             model::String,
                             rep::Int,
                             mechanism::String,
                             data_version::String="mis",
                             yvar::String="y",
                             idvar::String="id",
                             n_iter::Int=4000,
                             burnin::Int=1000,
                             thin::Int=5,
                             tau2_beta::Float64=1.0,
                             c0::Float64=2.0,
                             d0::Float64=1.0,
                             seed::Int=202605,
                             impute_mode::String="cart_every",
                             cart_every::Int=25,
                             cart_max_depth::Int=3,
                             cart_min_samples_leaf::Int=20,
                             cart_max_donors::Int=50,
                             store_imputations::Bool=false,
                             chib_reduced_iter::Int=2000,
                             chib_reduced_burnin::Int=500,
                             chib_reduced_thin::Int=2,
                             evidence_target::String="conditional_on_fixed_X",
                             medoid_draw_id::Union{Missing,Int}=missing,
                             verbose_every::Int=500)
    rng = MersenneTwister(seed)

    y, group, J = get_y_id(dat; yvar=yvar, idvar=idvar)
    group_indices = [findall(==(j), group) for j in 1:J]
    X, miss, colnames = make_design_matrix(dat, terms; yvar=yvar, idvar=idvar)
    miss = BitMatrix(miss)
    n, p = size(X)

    binary_cols = [is_binary_observed(X[:, j], miss[:, j]) for j in 1:p]
    initialize_missing!(X, miss; binary_cols=binary_cols)
    enforce_known_dummy_pairs!(X, colnames; rng=rng)

    if impute_mode == "cart_once" || impute_mode == "cart_every"
        # Original fast approximation: tree predictions + residual draw.
        cart_update!(rng, X, miss; binary_cols=binary_cols,
                     max_depth=cart_max_depth,
                     min_samples_leaf=cart_min_samples_leaf)
        enforce_known_dummy_pairs!(X, colnames; rng=rng)
    elseif impute_mode == "cartbb_once" || impute_mode == "cartbb_every"
        # Alternative 2 uses CART terminal-node donor pools.  It needs the current
        # latent z, therefore the first update is performed inside the Gibbs loop.
    elseif impute_mode == "mean_once"
        # mean/dummy initialization only
    elseif impute_mode == "none"
        if any(miss)
            error("impute_mode='none' but X contains missing values.")
        end
    else
        error("Unknown impute_mode = $impute_mode")
    end

    n_save = Int(floor((n_iter - burnin) / thin))
    beta_draws = Matrix{Float64}(undef, n_save, p)
    a_draws = Matrix{Float64}(undef, n_save, J)
    sigma2_draws = Vector{Float64}(undef, n_save)
    loglik_draws = Vector{Float64}(undef, n_save)
    loglik_matrix = Matrix{Float64}(undef, n_save, n)
    z_draws = Matrix{Float64}(undef, n_save, n)

    missing_indices = findall(miss)
    imputation_draws = store_imputations ?
        Matrix{Float64}(undef, n_save, length(missing_indices)) :
        Matrix{Float64}(undef, 0, 0)

    scale_by_column = ones(Float64, p)
    for j in 1:p
        binary_cols[j] && continue
        observed = X[.!miss[:, j], j]
        sj = length(observed) > 1 ? std(observed) : 1.0
        scale_by_column[j] = isfinite(sj) && sj > 1e-8 ? sj : 1.0
    end
    imputation_cell_scales = Float64[scale_by_column[idx[2]] for idx in missing_indices]

    β = zeros(p)
    a = zeros(J)
    sigma2_a = 1.0
    z = zeros(n)

    prior_prec = Matrix{Float64}(I, p, p) ./ tau2_beta
    prior_mean = zeros(p)

    save_id = 0
    elapsed = @elapsed begin
        for iter in 1:n_iter
            eta = X * β .+ a[group]
            sample_truncated_z!(rng, z, y, eta)

            # β | z,a,X
            precision = Symmetric(X'X + prior_prec)
            covariance = inv(precision)
            mean_beta = covariance * (X' * (z .- a[group]) + prior_prec * prior_mean)
            β = rand(rng, MvNormal(mean_beta, Symmetric(covariance)))

            # a | z,β,σ²_a
            Xβ = X * β
            update_random_effects!(rng, a, z, Xβ, group_indices, sigma2_a)

            # σ²_a | a
            shape = c0 + J / 2
            scale = d0 + 0.5 * dot(a, a)
            sigma2_a = rand(rng, InverseGamma(shape, scale))

            # Approximate updates of originally missing covariates.
            if impute_mode == "cart_every" && any(miss) && cart_every > 0 && iter % cart_every == 0
                # Original fast approximation: tree predictions + residual draw.
                cart_update!(rng, X, miss; binary_cols=binary_cols,
                             max_depth=cart_max_depth,
                             min_samples_leaf=cart_min_samples_leaf)
                enforce_known_dummy_pairs!(X, colnames; rng=rng)
            elseif impute_mode == "cartbb_once" && any(miss) && iter == 1
                # Alternative 2: terminal-node donor pools + Bayesian bootstrap,
                # deliberately without likelihood weighting.
                cartbb_update!(rng, X, miss, y, z, colnames;
                               binary_cols=binary_cols,
                               max_depth=cart_max_depth,
                               min_samples_leaf=cart_min_samples_leaf,
                               max_donors=cart_max_donors)
            elseif impute_mode == "cartbb_every" && any(miss) && cart_every > 0 && iter % cart_every == 0
                cartbb_update!(rng, X, miss, y, z, colnames;
                               binary_cols=binary_cols,
                               max_depth=cart_max_depth,
                               min_samples_leaf=cart_min_samples_leaf,
                               max_donors=cart_max_donors)
            end

            if iter > burnin && ((iter - burnin) % thin == 0)
                save_id += 1
                beta_draws[save_id, :] .= β
                a_draws[save_id, :] .= a
                sigma2_draws[save_id] = sigma2_a
                z_draws[save_id, :] .= z
                lli = pointwise_loglik_probit(y, X * β .+ a[group])
                loglik_matrix[save_id, :] .= lli
                loglik_draws[save_id] = sum(lli)
                if store_imputations
                    imputation_draws[save_id, :] .= X[missing_indices]
                end
            end

            if verbose_every > 0 && iter % verbose_every == 0
                @printf("Bayes | rep=%d | %s | %s | iter=%d/%d | sigma2_a=%.3f\n",
                        rep, mechanism, model, iter, n_iter, sigma2_a)
            end
        end
    end

    return Dict(
        "model" => model,
        "rep" => rep,
        "mechanism" => mechanism,
        "terms" => colnames,
        "beta_draws" => beta_draws,
        "a_draws" => a_draws,
        "sigma2_a_draws" => sigma2_draws,
        "z_draws" => z_draws,
        "loglik_draws" => loglik_draws,
        "loglik_matrix" => loglik_matrix,
        "imputation_draws" => imputation_draws,
        "imputation_cell_scales" => imputation_cell_scales,
        "missing_indices" => missing_indices,
        "missing_mask" => copy(miss),
        "X_final" => copy(X),
        "y" => y,
        "group" => group,
        "group_indices" => group_indices,
        "tau2_beta" => tau2_beta,
        "c0" => c0,
        "d0" => d0,
        "runtime_seconds" => elapsed,
        "nobs" => n,
        "n_groups" => J,
        "n_iter" => n_iter,
        "burnin" => burnin,
        "thin" => thin,
        "seed" => seed,
        "impute_mode" => impute_mode,
        "data_version" => data_version,
        "evidence_target" => evidence_target,
        "medoid_draw_id" => medoid_draw_id,
        "chib_reduced_iter" => chib_reduced_iter,
        "chib_reduced_burnin" => chib_reduced_burnin,
        "chib_reduced_thin" => chib_reduced_thin,
        "cart_every" => cart_every,
        "cart_max_depth" => cart_max_depth,
        "cart_min_samples_leaf" => cart_min_samples_leaf,
        "cart_max_donors" => cart_max_donors
    )
end

# Reference chain used only to construct a common medoid completion. The
# substantive probit model is defined by analysis_terms (M2_true by default),
# while imputation_terms contains every covariate needed by any candidate model.
# This keeps auxiliary/noise variables available as CART targets and predictors
# without adding them to the reference outcome model.
function run_medoid_reference_imputation(
        dat::DataFrame,
        analysis_terms::AbstractVector{<:AbstractString},
        imputation_terms::AbstractVector{<:AbstractString};
        rep::Int,
        mechanism::String,
        n_iter::Int=4000,
        burnin::Int=1000,
        thin::Int=10,
        tau2_beta::Float64=1.0,
        c0::Float64=2.0,
        d0::Float64=1.0,
        seed::Int=202605,
        impute_mode::String="cartbb_every",
        cart_every::Int=10,
        cart_max_depth::Int=5,
        cart_min_samples_leaf::Int=10,
        cart_max_donors::Int=50,
        verbose_every::Int=500)
    n_iter > burnin || error("Reference-imputation n_iter must exceed burnin.")
    thin > 0 || error("Reference-imputation thin must be positive.")
    impute_mode in ("cartbb_every", "cart_every") || error(
        "Reference medoid imputation requires cartbb_every or cart_every."
    )

    rng = MersenneTwister(seed)
    y, group, J = get_y_id(dat)
    group_indices = [findall(==(j), group) for j in 1:J]
    X, miss, colnames = make_design_matrix(dat, unique(clean_token.(imputation_terms)))
    miss = BitMatrix(miss)
    n, p_imp = size(X)

    analysis_names = ["(Intercept)"; clean_token.(analysis_terms)]
    analysis_idx = Int[]
    for term in analysis_names
        idx = findfirst(==(term), colnames)
        idx === nothing && error(
            "Reference analysis term '$term' is absent from imputation_terms."
        )
        push!(analysis_idx, idx)
    end

    binary_cols = [is_binary_observed(X[:, j], miss[:, j]) for j in 1:p_imp]
    initialize_missing!(X, miss; binary_cols=binary_cols)
    enforce_known_dummy_pairs!(X, colnames; rng=rng)

    missing_indices = findall(miss)
    n_save = Int(floor((n_iter - burnin) / thin))
    imputation_draws = Matrix{Float64}(undef, n_save, length(missing_indices))

    scale_by_column = ones(Float64, p_imp)
    for j in 1:p_imp
        binary_cols[j] && continue
        observed = X[.!miss[:, j], j]
        sj = length(observed) > 1 ? std(observed) : 1.0
        scale_by_column[j] = isfinite(sj) && sj > 1e-8 ? sj : 1.0
    end
    imputation_cell_scales = Float64[
        scale_by_column[idx[2]] for idx in missing_indices
    ]

    q = length(analysis_idx)
    beta = zeros(q)
    a = zeros(J)
    sigma2_a = 1.0
    z = zeros(n)
    prior_prec = Matrix{Float64}(I, q, q) ./ tau2_beta
    save_id = 0

    elapsed = @elapsed begin
        for iter in 1:n_iter
            X_analysis = X[:, analysis_idx]
            eta = X_analysis * beta .+ a[group]
            sample_truncated_z!(rng, z, y, eta)

            covariance = inv(Symmetric(X_analysis' * X_analysis + prior_prec))
            mean_beta = covariance * (X_analysis' * (z .- a[group]))
            beta = rand(rng, MvNormal(mean_beta, Symmetric(covariance)))
            Xbeta = X_analysis * beta
            update_random_effects!(rng, a, z, Xbeta, group_indices, sigma2_a)

            shape = c0 + J / 2
            scale = d0 + 0.5 * dot(a, a)
            sigma2_a = rand(rng, InverseGamma(shape, scale))

            if any(miss) && cart_every > 0 && iter % cart_every == 0
                if impute_mode == "cartbb_every"
                    cartbb_update!(rng, X, miss, y, z, colnames;
                                   binary_cols=binary_cols,
                                   max_depth=cart_max_depth,
                                   min_samples_leaf=cart_min_samples_leaf,
                                   max_donors=cart_max_donors)
                else
                    cart_update!(rng, X, miss;
                                 binary_cols=binary_cols,
                                 max_depth=cart_max_depth,
                                 min_samples_leaf=cart_min_samples_leaf)
                    enforce_known_dummy_pairs!(X, colnames; rng=rng)
                end
            end

            if iter > burnin && ((iter - burnin) % thin == 0)
                save_id += 1
                imputation_draws[save_id, :] .= X[missing_indices]
            end

            if verbose_every > 0 && iter % verbose_every == 0
                @printf("Medoid reference | rep=%d | %s | iter=%d/%d | sigma2_a=%.3f\n",
                        rep, mechanism, iter, n_iter, sigma2_a)
            end
        end
    end

    return Dict(
        "X_final" => copy(X),
        "imputation_draws" => imputation_draws,
        "imputation_cell_scales" => imputation_cell_scales,
        "missing_indices" => missing_indices,
        "missing_mask" => copy(miss),
        "analysis_terms" => analysis_names,
        "imputation_terms" => colnames,
        "runtime_seconds" => elapsed,
        "seed" => seed,
        "n_iter" => n_iter,
        "burnin" => burnin,
        "thin" => thin
    )
end

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------
function method_label_from_mode(mode::AbstractString, data_version::AbstractString)
    if data_version == "medoid"
        return "CART-BB medoid + conditional Bayesian MCMC (Julia)"
    elseif mode == "cartbb_once" || mode == "cartbb_every"
        return "Bayesian CART-BB MCMC (Julia)"
    elseif mode == "mean_once"
        return "Bayesian mean-initialized MCMC (Julia)"
    elseif mode == "none" && data_version == "cc"
        return "Complete case + Bayesian MCMC (Julia)"
    elseif mode == "none" && data_version == "full"
        return "Bayesian full-data MCMC (Julia)"
    elseif mode == "none"
        return "Bayesian no-imputation MCMC (Julia)"
    else
        return "Bayesian CART-assisted MCMC (Julia)"
    end
end

method_label(fit::Dict) = method_label_from_mode(
    String(fit["impute_mode"]),
    String(fit["data_version"])
)

# -----------------------------------------------------------------------------
# Common medoid imputation
# -----------------------------------------------------------------------------
# The medoid is one retained completed design, chosen from a single common
# reference-imputation run. Candidate models must all use this same completion.
# Distances are Euclidean after scaling continuous cells by the observed-column
# standard deviation; squared differences of binary cells equal Hamming distance.

function select_medoid_imputation(fit::Dict)
    draws = fit["imputation_draws"]
    size(draws, 1) > 0 || error("No imputation draws were stored for medoid selection.")

    if size(draws, 2) == 0
        return (X=copy(fit["X_final"]), draw_id=1, total_distance=0.0)
    end

    scales = fit["imputation_cell_scales"]
    Z = draws ./ reshape(scales, 1, :)
    Z .-= mean(Z, dims=1)
    gram = Z * transpose(Z)
    norms = diag(gram)
    d2 = max.(reshape(norms, :, 1) .+ reshape(norms, 1, :) .- 2.0 .* gram, 0.0)
    distances = sqrt.(d2)
    total_distance = vec(sum(distances, dims=2))
    draw_id = argmin(total_distance)

    X = copy(fit["X_final"])
    X[fit["missing_indices"]] .= vec(draws[draw_id, :])
    return (X=X, draw_id=draw_id, total_distance=total_distance[draw_id])
end

function completed_dataframe_from_design(dat::DataFrame,
                                         terms::AbstractVector{<:AbstractString},
                                         X_completed::Matrix{Float64})
    out = copy(dat)
    for (j, term) in enumerate(clean_token.(terms))
        col = find_column_name(out, term)
        missing_rows = ismissing.(out[!, col])
        out[missing_rows, col] = X_completed[missing_rows, j + 1]
    end
    return out
end

# -----------------------------------------------------------------------------
# Conditional marginal likelihood via Chib (1995) and reduced Gibbs sampling
# -----------------------------------------------------------------------------
# Target:
#   p(y | X_medoid, M) = p(y | theta*,X) p(theta*) / p(theta* | y,X).
# The posterior ordinate is factorized as
#   p(beta*|y,X) p(a*|beta*,y,X) p(sigma2*|a*,beta*,y,X).
# The first ordinate is averaged over the full Gibbs run. The second is averaged
# over a reduced run with beta fixed at beta*. The third is available exactly.

function colmean_vec(M::Matrix{Float64})
    return vec(mean(M, dims=1))
end

function colmedian_vec(M::Matrix{Float64})
    return [median(M[:, j]) for j in 1:size(M, 2)]
end

function log_prior_theta(beta::Vector{Float64},
                         a::Vector{Float64},
                         sigma2_a::Float64,
                         tau2_beta::Float64,
                         c0::Float64,
                         d0::Float64)
    lp = sum(logpdf.(Normal(0.0, sqrt(tau2_beta)), beta))
    lp += sum(logpdf.(Normal(0.0, sqrt(sigma2_a)), a))
    lp += logpdf(InverseGamma(c0, d0), sigma2_a)
    return lp
end

function beta_conditional_parameters(X::Matrix{Float64},
                                     z::AbstractVector{<:Real},
                                     a::AbstractVector{<:Real},
                                     group::Vector{Int},
                                     prior_prec::Matrix{Float64},
                                     prior_mean::Vector{Float64})
    precision = Symmetric(X'X + prior_prec)
    covariance = inv(precision)
    mean_beta = covariance * (X' * (z .- a[group]) + prior_prec * prior_mean)
    return mean_beta, Symmetric(covariance)
end

function log_beta_ordinate(beta_star::Vector{Float64}, fit::Dict)
    X = fit["X_final"]
    z_draws = fit["z_draws"]
    a_draws = fit["a_draws"]
    p = size(X, 2)
    prior_prec = Matrix{Float64}(I, p, p) ./ Float64(fit["tau2_beta"])
    prior_mean = zeros(p)
    covariance = Symmetric(inv(Symmetric(X'X + prior_prec)))
    logdens = Vector{Float64}(undef, size(z_draws, 1))

    for s in axes(z_draws, 1)
        mean_beta = covariance * (
            X' * (view(z_draws, s, :) .- view(a_draws, s, :)[fit["group"]]) +
            prior_prec * prior_mean
        )
        logdens[s] = logpdf(MvNormal(mean_beta, covariance), beta_star)
    end
    return logmeanexp(logdens)
end

function log_random_effect_conditional(a_star::Vector{Float64},
                                       z::Vector{Float64},
                                       Xbeta::Vector{Float64},
                                       group_indices::Vector{Vector{Int}},
                                       sigma2_a::Float64)
    value = 0.0
    residual = z .- Xbeta
    for j in eachindex(group_indices)
        idx = group_indices[j]
        vj = 1.0 / (length(idx) + 1.0 / sigma2_a)
        mj = vj * sum(residual[idx])
        value += logpdf(Normal(mj, sqrt(vj)), a_star[j])
    end
    return value
end

function log_a_ordinate_reduced(a_star::Vector{Float64},
                                beta_star::Vector{Float64},
                                sigma2_star::Float64,
                                fit::Dict;
                                fixed_point::String)
    n_iter = Int(fit["chib_reduced_iter"])
    burnin = Int(fit["chib_reduced_burnin"])
    thin = Int(fit["chib_reduced_thin"])
    n_iter > burnin || error("Reduced-Gibbs n_iter must exceed burnin.")
    thin > 0 || error("Reduced-Gibbs thin must be positive.")

    seed_offset = deterministic_offset("reduced_" * fixed_point)
    rng = MersenneTwister(Int(fit["seed"]) + seed_offset)
    X = fit["X_final"]
    y = fit["y"]
    group = fit["group"]
    group_indices = fit["group_indices"]
    J = Int(fit["n_groups"])
    c0 = Float64(fit["c0"])
    d0 = Float64(fit["d0"])

    a = copy(a_star)
    sigma2_a = sigma2_star
    z = ifelse.(y .== 1, 0.5, -0.5)
    Xbeta = X * beta_star
    logdens = Float64[]

    for iter in 1:n_iter
        eta = Xbeta .+ a[group]
        sample_truncated_z!(rng, z, y, eta)
        update_random_effects!(rng, a, z, Xbeta, group_indices, sigma2_a)
        shape = c0 + J / 2
        scale = d0 + 0.5 * dot(a, a)
        sigma2_a = rand(rng, InverseGamma(shape, scale))

        if iter > burnin && ((iter - burnin) % thin == 0)
            push!(logdens, log_random_effect_conditional(
                a_star, z, Xbeta, group_indices, sigma2_a
            ))
        end
    end
    isempty(logdens) && error("Reduced-Gibbs run retained no draws.")
    return logmeanexp(logdens), length(logdens)
end

function fixed_point_theta(fit::Dict, kind::String)
    beta_draws = fit["beta_draws"]
    a_draws = fit["a_draws"]
    sigma2_draws = fit["sigma2_a_draws"]

    if kind == "mean"
        beta = colmean_vec(beta_draws)
        a = colmean_vec(a_draws)
        sigma2_a = mean(sigma2_draws)
    elseif kind == "median"
        beta = colmedian_vec(beta_draws)
        a = colmedian_vec(a_draws)
        sigma2_a = median(sigma2_draws)
    elseif kind == "high_density"
        # Select the retained draw with the largest approximate log posterior
        # ordinate contribution log f(y|theta_s) + log pi(theta_s).
        ll = fit["loglik_draws"]
        tau2_beta = Float64(fit["tau2_beta"])
        c0 = Float64(fit["c0"])
        d0 = Float64(fit["d0"])
        scores = Vector{Float64}(undef, length(ll))
        @inbounds for s in eachindex(ll)
            scores[s] = ll[s] + log_prior_theta(vec(beta_draws[s, :]), vec(a_draws[s, :]),
                                                 sigma2_draws[s], tau2_beta, c0, d0)
        end
        idx = argmax(scores)
        beta = vec(beta_draws[idx, :])
        a = vec(a_draws[idx, :])
        sigma2_a = sigma2_draws[idx]
    else
        error("Unknown Chib fixed point '$kind'. Use mean, median, or high_density.")
    end
    return beta, a, sigma2_a
end

function chib_marginal_likelihood(fit::Dict; fixed_point::String="mean")
    any(fit["missing_mask"]) && error(
        "Reduced-Gibbs Chib evidence requires a fixed completed design. " *
        "Select one common medoid completion first and fit with impute_mode='none'."
    )
    beta, a, sigma2_a = fixed_point_theta(fit, fixed_point)
    X = fit["X_final"]
    y = fit["y"]
    group = fit["group"]
    tau2_beta = Float64(fit["tau2_beta"])
    c0 = Float64(fit["c0"])
    d0 = Float64(fit["d0"])

    ll = loglik_probit(y, X * beta .+ a[group])
    lp = log_prior_theta(beta, a, sigma2_a, tau2_beta, c0, d0)
    lp_beta = log_beta_ordinate(beta, fit)
    lp_a, reduced_draws = log_a_ordinate_reduced(a, beta, sigma2_a, fit;
                                                  fixed_point=fixed_point)
    lp_sigma = logpdf(
        InverseGamma(c0 + length(a) / 2, d0 + 0.5 * dot(a, a)),
        sigma2_a
    )
    lpost = lp_beta + lp_a + lp_sigma
    lm = ll + lp - lpost

    return (
        fixed_point = fixed_point,
        log_likelihood_at_theta = ll,
        log_prior_at_theta = lp,
        log_posterior_ordinate = lpost,
        log_ordinate_beta = lp_beta,
        log_ordinate_a_given_beta = lp_a,
        log_ordinate_sigma2_given_beta_a = lp_sigma,
        reduced_draws = reduced_draws,
        log_marginal_likelihood = lm,
        beta_intercept = beta[1],
        sigma2_a = sigma2_a
    )
end

function marginal_likelihood_dataframe(fit::Dict)
    rows = DataFrame()
    for fp in ("mean", "median", "high_density")
        cm = chib_marginal_likelihood(fit; fixed_point=fp)
        push!(rows, (
            method = method_label(fit),
            model = fit["model"],
            rep = fit["rep"],
            mechanism = fit["mechanism"],
            fixed_point = cm.fixed_point,
            log_likelihood_at_theta = cm.log_likelihood_at_theta,
            log_prior_at_theta = cm.log_prior_at_theta,
            log_posterior_ordinate = cm.log_posterior_ordinate,
            log_ordinate_beta = cm.log_ordinate_beta,
            log_ordinate_a_given_beta = cm.log_ordinate_a_given_beta,
            log_ordinate_sigma2_given_beta_a = cm.log_ordinate_sigma2_given_beta_a,
            reduced_draws = cm.reduced_draws,
            log_marginal_likelihood = cm.log_marginal_likelihood,
            beta_intercept = cm.beta_intercept,
            sigma2_a = cm.sigma2_a,
            impute_mode = fit["impute_mode"],
            data_version = fit["data_version"],
            evidence_target = fit["evidence_target"],
            medoid_draw_id = fit["medoid_draw_id"],
            tau2_beta = fit["tau2_beta"],
            nobs = fit["nobs"],
            n_groups = fit["n_groups"],
            n_iter = fit["n_iter"],
            burnin = fit["burnin"],
            thin = fit["thin"]
        ))
    end
    return rows
end

# -----------------------------------------------------------------------------
# Predictive criteria from the retained pointwise log-likelihood matrix
# -----------------------------------------------------------------------------
# The matrix has one row per retained posterior draw and one column per observation.
# WAIC is computed as -2 * sum_i (logmeanexp_s l_{si} - var_s(l_{si})).
# LOO is a raw importance-sampling approximation computed as
# elpd_loo_i = -logmeanexp_s(-l_{si}).  This is not PSIS-smoothed, but is useful
# as a lightweight diagnostic/model-ranking criterion for the simulation workflow.

function logmeanexp(v::AbstractVector{<:Real})
    m = maximum(v)
    if !isfinite(m)
        return m
    end
    return m + log(mean(exp.(v .- m)))
end

function predictive_criteria(fit::Dict)
    L = fit["loglik_matrix"]
    S, n = size(L)

    lppd = 0.0
    p_waic = 0.0
    elpd_loo = 0.0

    @inbounds for i in 1:n
        li = view(L, :, i)
        lppd_i = logmeanexp(li)
        v_i = S > 1 ? var(collect(li)) : 0.0
        loo_i = -logmeanexp(-collect(li))

        lppd += lppd_i
        p_waic += v_i
        elpd_loo += loo_i
    end

    elpd_waic = lppd - p_waic
    waic = -2.0 * elpd_waic
    looic = -2.0 * elpd_loo
    p_loo = lppd - elpd_loo

    return (
        lppd = lppd,
        p_waic = p_waic,
        elpd_waic = elpd_waic,
        WAIC = waic,
        elpd_loo = elpd_loo,
        p_loo = p_loo,
        LOO = looic,
        nobs = n,
        n_draws = S
    )
end

function predictive_criteria_dataframe(fit::Dict)
    pc = predictive_criteria(fit)
    return DataFrame(
        method = [method_label(fit)],
        model = [fit["model"]],
        rep = [fit["rep"]],
        mechanism = [fit["mechanism"]],
        data_version = [fit["data_version"]],
        impute_mode = [fit["impute_mode"]],
        lppd = [pc.lppd],
        p_waic = [pc.p_waic],
        elpd_waic = [pc.elpd_waic],
        WAIC = [pc.WAIC],
        elpd_loo = [pc.elpd_loo],
        p_loo = [pc.p_loo],
        LOO = [pc.LOO],
        nobs = [pc.nobs],
        n_draws = [pc.n_draws],
        n_iter = [fit["n_iter"]],
        burnin = [fit["burnin"]],
        thin = [fit["thin"]]
    )
end

function summarize_fit(fit::Dict;
                       marglik::Union{Nothing,DataFrame}=nothing,
                       compute_evidence::Bool=true)
    beta_draws = fit["beta_draws"]
    sigma2 = fit["sigma2_a_draws"]
    terms = fit["terms"]

    est = DataFrame()
    for j in 1:size(beta_draws, 2)
        draws = beta_draws[:, j]
        push!(est, (
            method = method_label(fit),
            model = fit["model"],
            rep = fit["rep"],
            mechanism = fit["mechanism"],
            term = terms[j],
            posterior_mean = mean(draws),
            posterior_median = median(draws),
            posterior_sd = std(draws),
            ci_lower = quantile(draws, 0.025),
            ci_upper = quantile(draws, 0.975),
            component = "fixed",
            impute_mode = fit["impute_mode"],
            data_version = fit["data_version"]
        ))
    end

    push!(est, (
        method = method_label(fit),
        model = fit["model"],
        rep = fit["rep"],
        mechanism = fit["mechanism"],
        term = "sigma2_a",
        posterior_mean = mean(sigma2),
        posterior_median = median(sigma2),
        posterior_sd = std(sigma2),
        ci_lower = quantile(sigma2, 0.025),
        ci_upper = quantile(sigma2, 0.975),
        component = "random_effect_variance",
        impute_mode = fit["impute_mode"],
        data_version = fit["data_version"]
    ))

    ll = fit["loglik_draws"]
    if compute_evidence && marglik === nothing
        marglik = marginal_likelihood_dataframe(fit)
    end
    mean_row = compute_evidence ? marglik[marglik.fixed_point .== "mean", :] : nothing
    if compute_evidence
        nrow(mean_row) == 1 || error("Expected exactly one mean fixed-point evidence row.")
    end
    pc = predictive_criteria(fit)
    modelres = DataFrame(
        method = [method_label(fit)],
        model = [fit["model"]],
        rep = [fit["rep"]],
        mechanism = [fit["mechanism"]],
        converged = [true],
        error = [missing],
        logLik = [median(ll)],
        AIC = [missing],
        BIC = [missing],
        nobs = [fit["nobs"]],
        runtime_seconds = [fit["runtime_seconds"]],
        logLik_mean = [mean(ll)],
        logLik_median = [median(ll)],
        log_marginal_likelihood = [compute_evidence ? mean_row.log_marginal_likelihood[1] : missing],
        chib_fixed_point = [compute_evidence ? "mean" : missing],
        log_likelihood_at_theta = [compute_evidence ? mean_row.log_likelihood_at_theta[1] : missing],
        log_prior_at_theta = [compute_evidence ? mean_row.log_prior_at_theta[1] : missing],
        log_posterior_ordinate = [compute_evidence ? mean_row.log_posterior_ordinate[1] : missing],
        evidence_target = [fit["evidence_target"]],
        medoid_draw_id = [fit["medoid_draw_id"]],
        lppd = [pc.lppd],
        p_waic = [pc.p_waic],
        elpd_waic = [pc.elpd_waic],
        WAIC = [pc.WAIC],
        elpd_loo = [pc.elpd_loo],
        p_loo = [pc.p_loo],
        LOO = [pc.LOO],
        n_iter = [fit["n_iter"]],
        burnin = [fit["burnin"]],
        thin = [fit["thin"]],
        impute_mode = [fit["impute_mode"]],
        data_version = [fit["data_version"]],
        cart_every = [fit["cart_every"]],
        cart_max_depth = [fit["cart_max_depth"]],
        cart_min_samples_leaf = [fit["cart_min_samples_leaf"]],
        cart_max_donors = [fit["cart_max_donors"]]
    )

    return modelres, est
end

function run_from_cli(args::Dict{String,String})
    index_file = get_arg(args, "index", "data/simulated/simulation_index.csv")
    model_file = get_arg(args, "models", "data/simulated/model_specs.csv")
    outdir = get_arg(args, "outdir", "output/julia_bayes")

    n_iter = parse(Int, get_arg(args, "n_iter", "4000"))
    burnin = parse(Int, get_arg(args, "burnin", "1000"))
    thin = parse(Int, get_arg(args, "thin", "5"))
    tau2_beta = parse(Float64, get_arg(args, "tau2_beta", "1.0"))
    c0 = parse(Float64, get_arg(args, "c0", "2.0"))
    d0 = parse(Float64, get_arg(args, "d0", "1.0"))
    seed = parse(Int, get_arg(args, "seed", "202605"))
    impute_mode = get_arg(args, "impute_mode", "cart_every")
    conditional_medoid = lowercase(get_arg(args, "conditional_medoid", "false")) in ("true", "1", "yes")
    medoid_reference_model = get_arg(args, "medoid_reference_model", "M2_true")
    medoid_n_iter = parse(Int, get_arg(args, "medoid_n_iter", string(n_iter)))
    medoid_burnin = parse(Int, get_arg(args, "medoid_burnin", string(burnin)))
    medoid_thin = parse(Int, get_arg(args, "medoid_thin", string(thin)))
    chib_reduced_iter = parse(Int, get_arg(args, "chib_reduced_iter", "2000"))
    chib_reduced_burnin = parse(Int, get_arg(args, "chib_reduced_burnin", "500"))
    chib_reduced_thin = parse(Int, get_arg(args, "chib_reduced_thin", "2"))
    compute_evidence = lowercase(get_arg(args, "compute_evidence", "true")) in ("true", "1", "yes")
    data_version = lowercase(get_arg(args, "data_version", "mis"))
    if !(data_version in ("mis", "full", "cc"))
        error("Unknown data_version = $data_version. Use 'mis', 'full', or 'cc'.")
    end
    cart_every = parse(Int, get_arg(args, "cart_every", "25"))
    cart_max_depth = parse(Int, get_arg(args, "cart_max_depth", "3"))
    cart_min_samples_leaf = parse(Int, get_arg(args, "cart_min_samples_leaf", "20"))
    cart_max_donors = parse(Int, get_arg(args, "cart_max_donors", "50"))
    verbose_every = parse(Int, get_arg(args, "verbose_every", "500"))
    only_models = parse_only_models(get_arg(args, "only_models", ""))
    only_reps = parse_only_reps(get_arg(args, "only_reps", ""))
    only_mechanisms = parse_only_mechanisms(get_arg(args, "only_mechanisms", ""))

    export_draws = lowercase(get_arg(args, "export_draws", "false")) in ("true", "1", "yes")
    diagnostics_dir = get_arg(args, "diagnostics_dir", joinpath(outdir, "mcmc_diagnostics"))

    draw_reps_arg = get_arg(args, "draw_reps", "1")
    draw_models_arg = get_arg(args, "draw_models", "M2_true")
    draw_mechanisms_arg = get_arg(args, "draw_mechanisms", "MCAR,MAR")

    draw_reps = parse_only_reps(draw_reps_arg)
    draw_models = parse_only_models(draw_models_arg)
    draw_mechanisms = parse_only_mechanisms(draw_mechanisms_arg)

    if export_draws
        mkpath(diagnostics_dir)
    end

    mkpath(outdir)

    index = read_csv_dataframe(index_file)
    models = read_csv_dataframe(model_file)

    all_model_results = DataFrame(
        method = String[], model = String[], rep = Int[], mechanism = String[],
        converged = Bool[], error = Union{Missing,String}[], logLik = Union{Missing,Float64}[],
        AIC = Union{Missing,Float64}[], BIC = Union{Missing,Float64}[], nobs = Int[],
        runtime_seconds = Union{Missing,Float64}[], logLik_mean = Union{Missing,Float64}[],
        logLik_median = Union{Missing,Float64}[], log_marginal_likelihood = Union{Missing,Float64}[],
        chib_fixed_point = Union{Missing,String}[], log_likelihood_at_theta = Union{Missing,Float64}[],
        log_prior_at_theta = Union{Missing,Float64}[], log_posterior_ordinate = Union{Missing,Float64}[],
        lppd = Union{Missing,Float64}[], p_waic = Union{Missing,Float64}[],
        elpd_waic = Union{Missing,Float64}[], WAIC = Union{Missing,Float64}[],
        elpd_loo = Union{Missing,Float64}[], p_loo = Union{Missing,Float64}[],
        LOO = Union{Missing,Float64}[],
        n_iter = Union{Missing,Int}[], burnin = Union{Missing,Int}[], thin = Union{Missing,Int}[],
        impute_mode = Union{Missing,String}[], data_version = Union{Missing,String}[],
        cart_every = Union{Missing,Int}[],
        cart_max_depth = Union{Missing,Int}[],
        cart_min_samples_leaf = Union{Missing,Int}[], cart_max_donors = Union{Missing,Int}[]
    )
    all_estimates = DataFrame(
        method = String[], model = String[], rep = Int[], mechanism = String[], term = String[],
        posterior_mean = Float64[], posterior_median = Float64[], posterior_sd = Float64[],
        ci_lower = Float64[], ci_upper = Float64[], component = String[],
        impute_mode = String[], data_version = String[]
    )
    all_marglik = DataFrame()
    all_predictive = DataFrame()
    all_medoid_metadata = DataFrame()

    for row in eachrow(index)
        rep_id_row = Int(row.rep_id)
        mechanism_row = clean_token(String(row.mechanism))

        if only_reps !== nothing && !(rep_id_row in only_reps)
            continue
        end
        if only_mechanisms !== nothing && !(mechanism_row in only_mechanisms)
            continue
        end

        if data_version == "full"
    data_file = String(row.file_full)
elseif data_version == "mis"
    data_file = String(row.file_mis)
elseif data_version == "cc"
    if !hasproperty(row, :file_cc) || ismissing(row.file_cc) || isempty(String(row.file_cc))
        error("data_version='cc' requested, but column file_cc is missing or empty in simulation_index.csv.")
    end
    data_file = String(row.file_cc)
else
    error("Unknown data_version: $data_version")
end

        dat = read_csv_dataframe(data_file)
        dat_for_models = dat
        fit_data_version = data_version
        fit_impute_mode = impute_mode
        medoid_draw_id = missing

        if conditional_medoid
            data_version == "mis" || error(
                "conditional_medoid=true requires data_version='mis'."
            )
            impute_mode in ("cartbb_every", "cart_every") || error(
                "conditional_medoid=true requires a repeated imputation mode."
            )

            reference_rows = models[
                clean_token.(String.(models.model)) .== medoid_reference_model, :
            ]
            nrow(reference_rows) == 1 || error(
                "medoid_reference_model='$medoid_reference_model' was not found uniquely."
            )
            reference_terms = parse_terms(reference_rows.terms[1])

            imputation_terms = String[]
            for model_row in eachrow(models)
                append!(imputation_terms, parse_terms(model_row.terms))
            end
            imputation_terms = unique(imputation_terms)

            println("Julia medoid imputation | rep=$rep_id_row | mechanism=$mechanism_row | reference=$medoid_reference_model")
            reference_seed = seed + 100000 * rep_id_row +
                1000 * deterministic_offset(mechanism_row) +
                deterministic_offset("medoid_" * medoid_reference_model)
            reference_fit = run_medoid_reference_imputation(
                dat,
                reference_terms,
                imputation_terms;
                rep = rep_id_row,
                mechanism = mechanism_row,
                n_iter = medoid_n_iter,
                burnin = medoid_burnin,
                thin = medoid_thin,
                tau2_beta = tau2_beta,
                c0 = c0,
                d0 = d0,
                seed = reference_seed,
                impute_mode = impute_mode,
                cart_every = cart_every,
                cart_max_depth = cart_max_depth,
                cart_min_samples_leaf = cart_min_samples_leaf,
                cart_max_donors = cart_max_donors,
                verbose_every = verbose_every
            )
            medoid = select_medoid_imputation(reference_fit)
            dat_for_models = completed_dataframe_from_design(
                dat, imputation_terms, medoid.X
            )
            medoid_draw_id = medoid.draw_id
            fit_data_version = "medoid"
            fit_impute_mode = "none"

            medoid_dir = joinpath(outdir, "medoid_imputations")
            mkpath(medoid_dir)
            medoid_file = joinpath(
                medoid_dir,
                "rep_$(rep_id_row)_$(mechanism_row)_medoid.csv"
            )
            CSV.write(medoid_file, dat_for_models)
            append!(all_medoid_metadata, DataFrame(
                rep = [rep_id_row],
                mechanism = [mechanism_row],
                reference_model = [medoid_reference_model],
                reference_analysis_terms = [join(reference_terms, ",")],
                imputation_terms = [join(imputation_terms, ",")],
                medoid_draw_id = [medoid.draw_id],
                total_scaled_distance = [medoid.total_distance],
                n_imputation_draws = [size(reference_fit["imputation_draws"], 1)],
                n_missing_cells = [size(reference_fit["imputation_draws"], 2)],
                medoid_file = [medoid_file],
                reference_seed = [reference_seed]
            ), cols=:union)
        end

        for m in eachrow(models)
            model_name = String(m.model)
            if only_models !== nothing && !(model_name in only_models)
                continue
            end

            terms = parse_terms(m.terms)
            rep_id = rep_id_row
            mechanism = mechanism_row

            println("Julia Bayes | rep=$rep_id | mechanism=$mechanism | model=$model_name")

            try
                fit = run_gibbs_re_probit(
                    dat_for_models,
                    terms;
                    model = model_name,
                    rep = rep_id,
                    mechanism = mechanism,
                    data_version = fit_data_version,
                    n_iter = n_iter,
                    burnin = burnin,
                    thin = thin,
                    tau2_beta = tau2_beta,
                    c0 = c0,
                    d0 = d0,
                    seed = seed + 100000 * rep_id + 1000 * deterministic_offset(mechanism) + deterministic_offset(model_name),
                    impute_mode = fit_impute_mode,
                    cart_every = cart_every,
                    cart_max_depth = cart_max_depth,
                    cart_min_samples_leaf = cart_min_samples_leaf,
                    cart_max_donors = cart_max_donors,
                    chib_reduced_iter = chib_reduced_iter,
                    chib_reduced_burnin = chib_reduced_burnin,
                    chib_reduced_thin = chib_reduced_thin,
                    evidence_target = conditional_medoid ?
                        "p(y|X_medoid,M): conditional marginal likelihood" :
                        "p(y|X,M): fixed-design marginal likelihood",
                    medoid_draw_id = medoid_draw_id,
                    verbose_every = verbose_every
                )

                write_diagnostic_draws_if_selected(
                    fit;
                    diagnostics_dir = diagnostics_dir,
                    rep_id = rep_id,
                    mechanism = mechanism,
                    workflow = method_label(fit),
                    model_name = model_name,
                    export_draws = export_draws,
                    draw_reps = draw_reps,
                    draw_models = draw_models,
                    draw_mechanisms = draw_mechanisms
                )

                marglik = compute_evidence ? marginal_likelihood_dataframe(fit) : DataFrame()
                modelres, est = summarize_fit(
                    fit; marglik=marglik, compute_evidence=compute_evidence
                )
                predcrit = predictive_criteria_dataframe(fit)
                append!(all_model_results, modelres, cols=:union)
                append!(all_estimates, est, cols=:union)
                if compute_evidence
                    append!(all_marglik, marglik, cols=:union)
                end
                append!(all_predictive, predcrit, cols=:union)
            catch err
                @printf("ERROR in rep=%d mechanism=%s model=%s: %s\n", rep_id, mechanism, model_name, sprint(showerror, err))
                errrow = DataFrame(
                    method = [method_label_from_mode(fit_impute_mode, fit_data_version)],
                    model = [model_name],
                    rep = [rep_id],
                    mechanism = [mechanism],
                    converged = [false],
                    error = [sprint(showerror, err)],
                    logLik = [missing],
                    AIC = [missing],
                    BIC = [missing],
                    nobs = [nrow(dat_for_models)],
                    runtime_seconds = [missing],
                    impute_mode = [fit_impute_mode],
                    data_version = [fit_data_version]
                )
                append!(all_model_results, errrow, cols=:union)
            end
        end
    end

    CSV.write(joinpath(outdir, "bayes_model_results.csv"), all_model_results)
    CSV.write(joinpath(outdir, "bayes_estimates.csv"), all_estimates)
    if ncol(all_marglik) > 0
        CSV.write(joinpath(outdir, "bayes_marginal_likelihood.csv"), all_marglik)
    end
    CSV.write(joinpath(outdir, "bayes_predictive_criteria.csv"), all_predictive)
    if nrow(all_medoid_metadata) > 0
        CSV.write(joinpath(outdir, "medoid_imputation_metadata.csv"), all_medoid_metadata)
    end

    println("Wrote Julia Bayesian results to $outdir")
end

end # module
