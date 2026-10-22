# -----------------------------------------------------------------------------
# diagnostic_draw_export.jl
# Helper functions to export retained MCMC draws for reviewer diagnostics.
#
# Include this file in your Julia run script, e.g.:
#   include("diagnostic_draw_export.jl")
#
# Expected fit object fields:
#   fit.beta_draws            :: Matrix{Float64}  (S x p)
#   fit.beta_names            :: Vector{String}
#   fit.sigma2_a_draws        :: Vector{Float64}
#   fit.loglik_draws          :: Vector{Float64}  (optional but recommended)
# Optional fit fields:
#   fit.imputation_summary_draws :: Matrix{Float64} (S x q)
#   fit.imputation_summary_names :: Vector{String}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# diagnostic_draw_export.jl
# Export retained MCMC draws for reviewer diagnostics.
# This file is included inside module BayesREProbit.
# -----------------------------------------------------------------------------

function _fit_has(fit, key::AbstractString)
    if fit isa AbstractDict
        return haskey(fit, key) || haskey(fit, Symbol(key))
    else
        return hasproperty(fit, Symbol(key))
    end
end

function _fit_get(fit, key::AbstractString)
    if fit isa AbstractDict
        if haskey(fit, key)
            return fit[key]
        elseif haskey(fit, Symbol(key))
            return fit[Symbol(key)]
        else
            error("Fit object does not contain key '$key'. Available keys: $(collect(keys(fit)))")
        end
    else
        return getproperty(fit, Symbol(key))
    end
end

function _as_string_set(x; upper::Bool = false)
    if x === nothing
        return nothing
    end
    vals = collect(x)
    vals = strip.(string.(vals))
    vals = upper ? uppercase.(vals) : vals
    return Set(vals)
end

function _selected_string(value, allowed; upper::Bool = false)
    allowed === nothing && return true
    allowed_set = _as_string_set(allowed; upper = upper)
    length(allowed_set) == 0 && return true
    val = strip(string(value))
    val = upper ? uppercase(val) : val
    return val in allowed_set
end

function _selected_rep(rep_id::Int, allowed)
    allowed === nothing && return true
    vals = collect(allowed)
    length(vals) == 0 && return true
    vals_int = Int.(vals)
    return rep_id in vals_int
end

function _safe_slug(x)
    s = lowercase(strip(string(x)))
    s = replace(s, r"[^a-z0-9]+" => "_")
    s = replace(s, r"^_+|_+$" => "")
    return isempty(s) ? "unknown" : s
end

function draws_long_from_fit(fit;
                             rep_id::Int,
                             mechanism,
                             workflow,
                             model_name,
                             chain::Int = 1)

    beta_draws = _fit_get(fit, "beta_draws")
    sigma2_a_draws = _fit_get(fit, "sigma2_a_draws")

    beta_names =
        _fit_has(fit, "beta_names") ? string.(_fit_get(fit, "beta_names")) :
        _fit_has(fit, "terms")      ? string.(_fit_get(fit, "terms")) :
        ["beta_$j" for j in 1:size(beta_draws, 2)]

    S, p = size(beta_draws)

    out = DataFrame(
        iter = Int[],
        chain = Int[],
        rep = Int[],
        mechanism = String[],
        workflow = String[],
        model = String[],
        term = String[],
        value = Float64[]
    )

    mech = uppercase(string(mechanism))
    wf = string(workflow)
    mod = string(model_name)

    for s in 1:S
        for j in 1:p
            push!(
                out,
                (
                    iter = s,
                    chain = chain,
                    rep = rep_id,
                    mechanism = mech,
                    workflow = wf,
                    model = mod,
                    term = string(beta_names[j]),
                    value = Float64(beta_draws[s, j])
                )
            )
        end

        push!(
            out,
            (
                iter = s,
                chain = chain,
                rep = rep_id,
                mechanism = mech,
                workflow = wf,
                model = mod,
                term = "sigma2_a",
                value = Float64(sigma2_a_draws[s])
            )
        )

        if _fit_has(fit, "loglik_draws")
            loglik_draws = _fit_get(fit, "loglik_draws")
            push!(
                out,
                (
                    iter = s,
                    chain = chain,
                    rep = rep_id,
                    mechanism = mech,
                    workflow = wf,
                    model = mod,
                    term = "loglik",
                    value = Float64(loglik_draws[s])
                )
            )
        end
    end

    return out
end

function write_diagnostic_draws_if_selected(fit;
                                            diagnostics_dir::AbstractString,
                                            rep_id::Int,
                                            mechanism,
                                            workflow,
                                            model_name,
                                            export_draws::Bool = false,
                                            draw_reps = nothing,
                                            draw_models = nothing,
                                            draw_mechanisms = nothing)

    if !export_draws
        return nothing
    end

    if !_selected_rep(rep_id, draw_reps)
        return nothing
    end

    if !_selected_string(model_name, draw_models; upper = false)
        return nothing
    end

    if !_selected_string(mechanism, draw_mechanisms; upper = true)
        return nothing
    end

    mkpath(diagnostics_dir)

    draws = draws_long_from_fit(
        fit;
        rep_id = rep_id,
        mechanism = mechanism,
        workflow = workflow,
        model_name = model_name,
        chain = 1
    )

    file_name = string(
        "draws_",
        _safe_slug(workflow), "_",
        _safe_slug(mechanism), "_",
        _safe_slug(model_name), "_",
        "rep", lpad(string(rep_id), 3, "0"),
        ".csv"
    )

    out_file = joinpath(diagnostics_dir, file_name)
    CSV.write(out_file, draws)

    println("Diagnostic draws written: $out_file")
    return out_file
end