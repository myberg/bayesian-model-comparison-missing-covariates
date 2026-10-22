#!/usr/bin/env julia

include(joinpath(@__DIR__, "src", "BayesREProbit.jl"))

using .BayesREProbit

args = BayesREProbit.parse_cli_args(ARGS)
BayesREProbit.run_from_cli(args)

