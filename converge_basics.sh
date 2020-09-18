#!/bin/bash

# Change "mean" to "min" or "max" to do solar min or max.
julia converge_new_file_250.jl temp 216 130 205 mean &
julia converge_new_file_250.jl temp 160 130 205 mean &
julia converge_new_file_250.jl temp 270 130 205 mean &
julia converge_new_file_250.jl temp 216 100 205 mean &
wait
julia converge_new_file_250.jl temp 216 160 205 mean &
julia converge_new_file_250.jl temp 216 130 150 mean &
julia converge_new_file_250.jl temp 216 130 250 mean &
julia converge_new_file_250.jl water 1.32e-5 mean & # 1 pr μm
wait
julia converge_new_file_250.jl water 1.38e-4 mean &  # 10 pr μm
julia converge_new_file_250.jl water 3.63e-4 mean &  # 25 pr μm 
julia converge_new_file_250.jl water 8.05e-4 mean &  # 50 pr μm
julia converge_new_file_250.jl water 2.76e-3 mean & # 100 pr μm
