julia check_eq.jl temp 216 130 205 &
julia check_eq.jl temp 160 130 205 &
julia check_eq.jl temp 270 130 205 &
julia check_eq.jl temp 216 100 205 &
wait
julia check_eq.jl temp 216 160 205 &
julia check_eq.jl temp 216 130 150 &
julia check_eq.jl temp 216 130 250 &
julia check_eq.jl water 1.32e-5 & # 1 pr μm
wait
julia check_eq.jl water 1.38e-4 &  # 10 pr μm
julia check_eq.jl water 3.63e-4 &  # 25 pr μm 
julia check_eq.jl water 8.05e-4 &  # 50 pr μm
julia check_eq.jl water 2.76e-3 & # 100 pr μm
wait
julia check_eq.jl temp 150 130 205 &
julia check_eq.jl temp 160 130 205 &
julia check_eq.jl temp 170 130 205 &
julia check_eq.jl temp 180 130 205 &
wait
julia check_eq.jl temp 190 130 205 &
julia check_eq.jl temp 200 130 205 &
julia check_eq.jl temp 210 130 205 &
julia check_eq.jl temp 220 130 205 &
wait
julia check_eq.jl temp 230 130 205 &
julia check_eq.jl temp 240 130 205 &
julia check_eq.jl temp 250 130 205 &
julia check_eq.jl temp 260 130 205 &
wait
julia check_eq.jl temp 270 130 205 &

# tropo
julia check_eq.jl temp 216 100 205 &
julia check_eq.jl temp 216 110 205 &
julia check_eq.jl temp 216 120 205 &
wait
julia check_eq.jl temp 216 130 205 &
julia check_eq.jl temp 216 140 205 &
julia check_eq.jl temp 216 150 205 &
julia check_eq.jl temp 216 160 205 &
wait

# exobase
julia check_eq.jl temp 216 130 150 &
julia check_eq.jl temp 216 130 175 &
julia check_eq.jl temp 216 130 200 &
julia check_eq.jl temp 216 130 225 &
wait
julia check_eq.jl temp 216 130 250 &
julia check_eq.jl temp 216 130 275 &
julia check_eq.jl temp 216 130 300 &
julia check_eq.jl temp 216 130 325 &
wait
julia check_eq.jl temp 216 130 350 &
