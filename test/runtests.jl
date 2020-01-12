
using SafeTestsets

    @time @safetestset "Constant Matrices Tests" begin include("matrices.jl") end
    @time @safetestset "SB Model tests" begin include("sbmodel.jl") end
