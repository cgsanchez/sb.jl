
using SafeTestsets

    @time @safetestset "Constant Matrices Tests" begin include("matrices.jl") end
    @time @safetestset "SB Model tests" begin include("sbmodel.jl") end
    @time @safetestset "Ehrenfest dynamics test" begin include("ehrenfest.jl") end
    @time @safetestset "ECEID dynamics test" begin include("ECEID.jl") end
