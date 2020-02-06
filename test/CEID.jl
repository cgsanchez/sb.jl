using sb, LinearAlgebra, Statistics, StaticArrays

struct storages
    time :: Vector{Float64}
    energy :: Vector{Float64}
    function storages()
        new(Vector{Float64}(),Vector{Float64}())
    end
end

function store!(storage, t, ops, sbm)
    append!(storage.time,t)
    append!(storage.energy  ,sb.CEIDenergy(ops,sbm))
end

function my_run(knondiag, thermostat)
    storage = storages()
    # Create a model
    sbm = SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 7.5, α = 0.1, No = 60, β = 5.0))
    # Create containers starting from a localized state
    ops = CEIDOps(SA[1.0 0.0; 0.0 0.0],sbm)
    # run dynamics
    RunCEID!(sbm,ops,2000,0.001,store!, storage, thermostat = thermostat, knondiag = knondiag)
    return storage.time,storage.energy
end

# test non cross K integration
ts, energies = my_run(false,false)
@test mean(energies .- energies[1]) < 0.0001

# test cross K integration
ts, energies = my_run(true,false)
@test mean(energies .- energies[1]) < 0.0001

ts, energies = my_run(true,true)
@test mean(energies[2000] .- energies[1]) < 0.0
