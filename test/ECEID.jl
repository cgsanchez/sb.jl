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
    append!(storage.energy  ,sb.ECEIDenergy(ops,sbm))
end

function my_run(thermostat)
    storage = storages()
    # Create a model
    sbm = SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 7.5, α = 0.1, No = 60, β = 5.0))
    # Create containers starting from a localized state
    ops = ECEIDOps(SA[1.0+0.0im 0.0im; 0.0im 0.0im],sbm)
    # run dynamics
    RunECEID!(sbm,ops,2000,0.001,store!, storage, thermostat = thermostat)
    return storage.time,storage.energy
end

ts, energies = my_run(false)
@test mean(energies .- energies[1]) < 1e-13

ts, energies = my_run(true)
@test (energies[2000] .- energies[1]) < 0.0
