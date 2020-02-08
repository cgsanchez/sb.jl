using sb, LinearAlgebra, Statistics

# observable storage structure
struct storages
    time :: Vector{Float64}
    energy :: Vector{Float64}
    function storages()
        new(Vector{Float64}(),Vector{Float64}())
    end
end

# storage callback
function store!(storage, t, ops, sbm)
    append!(storage.time,t)
    append!(storage.energy,sb.ehenergy(ops,sbm))
end

function integrate()
    storage = storages()
    # Create a model
    sbm = SBModel(OhmicSBParams(No = 100, α = 0.1))
    # Create containers starting from a localized state
    ops = EhrenfestOps([[1.0 0.0]; [0.0 0.0]],sbm)
    # run dynamics
    RunEhrenfest!(sbm,ops,2000,0.001,store!, storage)
    return storage.time,storage.energy
end

ts, energies = integrate()

@test mean(energies.-energies[1]) < 0.003
