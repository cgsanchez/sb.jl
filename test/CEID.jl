using sb, LinearAlgebra, Statistics

# Create a model
sbp = OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 7.5, α = 0.1, No = 60, β = 5.0)
sbm = SBModel(sbp)

function integrate(nsteps,dt)

    # Create containers starting from a localized state
    ops = CEIDOps([[1.0 0.0]; [0.0 0.0]],sbm)
    dotops = CEIDOps(sb.zm,sbm)
    oldops = CEIDOps(sb.zm,sbm)

    # Bootstrap the integrator
    sb.CEIDcalcdots!(dotops, ops , dt , sbm)
    sb.CEIDbootstrap!(ops,oldops,dotops,dt)

    # Create a vector to store energies along the dynamics
    t = 0.0
    ts = Vector{Float64}(undef,nsteps)
    ene = Vector{Float64}(undef,nsteps)

    # Store the value of the energy a t = 0
    ts[1] = t
    ene[1] = sb.CEIDenergy(ops,sbm)

    # integrate
    for i in 1:nsteps-1
        t = i*dt
        sb.CEIDcalcdots!(dotops,ops,dt,sbm)
        sb.CEIDforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        ts[i+1] = t
        ene[i+1] = sb.CEIDenergy(ops,sbm)
    end

    return ts, ene
end

ts, energies = integrate(100,0.001)

@test mean(energies .- energies[1]) < 0.0001
