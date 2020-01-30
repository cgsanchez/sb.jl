using sb, LinearAlgebra, Statistics

# Create a model
# This corresponds to figure 1 in HW's paper
sbp = OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 7.5, α = 0.1, No = 60, β = 5.0)
sbm = SBModel(sbp)

"""

This function integrates the ECEID EOMS

"""
function integrate(nsteps,dt)

    # Create containers starting from a localized state
    ops = ECEIDOps([[1.0 0.0]; [0.0 0.0]],sbm)
    dotops = ECEIDOps(sb.zm,sbm)
    oldops = ECEIDOps(sb.zm,sbm)

    # Bootstrap the integrator
    sb.ECEIDcalcdots!(dotops, ops , dt , sbm)
    sb.ECEIDbootstrap!(ops,oldops,dotops,dt)

    # Create a vector to store energies along the dynamics
    t = 0.0
    ts = Vector{Float64}(undef,nsteps)
    energies = Vector{Float64}(undef,nsteps)

    # Store the value of the energy a t = 0
    ts[1] = t
    energies[1] = sb.ECEIDenergy(ops,sbm)

    # integrate
    for i in 1:nsteps-1
        t = i*dt
        sb.ECEIDcalcdots!(dotops,ops,dt,sbm)
        sb.ECEIDforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        ts[i+1] = t
        energies[i+1]  = sb.ECEIDenergy(ops,sbm)
    end
    return ts,energies
end

ts, energies = integrate(200,0.001)

@test mean(energies .- energies[1]) < 1e-13
