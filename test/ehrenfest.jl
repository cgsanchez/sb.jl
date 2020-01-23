using sb, LinearAlgebra, Statistics

# Create a model
sbp = OhmicSBParams(No = 100, α = 0.1)
sbm = SBModel(sbp)

"""

This function integrates the Ehrenfest EOMs

"""
function integrate(nsteps,dt)

    # Create containers starting from a localized state
    ops = EhrenfestOps([[1.0 0.0]; [0.0 0.0]],sbm)
    dotops = EhrenfestOps(sb.zm,sbm)
    oldops = EhrenfestOps(sb.zm,sbm)

    # Bootstrap the integrator
    sb.ehcalcdots!(dotops, ops ,dt , sbm)
    sb.ehbootstrap!(ops,oldops,dotops,dt)

    # Create a vector to store energies along the dynamics
    t = 0.0
    ts = Vector{Float64}(undef,nsteps)
    energies = Vector{Float64}(undef,nsteps)

    # Store the value of the energy a t = 0
    ts[1] = t
    energies[1] = sb.ehenergy(ops,sbm)

    # integrate
    for i in 1:nsteps-1
        t = i*dt
        sb.ehcalcdots!(dotops,ops,dt,sbm)
        sb.ehforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        ts[i+1] = t
        energies[i+1] = sb.ehenergy(ops,sbm)
    end
    return ts,energies
end

ts, energies = integrate(2000,0.001)

@test mean(energies) < 0.003
