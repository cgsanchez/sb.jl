using sb

sbp = OhmicSBParams(No = 5000, ωm = 300)

sbm = SBModel(sbp)

omegarange = range(0,sbp.ωm*1.2,length=1001)

dos = sb.densityofstates(sbm.ωs,omegarange)

teodos = sb.ρosc(convert(Vector{Float64},omegarange),sbp)

using Plots

plot(omegarange,dos)

findmax(dos)

omegarange[19]

plot!(omegarange,teodos)

specden = sb.spectraldensity(sbm.cs,sbm.ωs,omegarange)
teospecden = sb.ohmicJ0(omegarange,sbp)

plot(omegarange,specden)
plot!(omegarange,teospecden)

findmax(teospecden)[1]
findmax(specden)

omegarange[105]

sum(abs.(specden[200:600]-teospecden[200:600]))
sum(abs.(dos[200:600]-teodos[200:600]))

plot(omegarange[200:600],specden[200:600]-teospecden[200:600])
plot!(omegarange[200:600],dos[200:600]-teodos[200:600])

# This is testing to construct Ehrenfest

Pkg.activate(".")

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

mean(energies) < 0.003

# This is testing to construct ECEID

Pkg.activate(".")

using sb, LinearAlgebra, Statistics

# Create a model
sbp = OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 7.5, α = 0.1, No = 60, β = 5.0)
sbm = SBModel(sbp)

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

ts, energies = integrate(2000,0.001)

mean(energies .- energies[1]) < 1e-13


















# Strange stuff

function this()
    a = 1
    b = 2
    for i in 1:5
        println(i," ",a,)
        (a,b) = (b,a)
    end
end

this()

for i in 1:5
    println(a)
    a = 5
end
