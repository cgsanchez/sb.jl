
using LinearAlgebra

"""

Structure to hold data for Ehrenfest dynamics

"""
struct EhrenfestOps
    No :: Int64
    ρ  :: Matrix{Complex{Float64}}
    p  :: Vector{Float64}
    q  :: Vector{Float64}
    function EhrenfestOps(ρ0, sbm :: SBModel)
        ρ = zeros(ComplexF64,2,2)
        ρ .= ρ0
        p = zeros(Float64,sbm.No)
        q = zeros(Float64,sbm.No)
        new(sbm.No,ρ,p,q)
    end
end

"""

Leapfrog bootstrap of data (Ehrenfest)

"""
function ehbootstrap!(ops :: EhrenfestOps, oldops :: EhrenfestOps,
                      dotops :: EhrenfestOps, dt :: Float64)
    oldops.ρ .= ops.ρ - dt * dotops.ρ
    oldops.p .= ops.p - dt * dotops.p
    oldops.q .= ops.q - dt * dotops.q
end

"""

Leapfrog step forward (Ehrenfest)

"""
function ehforward!(ops :: EhrenfestOps, oldops :: EhrenfestOps,
                    dotops :: EhrenfestOps, dt :: Float64)
    oldops.ρ .= oldops.ρ + 2 * dt * dotops.ρ
    oldops.p .= oldops.p + 2 * dt * dotops.p
    oldops.q .= oldops.q + 2 * dt * dotops.q
end

"""

Leapfrog Ehrenfest calculate EOM RHS

"""
function ehcalcdots!(dotops :: EhrenfestOps, ops :: EhrenfestOps,
                    dt :: Float64, sbm :: SBModel)
    dotops.ρ .= - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.ρ)
    for i in 1:ops.No
        dotops.p[i] = real(tr(F(sbm,ops.q,i)*ops.ρ))
        dotops.q[i] = ops.p[i]
    end
end

"""

Ehrenfest total energy

"""
function ehenergy(ops :: EhrenfestOps, sbm :: SBModel)
    return real(tr(H(sbm,ops.q,ops.p) * ops.ρ))
end

"""

This function integrates the Ehrenfest EOMs starting from the
initial condition provided in ops for nsteps of dt.
The store! callback function is called at every timestep.

"""
function RunEhrenfest!(sbm, ops, nsteps, dt, store!, storage)

    dotops = EhrenfestOps(zm,sbm)
    oldops = EhrenfestOps(zm,sbm)

    # Bootstrap the integrator
    ehcalcdots!(dotops, ops ,dt , sbm)
    ehbootstrap!(ops,oldops,dotops,dt)

    # call results storage function at t = 0
    t = 0.0
    store!(storage,t,ops,sbm)

    # integrate
    for i in 1:nsteps-1
        t = i*dt
        ehcalcdots!(dotops,ops,dt,sbm)
        ehforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        store!(storage,t,ops,sbm)
    end
    return 1
end
