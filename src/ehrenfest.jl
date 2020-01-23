
using LinearAlgebra

"""

Structure to hold data for Ehrenfest dynamics

"""
struct EhrenfestOps
    No :: Integer
    ρ  :: Matrix{Complex{Float64}}
    p  :: Vector{Float64}
    q  :: Vector{Float64}
    function EhrenfestOps(ρ0, sbm :: SBModel)
        p = zeros(Float64,sbm.No)
        q = zeros(Float64,sbm.No)
        new(sbm.No,ρ0,p,q)
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
