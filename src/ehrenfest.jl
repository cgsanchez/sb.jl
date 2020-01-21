
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
        ρ = ρ0
        p = zeros(Float64,sbm.No)
        q = zeros(Float64,sbm.No)
        new(sbm.No,ρ,p,q)
    end
end

"""

Leapfrog bootstrap of data

"""
function ehbootstrap!(ops :: EhrenfestOps, oldops :: EhrenfestOps,
                      dotops :: EhrenfestOps, dt :: Float64)
    oldops.ρ .= ops.ρ - dt * dotops.ρ
end

"""

Leapfrog step forward

"""
function ehforward!(ops :: EhrenfestOps, oldops :: EhrenfestOps,
                    dotops :: EhrenfestOps, dt :: Float64)
    oldops.ρ .= oldops.ρ + 2 * dt * dotops.ρ
end

"""

Leapfrog Ehrenfest claculate EOM RHS

"""
function ehcalcdots!(dotops :: EhrenfestOps, ops :: EhrenfestOps,
                    dt :: Float64, sbm :: SBModel)
    dotops.ρ .= - iohbar * comm(H(sbm,ops.q,ops.p),ops.ρ)
end

"""

Ehrenfest total energy

"""
function ehenergy(ops :: EhrenfestOps, sbm :: SBModel)
    return real(tr(H(sbm,ops.q,ops.p) * ops.ρ))
end