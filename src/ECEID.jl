using sb, LinearAlgebra

"""

Structure to hold data for ECEID dynamics

"""
struct ECEIDOps
    No :: Integer
    ρ  :: Array{Complex{Float64},2}
    Cs :: Vector{Array{Complex{Float64},2}}
    Cc :: Vector{Array{Complex{Float64},2}}
    As :: Vector{Array{Complex{Float64},2}}
    Ac :: Vector{Array{Complex{Float64},2}}
    N  :: Vector{Float64}
    p  :: Vector{Float64}
    q  :: Vector{Float64}
    function ECEIDOps(ρ0, sbm :: SBModel)
        Cs = Vector{Array{Complex{Float64},2}}(undef,sbm.No)
        Cc = Vector{Array{Complex{Float64},2}}(undef,sbm.No)
        As = Vector{Array{Complex{Float64},2}}(undef,sbm.No)
        Ac = Vector{Array{Complex{Float64},2}}(undef,sbm.No)
        N = zeros(Float64,sbm.No)
        p = zeros(Float64,sbm.No)
        q = zeros(Float64,sbm.No)
        for i in 1:sbm.No
            Cs[i] = zeros(ComplexF64,2,2)
            Cc[i] = zeros(ComplexF64,2,2)
            Ac[i] = zeros(ComplexF64,2,2)
            As[i] = zeros(ComplexF64,2,2)
            N[i] = 0.0
        end
        new(sbm.No,ρ0,Cs,Cc,As,Ac,N,p,q)
    end
end

"""

Leapfrog bootstrap of data (ECEID)

"""
function ECEIDbootstrap!(ops :: ECEIDOps, oldops :: ECEIDOps,
                      dotops :: ECEIDOps, dt :: Float64)
    oldops.ρ .= ops.ρ - dt * dotops.ρ
    for i in 1:ops.No
        oldops.Cs[i] .= ops.Cs[i] - dt * dotops.Cs[i]
        oldops.Cc[i] .= ops.Cc[i] - dt * dotops.Cc[i]
        oldops.As[i] .= ops.As[i] - dt * dotops.As[i]
        oldops.Ac[i] .= ops.Ac[i] - dt * dotops.Ac[i]
        oldops.N[i]  = ops.N[i] - dt * dotops.N[i]
    end
end

"""

Leapfrog step forward (ECEID)

"""
function ECEIDforward!(ops :: ECEIDOps, oldops :: ECEIDOps,
                    dotops :: ECEIDOps, dt :: Float64)
    oldops.ρ .= oldops.ρ + 2 * dt * dotops.ρ
    Threads.@threads for i in 1:ops.No
        oldops.Cs[i] .= oldops.Cs[i] + 2 * dt * dotops.Cs[i]
        oldops.Cc[i] .= oldops.Cc[i] + 2 * dt * dotops.Cc[i]
        oldops.As[i] .= oldops.As[i] + 2 * dt * dotops.As[i]
        oldops.Ac[i] .= oldops.Ac[i] + 2 * dt * dotops.Ac[i]
        oldops.N[i]  = oldops.N[i] + 2 * dt * dotops.N[i]
    end
end

"""

Leapfrog ECEID calculate EOM RHS

"""
function ECEIDcalcdots!(dotops :: ECEIDOps, ops :: ECEIDOps,
                    dt :: Float64, sbm :: SBModel; thermostat = false)

    tempρ = - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.ρ)
    for i in 1:ops.No
        μ = (im*ops.Cc[i]-ops.As[i])/sbm.ωs[i]
        tempρ += IOHBAR * comm(F(sbm,ops.q,i),μ)
    end
    # Using a temporary variable to set the values of ρ since ECEIDOps struct is inmutable
    dotops.ρ .= tempρ
    Threads.@threads for i in 1:ops.No
        dotops.Cc[i] .= - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.Cc[i]) + sbm.ωs[i] * ops.Cs[i] + (ops.N[i] + 0.5) * comm(F(sbm,ops.q,i),ops.ρ)
        dotops.Cs[i] .= - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.Cs[i]) - sbm.ωs[i] * ops.Cc[i]
        dotops.Ac[i] .= - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.Ac[i]) + sbm.ωs[i] * ops.As[i] + 0.5 * acomm(F(sbm,ops.q,i),ops.ρ)
        dotops.As[i] .= - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.As[i]) - sbm.ωs[i] * ops.Ac[i]
        λ = (im*ops.Cs[i]+ops.Ac[i])
        if !thermostat
            dotops.N[i] = tr(F(sbm,ops.q,i) * λ) / HBAR / sbm.ωs[i]
        else
            dotops.N[i] = 0.0
        end
    end
end

"""

ECEID total energy

"""
function ECEIDenergy(ops :: ECEIDOps, sbm :: SBModel)
    eenergy = real(tr(H(sbm,ops.q,ops.p) * ops.ρ))
    oenergy = 0.0
    intenergy = 0.0
    for i in 1:ops.No
        oenergy += HBAR * sbm.ωs[i] * (ops.N[i] + 0.5)
        μ = (im*ops.Cc[i]-ops.As[i])/sbm.ωs[i]
        intenergy -= real(tr(F(sbm,ops.q,i) * μ))
    end
    return eenergy+oenergy+intenergy
end
