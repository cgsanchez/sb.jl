using sb, LinearAlgebra

# This file implements mean field many body second moment CEID

"""

Structure to hold data for CEID dynamics

"""
struct CEIDOps
    No :: Integer
    ρ  :: Array{Complex{Float64},2}
    μ  :: Vector{Array{Complex{Float64},2}}
    λ  :: Vector{Array{Complex{Float64},2}}
    CRR :: Array{Float64,2}
    CPR :: Array{Float64,2}
    CPP :: Array{Float64,2}
    p  :: Vector{Float64}
    q  :: Vector{Float64}
    ke :: Vector{Float64}
    function CEIDOps(ρ0, sbm :: SBModel)
        μ = Vector{Array{Complex{Float64},2}}(undef,sbm.No)
        λ = Vector{Array{Complex{Float64},2}}(undef,sbm.No)
        CRR = zeros(Float64,sbm.No,sbm.No)
        CPR = zeros(Float64,sbm.No,sbm.No)
        CPP = zeros(Float64,sbm.No,sbm.No)
        p = zeros(Float64,sbm.No)
        q = zeros(Float64,sbm.No)
        ke = zeros(Float64,1)
         for ν in 1:sbm.No
             μ[ν] = zeros(ComplexF64,2,2)
             λ[ν] = zeros(ComplexF64,2,2)
             # Initialize CPP and CRR to expectation values at zero temperature
             # This must be changed to an arbitrary temperature initialization
             nphon = 0
             CPP[ν,ν] = (nphon+0.5) * HBAR * sbm.ωs[ν]
             CRR[ν,ν] = (nphon+0.5) * HBAR / sbm.ωs[ν]
         end
        new(sbm.No,ρ0,μ,λ,CRR,CPR,CPP,p,q,ke)
    end
end

"""

Leapfrog bootstrap of data (CEID)

"""
function CEIDbootstrap!(ops :: CEIDOps, oldops :: CEIDOps,
                      dotops :: CEIDOps, dt :: Float64)
    oldops.ρ .= ops.ρ - dt * dotops.ρ
    for i in 1:ops.No
        oldops.μ[i] .= ops.μ[i] - dt * dotops.μ[i]
        oldops.λ[i] .= ops.λ[i] - dt * dotops.λ[i]
    end
    oldops.CRR .= ops.CRR - dt * dotops.CRR
    oldops.CPR .= ops.CPR - dt * dotops.CPR
    oldops.CPP .= ops.CPP - dt * dotops.CPP
    oldops.p .= ops.p - dt * dotops.p
    oldops.q .= ops.q - dt * dotops.q
    oldops.ke .= ops.ke - dt * dotops.ke
end

"""

Leapfrog step forward (CEID)

"""
function CEIDforward!(ops :: CEIDOps, oldops :: CEIDOps,
                    dotops :: CEIDOps, dt :: Float64)
    oldops.ρ .= oldops.ρ + 2 * dt * dotops.ρ
    for i in 1:ops.No
        oldops.μ[i] .= oldops.μ[i] + 2 * dt * dotops.μ[i]
        oldops.λ[i] .= oldops.λ[i] + 2 * dt * dotops.λ[i]
    end
    oldops.CRR .= oldops.CRR + 2 * dt * dotops.CRR
    oldops.CPR .= oldops.CPR + 2 * dt * dotops.CPR
    oldops.CPP .= oldops.CPP + 2 * dt * dotops.CPP
    oldops.p .= oldops.p + 2 * dt * dotops.p
    oldops.q .= oldops.q + 2 * dt * dotops.q
    oldops.ke .= oldops.ke + 2 * dt * dotops.ke
end

"""

Leapfrog CEID calculate EOM RHS

"""
function CEIDcalcdots!(dotops :: CEIDOps, ops :: CEIDOps,
                    dt :: Float64, sbm :: SBModel; thermostat = false)

    # Time derivatives of generalized coordinates and momenta ------------------
    fbar = zeros(Float64,sbm.No)
    for ν in 1:ops.No
        fbar[ν] = real(tr(F(sbm,ops.q,ν)*ops.ρ))
        for νp in 1:ops.No
            fbar[ν] += - real(tr(K(sbm,ops.q,ν,νp)*ops.μ[νp]))
        end
        dotops.p[ν] = fbar[ν]
        dotops.q[ν] = ops.p[ν]
    end
    #---------------------------------------------------------------------------

    # Time derivative of ρ -----------------------------------------------------
    tempρ = IOHBAR * comm(H(sbm,ops.q,ops.p),ops.ρ)
    for ν in 1:ops.No
        tempρ += - IOHBAR * comm(F(sbm,ops.q,ν),ops.μ[ν])
    end
    for ν in 1:ops.No
        for νp in 1:ops.No
            tempρ += 0.5 * IOHBAR * ops.CRR[ν,νp] * comm(K(sbm,ops.q,ν,νp),ops.ρ)
        end
    end
    # Using a temporary variable to set the values of ρ since CEIDOps struct is inmutable
    dotops.ρ .= tempρ
    #---------------------------------------------------------------------------

    # Time derivative of μ -----------------------------------------------------
    for ν in 1:ops.No
        tempμ = ops.λ[ν] + IOHBAR * comm(H(sbm,ops.q,ops.p),ops.μ[ν])
        for νp in 1:ops.No
            tempμ += - IOHBAR * ops.CRR[ν,νp] * comm(F(sbm,ops.q,νp),ops.ρ)
        end
        dotops.μ[ν] .= tempμ
    end
    #---------------------------------------------------------------------------

    # Time derivative of λ -----------------------------------------------------
    for ν in 1:ops.No
        tempλ = IOHBAR * comm(H(sbm,ops.q,ops.p),ops.λ[ν])
#        tempλ += 0.5 * acomm(F(sbm,ops.q,ν),ops.ρ) - fbar[ν] * ops.ρ
        tempλ += 0.5 * acomm(F(sbm,ops.q,ν),ops.ρ) - tr(F(sbm,ops.q,ν)*ops.ρ)*ops.ρ
        for νp in 1:ops.No
            tempλ += - IOHBAR * ops.CPR[ν,νp] * comm(F(sbm,ops.q,νp),ops.ρ)
            tempλ += - 0.5 * acomm(K(sbm,ops.q,ν,νp),ops.μ[νp])
        end
        dotops.λ[ν] .= tempλ
    end
    #---------------------------------------------------------------------------

    # Time derivatives of mean field second order correlators ------------------
    if !thermostat
        for ν in 1:ops.No
            for νp in 1:ops.No
                dotops.CRR[ν,νp] = ops.CPR[ν,νp] + ops.CPR[νp,ν]
                dotops.CPR[ν,νp] = ops.CPP[ν,νp] + real(tr(F(sbm,ops.q,ν)*ops.μ[νp]))
                dotops.CPP[ν,νp] = real(tr(F(sbm,ops.q,ν)*ops.λ[νp])) +
                                   real(tr(F(sbm,ops.q,νp)*ops.λ[ν]))
                for νpp in 1:ops.No
                    dotops.CPR[ν,νp] += - ops.CRR[νpp,νp] * real(tr(K(sbm,ops.q,ν,νpp)*ops.ρ))
                    dotops.CPP[ν,νp] += - ops.CPR[ν,νpp] * real(tr(K(sbm,ops.q,νpp,νp)*ops.ρ)) -
                                          ops.CPR[νp,νpp] * real(tr(K(sbm,ops.q,νpp,ν)*ops.ρ))
                end
            end
        end
    else
        dotops.CRR .= 0.0
        dotops.CPR .= 0.0
        dotops.CPP .= 0.0
    end
    #---------------------------------------------------------------------------

    # Time derivative of quantum kinetic energy --------------------------------
    dotops.ke[1] = 0.0
    for ν in 1:ops.No
        for νp in 1:ops.No
            dotops.ke[1] += real(tr(K(sbm,ops.q,ν,νp)*ops.ρ)) * ops.CPR[ν,νp] + real(0.5 * IOHBAR * ops.CRR[νp,ν] * tr(K(sbm,ops.q,ν,νp)*comm(H(sbm,ops.q,ops.p),ops.ρ)))
        end
    end
    #---------------------------------------------------------------------------

end

"""

CEID total energy

"""
function CEIDenergy(ops :: CEIDOps, sbm :: SBModel)
    energy = real(tr(H(sbm,ops.q,ops.p)*ops.ρ)) + ops.ke[1]
    for ν in 1:ops.No
        energy += 0.5 * (ops.p[ν]^2 + ops.CPP[ν,ν])
        energy += - real(tr(F(sbm,ops.q,ν)*ops.μ[ν]))
    end
    return energy
end

"""

This function integrates the CEID EOMs starting from the
initial condition provided in ops for nsteps of dt.
The store! callback function is called at every timestep.

"""
function RunCEID!(sbm,ops,nsteps,dt,store!,storage; thermostat = false)

    dotops = CEIDOps(sb.zm(),sbm)
    oldops = CEIDOps(sb.zm(),sbm)

    # Bootstrap the integrator
    sb.CEIDcalcdots!(dotops, ops , dt , sbm; thermostat = thermostat)
    sb.CEIDbootstrap!(ops,oldops,dotops,dt)

    # call results storage function at t = 0
    t = 0.0
    store!(storage,t,ops,sbm)

    # integrate
    for i in 1:nsteps-1
        t = i*dt
        sb.CEIDcalcdots!(dotops,ops,dt,sbm; thermostat = thermostat)
        sb.CEIDforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        store!(storage,t,ops,sbm)
    end

    return 1
end
