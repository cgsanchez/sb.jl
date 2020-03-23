using sb, LinearAlgebra, StaticArrays

# This file implements mean field many body second moment CEID

"""

Structure to hold data for CEID dynamics

"""
mutable struct CEIDOps
    No :: Int64
    ρ  :: SArray{Tuple{2,2},Complex{Float64},2,4}
    μ  :: Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}
    λ  :: Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}
    CRR :: Array{Float64,2}
    CPR :: Array{Float64,2}
    CPP :: Array{Float64,2}
    p  :: Vector{Float64}
    q  :: Vector{Float64}
    ke :: Float64
    function CEIDOps(ρ0, sbm :: SBModel)
        μ = Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}(undef,sbm.No)
        λ = Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}(undef,sbm.No)
        CRR = zeros(Float64,sbm.No,sbm.No)
        CPR = zeros(Float64,sbm.No,sbm.No)
        CPP = zeros(Float64,sbm.No,sbm.No)
        p = zeros(Float64,sbm.No)
        q = zeros(Float64,sbm.No)
        ke = 0.0
         for ν in 1:sbm.No
             μ[ν] = SA[0.0im 0.0im; 0.0im 0.0im]
             λ[ν] = SA[0.0im 0.0im; 0.0im 0.0im]
             if sbm.β == Inf # zero temperature
                 CPP[ν,ν] = 0.5 * HBAR * sbm.ωs[ν]
                 CRR[ν,ν] = 0.5 * HBAR / sbm.ωs[ν]
             else
                 hav = HBAR*sbm.ωs[ν]*(0.5 + 1.0/(exp(HBAR*sbm.ωs[ν]*sbm.β)-1.0))
                 CPP[ν,ν] = hav
                 CRR[ν,ν] = hav / sbm.ωs[ν]^2
             end
         end
        new(sbm.No,ρ0,μ,λ,CRR,CPR,CPP,p,q,ke)
    end
end

"""

Leapfrog bootstrap of data (CEID)

"""
function CEIDbootstrap!(ops :: CEIDOps, oldops :: CEIDOps,
                      dotops :: CEIDOps, dt :: Float64)
    oldops.ρ = ops.ρ - dt * dotops.ρ
    for i in 1:ops.No
        oldops.μ[i] = ops.μ[i] - dt * dotops.μ[i]
        oldops.λ[i] = ops.λ[i] - dt * dotops.λ[i]
    end
    oldops.CRR .= ops.CRR .- dt .* dotops.CRR
    oldops.CPR .= ops.CPR .- dt .* dotops.CPR
    oldops.CPP .= ops.CPP .- dt .* dotops.CPP
    oldops.p .= ops.p .- dt .* dotops.p
    oldops.q .= ops.q .- dt .* dotops.q
    oldops.ke = ops.ke - dt * dotops.ke

    return nothing
end

"""

Leapfrog step forward (CEID)

"""
function CEIDforward!(ops :: CEIDOps, oldops :: CEIDOps,
                    dotops :: CEIDOps, dt :: Float64)
    oldops.ρ = oldops.ρ + 2 * dt * dotops.ρ
    for i in 1:ops.No
        oldops.μ[i] = oldops.μ[i] + 2 * dt * dotops.μ[i]
        oldops.λ[i]= oldops.λ[i] + 2 * dt * dotops.λ[i]
    end
    oldops.CRR .= oldops.CRR .+ (2 * dt) .* dotops.CRR
    oldops.CPR .= oldops.CPR .+ (2 * dt) .* dotops.CPR
    oldops.CPP .= oldops.CPP .+ (2 * dt) .* dotops.CPP
    oldops.p .= oldops.p .+ (2 * dt) .* dotops.p
    oldops.q .= oldops.q .+ (2 * dt) .* dotops.q
    oldops.ke = oldops.ke + 2 * dt * dotops.ke

    return nothing
end

"""

Calculate CEID EOM RHS for K with no coupling between oscillators

"""
function CEIDcalcdots!(dotops :: CEIDOps, ops :: CEIDOps,
                    dt :: Float64, sbm :: SBModel; thermostat = false)

    # Time derivatives of generalized coordinates and momenta ------------------
    fbar = zeros(Float64,sbm.No)
    for ν in 1:ops.No
        fbar[ν] = real(tr(F(sbm,ops.q,ν)*ops.ρ)) -
                  real(tr(K(sbm,ops.q,ν,ν)*ops.μ[ν]))
        dotops.p[ν] = fbar[ν]
        dotops.q[ν] = ops.p[ν]
    end
    #---------------------------------------------------------------------------

    # Time derivative of ρ -----------------------------------------------------
    dotops.ρ = IOHBAR * comm(H(sbm,ops.q,ops.p),ops.ρ)
    for ν in 1:ops.No
        dotops.ρ += - IOHBAR * comm(F(sbm,ops.q,ν),ops.μ[ν]) +
                    0.5 * IOHBAR * ops.CRR[ν,ν] * comm(K(sbm,ops.q,ν,ν),ops.ρ)
    end
    #---------------------------------------------------------------------------

    # Time derivative of μ -----------------------------------------------------
    Threads.@threads for ν in 1:ops.No
        dotops.μ[ν] = ops.λ[ν] + IOHBAR * comm(H(sbm,ops.q,ops.p),ops.μ[ν])
        for νp in 1:ops.No
            dotops.μ[ν] += - IOHBAR * ops.CRR[ν,νp] * comm(F(sbm,ops.q,νp),ops.ρ)
        end
    # Time derivative of λ -----------------------------------------------------
        dotops.λ[ν] = IOHBAR * comm(H(sbm,ops.q,ops.p),ops.λ[ν]) +
                      0.5 * acomm(F(sbm,ops.q,ν),ops.ρ) - tr(F(sbm,ops.q,ν)*ops.ρ)*ops.ρ -
                      0.5 * acomm(K(sbm,ops.q,ν,ν),ops.μ[ν])
        for νp in 1:ops.No
            dotops.λ[ν] += - IOHBAR * ops.CPR[ν,νp] * comm(F(sbm,ops.q,νp),ops.ρ)
        end
    # Time derivatives of mean field second order correlators ------------------
        if !thermostat
            for νp in 1:ops.No
                dotops.CRR[ν,νp] = ops.CPR[ν,νp] + ops.CPR[νp,ν]
                dotops.CPR[ν,νp] = ops.CPP[ν,νp] + real(tr(F(sbm,ops.q,ν)*ops.μ[νp])) -
                                   ops.CRR[ν,νp] * real(tr(K(sbm,ops.q,ν,ν)*ops.ρ))
                dotops.CPP[ν,νp] = real(tr(F(sbm,ops.q,ν)*ops.λ[νp])) +
                                   real(tr(F(sbm,ops.q,νp)*ops.λ[ν])) -
                                   ops.CPR[ν,νp] * real(tr(K(sbm,ops.q,νp,νp)*ops.ρ)) -
                                   ops.CPR[νp,ν] * real(tr(K(sbm,ops.q,ν,ν)*ops.ρ))
            end
        else
            dotops.CRR .= 0.0
            dotops.CPR .= 0.0
            dotops.CPP .= 0.0
        end
    end
    #---------------------------------------------------------------------------

    # Time derivative of quantum kinetic energy --------------------------------
    dotops.ke = 0.0
    for ν in 1:ops.No
        for νp in 1:ops.No
            dotops.ke += real(tr(K(sbm,ops.q,ν,νp)*ops.ρ)) * ops.CPR[ν,νp] +
                            real(0.5 * IOHBAR * ops.CRR[νp,ν] *
                            tr(K(sbm,ops.q,ν,νp)*comm(H(sbm,ops.q,ops.p),ops.ρ)))
        end
    end
    #---------------------------------------------------------------------------

    return nothing
end


"""

Calculate CEID EOM RHS

"""
function CEIDcalcdots_full!(dotops :: CEIDOps, ops :: CEIDOps,
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
    dotops.ρ = IOHBAR * comm(H(sbm,ops.q,ops.p),ops.ρ)
    for ν in 1:ops.No
        dotops.ρ += - IOHBAR * comm(F(sbm,ops.q,ν),ops.μ[ν])
    end
    for ν in 1:ops.No
        for νp in 1:ops.No
            dotops.ρ += 0.5 * IOHBAR * ops.CRR[ν,νp] * comm(K(sbm,ops.q,ν,νp),ops.ρ)
        end
    end
    #---------------------------------------------------------------------------

    # Time derivative of μ -----------------------------------------------------
    Threads.@threads for ν in 1:ops.No
        dotops.μ[ν] = ops.λ[ν] + IOHBAR * comm(H(sbm,ops.q,ops.p),ops.μ[ν])
        for νp in 1:ops.No
            dotops.μ[ν] += - IOHBAR * ops.CRR[ν,νp] * comm(F(sbm,ops.q,νp),ops.ρ)
        end
    # Time derivative of λ -----------------------------------------------------
        dotops.λ[ν] = IOHBAR * comm(H(sbm,ops.q,ops.p),ops.λ[ν])
        dotops.λ[ν] += 0.5 * acomm(F(sbm,ops.q,ν),ops.ρ) - tr(F(sbm,ops.q,ν)*ops.ρ)*ops.ρ
        for νp in 1:ops.No
            dotops.λ[ν] += - IOHBAR * ops.CPR[ν,νp] * comm(F(sbm,ops.q,νp),ops.ρ)
            dotops.λ[ν] += - 0.5 * acomm(K(sbm,ops.q,ν,νp),ops.μ[νp])
        end
    # Time derivatives of mean field second order correlators ------------------
        if !thermostat
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
        else
            dotops.CRR .= 0.0
            dotops.CPR .= 0.0
            dotops.CPP .= 0.0
        end
    end
    #---------------------------------------------------------------------------

    # Time derivative of quantum kinetic energy --------------------------------
    dotops.ke = 0.0
    for ν in 1:ops.No
        for νp in 1:ops.No
            dotops.ke += real(tr(K(sbm,ops.q,ν,νp)*ops.ρ)) * ops.CPR[ν,νp] +
                            real(0.5 * IOHBAR * ops.CRR[νp,ν] *
                            tr(K(sbm,ops.q,ν,νp)*comm(H(sbm,ops.q,ops.p),ops.ρ)))
        end
    end
    #---------------------------------------------------------------------------

    return nothing
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
function RunCEID!(sbm,ops,nsteps,dt,store!,storage; thermostat = false, knondiag = false)

    dotops = CEIDOps(SA[0.0im 0.0im; 0.0im 0.0im],sbm)
    oldops = CEIDOps(SA[0.0im 0.0im; 0.0im 0.0im],sbm)

    # Bootstrap the integrator
    if !knondiag
        sb.CEIDcalcdots!(dotops, ops, dt, sbm; thermostat = thermostat)
    else
        sb.CEIDcalcdots_full!(dotops, ops, dt, sbm; thermostat = thermostat)
    end
    sb.CEIDbootstrap!(ops,oldops,dotops,dt)

    # call results storage function at t = 0
    t = 0.0
    store!(storage,t,ops,dotops,sbm)

    # integrate
    for i in 1:nsteps-1
        t = i*dt
        if !knondiag
            sb.CEIDcalcdots!(dotops, ops, dt, sbm; thermostat = thermostat)
        else
            sb.CEIDcalcdots_full!(dotops, ops, dt, sbm; thermostat = thermostat)
        end
        sb.CEIDforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        store!(storage,t,ops,dotops,sbm)
    end

    return nothing
end

#
