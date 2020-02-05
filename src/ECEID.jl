using sb, LinearAlgebra, StaticArrays

"""

Structure to hold data for ECEID dynamics

"""
mutable struct ECEIDOps
    No :: Int64
    ρ  :: SArray{Tuple{2,2},Complex{Float64},2,4}
    Cs :: Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}
    Cc :: Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}
    As :: Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}
    Ac :: Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}
    N  :: Vector{Float64}
    p  :: Vector{Float64}
    q  :: Vector{Float64}
    function ECEIDOps(ρ0, sbm :: SBModel)
        Cs = Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}(undef,sbm.No)
        Cc = Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}(undef,sbm.No)
        As = Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}(undef,sbm.No)
        Ac = Vector{SArray{Tuple{2,2},Complex{Float64},2,4}}(undef,sbm.No)
        N = sb.thermalns(sbm)
        p = zeros(Float64,sbm.No)
        q = zeros(Float64,sbm.No)
        for i in 1:sbm.No
            Cs[i] = SA[0.0im 0.0im; 0.0im 0.0im]
            Cc[i] = SA[0.0im 0.0im; 0.0im 0.0im]
            Ac[i] = SA[0.0im 0.0im; 0.0im 0.0im]
            As[i] = SA[0.0im 0.0im; 0.0im 0.0im]
        end
        new(sbm.No,ρ0,Cs,Cc,As,Ac,N,p,q)
    end
end

"""

Leapfrog bootstrap of data (ECEID)

"""
function ECEIDbootstrap!(ops :: ECEIDOps, oldops :: ECEIDOps,
                      dotops :: ECEIDOps, dt :: Float64)
    oldops.ρ = ops.ρ - dt * dotops.ρ
    for i in 1:ops.No
        oldops.Cs[i] = ops.Cs[i] - dt * dotops.Cs[i]
        oldops.Cc[i] = ops.Cc[i] - dt * dotops.Cc[i]
        oldops.As[i] = ops.As[i] - dt * dotops.As[i]
        oldops.Ac[i] = ops.Ac[i] - dt * dotops.Ac[i]
        oldops.N[i]  = ops.N[i] - dt * dotops.N[i]
    end
end

"""

Leapfrog step forward (ECEID)

"""
function ECEIDforward!(ops :: ECEIDOps, oldops :: ECEIDOps,
                    dotops :: ECEIDOps, dt :: Float64)
    oldops.ρ = oldops.ρ + 2 * dt * dotops.ρ
    Threads.@threads for i in 1:ops.No
        oldops.Cs[i] = oldops.Cs[i] + 2 * dt * dotops.Cs[i]
        oldops.Cc[i] = oldops.Cc[i] + 2 * dt * dotops.Cc[i]
        oldops.As[i] = oldops.As[i] + 2 * dt * dotops.As[i]
        oldops.Ac[i] = oldops.Ac[i] + 2 * dt * dotops.Ac[i]
        oldops.N[i]  = oldops.N[i] + 2 * dt * dotops.N[i]
    end
end

"""

Leapfrog ECEID calculate EOM RHS

"""
function ECEIDcalcdots!(dotops :: ECEIDOps, ops :: ECEIDOps,
                    dt :: Float64, sbm :: SBModel; thermostat = false)

    dotops.ρ = - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.ρ)
    for i in 1:ops.No
        μ = (1.0im*ops.Cc[i] - ops.As[i]) / sbm.ωs[i]
        dotops.ρ += IOHBAR * comm(F(sbm,ops.q,i),μ)
    end
    for i in 1:ops.No
        dotops.Cc[i] = - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.Cc[i]) + sbm.ωs[i] * ops.Cs[i] + (ops.N[i] + 0.5) * comm(F(sbm,ops.q,i),ops.ρ)
        dotops.Cs[i] = - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.Cs[i]) - sbm.ωs[i] * ops.Cc[i]
        dotops.Ac[i] = - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.Ac[i]) + sbm.ωs[i] * ops.As[i] + 0.5 * acomm(F(sbm,ops.q,i),ops.ρ)
        dotops.As[i] = - IOHBAR * comm(H(sbm,ops.q,ops.p),ops.As[i]) - sbm.ωs[i] * ops.Ac[i]
        λ = (1.0im*ops.Cs[i] + ops.Ac[i])
        if !thermostat
            dotops.N[i] = real(tr(F(sbm,ops.q,i) * λ)) / HBAR / sbm.ωs[i]
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

"""

This function integrates the ECEID EOMs starting from the
initial condition provided in ops for nsteps of dt.
The store! callback function is called at every timestep.

"""
function RunECEID!(sbm,ops,nsteps,dt,store!,storage; thermostat = false)

    dotops = ECEIDOps(SA[0.0im 0.0im; 0.0im 0.0im],sbm)
    oldops = ECEIDOps(SA[0.0im 0.0im; 0.0im 0.0im],sbm)

    # Bootstrap the integrator
    sb.ECEIDcalcdots!(dotops, ops , dt , sbm; thermostat = thermostat)
    sb.ECEIDbootstrap!(ops,oldops,dotops,dt)

    # call results storage function at t = 0
    t = 0.0
    store!(storage,t,ops,sbm)

    # integrate
    for i in 1:nsteps-1
        t = i*dt
        sb.ECEIDcalcdots!(dotops,ops,dt,sbm; thermostat = thermostat)
        sb.ECEIDforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        store!(storage,t,ops,sbm)
    end

    return nothing
end
