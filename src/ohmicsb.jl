
using Parameters, LinearAlgebra

"""

Parameter container for an Ohmic Spin Boson Model

"""
# NOTE: Change the default type tp be Real to allow for simple precision eventually
@with_kw struct OhmicSBParams @deftype Float64
    Δ = 5.0
    ϵ = 0.0
    α = 0.1
    ωc = 7.5 * Δ
    ωm = 5 * ωc
    β = 1.0
    No :: Int64 = 50
end

"""

structure container for an Ohmic Spin Boson Model

"""
struct SBModel
    No :: Int64
    Δ :: Float64
    ϵ :: Float64
    β :: Float64
    ωs :: Vector{Float64}
    cs :: Vector{Float64}
    function SBModel(sbp::OhmicSBParams)
        @unpack No, Δ, ϵ, ωm, ωc = sbp
        Δ = sbp.Δ
        ϵ = sbp.ϵ
        β = sbp.β
        ωs = zeros(Float64,No)
        cs = zeros(Float64,No)
        for j=1:No
            ωs[j] = ωm - ωc * log((j - exp(ωm/ωc)*j + exp(ωm/ωc)*No)/No)
            cs[j] = sqrt((2/π)*ωs[j]*ohmicJ0(ωs[j],sbp)/ρosc(ωs[j],sbp))
        end
        new(No,Δ,ϵ,β,ωs,cs)
    end
end

"""

Ohmic spectral density

"""
function ohmicJ0(ω, sbp)
    @unpack α, ωc = sbp
    (π/2)*α*ω.*exp.(-ω/ωc)
end

"""

Oscilator theoretical density of states

"""
function ρosc(ω, sbp)
    @unpack No, ωc, ωm = sbp
    (No/ωc)*exp.(-ω/ωc)/(1-exp.(-ωm/ωc))
end

"""

Calculate thermal ocupation numbers for the oscillators

"""
function thermalns(sbm)
    N = Vector{Float64}(undef,sbm.No)
    for i in 1:sbm.No
        N[i] = 1.0 / (exp(HBAR * sbm.ωs[i] * sbm.β) - 1.0)
    end
    return N
end

"""

Gaussian distribution

"""
function gaussian(σ,x)
    exp(-((x/σ)^2)/2)/(σ*sqrt(2*π))
end

"""

Calculate the oscillator density of states for a model

"""
function densityofstates(sbm,range;σ=1.5)
    dos = zeros(size(range))
    for x in sbm.ωs
        dos += gaussian.(σ,range .- x)
    end
    return dos
end

"""

Calculate spectral density for a model

"""
function spectraldensity(sbm,range;σ=1.5)
    sd = zeros(size(range))
    for (i,ω) in enumerate(sbm.ωs)
        sd += sbm.cs[i]^2*gaussian.(σ,range .- ω)/sbm.ωs[i]
    end
    return pi*sd/2
end

"""

Hamiltonian for an Ohmic Spin Boson Model

"""
@inline function H(sbm, q , p)
    H = sbm.ϵ * σz + sbm.Δ * σx
    for nu = 1:sbm.No
        H += eye * 0.5 * sbm.ωs[nu]^2 * q[nu]^2 + σz * sbm.cs[nu] * q[nu]
    end
    return H
end

"""

Force operator for an Ohmic Spin Boson Model

"""
@inline function F(sbm, q, nu)
    F = - eye * sbm.ωs[nu]^2 * q[nu] - σz * sbm.cs[nu]
    return F
end

"""

Second order derivative operator for an Ohmic Spin Boson Model (scalar)

"""
@inline function K(sbm, q, nu, nup)
    if nu == nup
        K = sbm.ωs[nu]^2
    else
        K = 0.0
    end
end
