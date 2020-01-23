
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
    No :: Integer = 50
end

"""

structure container for an Ohmic Spin Boson Model

"""
struct SBModel
    No :: Integer
    Δ :: Float64
    ϵ :: Float64
    ωs :: Vector{Float64}
    cs :: Vector{Float64}
    function SBModel(sbp::OhmicSBParams)
        @unpack No, Δ, ϵ, ωm, ωc = sbp
        Δ = sbp.Δ
        ϵ = sbp.ϵ
        ω = zeros(Float64,No)
        c = zeros(Float64,No)
        for j=1:No
            ω[j] = ωm - ωc * log((j - exp(ωm/ωc)*j + exp(ωm/ωc)*No)/No)
            c[j] = sqrt((2/π)*ω[j]*ohmicJ0(ω[j],sbp)/ρosc(ω[j],sbp))
        end
        new(No,Δ,ϵ,ω,c)
    end
end

"""

Ohmic spectral density

"""
function ohmicJ0(ω, sbp :: OhmicSBParams)
    @unpack α, ωc = sbp
    (π/2)*α*ω.*exp.(-ω/ωc)
end

"""

Oscilator theoretical density of states

"""
function ρosc(ω, sbp :: OhmicSBParams)
    @unpack No, ωc, ωm = sbp
    (No/ωc)*exp.(-ω/ωc)/(1-exp.(-ωm/ωc))
end

"""

Calculate thermal ocupation numbers for the oscillators

"""
function thermalns(sbp :: OhmicSBParams, sbm :: SBModel)
    N = Vector{Float64}(undef,sbm.No)
    for i in 1:sbm.No
        N[i] = 1.0 / (exp(sbm.ωs[i] * sbp.β) - 1.0)
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
function densityofstates(values,range;σ=1.5)
    dos = zeros(size(range))
    for x in values
        dos += gaussian.(σ,range .- x)
    end
    return dos
end

"""

Calculate spectral density for a model

"""
function spectraldensity(cs,ωs,range;σ=1.5)
    sd = zeros(size(range))
    for (i,ω) in enumerate(ωs)
        sd += cs[i]^2*gaussian.(σ,range .- ω)/ωs[i]
    end
    return pi*sd/2
end

"""

Hamiltonian for an Ohmic Spin Boson Model

"""
function H(sbm :: SBModel, q :: Vector{Float64}, p :: Vector{Float64})
    H = sbm.ϵ * σz + sbm.Δ * σx
    for nu = 1:sbm.No
        H += eye * (p[nu]^2 + sbm.ωs[nu]^2 * q[nu]^2) / 2 + σz * sbm.cs[nu] * q[nu]
    end
    return H
end

"""

Force operator for an Ohmic Spin Boson Model

"""
function F(sbm :: SBModel, q :: Vector{Float64}, nu :: Integer)
    F = - eye * sbm.ωs[nu]^2 * q[nu] - σz * sbm.cs[nu]
    return F
end

"""

Second order derivative operator for an Ohmic Spin Boson Model

"""
function K(sbm :: SBModel, q :: Vector{Float64}, nu :: Integer, nup :: Integer)
    if nu == nup
        K = eye * sbm.ωs[nu]^2 + σz * sbm.cs[nu]
    else
        K = zm
    end
end
