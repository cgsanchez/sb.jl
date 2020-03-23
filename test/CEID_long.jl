using sb, LinearAlgebra, DelimitedFiles, StaticArrays

struct storages
    ts :: Vector{Float64}
    sz :: Vector{Float64}
    function storages()
        new(Vector{Float64}(),Vector{Float64}())
    end
end

function store!(storage, t, ops, sbm)
    append!(storage.ts, t)
    append!(storage.sz, real(tr(sb.σz*ops.ρ)))
end

const NSTEPS = 2500
const DT = 0.010
const TOL = 1.0e-11

function my_run(sbm)
    storage = storages()
    ops = CEIDOps(SA[1.0 0.0; 0.0 0.0],sbm)
    RunCEID!(sbm,ops,NSTEPS,DT,store!, storage)
    return storage.sz
end

# Figure 1

sz = my_run(SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 7.5, α = 0.1, No = 60, β = Inf)))
refdata = readdlm("CEID_ref_results/fig1_CEID.dat")[:,2]
@test TOL > sum(abs.(refdata-sz))

# Figure 2

sz = my_run(SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 7.5, α = 0.6, No = 60, β = Inf)))
refdata = readdlm("CEID_ref_results/fig2_CEID.dat")[:,2]
@test TOL > sum(abs.(refdata-sz))

# Figure 3

sz = my_run(SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 2.5, α = 1.2, No = 60, β = Inf)))
refdata = readdlm("CEID_ref_results/fig3_CEID.dat")[:,2]
@test TOL > sum(abs.(refdata-sz))

# Figure 4

sz = my_run(SBModel(OhmicSBParams(ϵ = 1.0, Δ = 1.0, ωc = 7.5, α = 0.1, No = 60, β = Inf)))
refdata = readdlm("CEID_ref_results/fig4_CEID.dat")[:,2]
@test TOL > sum(abs.(refdata-sz))

# Figure 5

sz = my_run(SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 2.5, α = 0.09, No = 60, β = Inf)))
refdata = readdlm("CEID_ref_results/fig5_CEID.dat")[:,2]
@test TOL > sum(abs.(refdata-sz))

# Figure 6

sz = my_run(SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 2.5, α = 0.5625, No = 60, β = Inf)))
refdata = readdlm("CEID_ref_results/fig6_CEID.dat")[:,2]
@test TOL > sum(abs.(refdata-sz))

# Figure 7

sz = my_run(SBModel(OhmicSBParams(ϵ = 0.0, Δ = 1.0, ωc = 2.5, α = 0.09, No = 60, β = Inf)))
refdata = readdlm("CEID_ref_results/fig7_CEID.dat")[:,2]
@test TOL > sum(abs.(refdata-sz))

# 
