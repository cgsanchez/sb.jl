using sb

sbp = OhmicSBParams(No = 5000, ωm = 300)
sbm = SBModel(sbp)

omegarange = range(0,length = 1001, stop=sbp.ωm*1.2)

specden = sb.spectraldensity(sbm,omegarange)
teospecden = sb.ohmicJ0(omegarange,sbp)

dos = sb.densityofstates(sbm,omegarange)
teodos = sb.ρosc(convert(Vector{Float64},omegarange),sbp)

@test abs(findmax(teospecden)[1] - findmax(specden)[1]) < 0.01
@test sum(abs.(specden[200:600]-teospecden[200:600])) < 0.1
@test sum(abs.(dos[200:600]-teodos[200:600])) < 2.0
