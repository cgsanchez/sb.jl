using sb

sbp = OhmicSBParams(No = 5000, ωm = 300)

sbm = SBModel(sbp)

omegarange = range(0,sbp.ωm*1.2,length=1001)

dos = sb.densityofstates(sbm.ωs,omegarange)

teodos = sb.ρosc(convert(Vector{Float64},omegarange),sbp)

using Plots

plot(omegarange,dos)

findmax(dos)

omegarange[19]

plot!(omegarange,teodos)

specden = sb.spectraldensity(sbm.cs,sbm.ωs,omegarange)
teospecden = sb.ohmicJ0(omegarange,sbp)

plot(omegarange,specden)
plot!(omegarange,teospecden)

findmax(teospecden)[1]
findmax(specden)

omegarange[105]

sum(abs.(specden[200:600]-teospecden[200:600]))
sum(abs.(dos[200:600]-teodos[200:600]))

plot(omegarange[200:600],specden[200:600]-teospecden[200:600])
plot!(omegarange[200:600],dos[200:600]-teodos[200:600])

# This is testing to construct Ehrenfest

Pkg.activate(".")

using sb, LinearAlgebra

sbp = OhmicSBParams(No = 100)

sbm = SBModel(sbp)

ops = EhrenfestOps([[1.0 0.0]; [0.0 0.0]],sbm)
dotops = EhrenfestOps(sb.zm,sbm)
oldops = EhrenfestOps(sb.zm,sbm)

dt = 0.001

sb.ehcalcdots!(dotops, ops ,dt , sbm)

sb.ehbootstrap!(ops,oldops,dotops,dt)


function integrate(nsteps,dt)
    t = 0.0
    ts = Vector{Float64}(undef,nsteps)
    evals = Vector{Float64}(undef,nsteps)
    ts[1] = t
    evals[1] = real(tr(sb.σz * ops.ρ))
    for i in 1:nsteps-1
        t = i*dt
        println(i)
        sb.ehcalcdots!(dotops,ops,dt,sbm)
        sb.ehforward!(ops,oldops,dotops,dt)
        (ops,oldops) = (oldops,ops)
        println(ops)
        ts[i+1] = t
        evals[i+1] = real(tr(sb.σz * ops.ρ))
    end
    return ts,evals
end


# Strange stuff

function this()
    a = 1
    b = 2
    for i in 1:5
        println(i," ",a,)
        (a,b) = (b,a)
    end
end

this()

for i in 1:5
    println(a)
    a = 5
end
