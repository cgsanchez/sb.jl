using sb

@test sum(abs.(3*sb.σy() - [[0 -3im]; [3im 0]])) < eps(1.0)
@test sum(abs.(sb.zm())) < eps(1.0)
@test sum(abs.(3*sb.σz() - [[3 0]; [0 -3]])) < eps(1.0)
@test sum(abs.(3*sb.eye() - [[3 0]; [0 3]])) < eps(1.0)
@test sum(abs.(3*sb.σx() - [[0 3]; [3 0]])) < eps(1.0)
