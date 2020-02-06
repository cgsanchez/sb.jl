using StaticArrays

const σz = SA[1.0+0.0im 0.0+0.0im; 0.0+0.0im -1.0+0.0im]
const σx = SA[0.0+0.0im 1.0+0.0im; 1.0+0.0im  0.0+0.0im]
const σy = SA[0.0+0.0im 0.0-1.0im; 0.0+1.0im  0.0+0.0im]
const zm = SA[0.0+0.0im 0.0+0.0im; 0.0+0.0im  0.0+0.0im]
const eye = SA[1.0+0.0im 0.0+0.0im; 0.0+0.0im  1.0+0.0im]

@inline function comm(A,B)
    A*B - B*A
end

@inline function acomm(A,B)
    A*B + B*A
end
