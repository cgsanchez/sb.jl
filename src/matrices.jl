
# define Pauli matrices needed for Hamiltonian and other operators

function σz()
    return [[1.0+0.0im 0.0+0.0im]; [0.0+0.0im -1.0+0.0im]]
end

function σx()
    return [[0.0+0.0im 1.0+0.0im]; [1.0+0.0im  0.0+0.0im]]
end

function σy()
    return [[0.0+0.0im 0.0-1.0im]; [0.0+1.0im  0.0+0.0im]]
end

function zm()
    return [[0.0+0.0im 0.0+0.0im]; [0.0+0.0im  0.0+0.0im]]
end

function eye()
     return [[1.0+0.0im 0.0+0.0im]; [0.0+0.0im  1.0+0.0im]]
end

function comm(A,B)
    A*B-B*A
end

function acomm(A,B)
    A*B+B*A
end
