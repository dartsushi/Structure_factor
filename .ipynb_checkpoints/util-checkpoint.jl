using MPSKit, TensorKit, KrylovKit
T = ComplexF64
Ising_Tc = 2.0/log(1.0+sqrt(2))
Potts_Tc = 1.0/log(1.0+sqrt(3))

function initialize_Ising(β::Float64)
    V = ℂ^2
    A_Ising = TensorMap(zeros,V⊗V←V⊗V)
    
    s(ind) = 2*ind-3
    ss(i,j) = s(i)*s(j)
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    E = -(ss(i,j)+ss(j,k)+ss(k,l)+ss(l,i))
                    A_Ising[i,j,k,l] = exp(-β*E)
                end
            end
        end
    end
    return A_Ising
end