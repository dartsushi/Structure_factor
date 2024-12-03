using MPSKit, TensorKit, KrylovKit
T = ComplexF64
Ising_Tc = 2.0/log(1.0+sqrt(2))
Potts_Tc = 1.0/log(1.0+sqrt(3))

function Ising_mpo(β::Float64)
    # MPSKit notaion: (left,top,bottom,right) 
    V = ℂ^2
    A_Ising = TensorMap(zeros,V⊗V←V⊗V)
    
    s(ind) = 2*ind-3
    ss(i,j) = s(i)*s(j)
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    E = -(ss(i,j)+ss(j,k)+ss(k,l)+ss(l,i))
                    A_Ising[i,j,l,k] = exp(-β*E)
                end
            end
        end
    end
    return A_Ising
end

# magnetic impurity
function Ising_impurity(β::Float64)
    # MPSKit notaion: (left,top,bottom,right) 
    V = ℂ^2
    A_Ising = TensorMap(zeros,V⊗V←V⊗V)
    
    s(ind) = 2*ind-3
    ss(i,j) = s(i)*s(j)
    for i=1:2
        for j=1:2
            for k=1:2
                for l=1:2
                    E = -(ss(i,j)+ss(j,k)+ss(k,l)+ss(l,i))
                    A_Ising[i,j,l,k] = (s(i)+s(j)+s(k)+s(l))/4*exp(-β*E)
                end
            end
        end
    end
    return A_Ising
end