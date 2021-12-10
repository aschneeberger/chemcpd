
import numpy as np
import matplotlib.pyplot as plt 

def test_func(V) :

    x = V[0]
    y = V[1]
    z = V[2]
    a = V[3]

    r0 = a + x - 3*y + z - 2
    r1 = -5*a + 3*x - 4*y + z 
    r2 = a + 2*y - z 
    r3 = a + 2*x - 12

    return np.array([r0,r1,r2,r3])




def Jv(fun,U,V) :

    norm_U = np.sqrt( np.sum(U*U) )
    norm_V = np.sqrt( np.sum(V*V) )
    
    if norm_V == 0.0 :
        return  np.zeros(np.shape(V))

    eps = np.sqrt((1+norm_U)*np.finfo(float).eps) / norm_V 

    return (fun(U+eps*V) - fun(U)) / eps


def Arnoldi_basis_construct(fun,U,v1,tol) :
    j = 1
    H = np.array([[0]])#Hessenberg matrix = jacobian in krylov space H=Vk* J Vk
    Vk = np.array([v1]).T  #Arnoldi basis 
    res = 1
    while res > tol  :
        n_Vk = np.zeros((len(U),j+1))
        n_Vk[:,:j] = Vk
        Vk = n_Vk
        Vk[:,j] = Jv(fun,U,Vk[:,j-1])

        for i in range(j):
            H[i,j-1] = np.dot(Vk[:,j],Vk[:,i])
            Vk[:,j] = Vk[:,j] - np.dot(H[i,j-1],Vk[:,i])
        
        n_H = np.zeros((j+1,j+1))
        n_H[:j,:j] = H
        H = n_H
        H[j,j-1] = np.sqrt(np.sum(Vk[:,j] * Vk[:,j]))
        Vk[:,j] = Vk[:,j] / H[j,j-1]
        res = H[j,j-1]
        j=j+1
    return Vk,H

U0 = np.array([1.0,1.0,1.0,1.0])

v1 = Jv(test_func,U0,0.0) - test_func(U0)
v1 = v1/np.sqrt(np.sum(v1*v1))
Vk = Arnoldi_basis_construct(test_func,U0,v1,1e-70)
print(Vk)
plt.imshow(np.log10(np.abs(Vk[1])))
plt.colorbar()
plt.show()
plt.close()