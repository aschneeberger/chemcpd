
import numpy as np
import matplotlib.pyplot as plt 

def test_func(V) :

    x = V[0]
    y = V[1]
    z = V[2]
    a = V[3]

    r0 = a + x - 3*y + z - 2
    r1 = -5*a + 3*x - 4*y + z 
    r2 = a + 2*y - z -1
    r3 = a + 2*x - 12

    return np.array([r0,r1,r2,r3])




def Jv(fun,U,V) :


    norm_U = np.sqrt( np.sum(U*U) )
    norm_V = np.sqrt( np.sum(V*V) )
    
    if norm_V == 0.0 :
        return  np.zeros(np.shape(V))

    eps = np.sqrt((1+norm_U)*np.finfo(float).eps) / norm_V 

    return (fun(U+eps*V) - fun(U)) / eps


def Arnoldi_basis_construct(func,U,v1,tol) :
    """
    Construct and Arnoldi vector basis in the Krylov space span(Jdu,J²du,J^3du...)
    return the Vector basis matrix and the 
    """
    j = 1
    H = np.array([[0]])#Hessenberg matrix = jacobian in krylov space H=Vk* J Vk
    Vk = np.array([v1]).T  #Arnoldi basis 
    res = 1
    while res > tol  :
        # Create a new Vk basis matrix 
        n_Vk = np.zeros((len(U),j+1))
        n_Vk[:,:j] = Vk
        Vk = n_Vk
        # Compute the new Vkj+1 guesse 
        # from the jacobian vec estimate
        Vk[:,j] = Jv(func,U,Vk[:,j-1])

        #Construction of the new Hessenberg matrix elements in the 
        # new column 
        for i in range(j):
            H[i,j-1] = np.dot(Vk[:,j],Vk[:,i])
            Vk[:,j] = Vk[:,j] - np.dot(H[i,j-1],Vk[:,i])
        
        # Add a new row and columns to the H matrix
        n_H = np.zeros((j+1,j+1))
        n_H[:j,:j] = H
        H = n_H
        H[j,j-1] = np.sqrt(np.sum(Vk[:,j] * Vk[:,j]))
        Vk[:,j] = Vk[:,j] / H[j,j-1]

        # The resolution is the last term of the H matrix
        res = H[j,j-1]
        j=j+1
    return Vk,H[:,:j-1]


def GMRES_naive(func,X0,tol) :
    """
    Naive GMRES algorithm based on the paper of Ayachour 2002. 
    It use an alternative minimization of the residue in Krylov space 
    without using Given rotation. This naive version is whithout optimization and 
    Imply two matrix inversion at each iteration
    """

    #Init the iteration
    X = X0
    r0_norm = 1.0

    while r0_norm  > tol :

        #Compute the current residue 
        r0 = Jv(test_func,X,np.zeros_like(X)) - test_func(X)
        r0_norm = np.sqrt(np.sum(r0*r0))
        # First vector of the Arnoldi basis 
        v1 = r0/r0_norm
        #Construct get the arnoldi basis and Hessenberg matrix 
        Vk,Hk_tilde = Arnoldi_basis_construct(func,X,v1,1e-30)

        # Construct the new du from the Ayachour paper 
        # Same notations are used, for more details see the paper 
        w = Hk_tilde[0,:]
        Hk = Hk_tilde[1:,:]

        # intermediate vectors construction
        u = np.dot(np.linalg.inv(Hk.T) ,w.T)

        # t_prime is the minimization of the residue in a var change framework
        c = 1/(1 + np.sum(u*u))
        t_prime = c*u

        # retreiving the z that minimize the real residue in krylov space 
        z = np.dot(np.linalg.inv(Hk),t_prime) 

        #Retreiving du from the krylov solution by basis change 
        du = np.dot(Vk[:,:-1],r0_norm*z)

        #Do the newton step
        X = X + du

        print("X",X)
        print("du",du)

    return X


X0 = np.array([10.0,1.0,1.0,1.0])

X = GMRES_naive(test_func,X0,1e-10) 