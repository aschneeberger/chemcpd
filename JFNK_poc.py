
import numpy as np
import matplotlib.pyplot as plt 
from scipy.special import erf

def test_func(V) :

    x = V[0]
    y = V[1]
    z = V[2]
    a = V[3]

    r0 = a + np.cos(x*y + np.sin(a)) - 3*y + z - 20*a/y
    r1 = -5*a  - 4*y + z 
    r2 = a*x + 2*y*y**5 - z -1
    r3 = erf(a) + 2*x - 12

    return np.array([r0,r1,r2,r3])

def test_exp_diff(V) :

    global dr
    W= np.zeros_like(V)
    dV = np.gradient(V,r,edge_order=2) 
    dV[0] = (1-V[1])/(dr*2)
    dV[-1] = (np.exp(10) - V[-2])/(dr*2)
    
    W = dV - 10*V
    return W

def test_heat_eq(V):
    W= np.zeros_like(V)
    dV = np.gradient(V,r,edge_order=2)
    ddV = np.gradient(dV,r,edge_order=2)

    V0 = 800
    Vn = 200

    dV[0] = (V[1] - V0)/(2*dr) 
    dV[-1] = (Vn - V[-2])/(2*dr) 

    ddV[0] = (V0 + V[1] - 2*V[0]) / (dr**2)
    ddV[-1] = (Vn + V[-2] - 2*V[-1]) / (dr**2)
    

    W = ddV + 1e4 * np.exp(-(r-0.5)**2/0.02) + 1e4*np.exp(-(r-0.1)**2/0.02) 
    return W

def Jv(fun,U,V) :
    """
    Compute an approximation of the dot product 
    between the system's jacobian at U and a vector V
    It is approximated as (F[U+epsV] - F[U]) / eps
    """

    #Computation of optimum epsilon from vec norms and 
    # Machine precision
    norm_U = np.sqrt( np.sum(U*U) )
    norm_V = np.sqrt( np.sum(V*V) )
    
    #if norm_V is zeros then the dot product us null
    if norm_V == 0.0 :
        return  np.zeros(np.shape(V))

    # Computation of the optimal eps from Knoll 2003
    eps = np.sqrt((1+norm_U)*np.finfo(float).eps) / norm_V 

    #Computation of Jv
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

        #Upadate alpha 
        # gammak = 
        # sintetak = H[j,j-1] * gammak
        # alpha = alpha *sintetak
        j=j+1
            
    return Vk,H[:,:j-1]


def GMRES_restart_naive(func,X0,tol,it_max) :
    """
    Naive GMRES algorithm based on the paper of Ayachour 2002. 
    It use an alternative minimization of the residue in Krylov space 
    without using Given rotation. This naive version is whithout optimization and 
    Imply two matrix inversion at each iteration
    """

    #Init the iteration
    X = X0
    r0_norm = 1.0
    it=0
    res = np.sqrt(np.sum(func(X0)**2))

    while res  > tol and it<it_max :

        #Compute the current residue 
        r0 = Jv(func,X,np.zeros_like(X)) - func(X)
        r0_norm = np.sqrt(np.sum(r0*r0))
        # First vector of the Arnoldi basis 
        v1 = r0/r0_norm
        #Construct get the arnoldi basis and Hessenberg matrix 
        Vk,Hk_tilde = Arnoldi_basis_construct(func,X,v1,1e-5)

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

        #Do the newton step by line search
        s = 1.0
        ratio = np.sum(func(X+s*du)**2) / np.sum(func(X)**2) 
        while ratio >= 1.0 :
            s =s/2
            ratio = np.sum(func(X+s*du)**2) / np.sum(func(X)**2) 
            print(s,ratio)

        X = X + s*du 

        res = np.sqrt(np.sum(func(X)**2))
        it+=1
        plt.clf()
        plt.cla()
        plt.plot(r,X,'+')
        plt.pause(0.1)
        print(it,res)
    plt.close()
    return X





N=1000
dr = 0.001
r = np.arange(dr,1,dr)

X0 = np.ones_like(r)*1e5

X = GMRES_restart_naive(test_heat_eq,X0,1e-6,5e3) 
plt.figure()
plt.plot(r,X)
plt.show()