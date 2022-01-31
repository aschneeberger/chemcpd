
import numpy as np
import matplotlib.pyplot as plt 
from scipy.special import erf

def test_func(V) :

    x = V[0]
    y = V[1]
    z = V[2]
    a = V[3]

    r0 = a - 3*y + z 
    r1 = -5*a  - 4*y + z 
    r2 = a*x + 2*y*y**5 - z -1
    r3 = a+ 2*x - 12

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
    H = np.array([[0]]) #Hessenberg matrix = jacobian in krylov space H= Vk* J Vk
    Vk = np.array([v1]).T  #Arnoldi basis 
    res = 1


    while res > tol  :
        # Create a new Vk basis matrix 
        n_Vk = np.zeros((len(U),j+1))
        n_Vk[:,:j] = Vk
        Vk = n_Vk
        # Compute the new Vkj+1 guess 
        # from the jacobian vec estimate == J.Vk(j-1) == A.Vk(j-1) in literature 
        Vk[:,j] = Jv(func,U,Vk[:,j-1])

        #Construction of the new Hessenberg matrix elements in the 
        # new column. 
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
    without using Given rotation. This naive version is without optimization and 
    Imply two matrix inversion at each iteration
    """

    #Init the iteration
    X = X0
    r0_norm = 1.0
    it=0
    res = np.sqrt(np.sum(func(X0)**2))
    du = 1.00
    while res  > tol and it<it_max :

        du_old = du
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
        while ratio >= 1.0 and s != 0.0 :
            s =s/2
            ratio = np.sum(func(X+s*du)**2) / np.sum(func(X)**2) 
            print(s,ratio)

        if s == 0.0 :
            s= 0.5 
        X = X + s*du 

        res = np.sqrt(np.sum(func(X)**2))
        it+=1
        # plt.clf()
        # plt.cla()
        # plt.plot(r,X,'+')
        # plt.pause(0.1)
        print(it,res,np.sqrt(np.sum(du**2)),np.sqrt(np.sum((du/du_old)**2)))
    plt.close()
    return X


def GMRES_given_jfnk(func,u,du0,tol,max_iter):
    """
    Given rotation GMRES version, it include the computation of the Arnoldi basis and 
    the reduction of the Vk vector basis matrix. It return the du that minimize the 
    residual (||fu - Jdu||₂) for a system of dimension N
    
    -------
    inputs:
    -------
    
    func : function : the function we want to minimize 
    u : array(N) : Current guess vector u from the JFNK step we are seaching the du for. 
    du0 : array(N) : the du step guess, in JFNK it is the null vector 
    tol : float : the desired tolerance in the residual.
    max_iter : int : if iterations stagnates and the tolerance is not reached, this is the maximum iteration performed

    -------
    return:
    -------

    du : array(N) : The step that minimize the residual ||fu - Jdu||₂ for the current  JNFK step 
    """

    # Results of the function at u 
    fu_init = func(u)
    
    #Initialize the residual 
    res = fu_init - Jv(func,u,du0) 
    res_norm = np.sqrt(np.sum(res*res))

    #Arnoldi basis 
    Vk = [0] * max_iter 
    Vk[0] = res/res_norm

    #Hessenberg matrix
    H = np.zeros((max_iter+1,max_iter))
    
    #Given rotation coefficients 
    Sn = np.zeros((max_iter,1))
    Cs = np.zeros((max_iter, 1))


    #Preparing the vector fu for given rotation 
    fu = np.zeros(max_iter+1)
    fu[0] = res_norm

    for k in range(max_iter+1): 
        
        # Estimation of a Krylov vector 
        Vk_estimate = Jv(func,u,np.array(Vk[k]))
        
        for j in range(k+1) :
            # Orthogonalisation of the vector 
            H[j,k] = np.dot(Vk[j],Vk_estimate)
            Vk_estimate = Vk_estimate - H[j,k] * Vk[j]
        
        #The next element of the Hessenberg matrix 
        Vk_estimate_norm = np.sqrt(np.sum(Vk_estimate*Vk_estimate))
        H[k+1,k] = Vk_estimate_norm 

        #We create the new vector if the value is not 0 
        if H[k+1,k] != 0 and k!= max_iter -1 :
            Vk[k+1] = Vk_estimate/H[k+1,k]
        else :
            break

        #Since H[k+1,k] is the equivalent to the residual, it is our break condition
        print(H[k+1,k])
        if np.abs(H[k+1,k]) < tol : 
            break 

        # We compute the Hessenberg matrix after k Given Rotations
        for i in range(k) : 
            # Since we need the value of H[i,k] to compute the value
            #  of the QR of  H[i+1,k]
            qr_temp = Cs[i] * H[i,k] + Sn[i] * H[i+1, k ]
            H[i+1,k] = -1.0 * Sn[i] * H[i,k] + Cs[i] * H[i+1,k]
            H[i,k] = qr_temp

        
        
        # Computation of the next Given Rotation factors 

        t = np.sqrt(H[k,k]**2 + H[k+1,k]**2)

        Cs[k] = H[k,k] / t 
        Sn[k] = H[k+1,k] / t 

        #Do the k+1 th Given Rotation transforming the Hessenberg Matrix into
        # pure upper Triangular matrix (easy to invert)
        H[k,k] = Cs[k] * H[k,k] + Sn[k] * H[k+1,k]
        
        # Since the transformation give an upper triangular matrix, this term is 
        # nullified by the k+1 th Given rotation
        H[k+1,k] = 0

        # We apply do the k+1 th rotation to fu to stay in the same basis
        fu[k+1] = -1.0 * Sn[k] * fu[k]
        fu[k] = Cs[k] * fu[k]
    # Here is the BackSubstitution of H by fu to get the 
    # the current guess solution in the Krylov space Vk
    # Will be done explicitly in a next iteration. 
    lmbd = np.linalg.lstsq(H,fu)[0]

    print(np.array(Vk))
    return np.dot(np.asarray(Vk).transpose(),lmbd)




def gmres_algorithm_givens (A , b , x0 , error , max_iter ):
    """
    Given rotation algorithm from github. Used to understand the implementation to write 
    An equivalent for the JFNK 
    """
    res = b - np.asarray(np.dot(A,x0)).reshape(-1) # residual error

    x_pred = []

    q_ = [0] * max_iter

    x_pred.append(res)

    q_[0] = res / np.linalg.norm(res)

    h_ = np.zeros((max_iter + 1, max_iter))

    sn = np.zeros((max_iter , 1))

    cs = np.zeros((max_iter , 1))

    b_ = np.zeros(max_iter + 1)

    b_[0] = np.linalg.norm(res)


    
    for k in range(min(max_iter , A.shape[0])) : 


        y_out = np.asarray(np.dot(A,q_[k])).reshape(-1)



        for j in range(k+1) : 
            h_[j , k] = np.dot(q_[j],y_out)
            y_out = y_out - h_[j , k] * q_[j]
        
        h_[k+1 , k] = np.linalg.norm(y_out) 

        if (h_[k + 1, k] != 0 and k != max_iter - 1):
            q_[k+1] = y_out / h_[k+1 , k]
        
        for i in range(k): 
            temp   =  cs[i] * h_[i , k] + sn[i] * h_[i+1  , k]
            h_[i+1 ,k] = -1*sn[i] * h_[i , k] + cs[i] * h_[i+1 , k]
            h_[i , k]   = temp
        
        t = np.sqrt(h_[k , k]**2 + h_[k+1 , k]**2)

        cs[k] = (h_[k , k]) /t
        sn[k] = (h_[k+1 , k]) /t

        h_[k , k] = cs[k] * h_[k , k] + sn[k] * h_[k + 1 , k]
        h_[k + 1 , k] = 0

        b_[k + 1] = -1 * sn[k] * b_[k]
        b_[k] = cs[k] * b_[k]

        #print(h_)
        #print(b_)

        c_ = np.linalg.lstsq(h_ , b_)[0] 

        prod_ = np.asarray(np.dot(np.asarray(q_).transpose() , c_))

        if (k == max_iter - 1) :
            print('q_ ' + str(np.asarray(q_).shape) + ' c_shape = ' + str(c_.shape) + ' prod_ = ' + str(prod_.shape))

        #print(prod_)

        x_pred.append(prod_ + x0)  

        x_temp_ = (np.linalg.norm(b - np.dot(A ,(prod_ + x0)).reshape(-1)) / np.linalg.norm(b))


        print(x_temp_)

        if (x_temp_ < error) :
            print("exit val : ", k)
            break

    return x_pred

N=1000
dr = 0.001
r = np.arange(dr,1,dr)

X0 = np.ones_like(r)*1e5

du0 = np.array([0,0,0,0])
u = np.array([0,0,0,0])
Uf = GMRES_restart_naive(test_func,u,1e-10,1000)

while True :
    du = GMRES_given_jfnk(test_func,u,du0,1e-30,30)
    print(du)
    u = u + du 
    plt.pause(1)
# X = GMRES_restart_naive(test_heat_eq,X0,1e-6,5e3) 
# plt.figure()
# plt.plot(r,X)
# plt.show()