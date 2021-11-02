import pytest
import numpy as np 
import scipy.optimize as opt
import matplotlib.pyplot as plt
def residue(X,N,r) : 
    x = X[:N]
    y = X[N:]

    dxdr = np.gradient(x,r,edge_order=2)

    res1 = dxdr - y
    res2 = y + x + 2

    return np.concatenate((res1,res2)) 



def test_krylov():
    N=1000
    r = np.linspace(1,5,N)
    sol_theory = 2*np.exp(-r) -2 

    X0 = np.concatenate((2*r,2*r))

    sol = opt.root(residue,X0,args=(N,r), method='krylov')

    plt.plot(r,sol.x[N:])
    plt.plot(r,sol_theory)
    plt.show()

    assert (sol.x[N:] == sol_theory).all() 



def test_derivative():
    N=1000
    r = np.logspace(1,5,N)

    f = r**0.5 + np.sin(r)

    dfdr_theory = 0.5*r**(-0.5) + np.cos(r)
    dfdr = np.gradient(f,r,edge_order=2)


    assert (dfdr_theory == dfdr).all()