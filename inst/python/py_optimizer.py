import autograd.numpy as np

from pymanopt.manifolds import Stiefel
from pymanopt import Problem
from pymanopt.solvers import SteepestDescent

## compute negative log-likelihood
## if r > 0 and r + s < p (meaning we model with isotropic variation)
def Fpy1(V, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A):
    ## relevant variation
    G = V[:, 0:s]

    ## irrelevant variation
    G0 = V[:, s:(s+r)]
    G0part = (n + r + v0 + 1)/2 * np.linalg.slogdet(G0.T @ S @ G0 + U0)[1]
    sig2part = (n*(p-s-r)/2 + alpha) * np.log(np.trace(S)/2 - np.trace(V.T @ S @ V)/2 + k)
    ## Minimize the negative log likelihood
    return (n + s + v1 + 1 - q)/2 * np.linalg.slogdet(G.T @ A @ G + U1)[1] + G0part + sig2part


## compute negative log-likelihood
## if r > 0 and r + s = p (no isotropic variation)
def Fpy2(V, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A):

    ## relevant variation
    G = V[:, 0:s]

    ## irrelevant variation
    G0 = V[:, s:(s+r)]
    G0part = (n + r + v0 + 1)/2 * np.linalg.slogdet(G0.T @ S @ G0 + U0)[1]

    ## Minimize the negative log likelihood
    return (n + s + v1 + 1 - q)/2 * np.linalg.slogdet(G.T @ A @ G + U1)[1] + G0part


## compute negative log-likelihood
## if r = 0 (all isotropic variation)
def Fpy3(V, n, p, s, q, U1, alpha, v1, S, k, A):

    ## relevant variation
    G = V[:, 0:s]

    ## irrelevant variation
    sig2part = (n*(p-s)/2 + alpha) * np.log(np.trace(S)/2 - np.trace(V.T @ S @ V)/2 + k)

    
    ## Minimize the negative log likelihood
    return (n + s + v1 + 1 - q)/2 * np.linalg.slogdet(G.T @ A @ G + U1)[1] + sig2part



def optStiefel_py(Vinit, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A,
                  maxtime=1000, maxiters=1000):
    
    if (r > 0) & (r + s < p):
        def Fx(X):
            return Fpy1(X, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A)
    elif (r > 0) & (r + s == p):
        def Fx(X):
            return Fpy2(X, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A)
    else:
        def Fx(X):
            return Fpy3(X, n, p, s, q, U1, alpha, v1, S, k, A)
        
    manifold = Stiefel(p, (s+r))
    problem = Problem(manifold=manifold, cost=Fx)

    # Instantiate a Pymanopt solver
    solver = SteepestDescent(maxtime=maxtime, maxiter=maxiters)

    # let Pymanopt do the rest
    Xopt = solver.solve(problem, x=Vinit)
    return Xopt

def evalStiefel_py(X, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A):
    
    if (r > 0) & (r + s < p):
        return Fpy1(X, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A)
    elif (r > 0) & (r + s == p):
        return Fpy2(X, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A)
    else:
        return Fpy3(X, n, p, s, q, U1, alpha, v1, S, k, A)

