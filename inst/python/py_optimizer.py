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
    G0part = (n + v0)/2 * np.linalg.slogdet(G0.T @ S @ G0 + U0)[1]
    sig2part = (n*(p-s-r)/2 + alpha) * np.log(np.trace(S)/2 - np.trace(V.T @ S @ V)/2 + k)
    ## Minimize the negative log likelihood
    return (n + v1 - q)/2 * np.linalg.slogdet(G.T @ A @ G + U1)[1] + G0part + sig2part


## compute negative log-likelihood
## if r > 0 and r + s = p (no isotropic variation)
def Fpy2(V, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A):

    ## relevant variation
    G = V[:, 0:s]

    ## irrelevant variation
    G0 = V[:, s:(s+r)]
    G0part = (n + v0)/2 * np.linalg.slogdet(G0.T @ S @ G0 + U0)[1]

    ## Minimize the negative log likelihood
    return (n + v1 - q)/2 * np.linalg.slogdet(G.T @ A @ G + U1)[1] + G0part


## compute negative log-likelihood
## if r = 0 (all isotropic variation)
def Fpy3(V, n, p, s, q, U1, alpha, v1, S, k, A):

    ## relevant variation
    G = V[:, 0:s]

    ## irrelevant variation
    sig2part = (n*(p-s)/2 + alpha) * np.log(np.trace(S)/2 - np.trace(V.T @ S @ V)/2 + k)

    
    ## Minimize the negative log likelihood
    return (n + v1 - q)/2 * np.linalg.slogdet(G.T @ A @ G + U1)[1] + sig2part



def create_envelope_problem_py(Vinit, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A, L=0):

    
    if (r > 0) & (r + s < p):
        def Fx(X):
            return Fpy1(X, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A)
    elif (r > 0) & (r + s == p):
        def Fx(X):
            return Fpy2(X, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A)
    else:
        def Fx(X):
            return Fpy3(X, n, p, s, q, U1, alpha, v1, S, k, A)

    def Gx(X):
        if L > 0:
            return Fx(X) + L*np.sum(np.sqrt(np.diag(X[:, 0:s] @ X[:, 0:s].T)))
        else:
            return Fx(X)
            
    manifold = Stiefel(p, (s+r))
    problem = Problem(manifold=manifold, cost=Gx)
    return problem

def optim_envelope_py(Vinit, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A, L=0, maxtime=1000, maxiters=1000):

    problem = create_envelope_problem_py(Vinit, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A, L)
    
    # Instantiate a Pymanopt solver
    solver = SteepestDescent(maxtime=maxtime, maxiter=maxiters)

    # let Pymanopt do the rest
    Xopt = solver.solve(problem, x=Vinit)
    return Xopt

def envelope_cost_py(V, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A, L=0):
    problem = create_envelope_problem_py(V, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A, L)
    return problem.cost(V)

def envelope_grad_py(V, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A, L=0):
    problem = create_envelope_problem_py(V, n, p, s, r, q, U0, U1, alpha, v0, v1, S, k, A, L)
    return problem.grad(V)
