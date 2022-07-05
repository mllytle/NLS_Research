import numpy as np
import matplotlib.pyplot as plt

A = 1
GAM = -1
MU = -(A**2)/2

def jac(v, stepsize):
    dim = len(v)
    J = np.zeros((dim, dim))
    
    for ii in range(dim):
        J[ii][ii] = -2/(stepsize**2) - 6*GAM*v[ii]**2 + 2*MU

    for ii in range(dim-1):
        J[ii+1][ii] = 1/(stepsize**2)
        J[ii][ii+1] = 1/(stepsize**2)

    return J

def f(v, boundaries, stepsize):
    out = np.zeros(len(v))
    out[0] = (v[1] - 2*v[0] + boundaries[0])/(stepsize**2) - 2*GAM*v[0]**3 + 2*MU*v[0]
    out[-1] = (boundaries[1] - 2*v[-1] + v[-2])/(stepsize**2) - 2*GAM*v[-1]**3 + 2*MU*v[-1]
    for ii in range(1,len(v)-1):
        out[ii] = (v[ii+1] - 2*v[ii] + v[ii-1])/(stepsize**2) - 2*GAM*v[ii]**3 + 2*MU*v[ii]

    return out

def sol(bounds, steps):
    x = np.linspace(bounds[0], bounds[1], steps+1)
    return A*1/np.cosh(A*x)

def newtons(v0, niters, boundary_conditions, stepsize):
    v = v0
    residual = np.zeros(niters)
    for ii in range(niters):
        v = v - 0.5*np.linalg.inv(jac(v,stepsize))@f(v, boundary_conditions, stepsize)
        residual[ii] = np.linalg.norm(f(v, boundary_conditions, stepsize))
    return v, residual

m = 1000
bounds = np.array([-20,20])
h = (bounds[1] - bounds[0])/(m+1)
print(h)
bc = np.array([0,0])
iterations = 100

s = sol(bounds, m)
x = f(s, bc, h)
print(np.linalg.norm(x))
print(x)

y0 = s*2
y, res = newtons(y0, iterations, bc, h)

plt.figure()
plt.plot(s)
plt.plot(y)

plt.figure()
plt.plot(res)
plt.show()


