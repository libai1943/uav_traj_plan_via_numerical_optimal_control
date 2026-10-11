"""Seven-state coordinated-flight point-mass dynamics, SI units, angles in rad."""
import numpy as np


def motion(x, u, jacobian=False):
    x, u = np.asarray(x), np.asarray(u)
    V, ga, ch, mu = (x[..., i] for i in (3, 4, 5, 6))
    at, n, om = (u[..., i] for i in (0, 1, 2))
    g = 9.81
    f = np.stack((V*np.cos(ga)*np.cos(ch), V*np.cos(ga)*np.sin(ch), V*np.sin(ga),
                  at-g*np.sin(ga), g/V*(n*np.cos(mu)-np.cos(ga)),
                  g*n*np.sin(mu)/(V*np.cos(ga)), om), axis=-1)
    if not jacobian:
        return f
    A = np.zeros(x.shape[:-1]+(7, 7)); B = np.zeros(x.shape[:-1]+(7, 3))
    A[..., 0, 3] = np.cos(ga)*np.cos(ch); A[..., 0, 4] = -V*np.sin(ga)*np.cos(ch)
    A[..., 0, 5] = -V*np.cos(ga)*np.sin(ch)
    A[..., 1, 3] = np.cos(ga)*np.sin(ch); A[..., 1, 4] = -V*np.sin(ga)*np.sin(ch)
    A[..., 1, 5] = V*np.cos(ga)*np.cos(ch)
    A[..., 2, 3] = np.sin(ga); A[..., 2, 4] = V*np.cos(ga)
    A[..., 3, 4] = -g*np.cos(ga); B[..., 3, 0] = 1
    A[..., 4, 3] = -f[..., 4]/V; A[..., 4, 4] = g/V*np.sin(ga)
    A[..., 4, 6] = -g*n/V*np.sin(mu); B[..., 4, 1] = g/V*np.cos(mu)
    A[..., 5, 3] = -f[..., 5]/V; A[..., 5, 4] = f[..., 5]*np.tan(ga)
    A[..., 5, 6] = g*n*np.cos(mu)/(V*np.cos(ga)); B[..., 5, 1] = g*np.sin(mu)/(V*np.cos(ga))
    B[..., 6, 2] = 1
    return f, A, B


def rk4(x, u0, u1, h):
    um = (u0+u1)/2
    k1 = motion(x, u0); k2 = motion(x+h*k1/2, um)
    k3 = motion(x+h*k2/2, um); k4 = motion(x+h*k3, u1)
    return x+h*(k1+2*k2+2*k3+k4)/6
