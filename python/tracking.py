"""Lecture supplement: reference-linearized, receding-horizon tracking MPC.

Input bounds only; no obstacle/state constraints or estimator are designed here.
Uses the full nonlinear plant and an initial position/heading error, no wind.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.linalg import expm
from .model import motion, rk4


def box_qp(H, f, lo, hi, guess):
    """Strictly convex box QP via damped Newton on a logarithmic barrier."""
    x = np.clip(guess, lo+.001*(hi-lo), hi-.001*(hi-lo))
    for mu in (1., .1, .01, .001, 1e-4, 1e-5, 1e-6, 1e-7):
        for _ in range(70):
            dl, du = x-lo, hi-x
            grad = H@x+f-mu/dl+mu/du
            step = np.linalg.solve(H+np.diag(mu/dl**2+mu/du**2), -grad)
            if -grad@step < 1e-12: break
            alpha = 1.
            if np.any(step < 0): alpha = min(alpha, .99*np.min(-dl[step < 0]/step[step < 0]))
            if np.any(step > 0): alpha = min(alpha, .99*np.min(du[step > 0]/step[step > 0]))
            old = .5*x@H@x+f@x-mu*np.sum(np.log(dl)+np.log(du))
            while alpha > 1e-12:
                xx = x+alpha*step
                value = .5*xx@H@xx+f@xx-mu*np.sum(np.log(xx-lo)+np.log(hi-xx))
                if value <= old+1e-4*alpha*(grad@step): break
                alpha /= 2
            x += alpha*step
    residual = float(np.max(np.abs(x-np.clip(x-(H@x+f), lo, hi))))
    if residual > 2e-3: raise RuntimeError(f'MPC QP projected-gradient residual: {residual}')
    return x, residual


def track(result, cfg):
    h, horizon = .2, 15
    ts = np.arange(round(cfg['duration_s']/h)+1)*h
    ref = PchipInterpolator(result['time'], result['X'])(ts)
    uu = np.column_stack([np.interp(ts, result['time'], result['U'][:, d]) for d in range(3)])
    # A fixed wing continues flying beyond the destination instead of hovering.
    for j in range(1, horizon+1):
        ref = np.vstack((ref, ref[len(ts)-1]+[20*j*h, 0, 0, 0, 0, 0, 0]))
        uu = np.vstack((uu, [0., 1., 0.]))
    actual = ref[0]+[0, 5, 3, 0, 0, np.pi/60, 0]
    opened = actual.copy(); guess = np.zeros(3*horizon)
    rows, applied, residuals = [], [], []
    # Reference linearizations can be precomputed; feedback changes each QP's linear term.
    mats = []
    for i in range(len(ref)-1):
        _, A, B = motion(ref[i], uu[i], True)
        aug = np.zeros((10, 10)); aug[:7, :7], aug[:7, 7:] = A, B
        disc = expm(aug*h)
        _, Av, _ = motion(ref[i+1], uu[i], True)
        C = np.vstack((np.c_[np.eye(3), np.zeros((3, 4))], Av[:3]))
        Q = C.T@np.diag([1/25]*3+[1/4]*3)@C+np.diag([0, 0, 0, .01, .03, .03, .03])
        mats.append((disc[:7, :7], disc[:7, 7:], Q, ref[i+1]-rk4(ref[i], uu[i], uu[i], h)))
    for k in range(len(ts)-1):
        P = np.eye(7); G = np.zeros((7, 3*horizon)); offset = np.zeros(7)
        H = np.zeros((3*horizon, 3*horizon)); f = np.zeros(3*horizon)
        e = ref[k]-actual
        for j in range(horizon):
            A, B, Q, defect = mats[k+j]
            P, G, offset = A@P, A@G, A@offset+defect
            G[:, 3*j:3*j+3] -= B
            H += 2*G.T@Q@G; f += 2*G.T@Q@(P@e+offset)
        H += 2*np.kron(np.eye(horizon), np.diag([.02, .15, .05]))+1e-9*np.eye(3*horizon)
        lo = (np.array(cfg['input_lower'])-uu[k:k+horizon]).ravel()
        hi = (np.array(cfg['input_upper'])-uu[k:k+horizon]).ravel()
        du, residual = box_qp((H+H.T)/2, f, lo, hi, guess)
        guess = np.r_[du[3:], [0., 0., 0.]]  # Shift the previous solution: warm start.
        control = uu[k]+du[:3]
        rows.append([ts[k], *actual, *opened, *ref[k]])
        applied.append(control); residuals.append(residual)
        for _ in range(5):
            actual = rk4(actual, control, control, h/5)
            opened = rk4(opened, uu[k], uu[k], h/5)
    rows.append([ts[-1], *actual, *opened, *ref[len(ts)-1]])
    rows = np.array(rows)
    mpc_error = np.linalg.norm(rows[:, 1:4]-rows[:, 15:18], axis=1)
    open_error = np.linalg.norm(rows[:, 8:11]-rows[:, 15:18], axis=1)
    min_clear = min(float(np.min(np.linalg.norm(rows[:, 1:3]-o[:2], axis=1)-o[2])) for o in result['columns'])
    return dict(time=ts, rows=rows, inputs=np.array(applied), qp_residuals=np.array(residuals),
                sample_time_s=h, prediction_steps=horizon, mpc_position_error=mpc_error,
                open_loop_position_error=open_error, min_sampled_mpc_clearance_m=min_clear,
                rms_mpc_position_error_m=float(np.sqrt(np.mean(mpc_error**2))),
                rms_open_loop_position_error_m=float(np.sqrt(np.mean(open_error**2))))


def plot_tracking(result, out):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    rows = result['rows']; t = result['time']
    axes[0].plot(rows[:, 15], rows[:, 16], '--', color='.45', label='Planned reference')
    axes[0].plot(rows[:, 8], rows[:, 9], color='#ac681f', label='Open loop')
    axes[0].plot(rows[:, 1], rows[:, 2], color='#0066ad', label='MPC feedback')
    axes[0].set(xlabel='X (m)', ylabel='Y (m)', title='Same initial position and heading error')
    axes[0].axis('equal'); axes[0].legend(frameon=False)
    axes[1].plot(t, result['open_loop_position_error'], color='#ac681f', label='Open loop')
    axes[1].plot(t, result['mpc_position_error'], color='#0066ad', label='MPC feedback')
    axes[1].set(xlabel='Time (s)', ylabel='Position error norm (m)', title='Position / velocity tracking MPC')
    axes[1].legend(frameon=False)
    for ax in axes: ax.grid(alpha=.15)
    fig.tight_layout(); fig.savefig(out/'mpc_tracking.png', dpi=240)
