"""Scaled direct transcription (trapezoidal collocation) solved with SLSQP."""
import time
import numpy as np
from scipy.optimize import minimize, Bounds
from .model import motion, rk4
from .geometry import astar_3d, spline_seed, build_corridors, box_clearance


class PlanningProblem:
    def __init__(self, cfg, case):
        self.cfg, self.case = cfg, case
        self.N, self.dt = cfg['intervals'], cfg['duration_s']/cfg['intervals']
        self.M, self.obs = self.N+1, np.array(case['columns'], float)
        self.sx, self.su = np.array([600, 150, 150, 20, .3, .6, .6]), np.array([2, 1, .4])
        self.scale = np.r_[np.tile(self.sx, self.M), np.tile(self.su, self.M)]
        start, goal = np.array(cfg['start_state']), np.array(case['goal_state'])
        self.raw, self.expanded = astar_3d(start[:3], goal[:3], self.obs, cfg)
        seed, self.spline, self.knots, self.spline_clearance = spline_seed(self.raw, self.obs, cfg)
        self.boxes, self.owner, plo, phi = build_corridors(seed, self.obs, cfg)
        X = np.zeros((self.M, 7)); X[:, :3] = seed
        dp = np.gradient(seed, self.dt, axis=0)
        X[:, 3] = np.linalg.norm(dp, axis=1)
        X[:, 4] = np.arctan2(dp[:, 2], np.linalg.norm(dp[:, :2], axis=1))
        X[:, 5] = np.unwrap(np.arctan2(dp[:, 1], dp[:, 0]))
        X[:, 6] = np.arctan(X[:, 3]/9.81*np.gradient(X[:, 5], self.dt))
        X[0], X[-1] = start, goal
        U = np.column_stack((np.gradient(X[:, 3], self.dt)+9.81*np.sin(X[:, 4]),
                             np.cos(X[:, 4])/np.cos(X[:, 6]), np.gradient(X[:, 6], self.dt)))
        lx, hx = np.tile(cfg['state_lower'], (self.M, 1)), np.tile(cfg['state_upper'], (self.M, 1))
        lx[:, :3], hx[:, :3] = plo, phi
        if np.any(start < lx[0]) or np.any(start > hx[0]) or np.any(goal < lx[-1]) or np.any(goal > hx[-1]):
            raise RuntimeError('Endpoint lies outside its corridor or physical bounds.')
        lx[0] = hx[0] = start; lx[-1] = hx[-1] = goal
        self.lb = np.r_[lx.ravel(), np.tile(cfg['input_lower'], self.M)]/self.scale
        self.ub = np.r_[hx.ravel(), np.tile(cfg['input_upper'], self.M)]/self.scale
        self.z0 = np.clip(np.r_[X.ravel(), U.ravel()]/self.scale, self.lb, self.ub)
        self.W = np.ones(self.M); self.W[[0, -1]] = .5
        self._cache_z, self._cache = None, None

    def unpack(self, z):
        w = z*self.scale
        return w[:7*self.M].reshape(self.M, 7), w[7*self.M:].reshape(self.M, 3)

    def objective(self, z):
        X, U = self.unpack(z); weights = np.array([.15, 2, .8])
        du = U-[0, 1, 0]
        J = self.dt*np.sum(weights*(du/self.su)**2*self.W[:, None])
        J += .02*self.dt*np.sum((X[:, 1]/100)**2*self.W)
        gx = np.zeros_like(X); gx[:, 1] = .04*self.dt*X[:, 1]/1e4*self.W
        gu = 2*self.dt*weights*du/self.su**2*self.W[:, None]
        return float(J), np.r_[gx.ravel(), gu.ravel()]*self.scale

    def constraints(self, z):
        if self._cache_z is not None and np.array_equal(z, self._cache_z):
            return self._cache
        X, U = self.unpack(z); F, A, B = motion(X, U, True)
        eq = ((X[1:]-X[:-1]-self.dt/2*(F[:-1]+F[1:]))/self.sx).ravel()
        # SciPy ineq convention is >= 0 (opposite to MATLAB fmincon).
        P = X[None, :, :2]-self.obs[:, None, :2]
        rr = self.obs[:, 2]+self.cfg['safety_margin_m']
        ine = (np.sum(P**2, axis=2)/rr[:, None]**2-1).ravel()
        ji = np.zeros((len(ine), len(z))); je = np.zeros((len(eq), len(z)))
        for b in range(len(self.obs)):
            for k in range(self.M):
                ji[b*self.M+k, 7*k:7*k+2] = 2*P[b, k]*self.sx[:2]/rr[b]**2
        for k in range(self.N):
            rows, ix, iu = slice(7*k, 7*k+7), 7*k, 7*self.M+3*k
            je[rows, ix:ix+7] = (-np.eye(7)-self.dt/2*A[k])*self.sx/self.sx[:, None]
            je[rows, ix+7:ix+14] = (np.eye(7)-self.dt/2*A[k+1])*self.sx/self.sx[:, None]
            je[rows, iu:iu+3] = -self.dt/2*B[k]*self.su/self.sx[:, None]
            je[rows, iu+3:iu+6] = -self.dt/2*B[k+1]*self.su/self.sx[:, None]
        self._cache_z, self._cache = z.copy(), (eq, ine, je, ji)
        return self._cache

    def solve(self):
        history = []
        def record(z):
            eq, ine, _, _ = self.constraints(z)
            history.append([self.objective(z)[0], max(np.max(np.abs(eq)), np.max(-ine), 0)])
        begin = time.perf_counter()
        result = minimize(self.objective, self.z0, jac=True, method='SLSQP',
                          bounds=Bounds(self.lb, self.ub),
                          constraints=[{'type': 'eq', 'fun': lambda z: self.constraints(z)[0],
                                        'jac': lambda z: self.constraints(z)[2]},
                                       {'type': 'ineq', 'fun': lambda z: self.constraints(z)[1],
                                        'jac': lambda z: self.constraints(z)[3]}],
                          callback=record, options={'maxiter': 350, 'ftol': 1e-8, 'disp': True})
        elapsed = time.perf_counter()-begin
        if not result.success:
            raise RuntimeError(f"{self.case['name']}: SLSQP failed ({result.message}); no successful result exported.")
        X, U = self.unpack(result.x); X0, U0 = self.unpack(self.z0)
        checks, dense, rollout = self.validate(result.x)
        return dict(name=self.case['name'], X=X, U=U, initial_X=X0, initial_U=U0,
                    astar_path=self.raw, spline_path=self.spline, spline_knots=self.knots,
                    boxes=self.boxes, interval_box=self.owner, columns=self.obs,
                    time=np.arange(self.M)*self.dt, dense=dense, open_loop_rollout=rollout,
                    objective=float(result.fun), iterations=int(result.nit), seconds=elapsed,
                    solver_message=str(result.message), astar_expanded_nodes=self.expanded,
                    history=np.array(history), checks=checks)

    def validate(self, z):
        X, U = self.unpack(z); eq, ine, _, _ = self.constraints(z)
        residual = max(0., float(np.max(np.abs(eq))), float(np.max(-ine)),
                       float(np.max(self.lb-z)), float(np.max(z-self.ub)))
        dense, rollout = [], [[0., *X[0]]]
        defect, violation, forward = 0., 0., X[0].copy()
        for k in range(self.N):
            local = X[k].copy(); box = self.boxes[self.owner[k]]
            # Dense checks reset at EACH collocation node. Full rollout is reported separately.
            for j in range(21):
                dense.append([k*self.dt+j*self.dt/20, *local])
                violation = max(violation, float(np.max(box[:3]-local[:3])), float(np.max(local[:3]-box[3:])))
                if j == 20: break
                u0 = U[k]+j/20*(U[k+1]-U[k]); u1 = U[k]+(j+1)/20*(U[k+1]-U[k])
                local = rk4(local, u0, u1, self.dt/20)
                forward = rk4(forward, u0, u1, self.dt/20)
                rollout.append([k*self.dt+(j+1)*self.dt/20, *forward])
            defect = max(defect, float(np.linalg.norm(local[:3]-X[k+1, :3])))
        dense, rollout = np.array(dense), np.array(rollout)
        clearance = min(float(np.min(np.linalg.norm(dense[:, 1:3]-o[:2], axis=1)-o[2])) for o in self.obs)
        global_clearance = min(float(np.min(np.linalg.norm(rollout[:, 1:3]-o[:2], axis=1)-o[2])) for o in self.obs)
        whole_box_clearance = min(box_clearance(b[:3], b[3:], self.obs) for b in self.boxes)
        checks = dict(max_scaled_constraint_violation=residual, min_dense_clearance_m=clearance,
                      min_whole_box_clearance_m=whole_box_clearance,
                      max_dense_corridor_violation_m=violation, max_interval_position_defect_m=defect,
                      spline_certified_clearance_m=self.spline_clearance,
                      open_loop_terminal_position_error_m=float(np.linalg.norm(forward[:3]-X[-1, :3])),
                      open_loop_min_clearance_m=global_clearance)
        if residual > 1e-5 or violation > 1e-5 or clearance < self.cfg['safety_margin_m'] or whole_box_clearance <= self.cfg['safety_margin_m'] or defect > .5:
            raise RuntimeError(f'Post-solve validation failed: {checks}. Refine mesh or initialization.')
        return checks, dense, rollout


def check_derivatives(problem):
    z = problem.z0.copy(); _, g = problem.objective(z)
    _, _, je, ji = problem.constraints(z)
    error = 0.
    for j in np.linspace(7, len(z)-4, 18).astype(int):
        dz = np.zeros_like(z); dz[j] = 1e-6
        ep, ip, _, _ = problem.constraints(z+dz); em, im, _, _ = problem.constraints(z-dz)
        error = max(error, abs((problem.objective(z+dz)[0]-problem.objective(z-dz)[0])/2e-6-g[j]),
                    np.max(np.abs((ep-em)/2e-6-je[:, j])), np.max(np.abs((ip-im)/2e-6-ji[:, j])))
    if error > 1e-5: raise AssertionError(f'Analytic derivatives failed: {error}')
    return float(error)
