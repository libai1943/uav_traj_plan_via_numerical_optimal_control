"""3-D A*, checked cubic-spline initialization and axis-aligned safe corridors.

Columns span every admissible altitude. Each row is (center X, center Y, radius).
The search really has 26 neighbors in XYZ; it is NOT a fixed-wing kinodynamic A*.
"""
import heapq
import itertools
import numpy as np
from scipy.interpolate import CubicSpline


def segment_clearance(a, b, columns):
    """Exact minimum distance of an XY segment to the physical column surfaces."""
    v = np.asarray(b)[:2] - np.asarray(a)[:2]
    w = columns[:, :2] - np.asarray(a)[:2]
    t = np.clip(w @ v / max(v @ v, 1e-30), 0, 1)
    return float(np.min(np.linalg.norm(w - t[:, None] * v, axis=1) - columns[:, 2]))


def box_clearance(lo, hi, columns):
    """Exact minimum distance from a whole axis-aligned box to column surfaces."""
    d = np.maximum(np.maximum(lo[:2] - columns[:, :2], columns[:, :2] - hi[:2]), 0)
    return float(np.min(np.linalg.norm(d, axis=1) - columns[:, 2]))


def astar_3d(start, goal, columns, cfg):
    step = np.array(cfg['grid_step_m'], float)
    lo, hi = np.array(cfg['search_lower_m']), np.array(cfg['search_upper_m'])
    dims = np.rint((hi - lo) / step).astype(int) + 1
    s, g = tuple(np.rint((start-lo)/step).astype(int)), tuple(np.rint((goal-lo)/step).astype(int))
    if not (np.allclose(lo+step*s, start) and np.allclose(lo+step*g, goal)):
        raise ValueError('Start and goal must lie on the configured 3-D grid.')
    margin = cfg['safety_margin_m'] + cfg['astar_extra_clearance_m']
    if segment_clearance(start, start, columns) <= margin or segment_clearance(goal, goal, columns) <= margin:
        raise ValueError('Start/goal violates A* clearance.')
    # Identical ordering in MATLAB; deterministic ties by insertion order.
    neighbors = list(itertools.product((-1, 0, 1), repeat=3))
    neighbors.remove((0, 0, 0))
    costs, parent, closed = {s: 0.0}, {}, set()
    serial, expanded = 0, 0
    heap = [(float(np.linalg.norm(goal-start)), serial, s)]
    while heap:
        _, _, node = heapq.heappop(heap)
        if node in closed:
            continue
        if node == g:
            nodes = [g]
            while nodes[-1] != s:
                nodes.append(parent[nodes[-1]])
            path = lo + np.array(nodes[::-1]) * step
            return path.astype(float), expanded
        closed.add(node)
        expanded += 1
        a = lo + np.array(node)*step
        for delta in neighbors:
            nxt = tuple(np.array(node) + delta)
            if nxt in closed or any(nxt[d] < 0 or nxt[d] >= dims[d] for d in range(3)):
                continue
            b = lo + np.array(nxt)*step
            # Edge check, not only endpoint occupancy: diagonal corner-cutting is disallowed.
            if segment_clearance(a, b, columns) <= margin:
                continue
            candidate = costs[node] + float(np.linalg.norm(b-a))
            if candidate < costs.get(nxt, np.inf)-1e-10:
                costs[nxt], parent[nxt] = candidate, node
                serial += 1
                heapq.heappush(heap, (candidate+float(np.linalg.norm(goal-b)), serial, nxt))
    raise RuntimeError('3-D A* found no path. Check grid, bounds and inflation.')


def spline_seed(raw, columns, cfg):
    """Visibility shortcut -> clamped cubic spline -> arc-length resampling.

    A conservative derivative bound certifies each sampled spline interval:
    chord clearance minus a second-derivative interpolation error bound.
    """
    margin = cfg['safety_margin_m']
    knots, i = [raw[0]], 0
    while i < len(raw)-1:
        j = len(raw)-1
        while j > i+1 and segment_clearance(raw[i], raw[j], columns) <= margin+cfg['astar_extra_clearance_m']:
            j -= 1
        knots.append(raw[j]); i = j
    knots = np.array(knots)
    # Try the shortcut first, then increasingly many original A* knots if needed.
    candidates = [knots, raw[::2] if (len(raw)-1) % 2 == 0 else np.vstack((raw[::2], raw[-1])), raw]
    for knots in candidates:
        arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(knots, axis=0), axis=1))]
        curve = CubicSpline(arc, knots, bc_type=((1, [1., 0., 0.]), (1, [1., 0., 0.])))
        # Include every polynomial break; one interval never crosses two cubics.
        samples = np.unique(np.r_[np.linspace(0, arc[-1], 1601), arc])
        dense = curve(samples)
        second = np.maximum(np.linalg.norm(curve(samples[:-1], 2), axis=1),
                            np.linalg.norm(curve(samples[1:], 2), axis=1))
        errors = second * np.diff(samples)**2 / 8
        lower, upper = np.array(cfg['state_lower'][:3]), np.array(cfg['state_upper'][:3])
        clear = min(segment_clearance(a, b, columns)-err
                    for a, b, err in zip(dense[:-1], dense[1:], errors))
        if clear > margin+2 and np.all(dense >= lower+errors.max()) and np.all(dense <= upper-errors.max()):
            length = np.r_[0, np.cumsum(np.linalg.norm(np.diff(dense, axis=0), axis=1))]
            query = np.interp(np.linspace(0, length[-1], cfg['intervals']+1), length, samples)
            return curve(query), dense, knots, float(clear)
    raise RuntimeError('The cubic spline is not safely clear. Increase A* clearance or refine the grid.')


def build_corridors(seed, columns, cfg):
    """Two collocation intervals per box; shared nodes lie in both neighboring boxes."""
    N = cfg['intervals']
    if N % 2:
        raise ValueError('intervals must be even (two intervals per fragment).')
    lower, upper = np.array(cfg['state_lower'][:3]), np.array(cfg['state_upper'][:3])
    boxes, owner = [], np.repeat(np.arange(N//2), 2)
    node_lo = np.tile(lower, (N+1, 1)); node_hi = np.tile(upper, (N+1, 1))
    margin = cfg['safety_margin_m']
    for i in range(N//2):
        ids = slice(2*i, 2*i+3)
        lo = np.maximum(seed[ids].min(axis=0)-[2, 2, 8], lower)
        hi = np.minimum(seed[ids].max(axis=0)+[2, 2, 8], upper)
        if box_clearance(lo, hi, columns) <= margin:
            raise RuntimeError('Seed box intersects an inflated column; refine intervals or increase seed clearance.')
        base_lo, base_hi = lo.copy(), hi.copy()
        for _ in range(26):
            for face in range(6):
                d, cl, ch = face//2, lo.copy(), hi.copy()
                if face % 2 == 0: cl[d] -= 2
                else: ch[d] += 2
                if cl[d] < base_lo[d]-[34, 40, 18][d] or ch[d] > base_hi[d]+[34, 40, 18][d]:
                    continue
                cl, ch = np.maximum(cl, lower), np.minimum(ch, upper)
                if box_clearance(cl, ch, columns) > margin+.05:
                    lo, hi = cl, ch
        boxes.append(np.r_[lo, hi])
        # 1 m interior buffer, additionally checked by dense RK4 integration.
        node_lo[ids] = np.maximum(node_lo[ids], lo+1)
        node_hi[ids] = np.minimum(node_hi[ids], hi-1)
    if np.any(node_lo > node_hi):
        raise RuntimeError('Consecutive corridor fragments have an empty node intersection.')
    return np.array(boxes), owner, node_lo, node_hi
