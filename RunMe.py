"""One entry point: python RunMe.py (add --no-show for a headless run)."""
import argparse
import json
import os
from pathlib import Path
import sys

# Avoid excessive BLAS threads in the small dense SQP subproblems.
# Explicit user settings always take precedence.
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ.setdefault(key, '1')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', default='all', choices=['all', 'two_columns', 'four_columns', 'five_columns'])
    parser.add_argument('--no-show', action='store_true')
    parser.add_argument('--skip-mpc', action='store_true', help='Skip the two-column tracking supplement')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        import numpy as np
        import scipy
        import matplotlib
        if args.no_show: matplotlib.use('Agg')
        from python.planner import PlanningProblem, check_derivatives
        from python.plotting import plot_result
    except ImportError as exc:
        sys.exit(f'Missing dependency: {exc}\nInstall first: python -m pip install -r requirements.txt')
    root = Path(__file__).resolve().parent
    cfg = json.loads((root/'scenes.json').read_text(encoding='utf-8'))
    out = args.output or root/'outputs'/'python'
    out.mkdir(parents=True, exist_ok=True)
    results = []
    for case in cfg['cases']:
        if args.case != 'all' and args.case != case['name']: continue
        print(f"\n=== {case['name']}: 3-D A* -> cubic spline -> corridors -> NLP ===", flush=True)
        problem = PlanningProblem(cfg, case)
        err = check_derivatives(problem)
        print(f'A*: {len(problem.raw)} points; {problem.expanded} expanded nodes; derivative error {err:.3g}', flush=True)
        result = problem.solve()
        result['derivative_check_max_error'] = err
        folder = out/case['name']; folder.mkdir(exist_ok=True)
        serial = lambda obj: obj.tolist() if isinstance(obj, np.ndarray) else obj
        (folder/'result.json').write_text(json.dumps(result, default=serial, indent=2), encoding='utf-8')
        np.savez_compressed(folder/'result.npz', **{k: v for k, v in result.items() if isinstance(v, np.ndarray)})
        plot_result(result, cfg, folder)
        if case['name'] == 'two_columns' and not args.skip_mpc:
            from python.tracking import track, plot_tracking
            print('Tracking supplement: local linear model -> QP -> first input -> feedback', flush=True)
            tracking = track(result, cfg)
            (folder/'mpc_tracking.json').write_text(json.dumps(tracking, default=serial, indent=2), encoding='utf-8')
            plot_tracking(tracking, folder)
            print(f"MPC/open-loop RMS position error: {tracking['rms_mpc_position_error_m']:.3f} / {tracking['rms_open_loop_position_error_m']:.3f} m", flush=True)
        print(json.dumps(result['checks'], indent=2), flush=True)
        results.append({k: result[k] for k in ('name', 'objective', 'iterations', 'seconds', 'checks')})
    report = dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                  matplotlib=matplotlib.__version__, config=cfg, results=results)
    (out/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(f'\nAll selected cases passed. Figures/data saved to {out}', flush=True)
    if not args.no_show:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == '__main__':
    main()
