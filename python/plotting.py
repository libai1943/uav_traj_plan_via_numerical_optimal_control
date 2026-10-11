"""Equal-scale 3-D figures; translucent columns and safe corridor fragments."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Circle, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

BLUE, ORANGE, TEAL = '#0066ad', '#ac681f', '#478e99'


def box_faces(box):
    vertices = np.array([[box[d+3] if (q >> d) & 1 else box[d] for d in range(3)] for q in range(8)])
    return [vertices[f] for f in ([0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4],
                                  [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5])]


def plot_result(result, cfg, folder):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'savefig.facecolor': 'white'})
    r = result; raw, seed, X, U = (r[k] for k in ('astar_path', 'spline_path', 'X', 'U'))
    legend = [Line2D([], [], color='.55', lw=1.5, label='3-D A* path'),
              Line2D([], [], color=ORANGE, lw=1.7, ls='--', label='Cubic-spline initial guess'),
              Line2D([], [], color=BLUE, lw=2.8, label='Optimized trajectory'),
              Patch(facecolor=TEAL, alpha=.22, label='Safe corridor fragments'),
              Patch(facecolor='.5', alpha=.24, label='Column obstacles')]
    fig = plt.figure(figsize=(13.5, 7.2)); ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
    theta, h = np.meshgrid(np.linspace(0, 2*np.pi, 80), [0, 200])
    for cx, cy, radius in r['columns']:
        ax.plot_surface(cx+radius*np.cos(theta), cy+radius*np.sin(theta), h,
                        color='.5', alpha=.22, linewidth=0, shade=False)
        for z in (0, 200):
            ax.plot(cx+radius*np.cos(theta[0]), cy+radius*np.sin(theta[0]), z,
                    color='.55', lw=.7, alpha=.6)
    for box in r['boxes']:
        ax.add_collection3d(Poly3DCollection(box_faces(box), facecolors=(.28, .56, .60, .025),
                                            edgecolors=(.28, .56, .60, .45), linewidths=.55))
    ax.plot(*raw.T, color='.55', lw=1.1, alpha=.8)
    ax.plot(*seed.T, color=ORANGE, ls='--', lw=1.8)
    ax.plot(*X[:, :3].T, color=BLUE, lw=3)
    ax.scatter(*X[0, :3], color=BLUE, s=45); ax.scatter(*X[-1, :3], color=BLUE, marker='*', s=110)
    ax.set(xlim=(-20, 640), ylim=(-220, 220), zlim=(0, 205), xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
    ax.set_box_aspect((660, 440, 205), zoom=1.12); ax.view_init(elev=28, azim=-63)
    ax.set_xticks([0, 200, 400, 600]); ax.set_yticks([-200, 0, 200]); ax.set_zticks([0, 100, 200])
    ax.xaxis.labelpad = 12; ax.yaxis.labelpad = 12; ax.zaxis.labelpad = 10
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis): axis.pane.fill = False
    ax.set_title(r['name'].replace('_', ' ').title()+' | fixed-wing trajectory and safe corridors', pad=15)
    fig.legend(handles=legend, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(.5, .025))
    fig.subplots_adjust(left=0, right=1, bottom=.12, top=.92)
    fig.savefig(folder/'trajectory_3d.png', dpi=240)

    fig2, ax2 = plt.subplots(figsize=(12, 6.6))
    for box in r['boxes']:
        ax2.add_patch(Rectangle(box[:2], box[3]-box[0], box[4]-box[1], facecolor=TEAL, alpha=.045, edgecolor=TEAL))
        ax2.add_patch(Rectangle(box[:2], box[3]-box[0], box[4]-box[1], fill=False, ec=TEAL, lw=.55, alpha=.45))
    for cx, cy, radius in r['columns']:
        ax2.add_patch(Circle((cx, cy), radius, facecolor='.82', edgecolor='.4'))
        ax2.add_patch(Circle((cx, cy), radius+cfg['safety_margin_m'], fill=False, ls=':', ec='.4'))
    ax2.plot(*raw[:, :2].T, color='.55', lw=1.3)
    ax2.plot(*seed[:, :2].T, '--', color=ORANGE, lw=1.7)
    ax2.plot(*X[:, :2].T, color=BLUE, lw=2.7)
    ax2.scatter(*X[0, :2], color=BLUE, s=40); ax2.scatter(*X[-1, :2], color=BLUE, marker='*', s=100)
    ax2.set(xlim=(-20, 640), ylim=(-220, 220), xlabel='X (m)', ylabel='Y (m)',
            title='Plan view | dotted circles include the 8 m safety margin')
    ax2.set_aspect('equal', adjustable='box'); ax2.grid(alpha=.15)
    fig2.legend(handles=legend, loc='lower center', ncol=3, frameon=False)
    fig2.tight_layout(rect=(0, .1, 1, 1)); fig2.savefig(folder/'trajectory_top.png', dpi=240)

    fig3, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    quantities = [X[:, 3], np.rad2deg(X[:, 4]), np.rad2deg(X[:, 6]), U[:, 0], U[:, 1], np.rad2deg(U[:, 2])]
    labels = ['Airspeed V (m/s)', 'Flight-path angle gamma (deg)', 'Bank angle mu (deg)',
              'Tangential input a_t (m/s^2)', 'Load factor n', 'Bank rate omega_mu (deg/s)']
    bounds = [(cfg['state_lower'][3], cfg['state_upper'][3]), np.rad2deg([cfg['state_lower'][4], cfg['state_upper'][4]]),
              np.rad2deg([cfg['state_lower'][6], cfg['state_upper'][6]]),
              (cfg['input_lower'][0], cfg['input_upper'][0]), (cfg['input_lower'][1], cfg['input_upper'][1]),
              np.rad2deg([cfg['input_lower'][2], cfg['input_upper'][2]])]
    for ax, q, label, bound in zip(axes.ravel(), quantities, labels, bounds):
        ax.plot(r['time'], q, color=BLUE, lw=2)
        for b in bound: ax.axhline(b, ls=':', color='.5', lw=1)
        ax.set_ylabel(label); ax.grid(alpha=.16)
    for ax in axes[-1]: ax.set_xlabel('Time (s)')
    fig3.suptitle('Optimized states and inputs | dotted lines: prescribed bounds')
    fig3.tight_layout(); fig3.savefig(folder/'states_inputs.png', dpi=200)

    fig4, axes = plt.subplots(1, 2, figsize=(11, 4))
    hist = r['history']
    axes[0].plot(np.arange(1, len(hist)+1), hist[:, 0], color=BLUE)
    axes[1].semilogy(np.arange(1, len(hist)+1), np.maximum(hist[:, 1], 1e-16), color=BLUE)
    for ax, label in zip(axes, ['Objective J', 'Maximum scaled constraint violation']):
        ax.set(xlabel='SQP iteration', ylabel=label); ax.grid(alpha=.2)
    fig4.tight_layout(); fig4.savefig(folder/'convergence.png', dpi=200)
