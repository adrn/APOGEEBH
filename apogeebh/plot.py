# Third-party
from astropy.utils.misc import isiterable
import astropy.units as u
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np

from twoface.mass import get_m2_min, mf
from twoface.plot import plot_data_orbits, plot_phase_fold_residual
from twoface.samples_analysis import MAP_sample

__all__ = ['plot_diag_panel']

def plot_diag_panel(data, samples, joker_params, M1, data_src=None):
    # Compute things we'll need in plots:
    mfs = mf(samples['P'], samples['K'], samples['e'])
    if not isiterable(M1):
        M1 = [M1.value] * len(mfs) * M1.unit
    M2_mins = get_m2_min(M1, mfs)

    # Now make the plot!
    fig = plt.figure(figsize=(9, 10))
    gs = GridSpec(3, 5)

    ax1 = fig.add_subplot(gs[0, :])
    colormap = dict(DR15='k', LAMOST='tab:orange', Keck='tab:green')
    plot_data_orbits_kw = dict()
    plot_data_orbits_kw.setdefault('xlim_choice', 'tight')
    plot_data_orbits_kw.setdefault('highlight_P_extrema', False)
    plot_data_orbits(data, samples, ax=ax1, n_times=16384,
                     **plot_data_orbits_kw)

    if data_src is not None:
        for src in np.unique(data_src):
            sub_data = data[data_src == src]
            sub_data.plot(ax=ax1, markerfacecolor='none', markeredgewidth=1,
                          markeredgecolor=colormap[src], label=src)
        ax1.legend(loc='best')

    # ---
    sample = MAP_sample(data, samples, joker_params)
    axes2 = [fig.add_subplot(gs[1, 0:3]),
             fig.add_subplot(gs[2, 0:3])]
    plot_phase_fold_residual(data, sample, axes=axes2)
    # plt.setp(axes2[0].get_xticklabels(), visible=False)

    # ---
    ax3 = fig.add_subplot(gs[1, 3:])
    ax3.scatter(samples['P'].to(u.day).value, samples['e'],
                marker='o', linewidth=0, alpha=0.5, s=10)
    ax3.set_xlim(0.8, 2500)
    ax3.set_xscale('log')
    ax3.set_ylim(0, 1)
    ax3.set_ylabel('$e$')
    # plt.setp(ax3.get_xticklabels(), visible=False)

    # ---
    ax4 = fig.add_subplot(gs[2, 3:])
    ax4.scatter(samples['P'].to(u.day).value[:len(M2_mins)], M2_mins.value,
                marker='o', linewidth=0, alpha=0.5, s=10)
    ax4.set_xlim(0.8, 1024)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_ylim(1E-2, 1E3)
    ax4.axhline(1.4, zorder=-10, color='#cccccc')
    ax4.axhline(M1.value[0], zorder=-1, color='tab:green',
                alpha=0.5, label='$M_1$', marker='')
    ax4.legend(loc='best')
    ax4.set_ylabel(r'$M_{2,{\rm min}}$' + ' [{0:latex_inline}]'.format(u.Msun))
    ax4.set_xlabel('$P$ [day]')

    for ax in [ax3, ax4]:
        ax.xaxis.set_ticks(10**np.arange(0, 3+0.1, 1))
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

    return fig
