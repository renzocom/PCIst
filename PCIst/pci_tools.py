import os
import glob
import scipy.io as sio
import sys
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

from collections import OrderedDict
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from PCIst import pci_st as pci


# Test get_dataset_info function
# project_name = 'storm'
# dataset_name = 'test'
# project_dir = Path.home()/'Documents/research'/project_name
# dataset_dir = project_dir/'data'/dataset_name
# filename = '*.rtf'

# data_paths, data_info = get_dataset_info(dataset_dir, filename, n_levels=2)
# data_info, data_paths

def get_par_id(par, name='PCIst'):

    s = [name,par['baseline_window'][0],par['baseline_window'][1],par['response_window'][0],par['response_window'][1],
                 'k',round(par['k'],1),'maxVar',par['max_var'],'nsteps',par['n_steps'],'min_snr',par['min_snr']]
    s = [str(x) for x in s]
    if par['embed']:
        embedding = '_'.join(['m',str(par['m']),'tau',str(par['tau'])])
        s.append(embedding)
    else:
        s.append('m_1')
    s = '_'.join(s)
    return s

def normalize_byvar(signal, times, window):
    ''' Normalize signal (dim x times) baseline to var==1. '''
    d,n = signal.shape
    norm_signal = np.zeros((d,n))
    ix_ini = pci.get_time_index(times, onset=window[0])
    ix_end = pci.get_time_index(times, onset=window[1])
    for i in range(d):
        norm_signal[i,:] = signal[i,:]/np.std(signal[i,ix_ini:ix_end]) # normalize by variance of BASELINE
    return norm_signal

def view_PCIst(signal_evk, times, params, ax=None):
    r = pci.calc_PCIst(signal_evk, times, full_return=True, **params)
    var_exp = r['var_exp']
    signal_svd = normalize_byvar(r['signal_svd'], r['times'], params['baseline_window'])

    if ax is None:
        fig, ax = plt.subplots(figsize=(10,4))

    ax.plot(r['times'], signal_svd.T)
    
    legend_str = [f'$∆NST_{i+1}$: {s[0]:<6.1f} ({s[1]:.1f}%, snr={s[2]:.1f})' for i, s in enumerate(zip(r['dNST'][:10], r['var_exp'][:10], r['snrs'][:10]))]
    plt.xlim(params['baseline_window'][0], params['response_window'][1])
    leg = plt.legend(legend_str, loc=3, fontsize=10)
    title_str = 'PCI = {} x {:.1f} = {:.1f}'.format(r['n_dims'], np.mean(r['dNST']), r['PCI'])
    leg.set_title(title_str,prop={'size':12})

def plot_PCIscatter_bysubject(df, par_id, condition2color, threshold=None, condition_label='Condition', subject_label='subject'):
    for i, (subject, df2) in enumerate(df.groupby(subject_label, sort=False)):
        for cond, df3 in df2.groupby(condition_label):
            plt.plot([i]*len(df3[par_id]), df3[par_id].values, 'o', color='white', markeredgecolor=condition2color[cond])
            plt.plot([i], [df3[par_id].max()], 'o', color=condition2color[cond], markeredgecolor='black')

    plt.xlim(-1,i+1)

    if threshold is not None:
        plt.axhline(threshold, linestyle='--', color='grey')

    # LEGEND
    legendmap = OrderedDict(condition2color)
    legend_elements = [Line2D([0], [0], marker='o',color='white',markerfacecolor=legendmap[cond],label=str(cond), markersize=8) for cond in legendmap.keys()]
    plt.legend(handles=legend_elements, loc=1)

def plot_PCIscatter(df, par_id, condition2color, threshold=None, condition_label='Condition', subject_label='subject'):

    conditions = df[condition_label].unique()
    k = 1
    for i, (condition, df2) in enumerate(df.groupby(condition_label, sort=False)):
        for j, (subject, df3) in enumerate(df2.groupby(subject_label)):
            plt.plot([k]*len(df3[par_id]), df3[par_id].values, 'o', color='white', markeredgecolor='black')
            plt.plot([k], [df3[par_id].max()], 'o', color=condition2color[condition], markeredgecolor='black')
            k += 1
    plt.xlim(-1,k+10)

    if threshold is not None:
        plt.axhline(threshold, linestyle='--', color='grey')

    # LEGEND
    legendmap = OrderedDict(condition2color)
    legend_elements = [Line2D([0], [0], marker='o',color='white',markerfacecolor=legendmap[cond],label=str(cond), markersize=8) for cond in legendmap.keys()]
    plt.legend(handles=legend_elements, loc=1)
