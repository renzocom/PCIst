import pytest
from pytest import approx
import numpy as np
import scipy.io as sio

from PCIst import pci_st

def test_time2ix():
    times = [-2, -1, 0, 1, 2]
    assert pci_st.time2ix(times, 0) == 2
    assert pci_st.time2ix(times, 0.1) == 2
    assert pci_st.time2ix(times, -0.1) == 2
    assert pci_st.time2ix(times, -4) == 0
    assert pci_st.time2ix(times, 4) == 4

def test_recurrence_matrix():
    d = sio.loadmat("sample_tms.mat")
    evoked, times = d['evoked'], d['times']
    resp_evoked = pci_st.trim_signal(evoked, times, (0,300))
    signal = resp_evoked[0, :]

    D = pci_st.recurrence_matrix(signal, 'distance')
    D2 = pci_st.recurrence_matrix(signal, 'distance', thr=2)
    R = pci_st.recurrence_matrix(signal, 'recurrence', thr=2)
    T = pci_st.recurrence_matrix(signal, 'transition', thr=2)
    T2 = pci_st.diff_matrix(R, symmetric=True) # Transition matrix with derivative in 2D

    assert approx(np.sum(D)) == 85133.84059843828
    assert approx(np.sum(D2)) == 20379.882743126567
    assert approx(np.sum(R)) == 81422
    assert approx(np.sum(T)) == 1176.0
    assert approx(np.sum(T2)) == 2020.0

def test_state_transition_quantification():
    d = sio.loadmat("sample_tms.mat")
    evoked, times = d['evoked'], np.squeeze(d['times'])
    signal = evoked[[0,1], :]
    dNST = pci_st.state_transition_quantification(signal, times, k=1.2, baseline_window=(-400, -50), response_window=(0,300), n_steps=100)['dNST']

    assert len(dNST)==2
    assert approx(np.mean(dNST)) == 10.827285714285713
