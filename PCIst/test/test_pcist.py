import pytest
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
    R = pci_st.recurrence_matrix(signal, 'recurrence', thr=2)
    T = pci_st.recurrence_matrix(signal, 'transition', thr=2)
    T2 = pci_st.diff_matrix(R, symmetric=True) # Transition matrix with derivative in 2D

    assert np.sum(D) == 85133.84059843828
    assert np.sum(R) == 81422
    assert np.sum(T) == 1176.0
    assert np.sum(T2) == 2020.0