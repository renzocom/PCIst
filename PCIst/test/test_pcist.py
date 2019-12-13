import pytest

from PCIst import pci_st

def test_time2ix():
    times = [-2, -1, 0, 1, 2]

    assert pci_st.time2ix(times, 0) == 2
    assert pci_st.time2ix(times, 0.1) == 2
    assert pci_st.time2ix(times, -0.1) == 2
    assert pci_st.time2ix(times, -4) == 0
    assert pci_st.time2ix(times, 4) == 4