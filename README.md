# PCI<sup>ST</sup>
# NEW CHANGE
A short library for calculating the state transitions Perturbational Complexity Index (PCI<sup>ST</sup>).

The main function of the library is `calc_PCIst()`, which  is composed of two functions corresponding to the two steps involved is the computation of PCI<sup>ST</sup>: `dimensionality reduction()` and `state_transition_quantification()`. The parameters of calc_PCIst() are the inputs of these two functions. The output of `calc_PCIst()` is the PCI<sup>ST</sup> value and a list with the component wise PCI<sup>ST</sup> (∆NST<sub>n</sub>).

## Basic Usage
**TMS/EEG**
```python
from pci_st import **
par = {'baseline_window':(-400,-50), 'response_window':(0,300), 'k':1.2, 'min_snr':1.1, 'max_var':99, 'embed':False,'n_steps':100}
pci = calc_PCIst(signal_evoked, times, **par)
```
**SPES/SEEG**
```python
from pci_st import **
par = {'baseline_window':(-250,-50), 'response_window':(10,600), 'k':1.2, 'min_snr':1.1, 'max_var':99, 'embed':False,'n_steps':100, 'avgref': False}
pci = calc_PCIst(signal_evoked, times, **par)
```

## Credit
**Please cite this paper if you use this code:**

Comolatti et al., "A fast and general method to empirically estimate the complexity of distributed causal interactions in the brain" (to be submitted)

Correspondance regarding the code can be directed to renzo.com@gmail.com 
