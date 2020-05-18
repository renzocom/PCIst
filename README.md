# PCI<sup>ST</sup> (Python and Matlab)
A short library for calculating the state transitions Perturbational Complexity Index (PCI<sup>ST</sup>).

The main function of the `python` library is `calc_PCIst()`, which  is composed of two functions corresponding to the two steps involved is the computation of PCI<sup>ST</sup>: `dimensionality reduction()` and `state_transition_quantification()`. The parameters of calc_PCIst() are the inputs of these two functions. The output of `calc_PCIst()` is the PCI<sup>ST</sup> value and a list with the component wise PCI<sup>ST</sup> (∆NST<sub>n</sub>).

## Basic Usage
### Python
**TMS/EEG**
```python
from pci_st import **
par = {'baseline_window':(-400,-50), 'response_window':(0,300), 'k':1.2, 'min_snr':1.1, 'max_var':99, 'embed':False,'n_steps':100} # 
pci = calc_PCIst(evoked, times, **par)
```
**SPES/SEEG**
```python
from pci_st import **
par = {'baseline_window':(-250,-50), 'response_window':(10,600), 'k':1.2, 'min_snr':1.1, 'max_var':99, 'embed':False,'n_steps':100, 'avgref': False}
pci = calc_PCIst(evoked, times, **par)
```

### Matlab
**TMS/EEG**
```matlab
parameters = []; % Default parameters for TMS/EEG used in (Brain Stim, 2019)
[pci,dNST] = PCIst(evoked, times, parameters)

% Equivalent to:

par=struct('baseline',[-400 -50],'response',[0 300],'k',1.2,'min_snr',1.1,'max_var',99,'l',1,'nsteps',100);
[pci,dNST] = PCIst(evoked, times, par)
```
**SPES/SEEG**
```matlab
par=struct('baseline',[-250 -50],'response',[10 600],'k',1.2,'min_snr',1.1,'max_var',99,'l',1,'nsteps',100);
[pci,dNST] = PCIst(evoked, times, par)
```

## Credit
**Please cite this paper if you use this code:**

Comolatti R et al., "A fast and general method to empirically estimate the complexity of brain responses to transcranial and intracranial stimulations" Brain Stimulation (2019) https://doi.org/10.1016/j.brs.2019.05.013

Correspondance regarding the code can be directed to renzo.com@gmail.com 
