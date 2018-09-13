# PCIst
A short library for calculating the state transitions Perturbational Complexity Index (PCIst).

The main function of the library is calc_PCIst(), which  is composed of two auxiliary functions that correspond to the two steps involved is the computation of PCIST: dimensionality reduction() and state_transition_quantification(). The parameters of calc_PCIst() are the inputs of these two functions (Table 1 and Table 2). The output of calc_PCIst() is the PCIST value and a list with the component wise PCIST (∆NSTn).

## Basic Usage
`par = {'baseline_window':(-400,-50), 'response_window':(0,300), 'k':1.2, 'min_snr':1.1, 'max_var':99, 'embed':False,'n_steps':100}
PCIst, PCIst_bydim = calc_PCIst(signal_evoked, times, **par)`
