#
# Renzo Comolatti (renzo.com@gmail.com) and Adenauer G. Casali
#
# Please cite this paper if you use this code:
# Comolatti R et al., "A fast and general method to empirically estimate the complexity of brain responses
# to transcranial and intracranial stimulations" Brain Stimulation (in press)
# https://doi.org/10.1016/j.brs.2019.05.013
#
# started: 26/10/2017
# last update: 27/05/2019

import numpy as np
from numpy import linalg
import scipy.signal


def calc_PCIst(signal_evk, times, full_return=False, **par):
    ''' Calculates PCIst (Perturbational Complexity Index based on State transitions) of a signal.
    Parameters
    ----------
    signal_evk : ndarray
        2D array (ch, times) containing signal.
    times : ndarray
        1D array (time,) containing timepoints (negative values are baseline).
    full_return : bool
        Returns multiple variables involved in PCI computation.
    **pars : dictionary
        Dictionary containing parameters (see dimensionality_reduction(),
        state_transition_quantification()and preprocess_signal() documentation).
        Example:
        >> par = {'baseline_window':(-400,-50), 'response_window':(0,300), 'k':1.2, 'min_snr':1.1,
        'max_var':99, 'embed':False,'n_steps':100}
        >> PCIst, PCIst_bydim = calc_PCIst(signal_evoked, times, **par)

    Returns
    -------
    float
        PCIst value
    OR (if full_return==True)
    dict
        Dictionary containing all variables from calculation including array 'dNSTn' with PCIst decomposition.
    '''
    if np.any(np.isnan(signal_evk)):
        print('Data contains nan values.')
        return 0

    signal_evk, times = preprocess_signal(signal_evk, times, (par['baseline_window'][0],
                                                              par['response_window'][1]), **par)
    signal_svd, var_exp, eigenvalues, snrs = dimensionality_reduction(signal_evk, times, **par)
    STQ = state_transition_quantification(signal_svd, times, **par)

    PCI = np.sum(STQ['dNST'])

    if full_return:
        return {'PCI':PCI, **STQ, 'signal_evk':signal_evk, 'times':times, 'signal_svd':signal_svd,
                'eigenvalues':eigenvalues, 'var_exp':var_exp, 'snrs':snrs}
    return PCI

## DIMENSIONALITY REDUCTION
def dimensionality_reduction(signal, times, response_window, max_var=0.99, min_snr=1.1,
                             n_components=None, **kwargs):
    '''Returns principal components of signal according to SVD of the response.

    Calculates SVD at a given time interval (response_window) and uses the new basis to transform
    the whole signal yielding `n_components` principal components. The principal components are
    then selected to account for at least `max_var`% of the variance basesent in the signal's
    response.

    Parameters
    ----------
    signal : ndarray
        2D array (ch,time) containing signal.
    times : ndarray
        1D array (time,) containing timepoints
    response_window : tuple
        Signal's response time interval (ini,end).
    max_var: 0 < float <= 100
        Percentage of variance accounted for by the selected principal components.
    min_snr : float, optional
        Selects principal components with a signal-to-noise ratio (SNR) > min_snr.
    n_components : int, optional
        Number of principal components calculated (before selection).


    Returns
    -------
    np.ndarray
        2D array (ch,time) with selected principal components.
    np.ndarray
        1D array (n_components,) with `n_components` SVD eigenvalues of the signal's response.
    '''

    if not n_components:
        n_components = signal.shape[0]

    Vk, eigenvalues = get_svd(signal, times, response_window, n_components)
    var_exp = 100 * eigenvalues**2/np.sum(eigenvalues**2)

    signal_svd = apply_svd(signal, Vk)

    max_dim = calc_maxdim(eigenvalues, max_var)

    signal_svd = signal_svd[:max_dim, :]

    # if min_snr:
        # base_ini_ix = get_time_index(times, kwargs['baseline_window'][0])
        # base_end_ix = get_time_index(times, kwargs['baseline_window'][1])
        # resp_ini_ix = get_time_index(times, response_window[0])
        # resp_end_ix = get_time_index(times, response_window[1])
        # n_dims = np.size(signal_svd, 0)
        # snrs = np.zeros(n_dims)
        # for c in range(n_dims):
        #     resp_power = np.mean(np.square(signal_svd[c, resp_ini_ix:resp_end_ix]))
        #     base_power = np.mean(np.square(signal_svd[c, base_ini_ix:base_end_ix]))
        #     snrs[c] = np.sqrt(np.divide(resp_power, base_power))
    snrs = calc_snr(signal_svd, times, kwargs['baseline_window'], response_window)
    signal_svd = signal_svd[snrs > min_snr, :]
    snrs = snrs[snrs > min_snr]

    Nc = signal_svd.shape[0]

    return signal_svd, var_exp[:Nc], eigenvalues, snrs

def calc_snr(signal_svd, times, baseline_window, response_window):

    base_ini_ix = get_time_index(times, baseline_window[0])
    base_end_ix = get_time_index(times, baseline_window[1])
    resp_ini_ix = get_time_index(times, response_window[0])
    resp_end_ix = get_time_index(times, response_window[1])

    resp_power = np.mean(np.square(signal_svd[:,resp_ini_ix:resp_end_ix]), axis=1)
    base_power = np.mean(np.square(signal_svd[:,base_ini_ix:base_end_ix]), axis=1)
    snrs = np.sqrt(resp_power / base_power)
    return snrs

def get_svd(signal_evk, times, response_window, n_components):
    ini_t, end_t = response_window
    ini_ix = get_time_index(times, onset=ini_t)
    end_ix = get_time_index(times, onset=end_t)
    signal_resp = signal_evk[:, ini_ix:end_ix].T
    U, S, V = linalg.svd(signal_resp, full_matrices=False)
    V = V.T
    Vk = V[:, :n_components]
    eigenvalues = S[:n_components]
    return Vk, eigenvalues

def apply_svd(signal, V):
    '''Transforms signal according to SVD basis.'''
    return signal.T.dot(V).T

## STATE TRANSITION QUANTIFICATION
def state_transition_quantification(signal, times, k, baseline_window, response_window, embed=False,
                                    L=None, tau=None, n_steps=100, max_thr_p=1.0, **kwargs):
    ''' Receives selected principal components of perturbational signal and
    performs state transition quantification.

    Parameters
    ----------
    signal : ndarray
        2D array (component,time) containing signal (typically, the selected
        principal components).
    times : ndarray
        1D array (time,) containing timepoints
    k : float > 1
        Noise control parameter.
    baseline_window : tuple
        Signal's baseline time interval (ini,end).
    response_window : tuple
        Signal's response time interval (ini,end).
    embed : bool, optional
        Perform time-delay embedding.
    L : int
        Number of embedding dimensions.
    tau : int
        Number of timesamples of embedding delay
    n_steps : int, optional
        Number of steps used to search for the threshold that maximizes ∆NST.
        Search is performed by partitioning  the interval (defined from the median
        of the baseline’s distance matrix to the maximum of the response’s
        distance matrix) into ‘n_steps’ equal lengths.

    Returns
    -------
    float
        PCIst value.
    ndarray
        List containing component wise PCIst value (∆NSTn).
    '''

    n_dims = signal.shape[0]
    if n_dims == 0:
        print('No components --> PCIst=0')
        return {'dNST':np.array([]), 'n_dims':0}

    # EMBEDDING
    if embed:
        cut = (L-1)*tau
        times = times[cut:]
        temp_signal = np.zeros((n_dims, L, len(times)))
        for i in range(n_dims):
            temp_signal[i, :, :] = dimension_embedding(signal[i, :], L, tau)
        signal = temp_signal

    else:
        signal = signal[:, np.newaxis, :]

    # BASELINE AND RESPONSE DEFINITION
    base_ini_ix = get_time_index(times, baseline_window[0])
    base_end_ix = get_time_index(times, baseline_window[1])
    resp_ini_ix = get_time_index(times, response_window[0])
    resp_end_ix = get_time_index(times, response_window[1])
    n_baseline = len(times[base_ini_ix:base_end_ix])
    n_response = len(times[resp_ini_ix:resp_end_ix])

    if n_response <= 1 or n_baseline <= 1:
        print('Warning: Bad time interval defined.')

    baseline = signal[:, :, base_ini_ix:base_end_ix]
    response = signal[:, :, resp_ini_ix:resp_end_ix]

    # NST CALCULATION
        # Distance matrix
    D_base = np.zeros((n_dims, n_baseline, n_baseline))
    D_resp = np.zeros((n_dims, n_response, n_response))
        # Transition matrix
    T_base = np.zeros((n_steps, n_dims, n_baseline, n_baseline))
    T_resp = np.zeros((n_steps, n_dims, n_response, n_response))
        # Number of state transitions
    NST_base = np.zeros((n_steps, n_dims))
    NST_resp = np.zeros((n_steps, n_dims))
    thresholds = np.zeros((n_steps, n_dims))
    for i in range(n_dims):
        D_base[i, :, :] = recurrence_matrix(baseline[i, :, :], thr=None, mode='distance')
        D_resp[i, :, :] = recurrence_matrix(response[i, :, :], thr=None, mode='distance')
        min_thr = np.median(D_base[i, :, :].flatten())
        max_thr = np.max(D_resp[i, :, :].flatten()) * max_thr_p
        thresholds[:, i] = np.linspace(min_thr, max_thr, n_steps)
    for i in range(n_steps):
        for j in range(n_dims):
            T_base[i, j, :, :] = distance2transition(D_base[j, :, :], thresholds[i, j])
            T_resp[i, j, :, :] = distance2transition(D_resp[j, :, :], thresholds[i, j])

            NST_base[i, j] = np.sum(T_base[i, j, :, :])/n_baseline**2
            NST_resp[i, j] = np.sum(T_resp[i, j, :, :])/n_response**2

    # PCIST
    NST_diff = NST_resp - k * NST_base
    ixs = np.argmax(NST_diff, axis=0)
    max_thresholds = np.array([thresholds[ix, i] for ix, i in zip(ixs, range(n_dims))])
    dNST = np.array([NST_diff[ix, i] for ix, i in zip(ixs, range(n_dims))]) * n_response
    dNST = [x if x>0 else 0 for x in dNST]

    temp = np.zeros((n_dims, n_response, n_response))
    temp2 = np.zeros((n_dims, n_baseline, n_baseline))
    for i in range(n_dims):
        temp[i, :, :] = T_resp[ixs[i], i, :, :]
        temp2[i, :, :] = T_base[ixs[i], i, :, :]
    T_resp = temp
    T_base = temp2

    return {'dNST':dNST, 'n_dims':n_dims,
    'D_base':D_base, 'D_resp':D_resp, 'T_base':T_base,'T_resp':T_resp,
    'thresholds':thresholds, 'NST_diff':NST_diff, 'NST_resp':NST_resp, 'NST_base':NST_base,'max_thresholds':max_thresholds}


def recurrence_matrix(signal, mode, thr=None):
    ''' Calculates distance, recurrence or transition matrix. Signal can be
    embedded (m, n_times) or not (, n_times).

    Parameters
    ----------
    signal : ndarray
        Time-series; may be a 1D (time,) or a m-dimensional array (m, time) for
        time-delay embeddeding.
    mode : str
        Specifies calculated matrix: 'distance', 'recurrence' or 'transition'
    thr : float, optional
        If transition matrix is chosen (`mode`=='transition'), specifies threshold value.

    Returns
    -------
    ndarray
        2D array containing specified matrix.
    '''
    if len(signal.shape) == 1:
        signal = signal[np.newaxis, :]
    n_dims = signal.shape[0]
    n_times = signal.shape[1]

    R = np.zeros((n_dims, n_times, n_times))
    for i in range(n_dims):
        D = np.tile(signal[i, :], (n_times, 1))
        D = D - D.T
        R[i, :, :] = D
    R = np.linalg.norm(R, ord=2, axis=0)

    mask = (R <= thr) if thr else np.zeros(R.shape).astype(bool)
    if mode == 'distance':
        R[mask] = 0
        return R
    if mode == 'recurrence':
        return mask.astype(int)
    if mode == 'transition':
        return diff_matrix(mask.astype(int), symmetric=False)
    return 0

def distance2transition(dist_R, thr):
    ''' Receives 2D distance matrix and calculates transition matrix. '''
    mask = dist_R <= thr
    R = diff_matrix(mask.astype(int), symmetric=False)
    return R

def distance2recurrence(dist_R, thr):
    ''' Receives 2D distance matrix and calculates recurrence matrix. '''
    mask = dist_R <= thr
    return mask.astype(int)

def diff_matrix(A, symmetric=False):
    B = np.abs(np.diff(A))
    if B.shape[1] != B.shape[0]:
        B2 = np.zeros((B.shape[0], B.shape[1]+1))
        B2[:, :-1] = B
        B = B2
    if symmetric:
        B = (B + B.T)
        B[B > 0] = 1
    return B

def calc_maxdim(eigenvalues, max_var):
    ''' Get number of dimensions that accumulates at least `max_var`% of total variance'''
    if max_var == 100:
        return len(eigenvalues)
    eigenvalues = np.sort(eigenvalues)[::-1] # Sort in descending order
    var = eigenvalues ** 2
    var_p = 100 * var/np.sum(var)
    var_cum = np.cumsum(var_p)
    max_dim = len(eigenvalues) - np.sum(var_cum >= max_var) + 1
    return max_dim

def dimension_embedding(x, L, tau):
    '''
    Returns time-delay embedding of vector.
    Parameters
    ----------
    x : ndarray
        1D array time series.
    L : int
        Number of dimensions in the embedding.
    tau : int
        Number of samples in delay.
    Returns
    -------
    ndarray
        2D array containing embedded signal (L, time)

    '''
    assert len(x.shape) == 1, "x must be one-dimensional array (n_times,)"
    n_times = x.shape[0]
    s = np.zeros((L, n_times - (L-1) * tau))
    ini = (L-1) * tau if L > 1 else None
    s[0, :] = x[ini:]
    for i in range(1, L):
        ini = (L-i-1) * tau
        end = -i * tau
        s[i, :] = x[ini:end]
    return s

## PREPROCESS
def preprocess_signal(signal_evk, times, time_window, baseline_corr=False, resample=None,
                      avgref=False, **kwargs):
    assert signal_evk.shape[1] == len(times), 'Signal and Time arrays must be of the same size.'
    if avgref:
        signal_evk = avgreference(signal_evk)
    if baseline_corr:
        signal_evk = baseline_correct(signal_evk, times, delta=-50)
    t_ini, t_end = time_window
    ini_ix = get_time_index(times, t_ini)
    end_ix = get_time_index(times, t_end)
    signal_evk = signal_evk[:, ini_ix:end_ix]
    times = times[ini_ix:end_ix]
    if resample:
        signal_evk, times = undersample_signal(signal_evk, times, new_fs=resample)
    return signal_evk, times

def avgreference(signal):
    ''' Performs average reference to signal. '''
    new_signal = np.zeros(signal.shape)
    channels_mean = np.mean(signal, axis=0)[np.newaxis]
    new_signal = signal - channels_mean
    return new_signal

def undersample_signal(signal, times, new_fs):
    '''
    signal : (ch x times)
    times : (times,) [ms]
    new_fs : [hz]
    '''
    n_samples = int((times[-1]-times[0])/1000 * new_fs)
    new_signal_evk, new_times = scipy.signal.resample(signal, n_samples, t=times, axis=1)
    return new_signal_evk, new_times

def baseline_correct(Y, times, delta=0):
    ''' Baseline correct signal using times < delta '''
    newY = np.zeros(Y.shape)
    onset_ix = get_time_index(times, delta)
    baseline_mean = np.mean(Y[:, :onset_ix], axis=1)[np.newaxis]
    newY = Y - baseline_mean.T
    close_enough = np.all(np.isclose(np.mean(newY[:, :onset_ix], axis=1), 0, atol=1e-08))
    assert close_enough, "Baseline mean is not zero"
    return newY

def get_time_index(times, onset=0):
    ''' Returns index of first time greater then delta. For delta=0 gets index of
    first non-negative time.
    '''
    return np.sum(times < onset)
