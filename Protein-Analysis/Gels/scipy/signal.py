"""Tiny shim providing only find_peaks (distance + scalar prominence), for
gel_annotator.py in this sandbox where real scipy isn't installable.
Reimplements scipy's local-maxima / distance-suppression / prominence
algorithm closely enough for 1D intensity-profile peak finding."""
import numpy as np


def _local_maxima_1d(x):
    midpoints = []
    i = 1
    i_max = len(x) - 1
    while i < i_max:
        if x[i - 1] < x[i]:
            i_ahead = i + 1
            while i_ahead < i_max and x[i_ahead] == x[i]:
                i_ahead += 1
            if x[i_ahead] < x[i]:
                left_edge = i
                right_edge = i_ahead - 1
                midpoints.append((left_edge + right_edge) // 2)
                i = i_ahead
        i += 1
    return np.array(midpoints, dtype=int)


def _select_by_peak_distance(peaks, priority, distance):
    n = len(peaks)
    distance_ = int(np.ceil(distance))
    keep = np.ones(n, dtype=bool)
    order = np.argsort(priority)
    for idx in range(n - 1, -1, -1):
        j = order[idx]
        if not keep[j]:
            continue
        k = j - 1
        while k >= 0 and peaks[j] - peaks[k] < distance_:
            keep[k] = False
            k -= 1
        k = j + 1
        while k < n and peaks[k] - peaks[j] < distance_:
            keep[k] = False
            k += 1
    return keep


def _peak_prominences(x, peaks):
    n = len(x)
    prominences = np.empty(len(peaks), dtype=float)
    left_bases = np.empty(len(peaks), dtype=int)
    right_bases = np.empty(len(peaks), dtype=int)
    for idx, peak in enumerate(peaks):
        i = peak
        left_min = x[peak]
        while i > 0:
            i -= 1
            if x[i] > x[peak]:
                break
            if x[i] < left_min:
                left_min = x[i]
        left_base = i
        j = peak
        right_min = x[peak]
        while j < n - 1:
            j += 1
            if x[j] > x[peak]:
                break
            if x[j] < right_min:
                right_min = x[j]
        right_base = j
        prominences[idx] = x[peak] - max(left_min, right_min)
        left_bases[idx] = left_base
        right_bases[idx] = right_base
    return prominences, left_bases, right_bases


def find_peaks(x, distance=None, prominence=None):
    x = np.asarray(x, dtype=float)
    peaks = _local_maxima_1d(x)
    properties = {}
    if distance is not None and len(peaks) > 1:
        keep = _select_by_peak_distance(peaks, x[peaks], distance)
        peaks = peaks[keep]
    prominences, left_bases, right_bases = _peak_prominences(x, peaks)
    properties["prominences"] = prominences
    properties["left_bases"] = left_bases
    properties["right_bases"] = right_bases
    if prominence is not None:
        if np.isscalar(prominence):
            pmin, pmax = prominence, None
        else:
            pmin, pmax = prominence
        keep = np.ones(len(peaks), dtype=bool)
        if pmin is not None:
            keep &= prominences >= pmin
        if pmax is not None:
            keep &= prominences <= pmax
        peaks = peaks[keep]
        for k in properties:
            properties[k] = properties[k][keep]
    return peaks, properties
