import csv

import numpy as np

from ms_deisotope.data_source.scan.base import RawDataArrays, ChargeNotProvided, PrecursorInformation
from ms_deisotope.data_source.memory import make_scan
from ms_deisotope.data_source._compression import get_opener


def scan_from_csv(file_handle, delimiter=',', ms_level=2, is_profile=True, polarity=1,
                  precursor_mz=None, precursor_charge=None):
    """Read an m/z-intensity point list from a text file stream.

    Parameters
    ----------
    file_handle : file
        The file to read from.
    delimiter : str, optional
        The field separator between m/z and intensity. (the default is ',')
    ms_level : int, optional
        The MS level of the read scan. (the default is 2)
    is_profile : bool, optional
        Whether the scan is in profile mode (the default is True)
    polarity : int, optional
        Whether the scan is positive mode or negative (the default is 1, which means positive)
    precursor_mz : float, optional
        The m/z of the precursor ion to record. If provided, it will be
        provided in :attr:`~.Scan.precursor_information`.
    precursor_charge : int, optional
        If provided, and `precursor_mz` is not :const:`None`, then it will be provided
        in :attr:`~.Scan.precursor_information`.

    Returns
    -------
    :class:`~.Scan`
    """
    file_handle = get_opener(file_handle)
    reader = csv.reader(file_handle, delimiter=delimiter)
    mzs = []
    intensities = []
    for i, line in enumerate(reader):
        if i == 0:
            try:
                float(line[0])
            except Exception:
                continue
        mzs.append(line[0])
        intensities.append(line[1])
    mzs = np.array(mzs, dtype=float)
    intensities = np.array(intensities, dtype=float)
    signal = RawDataArrays(mzs, intensities)
    pinfo = None
    if precursor_mz is not None:
        pinfo = PrecursorInformation(
            precursor_mz,
            0,
            ChargeNotProvided if precursor_charge is None else precursor_charge)
    scan = make_scan(
        signal, ms_level, "index=1", 0, 0, is_profile,
        polarity,
        precursor_information=pinfo)
    return scan
