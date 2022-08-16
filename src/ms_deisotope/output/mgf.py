"""Code for writing Mascot Generic Format files and for reading deconvoluted
generic format files.

"""

from ms_deisotope.averagine import neutral_mass, mass_charge_ratio
from ms_deisotope.peak_set import DeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.data_source.mgf import MGFLoader, mgf as pymgf, _MGFParser

from .text import HeaderedDelimitedWriter


def _format_parameter(key, value):
    return "{0}={1}\n".format(str(key).upper(), str(value))


class MGFSerializer(HeaderedDelimitedWriter):
    """A MASCOT Generic Format writer which is deconvolution aware.

    Implements :class:`~.HeaderedDelimitedWriter, as well as allowing
    global parameters to be included at the beginning of the file.

    Attributes
    ----------
    sample_name: str
        The name of the sample. This value is currently ignored.
    started: bool
        Whether scans have been written out yet or not. Once :attr:`started` is
        true, it is no longer possible to add global parameters using :meth:`add_global_parameter`.

    """

    file_extensions = {
        "mgf",
        "mgf.gz"
    }

    def __init__(self, stream, sample_name=None, deconvoluted=True):
        super(MGFSerializer, self).__init__(stream, deconvoluted)
        self.sample_name = sample_name
        self.started = False

    def add_global_parameter(self, name, value):
        """Add a global parameter at the beginning of the file, before scans
        are written.

        Global parameters follow the same ``KEY=value`` format as scan-level
        parameters. New global parameters may not be added once a scan has been
        written.

        Parameters
        ----------
        name : str
            The parameter name. Will be made upper-case.
        value : object
            The parameter's value. Will be converted to a string.

        Raises
        ------
        ValueError:
            If scans have already been written, this method raises a :class:`ValueError`
        """
        if self.started:
            raise ValueError("Cannot add global parameter if scan data has begun being written")
        self._add_parameter(name, value)

    def _add_parameter(self, name, value):
        self.stream.write(_format_parameter(name, value).encode('utf-8'))

    def add_parameter(self, name, value):
        """Add a parameter to the current block.

        Parameters are written ``KEY=value`` format.

        Parameters
        ----------
        name : str
            The parameter name. Will be made upper-case.
        value : object
            The parameter's value. Will be converted to a string.
        """
        self._add_parameter(name, value)

    def save_scan_bunch(self, bunch, **kwargs):
        for scan in bunch.products:
            self.save_scan(scan, **kwargs)

    def format_peak_vectors(self, scan):
        """As in :class:`~.HeaderedDelimitedTextWriter` but always writes m/z, even when charge
        is known.
        """
        if self.deconvoluted:
            (neutral_mass_array, intensity_array, charge_array) = super(
                MGFSerializer, self).format_peak_vectors(scan)
            mz_array = [mass_charge_ratio(
                neutral_mass_array[i], charge_array[i]) for i in range(len(charge_array))]
        else:
            (mz_array, intensity_array, charge_array) = super(
                MGFSerializer, self).format_peak_vectors(scan)
        return (mz_array, intensity_array, charge_array)

    def write_header(self, header_dict):
        pepmass = header_dict['precursor_mz']
        charge = header_dict['precursor_charge']
        intensity = header_dict['precursor_intensity']
        self.add_parameter("pepmass", "%f %f" % (pepmass, intensity))
        polarity = header_dict['polarity']
        if polarity is None:
            polarity = 1
        try:
            self.add_parameter("charge", "%d%s" % (charge, "+" if polarity > 0 else '-'))
        except TypeError:
            pass
        self.add_parameter("title", header_dict['title'])
        self.add_parameter("rtinseconds", header_dict['scan_time'] * 60.0)

    def write_scan(self, scan_header, data_vectors):
        self.stream.write(b'BEGIN IONS\n')
        self.write_header(scan_header)
        self.write_vectors(data_vectors)
        self.stream.write(b'END IONS\n')


class ProcessedMGFLoader(MGFLoader):
    """A variant of the :class:`~.MGFLoader` that reads deconvoluted mass spectra
    with a charge value for each peak, and looks for additional annotation of the
    precursor ion.

    """

    file_extensions = {
        "mgf",
        "mgf.gz"
    }

    def __init__(self, source_file, encoding='ascii', use_index=True, ** kwargs):
        super(ProcessedMGFDeserializer, self).__init__(
            source_file, encoding, use_index, **kwargs)

    def _create_parser(self):
        if self._use_index:
            return _MGFParser(self.source_file, read_charges=True,
                              convert_arrays=1, encoding=self.encoding)
        else:
            return pymgf.MGF(self.source_file, read_charges=True,
                             convert_arrays=1, encoding=self.encoding)

    def _build_peaks(self, scan):
        mz_array = scan['m/z array']
        intensity_array = scan["intensity array"]
        charge_array = scan['charge array']
        return build_deconvoluted_peak_set_from_arrays(mz_array, intensity_array, charge_array)

    def _make_scan(self, data):
        scan = super(ProcessedMGFDeserializer, self)._make_scan(data)
        scan.peak_set = None
        scan.deconvoluted_peak_set = self._build_peaks(scan._data)
        return scan.pack(bind=True)

    def _precursor_information(self, scan):
        pinfo = super(ProcessedMGFDeserializer, self)._precursor_information(scan)
        defaulted = pinfo.defaulted
        orphan = pinfo.orphan
        pinfo.default()
        pinfo.defaulted = defaulted
        pinfo.orphan = orphan
        return pinfo


ProcessedMGFDeserializer = ProcessedMGFLoader


def build_deconvoluted_peak_set_from_arrays(mz_array, intensity_array, charge_array):
    peaks = []
    for i in range(len(mz_array)):
        peak = DeconvolutedPeak(
            neutral_mass(mz_array[i], charge_array[i]), intensity_array[i], charge_array[i],
            intensity_array[i], i, 0)
        peaks.append(peak)
    peak_set = DeconvolutedPeakSet(peaks)
    peak_set.reindex()
    return peak_set


try:
    _build_deconvoluted_peak_set_from_arrays = build_deconvoluted_peak_set_from_arrays
    from ms_deisotope._c.utils import build_deconvoluted_peak_set_from_arrays
except ImportError:
    pass
