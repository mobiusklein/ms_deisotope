
from ms_deisotope.averagine import neutral_mass, mass_charge_ratio
from ms_deisotope.peak_set import DeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.data_source.mgf import MGFLoader, mgf as pymgf

from .text import HeaderedDelimitedWriter


def _format_parameter(key, value):
    return "{0}={1}\n".format(str(key).upper(), str(value))


class MGFSerializer(HeaderedDelimitedWriter):
    def __init__(self, stream, sample_name=None, deconvoluted=True):
        super(MGFSerializer, self).__init__(stream, deconvoluted)
        self.sample_name = sample_name
        self.started = False

    def add_global_parameter(self, name, value):
        if self.started:
            raise ValueError("Cannot add global parameter if scan data has begun being written")
        self._add_parameter(name, value)

    def _add_parameter(self, name, value):
        self.stream.write(_format_parameter(name, value))

    def add_parameter(self, name, value):
        self._add_parameter(name, value)

    def save_scan_bunch(self, bunch):
        for scan in bunch.products:
            self.write_scan(*self.prepare_scan_data(scan))

    def format_peak_vectors(self, scan):
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
        pepmass = header_dict['precursor_neutral_mass']
        charge = header_dict['precursor_charge']
        intensity = header_dict['precursor_intensity']
        self.add_parameter("pepmass", "%f %f" % (pepmass, intensity))
        try:
            self.add_parameter("charge", "%d%s" % (charge, "+" if header_dict['polarity'] > 0 else '-'))
        except TypeError:
            pass
        self.add_parameter("title", header_dict['title'])
        self.add_parameter("rtinseconds", header_dict['scan_time'] * 60.0)

    def write_scan(self, scan_header, data_vectors):
        self.stream.write('BEGIN IONS\n')
        self.write_header(scan_header)
        self.write_vectors(data_vectors)
        self.stream.write('END IONS\n')


class ProcessedMGFDeserializer(MGFLoader):
    def _create_parser(self):
        return pymgf.read(self.source_file, read_charges=True,
                          convert_arrays=1, encoding=self.encoding)

    def _build_peaks(self, scan):
        mz_array = scan['m/z array']
        intensity_array = scan["intensity array"]
        charge_array = scan['charge array']
        peaks = []
        for i in range(len(mz_array)):
            peak = DeconvolutedPeak(
                neutral_mass(mz_array[i], charge_array[i]), intensity_array[i], charge_array[i],
                intensity_array[i], i, 0)
            peaks.append(peak)
        peak_set = DeconvolutedPeakSet(peaks)
        peak_set.reindex()
        return peak_set

    def _make_scan(self, scan):
        scan = super(ProcessedMGFDeserializer, self)._make_scan(scan)
        scan.peak_set = None
        scan.deconvoluted_peak_set = self._build_peaks(scan._data)
        return scan.pack()
