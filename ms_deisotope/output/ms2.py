# pragma: no cover

from ms_deisotope import mass_charge_ratio

from .text import TextScanSerializerBase


class MS2Serializer(TextScanSerializerBase):
    def __init__(self, stream, sample_name=None, deconvoluted=True):
            super(MS2Serializer, self).__init__(stream, deconvoluted)
            self.sample_name = sample_name
            self.started = False

    def add_global_parameter(self, name, value):
        if self.started:
            raise ValueError("Cannot add global parameter if scan data has begun being written")
        self.stream.write("H\t{0}\t{1}\n".format(name, value))

    def _format_parameter(self, key, value, symbol):
        return "{2}\t{0}\t{1}\n".format(str(key).upper(), str(value), symbol)

    def _add_parameter(self, name, value, symbol):
        self.stream.write(self._format_parameter(name, value, symbol))

    def add_parameter(self, name, value, symbol="I"):
        self._add_parameter(name, value, symbol)

    def save_scan_bunch(self, bunch, **kwargs):
        for scan in bunch.products:
            self.write_scan(*self.prepare_scan_data(scan))

    def format_peak_vectors(self, scan):
        if self.deconvoluted:
            mz = [i.mz for p in scan.deconvoluted_peak_set for i in p.envelope]
            intensity = [i.intensity for p in scan.deconvoluted_peak_set for i in p.envelope]
            return (mz, intensity, None)
        else:
            mz = [p.mz for p in scan.peak_set]
            intensity = [p.intensity for p in scan.peak_set]
            return mz, intensity, None

    def write_header(self, header_dict):
        mass = header_dict['precursor_neutral_mass']
        charge = header_dict['precursor_charge']
        mz = mass_charge_ratio(mass, charge)
        second_index = header_dict.scan.index
        try:
            first_index = header_dict.scan.precursor_information.precursor.index
        except AttributeError:
            first_index = second_index
        self.stream.write("S\t%r\t%r\t%r\n" % (first_index, second_index, mz))
        self.stream.write("Z\t%r\t%r\n" % (charge, mass))
