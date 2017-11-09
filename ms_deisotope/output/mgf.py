from .text import HeaderedDelimitedWriter


def _format_parameter(key, value):
    return "{0}={1}\n".format(str(key).upper(), str(value))


class MGFSerializer(HeaderedDelimitedWriter):
    def __init__(self, stream, sample_name=None, deconvoluted=True):
        super(MGFSerializer, self).__init__(stream)
        self.sample_name = sample_name
        self.started = False

    def add_global_parameter(self, name, value):
        if self.started:
            raise ValueError("Cannot add global parameter if scan data has begun being written")
        self._add_parameter(name, value)

    def _add_parameter(self, name, value):
        self.stream.write(_format_parameter(name, value))

    def save_scan_bunch(self, bunch):
        for scan in bunch.products:
            self.write_scan(*self.prepare_scan_data(scan))

    def write_header(self, header_dict):
        pepmass = header_dict['precursor_neutral_mass']
        charge = header_dict['precursor_charge']
        intensity = header_dict['precursor_intensity']
        self._add_parameter("pepmass", "%f %f" % (pepmass, intensity))
        try:
            self._add_parameter("charge", "%d%s" % (charge, "+" if header_dict['polarity'] > 0 else '-'))
        except TypeError:
            pass
        self._add_parameter("title", header_dict['title'])
        self._add_parameter("rtinseconds", header_dict['scan_time'] * 60.0)

    def write_scan(self, scan_header, data_vectors):
        self.stream.write('BEGIN IONS\n')
        self.write_header(scan_header)
        self.write_vectors(data_vectors)
        self.stream.write('END IONS\n')
