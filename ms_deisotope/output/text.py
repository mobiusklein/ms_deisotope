from .common import ScanSerializerBase


class HeaderInformation(dict):
    def __init__(self, scan):
        dict.__init__(self)
        self.scan = scan


class TextScanSerializerBase(ScanSerializerBase):
    def __init__(self, stream, deconvoluted=True):
        self.stream = stream
        self.deconvoluted = deconvoluted

    def close(self):
        self.stream.close()

    def save_scan_bunch(self, bunch, **kwargs):
        self.write_scan(*self.prepare_scan_data(bunch.precursor))
        for scan in bunch.products:
            self.write_scan(*self.prepare_scan_data(scan))

    def construct_header(self, scan):
        prec_info = scan.precursor_information
        header_dict = HeaderInformation(scan)
        if prec_info is not None:
            header_dict['precursor_neutral_mass'] = prec_info.extracted_neutral_mass
            header_dict['precursor_charge'] = prec_info.extracted_charge
            header_dict['precursor_intensity'] = prec_info.extracted_intensity
            header_dict['precursor_scan_id'] = prec_info.precursor_scan_id
            header_dict['defaulted'] = prec_info.defaulted
            activation = scan.activation
            if activation:
                header_dict['precursor_activation_method'] = activation.method
                header_dict['precursor_activation_energy'] = activation.energy
                if activation.is_multiple_dissociation():
                    header_dict['precursor_activation_method'] = ';'.join([str(m) for m in activation.methods])
                    header_dict['precursor_activation_energy'] = ';'.join([str(m) for m in activation.energies])
                for key, value in activation.data.items():
                    header_dict['precursor_activation_' + key.lower()] = value
        header_dict['scan_time'] = scan.scan_time
        header_dict['id'] = scan.id
        header_dict['title'] = scan.title
        header_dict['ms_level'] = scan.ms_level
        header_dict['polarity'] = scan.polarity
        header_dict['index'] = scan.index
        instrument_config = scan.instrument_configuration
        if instrument_config:
            header_dict['analyzers'] = ';'.join(a.name for a in instrument_config.analyzers)
        if scan.has_ion_mobility():
            header_dict['drift_time'] = scan.drift_time
        header_dict['annotations'] = scan.annotations
        return header_dict

    def format_peak_vectors(self, scan):
        if self.deconvoluted:
            neutral_mass = [p.neutral_mass for p in scan.deconvoluted_peak_set]
            intensity = [p.intensity for p in scan.deconvoluted_peak_set]
            charge = [p.charge for p in scan.deconvoluted_peak_set]
            return (neutral_mass, intensity, charge)
        else:
            mz = [p.mz for p in scan.peak_set]
            intensity = [p.intensity for p in scan.peak_set]
            return mz, intensity, None

    def prepare_scan_data(self, scan):
        header_dict = self.construct_header(scan)
        vectors = self.format_peak_vectors(scan)
        return (header_dict, vectors)

    def write_scan(self, scan_header, data_vectors):
        raise NotImplementedError()


class HeaderedDelimitedWriter(TextScanSerializerBase):
    def __init__(self, stream, deconvoluted=True):
        super(HeaderedDelimitedWriter, self).__init__(stream, deconvoluted)

    def write_header(self, header_dict):
        for key, value in header_dict.items():
            self.stream.write("#%s=%s\n" % (key, value))

    def write_vectors(self, vectors):
        vectors = [v for v in vectors if v is not None]
        if sum(map(len, vectors)) != len(vectors[0]) * len(vectors):
            raise ValueError("All data vectors must be equal!")
        for i in range(len(vectors[0])):
            line = []
            for j in range(len(vectors)):
                if j != 0:
                    line.append(' ')
                line.append(str(vectors[j][i]))
            line.append("\n")
            self.stream.write(''.join(line))

    def write_scan(self, scan_header, data_vectors):
        self.write_header(scan_header)
        self.write_vectors(data_vectors)
