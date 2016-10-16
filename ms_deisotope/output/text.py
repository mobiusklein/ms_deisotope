from .common import ScanSerializerBase


class TextScanSerializerBase(ScanSerializerBase):
    def __init__(self, *args, **kwargs):
        pass

    def save_scan_bunch(self, bunch):
        self.write_scan(*self.prepare_scan_data(bunch.precursor))
        for scan in bunch.products:
            self.write_scan(*self.prepare_scan_data(scan))

    def construct_header(self, scan):
        prec_info = scan.precursor_information
        header_dict = {}
        if prec_info is not None:
            header_dict['precursor_neutral_mass'] = prec_info.extracted_neutral_mass
            header_dict['precursor_charge'] = prec_info.extracted_charge
            header_dict['precursor_intensity'] = prec_info.extracted_intensity
            header_dict['precursor_scan_id'] = prec_info.precursor_scan_id
        header_dict['scan_time'] = scan.scan_time
        header_dict['id'] = scan.id
        header_dict['title'] = scan.title
        header_dict['ms_level'] = scan.ms_level
        return header_dict

    def format_peak_vectors(self, scan):
        mz = [p.mz for p in scan.deconvoluted_peak_set]
        neutral_mass = [p.neutral_mass for p in scan.deconvoluted_peak_set]
        intensity = [p.intensity for p in scan.deconvoluted_peak_set]
        charge = [p.charge for p in scan.deconvoluted_peak_set]
        return (mz, neutral_mass, intensity, charge)

    def prepare_scan_data(self, scan):
        header_dict = self.construct_header(scan)
        vectors = self.format_peak_vectors(scan)
        return (header_dict, vectors)

    def write_scan(self, scan_header, data_vectors):
        raise NotImplementedError()


class HeaderedDelimitedWriter(TextScanSerializerBase):
    def __init__(self, stream):
        self.stream = stream

    def write_header(self, header_dict):
        for key, value in header_dict.items():
            self.stream.write("#%s=%s\n" % (key, value))

    def write_vectors(self, vectors):
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
        self.stream.write('\n')

    def write_scan(self, scan_header, data_vectors):
        self.write_header(scan_header)
        self.write_vectors(data_vectors)
