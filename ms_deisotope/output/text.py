from .common import ScanSerializerBase


class _HeaderInformation(dict):
    """A simple key-value store that also carries a reference to a :class:`~.Scan` object.

    Intended to be used internally by the :class:`TextScanSerializerBase` class hierarchy,
    but should not be needed elsewhere either.

    """
    def __init__(self, scan):
        dict.__init__(self)
        self.scan = scan


class TextScanSerializerBase(ScanSerializerBase):
    """A base class for text file-based scan serialization formats like
    simple custom tabular formats, MS1, MS2, or MGF.

    This class tries to flatten out the hierarchical structure of a :class:`~.Scan`
    into a dictionary that is more convenient to project onto a text file.

    Attributes
    ----------
    stream: writable file-like
        The stream to write content to. Should support basic IO operations like
        write, flush, and close.
    deconvoluted: bool
        Whether the file is storing spectra in deconvoluted mode or not.

    """
    def __init__(self, stream, deconvoluted=True):
        super(TextScanSerializerBase, self).__init__(stream, deconvoluted=deconvoluted)
        self.stream = stream
        self.deconvoluted = deconvoluted

    def close(self):
        """Close :attr:`stream`.
        """
        self.stream.close()

    def save_scan(self, scan, **kwargs):
        self.write_scan(*self.prepare_scan_data(scan))

    def construct_header(self, scan):
        """Build a dictionary of properties describing `scan` specifically not
        including the peak list or signal information.

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The scan to describe

        Returns
        -------
        :class:`_HeaderInformation`
        """
        prec_info = scan.precursor_information
        header_dict = _HeaderInformation(scan)
        if prec_info is not None:
            if prec_info.extracted_charge:
                header_dict['precursor_neutral_mass'] = prec_info.extracted_neutral_mass
                header_dict['precursor_mz'] = prec_info.extracted_mz
                header_dict['precursor_charge'] = prec_info.extracted_charge
                header_dict['precursor_intensity'] = prec_info.extracted_intensity
            else:
                header_dict['precursor_neutral_mass'] = prec_info.neutral_mass
                header_dict['precursor_mz'] = prec_info.mz
                header_dict['precursor_charge'] = prec_info.charge
                header_dict['precursor_intensity'] = prec_info.intensity
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
        """Build a simple peak attribute list format commonly used by most
        text formats.

        Most text based formats represent peaks as a list of delimited values
        per line. This method unpacks peak lists into parallel lists for each
        peak attribute.

        This method may behave differently if :attr:`deconvoluted` is true or
        not.

        Parameters
        ----------
        scan : :class:`~.ScanBase`
            The scan to extract peak attributes for.

        Returns
        -------
        mz: list of float
            The m/z of each peak. If :attr:`deconvoluted` is true, this may be
            the peak's neutral mass.
        intensity: list of float
            The intensity of each peak
        charge: list of int, optional
            If :attr:`deconvoluted` is true, this may be a list of each peak's charge
            state. :const:`None` otherwise.
        """
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
        """Prepare data from `scan` for formatting as text.

        Parameters
        ----------
        scan: :class:`~.ScanBase`
            The scan to prepare

        Returns
        -------
        header: :class:`_HeaderInformation`
            The scan header metadata
        vectors: :class:`tuple` of :class:`list`
            The peak attribute lists for all peaks in `scan`
        """
        header_dict = self.construct_header(scan)
        vectors = self.format_peak_vectors(scan)
        return (header_dict, vectors)

    def write_scan(self, scan_header, data_vectors):
        """Format the given scan data and write it to :attr:`stream`

        Should receive its parameter information from :meth:`prepare_scan_data`

        Parameters
        ----------
        scan_header : :class:`_HeaderInformation`
            The scan header metadata
        data_vectors : :class:`Sequence`
            The peak attribute lists

        """
        raise NotImplementedError()


class HeaderedDelimitedWriter(TextScanSerializerBase):
    """A simple but usable implementation of :class:`TextScanSerializerBase`

    Writes all metadata headers as "#{key}={value}\\n" immediately preceding
    the peak attribute list.

    """
    def __init__(self, stream, deconvoluted=True):
        super(HeaderedDelimitedWriter, self).__init__(stream, deconvoluted)

    def write_header(self, header_dict):
        """Write the header information in the form "#{key}={value}" with one
        key-value pair per line.

        Parameters
        ----------
        header_dict : :class:`_HeaderInformation`
        """
        for key, value in header_dict.items():
            self.stream.write(("#%s=%s\n" % (key, value)).encode("utf8"))

    def write_vectors(self, vectors):
        """Write out the peak attribute lists along peak index.

        Parameters
        ----------
        vectors : :class:`Sequence`
            The peak attribute lists to write.
        """
        vectors = [v for v in vectors if v is not None]
        if sum(map(len, vectors)) != len(vectors[0]) * len(vectors):
            raise ValueError("All data vectors must be equal!")
        for i in range(len(vectors[0])):
            line = []
            for j, vector in enumerate(vectors):
                if j != 0:
                    line.append(' ')
                line.append(str(vector[i]))
            line.append("\n")
            self.stream.write(''.join(line).encode('utf8'))

    def write_scan(self, scan_header, data_vectors):
        self.write_header(scan_header)
        self.write_vectors(data_vectors)
