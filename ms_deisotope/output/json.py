from __future__ import absolute_import

import json

from ms_deisotope.data_source import ChargeNotProvided

class JSONScanFormatter(object):

    def _pack_activation(self, activation_information):
        """Pack :class:`~.ActivationInformation` into a :class:`dict` structure
        which that :class:`~psims.mzml.writer.MzMLWriter` expects.

        Parameters
        ----------
        activation_information: :class:`~.ActivationInformation`

        Returns
        -------
        :class:`dict`
        """
        params = []
        params.append({
            "name": str(activation_information.method),
        })
        if activation_information.is_multiple_dissociation():
            for method in activation_information.methods[1:]:
                params.append({"name": str(method)})
        # NOTE: Only correct for CID/HCD spectra with absolute collision energies, but that is all I have
        # to test with.
        params.append({
            "name": "collision energy",
            "value": activation_information.energy,
            "unitName": "electronvolt"
        })
        if activation_information.is_multiple_dissociation():
            energies = activation_information.energies[1:]
            supplemental_energy = None
            if activation_information.has_supplemental_dissociation():
                supplemental_energy = energies[-1]
                energies = energies[:-1]
            for energy in energies:
                params.append({
                    "name": "collision energy",
                    "value": energy,
                    "unitName": "electronvolt"
                })
            if supplemental_energy is not None:
                params.append({
                    "name": 'supplemental collision energy',
                    "value": supplemental_energy,
                    "unitName": "electronvolt"
                })

        for key, val in activation_information.data.items():
            arg = {
                "name": key,
                "value": val
            }
            try:
                arg['unitName'] = val.unit_info
            except AttributeError:
                pass
            params.append(arg)
        return params

    def _pack_precursor_information(self, precursor_information, activation_information=None,
                                    isolation_window=None):
        """Repackage the :class:`~.PrecursorInformation`, :class:`~.ActivationInformation`,
        and :class:~.IsolationWindow` into the nested :class:`dict` structure that
        :class:`~psims.mzml.writer.MzMLWriter` expects.

        Parameters
        ----------
        precursor_information : :class:`~.PrecursorInformation`
        activation_information : :class:`~.ActivationInformation`, optional
        isolation_window : :class:`~.IsolationWindow`, optional

        Returns
        -------
        :class:`dict`
        """
        package = {}
        # If the scan bunch has been fully deconvoluted and it's PrecursorInformation
        # filled in, its extracted fields will be populated and should be used, otherwise
        # use the default read values.
        if precursor_information is not None:
            extracted_neutral_mass = precursor_information.extracted_neutral_mass
            if (extracted_neutral_mass != 0):
                package = {
                    "mz": precursor_information.extracted_mz,
                    "intensity": precursor_information.extracted_intensity,
                    "charge": precursor_information.extracted_charge,
                    "scan_id": precursor_information.precursor_scan_id,
                    "params": [
                        {"ms_deisotope:defaulted": precursor_information.defaulted},
                        {"ms_deisotope:orphan": precursor_information.orphan}
                    ]
                }
                if precursor_information.coisolation:
                    for p in precursor_information.coisolation:
                        package['params'].append({
                            "name": "ms_deisotope:coisolation",
                            "value": "%f %f %d" % (p.neutral_mass, p.intensity, p.charge)
                        })
            else:
                package = {
                    "mz": precursor_information.mz,
                    "intensity": precursor_information.intensity,
                    "charge": precursor_information.charge,
                    "scan_id": precursor_information.precursor_scan_id
                }
        else:
            package['mz'] = None
            package["charge"] = None
        if package['charge'] == ChargeNotProvided:
            package["charge"] = None
        if activation_information is not None:
            package['activation'] = self._pack_activation(
                activation_information)
        if isolation_window is not None:
            package['isolation_window'] = {
                "lower": isolation_window.lower,
                "target": isolation_window.target,
                "upper": isolation_window.upper
            }
        return package

    def peak_set_to_json(self, peak_set):
        points = []
        for peak in peak_set:
            points.append({
                "mz": peak.mz,
                "intensity": peak.intensity
            })
        return points

    def deconvoluted_peak_set_to_json(self, peak_set):
        points = []
        for peak in peak_set:
            points.append({
                "mz": peak.mz,
                "neutral_mass": peak.neutral_mass,
                "charge": peak.charge,
                "envelope": [
                    {'mz': p.mz, 'intensity': p.intensity} for p in peak.envelope
                ],
                "intensity": peak.intensity,
            })
        return points

    def raw_data_arrays_to_json(self, arrays):
        data = {
            "mz": arrays.mz.tolist(),
            "intensity": arrays.intensity.tolist()
        }
        return data

    def format_scan(self, scan):
        data = {}
        data['ms_level'] = scan.ms_level
        data['is_profile'] = scan.is_profile
        data['polarity'] = scan.polarity
        data['scan_time'] = scan.scan_time
        data['id'] = scan.id
        if scan.is_profile:
            data['arrays'] = self.raw_data_arrays_to_json(scan.arrays)
        if scan.peak_set is not None:
            data['peak_set'] = self.peak_set_to_json(scan.peak_set)
        if scan.deconvoluted_peak_set is not None:
            data['deconvoluted_peak_set'] = self.deconvoluted_peak_set_to_json(scan.deconvoluted_peak_set)
        if scan.precursor_information is not None:
            pinfo = data['precursor_information'] = self._pack_precursor_information(
                scan.precursor_information, scan.activation, scan.isolation_window)
            data['activation'] = pinfo.pop("activation", None)
            data['isolation_window'] = pinfo.pop("isolation_window")
        return data

    def dump(self, scan, fh, **kwargs):
        data = self.format_scan(scan)
        return json.dump(data, fh, **kwargs)

    def dumps(self, scan, **kwargs):
        data = self.format_scan(scan)
        return json.dumps(data, **kwargs)


formatter = JSONScanFormatter()


dump = formatter.dump
dumps = formatter.dumps
format = formatter.format_scan
