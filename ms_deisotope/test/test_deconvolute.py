import unittest
import logging
import numpy as np

import brainpy

from ms_deisotope.data_source import common, mzml, MSFileLoader
from ms_deisotope.averagine import peptide, glycopeptide, TheoreticalIsotopicPattern, AveragineCache
from ms_peak_picker import reprofile
from ms_deisotope.deconvolution import (
    deconvolute_peaks, AveragineDeconvoluter,
    AveraginePeakDependenceGraphDeconvoluter,
    CompositionListDeconvoluter,
    CompositionListPeakDependenceGraphDeconvoluter,
    count_placeholders, drop_placeholders,
    ChargeIterator)
from ms_deisotope.scoring import PenalizedMSDeconVFitter, MSDeconVFitter
from brainpy import neutral_mass

from ms_deisotope.test.test_scan import (
    make_profile, points, fwhm,
    FittedPeak)
from ms_deisotope.test.common import datafile

import pickle


logger = logging.getLogger('ms_deisotope.test.test_deconvolution')


class TestAveragineDeconvolution(unittest.TestCase):
    def make_scan(self):
        mz, intensity = make_profile(points, fwhm)
        scan = common.Scan(
            {
                "id": "test-scan",
                "index": 0,
                "m/z array": mz,
                "ms level": 1,
                "scan time": 0.0,
                "intensity array": intensity,
                "profile spectrum": "",
                "positive scan": "",
            },
            mzml.MzMLDataInterface())
        scan.pick_peaks()
        for point in points:
            self.assertIsNotNone(scan.peak_set.has_peak(point[0]))
        return scan

    def test_quick_charge(self):
        scan = self.make_scan()
        peaks = scan.peak_set
        peak = peaks[0]
        charge_states = ChargeIterator(1, 8)
        charge_states.sequence_from_quickcharge(peaks, peak)
        states = list(charge_states)
        self.assertEqual(states, [1, 3])

    def test_deconvolution(self):
        scan = self.make_scan()
        algorithm_type = AveragineDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, {
                "averagine": peptide,
                "scorer": PenalizedMSDeconVFitter(5., 1.),
                "use_subtraction": False,
            }, left_search_limit=3, deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        assert len(dpeaks) == 6
        for point in points:
            peak = dpeaks.has_peak(neutral_mass(point[0], point[1]))
            self.assertIsNotNone(peak)

    def test_graph_deconvolution(self):
        scan = self.make_scan()
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        algorithm_type = AveraginePeakDependenceGraphDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, {
                "averagine": peptide,
                "scorer": PenalizedMSDeconVFitter(5., 1.)
            }, deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        assert len(dpeaks) == 2
        deconvoluter = deconresult.deconvoluter
        for point in points:
            peak = dpeaks.has_peak(neutral_mass(point[0], point[1]))
            self.assertIsNotNone(peak)
            fp = scan.has_peak(peak.mz)
            self.assertAlmostEqual(
                deconvoluter.peak_dependency_network.find_solution_for(fp).mz,
                peak.mz, 3)


class TestCompositionListDeconvolution(unittest.TestCase):
    compositions = [
        brainpy.parse_formula('C84H138N6O62'),
        brainpy.parse_formula('C90H148N6O66'),
        brainpy.parse_formula('C90H148N6O69S1'),
        brainpy.parse_formula('C62H94N4O69S6'),
    ]

    charges = [
        ((2, 100), (3, 250),),
        ((2, 75), (3, 200),),
        ((2, 90), (3, 270),),
        ((3, 20), (4, 90), (5, 45))
    ]

    def make_tids(self):
        tids = []
        ions = []
        for comp, charges in zip(self.compositions, self.charges):
            for charge, abundance in charges:
                tid = brainpy.isotopic_variants(comp, charge=-charge)
                tid = TheoreticalIsotopicPattern(tid, tid[0].mz)
                tid.scale_raw(abundance * 100)
                tids.append(tid)
                ions.append((comp, -charge))
        return tids, ions

    @staticmethod
    def get_nearest_index(query_mz, tid_list):
        best_index = None
        best_error = float('inf')

        for i, tid in enumerate(tid_list):
            error = abs(tid.monoisotopic_mz - query_mz)
            if error < best_error:
                best_error = error
                best_index = i
        return best_index

    def make_scan(self):
        peaks = []
        tids, ions = self.make_tids()
        list(map(peaks.extend, tids))
        peaks.sort(key=lambda x: x.mz)

        mz = np.array([0])
        intensity = np.array([0])

        fpeaks = []
        for p in peaks:
            fpeaks.append(FittedPeak(
                mz=p.mz, intensity=p.intensity, signal_to_noise=p.intensity, peak_count=-1, index=-1,
                full_width_at_half_max=fwhm, area=p.intensity))
        mz, intensity = reprofile(fpeaks)
        scan = common.Scan(
            {
                "id": "test-scan",
                "index": 0,
                "m/z array": mz,
                "intensity array": intensity,
                "ms level": 1,
                "scan time": 0.0,
                "profile spectrum": "",
                "negative scan": "",
            },
            mzml.MzMLDataInterface())
        return scan

    def test_deconvolution(self):
        scan = self.make_scan()
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        algorithm_type = CompositionListDeconvoluter
        decon_config = {
            "composition_list": self.compositions,
            "scorer": PenalizedMSDeconVFitter(5., 2.),
            # "scorer": MSDeconVFitter(5.),
            "use_subtraction": False
        }
        deconresult = deconvolute_peaks(
            scan.peak_set, decon_config, charge_range=(-1, -8), deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        n_cases = sum(map(len, self.charges))
        if not (len(dpeaks) == n_cases):
            tids, ions = self.make_tids()
            tids, ions = zip(*sorted(zip(tids, ions), key=lambda x: x[0].monoisotopic_mz))
            seen = set()
            for i, dp in enumerate(sorted(dpeaks, key=lambda x: x.mz)):
                ix = self.get_nearest_index(dp.mz, tids)
                logger.warning("%0.3f %d %0.3f %r (Matched %d)", dp.neutral_mass, dp.charge, dp.score, dp.solution, ix)
                seen.add(ix)
            indices = set(range(len(ions)))
            missed = list(indices - seen)
            deconvoluter = algorithm_type(scan.peak_set.clone(), **decon_config)
            for ix in missed:
                tid = deconvoluter.generate_theoretical_isotopic_cluster(*ions[ix])
                assert np.isclose(sum(p.intensity for p in tid), 1.0)
                monoisotopic_peak = deconvoluter.peaklist.has_peak(tid[0].mz, 2e-5)
                if monoisotopic_peak is not None:
                    tid = deconvoluter.recalibrate_theoretical_mz(tid, monoisotopic_peak.mz)
                eid = deconvoluter.match_theoretical_isotopic_distribution(
                    tid.peaklist, 2e-5)
                missed_peaks = count_placeholders(eid)
                deconvoluter.scale_theoretical_distribution(tid, eid)
                score = deconvoluter.scorer.evaluate(deconvoluter.peaklist, eid, tid.peaklist)
                fit_record = deconvoluter.fit_composition_at_charge(*ions[ix])
                eid = fit_record.experimental
                tid = fit_record.theoretical
                rep_eid = drop_placeholders(eid)
                validation = (len(rep_eid) < 2), (len(rep_eid) < len(tid) / 2.), (
                    len(rep_eid) == 1 and fit_record.charge > 1)
                composition, charge = ions[ix]
                logger.warning("Missed %r %d (%d missed peaks, score = %0.3f, record = %r, validation = %r)" % (
                    composition, charge, missed_peaks, score, fit_record, validation))
            assert not missed

    def test_graph_deconvolution(self):
        scan = self.make_scan()
        scan.pick_peaks()
        self.assertIsNotNone(scan.peak_set)
        algorithm_type = CompositionListPeakDependenceGraphDeconvoluter
        decon_config = {
            "composition_list": self.compositions,
            "scorer": PenalizedMSDeconVFitter(5., 2.),
            "use_subtraction": True
        }
        deconresult = deconvolute_peaks(
            scan.peak_set, decon_config, charge_range=(-1, -8), deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        n_cases = sum(map(len, self.charges))
        # assert len(dpeaks) == n_cases
        if not (len(dpeaks) == n_cases):
            tids, ions = self.make_tids()
            tids, ions = zip(*sorted(zip(tids, ions), key=lambda x: x[0].monoisotopic_mz))
            seen = set()
            for i, dp in enumerate(sorted(dpeaks, key=lambda x: x.mz)):
                ix = self.get_nearest_index(dp.mz, tids)
                logger.warning("%0.3f %d %0.3f %r (Matched %d)", dp.neutral_mass, dp.charge, dp.score, dp.solution, ix)
                seen.add(ix)
            indices = set(range(len(ions)))
            missed = list(indices - seen)
            deconvoluter = algorithm_type(scan.peak_set.clone(), **decon_config)
            for ix in missed:
                tid = deconvoluter.generate_theoretical_isotopic_cluster(*ions[ix])
                assert np.isclose(sum(p.intensity for p in tid), 1.0)
                monoisotopic_peak = deconvoluter.peaklist.has_peak(tid[0].mz, 2e-5)
                if monoisotopic_peak is not None:
                    tid = deconvoluter.recalibrate_theoretical_mz(tid, monoisotopic_peak.mz)
                eid = deconvoluter.match_theoretical_isotopic_distribution(
                    tid.peaklist, 2e-5)
                missed_peaks = count_placeholders(eid)
                deconvoluter.scale_theoretical_distribution(tid, eid)
                score = deconvoluter.scorer.evaluate(deconvoluter.peaklist, eid, tid.peaklist)
                fit_record = deconvoluter.fit_composition_at_charge(*ions[ix])
                eid = fit_record.experimental
                tid = fit_record.theoretical
                rep_eid = drop_placeholders(eid)
                validation = (len(rep_eid) < 2), (len(rep_eid) < len(tid) / 2.), (
                    len(rep_eid) == 1 and fit_record.charge > 1)
                composition, charge = ions[ix]
                logger.warning("Missed %r %d (%d missed peaks, score = %0.3f, record = %r, validation = %r)" % (
                    composition, charge, missed_peaks, score, fit_record, validation))
            assert not missed


class TestSolutionRetrieval(unittest.TestCase):
    def make_scan(self):
        complex_compressed_mzml = datafile("20150710_3um_AGP_001_29_30.mzML.gz")
        reader = MSFileLoader(complex_compressed_mzml)
        bunch = next(reader)
        return bunch

    def target_envelopes(self):
        envelopes = [
            [1161.00927734375,
             1161.25830078125,
             1161.508544921875,
             1161.7586669921875,
             1162.0111083984375,
             1162.2607421875,
             1162.5128173828125],
            [1197.52197265625,
             1197.7735595703125,
             1198.024169921875,
             1198.2742919921875,
             1198.524658203125,
             1198.774169921875,
             1199.02490234375],
            [929.0088500976562,
             929.2091674804688,
             929.4107055664062,
             929.6107788085938,
             929.8109741210938,
             930.01171875,
             930.213134765625],
            [1088.234619140625,
             1088.4853515625,
             1088.7354736328125,
             1088.987060546875,
             1089.237060546875,
             1089.4873046875,
             1089.7379150390625],
            [1002.0338745117188,
             1002.2335815429688,
             1002.4354858398438,
             1002.63671875,
             1002.837646484375,
             1003.0364990234375,
             1003.2337646484375],
        ]
        return envelopes

    def test_retrieve_deconvolution_solution(self):
        bunch = self.make_scan()
        scan = bunch.precursor
        scan.pick_peaks()
        ms1_deconvolution_args = {
            "averagine": glycopeptide,
            "scorer": PenalizedMSDeconVFitter(20., 2.),
        }
        priorities = []
        for product in bunch.products:
            priorities.append(scan.has_peak(product.precursor_information.mz))
        algorithm_type = AveraginePeakDependenceGraphDeconvoluter
        deconresult = deconvolute_peaks(
            scan.peak_set, ms1_deconvolution_args,
            priority_list=priorities,
            deconvoluter_type=algorithm_type)
        dpeaks = deconresult.peak_set
        deconvoluter = deconresult.deconvoluter
        priority_results = deconresult.priorities
        reference_deconvoluter = algorithm_type(scan.peak_set.clone(), **ms1_deconvolution_args)
        for i, result in enumerate(priority_results):
            query = priorities[i].mz
            if result is None:
                logger.warn("Query %d (%f) had no result", i, query)
                raw_peaks = scan.peak_set.between(query - 2, query + 3)
                anchor_peak = scan.peak_set.has_peak(query)
                deconvoluted_peaks = dpeaks.between(query - 2, query + 3, use_mz=True)
                possible_solutions = reference_deconvoluter._fit_all_charge_states(anchor_peak)
                sols = []
                logger.warn("Possible Solutions %r", possible_solutions)
                if not possible_solutions:
                    for charge in [3, 4, 5]:
                        tid = reference_deconvoluter.averagine.isotopic_cluster(
                            anchor_peak.mz, charge)
                        assert np.isclose(tid.monoisotopic_mz, anchor_peak.mz)
                        assert np.isclose(sum(p.intensity for p in tid), 1.0)
                        eid = reference_deconvoluter.match_theoretical_isotopic_distribution(
                            tid.peaklist, error_tolerance=2e-5)
                        assert len(eid) == len(tid)
                        record = reference_deconvoluter._evaluate_theoretical_distribution(
                            eid, tid, anchor_peak, charge)
                        sols.append(record)
                    logger.warn("Manually Generated Solutions %r", sols)
                assert anchor_peak is not None and raw_peaks and (possible_solutions or sols) and not deconvoluted_peaks
                assert deconvoluter.peak_dependency_network.find_solution_for(anchor_peak) is not None
                assert dpeaks.has_peak(query, use_mz=True)
                # error out
                assert result is not None
            else:
                assert 0 <= abs(result.mz - query) < 1
                anchor_peak = scan.peak_set.has_peak(query)
                assert deconvoluter.peak_dependency_network.find_solution_for(anchor_peak) is not None


class TestIncrementalGraphExtraction(unittest.TestCase):
    @staticmethod
    def make_scan():
        reader = MSFileLoader(datafile("20150710_3um_AGP_001_29_30.mzML.gz"))
        scan = reader.get_scan_by_id("scanId=1740086")
        return scan

    @classmethod
    def setUpClass(cls):
        cls.scan = cls.make_scan()
        cls.scan.pick_peaks()

    @classmethod
    def tearDownClass(cls):
        cls.scan.clear()
        cls.scan.source.close()

    def build_deconvoluter(self, scan, averagine, **kwargs):
        deconvoluter = AveraginePeakDependenceGraphDeconvoluter(
            scan.peak_set, averagine=averagine, scorer=MSDeconVFitter(10),
            **kwargs)
        peak = deconvoluter.has_peak(138.5309, 2e-5)
        return peak, deconvoluter

    def test_extraction(self):
        scan = self.scan

        peak, deconvoluter = self.build_deconvoluter(scan, peptide)
        deconvoluter._deconvolution_step(0, truncate_after=0.8, charge_range=(1, 4))

        with open(datafile("extraction_base_averagine.pkl"), 'rb') as fh:
            reference_averagine = pickle.load(fh)

        diff = set(deconvoluter.averagine.backend) - set(reference_averagine.backend)
        assert len(diff) == 0
        assert len(deconvoluter.averagine.backend) == 8868
        assert reference_averagine == deconvoluter.averagine

        cluster = deconvoluter.peak_dependency_network.find_cluster_for(peak)
        spanned = cluster.fits_using_mz(peak.mz)
        assert len(cluster) == 4
        assert len(spanned) == 2
        assert np.isclose(cluster.best_fit.monoisotopic_peak.mz, 138.19520)
        assert cluster.best_fit.charge == 3

    def test_extraction_cached_averagine(self):
        scan = self.scan
        cache = AveragineCache(peptide)
        cache.populate(truncate_after=0.8)
        peak2, deconvoluter = self.build_deconvoluter(scan, cache)
        deconvoluter._deconvolution_step(
            0, truncate_after=0.8, charge_range=(1, 4))

        with open(datafile("extraction_cached_averagine.pkl"), 'rb') as fh:
            reference_averagine = pickle.load(fh)

        diff = set(deconvoluter.averagine.backend) - set(reference_averagine.backend)
        assert len(diff) == 0
        assert len(deconvoluter.averagine.backend) == 23960
        assert reference_averagine == deconvoluter.averagine

        cluster2 = deconvoluter.peak_dependency_network.find_cluster_for(peak2)
        spanned2 = cluster2.fits_using_mz(peak2.mz)
        assert len(cluster2) == 4
        assert len(spanned2) == 2
        assert np.isclose(cluster2.best_fit.monoisotopic_peak.mz, 138.19520)
        assert cluster2.best_fit.charge == 3

    def test_extraction_quick_charge(self):
        scan = self.scan

        peak3, deconvoluter = self.build_deconvoluter(scan, peptide, use_quick_charge=True)
        deconvoluter._deconvolution_step(0, truncate_after=0.8, charge_range=(1, 4))

        with open(datafile("extraction_quick_charge_averagine.pkl"), 'rb') as fh:
            reference_averagine = pickle.load(fh)

        diff = set(deconvoluter.averagine.backend) - set(reference_averagine.backend)
        assert len(diff) == 0
        assert len(deconvoluter.averagine.backend) == 5136
        assert reference_averagine == deconvoluter.averagine

        cluster3 = deconvoluter.peak_dependency_network.find_cluster_for(peak3)
        spanned3 = cluster3.fits_using_mz(peak3.mz)
        assert len(cluster3) == 1
        assert len(spanned3) == 1
        assert np.isclose(cluster3.best_fit.monoisotopic_peak.mz, 138.53090)
        assert cluster3.best_fit.charge == 4


if __name__ == '__main__':
    unittest.main()
