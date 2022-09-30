import sys

from ms_deisotope.output import ProcessedMzMLDeserializer
from ms_deisotope.data_source._compression import get_opener

def compare_peaks(peaks_a, peaks_b):
    missing = []
    for peak in peaks_a:
        peaks = peaks_b.all_peaks_for(
            peak.neutral_mass, 1e-6)
        for ref in peaks:
            if peak == ref:
                break
        else:
            missing.append(peak)
    return missing


def diff_deconvoluted_peak_set(peaks_a, peaks_b):
    a_missing = compare_peaks(peaks_a, peaks_b)
    b_missing = compare_peaks(peaks_b, peaks_a)
    return a_missing, b_missing


def compare_readers(result_reader, reference_reader):
    assert len(result_reader) == len(reference_reader)
    for a_bunch, b_bunch in zip(result_reader, reference_reader):
        assert len(a_bunch.products) == len(b_bunch.products)
        aprec = a_bunch.precursor
        bprec = b_bunch.precursor
        assert aprec.id == bprec.id
        diffa, diffb = diff_deconvoluted_peak_set(
            aprec.deconvoluted_peak_set, bprec.deconvoluted_peak_set)
        assert len(aprec.deconvoluted_peak_set) == len(
            bprec.deconvoluted_peak_set), "Peak Counts Diff On %r, (%r, %r)" % (aprec.id, diffa, diffb)
        assert aprec.deconvoluted_peak_set == bprec.deconvoluted_peak_set, "Peaks Diff On %r, (%r, %r)" % (
            aprec.id, diffa, diffb)

        for aprod, bprod in zip(a_bunch.products, b_bunch.products):
            assert aprod.id == bprod.id
            diffa, diffb = diff_deconvoluted_peak_set(aprod.deconvoluted_peak_set, bprod.deconvoluted_peak_set)
            assert len(aprod.deconvoluted_peak_set) == len(
                bprod.deconvoluted_peak_set), "Peak Counts Diff On %r, (%r, %r)" % (aprod.id, diffa, diffb)
            assert aprod.deconvoluted_peak_set == bprod.deconvoluted_peak_set, "Peaks Diff On %r" % (
                aprod.id, diffa, diffb)


def main(path, reference_path):
    reader = ProcessedMzMLDeserializer(get_opener(path))
    reference_reader = ProcessedMzMLDeserializer(get_opener(reference_path))
    compare_readers(reader, reference_reader)
    print("Processed Files Appear to Match Perfectly.")


if __name__ == "__main__":
    path = sys.argv[1]
    reference_path = sys.argv[2]
    main(path, reference_path)