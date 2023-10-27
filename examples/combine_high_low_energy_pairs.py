from typing import Deque, List, Dict, Tuple, Optional, DefaultDict, Any

import click

from ms_deisotope.data_source.scan.loader import ScanIterator

from ms_deisotope.peak_set import (DeconvolutedPeak, DeconvolutedPeakSet, DeconvolutedPeakSolution, ppm_error)
from ms_deisotope.output import ProcessedMSFileLoader, MzMLbSerializer
from ms_deisotope.feature_map import NeutralMassIndex
from ms_deisotope.data_source.scan import ProcessedScan, ScanBunch, PrecursorInformation


class DeconvolutedPeakMerger:
    charge_buckets: DefaultDict[int, List[DeconvolutedPeak]]

    def __init__(self, *args, **kwargs):
        self.charge_buckets = DefaultDict(list)

    def add(self, peaks: DeconvolutedPeakSet):
        for peak in peaks:
            self.charge_buckets[peak.charge].append(peak)

    def sort(self):
        for z, peaks in self.charge_buckets.items():
            peaks.sort(key=lambda x: x.neutral_mass)

    def merge_peaks_in_bucket(self, bucket: List[DeconvolutedPeak],
                              error_tolerance: float=5e-6) -> List[DeconvolutedPeak]:
        n = len(bucket)
        if not n:
            return bucket
        elif n == 1:
            return [bucket[0].clone()]
        acc = []
        current = bucket[0].clone()
        for i in range(1, n):
            pnext = bucket[i]
            if error_tolerance >= abs(ppm_error(current.neutral_mass, pnext.neutral_mass)):
                total_intensity = current.intensity + pnext.intensity
                current.neutral_mass = (current.intensity * current.neutral_mass) + (
                    pnext.neutral_mass * pnext.intensity) / total_intensity
                current.intensity = total_intensity
            else:
                acc.append(current)
                current = pnext.clone()
        acc.append(current)
        return acc

    def merge(self, error_tolerance: float=5e-6) -> DeconvolutedPeakSet:
        self.sort()
        acc = []
        for z, bucket in self.charge_buckets.items():
            acc.extend(self.merge_peaks_in_bucket(bucket, error_tolerance))
        peaks = DeconvolutedPeakSet(acc)
        peaks.reindex()
        return peaks


class PrecursorProxy:
    scan: ProcessedScan
    neutral_mass: float
    charge: int
    index: int

    def __init__(self, scan: ProcessedScan):
        self.scan = scan
        self.neutral_mass = scan.precursor_information.neutral_mass
        self.charge = scan.precursor_information.charge
        self.index = scan.index


def _index_precursors(product_scans: List[ProcessedScan]) -> NeutralMassIndex[PrecursorProxy]:
    index = NeutralMassIndex(list(map(PrecursorProxy, product_scans)))
    return index


def match_precursors(batch1: List[ProcessedScan],
                     batch2: List[ProcessedScan],
                     error_tolerance: float = 1e-6) -> Tuple[List[Tuple[ProcessedScan, ProcessedScan]],
                                                             List[ProcessedScan],
                                                             List[ProcessedScan]]:
    unmatched1 = []
    unmatched2 = []
    shared = []

    indexed = _index_precursors([p for p in batch2])
    mask = set()
    scan1: ProcessedScan
    for scan1 in batch1:
        pinfo = scan1.precursor_information
        for proxy in indexed.find_all(pinfo.neutral_mass, error_tolerance):
            if proxy.charge == pinfo.charge:
                key = proxy.index
                if key in mask:
                    continue
                else:
                    mask.add(key)
                    shared.append((scan1, proxy.scan))
                    break
        else:
            unmatched1.append(scan1)
    for proxy in indexed:
        if proxy.index not in mask:
            unmatched2.append(proxy.scan)
    return shared, unmatched1, unmatched2


class BlindScanMerger:
    scans: List[ProcessedScan]

    def __init__(self, scans=None):
        if scans is None:
            scans = []
        self.scans = []
        self.peaks = DeconvolutedPeakMerger()
        for scan in scans:
            self.add(scan)

    def add(self, scan: ProcessedScan):
        self.scans.append(scan)
        self.peaks.add(scan.deconvoluted_peak_set)

    def merge(self, error_tolerance: float=5e-6) -> Optional[ProcessedScan]:
        if not self.scans:
            return None
        if len(self.scans) == 1:
            return self.scans[0]
        scan = self.scans[-1]
        peaks = self.peaks.merge(error_tolerance)
        self.scans.sort(key=lambda x: x.tic())
        scan = scan.clone(deep=True)
        scan.deconvoluted_peak_set = peaks
        return scan


_sentinel = object()

class ScanMergingIterator:
    source: ScanIterator[Any, ProcessedScan]

    buffer_size: int
    buffers: Deque[List[ProcessedScan]]

    error_tolerance: float

    def __init__(self, source: ScanIterator, buffer_size: int=5, error_tolerance: float=5e-6) -> None:
        self.source = source
        self.error_tolerance = error_tolerance

        self.buffers = Deque([] for i in range(buffer_size))
        self.buffer_size = buffer_size

    def add_batch(self, batch: ScanBunch[ProcessedScan]):
        scans = batch.products
        paired = []
        for i, generation in enumerate(self.buffers):
            shared, unmatched1, scans = match_precursors(generation, scans)
            paired.extend(shared)
            self.buffers[i] = unmatched1
        unmatched = self.buffers.popleft()
        self.buffers.append(scans)
        return unmatched, paired

    def merge_pairs(self, pairs: List[Tuple[ProcessedScan, ProcessedScan]]) -> List[ProcessedScan]:
        scans = []
        for pair in pairs:
            scans.append(BlindScanMerger(pair).merge(self.error_tolerance))
        return scans

    def __iter__(self):
        return self

    def __next__(self):
        batch: ScanBunch[ProcessedScan] = next(self.source)
        products, paired = self.add_batch(batch)
        products.extend(self.merge_pairs(paired))
        return ScanBunch(batch.precursor, products)




@click.command()
@click.argument("inpath", type=click.Path(readable=True))
@click.argument("outpath", type=click.Path(writable=True))
def main(inpath, outpath):
    import faulthandler
    faulthandler.enable()
    reader = ProcessedMSFileLoader(inpath)
    writer = MzMLbSerializer(outpath, len(reader), include_software_entry=False)

    with writer:
        writer.copy_metadata_from(reader)
        iterator = ScanMergingIterator(reader)
        t = 0
        n = len(reader)
        for batch in iterator:
            ts = batch.precursor.scan_time
            if (ts - t) > 1:
                print(f"Processed {batch.precursor.id} {batch.precursor.index} ({batch.precursor.index / n * 100.0:0.2f}%)")
                t = ts
            writer.save_scan_bunch(batch)

if __name__ == "__main__":
    main.main()
