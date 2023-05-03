import os

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, List, Tuple, Dict, Optional, Union, Set


import click
import numpy as np

from ms_deisotope.output import ProcessedMSFileLoader, get_writer
from ms_deisotope.task import TaskBase
from ms_deisotope.data_source import ProcessedScan, ScanBunch, ProcessedRandomAccessScanSource
from ms_deisotope.output.mzml import MzMLSerializer
from ms_deisotope.tools.utils import progress


@dataclass
class PrecursorShiftBase:
    def permute(self, scan: ProcessedScan) -> ProcessedScan:
        raise NotImplementedError()

    def __call__(self, scan: ProcessedScan) -> ProcessedScan:
        return self.permute(scan)


@dataclass
class ShiftByConstant(PrecursorShiftBase):
    constant: float

    def permute(self, scan: ProcessedScan) -> ProcessedScan:
        scan = scan.copy()
        scan.precursor_information.extracted_neutral_mass += self.constant
        scan.precursor_information.mz += self.constant / scan.precursor_information.charge
        return scan


@dataclass
class ShiftByRandomBase(PrecursorShiftBase):
    seed: int
    rng: np.random.BitGenerator = field(default=None, init=False)

    def __post_init__(self):
        self.rng = np.random.default_rng(self.seed)

    def sample_offset(self) -> float:
        raise NotImplementedError()

    def permute(self, scan: ProcessedScan) -> ProcessedScan:
        scan = scan.copy()
        offset = self.sample_offset()
        scan.precursor_information.extracted_neutral_mass += offset
        scan.precursor_information.mz += offset / scan.precursor_information.charge
        return scan


class PermutePrecursorTask(TaskBase):
    inpath: Path
    outpath: Path

    permuter: PrecursorShiftBase
    reader: ProcessedRandomAccessScanSource
    writer: MzMLSerializer

    def __init__(self, inpath: Path, outpath: Path, permuter: PrecursorShiftBase):
        self.inpath = inpath
        self.outpath = outpath
        self.permuter = permuter
        self.reader = None
        self.writer = None

    def run(self):
        self.reader = ProcessedMSFileLoader(str(self.inpath))
        self.writer = get_writer(str(self.outpath), len(self.reader))
        self.writer.copy_metadata_from(self.reader)

        i = 0
        bar = progress(
            length=len(self.reader), label='Scans Processed',
            item_show_func=lambda x: f"{x.id}" if x else '-',
            color=True, fill_char=click.style('-', 'green'))
        with bar, self.writer:
            for precursor, products in self.reader.make_iterator(grouped=True):
                products = [self.permuter(product) for product in products]
                batch = ScanBunch(precursor, products)
                self.writer.save_scan_bunch(batch)
                j = len(products) + 1
                bar.update(j, precursor if i % 10 == 0 else None)
                i += 1


@click.command()
@click.option("-i", "--inpath", type=click.Path(exists=True, readable=True), help="The processed MS data file to permute")
@click.option("-o", "--outpath", type=click.Path(writable=True), default=None, help="Where to write the permuted MS data file to")
@click.option("-f", "--shift", type=float, help="How much to shift the precursor by")
def main(inpath: Path, outpath: Path, shift: float):
    task = PermutePrecursorTask(Path(inpath), Path(outpath), ShiftByConstant(shift))
    task.run()

if __name__ == "__main__":
    main.main()
