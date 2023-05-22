import sys
import os

from ms_deisotope.data_source.scan.scan import ProcessedScan
from ms_deisotope.output import MzMLSerializer, ProcessedMzMLLoader
from ms_deisotope.task.log_utils import init_logging, logger


def main(inpath: os.PathLike, outpath: os.PathLike):
    reader = ProcessedMzMLLoader(inpath)
    init_logging()
    scan_count = 0
    split_count = 0
    logger.info("Reading spectra from %r", inpath)
    logger.info("Writing spectra to %r", outpath)
    with MzMLSerializer(outpath, int(len(reader) * 1.5), sample_name=reader.sample_run.name) as writer:
        writer.copy_metadata_from(reader)
        for bunch_i, (precursor, products) in enumerate(reader):
            if bunch_i % 1000 == 0 and bunch_i:
                logger.info("... Handled %d/%d spectra (%0.2f%%); Wrote %d spectra (%d split)", precursor.index,
                            len(reader), precursor.index / len(reader) * 100.0, scan_count, split_count)
            writer.save(precursor)
            scan_count += 1
            product: ProcessedScan
            for product in products:
                writer.save(product)
                scan_count += 1
                for isolation_j, coisolation in enumerate(product.precursor_information.split_coisolations()[1:], 1):
                    dup = product.clone()
                    dup.precursor_information = coisolation
                    dup.id += f".{isolation_j}"
                    writer.save(dup)
                    scan_count += 1
                    split_count += 1


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
