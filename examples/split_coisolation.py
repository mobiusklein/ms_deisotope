import sys
import os

from ms_deisotope.data_source.scan.scan import ProcessedScan
from ms_deisotope.output import MzMLSerializer, ProcessedMzMLLoader
from ms_deisotope.task.log_utils import init_logging, logger


def main(inpath: os.PathLike, outpath: os.PathLike):
    reader = ProcessedMzMLLoader(inpath)
    init_logging()
    j = 0
    with MzMLSerializer(outpath, int(len(reader) * 1.5)) as writer:
        writer.copy_metadata_from(reader)
        for bunch_i, (precursor, products) in enumerate(reader):
            if bunch_i % 5000 == 0 and bunch_i:
                logger.info("... Handled %d/%d spectra (%0.2f%%); Wrote %d spectra", precursor.index, len(reader), precursor.index / len(reader) * 100.0, j)
            writer.save(precursor)
            j += 1
            product: ProcessedScan
            for product in products:
                writer.save(product)
                j += 1
                for isolation_j, coisolation in enumerate(product.precursor_information.split_coisolations()[1:], 1):
                    dup = product.clone()
                    dup.precursor_information = coisolation
                    dup.id += f".{isolation_j}"
                    writer.save(dup)
                    j += 1


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
