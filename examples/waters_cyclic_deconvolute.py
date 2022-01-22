import os
import logging
import sys
import psutil

from itertools import chain

import numpy as np

import ms_deisotope
from ms_deisotope.peak_dependency_network.intervals import Interval, IntervalTreeNode
from ms_deisotope.peak_set import IonMobilityDeconvolutedPeak, DeconvolutedPeakSet
from ms_deisotope.data_source.scan import ProcessedScan, PrecursorInformation, ScanBunch
from ms_deisotope.data_source.scan.scan_iterator import MSEIterator
from ms_deisotope.output.mzml import MzMLSerializer, IonMobilityAware3DMzMLSerializer, ProcessedGeneric3DIonMobilityFrameSource
from ms_deisotope.data_source.scan.mobility_frame import IonMobilityFrame, Generic3DIonMobilityFrameSource

from ms_deisotope.tools.utils import register_debug_hook


import click


logger = logging.getLogger("mse_deconvolute")
logger.addHandler(logging.NullHandler())


def make_iterator(reader, start_index, stop_index=float('inf'), low_energy_function=1, lock_mass_function=3):
    iterator = MSEIterator(
        reader.start_from_frame(
            index=start_index, require_ms1=False, grouped=False),
        lambda x: x, low_energy_function, lock_mass_function)
    for bunch in iterator:
        if bunch.precursor:
            i = bunch.precursor.index
        elif bunch.products:
            i = bunch.products[0].index
        else:
            i = 0
        if i >= stop_index:
            break
        yield bunch


def connect_products_to_precursors(parent, child):
    cft = IntervalTreeNode.build(
        [Interval(f.start_time, f.end_time, [f]) for f in child.deconvoluted_features])

    prec_prod_rels = []
    for f in parent.deconvoluted_features:
        if len(f) < 3:
            continue
        prods = []
        for prod in (chain.from_iterable(cft.overlaps(f.start_time, f.end_time))):
            if prod.charge > f.charge:
                continue
            if prod.neutral_mass > f.neutral_mass:
                continue
            chunk = prod._copy_chunk(
                prod[prod.find_time(f.start_time)[1]:
                     prod.find_time(f.end_time)[1] + 1])
            prods.append(chunk)
        prods.sort(key=lambda x: x.neutral_mass)
        prec_prod_rels.append((f, prods))
    prec_prod_rels.sort(key=lambda x: (
        x[0].start_time, x[0].neutral_mass, len(x[1])))
    return prec_prod_rels


def weighted_centroid(feature):
    total = 0
    normalizer = 0
    for node in feature:
        weight = node.total_intensity()
        total += node.time * weight
        normalizer += weight
    return total / normalizer


def merge_envelopes(envelopes):
    base = envelopes[0].clone()
    for env in envelopes[1:]:
        for i, p in enumerate(env):
            base[i].intensity += p.intensity

    return base


def feature_to_peak(feature):
    peak_cluster = feature.peaks
    peak_cluster = [pi for p in peak_cluster for pi in p]
    total_intensity = sum(p.intensity for p in peak_cluster)
    mz = sum(p.mz * p.intensity for p in peak_cluster) / total_intensity
    neutral_mass = sum(
        p.neutral_mass * p.intensity for p in peak_cluster) / total_intensity
    most_abundant_mass = sum(
        p.most_abundant_mass * p.intensity for p in peak_cluster) / total_intensity
    a_to_a2_ratio = sum(
        p.a_to_a2_ratio * p.intensity for p in peak_cluster) / total_intensity
    average_mass = sum(
        p.average_mass * p.intensity for p in peak_cluster) / total_intensity
    signal_to_noise = sum(p.signal_to_noise *
                          p.intensity for p in peak_cluster) / total_intensity
    fwhm = sum(p.full_width_at_half_max *
               p.intensity for p in peak_cluster) / total_intensity
    area = sum(p.area * p.intensity for p in peak_cluster) / total_intensity
    score = sum(p.score * p.intensity for p in peak_cluster) / total_intensity
    charge = peak_cluster[0].charge
    envelope = merge_envelopes([p.envelope for p in peak_cluster])
    return IonMobilityDeconvolutedPeak(
        neutral_mass=neutral_mass, intensity=total_intensity, charge=charge, signal_to_noise=signal_to_noise,
        full_width_at_half_max=fwhm, index=-1, a_to_a2_ratio=a_to_a2_ratio, most_abundant_mass=most_abundant_mass,
        average_mass=average_mass, score=score, envelope=envelope, mz=mz, fit=None, chosen_for_msms=False,
        area=area, drift_time=weighted_centroid(feature))


def frame_pair_to_scan_bunch(parent, child, prec_prod_rels, scan_index):
    precursor = ProcessedScan(
        parent.id, parent.id, None, 1, parent.time, scan_index, None,
        DeconvolutedPeakSet([feature_to_peak(f)
                             for f in parent.deconvoluted_features]).reindex(),
        # Hard-coded positive mode for the moment while IonMobilityFrame-Source-like does not have a polarity
        1,
    )

    scan_index += 1

    pseudo_ms2s = []

    for prec_feat, prod_feats in prec_prod_rels:
        if abs(prec_feat.charge) == 1:
            continue
        pinfo = PrecursorInformation(
            prec_feat.mz, prec_feat.total_signal, prec_feat.charge, precursor.id, None, annotations={
                "ion mobility drift time": weighted_centroid(prec_feat),
                "source frame id": child.id,
            })
        assert pinfo.has_ion_mobility()
        # Eventually, create an artificial ScanAcquisitionInformation instance here to let us set the
        # drift time of the product scan to match the precursor drift time too?
        prod = ProcessedScan(
            'merged=%d' % scan_index, child.id + '.%d' % scan_index, pinfo, 2, child.time,
            scan_index, None,
            DeconvolutedPeakSet([feature_to_peak(f)
                                 for f in prod_feats]).reindex(),
            # Hard-coded positive mode for the moment while IonMobilityFrame-Source-like does not have a polarity
            1, activation=child.activation)
        prod.annotations['spectrum title'] = child.id + '.%d' % scan_index
        pseudo_ms2s.append(prod)
        scan_index += 1

    return ScanBunch(precursor, pseudo_ms2s), scan_index


class IonMobilityAwareMzMLSerializer(MzMLSerializer):
    def _prepare_extra_arrays(self, scan, **kwargs):
        extra_arrays = super(IonMobilityAwareMzMLSerializer,
                             self)._prepare_extra_arrays(scan, **kwargs)
        if scan.deconvoluted_peak_set is not None:
            # This is sensitive to units used? Shouldn't there be a unit key?
            # Waters uses milliseconds
            extra_arrays.append(("mean drift time array", [
                p.drift_time for p in scan.deconvoluted_peak_set
            ]))
        return extra_arrays


averagine_map = {
    "glycopeptide": ms_deisotope.glycopeptide,
    "heparin": ms_deisotope.heparin,
    "peptide": ms_deisotope.peptide,
    "glycan": ms_deisotope.glycan,
    "heparan_sulfate": ms_deisotope.heparan_sulfate,
    "permethylated_glycan": ms_deisotope.permethylated_glycan,
}


@click.command("cyclic_deconvolute")
@click.argument("input_path", type=click.Path())
@click.argument("output_path", type=click.Path(writable=True))
@click.option("-m", "--lockmass-config", type=float, help="The lock mass used", default=785.8421)
@click.option("-s", "--start-time", type=float, help="The time to start processing cycles from", default=0)
@click.option("-e", "--end-time", type=float, help="The time to stop processing cycles at", default=None)
@click.option("-a", "--averagine", type=click.Choice(list(averagine_map)), default='glycopeptide',
              help='The isotopic model to use. Defaults to the glycopeptide averagine.')
@click.option("-i", "--minimum-intensity", type=float, default=10.0, help="The minimum intensity to accept a peak")
@click.option("-n", "--no-product-splitting", is_flag=True, default=False,
              help=(
                  'If true, do not split high energy frames into pseudo-MS2 spectra,'
                  ' creating continuous ion mobility features'))
@click.option("-k", "--lockmass-function", type=int, default=3, help="The number of the lock mass function. For normal low-high MSE this is 3.")
def main(input_path, output_path, lockmass_config, start_time=0, end_time=None, averagine='glycopeptide',
         minimum_intensity=10.0, no_product_splitting=False, lockmass_function=3):
    logging.basicConfig(
        level="INFO", format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p',
        filename="cyclic_deconvolute_%s_%s.log" % (os.path.basename(input_path).rsplit(".", 1)[0], start_time))
    logging.getLogger().addHandler(logging.StreamHandler(sys.stderr))
    input_path = str(input_path)
    is_mzml = input_path.lower().endswith(
        "mzml") or input_path.lower().endswith(".mzml.gz")

    reader = ms_deisotope.MSFileLoader(
        input_path, lockmass_config=lockmass_config)

    if no_product_splitting:
        writer = IonMobilityAware3DMzMLSerializer(
            open(output_path, 'wb'), len(reader))
    else:
        writer = IonMobilityAwareMzMLSerializer(
            open(output_path, 'wb'), len(reader))

    if is_mzml:
        reader = Generic3DIonMobilityFrameSource(reader)

    averagine_model: ms_deisotope.AveragineCache = ms_deisotope.AveragineCache(
        averagine_map[averagine])
    averagine_model.populate(truncate_after=0.95)
    scan_index = 0

    if start_time is not None:
        start_index = reader.get_frame_by_time(start_time).index
    else:
        start_index = 0
    if end_time is not None:
        end_index = reader.get_frame_by_time(end_time).index
    else:
        end_index = float('inf')

    strategy = make_iterator(reader, start_index=start_index,
                             stop_index=end_index, lock_mass_function=lockmass_function)
    writer.add_file_contents("MS1 spectrum")
    writer.add_file_contents("MSn spectrum")
    writer.add_data_processing(writer.build_processing_method())
    parent: IonMobilityFrame
    child: IonMobilityFrame

    pid = os.getpid()
    logger.info("Running on PID %d", pid)
    proc = psutil.Process(pid)
    with writer:
        for parent, children in strategy:
            if children:
                child = children[0]
            else:
                break
            if parent:
                logger.info("Begin Memory RSS: %s", proc.memory_percent())
                logger.info("Parent Frame: %s", parent)
                delta_im = np.median(np.diff(parent.ion_mobilities))
                max_gap_size = (delta_im * 3) + (delta_im / 2)
                parent.extract_features(max_gap_size=max_gap_size)
                logger.info("Parent Frame Contains %d raw features",
                            len(parent.features))
                parent.deconvolute_features(averagine_model,
                                            scorer=ms_deisotope.PenalizedMSDeconVFitter(
                                                5, 2),
                                            truncate_after=0.95, minimum_intensity=minimum_intensity,
                                            max_gap_size=max_gap_size)
                logger.info("Parent Frame Contains %d deconvoluted features",
                            len(parent.deconvoluted_features))
                if len(parent.deconvoluted_features) == 0:
                    logger.info(
                        "Skipping child frames, no viable features extracted from the parent.")
                    continue
                logger.info("Child Frame: %s", child)

            for child in children:
                delta_im = np.median(np.diff(parent.ion_mobilities))
                max_gap_size = (delta_im * 3) + (delta_im / 2)
                child.extract_features(max_gap_size=max_gap_size)
                logger.info("Child Frame Contains %d raw features",
                            len(child.features))
                child.deconvolute_features(
                    averagine_model, scorer=ms_deisotope.MSDeconVFitter(1), truncate_after=0.95,
                    minimum_intensity=minimum_intensity, max_gap_size=max_gap_size)
                logger.info("Child Frame Contains %d deconvoluted features",
                            len(child.deconvoluted_features))

            if no_product_splitting:
                if parent:
                    writer.save_scan(parent)
                for child in children:
                    writer.save_scan(child)
            else:
                if parent:
                    for child in children:
                        prec_prod_rels = connect_products_to_precursors(
                            parent, child)
                        scan_bunch, scan_index = frame_pair_to_scan_bunch(
                            parent, child, prec_prod_rels, scan_index)
                        logger.info("Produced %d Pseudo-MS2s",
                                    len(scan_bunch.products))
                        writer.save(scan_bunch)
            logger.info("End Memory RSS: %s", proc.memory_percent())


if __name__ == "__main__":
    register_debug_hook()
    main.main()
