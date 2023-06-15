
import logging
import os
from typing import Any, Callable, Dict, List, Optional, Tuple, Union
import warnings

from six import string_types as basestring

import numpy as np

from ms_deisotope import constants, MSFileLoader
from ms_deisotope.data_source.scan.base import ScanBunch
from ms_deisotope.data_source.scan.loader import ScanIterator
from ms_deisotope.task import LogUtilsMixin
from ms_deisotope.utils import Base
from ms_deisotope.averagine import AveragineCache, PROTON

from ms_deisotope.data_source.scan.mobility_frame import IonMobilityFrame, IonMobilitySourceRandomAccessFrameSource, Generic3DIonMobilityFrameSource

logger = logging.getLogger("ms_deisotope.frame_processor")
logger.addHandler(logging.NullHandler())


def _loader_creator(specification, **kwargs):
    if isinstance(specification, (tuple, list)):
        specification, options = specification
        if options:
            kwargs.update(options)
    if isinstance(specification, (basestring, os.PathLike)):
        specification = MSFileLoader(specification, **kwargs)
    if isinstance(specification, IonMobilitySourceRandomAccessFrameSource):
        return specification
    elif isinstance(specification, ScanIterator):
        return Generic3DIonMobilityFrameSource(specification, **kwargs)
    else:
        raise ValueError(
            "Cannot determine how to get a frame iterator from %r" % (specification,))


class IonMobilityFrameProcessor(Base, LogUtilsMixin):
    data_source: Union[os.PathLike, IonMobilitySourceRandomAccessFrameSource]
    loader_type: Callable

    ms1_averaging: int
    msn_averaging: int
    num_threads: int

    ms1_peak_picking_args: Dict[str, Any]
    msn_peak_picking_args: Dict[str, Any]

    ms1_deconvolution_args: Dict[str, Any]
    msn_deconvolution_args: Dict[str, Any]


    def __init__(self, data_source, ms1_peak_picking_args=None,
                 msn_peak_picking_args=None,
                 ms1_deconvolution_args=None,
                 msn_deconvolution_args=None,
                 loader_type=None,
                 terminate_on_error=True,
                 ms1_averaging=0,
                 msn_averaging=0, num_threads=3):
        if loader_type is None:
            loader_type = _loader_creator

        self.data_source = data_source
        self.ms1_peak_picking_args = ms1_peak_picking_args or {}
        self.msn_peak_picking_args = msn_peak_picking_args or ms1_peak_picking_args or {}
        self.ms1_deconvolution_args = ms1_deconvolution_args or {}
        self.ms1_deconvolution_args.setdefault("charge_range", (1, 8))
        self.msn_deconvolution_args = msn_deconvolution_args or {}
        self.msn_deconvolution_args.setdefault("charge_range", (1, 8))

        self.ms1_averaging = int(ms1_averaging) if ms1_averaging else 0
        self.msn_averaging = int(msn_averaging) if msn_averaging else 0
        self.num_threads = num_threads
        self.loader_type = loader_type

        self._signal_source = self.loader_type(data_source)
        self.terminate_on_error = terminate_on_error
        self._prepopulate_averagine_cache()

    @property
    def reader(self) -> IonMobilitySourceRandomAccessFrameSource:
        """
        The :class:`~.IonMobilitySourceRandomAccessFrameSource` which generates the raw scans that will
        be processed.

        Returns
        -------
        :class:`~.IonMobilitySourceRandomAccessFrameSource`
        """
        return self._signal_source

    def _prepopulate_averagine_cache(self):
        if 'averagine' in self.ms1_deconvolution_args:
            averagine = self.ms1_deconvolution_args['averagine']
            ms1_truncate_after = self.ms1_deconvolution_args.get(
                'truncate_after', constants.TRUNCATE_AFTER)
            ms1_ignore_below = self.ms1_deconvolution_args.get(
                'ignore_below', constants.IGNORE_BELOW)
            ms1_charge_range = self.ms1_deconvolution_args.get(
                'charge_range', (1, 8))
            ms1_charge_carrier = self.ms1_deconvolution_args.get(
                'charge_carrier', PROTON)
            if isinstance(averagine, (list, tuple)):
                averagine = [
                    AveragineCache(a).populate(
                        truncate_after=ms1_truncate_after,
                        ignore_below=ms1_ignore_below,
                        min_charge=ms1_charge_range[0],
                        max_charge=ms1_charge_range[1],
                        charge_carrier=ms1_charge_carrier)
                    for a in averagine]
            else:
                averagine = AveragineCache(averagine).populate(
                    truncate_after=ms1_truncate_after,
                    ignore_below=ms1_ignore_below,
                    min_charge=ms1_charge_range[0],
                    max_charge=ms1_charge_range[1],
                    charge_carrier=ms1_charge_carrier)
            self.ms1_deconvolution_args['averagine'] = averagine
        if 'averagine' in self.msn_deconvolution_args:
            averagine = self.msn_deconvolution_args['averagine']
            msn_truncate_after = self.msn_deconvolution_args.get(
                'truncate_after', constants.TRUNCATE_AFTER)
            msn_ignore_below = self.msn_deconvolution_args.get(
                'ignore_below', constants.IGNORE_BELOW)
            msn_charge_range = self.msn_deconvolution_args.get(
                'charge_range', (1, 8))
            msn_charge_carrier = self.msn_deconvolution_args.get(
                'charge_carrier', PROTON)
            if isinstance(averagine, (list, tuple)):
                averagine = [
                    AveragineCache(a).populate(
                        truncate_after=msn_truncate_after,
                        ignore_below=msn_ignore_below,
                        min_charge=msn_charge_range[0],
                        max_charge=msn_charge_range[1],
                        charge_carrier=msn_charge_carrier
                    ) for a in averagine]
            else:
                averagine = AveragineCache(averagine).populate(
                    truncate_after=msn_truncate_after,
                    ignore_below=msn_ignore_below,
                    min_charge=msn_charge_range[0],
                    max_charge=msn_charge_range[1],
                    charge_carrier=msn_charge_carrier)
            self.msn_deconvolution_args['averagine'] = averagine

    def process_frame_group(self, precursor_frame: IonMobilityFrame, product_frames: List[IonMobilityFrame]) -> Tuple[IonMobilityFrame, None, List[IonMobilityFrame]]:
        if precursor_frame is not None:
            self.extract_precursor_features(precursor_frame)
        return precursor_frame, None, product_frames

    def _default_max_gap_size(self, frame: IonMobilityFrame) -> float:
        delta_im = np.median(np.diff(frame.ion_mobilities))
        max_gap_size = (delta_im * 3) + (delta_im / 2)
        return max_gap_size

    def extract_precursor_features(self, precursor_frame: IonMobilityFrame):
        self.log(f"Extracting MS1 raw features from {precursor_frame.id}")
        options = self.ms1_peak_picking_args.copy()
        options.setdefault("average_within", self.ms1_averaging)
        options.setdefault("num_threads", self.num_threads)

        if "max_gap_size" not in options:
            options["max_gap_size"] = self._default_max_gap_size(precursor_frame)

        precursor_frame.extract_features(**options)

    def deconvolute_precursor_features(self, precursor_frame: IonMobilityFrame):
        self.log(f"Deconvolving MS1 IMS features from {precursor_frame.id}")
        options = self.ms1_deconvolution_args.copy()

        if "max_gap_size" not in options:
            options["max_gap_size"] = self._default_max_gap_size(
                precursor_frame)

        precursor_frame.deconvolute_features(**options)

    def extract_product_features(self, product_frame: IonMobilityFrame):
        self.log(f"Extracting MSn raw features from {product_frame.id}")
        options = self.msn_peak_picking_args.copy()
        options.setdefault("average_within", self.msn_averaging)
        options.setdefault("num_threads", self.num_threads)

        if "max_gap_size" not in options:
            options["max_gap_size"] = self._default_max_gap_size(
                product_frame)

        product_frame.extract_features(**options)

    def deconvolute_product_features(self, product_frame: IonMobilityFrame):
        self.log(f"Deconvolving MSn IMS features from {product_frame.id}")
        options = self.msn_deconvolution_args.copy()

        if "max_gap_size" not in options:
            options["max_gap_size"] = self._default_max_gap_size(
                product_frame)

        product_frame.deconvolute_features(**options)

    def process(self, precursor: IonMobilityFrame, products: List[IonMobilityFrame]) -> ScanBunch:
        precursor, _, products = self.process_frame_group(precursor, products)
        if precursor is not None:
            self.deconvolute_precursor_features(precursor)
        for product in products:
            self.extract_product_features(product)
            self.deconvolute_product_features(product)
        return ScanBunch(precursor, products)

    def __next__(self):
        precursor, products = next(self.reader)
        return self.process(precursor, products)

    def __iter__(self) -> 'IonMobilityFrameProcessor':
        return self

    def reset(self):
        self.reader.reset()

    def start_from_frame(self, *args, **kwargs) -> 'IonMobilityFrameProcessor':
        """
        A wrapper around :meth:`~.IonMobilitySourceRandomAccessFrameSource.start_from_frame` provided by
        :attr:`reader`, if available.

        Returns
        -------
        self

        See Also
        --------
        :meth:`~.IonMobilitySourceRandomAccessFrameSource.start_from_frame`
        """
        self.reader.start_from_frame(*args, **kwargs)
        return self


IonMobilityFrameProcessor.log_with_logger(logger)
