'''A Tk-based spectrum viewer that integrates :mod:`ms_deisotope`'s spectrum processing
features.
'''
import sys
import os

import array

import matplotlib as mpl

import click

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,)
from matplotlib.widgets import SpanSelector
import numpy as np

try:
    import Tkinter as tk
    from Tkinter import Tk
    import ttk
except ImportError:
    import tkinter as tk
    from tkinter import Tk
    from tkinter import ttk


try:
    import tkFileDialog as tkfiledialog
except ImportError:
    import tkinter.filedialog as tkfiledialog

import ms_deisotope
from ms_deisotope.data_source import ScanBunch, ChargeNotProvided
from ms_deisotope.peak_set import EnvelopePair
from ms_deisotope.plot import (draw_raw, draw_peaklist, annotate_isotopic_peaks)

averagine_label_map = {
    "peptide": ms_deisotope.peptide,
    "glycan": ms_deisotope.glycan,
    "glycopeptide": ms_deisotope.glycopeptide,
    "heparan sulfate": ms_deisotope.heparan_sulfate
}


def interactive_console(**kwargs):
    import code
    code.interact(local=kwargs)


class Cursor(object):
    def __init__(self, ax, binding):
        self.ax = ax
        self.binding = binding

    def mouse_move(self, event):
        if not event.inaxes:
            return
        x, y = event.xdata, event.ydata
        if y > 1e4:
            text = "m/z=%1.2f, inten=%1.2g" % (x, y)
        else:
            text = "m/z=%1.2f, inten=%1.2f" % (x, y)
        self.binding.set(text)


class SummaryChromatogramBuilder(object):
    def __init__(self):
        self.tic = array.array('d')
        self.bpc = array.array('d')
        self.time = array.array('d')

    def add(self, scan):
        self.time.append(scan.scan_time)
        self.tic.append(scan.arrays.intensity.sum())
        self.bpc.append(scan.arrays.intensity.max())

    def reset(self):
        self.tic = array.array('d')
        self.bpc = array.array('d')
        self.time = array.array('d')


if issubclass(ttk.Frame, object):
    _BaseFrame = ttk.Frame
else:
    class _BaseFrame(object, ttk.Frame):
        pass


class SpectrumViewer(_BaseFrame):
    def __init__(self, master, denoise=4, ms1_score_threshold=20.0):
        ttk.Frame.__init__(self, master)
        self.root = master
        self._ms_file_name = None
        self.reader = None
        self.scan = None
        self.denoise = denoise
        self.ms1_score_threshold = ms1_score_threshold
        self.draw_isotopic_patterns = True
        self.configure_toolbar()
        self.configure_canvas()
        self.configure_treeview()
        self.configure_display_row()
        self.populate()
        self.draw_plot()

    @property
    def ms_file_name(self):
        return self._ms_file_name

    @ms_file_name.setter
    def ms_file_name(self, value):
        if value not in (None, '') and os.path.exists(value):
            print("Loading %r" % value)
            self._ms_file_name = value
            self.reader = ms_deisotope.MSFileLoader(self.ms_file_name)
            self.populate()

    def set_ms_file(self, value, populate=True):
        self._ms_file_name = value
        self.reader = ms_deisotope.MSFileLoader(self.ms_file_name)
        if populate:
            self.populate()

    def select_ms_file(self):
        file_name = tkfiledialog.askopenfilename()
        if file_name is not None and file_name != '':
            self.ms_file_name = file_name

    def do_layout(self):
        self.grid(sticky=tk.N + tk.W + tk.E + tk.S)
        tk.Grid.rowconfigure(self, 0, weight=1)
        tk.Grid.columnconfigure(self, 0, weight=1)

    def configure_canvas(self):
        self.figure = Figure(dpi=100)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.axis = self.figure.add_subplot(111)
        self.canvas.draw()
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.grid(row=0, column=0, sticky=tk.N + tk.W + tk.E + tk.S)
        self.canvas_cursor = Cursor(self.axis, tk.StringVar(master=self.root))
        self.canvas.mpl_connect('motion_notify_event', self.canvas_cursor.mouse_move)
        self.span = SpanSelector(
            self.axis, self.zoom, 'horizontal', useblit=True,
            rectprops=dict(alpha=0.5, facecolor='red'))
        self.mz_span = None
        self.scan = None
        self.annotations = []
        self.canvas.draw()

    def configure_toolbar(self):
        self.toolbar = tk.Menu(self)
        self.toolbar.add_command(label='Open', command=self.select_ms_file)
        self.toolbar.add_command(label='Interact', command=self._interact)
        self.root.config(menu=self.toolbar)

    def _interact(self):
        local_vars = {
            "self": self,
            "scan": self.scan,
        }
        return interactive_console(**local_vars)

    def _zoom_in(self, min_mz, max_mz):
        arrays = self.scan.arrays
        subset = arrays.between_mz(min_mz, max_mz).intensity
        if len(subset) > 0 and subset.max() > 0 and subset.min() < subset.max():
            most_intense = np.max(subset)
            y_cap = most_intense * 1.2
            delta = max_mz - min_mz
            if delta > 20:
                max_mz += min(delta * 0.1, 100)
                min_mz -= min(delta * 0.1, 100)
            self.mz_span = (min_mz, max_mz)
            self.axis.set_xlim(min_mz, max_mz)
            self.axis.set_ylim(0, y_cap)
            self.annotate_plot(min_mz, max_mz)
            self.figure.canvas.draw()

    def zoom(self, xmin, xmax):
        # If the bounds are not close, we're zooming in
        if not self.scan:
            return
        arrays = self.scan.arrays
        # Default zoom
        if (xmax - xmin) <= .1:
            min_peak = 0
            max_peak = len(arrays[0]) - 1
            self.mz_span = None
            xmin, xmax = arrays.mz[min_peak], arrays.mz[max_peak]
            xmin = max(xmin - 100, 0)
        if (xmin - 50) < 0:
            if xmin > -20:
                xmin = -20
        self._zoom_in(xmin, xmax)

    def clear_annotations(self):
        for anno in self.annotations:
            try:
                anno.remove()
            except ValueError:
                break
        self.annotations = []

    def annotate_plot(self, min_mz=None, max_mz=None):
        self.clear_annotations()
        if min_mz is None:
            min_mz = 0
        if max_mz is None:
            max_mz = self.scan.arrays.mz[-1]
        subset = self.scan.deconvoluted_peak_set.between(min_mz, max_mz, use_mz=True)
        if not subset:
            return
        threshold = 0.0
        threshold_list = ([max(i.intensity for i in p.envelope) for p in subset])
        if threshold_list:
            threshold = np.mean(threshold_list)
        threshold_list = ([max(i.intensity for i in p.envelope) for p in subset
                           if max(i.intensity for i in p.envelope) > threshold])
        if threshold_list:
            threshold = np.mean(threshold_list)
        for peak in subset:
            if peak.intensity > threshold:
                label = '%0.2f (%d)' % (peak.neutral_mass, peak.charge)
                pt = max(peak.envelope, key=lambda x: x.intensity)
                y = pt.intensity * 1.05
                x = np.average(
                    [p.mz for p in peak.envelope],
                    weights=[p.intensity for p in peak.envelope])
                self.annotations.append(
                    self.axis.text(x, y, label, ha='center', clip_on=True, fontsize=10))

    def _process_scan(self):
        self.scan.pick_peaks(signal_to_noise_threshold=1.5)
        if self.scan.ms_level == 1:
            averagine_value = self.ms1_averagine_combobox.get()
            averagine_value = averagine_label_map[averagine_value]
            truncate_to = 0.95
            scorer = ms_deisotope.PenalizedMSDeconVFitter(self.ms1_score_threshold, 2.)
        else:
            averagine_value = self.msn_averagine_combobox.get()
            averagine_value = averagine_label_map[averagine_value]
            truncate_to = 0.8
            scorer = ms_deisotope.MSDeconVFitter(10.)
        max_charge_state = self.max_charge_state_var.get()
        min_charge_state = self.min_charge_state_var.get()
        if self.scan.ms_level > 1:
            prec_charge = self.scan.precursor_information.charge
            if prec_charge and prec_charge != ChargeNotProvided:
                max_charge_state = prec_charge
                if prec_charge < 0:
                    min_charge_state = -abs(min_charge_state)
        self.scan.deconvolute(
            averagine=averagine_value,
            scorer=scorer,
            truncate_after=1 - 1e-4,
            incremental_truncation=truncate_to,
            charge_range=(min_charge_state, max_charge_state))

    def draw_plot(self, scan=None, children=None):
        if children is None:
            children = []
        if scan is None or scan.arrays.mz.shape[0] == 0:
            return
        self.axis.clear()
        self.scan = scan
        self.children_scans = children
        if scan.ms_level == 1:
            averaging_level = self.ms1_scan_averaging_var.get()
            if averaging_level > 0:
                print("Averaging")
                if scan.arrays.mz.shape[0] > 1 :
                    self.scan = scan = scan.average(averaging_level)
            print("Denoising")
            self.scan = scan.denoise(self.denoise) if self.denoise else scan
        print("Processing Peaks")
        self._process_scan()
        print("Begin Drawing")
        scan = self.scan
        if scan.is_profile:
            print("Drawing Profile")
            draw_raw(*scan.arrays, ax=self.axis, color='black', lw=0.75)
            # lock in axis dimensions for profile spectra according to the raw profile
            self.axis.set_xlim(0, max(self.axis.get_xlim()))
        if scan.peak_set:
            print("Drawing Centroids")
            draw_peaklist(scan.peak_set, ax=self.axis, alpha=0.5, lw=1, color='grey', linestyle='--')
            if not scan.is_profile:
                # lock in axis dimensions for centroid spectra according to the picked peaks
                self.axis.set_xlim(0, max(self.axis.get_xlim()) + 2)
        if scan.deconvoluted_peak_set:
            print("Drawing Isotopic Envelope")
            draw_peaklist(
                [i for p in scan.deconvoluted_peak_set for i in p.envelope],
                ax=self.axis, alpha=0.6, lw=0.5, color='orange')
            draw_peaklist(
                [p.envelope[0] for p in scan.deconvoluted_peak_set if p.envelope[0].intensity > 0],
                ax=self.axis, alpha=0.6, lw=1,
                color='red')

        # draw isolation window and instrument reported precursor
        if children:
            ylim = scan.arrays.intensity.max()
            for child in children:
                if child.isolation_window and not child.isolation_window.is_empty():
                    self.axis.vlines(
                        [child.isolation_window.lower_bound, child.isolation_window.upper_bound],
                        0, ylim, linestyle='--', color='skyblue', lw=1)
                self.axis.vlines(child.precursor_information.mz, 0, ylim,
                                 linestyle='--', color='black', lw=1)
        if self.draw_isotopic_patterns:
            annotate_isotopic_peaks(self.scan, self.axis)
        if self.scan.precursor_information:
            ylim = scan.arrays.intensity.max()
            self.axis.vlines(self.scan.precursor_information.mz, 0, ylim, linestyle='-.', color='orange')
        if self.mz_span is not None:
            self._zoom_in(*self.mz_span)
        else:
            self.annotate_plot(None, None)
        self.canvas.draw()

    def on_row_click(self, event):
        selection = self.treeview.focus()
        item = self.treeview.item(selection)
        index = item['text']

        if self.reader is not None:
            scan = self.reader.get_scan_by_index(index)
            if scan.ms_level == 1:
                bunch = next(self.reader.start_from_scan(scan.id))
                if not isinstance(bunch, ScanBunch):
                    children = []
                else:
                    children = bunch.products
            else:
                children = []
            self.draw_plot(scan, children)
        else:
            self.draw_plot(None)

    def configure_display_row(self):
        self.display_row = ttk.Frame(self)
        self.display_row.grid(row=1, column=0, sticky=tk.W + tk.S + tk.E)
        self.cursor_label = ttk.Label(self.display_row, text=" " * 20)
        self.cursor_label.grid(row=0, padx=(10, 10))

        def update_label(*args, **kwargs):
            self.cursor_label['text'] = self.canvas_cursor.binding.get()

        self.canvas_cursor.binding.trace('w', update_label)

        self.ms1_averagine_combobox = ttk.Combobox(self.display_row, values=[
            "peptide",
            "glycan",
            "glycopeptide",
            "heparan sulfate",
        ], width=15)
        self.ms1_averagine_combobox.set("glycopeptide")
        self.ms1_averagine_combobox_label = ttk.Label(self.display_row, text="MS1 Averagine:")
        self.ms1_averagine_combobox_label.grid(row=0, column=1, padx=(10, 0))
        self.ms1_averagine_combobox.grid(row=0, column=2, padx=(1, 10))

        self.msn_averagine_combobox = ttk.Combobox(self.display_row, values=[
            "peptide",
            "glycan",
            "glycopeptide",
            "heparan sulfate",
        ], width=15)
        self.msn_averagine_combobox.set("peptide")
        self.msn_averagine_combobox_label = ttk.Label(self.display_row, text="MSn Averagine:")
        self.msn_averagine_combobox_label.grid(row=0, column=3, padx=(10, 0))
        self.msn_averagine_combobox.grid(row=0, column=4, padx=(1, 10))

        self.ms1_scan_averaging_label = ttk.Label(
            self.display_row, text="MS1 Signal Averaging:")
        self.ms1_scan_averaging_label.grid(row=0, column=5, padx=(10, 0))
        self.ms1_scan_averaging_var = tk.IntVar(self, value=2)
        self.ms1_scan_averaging = ttk.Entry(self.display_row, width=3)
        self.ms1_scan_averaging['textvariable'] = self.ms1_scan_averaging_var
        self.ms1_scan_averaging.grid(row=0, column=6, padx=(1, 3))

        self.min_charge_state_var = tk.IntVar(self, value=1)
        self.max_charge_state_var = tk.IntVar(self, value=12)
        self.min_charge_state_label = ttk.Label(
            self.display_row, text="Min Charge:")
        self.min_charge_state_label.grid(row=0, column=7, padx=(5, 0))
        self.min_charge_state = ttk.Entry(self.display_row, width=3)
        self.min_charge_state['textvariable'] = self.min_charge_state_var
        self.min_charge_state.grid(row=0, column=8, padx=(1, 3))

        self.max_charge_state_label = ttk.Label(
            self.display_row, text="Max Charge:")
        self.max_charge_state_label.grid(row=0, column=9, padx=(5, 0))
        self.max_charge_state = ttk.Entry(self.display_row, width=3)
        self.max_charge_state['textvariable'] = self.max_charge_state_var
        self.max_charge_state.grid(row=0, column=10, padx=(1, 3))

    def configure_treeview(self):
        self.treeview = ttk.Treeview(self)
        self.treeview['columns'] = ["id", "time", 'ms_level', 'precursor_mz', 'precursor_charge', 'activation']
        self.treeview.grid(row=2, column=0, sticky=tk.S + tk.W + tk.E + tk.N)

        self.treeview_scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.treeview.yview)
        self.treeview_scrollbar.grid(row=2, column=0, sticky=tk.S + tk.E + tk.N)
        self.treeview.configure(yscrollcommand=self.treeview_scrollbar.set)

        self.treeview.heading('id', text="Scan ID")
        self.treeview.heading('#0', text='Index')
        self.treeview.heading("time", text='Time (min)')
        self.treeview.heading("ms_level", text='MS Level')
        self.treeview.heading("precursor_mz", text='Precursor M/Z')
        self.treeview.heading("precursor_charge", text='Precursor Z')
        self.treeview.heading("activation", text='Activation')
        self.treeview.column("#0", width=75)
        self.treeview.column("ms_level", width=75)
        self.treeview.column("time", width=75)
        self.treeview.column("precursor_mz", width=100)
        self.treeview.column("precursor_charge", width=100)
        self.treeview.bind("<<TreeviewSelect>>", self.on_row_click)

    def clear_treeview(self):
        children = self.treeview.get_children()
        if children:
            self.treeview.delete(*children)

    def _populate_range(self, start, stop):
        scans = []
        ended = False
        for i in range(start, stop):
            try:
                scan = self.reader[i]
            except Exception:
                ended = True
                break
            if scan.index % 5000 == 0:
                print(scan)
            i = scan.index
            values = [scan.id, "%0.4f" % scan.scan_time, scan.ms_level]
            if scan.ms_level > 1:
                values.extend([scan.precursor_information.mz, scan.precursor_information.charge,
                               str(scan.activation)])
            else:
                values.extend(['-', '-', '-'])
            scans.append(values)
        i = start
        for values in scans:
            self.treeview.insert('', 'end', values=values, text=i)
            i += 1
        if not ended:
            self.after(100, self._populate_range, stop, stop + 500)

    def populate(self, clear=True):
        if clear:
            self.clear_treeview()
        if self.reader is not None:
            self.reader.make_iterator(grouped=False)
            self.after(10, self._populate_range, 0, 500)


@click.command("ms-view")
@click.argument("filename", type=click.Path(exists=True, readable=True), required=False)
@click.option("-d", "--denoise", type=float, default=5.0)
def main(filename=None, denoise=5, ms1_score=20.0, msn_score=10.0, ms1_truncation=0.95, msn_truncation=0.8):
    base = Tk()
    base.title("ms_deisotope Spectrum Viewer")
    tk.Grid.rowconfigure(base, 0, weight=1)
    tk.Grid.columnconfigure(base, 0, weight=1)
    app = SpectrumViewer(base, denoise=denoise)
    app.do_layout()
    app.ms_file_name = filename

    app.mainloop()


if __name__ == '__main__':
    main.main()
