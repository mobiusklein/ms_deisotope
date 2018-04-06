# pragma: no cover
try:
    from matplotlib import font_manager
    from matplotlib import pyplot as plt
    from matplotlib import patches as mpatches
    from matplotlib.textpath import TextPath
except (RuntimeError, ImportError):
    pass

from ms_peak_picker.utils import draw_peaklist

from .subgraph import (layout_layers, peak_overlap)

font_options = font_manager.FontProperties(family='monospace')


def draw_envelope_subgraph(envelopes, scale_factor=1.0, overlap_fn=peak_overlap, ax=None, **kwargs):
    layers = layout_layers(envelopes)
    max_score = max(e.score for e in envelopes)
    peaks = set()

    for e in envelopes:
        peaks.update(e.fit.experimental)

    peaks = sorted(peaks, key=lambda x: x.mz)
    peaks = [p.clone() for p in peaks]
    total_intensity = sum(p.intensity for p in peaks)

    start = peaks[0].mz
    end = peaks[-1].mz

    if ax is None:
        figure, ax = plt.subplots(1, 1)

    row_width = float('inf')
    annotation_text_size = 3. * scale_factor
    layer_height = 0.56 * scale_factor
    y_step = (layer_height + 0.05) * - scale_factor
    origin_y = cur_y = -layer_height - 0.075

    cur_position = peaks[0].mz

    for layer in layers:
        layer.sort(key=lambda x: x.start)

    while cur_position < end:
        next_row = cur_position + row_width

        for layer in layers:
            c = 0
            for envelope in layer:
                if envelope.start < cur_position:
                    continue
                elif envelope.start > next_row:
                    break
                c += 1
                rect = mpatches.Rectangle(
                    (envelope.start - 0.01, cur_y),
                    width=0.01 + envelope.end - envelope.start, height=layer_height,
                    facecolor='lightblue', edgecolor='black', linewidth=0.15,
                    alpha=min(max(envelope.score / max_score, 0.2), 0.8))
                ax.add_patch(rect)
                text_path = TextPath(
                    (envelope.start + 0.1, cur_y + 0.2),
                    "%0.2f, %d" % (envelope.score, envelope.fit.charge), size=annotation_text_size / 14.5,
                    prop=font_options, stretch=200)
                patch = mpatches.PathPatch(text_path, facecolor='grey', lw=0.04)
                ax.add_patch(patch)

            if c > 0:
                cur_y += y_step
        cur_y += y_step / 5
        cur_position = next_row

    for p in peaks:
        p.intensity = (p.intensity / total_intensity) * abs(origin_y - cur_y) * 8

    draw_peaklist(peaks, ax=ax)

    ax.set_ylim(cur_y, max(p.intensity for p in peaks) + 0.2)
    ax.set_xlim(start - 0.2, end + 0.2)
    ax.axes.get_yaxis().set_visible(False)
    return ax
