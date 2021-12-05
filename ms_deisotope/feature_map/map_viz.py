from itertools import cycle
from collections import defaultdict

try:
    from matplotlib import pyplot as plt
    from matplotlib import patches as mpatches, colors as mcolors, collections as mcollections
    from matplotlib.colors import Normalize
except (RuntimeError, ImportError):
    pass
import numpy as np


from .profile_transform import interpolate, smooth_leveled


Ellipse = mpatches.Ellipse
FancyBboxPatch = mpatches.FancyBboxPatch
Rectangle = mpatches.Rectangle
LineCollection = mcollections.LineCollection


def compute_intensity_range(features):
    min_int = float('inf')
    max_int = 0
    for feature in features:
        for node in feature:
            t = node.total_intensity()
            min_int = min(min_int, t)
            max_int = max(max_int, t)
    return min_int, max_int


def feature_to_segments(feature):
    mzs = []
    times = []
    intensities = []
    for node in feature:
        mzs.append(node.mz)
        times.append(node.time)
        intensities.append(node.total_intensity())
    points = np.array([mzs, times]).T.reshape(-1, 1, 2)
    return np.concatenate([points[:-1], points[1:]], axis=1), np.array(intensities)


def draw_features(features, ax=None, alpha=0.65, norm=None, cmap=None, **kwargs):
    if norm is None:
        norm = Normalize(*compute_intensity_range(features))
    if ax is None:
        fig, ax = plt.subplots(1)
    if not features:
        return ax
    lines = []
    kwargs.setdefault("lw", 0.05)
    lw = kwargs.get("linewidth", kwargs.get("lw"))
    for feat in features:
        segments, intensities = feature_to_segments(feat)
        art = LineCollection(segments, linewidths=lw, cmap=cmap, alpha=alpha)
        art.set_array(intensities)
        lines.append(art)

    for line in lines:
        ax.add_collection(line)

    ax.set_xlim(
        min(features, key=lambda x: x.mz if x is not None else float('inf')).mz - 1,
        max(features, key=lambda x: x.mz if x is not None else -float('inf')).mz + 1)
    ax.set_ylim(
        min(features, key=lambda x: x.start_time if x is not None else float(
            'inf')).start_time - 1,
        max(features, key=lambda x: x.end_time if x is not None else -float('inf')).end_time + 1)
    return ax


def draw_feature_sets(feature_sets, ax=None, alpha=0.65, width=2e-5, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1)

    kwargs.setdefault("lw", 0.05)
    lw = kwargs.get("linewidth", kwargs.get("lw"))

    features = []
    ellipse_sets = []
    for feature_set in feature_sets:
        ellipses = []
        for feat in feature_set:
            if feat is None:
                continue
            center = (feat.end_time + feat.start_time) / 2.
            height = feat.end_time - feat.start_time

            center_mz = feat.mz
            mz_width = center_mz * width

            ellipses.append(
                FancyBboxPatch((feat.mz - mz_width / 4., center - height / 2.),
                               width=mz_width / 2., height=height,
                               boxstyle=mpatches.BoxStyle.Round(pad=mz_width / 2.)))
            features.append(feat)
        ellipse_sets.append(ellipses)

    for ellipses in ellipse_sets:
        color = np.random.rand(3)
        for ell in ellipses:
            ell.set_alpha(alpha)
            ell.set_facecolor(color)
            ell.set_edgecolor(color)
            ell.set_linewidth(lw)
            ax.add_artist(ell)

    ax.set_xlim(
        min(features, key=lambda x: x.mz if x is not None else float('inf')).mz - 1,
        max(features, key=lambda x: x.mz if x is not None else -float('inf')).mz + 1)
    ax.set_ylim(
        min(features, key=lambda x: x.start_time if x is not None else float('inf')).start_time - 1,
        max(features, key=lambda x: x.end_time if x is not None else -float('inf')).end_time + 1)
    return ax


_nice_color_cycle = ([
    'blue', 'red', 'green', 'pink', 'orange', "grey", "purple",
    "seagreen", "darkgoldenrod", "darkcyan", "skyblue",
    "maroon", "darkgreen", "slategrey", "darkslateblue"
])


nice_color_cycle = cycle(_nice_color_cycle)


def random_colorizer(profile, *args, **kwargs):
    return mcolors.rgb2hex(np.random.rand(3))


def nice_colorizer(profile, *args, **kwargs):
    return next(nice_color_cycle)


def labeler(profile, *args, **kwargs):
    label = "%0.4f" % profile.mz
    if hasattr(profile, 'charge'):
        label += ", %d" % profile.charge
    return label


def draw_profiles(profiles, ax=None, smooth=False, interp=False, label_font_size=10,
                  axis_label_font_size=16, axis_font_size=16, label=True,
                  colorizer=random_colorizer, label_function=labeler):
    if ax is None:
        _fig, ax = plt.subplots(1)

    if not profiles:
        return ax
    minimum_ident_time = float("inf")
    maximum_ident_time = 0
    maximum_intensity = 0

    _label_apexes = label and (label_function is not None)

    for profile in profiles:
        if profile is None:
            continue
        rt, heights = profile.as_arrays()

        if smooth:
            heights = smooth_leveled(rt, heights, smooth)
        if interp:
            rt, heights = interpolate(rt, heights)

        maximum_ident_time = max(max(rt), maximum_ident_time)
        minimum_ident_time = min(min(rt), minimum_ident_time)

        label = label_function(profile) if label_function is not None else None

        maximum_intensity = max(max(heights), maximum_intensity)
        color = colorizer(profile)
        ax.scatter(rt, heights, color=color)
        ax.fill_between(
            rt,
            heights,
            alpha=0.25,
            color=color,
            label=label
        )
        apex = max(heights)
        apex_ind = np.argmax(heights)
        rt_apex = rt[apex_ind]
        if _label_apexes:
            ax.text(rt_apex, apex + min(apex * (
                .1), 1200), label,
                ha='center', fontsize=label_font_size, alpha=0.75)
    ax.set_xlim(minimum_ident_time - 0.02,
                maximum_ident_time + 0.02)
    ax.set_ylim(0, maximum_intensity * 1.1)
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.), ncol=2, fontsize=10)
    ax.axes.spines['right'].set_visible(False)
    ax.axes.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.set_xlabel("Time", fontsize=axis_label_font_size)
    ax.set_ylabel("Relative Abundance", fontsize=axis_label_font_size)
    [t.set(fontsize=axis_font_size) for t in ax.get_xticklabels()]
    [t.set(fontsize=axis_font_size) for t in ax.get_yticklabels()]
    return ax


def features_to_peak_time_pairs(features):
    pairs = []
    for feature in features:
        for node in feature:
            for peak in node.members:
                pairs.append((peak, node.time))
    pairs.sort(key=lambda x: x[1])
    return pairs


def extract_intensity_array(peaks):
    mzs = []
    intensities = []
    rts = []
    current_mzs = []
    current_intensities = []
    last_time = None
    time = None
    for peak, time in peaks:
        if time != last_time:
            if last_time is not None:
                mzs.append(current_mzs)
                intensities.append(current_intensities)
                rts.append(time)
            last_time = time
            current_mzs = []
            current_intensities = []
        current_mzs.append(peak.mz)
        current_intensities.append(peak.intensity)
    if time is not None:
        mzs.append(current_mzs)
        intensities.append(current_intensities)
        rts.append(time)
    return mzs, intensities, rts


def binner(x):
    return np.floor(np.array(x) / 10.) * 10


def make_map(mzs, intensities):
    binned_mzs = [
        binner(mz_row) for mz_row in mzs
    ]
    unique_mzs = set()
    map(unique_mzs.update, binned_mzs)
    unique_mzs = np.array(sorted(unique_mzs))
    assigned_bins = []
    j = 0
    for mz, inten in zip(mzs, intensities):
        j += 1
        mz = binner(mz)
        bin_row = defaultdict(float)
        for i in range(len(mz)):
            k = mz[i]
            v = inten[i]
            bin_row[k] += v
        array_row = np.array([bin_row[m] for m in unique_mzs])
        assigned_bins.append(array_row)

    assigned_bins = np.vstack(assigned_bins)
    assigned_bins[assigned_bins == 0] = 1
    return assigned_bins, unique_mzs


def render_map(assigned_bins, rts, unique_mzs, ax=None, color_map=None, scaler=np.sqrt, colorbar=True):
    if ax is None:
        _fig, ax = plt.subplots(1)
    if color_map is None:
        color_map = plt.cm.viridis
    data = np.array(scaler(assigned_bins.T))
    img = ax.pcolormesh(data, cmap=color_map)

    xticks = ax.get_xticks()
    newlabels = np.array(ax.get_xticklabels())[np.arange(0, len(xticks))]
    n = len(newlabels)
    interp = np.linspace(rts[0], rts[-1], n)
    for i, label in enumerate(newlabels):
        num = round(interp[i], 1)
        label.set_rotation(90)
        label.set_text(str((num)))

    # ax.set_xticklabels(newlabels)

    yticks = ax.get_yticks()
    n = len(yticks)
    # print n, yticks
    newlabels = np.array(ax.get_yticklabels())[np.arange(0, len(yticks))]

    va = unique_mzs
    # Hacky and doesn't align quite right
    n = len(newlabels)
    # if n % 2 != 0:
    #     n -= 1
    interp = np.linspace(va[0], va[-1], n)
    newlabels = []
    # for i, label in enumerate(newlabels):
    for i in range(n):
        num = interp[i]
        newlabels.append(str(round(num, 2)))
    ax.set_yticklabels((newlabels))
    ax.set_xlabel("Time")
    ax.set_ylabel("m/z")
    if colorbar:
        _cbar = ax.figure.colorbar(img, ax=ax)
    return ax


def heatmap(features, ax=None, color_map=None, scaler=np.sqrt, colorbar=True):
    mzs, intensities, rts = extract_intensity_array(
        features_to_peak_time_pairs(features))
    binned_intensities, unique_mzs = make_map(mzs, intensities)
    return render_map(
        binned_intensities, rts, unique_mzs, ax=ax, color_map=color_map,
        scaler=scaler, colorbar=colorbar)
