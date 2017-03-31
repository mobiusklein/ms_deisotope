from itertools import cycle
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import numpy as np


from .profile_transform import sliding_mean, sliding_median, gaussian_smooth, interpolate, smooth_leveled


Ellipse = mpatches.Ellipse


def draw_features(features, ax=None, alpha=0.65, width=0.025, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(1)

    ellipses = []
    kwargs.setdefault("lw", 0.05)
    lw = kwargs.get("linewidth", kwargs.get("lw"))
    for feat in features:
        if feat is None:
            continue
        center = (feat.end_time + feat.start_time) / 2.
        height = feat.end_time - feat.start_time
        ellipses.append(
            Ellipse((feat.mz, center), width=width, height=height, angle=0))

    for ell in ellipses:
        ell.set_alpha(alpha)
        ell.set_facecolor("blue")
        ell.set_edgecolor("blue")
        ell.set_linewidth(lw)
        ax.add_artist(ell)

    ax.set_xlim(
        min(features, key=lambda x: x.mz if x is not None else float('inf')).mz - 1,
        max(features, key=lambda x: x.mz if x is not None else -float('inf')).mz + 1)
    ax.set_ylim(
        min(features, key=lambda x: x.start_time if x is not None else float('inf')).start_time - 1,
        max(features, key=lambda x: x.end_time if x is not None else -float('inf')).end_time + 1)
    return ax


def draw_feature_sets(feature_sets, ax=None, alpha=0.65, width=0.025, **kwargs):
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
            ellipses.append(
                Ellipse((feat.mz, center), width=width, height=height, angle=0))
            features.append(feat)
        ellipse_sets.append(ellipses)

    for ellipses in ellipse_sets:
        color = np.random.rand(3.)
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
    return np.random.rand(3.)


def nice_colorizer(profile, *args, **kwargs):
    return next(nice_color_cycle)


def labeler(profile, *args, **kwargs):
    label = "%0.4f" % profile.mz
    if hasattr(profile, 'charge'):
        label += ", %d" % profile.charge
    return label


def draw_profiles(profiles, ax=None, smooth=False, interp=False, label_font_size=10,
                  axis_label_font_size=20, axis_font_size=16, label=True,
                  colorizer=random_colorizer, label_function=labeler):
    if ax is None:
        fig, ax = plt.subplots(1)
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
    ax.legend(bbox_to_anchor=(1.7, 1.), ncol=2, fontsize=10)
    ax.axes.spines['right'].set_visible(False)
    ax.axes.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.set_xlabel("Retention Time", fontsize=axis_label_font_size)
    ax.set_ylabel("Relative Abundance", fontsize=axis_label_font_size)
    [t.set(fontsize=axis_font_size) for t in ax.get_xticklabels()]
    [t.set(fontsize=axis_font_size) for t in ax.get_yticklabels()]
    return ax
