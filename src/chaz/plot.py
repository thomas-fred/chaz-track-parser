import random

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.figure import Figure


COASTLINE_URL = "https://naturalearth.s3.amazonaws.com/50m_physical/ne_50m_coastline.zip"


# For example: {"max_wind_speed_ms": ((10, 120), "log")}
type PlotMetadata = dict[str, tuple[tuple[float, float], str]]


def pairplot(df: pd.DataFrame, metadata: PlotMetadata, **kwargs: dict) -> sns.PairGrid:
    """
    Extend seaborn's pairplot to enforce consistent axis limits and optional log scaling.

    Args:
        df: Table of data to plot, must include columns referenced in metadata.
        metadata: Mapping from column names to tuple of limits (min, max) and
            scaling, i.e. "linear" or "log"
        kwargs: Keyword args to pass through to seaborn.pairplot

    Returns:
        Figure object.
    """
    # We don't use the index, but 10 functions down pandas will
    # Reindex to avoid: ValueError("cannot reindex on an axis with duplicate labels")
    g = sns.pairplot(df.reset_index(drop=True), **kwargs)
    
    for i, row_var in enumerate(metadata.keys()):
        for j, col_var in enumerate(metadata.keys()):
            if not g.axes[i, j]:
                continue
            elif i == j:
                limits, scale = metadata[row_var]
                g.axes[i, j].set_xlim(*limits)
                g.axes[i, j].set_xscale(scale)
            else:
                limits, scale = metadata[col_var]
                g.axes[i, j].set_xlim(*limits)
                g.axes[i, j].set_xscale(scale)
                limits, scale = metadata[row_var]
                g.axes[i, j].set_ylim(*limits)
                g.axes[i, j].set_yscale(scale)
                
    return g


def plot_joint_distributions(
    df: gpd.GeoDataFrame,
    metadata: PlotMetadata,
    title: str,
    stride: int = 5000,
    **kwargs: dict
) -> Figure:
    to_plot = df.iloc[::stride].loc[:, ["basin_id"] + list(metadata.keys())]
    default_kwargs = dict(
        hue="basin_id",
        plot_kws={"alpha": 0.5},
        diag_kws={"fill": False, "alpha": 0.5},
        corner=True,
        kind="kde",
    )
    g = pairplot(to_plot, metadata, **dict(default_kwargs, **kwargs))
    g.figure.suptitle(title)
    return g.figure


def plot_tc_frequency(df: gpd.GeoDataFrame, title: str, **kwargs: dict) -> Figure:
    f, ax = plt.subplots()
    to_plot = df.loc[:, ["year", "track_id"]].groupby("year").nunique()
    ax = sns.histplot(to_plot, legend=False, discrete=True, **kwargs)
    ax.grid(alpha=0.2)
    ax.set_xlabel("Annual TC frequency")
    ax.set_title(title)
    return f


def plot_tc_frequency_per_basin(df: gpd.GeoDataFrame, title: str, **kwargs: dict) -> Figure:
    f, ax = plt.subplots()
    to_plot = (
        df.loc[:, ["basin_id", "year", "track_id"]]
            .groupby(["basin_id", "year"])["track_id"]
            .nunique().reset_index()
    )
    ax = sns.histplot(to_plot, x="track_id", ax=ax, discrete=True, hue="basin_id", **kwargs)
    ax.grid(alpha=0.2)
    ax.set_xlabel("Annual TC frequency")
    ax.set_title(title)
    return f


def plot_scatter_map(df: gpd.GeoDataFrame, var: str, label: str, title: str, max_points: int = 200_000) -> Figure:
    f, ax = plt.subplots(figsize=(12, 5))

    # Subsample the input to show approximately max_points (of complete tracks)
    mean_points_per_track = df.loc[:, ["track_id", "timestep"]].groupby("track_id").count().mean()
    n_tracks_to_plot = int(round(max_points / mean_points_per_track, 0))
    tracks_to_plot = random.sample(list(df.track_id.unique()), n_tracks_to_plot)
    to_plot = df[df.track_id.isin(tracks_to_plot)].copy()

    gpd.read_file(COASTLINE_URL).plot(ax=ax, alpha=0.4, lw=1, color="k")

    long = to_plot.geometry.x
    to_plot["geometry"] = gpd.points_from_xy(np.where(long > 180, long - 360, long), to_plot.geometry.y)

    to_plot = to_plot.sort_values(var, ascending=True)
    to_plot.plot(
        var,
        vmin=40,
        vmax=120,
        ax=ax,
        s=4,
        alpha=0.5,
        legend=True,
        cmap="inferno_r",
        legend_kwds={"shrink":0.56, "label": label}
    )

    ax.set_xticks(np.linspace(-180, 180, 13))
    ax.set_yticks(np.linspace(-60, 60, 5))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-65, 65)
    ax.grid(alpha=0.2)
    ax.set_title(title)
    plt.subplots_adjust(left=0.15, right=0.85)

    return f

