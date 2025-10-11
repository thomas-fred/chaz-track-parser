import random

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.figure import Figure


COASTLINE_URL = "https://naturalearth.s3.amazonaws.com/50m_physical/ne_50m_coastline.zip"


def plot_joint_distributions(
    df: gpd.GeoDataFrame,
    variables: list[str],
    title: str,
    stride: int = 5000,
    **kwargs: dict
) -> Figure:
    to_plot = df.iloc[::stride].loc[:, ["basin_id"] + variables]
    default_kwargs = dict(
        hue="basin_id",
        plot_kws={"alpha": 0.5},
        diag_kws={"fill": False, "alpha": 0.5},
        corner=True,
        kind="kde",
    )
    g = sns.pairplot(to_plot, **dict(default_kwargs, **kwargs))
    g.fig.suptitle(title)
    return g.fig


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

