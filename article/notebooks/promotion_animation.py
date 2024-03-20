import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines
import seaborn as sns
import numpy as np
import pandas as pd
from patsy import dmatrices
import statsmodels.api as sm
from matplotlib import animation
from functools import partial
from typing import Tuple, List

from pseudobatch import pseudobatch_transform, pseudobatch_transform_pandas
from pseudobatch.datasets._dataloaders import _prepare_simulated_dataset


def animate_line(i, line, xdata: pd.Series, ydata: pd.Series):
    line.set_data(
        xdata.iloc[:i],
        ydata.iloc[:i],
    )

    return (line,)


###############################################################
##### Prepare the data and fit the growth rate ################
###############################################################
def prepare_data():
    # Here we do not use the datasets stored in the package.
    # Instead we import the dataset from the article/data folder,
    # this makes sure if the simulations are rerun the new data is used
    data_path = pathlib.Path("article/data/standard_fed-batch_process.csv")
    fedbatch_df = _prepare_simulated_dataset(data_path)

    FIGURES_DIR = pathlib.Path("article/figures")

    # pseudo batch transform
    fedbatch_df["c_Biomass_pseudo"] = pseudobatch_transform(
        measured_concentration=fedbatch_df["c_Biomass"].to_numpy(),
        reactor_volume=fedbatch_df["v_Volume"].to_numpy(),
        accumulated_feed=fedbatch_df["v_Feed_accum"].to_numpy(),
        concentration_in_feed=0,
        sample_volume=fedbatch_df["sample_volume"].fillna(0).to_numpy(),
    )
    fedbatch_df["m_Biomass_pseudo"] = (
        fedbatch_df["c_Biomass_pseudo"] * fedbatch_df["v_Volume"].iloc[0]
    )

    return fedbatch_df.copy().iloc[::10]


def init(lines: List[mpl.lines.Line2D], ax) -> List[mpl.lines.Line2D]:
    for line in lines:
        line.set_data([], [])
    # remove the text
    ax.texts[-1].set_text("")
    return lines


def main():
    plot_data = prepare_data()

    ## Fit growth rate
    y = plot_data["m_Biomass"].transform(np.log)
    X_noncorrected = sm.add_constant(plot_data["timestamp"])
    model_noncorrected = sm.OLS(endog=y, exog=X_noncorrected)
    res_noncorrected = model_noncorrected.fit()
    x_pred = X_noncorrected.iloc[::10]
    y_pred_noncorrected = res_noncorrected.predict(x_pred).transform(np.exp)

    # create a figure and axis
    fig, ax = plt.subplots()
    ax.set_xlim(plot_data["timestamp"].min(), plot_data["timestamp"].max())
    ax.set_ylim(
        plot_data["m_Biomass_pseudo"].min(), plot_data["m_Biomass_pseudo"].max()
    )
    ax.semilogy()
    # ax.set_title("Fed-batch culture with sampling")

    # create a line for the raw data
    (line_raw,) = ax.plot(
        [], [], color="C0", marker=None, label="Total biomass"
    )
    (line_pred_noncorrected,) = ax.plot(
        [],
        [],
        color="grey",
        linestyle="--",
        marker=None,
        label="Exponential \ngrowth rate fit",
    )
    (line_pseudo,) = ax.plot(
        [],
        [],
        color="C1",
        marker=None,
        label="Pseudo batch \ntranformed \nbiomass",
    )

    ## Setup explanatory text
    ax.text(
        0.01,
        1.05,
        "" + "",
        fontdict={"fontsize": 14},
        horizontalalignment="left",
        verticalalignment="bottom",
        transform=ax.transAxes,
    )

    ### Configure legend
    # shrink the axis to make space for the legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height * 0.8])
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

    def animate(i, scene_splits: list[int]):
        if i <= scene_splits[0]:
            animate_line(
                i, line_raw, plot_data["timestamp"], plot_data["m_Biomass"]
            )
        elif i <= scene_splits[1]:
            # waiting a bit
            pass
        elif i <= scene_splits[2]:
            animate_line(
                i - scene_splits[0],
                line_pred_noncorrected,
                x_pred.iloc[:, 1],
                y_pred_noncorrected,
            )
        elif i < scene_splits[3]:
            pass
        elif i < scene_splits[4]:
            line_pseudo.set_data(
                plot_data["timestamp"], plot_data["m_Biomass_pseudo"]
            )
        if i == 0:
            ax.texts[-1].set_text(
                "When a fedbatch culture is sampled the \n"
                + "biomass curve becomes discontinuous",
            )
        if i == scene_splits[2]:
            ax.texts[-1].set_text(
                "This leads to an INCORRECT estimate of\n"
                + "the exponential growth rate"
            )
        return None

    # call the animator.  blit=True means only re-draw the parts that have changed.
    scenes = [len(plot_data["timestamp"]), 5, len(x_pred), 50, 100]
    n_frames = sum(scenes)
    scene_splits_list = list(np.cumsum(scenes))

    animate = partial(animate, scene_splits=scene_splits_list)
    init_anim = partial(
        init,
        lines=[line_raw, line_pseudo, line_pred_noncorrected],
        ax=ax,
    )
    anim = animation.FuncAnimation(
        fig,
        animate,
        frames=n_frames,
        init_func=init_anim,
        interval=10,
    )
    print("Saving animation")
    plt.show()
    # anim.save(FIGURES_DIR / "animation.gif", fps=10)


if __name__ == "__main__":
    main()
