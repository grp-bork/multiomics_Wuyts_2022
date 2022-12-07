#!/usr/bin/env python
# coding: utf-8

# from scipy import stats
# from statsmodels.stats.multitest import multipletests

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

fileFolder = "../data/ancom/"

meta16Stable = pd.read_csv(
    fileFolder + "ancom_FC_16S_time_conditions.txt.gz", delimiter=" ", index_col=1
)
meta16Stable.columns = ["meta16S_" + item for item in meta16Stable.columns]

metaGtable = pd.read_csv(
    fileFolder + "ancom_FC_metaG_time_conditions.txt.gz", delimiter=" ", index_col=1
)
metaGtable.columns = ["metaG_" + item for item in metaGtable.columns]

metaTtable = pd.read_csv(
    fileFolder + "ancom_FC_metaT_time_conditions.txt.gz", delimiter=" ", index_col=1
)
metaTtable.columns = ["metaT_" + item for item in metaTtable.columns]


metaPtable = pd.read_csv(
    fileFolder + "ancom_FC_metaP_time_conditions.txt.gz", delimiter=" ", index_col=1
)
metaPtable.columns = ["metaP_" + item for item in metaPtable.columns]


plotANCOM = pd.concat(
    [meta16Stable, metaGtable, metaTtable, metaPtable], axis=1, join="outer"
)

# add species names
fileFolder = "../data/relative_abundance_tables/"
infilename = fileFolder + "metagenomics_relabun.tsv.gz"
nameDF = pd.read_csv(infilename, delimiter="\t", index_col=0)
relabMetaGDF = nameDF[
    [
        item
        for item in nameDF.columns
        if not (("Inoculum" in item) or ("Transfer" in item) or ("species" in item))
    ]
]
nameDF = nameDF["species_name"]

# get species order in phylogenetic tree
treeDF = pd.read_csv(
    "../data/pathway_enrichment_analysis/species_tree_order_updated_extendedcolor.txt.gz",
    delimiter=",",
    index_col=0,
)
speciesColor_dict = dict(zip(treeDF.index, treeDF.color))

# get updated species order in phylogenetic tree
fileFolder = "../data/colours/"
fileName = "colour_scheme.tsv.gz"
colorDF = pd.read_csv(fileFolder + fileName, delimiter="\t", index_col=0)

speciesColor_dict = dict(zip(colorDF.index, colorDF.colour))

nameDF = pd.concat([treeDF, plotANCOM, relabMetaGDF], axis=1, join="outer")
# try to plot in the order of the tree
plotSORT = list(range(np.shape(nameDF)[0] - 1, -1, -1))

#################################################
# select data and sort according to metaG abundance
replic = ["A", "B"]

for repl in replic:
    sign_threshold = "detected_0.7"
    curcolumns = [
        item
        for item in nameDF.columns
        if (sign_threshold in item)
        and (
            (("_" + repl) in item[-len("_" + repl):])
            or (("_ " + repl) in item[-len("_ " + repl):])
        )
    ]

    # Make asterisk matrix to mark significantly different changes
    stardf = nameDF[curcolumns].copy()
    stardf.columns = [
        item.replace(sign_threshold, "signstar") for item in stardf.columns
    ]

    for col in range(np.shape(stardf)[1]):
        stardf.iloc[:, col] = [
            "*" if nameDF[curcolumns[col]].iloc[i] else " "
            for i in range(np.shape(stardf)[0])
        ]

    # add double and triple stars
    sign_threshold = "detected_0.8"
    curcolumns = [
        item
        for item in nameDF.columns
        if (sign_threshold in item)
        and (
            (("_" + repl) in item[-len("_" + repl):])
            or (("_ " + repl) in item[-len("_ " + repl):])
        )
    ]
    for col in range(np.shape(stardf)[1]):
        stardf.iloc[:, col] = [
            "**" if nameDF[curcolumns[col]].iloc[i] else stardf.iloc[i, col]
            for i in range(np.shape(stardf)[0])
        ]
    sign_threshold = "detected_0.9"
    curcolumns = [
        item
        for item in nameDF.columns
        if (sign_threshold in item)
        and (
            (("_" + repl) in item[-len("_" + repl):])
            or (("_ " + repl) in item[-len("_ " + repl):])
        )
    ]
    for col in range(np.shape(stardf)[1]):
        stardf.iloc[:, col] = [
            "***" if nameDF[curcolumns[col]].iloc[i] else stardf.iloc[i, col]
            for i in range(np.shape(stardf)[0])
        ]
    for col in range(np.shape(stardf)[1]):
        stardf.iloc[:, col] = [
            "nan"
            if type(nameDF[curcolumns[col]].iloc[i]) == float
            else stardf.iloc[i, col]
            for i in range(np.shape(stardf)[0])
        ]

    # Add star dataframe to nameDF
    nameDF = pd.concat([nameDF, stardf], axis=1)

    selectTime = "T48h"
    curcolumnsFC = [
        item
        for item in nameDF.columns
        if ("FCcondctr" in item)
        and (selectTime in item)
        and (
            (("_" + repl) in item[-len("_" + repl):])
            or (("_ " + repl) in item[-len("_ " + repl):])
        )
    ]

    curcolumnsFC = [
        item
        for item in nameDF.columns
        if ("FCcondctr" in item)
        and ("metaP" in item)
        and ("Niclo" in item)
        and (
            (("_" + repl) in item[-len("_" + repl):])
            or (("_ " + repl) in item[-len("_ " + repl):])
        )
    ]

    # calculate average FC per omic, replicate and condition
    replicates = ["A", "B", ""]
    conditions = ["Niclo", "Metfo", "Chlor"]
    omics = ["meta16S", "metaG", "metaT", "metaP"]
    for omic in omics:
        for cond in conditions:
            for currepl in replicates:
                curcolumnsFC = [
                    item
                    for item in nameDF.columns
                    if ("FCcondctr" in item)
                    and (omic in item)
                    and (cond in item)
                    and (
                        (("_" + currepl) in item[-len("_" + currepl):])
                        or (("_ " + currepl) in item[-len("_ " + currepl):])
                    )
                ]

                plotANCOMFC = nameDF[curcolumnsFC]
                nameDF[
                    omic + "_MeanFCcondctr_" + cond + "_" + currepl
                ] = plotANCOMFC.mean(axis=1)

    species_names_plot = nameDF.index + "-" + nameDF.species

    plt.rcParams.update({"font.size": 14})
    #################################################################
    fig = plt.figure(figsize=(15, 10))
    # set up subplot grid
    gridcol = 6
    gridspec.GridSpec(1, gridcol)
    #################################################################
    metaOdata = nameDF.copy()
    # select data and sort according to metaG abundance
    curcolumns = [
        item
        for item in metaOdata.columns
        if (("_r1" in item) or ("_r2" in item)) and ((repl + "_") in item)
    ]

    curdata = np.log10(metaOdata[curcolumns] * 100)
    curdata = curdata.clip(-5)
    # sort
    curdata = curdata.iloc[plotSORT, :]
    # linearize for violin plot
    res = pd.melt(curdata.assign(index=curdata.index), id_vars=["index"])

    spidx = 0
    axs = [0] * gridcol
    # equal subplots
    axs[spidx] = plt.subplot2grid((1, gridcol), (0, spidx), colspan=2, rowspan=1)
    sns.violinplot(
        x="value",
        y="index",
        data=res,
        scale="count",
        palette=speciesColor_dict,  # "Set3",
        ax=axs[spidx],
    )
    axs[spidx].plot([0, 0], [-1, np.shape(curdata)[0]], "--k")
    axs[spidx].set_ylim([-1, np.shape(curdata)[0]])
    axs[spidx].set_xlabel("Abundance across samples, log10")
    axs[spidx].set_ylabel("Species")
    axs[spidx].set_title("MetaG")
    curYlabels = species_names_plot[plotSORT]
    axs[spidx].set_yticks(np.arange(0, len(curYlabels), 1))
    axs[spidx].set_yticklabels(curYlabels)

    # plot per condition
    dataCondition_unique = ["Niclo", "Metfo", "Chlor"]
    spidx = spidx + 1
    for curcond in dataCondition_unique:
        curcolumnsFC = [
            item
            for item in nameDF.columns
            if ("MeanFCcondctr" in item)
            and (
                (("_" + repl) in item[-len("_" + repl):])
                or (("_ " + repl) in item[-len("_ " + repl):])
            )
            and (curcond in item)
        ]
        curcolumns = [
            item
            for item in nameDF.columns
            if ("signstar" in item)
            and (
                (("_" + repl) in item[-len("_" + repl):])
                or (("_ " + repl) in item[-len("_ " + repl):])
            )
            and (curcond in item)
        ]

        plotANCOMFC = nameDF[curcolumnsFC]
        plotANCOMFC = np.asarray(plotANCOMFC)

        plotANCOMstar = nameDF[curcolumns]

        plotANCOM = np.asarray(plotANCOMFC)

        omics_names = [item[0: item.find("_")] for item in curcolumns]

        axs[spidx] = plt.subplot2grid(
            (1, gridcol), (0, spidx + 1), colspan=1, rowspan=1
        )
        c = axs[spidx].pcolor(plotANCOM[plotSORT, :], cmap="RdBu_r")
        c.set_clim(-3, 3)
        axs[spidx].set_yticks(np.arange(0.5, plotANCOM.shape[0], 1))
        axs[spidx].set_yticklabels("")
        # axs[spidx].set_yticklabels(species_names_plot);
        axs[spidx].set_ylim([-0.5, np.shape(plotANCOM)[0] + 0.5])
        axs[spidx].set_xticks(np.arange(0.5, len(curcolumns), 1))
        axs[spidx].set_xticklabels(omics_names, rotation="vertical")
        axs[spidx].set_title(curcond)

        ###############################
        # add stars
        for i in range(plotANCOMstar.shape[0]):
            for j in range(plotANCOMstar.shape[1]):
                plt.text(j, plotANCOMstar.shape[0] - i - 1, plotANCOMstar.iloc[i, j])
        spidx = spidx + 1

    cbar = axs[spidx - 1].figure.colorbar(c)
    cbar.ax.set_ylabel("FC to control at T=48h", rotation=-90, va="bottom")
    fig.suptitle("ANCOM analysis for replicate " + repl)
    # fig.tight_layout()

    # save figure to file
    plt.rcParams["svg.fonttype"] = "none"
    plt.savefig(
        "Figure_3_B" + str(replic.index(repl) + 1) + ".svg",
        format="svg",
        bbox_inches="tight",
    )
