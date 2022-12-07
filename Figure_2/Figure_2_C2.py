#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns


# Plot in the order of the fraction of covered genes
repl = ""  # 'B'#'A'#

treeDF = pd.read_csv(
    "../data/pathway_enrichment_analysis/species_coverage_order_updated_extendedcolor.txt.gz",
    delimiter=",",
    index_col=0,
)
speciesColor_dict = dict(zip(treeDF.index, treeDF.color))

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

nameDF = pd.concat([treeDF, relabMetaGDF], axis=1, join="outer")

species_names_plot = nameDF.index + "-" + nameDF.species

# try to plot in the order of the tree
plotSORT = list(range(np.shape(nameDF)[0] - 1, -1, -1))
# plotSORT = list(range(np.shape(nameDF)[0]))

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

curdata = np.log10(metaOdata[curcolumns])
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
    color="#d4d4d4",  # palette=speciesColor_dict, #"Set3",
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

# save figure to file
figfileName = "Figure_2_C2"
plt.rcParams["svg.fonttype"] = "none"
plt.savefig(figfileName + repl + ".svg", format="svg", bbox_inches="tight")
