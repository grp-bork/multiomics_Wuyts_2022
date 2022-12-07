#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

fileFolder = "../data/differential_expression/"
metaTfile = "DE_genes_and_proteins_unfiltered.tsv.gz"
metaTdata = pd.read_csv(fileFolder + metaTfile, sep="\t", index_col=0, nrows=10)

print("Read sample DE data")

metacolumns = [
    item
    for item in metaTdata.columns
    if (item == "omic")
    or ("drug" in item)
    or ("timepoint" in item)
    or ("gene" in item)
    or ("log2FoldChange" in item)
    or ("pvalue" in item)
    or ("padj" in item)
    or ("accession_number" in item)
    or ("EC" in item)
    or ("KEGG_ko" in item)
    or ("KEGG_Pathway" in item)
    or ("KEGG_Reaction" in item)
    or ("eggNOG_free_text_desc" in item)
    or ("description" in item)
]

metaTdata = pd.read_csv(fileFolder + metaTfile, sep="\t")  # ,index_col=0)

print("Read all DE data. Beginning data transformations")

metaTdata[
    (metaTdata["product"].str.contains("AcrB"))
    & (metaTdata["drug"].str.contains("chlor"))
    & (metaTdata["timepoint"].str.contains("180"))
]

metaTdata[
    (metaTdata["product"].str.contains("efflux pump"))
    & (metaTdata["drug"].str.contains("chlor"))
    & (metaTdata["timepoint"].str.contains("180"))
]

# calculate number of changes per drug per timepoint per omic
omic_unique = list(set(metaTdata["omic"]))
drug_unique = list(set(metaTdata["drug"]))
time_unique = list(set(metaTdata["timepoint"]))

time_unique = [
    "timepoint15",
    "timepoint30",
    "timepoint60",
    "timepoint180",
    "timepoint2880",
    "timepoint5760",
]

# read gene and protein DE results
fileFolder = "../data/ranova/"
metaPfile = "FC.metaP.tsv.gz"
metaPdata = pd.read_csv(fileFolder + metaPfile, sep="\t", index_col=0)


# get p-values
fileFolder = "../data/ranova/"
metaPfile = "ranova.metaP.runAB.post_hoc_ttest.tsv.gz"  # 'FC.metaP.tsv.gz'
metaPPvaluesdata = pd.read_csv(fileFolder + metaPfile, sep="\t", index_col=0)

# add p-values
metaPdata["prot_drug_time"] = (
    metaPdata["protein"] + "_" + metaPdata["drug"] + "_" + metaPdata["time"]
)

metaPPvaluesdata["prot_drug_time"] = (
    metaPPvaluesdata["protein"]
    + "_"
    + metaPPvaluesdata["drug"]
    + "_"
    + metaPPvaluesdata["time"]
)
metaPPvaluesdata.iloc[0:1, :]

metaPPvaluesdata = metaPPvaluesdata.set_index("prot_drug_time")
metaPdata = metaPdata.set_index("prot_drug_time")

metaPdata["prot_padj"] = 1
metaPdata.loc[metaPPvaluesdata.index, "prot_padj"] = metaPPvaluesdata["p_adj"]

metaPdata = metaPdata.reset_index(drop=True)

time_unique = list(set(metaPdata["time"]))
time_unique = ["T15", "T30", "T1h", "T3h", "T48h", "T96h"]

metaGdata = metaTdata[metaTdata["omic"] == "metagenomics"].copy()
metaTdata = metaTdata[metaTdata["omic"] == "metatranscriptomics"]

timeConversion = {
    "timepoint15": "T15",
    "timepoint30": "T30",
    "timepoint60": "T1h",
    "timepoint180": "T3h",
    "timepoint2880": "T48h",
    "timepoint5760": "T96h",
}

timeSort = {"T15": 1, "T30": 2, "T1h": 3, "T3h": 4, "T48h": 5, "T96h": 6}

drugConversion = {
    "chlorpromazine": "Chlor",
    "metformin": "Metfo",
    "niclosamide": "Niclo",
}

metaTdata_time = [timeConversion[item] for item in metaTdata["timepoint"]]

metaTdata["timepoint"] = metaTdata_time

metaGdata_time = [timeConversion[item] for item in metaGdata["timepoint"]]
metaGdata["timepoint"] = metaGdata_time

metaTdata_drug = [drugConversion[item] for item in metaTdata["drug"]]

metaTdata["drug"] = metaTdata_drug


metaGdata_drug = [drugConversion[item] for item in metaGdata["drug"]]
metaGdata["drug"] = metaGdata_drug

metaTdata["gene_drug_time"] = (
    metaTdata["gene"] + "_" + metaTdata["drug"] + "_" + metaTdata["timepoint"]
)

metaGdata["gene_drug_time"] = (
    metaGdata["gene"] + "_" + metaGdata["drug"] + "_" + metaGdata["timepoint"]
)

# join metaP
metaPdata["gene_drug_time"] = (
    metaPdata["protein"] + "_" + metaPdata["drug"] + "_" + metaPdata["time"]
)

metaPTdata = metaPdata.copy()

metaPTdata = metaPTdata.merge(metaTdata, on="gene_drug_time", how="inner")

metaPTdata["timesort"] = [timeSort[x] for x in metaPTdata["time"]]

# leave only joint comparison
metaPTdata = metaPTdata[metaPTdata["run"] == "runAB"].copy()

plotdata_colorby = metaPTdata["timepoint"]

print("Finished data transformations")

categories = ["T15", "T30", "T1h", "T3h", "T48h", "T96h"]
pthres = 0.1
fcthres = np.log2(1.5)

fig = plt.figure(figsize=(30, 30))
plt.rcParams.update({"font.size": 10})
pltidx = 1
PCCcorrmat = np.zeros([len(categories), len(categories)])
SCCcorrmat = np.zeros([len(categories), len(categories)])

for i in range(len(categories)):
    for j in range(len(categories)):
        curtime1 = categories[i]
        curtime2 = categories[j]
        plotdataX = np.asarray(
            metaPTdata[metaPTdata["time"] == curtime2]["log2FoldChange"]
        )
        plotdataY = np.asarray(metaPTdata[metaPTdata["time"] == curtime1]["log2fc"])
        # remove nan and also remove zeros
        plotdatanan = np.asarray(
            np.isfinite(plotdataX)
            & np.isfinite(plotdataY)
            & (plotdataX != 0)
            & (plotdataY != 0)
        )
        plotdataX = plotdataX[plotdatanan.nonzero()]
        plotdataY = plotdataY[plotdatanan.nonzero()]

        # pearson
        curcorr = np.corrcoef(plotdataX, plotdataY)
        curcorr = curcorr[0, 1]
        PCCcorrmat[i, j] = curcorr

        # spearman
        curcorrSCC = stats.spearmanr(plotdataX, plotdataY)
        curcorrSCC = curcorrSCC[0]
        SCCcorrmat[i, j] = curcorrSCC

        pltidx = (len(categories) - 1 - i) * len(categories) + j + 1
        plt.subplot(len(categories), len(categories), pltidx)
        scatter = plt.scatter(plotdataX, plotdataY, c="k")
        # produce a legend with the unique colors from the scatter
        plt.xlabel(curtime2)  # + ' metaT')
        plt.ylabel(curtime1)  # + ' metaP')
        plt.title(
            "PCC=" + "{:.2f}".format(curcorr) + " SCC=" + "{:.2f}".format(curcorrSCC)
        )

        # add significant
        plotdataXp = np.asarray(metaPTdata[(metaPTdata["time"] == curtime2)]["padj"])
        plotdataYp = np.asarray(
            metaPTdata[(metaPTdata["time"] == curtime1)]["prot_padj"]
        )
        plotdataXp = plotdataXp[plotdatanan.nonzero()]
        plotdataYp = plotdataYp[plotdatanan.nonzero()]

        plot_significant = (
            (plotdataXp <= pthres)
            & (plotdataYp <= pthres)
            & (abs(plotdataX) >= fcthres)
            & (abs(plotdataY) >= fcthres)
        )

        scatter = plt.scatter(
            plotdataX[plot_significant], plotdataY[plot_significant], c="#009E73"
        )

plt.rcParams["svg.fonttype"] = "none"

plt.savefig("Supp_Figure_11.svg", format="svg", bbox_inches="tight")
