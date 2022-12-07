#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# read gene and protein DE results
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

metaGTdata = metaGdata.copy()
metaGTdata = metaGTdata.merge(metaTdata, on="gene_drug_time", how="inner")
metaGTdata["timesort"] = [timeSort[x] for x in metaGTdata["timepoint_x"]]

metaGTdata = metaGTdata[metaGTdata["log2FoldChange_x"].notna()]
metaGTdata = metaGTdata[metaGTdata["log2FoldChange_y"].notna()]

plotdataX = metaPTdata["log2FoldChange"]
plotdataY = metaPTdata["log2fc"]

plotdataX = plotdataX.fillna(0)
plotdataY = plotdataY.fillna(0)

curcorr = np.corrcoef(np.ma.masked_invalid(plotdataX), np.ma.masked_invalid(plotdataY))
curcorr[0, 1]

plotdata_colorby = metaPTdata["timepoint"]

print("Finished data transformations")

categories = np.unique(plotdata_colorby)
categories = ["T15", "T30", "T1h", "T3h", "T48h", "T96h"]
colors = np.linspace(0, 1, len(categories))
colordict = dict(zip(categories, colors))

# High positive correlation between metaG and metaT fold changes, but we only consider times 3h, 48h and 96h, where majority of changes are driven in late time points where there are big changes in species abundances.
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
        # print(curcorr)
        curcorr = curcorr[0, 1]
        PCCcorrmat[i, j] = curcorr

        # spearman
        curcorrSCC = stats.spearmanr(plotdataX, plotdataY)
        # print(curcorrSCC)
        curcorrSCC = curcorrSCC[0]
        SCCcorrmat[i, j] = curcorrSCC


# Calculate number of changing transcripts and proteins per species, per time and per drug
species_names_unique = list(set(metaPTdata.species_name))

drug_names_unique = list(set(metaPTdata.drug_x))

column_names_species = []
for drug in drug_names_unique:
    for time_i in categories:
        column_names_species.append(drug + "_" + time_i + "_metaT_up")
        column_names_species.append(drug + "_" + time_i + "_metaT_down")
        column_names_species.append(drug + "_" + time_i + "_metaP_up")
        column_names_species.append(drug + "_" + time_i + "_metaP_down")

species_changes_df = pd.DataFrame(
    index=species_names_unique, columns=column_names_species
)

metaT_fcthres = np.log2(2)
metaP_fcthres = np.log2(1.5)
padj_thres = 0.05

for species in species_names_unique:
    for drug in drug_names_unique:
        for time_i in categories:
            metaTPtemp = metaPTdata[
                (metaPTdata["species_name"] == species)
                & (metaPTdata["drug_x"] == drug)
                & (metaPTdata["time"] == time_i)
            ]
            species_changes_df.loc[species, drug + "_" + time_i + "_metaT_up"] = sum(
                (metaTPtemp["log2FoldChange"] >= metaT_fcthres)
                & (metaTPtemp["padj"] <= padj_thres)
            )
            species_changes_df.loc[species, drug + "_" + time_i + "_metaT_down"] = sum(
                (metaTPtemp["log2FoldChange"] <= -metaT_fcthres)
                & (metaTPtemp["padj"] <= padj_thres)
            )
            species_changes_df.loc[species, drug + "_" + time_i + "_metaP_up"] = sum(
                (metaTPtemp["log2fc"] >= metaP_fcthres)
                & (metaTPtemp["prot_padj"] <= padj_thres)
            )
            species_changes_df.loc[species, drug + "_" + time_i + "_metaP_down"] = sum(
                (metaTPtemp["log2fc"] <= -metaP_fcthres)
                & (metaTPtemp["prot_padj"] <= padj_thres)
            )


fig = plt.figure(figsize=(8, 5))
plotdata = SCCcorrmat
# Plot corr matrix
ax = sns.heatmap(plotdata, vmin=-0.5, vmax=0.5, cmap="vlag", annot=True)
ax.invert_yaxis()
ax.set_xticklabels(categories)
ax.set_yticklabels(categories)
ax.set_xlabel("metaT log2FC")
ax.set_ylabel("metaP log2FC")
ax.set_title("Spearman corr")

filename = "Figure_4_C"

plt.rcParams["svg.fonttype"] = "none"

plt.savefig(filename + ".svg", format="svg", bbox_inches="tight")

# plot with matpltlib for export
fig = plt.figure(figsize=(10, 10))
plt.rcParams.update({"font.size": 10})

curtime1 = categories[0]
curtime2 = categories[2]
dfP = metaPTdata[
    (metaPTdata["time"] == curtime2) & (metaPTdata["drug_x"] == "Chlor")
].copy()
dfP = dfP.set_index("protein")
# get metaT of the correct time
dfT = metaPTdata[
    (metaPTdata["time"] == curtime1) & (metaPTdata["drug_x"] == "Chlor")
].copy()
dfT = dfT.set_index("protein")
dfP.loc[dfP.index, "log2FoldChange"] = dfT.loc[dfP.index, "log2FoldChange"]
dfP.loc[dfP.index, "padj"] = dfT.loc[dfP.index, "padj"]

df = dfP.copy()
df["protein"] = df.index
df["prot_product"] = df["protein"] + " " + df["product"]

df["signif"] = "none"
df.loc[df["prot_padj"] <= 0.05, "signif"] = "metaP"
df.loc[df["padj"] <= 0.05, "signif"] = "metaT"
df.loc[(df["prot_padj"] <= 0.05) & (df["padj"] <= 0.05), "signif"] = "both"

df.loc[df["log2FoldChange"] > 10, "log2FoldChange"] = 8
df.loc[df["log2FoldChange"] < -10, "log2FoldChange"] = -8


plotdataX = np.asarray(df["log2FoldChange"])
plotdataY = np.asarray(df["log2fc"])

plotdatanan = np.asarray(
    np.isfinite(plotdataX)
    & np.isfinite(plotdataY)
    & (plotdataX != 0)
    & (plotdataY != 0)
)
plotdataX = plotdataX[plotdatanan.nonzero()]
plotdataY = plotdataY[plotdatanan.nonzero()]

scatter = plt.scatter(plotdataX, plotdataY, c="k")
plt.ylim([-8, 8])
plt.xlim([-8, 8])

plt.xlabel("MetaT log2FC " + curtime1)
plt.ylabel("MetaP log2FC " + curtime2)

curcorrSCC = stats.spearmanr(plotdataX, plotdataY)
# print(curcorrSCC)
curcorrSCC = curcorrSCC[0]

plt.text(4, -6, "rho={:.2f}".format(curcorrSCC))

# scatter = plt.scatter(df.loc[df['Bacteroides']=='bact','log2FoldChange'],
#                       df.loc[df['Bacteroides']=='bact','log2fc'],
#                       c='skyblue')#, c=colordict[curtime])

plt.rcParams["svg.fonttype"] = "none"

filename = "Figure_4_D"
plt.savefig(filename + ".svg", format="svg", bbox_inches="tight")
