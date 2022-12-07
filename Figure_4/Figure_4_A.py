import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


fileFolder = "../data/ranova/"
fileName = "ranova.separate.metaB.log10.runAB.tsv.gz"
metaB_anova_df = pd.read_csv(fileFolder + fileName, sep="\t", index_col=0)

# get fold change per time point
fileName = "ranova.metaB.runAB.post_hoc_ttest.tsv.gz"
metaB_posthoc_df = pd.read_csv(fileFolder + fileName, sep="\t", index_col=0)

# save max fold change per ion per drug
metaB_anova_df["maxlog2FCtime"] = 0
for i in range(len(metaB_anova_df)):
    curfc = metaB_posthoc_df[
        (metaB_posthoc_df["m/z"] == metaB_anova_df.index[i])
        & (metaB_posthoc_df["drug"] == metaB_anova_df["drug"].iloc[i])
    ]
    curfc = curfc["log2fc"].values
    metaB_anova_df.iloc[i, list(metaB_anova_df.columns).index("maxlog2FCtime")] = max(
        np.nanmax(curfc), np.nanmin(curfc), key=abs
    )

drug_unique = list(set(metaB_posthoc_df["drug"]))
time_unique = list(set(metaB_posthoc_df["time"]))

metaBnumDOWN = pd.DataFrame(columns=drug_unique, index=time_unique)
metaBnumUP = pd.DataFrame(columns=drug_unique, index=time_unique)

fcThreshold = np.log2(1.5)
pThreshold = 0.05

for curdrug in drug_unique:
    for curtime in time_unique:
        cursubset_df = set(
            metaB_anova_df[
                (metaB_anova_df["drug"] == curdrug)
                & (metaB_anova_df["ihw.fdr"] <= pThreshold)
            ].index
        )
        curfoldchange = metaB_posthoc_df[
            (metaB_posthoc_df["drug"] == curdrug)
            & (metaB_posthoc_df["time"] == curtime)
            & (metaB_posthoc_df["log2fc"] >= fcThreshold)
        ]

        metaBnumUP.loc[curtime, curdrug] = len(
            cursubset_df.intersection(curfoldchange["m/z"])
        )
        curfoldchange = metaB_posthoc_df[
            (metaB_posthoc_df["drug"] == curdrug)
            & (metaB_posthoc_df["time"] == curtime)
            & (metaB_posthoc_df["log2fc"] <= -fcThreshold)
        ]
        metaBnumDOWN.loc[curtime, curdrug] = len(
            cursubset_df.intersection(curfoldchange["m/z"])
        )

metaPnumDOWN = pd.DataFrame(columns=drug_unique, index=time_unique)
metaPnumUP = pd.DataFrame(columns=drug_unique, index=time_unique)

fileFolder = "../data/ranova/"
fileName = "ranova.separate.metaP.log10.runAB.tsv.gz"
metaP_anova_df = pd.read_csv(fileFolder + fileName, sep="\t", index_col=0)

# get fold change per time point
fileName = "ranova.metaP.runAB.post_hoc_ttest.tsv.gz"
metaP_posthoc_df = pd.read_csv(fileFolder + fileName, sep="\t", index_col=0)

# save max fold change per ion per drug
metaP_anova_df["maxlog2FCtime"] = 0
for i in range(len(metaP_anova_df)):
    curfc = metaP_posthoc_df[
        (metaP_posthoc_df["protein"] == metaP_anova_df.index[i])
        & (metaP_posthoc_df["drug"] == metaP_anova_df["drug"].iloc[i])
    ]
    curfc = curfc["log2fc"].values
    metaP_anova_df.iloc[i, list(metaP_anova_df.columns).index("maxlog2FCtime")] = max(
        np.nanmax(curfc), np.nanmin(curfc), key=abs
    )

fcThreshold = np.log2(1.5)
pThreshold = 0.05

for curdrug in drug_unique:
    for curtime in time_unique:
        cursubset_df = set(
            metaP_anova_df[
                (metaP_anova_df["drug"] == curdrug)
                & (metaP_anova_df["ihw.fdr"] <= pThreshold)
            ].index
        )
        curfoldchange = metaP_posthoc_df[
            (metaP_posthoc_df["drug"] == curdrug)
            & (metaP_posthoc_df["time"] == curtime)
            & (metaP_posthoc_df["log2fc"] >= fcThreshold)
        ]

        metaPnumUP.loc[curtime, curdrug] = len(
            cursubset_df.intersection(curfoldchange["protein"])
        )
        curfoldchange = metaP_posthoc_df[
            (metaP_posthoc_df["drug"] == curdrug)
            & (metaP_posthoc_df["time"] == curtime)
            & (metaP_posthoc_df["log2fc"] <= -fcThreshold)
        ]
        metaPnumDOWN.loc[curtime, curdrug] = len(
            cursubset_df.intersection(curfoldchange["protein"])
        )

fileFolder = "../data/differential_expression/"
fileName = "DE_genes_LRT_unfiltered.tsv.gz"
metaT_anova_df = pd.read_csv(fileFolder + fileName, sep="\t", index_col=0)

metaT_anova_df = metaT_anova_df[metaT_anova_df["omic"] == "metatranscriptomics"]

fileName = "DE_genes_and_proteins_unfiltered.tsv.gz"
metaT_posthoc_df = pd.read_csv(fileFolder + fileName, sep="\t", index_col=0)

metaT_posthoc_df = metaT_posthoc_df.loc[
    "metatranscriptomics",
]

drug_unique_t = list(set(metaT_anova_df.index))
time_unique_t = list(set(metaT_posthoc_df["timepoint"]))

# calculate number of changes
metaTnumDOWN = pd.DataFrame(columns=drug_unique_t, index=time_unique_t)
metaTnumUP = pd.DataFrame(columns=drug_unique_t, index=time_unique_t)

fcThreshold = np.log2(4)
pThreshold = 0.001

for curdrug in drug_unique_t:
    for curtime in time_unique_t:
        cursubset_df = set(
            metaT_anova_df[
                (metaT_anova_df.index == curdrug)
                & (metaT_anova_df["padj"] <= pThreshold)
            ]["gene"]
        )
        curfoldchange = metaT_posthoc_df[
            (metaT_posthoc_df["drug"] == curdrug)
            & (metaT_posthoc_df["timepoint"] == curtime)
            & (metaT_posthoc_df["log2FoldChange"] >= fcThreshold)
        ]

        metaTnumUP.loc[curtime, curdrug] = len(
            cursubset_df.intersection(curfoldchange["gene"])
        )
        curfoldchange = metaT_posthoc_df[
            (metaT_posthoc_df["drug"] == curdrug)
            & (metaT_posthoc_df["timepoint"] == curtime)
            & (metaT_posthoc_df["log2FoldChange"] <= -fcThreshold)
        ]
        metaTnumDOWN.loc[curtime, curdrug] = len(
            cursubset_df.intersection(curfoldchange["gene"])
        )

# rename and sort drugs and time points
time_dict = {
    "timepoint2880": "T48h",
    "timepoint180": "T3h",
    "timepoint5760": "T96h",
    "timepoint15": "T15",
    "timepoint60": "T1h",
    "timepoint30": "T30",
}

drug_dict = {"chlorpromazine": "Chlor", "niclosamide": "Niclo", "metformin": "Metfo"}

metaTnumDOWN.columns = [drug_dict[item] for item in metaTnumDOWN.columns]
metaTnumDOWN.index = [time_dict[item] for item in metaTnumDOWN.index]

metaTnumUP.columns = [drug_dict[item] for item in metaTnumUP.columns]
metaTnumUP.index = [time_dict[item] for item in metaTnumUP.index]

drug_unique = ["Chlor", "Metfo", "Niclo"]
time_unique = ["T15", "T30", "T1h", "T3h", "T48h", "T96h"]

fig = plt.figure(figsize=(20, 20))
plt.rcParams.update({"font.size": 20})

spidx = 1
for curdrug in drug_unique:
    curchangesDOWN = metaTnumDOWN.loc[time_unique, curdrug]
    curchangesUP = metaTnumUP.loc[time_unique, curdrug]
    ax = plt.subplot(3, 3, spidx)
    plt.plot(range(-2, len(time_unique) + 1), np.zeros([len(time_unique) + 3, 1]), "k")

    plt.bar(range(len(time_unique)), [-x for x in curchangesDOWN], color="#B30000")

    # time_unique
    plt.bar(range(len(time_unique)), curchangesUP, color="#B30000")
    plt.gca().invert_yaxis()
    plt.title(curdrug)
    plt.ylim([-15000, 15000])
    plt.xticks(range(-1, len(time_unique)), ["T0"] + time_unique)
    if spidx == 1:
        plt.ylabel("Metatranscriptomics\nDown     |     Up   ")
    else:
        plt.yticks([], [])
    plt.xlim([-1.5, len(time_unique) - 0.5])
    spidx = spidx + 1

for curdrug in drug_unique:
    curchangesDOWN = metaPnumDOWN.loc[time_unique, curdrug]
    curchangesUP = metaPnumUP.loc[time_unique, curdrug]
    ax = plt.subplot(3, 3, spidx)
    plt.plot(range(-2, len(time_unique) + 1), np.zeros([len(time_unique) + 3, 1]), "k")
    plt.bar(range(len(time_unique)), [-x for x in curchangesDOWN], color="#006D2C")
    # time_unique
    plt.bar(range(len(time_unique)), curchangesUP, color="#006D2C")
    plt.gca().invert_yaxis()
    plt.title(curdrug)
    plt.ylim([-200, 200])
    plt.xticks(range(-1, len(time_unique)), ["T0"] + time_unique)
    if spidx == 4:
        plt.ylabel("Metaproteomics\nDown     |     Up   ")
    else:
        plt.yticks([], [])
    plt.xlim([-1.5, len(time_unique) - 0.5])
    spidx = spidx + 1

for curdrug in drug_unique:
    curchangesDOWN = metaBnumDOWN.loc[time_unique, curdrug]
    curchangesUP = metaBnumUP.loc[time_unique, curdrug]
    ax = plt.subplot(3, 3, spidx)
    plt.plot(range(-2, len(time_unique) + 1), np.zeros([len(time_unique) + 3, 1]), "k")

    plt.bar(range(len(time_unique)), [-x for x in curchangesDOWN], color="#984EA3")
    # time_unique
    plt.bar(range(len(time_unique)), curchangesUP, color="#984EA3")
    plt.gca().invert_yaxis()
    plt.title(curdrug)
    plt.ylim([-200, 200])
    plt.xticks(range(-1, len(time_unique)), ["T0"] + time_unique)
    if spidx == 7:
        plt.ylabel("Metabolomics\nDown     |     Up   ")
    else:
        plt.yticks([], [])
    plt.xlim([-1.5, len(time_unique) - 0.5])
    spidx = spidx + 1

filename = "Figure_4_A.svg"
plt.savefig(filename, format="svg", bbox_inches="tight")
