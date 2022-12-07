import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fileFolder = "../data/metabolomics/"

# metaB neg
metaBnegfile = "metabolomics_negative_ions_abundance.tsv.gz"
metaBnegdata = pd.read_csv(fileFolder + metaBnegfile, sep="\t", index_col=0)

# metaB pos
metaBposfile = "metabolomics_positive_ions_abundance.tsv.gz"
metaBposdata = pd.read_csv(fileFolder + metaBposfile, sep="\t", index_col=0)

# metaB pos annotation
metaBposfile = "metabolomics_positive_annotations.tsv.gz"
metaBposann = pd.read_csv(fileFolder + metaBposfile, sep="\t", index_col=0)

# metaB neg annotation
metaBnegfile = "metabolomics_negative_annotations.tsv.gz"
metaBnegann = pd.read_csv(fileFolder + metaBnegfile, sep="\t", index_col=0)

# read no-bacteria control samples
metaBnegfile = "metabolomics_negative_ions_abundance.without_bacteria.tsv.gz"
metaBnegdata_control = pd.read_csv(fileFolder + metaBnegfile, sep="\t", index_col=0)

# metaB pos
metaBposfile = "metabolomics_positive_ions_abundance.without_bacteria.tsv.gz"
metaBposdata_control = pd.read_csv(fileFolder + metaBposfile, sep="\t", index_col=0)

print("Read all data")

metaBposdata_control.columns = [
    item.replace("_r", "_nobact_r") for item in metaBposdata_control.columns
]
metaBposdata = pd.concat([metaBposdata, metaBposdata_control], axis=1)

metaBnegdata_control.columns = [
    item.replace("_r", "_nobact_r") for item in metaBnegdata_control.columns
]
metaBnegdata = pd.concat([metaBnegdata, metaBnegdata_control], axis=1)

print("Concat main data")

def quantileNormalize(df_input):
    df = df_input.copy()
    # compute rank
    dic = {}
    for col in df:
        dic[col] = df[col].sort_values(na_position="first").values
    sorted_df = pd.DataFrame(dic)
    # rank = sorted_df.mean(axis = 1).tolist()
    rank = sorted_df.median(axis=1).tolist()
    # sort
    for col in df:
        # compute percentile rank [0,1] for each score in column
        t = df[col].rank(pct=True, method="max").values
        # replace percentile values in column with quantile normalized score
        # retrieve q_norm score using calling rank with percentile value
        df[col] = [
            np.nanpercentile(rank, i * 100) if ~np.isnan(i) else np.nan for i in t
        ]
    return df


metaBnegdata[metaBnegdata == 0] = np.nan
metaBposdata[metaBposdata == 0] = np.nan

print("Normalizing quantiles")

metaBnegdata_quantnorm = quantileNormalize(metaBnegdata)
metaBposdata_quantnorm = quantileNormalize(metaBposdata)

print("Normalized quantiles")

metaBnegdata_quantnorm = metaBnegdata_quantnorm.fillna(0)
metaBposdata_quantnorm = metaBposdata_quantnorm.fillna(0)

metaBdata_quantnorm = metaBposdata_quantnorm.copy()
metaBdata_quantnorm["mode"] = "pos"

test = metaBnegdata_quantnorm.copy()
test["mode"] = "neg"
metaBdata_quantnorm = pd.concat([metaBdata_quantnorm, test])
del test

print("Finalized plotting data... plotting")

def make_plot(index):
    plotxA = []
    plotxB = []
    plottimeA = []
    plottimeB = []
    plotxNobact = []
    plottimeNobact = []
    time_unique = ["T0", "T15", "T30", "T1h", "T3h", "T48h", "T96h"]
    for time_i in time_unique:
        curcol = [
            item
            for item in curdata.columns
            if ("A" in item) and (time_i in item) and not ("nobact" in item)
        ]
        plotxA.extend(curdata[curcol].values)
        plottimeA.extend([time_unique.index(time_i)] * len(curcol))
        curcol = [
            item
            for item in curdata.columns
            if ("B" in item) and (time_i in item) and not ("nobact" in item)
        ]
        plotxB.extend(curdata[curcol].values)
        plottimeB.extend([time_unique.index(time_i)] * len(curcol))

        curcol = [
            item for item in curdata.columns if (time_i in item) and ("nobact" in item)
        ]
        plotxNobact.extend(curdata[curcol].values)
        plottimeNobact.extend([time_unique.index(time_i)] * len(curcol))

    plt.clf()
    plt.figure(figsize=(5, 5))
    plt.rcParams.update({"font.size": 8})
    plt.scatter(plottimeA, plotxA)
    plt.scatter(plottimeB, plotxB)
    plt.scatter(plottimeNobact, plotxNobact, color="gray")
    plt.legend(["Run A", "Run B", "No bacteria"])
    plt.ylim([0, 10])
    plt.xticks(list(set(plottimeA)), time_unique)
    plt.title(plottitle)

    # # save figure to file
    plt.rcParams["svg.fonttype"] = "none"
    filename = f"Supp_Figure_10_{index}.svg"
    plt.savefig(filename, format="svg", bbox_inches="tight")


drugmzThreshold = 0.001

drugmz = 319.10093
annotated_indeces = [
    i for i in metaBdata_quantnorm.index if abs(i - drugmz) <= drugmzThreshold
]

samplecolumnsB = [
    item
    for item in metaBdata_quantnorm.columns
    if (("_r1" in item) or ("r2" in item))
    and ("Chlor" in item) and not ("r3" in item)
]
curdata = metaBdata_quantnorm.loc[annotated_indeces, samplecolumnsB].copy()
plottitle = "Chlorpromazine"

print("Plotting Chlorpromazine")
make_plot(1)

drugmz = 130.10812  # ,  319.10093]324.9787 #
annotated_indeces = [
    i for i in metaBdata_quantnorm.index if abs(i - drugmz) <= drugmzThreshold
]

samplecolumnsB = [
    item
    for item in metaBdata_quantnorm.columns
    if (("_r1" in item) or ("r2" in item))
    and ("Metfo" in item) and not ("r3" in item)
]
curdata = metaBdata_quantnorm.loc[annotated_indeces, samplecolumnsB].copy()
plottitle = "Metformin"

print("Plotting Metformin")
make_plot(2)

# check for compounds that correlate with drug (massspec drug artifacts)
test = metaBnegdata.copy()

drugmz = 324.9787  #
annotated_indeces = [i for i in test.index if abs(i - drugmz) <= drugmzThreshold]


samplecolumnsB = [
    item
    for item in metaBdata_quantnorm.columns
    if (("_r1" in item) or ("r2" in item))
    and ("Niclo" in item) and not ("r3" in item)
]
curdata = test.loc[annotated_indeces, samplecolumnsB].copy()
plottitle = "Niclosamide"

print("Plotting Niclosamide")
make_plot(3)
