from ete3 import Tree
import openpyxl


wb_obj = openpyxl.load_workbook("../data/phylogeny/Lineage_of_32_sp.xlsx")
sheet = wb_obj.active

for i, cols in enumerate(sheet.columns):
    if i == 0:
        ntid = tuple(map(lambda x: x.value, cols[2:]))
    elif i == 7:
        species = tuple(map(lambda x: x.value, cols[2:]))
    elif i == 8:
        labels = tuple(map(lambda x: x.value, cols[2:]))

prefixed_names = map(lambda x: "-".join(x), zip(ntid, species))
names = dict(zip(labels, prefixed_names))

t = Tree("../data/phylogeny/bac120_r95.tree", format=1, quoted_node_names=True)
print("Pruning...")
t.prune(labels)
print("Done pruning")

# Replace IDs with species names
for n in t.get_leaves():
    n.name = names[n.name]

t.write(format=9, quoted_node_names=True, outfile="Figure-1-B.tree")
