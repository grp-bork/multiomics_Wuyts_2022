.PHONY: all clean

all: Supp_Figure_2_AB.pdf Supp_Figure_2_CD.pdf

Supp_Figure_2_CD.pdf:
	Rscript Supp_Figure_2_CD.R

Supp_Figure_2_AB.pdf:
	Rscript Supp_Figure_2_AB.R
	
clean:
	rm -f *.svg *.pdf
