.PHONY: all clean

all: Figure_3_A.pdf Figure_3_B2.svg

Figure_3_A.pdf:
	Rscript Figure_3_A.R

Figure_3_B2.svg: ../.venv/bin/activate
	source ../.venv/bin/activate && python Figure_3_B.py

clean:
	rm -f Figure_3_A.pdf Figure_3_B1.svg Figure_3_B2.svg
