.PHONY: all clean

all: Figure_4_A.svg Figure_4_B.pdf Figure_4_D.svg Figure_4_E.pdf Figure_4_F.pdf

Figure_4_A.svg: ../.venv/bin/activate
	source ../.venv/bin/activate && python Figure_4_A.py

Figure_4_B.pdf:
	Rscript Figure_4_B.R

Figure_4_D.svg: ../.venv/bin/activate
	source ../.venv/bin/activate && python Figure_4_CD.py

Figure_4_E.pdf:
	Rscript Figure_4_E.R

Figure_4_F.pdf:
	Rscript Figure_4_F.R

clean:
	rm -f *.pdf *.svg
