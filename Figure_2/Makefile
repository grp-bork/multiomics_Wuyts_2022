.PHONY: all clean

all: Figure_2_E2.svg Figure_2_D.svg Figure_2_C1.svg Figure_2_C2.svg Figure_2_B.pdf Figure_2_A.pdf

Figure_2_E2.svg:
	Rscript Figure_2_E.R

Figure_2_D.svg:
	Rscript Figure_2_D.R
	
Figure_2_C1.svg:
	Rscript Figure_2_C1.R

Figure_2_C2.svg: ../.venv/bin/activate
	source ../.venv/bin/activate && python Figure_2_C2.py
	
Figure_2_B.pdf:
	Rscript Figure_2_B.R
	
Figure_2_A.pdf:
	Rscript Figure_2_A.R
	
clean:
	rm -f *.svg *.pdf
