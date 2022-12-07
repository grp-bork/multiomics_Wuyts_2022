.PHONY: all download

export R_LIBS := $(shell pwd)/.R_packages

all: .venv .data_downloaded 
	mkdir -p $(R_LIBS)
	Rscript init.R && \
	cd Figure_1 && make && \
	cd ../Figure_2 && make && \
	cd ../Figure_3 && make && \
	cd ../Figure_4 && make && \
	cd ../Supp_Figure_01 && make && \
	cd ../Supp_Figure_02 && make && \
	cd ../Supp_Figure_03 && make && \
	cd ../Supp_Figure_04 && make && \
	cd ../Supp_Figure_05 && make && \
	cd ../Supp_Figure_06 && make && \
	cd ../Supp_Figure_07 && make && \
	cd ../Supp_Figure_08 && make && \
	cd ../Supp_Figure_09 && make && \
	cd ../Supp_Figure_10 && make && \
	cd ../Supp_Figure_11 && make && \
	cd ../Supp_Figure_12 && make && \
	cd ../Supp_Figure_13 && make

.data_downloaded:
	cd data && make && cd .. && touch $@

clean:
	cd Figure_1 && make clean
	cd Figure_2 && make clean
	cd Figure_3 && make clean
	cd Figure_4 && make clean
	cd Supp_Figure_01 && make clean
	cd Supp_Figure_02 && make clean
	cd Supp_Figure_03 && make clean
	cd Supp_Figure_04 && make clean
	cd Supp_Figure_05 && make clean
	cd Supp_Figure_06 && make clean
	cd Supp_Figure_07 && make clean
	cd Supp_Figure_08 && make clean
	cd Supp_Figure_09 && make clean
	cd Supp_Figure_10 && make clean
	cd Supp_Figure_11 && make clean
	cd Supp_Figure_12 && make clean
	cd Supp_Figure_13 && make clean

.venv:
	python3 -m venv .venv
	source .venv/bin/activate && pip install ete3 six numpy openpyxl pandas scipy matplotlib seaborn

.DELETE_ON_ERROR:
