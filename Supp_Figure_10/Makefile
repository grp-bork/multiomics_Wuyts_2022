.PHONY: all clean

all: Supp_Figure_10_3.svg

Supp_Figure_10_3.svg: ../.venv/bin/activate
	source ../.venv/bin/activate && python Supp_Figure_10.py
	
clean:
	rm -f *.svg
