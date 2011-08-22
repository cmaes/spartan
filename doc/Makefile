FILES = $(wildcard *.tex)
PDF = $(FILES:.tex=.pdf)
BASE = $(firstword $(notdir $(basename $(FILES))))

all: $(PDF)

%.pdf: %.tex references.bib
	pdflatex --file-line-error-style $< </dev/null
	#uncomment for bibtex
	bibtex $(basename $<)
	pdflatex --file-line-error-style $< </dev/null
	pdflatex --file-line-error-style $< </dev/null

run:
	kpdf $(PDF)
clean: 
	rm -f *.log *.toc *.aux 
