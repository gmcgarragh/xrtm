include ../make.inc

PRODUCTS = xrtm.aux xrtm.bbl xrtm.blg xrtm.brf xrtm.log xrtm.out xrtm.toc # xrtm.pdf *.html

all: xrtm.pdf index.html

xrtm.pdf: xrtm.tex common.tex
	pdflatex -halt-on-error xrtm.tex
	bibtex xrtm
	pdflatex -halt-on-error xrtm.tex
	pdflatex -halt-on-error xrtm.tex

index.html: common.tex xrtm.tex xrtm.bib xrtm.pdf
	$(HOME)/src/tth_C-4.03/tth -Lxrtm -s < xrtm.tex | sed 's/\&times;</\&times; </g' | $(HOME)/src/tthgold/tthsplit
	sed 's|href="chap|href="sec|g' index.html > temp.out
	sed 's|Callxrtm input format</a><br />|Callxrtm input format</a><br /><a href="refs.html">Bibliography</a><br />|g' temp.out > index.html
	rm -f temp.out

include dep.inc
