include ../make.inc
include ../common.inc

SUBDIRS =

OBJECTS =

BINARIES =

PRODUCTS = xrtm.aux xrtm.bbl xrtm.blg xrtm.brf xrtm.log xrtm.out xrtm.toc \
           xrtm_equations.aux xrtm_equations.log xrtm_equations.out xrtm_equations.toc \
           # *.html xrtm.pdf xrtm_equations.pdf

all: index.html xrtm.pdf xrtm_equations.pdf

xrtm.pdf: common.tex xrtm.tex xrtm.bib
	pdflatex -halt-on-error xrtm.tex
	bibtex xrtm
	pdflatex -halt-on-error xrtm.tex
	pdflatex -halt-on-error xrtm.tex

index.html: common.tex xrtm.tex xrtm.bib xrtm.pdf
	$(HOME)/src/tth_C-4.03/tth -Lxrtm -s < xrtm.tex | sed 's/\&times;</\&times; </g' | $(HOME)/src/tthgold/tthsplit
	sed 's|href="chap|href="sec|g' index.html > temp.out
	sed 's|Callxrtm input format</a><br />|Callxrtm input format</a><br /><a href="refs.html">Bibliography</a><br />|g' temp.out > index.html
	rm -f temp.out

xrtm_equations.pdf: xrtm_equations.tex equations/*.tex
#	pdflatex -halt-on-error xrtm_equations.tex
#	bibtex xrtm_equations
	pdflatex -halt-on-error xrtm_equations.tex
	pdflatex -halt-on-error xrtm_equations.tex

spell_check:
	aspell --disable-backup -c xrtm_equations.tex

xrtm_spell_check:
	aspell --disable-backup \
               --add-tex-command="bibliographystyle p" \
               --add-tex-command="citep p" \
               --add-tex-command="createarg p" \
               --add-tex-command="sourcename p" \
               --add-tex-command="sourcenamearg p" \
               --add-tex-command="textcolor p" \
               -c xrtm.tex

include depend.inc
