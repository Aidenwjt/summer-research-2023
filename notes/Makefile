all: compile

compile:
	pdflatex notes.tex
	bibtex notes.aux
	pdflatex notes.tex
	pdflatex notes.tex
	rm references.bib
	xdg-open notes.pdf
