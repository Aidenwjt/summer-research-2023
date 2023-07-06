all: compile

compile:
	pdflatex notes.tex
	bibtex notes.aux
	pdflatex notes.tex
	pdflatex notes.tex
	xdg-open notes.pdf
