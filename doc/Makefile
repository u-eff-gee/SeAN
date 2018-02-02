BIBTEX=bibtex
LATEX=pdflatex

OUTPUT = theory.aux theory.bbl theory.blg theory.log theory.out theory.pdf theory.toc

theory.pdf: theory.tex
	$(LATEX) $<
	$(BIBTEX) theory
	$(LATEX) $<
	$(LATEX) $<

clean:
	rm -f $(OUTPUT)
