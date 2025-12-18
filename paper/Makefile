all: pdf code clean
.PHONY: all tex pdf pdf2 code clean cleanall

# names of .Rnw .bib. and .bst files
FILE=BFGSD
BIB=bibliography

## generate .tex file from .Rnw file with knitr
tex: $(FILE).Rnw $(BIB).bib
	Rscript -e "knitr::knit('$(FILE).Rnw')" --vanilla

## generate .pdf file from .tex file
pdf: tex
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## generate .pdf file (without first knitting tex)
pdf2:
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

## extract R code
code: $(FILE).Rnw
	Rscript -e "knitr::purl('$(FILE).Rnw')" --vanilla

## clean
clean:  
	-rm $(FILE).aux  $(FILE).blg  $(FILE).log  $(FILE).tex  $(FILE).bbl \
	$(FILE).out $(FILE).brf $(FILE).tex $(FILE).synctex.gz

cleanall: clean
	-rm -r ./figure/
	-rm -r ./cache/
