
#This will overwrite, but assume everything in Rcode/BAMMtools is current
rm BAMMtools.R
cat ../../Rcode/BAMMtools/*.R >> BAMMtools.R

 
R CMD batch bammdoc_makefigs.R

R CMD sweave bammdoc1.Rnw

# Run this three times for good measure!

pdflatex bammdoc1.tex
pdflatex bammdoc1.tex
pdflatex bammdoc1.tex

