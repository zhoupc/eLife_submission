#!/bin/bash 
pdflatex Fig_introduction.tex
pdflatex Fig_BG.tex 
pdflatex Fig_initialization.tex
pdflatex example_neurons.tex
pdflatex Fig_SIM.tex

# clean extra files 
rm *.aux *.log
