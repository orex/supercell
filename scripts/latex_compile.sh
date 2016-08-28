#!/bin/sh

pdflatex -interaction=nonstopmode -shell-escape  $1 
pdflatex -interaction=nonstopmode -shell-escape  $1

exit $?
