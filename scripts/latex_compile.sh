#!/bin/sh

pdflatex -interaction=batchmode -shell-escape  $1 
pdflatex -interaction=batchmode -shell-escape  $1

exit $?
