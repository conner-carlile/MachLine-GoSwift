#!/bin/bash
# Saves Contour plot as .svg
# Opens .svg in inkscape
# Resizes inkscape document
# Saves inkscape document as .emf for windows applications

#Contour Plot Script
SCRIPT='contour_plot.py'
TEMPFILE='tempfile.svg'

echo 'Please Enter Input Filename (without Extension): '
read INPUTFILE

python $SCRIPT $INPUTFILE

if [ -e $INPUTFILE$TEMPFILE ]
then
    #Ask for Final Filename
    echo ' '
    echo 'Please Enter .emf Filename (without Extension): '
    #Read user input to the EMFFILE variable
    read EMFFILE
    #open inkscape
    inkscape $INPUTFILE$TEMPFILE --export-emf=$EMFFILE'.emf'
    #delete temporary file
    rm -r $INPUTFILE$TEMPFILE
    echo ' '
    echo 'File '$EMFFILE'.emf was saved successfully!'
else
    echo ' '
    echo 'File was not saved!'
    exit 1
fi




