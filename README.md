# v2 Analysis port from C++ to python

## Name
Python and ROOT (PyRoot) based v2 analysis.

## Description
This repository aims to compare the same analysis in C++, pure python and python with C++ shared libraries loaded enhanced with ROOT/PyROOT. The comparison is based on a toy version of PHENIX, RHIC data of 200 GeV/c Au+Au collisions. 

## Dependencies
* ROOT (version compatible with your python) root-config --python-version
* Python 3.x
* C++ (clang 14.0.0 was used)

## Usage
Since this repository stores python and C++ codes, the usage is slightly different.
To compile C++ and C++ shared library, simply say

make all

To run C++, simply say
./exe/v2\_analysis.exe \<input file name\> \<output file name\> \<max events=-1\>

To run Python, simply say
python v2\_analysis\_cppLibraries.py \<input file name\> \<output file name\> \<max events=-1\> 


## Project status
Ongoing

