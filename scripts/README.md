October 2017
CONVERTER reads a matrix in CSR format and converts it to MTX format.

HQP outputs a matrix in CSR format in a single file, e.g., FullUMF_A.txt. It needs two steps to transform this matrix
file into Matrix Market format:
    1. Split the CSR single file into three files containing v, iv, and jv, respectively. (convert.sh)
    2. Convert the CSR matrix into Matrix Market format using Python Scipy library. (convert.py)


# Authors

N. Marcus      nora_lynn.marcus@tu-dresden.de
M. Schroschk   martin.schroschk@tu-dresden.de

# Content

README.md
converter.sh            reads the matrix in CSR format, invokes converter.py
converter.py            assembles the CSR matrix, converts and writes the MTX file
converter_multiple.sh   example on how to use converter for multiple files

# Installation
No compilation necessary. You need a python and bash interpreter.

# Usage
## Input from HQP
bash converter.sh <path-to-input-file>
Be aware that sh will usually not work.

## Usual CSR input files
Usually, a CSR format should consist of three files. If you name them v, jv and pv and put them into ./tmp/, you can 
directly invoke python convert.py to convert the matrix.

## Output
Output is always A.mtx and b.vec.

## Faulty matrices
Faulty matrices will be recognized during the execution of convert.py.


# Known Issues

## Base
The MTX format is 1-based. Some programs want a 0-based MTX file. For base 0, uncomment the line during output in convert.py.

## Converting multiple matrices
The output files will be overwritten without prompt. If you want to convert multiple files at once, write your own bash
script. You can find an example in convert_multiple.sh.

## ValueError: could not convert string to float
If you encounter this error message from Python

    Traceback (most recent call last):
      File "convert.py", line 32, in <module>
      v.append(float(line))
      ValueError: could not convert string to float: Vector:

than, check the input file for lines containing no '#' in lines containing no coordinate value or matrix entry.
