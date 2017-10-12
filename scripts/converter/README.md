October 2017
CONVERTER reads a matrix in CSR format and converts it to MTX format.

# Author

N. Marcus 	nora_lynn.marcus@tu-dresden.de

# Content

README.md	
converter.sh		reads the matrix in CSR format, invokes converter.py
converter.py		assembles the CSR matrix, converts and writes the MTX file
converter_multiple.sh	example on how to use converter for multiple files

# Installation
No compilation necessary. You need a python and bash interpreter.

# Usage
## Input from HQP
bash converter.sh <path-to-input-file>
Be aware that sh will usually not work.
## Usual CSR input files
Usually, a CSR format should consist of three files. If you name them v, jv and pv and put them into ./tmp/, you can invoke 
python convert.py 
to convert the matrix.
## Output
Output is always A.mtx and b.vec.
## Faulty matrices
Faulty matrices will be recognized during the execution of convert.py.

# Known Issues
## Base
The MTX format is 1-based. Some programs want a 0-based MTX file. For base 0, uncomment the line during output in convert.py.
## Converting multiple matrices
The output files will be overwritten without prompt. If you want to convert multiple files at once, write your own bash script. You can find an example in convert_multiple.sh.
