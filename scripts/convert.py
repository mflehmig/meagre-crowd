"""
Synopsis:
    This script reads 3 files that make a sparse matrix in CSR format and converts it to Matrix Market format which is
    basicaly coordinate format (COO).
    We expect
        ./tmp/v
        ./tmp/iv
        ./tmp/jv
    to be present and only containing real/float values (no comments!).
    Furthermore, we expect the input CSR matrix to be 0-based, i.e., indices starting with 0.
    The Matrix Market output format will be 1-based so that it can be processed by meagre-crowd.

    For further information see README.md.

Usage:
    python createModuleFiles.py --sw /PATH/TO/OMC/DEPENDENCIES/build/clang3.7/ \
                                --inst /PATH/TO/OMC/DEPENDENCIES/modules/ \
                                --compiler=clang --version=3.7 \
                                --templates /PATH/TO/OMC/DEPENDENCIES/templates/

Author:
    nora_lynn.marcus@tu-dresden.de
    martin.schroschk@tu-dresden.de

"""

from scipy import sparse


# read matrix in CSR format
file_object=open("./tmp/v", "r")
v=[]
for line in file_object:
  v.append(float(line))
file_object.close()

file_object=open("./tmp/jv", "r")
jv=[]
for line in file_object:
  jv.append(int(line))
file_object.close()

file_object=open("./tmp/iv", "r")
iv=[]
for line in file_object:
  iv.append(int(line))
file_object.close()

# assemble matrix in scipy CSR format
A = sparse.csr_matrix((v,jv,iv))

# convert CSR to COO
A = A.tocoo()

file_object=open("./A.mtx", "w")

# write MTX file
## write header
file_object.write("%%MatrixMarket matrix coordinate real general\n")
file_object.write("%-------------------------------------------------------------------------------\n")
## add matrix information here (name, usage, creator, etc)
## file_object.write("")
file_object.write("%-------------------------------------------------------------------------------\n")

## maximum row, column and number of nonzeros
### if, for some strange reason, you need a zero-based mtx file:
# file_object.write(str(max(A.row)) + " " + str(max(A.col)) + " " + str(len(A.data)) + "\n")
### otherwise, use this:
file_object.write(str(max(A.row)+1) + " " + str(max(A.col)+1) + " " + str(len(A.data)) + "\n")

## write matrix entries
## Add value 1 to row and col indices in order to convert from 0-based to 1-based indexing!
for i in range(0, len(A.data)):
  value=A.data[i]
  row=A.row[i]+1
  column=A.col[i]+1

  file_object.write(str(row) + " " + str(column) + " " + str(value) + "\n")

file_object.close()





