# This script reads 3 files that make a sparse matrix in CSR format and converts it to CSC format.

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
#file_object.write(str(max(A.row)+1) + " " + str(max(A.col)+1) + " " + str(len(A.data)) + "\n")
### otherwise, use this:
file_object.write(str(max(A.row)) + " " + str(max(A.col)) + " " + str(len(A.data)) + "\n")

## write matrix entries
for i in range(1, len(A.data)):
  value=A.data[i]
  row=A.row[i]
  column=A.col[i]

  file_object.write(str(row) + " " + str(column) + " " + str(value) + "\n")

file_object.close()





