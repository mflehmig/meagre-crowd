"""
Synopsis:
    This script transforms a linear system Ax=b from HQP into Matrix Market format.
    
Usage:
    python hqp2mtx.py -b right-hand-side.txt -A matrix.txt -x solution.txt
    
        --> Output: right-hand-side.mtx, matrix.mtx and solution.mtx


Command Line Parameters:
    Mandatory:
        -A FILE        File (incl. path) containing the matrix A from HQP.
        -b FILE        File (incl. path) containing the right hand side from HQP.
        -x FILE        File (incl. path) containing the solution from HQP.

Author:
    martin.schroschk@tu-dresden
    nora_lynn.marcus@tu-dresden.de
"""

import sys
import os
import re
import argparse
from scipy import sparse


def read_A(file):
    """
    Read matrix A from file and create scipy.sparse.csr_matrix A.
    
    Note: The order of v, iv and jv in the file A.txt does not matter. 

    :param file: The file (incl. path) to the matrix file.
    :return: Matrix A in scipy.sparse.csr_matrix format.
    """
    # The three arrays representing a matrix in CSR format.
    iv=[]
    jv=[]
    v=[]
    
    file_object=open(file, "r")
    lines = file_object.read().split('#')
    
    for line in lines:
        #match = re.search('.*[a-zA-Z]*Vector.*\n', line)
        #if match is not None:
            #i = 1 +3
        id = re.search('(\S*v):\n', line)
        if id:
            csr_id = id.group(1)
        tmp = re.split('.*[a-zA-Z]*Vector.*\n', line)
        if len(tmp) > 1:
            if csr_id == "iv":
                floats = [float(x) for x in tmp[1].split()]
                iv = iv + floats
            elif csr_id == "jv":
                floats = [float(x) for x in tmp[1].split()]
                jv = jv + floats
            elif csr_id == "v":
                floats = [float(x) for x in tmp[1].split()]
                v = v + floats
            
    return sparse.csr_matrix((v, jv, iv))


def write_A(file, A):
    """
    Write the matrix A in Matrix Market format to file. If the original file is
    A.txt than we create the file A.mtx.

    :param file: The file (incl. path) to the original matrix A file.
    :param A: Matrix A in coordinate format.
    :return: Name of MTX file.
    """

    mtx_file = str.replace(file, '.txt', '.mtx')
    file_object=open(mtx_file, "w")

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


def read_vec(file):
    """
    Read a vector (,e.g., b and x) from file and write the values into the list vec.

    :param file: The file (incl. path) to the vector file.
    :return: List vec containg the values of vector.
    """ 
    file_object=open(file, "r")
    vec=[]
    for line in file_object:
        if not line.startswith('#'):
            floats = [float(x) for x in line.split()]
            vec = vec + floats
    file_object.close()

    return vec


def write_vec(file, vec):
    """
    Write a vector in Matrix Market format to file. If the original file is vec.txt
    than we create the file vec.mtx.

    :param file: The file (incl. path) to the original vector file.
    :param vec: Vector in coordinate format.
    :return: Name of MTX file.
    """

    mtx_file = str.replace(file, '.txt', '.mtx')
    file_object=open(mtx_file, "w")

    # write MTX file
    ## write header
    file_object.write("%%MatrixMarket matrix coordinate real general\n")
    file_object.write("%-------------------------------------------------------------------------------\n")
    ## add matrix information here (name, usage, creator, etc)
    ## file_object.write("")
    file_object.write("%-------------------------------------------------------------------------------\n")
    
    num_elems = len(vec)
    file_object.write("{} {} {}\n".format(num_elems, '1', num_elems))
    cnt = 1
    for item in vec:
        #file_object.write("%s 1 %s\n" % (cnt item))
        file_object.write("{} {} {}\n".format(cnt, '1', item))
        cnt += 1


def create_argument_parser():
    """
    """
    parser = argparse.ArgumentParser(prog='hqp2mtx', description='Convert right-hand side b into Matrix Market format.')
    parser.add_argument('-A', metavar='FILE', type=str, required=False,
                        help='Matrix A from HQP.')
    parser.add_argument('-b', metavar='FILE', type=str, required=False,
                        help='Right-hand side b from HQP.')
    parser.add_argument('-x', metavar='FILE', type=str, required=False,
                        help='Solution of linear system.')

    return parser


def main():
    # Parse command line arguments
    parser = create_argument_parser()
    args = parser.parse_args()  # will raise an error if the arguments are invalid and terminate the
#     if not (args.A or args.b or args.x):
#         parser.error('No action requested.')
    
    try:
        if args.A:
            file_A = os.path.abspath(args.A)
            print("  Convert matrix A from HQP-CSR format into Matrix Market format ...")
            A = read_A(file_A);
            # convert CSR to COO
            A = A.tocoo()
            write_A(file_A, A);
            print("  Done.")

        if args.b:
            file_b = os.path.abspath(args.b)
            print("\n  Convert right-hand side b from HQP format into Matrix Market format ...")
            b = read_vec(file_b);
            write_vec(file_b, b);
            print("  Done.")

        if args.x:
            file_x = os.path.abspath(args.x)
            print("\n  Convert solution x from HQP format into Matrix Market format ...")
            x = read_vec(file_x);
            write_vec(file_x, x);
            print("  Done.")  

        print("Finished.\n")
    except:
        e = sys.exc_info()[0]
        print("Error: %s" % e )


if __name__ == '__main__':
    main()