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
        --dir PATH     All ".txt" files in the provided directory (not recursive!).

Author:
    martin.schroschk@tu-dresden
    nora_lynn.marcus@tu-dresden.de
"""

import sys
import os
import re
import argparse
from scipy import sparse
#from lxml import _elementpath
#from pyanaconda.exception import list_addons_callback


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


def write_A(file, A, zerobased):
    """
    Write the matrix A in Matrix Market format to file. If the original file is
    A.txt than we create the file A.mtx.

    :param file: The file (incl. path) to the original matrix A file.
    :param A: Matrix A in coordinate format.
    :param zerobased: If true, create zero-based mtx file, otherwise one-based.
    :return: Name of MTX file.
    """

    ### if, for some reason, you need a zero-based mtx file:
    if zerobased:
      offset=0
    # Default: one-index mtx file
    else:
      offset=1

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
    file_object.write(str(max(A.row)+offset) + " " + str(max(A.col)+offset) + " " + str(len(A.data)) + "\n")

    ## write matrix entries
    ## Add value 1 to row and col indices in order to convert from 0-based to 1-based indexing!
    for i in range(0, len(A.data)):
        value=A.data[i]
        row=A.row[i]+offset
        column=A.col[i]+offset
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


def write_vec(file, vec, zerobased):
    """
    Write a vector in Matrix Market format to file. If the original file is vec.txt
    than we create the file vec.mtx.

    :param file: The file (incl. path) to the original vector file.
    :param vec: Vector in coordinate format.
    :param zerobased: If true, create zero-based mtx file, otherwise one-based.
    :return: Name of MTX file.
    """

    ### if, for some reason, you need a zero-based mtx file:
    if zerobased:
      offset=1
    # Default: one-index mtx file
    else:
      offset=0

    mtx_file = str.replace(file, '.txt', '.mtx')
    file_object=open(mtx_file, "w")

    # write MTX file
    ## write header
    file_object.write("%%MatrixMarket matrix coordinate real general\n")
    #file_object.write("%-------------------------------------------------------------------------------\n")
    ## add matrix information here (name, usage, creator, etc)
    ## file_object.write("")
    #file_object.write("%-------------------------------------------------------------------------------\n")

    num_elems = len(vec) - offset
    #file_object.write("{} {} {}\n".format(num_elems, '1', num_elems))
    file_object.write("{} {} {}\n".format(num_elems, 1-offset, num_elems))
    cnt = 1
    for item in vec:
        #file_object.write("%s 1 %s\n" % (cnt item))
        file_object.write("{} {} {}\n".format(cnt - offset, 1-offset, item))
        cnt += 1


def create_argument_parser():
    """
    Define command line arguments and provide help.
    """
    parser = argparse.ArgumentParser(prog='hqp2mtx',
                                     description='Convert linear system Ax=b from HQP-dump format (CSR) into Matrix Market format.')
    parser.add_argument('-A', metavar='FILE', type=str, required=False,
                        help='Matrix A from HQP.')
    parser.add_argument('-b', metavar='FILE', type=str, required=False,
                        help='Right-hand side b from HQP.')
    parser.add_argument('-x', metavar='FILE', type=str, required=False,
                        help='Solution of linear system.')
    parser.add_argument('--zerobased',  action='store_true',
                        help='Set to true, if you want zero-based mtx files.')
    parser.add_argument('--dir', metavar='PATH', type=str, required=False,
                        help='Convert all .txt files in provided directory.')

    return parser


def get_filenames(path):
    """
    Run over all files in the given directory. Files that end with "A.txt", "b.txt" or "x.txt" are saved to three
    separate lists. Return the lists containg the matrix, right-hand and solution files.
    """
    abspath = os.path.abspath(path)
    if not os.path.exists(abspath):
        exit

    list_A = []
    list_b = []
    list_x = []
    # Do not search in subdirectories, i.e., do not work recursive
    for root, dirs, files in os.walk(abspath):
        for file in files:
            if file.endswith("A.txt"):
                list_A.append(os.path.join(abspath, file))
            if file.endswith("b.txt"):
                list_b.append(os.path.join(abspath, file))
            if file.endswith("x.txt"):
                list_x.append(os.path.join(abspath, file))

    return (list_A, list_b, list_x)

def main():
    # Parse command line arguments
    parser = create_argument_parser()
    args = parser.parse_args()  # will raise an error if the arguments are invalid and terminate the
    if not (args.A or args.b or args.x or args.dir):
        parser.error('No action requested.')

    print("  zerobased: ", args.zerobased)

    if args.dir:
        (list_A, list_b, list_x) = get_filenames(args.dir)

        if list_A:
            print("Convert matrix files: ")
            for file in list_A:
                #print("  Convert matrix A from HQP-CSR format into Matrix Market format ...")
                print("  ", file)
                A = read_A(file);
                # convert CSR to COO
                A = A.tocoo()
                write_A(file, A, args.zerobased);
            print("Done.\n")

        if list_b:
            print("Convert vector b files: ")
            for file in list_b:
                #print("  Convert matrix A from HQP-CSR format into Matrix Market format ...")
                print("  ",  file)
                b = read_vec(file);
                write_vec(file, b, args.zerobased);
            print("Done.\n")

        if list_x:
            print("Convert vector x files: ")
            for file in list_x:
                #print("  Convert matrix A from HQP-CSR format into Matrix Market format ...")
                print("  ",  file)
                x = read_vec(file);
                write_vec(file, x, args.zerobased);
            print("Done.\n")

    try:
        if args.A:
            file_A = os.path.abspath(args.A)
            print("  Convert matrix A from HQP-CSR format into Matrix Market format ...")
            A = read_A(file_A);
            # convert CSR to COO
            A = A.tocoo()
            write_A(file_A, A, args.zerobased);
            print("  Done.")

        if args.b:
            file_b = os.path.abspath(args.b)
            print("\n  Convert right-hand side b from HQP format into Matrix Market format ...")
            b = read_vec(file_b);
            write_vec(file_b, b, args.zerobased);
            print("  Done.")

        if args.x:
            file_x = os.path.abspath(args.x)
            print("\n  Convert solution x from HQP format into Matrix Market format ...")
            x = read_vec(file_x);
            write_vec(file_x, x, args.zerobased);
            print("  Done.")  

        print("Finished.\n")
    except:
        e = sys.exc_info()[0]
        print("Error: %s" % e )


if __name__ == '__main__':
    main()
