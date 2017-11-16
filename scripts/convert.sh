#!/bin/bash

# This is a short script to read a matrix in CSR format A and a vector b and save it in temporary files.
# A python program is invoked to reorder the matrix to Matrix Market Format. The Matrix is saved in an
# MTX file and the temporary files deleted.

# For further information see README.md.

# Directory of this script (to call convert.py)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# VARIABLES #
## input file is the first argument
input_file=$1
echo "parsing $input_file"

## output files
### temporary
mkdir -p ./tmp
V=./tmp/v
echo -e -n "" > $V
JV=./tmp/jv
echo -e -n "" > $JV
IV=./tmp/iv
echo -e -n "" > $IV
B=./b.vec
echo -e -n "" > $B
### for the python script
echo -e -n "" > ./A.mtx
echo "saving temporary files in ./tmp"

## vectors to hold the matrix in CSR format
declare -a v
declare -a jv
declare -a iv

## vector to hold the b vector
declare -a b

# CODE for CSR matrix #
# ## V ##
# ### read from '# v:' to '# jv:', delete all lines starting with #
# v=$( sed -n '/# v:/,/# jv:/p' $input_file | sed '/^#/d' )
# ### echo deletes multiple blanks, tr replaces blanks with \n
# v=$( echo $v | tr ' ' '\n' )
# ### output to file
# for element in $v; do
#   echo $element >> $V
# done
# 
# ## JV ##
# jv=$( sed -n '/ jv:/,/# iv:/p' $input_file | sed '/^#/d' )
# jv=$( echo $jv | tr ' ' '\n' )
# for element in $jv; do
#   echo $element >> $JV
# done
# 
# ## IV ##
# iv=$( sed -n '/# iv:/,/# b/p' $input_file | sed '/^#/d' )
# iv=$( echo $iv | tr ' ' '\n' )
# for element in $iv; do
#   echo $element >> $IV
# done
# 
# # B #
# ## read from '# b' until EOF, delete lines starting with #
# b=$( sed -n '/# b/,$p' $input_file | sed '/^#/d' )
# b=$( echo $b | tr ' ' '\n' )
# for element in $b; do
#   echo $element >> $B
# done

## IV ##
### read from '# iv:' to '# jv:', delete all lines starting with #
iv=$( sed -n '/# iv:/,/# jv:/p' $input_file | sed '/^#/d' )
iv=$( echo $iv | tr ' ' '\n' )
for element in $iv; do
  echo $element >> $IV
done

## JV ##
jv=$( sed -n '/ jv:/,/# v:/p' $input_file | sed '/^#/d' )
jv=$( echo $jv | tr ' ' '\n' )
for element in $jv; do
  echo $element >> $JV
done

## V ##
## read from '# v:' until EOF, delete lines starting with #
v=$( sed -n '/# v:/,$p' $input_file | sed '/^#/d' )
### echo deletes multiple blanks, tr replaces blanks with \n
v=$( echo $v | tr ' ' '\n' )
### output to file
for element in $v; do
  echo $element >> $V
done



# INVOKE python #
echo "converting matrix"
echo "writing MTX file"
python $DIR/convert.py

# DELETE tmp files #
#rm -rf ./tmp

echo "Finished."
