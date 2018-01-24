# hqp2mtx

`hqp2mtx` converts a linear system Ax=b dumped by [HQP](https://github.com/omuses/hqp/) into Matrix Market format.

HQP dumps the linear system in Compressed Sparse Row format (CSR) with some small extensions. To solve this linear
systems using Meagre-Crowd, we need to convert them into Matrix Market format.


# Authors

N. Marcus      nora_lynn.marcus@tu-dresden.de
M. Schroschk   martin.schroschk@tu-dresden.de


# Usage
See `python3.6 hqp2mtx.py -h` for usage.
