#!/usr/bin/env python3

import numpy as np
from numpy.random import default_rng

# import sympy as sp
# from sympy import pprint


LINE_WIDTH = 200
DATA_TYPE = np.float64
PRECISION = 8
VERBOSE = False 

np.set_printoptions(linewidth=LINE_WIDTH, suppress=True, precision=PRECISION)

def verbose_print(text = ''):
    if VERBOSE:
        print(text)

def gaussian_elimination_partial_pivoting(A):
    m = len(A)      # Rows
    n = len(A[0])   # Columns
    A = np.asarray(A, dtype=DATA_TYPE)

    FLOPs_Counter = 0

    print("ORIGINAL MATRIX:")
    print(f"{m} x {n} matrix | DATA_TYPE = {DATA_TYPE}")
    print(A)
    print()

    # Loop through the cols of A
    for col in range(n):

        verbose_print(A)
        verbose_print()

        row_swapped = False 
        # Loop through the rows of the current column starting at the diagonal
        for i in range(col, m):
            # Check if the abs value of the current row is > than the diagonal "default pivot"
            if abs(A[i, col]) > abs(A[col, col]):
                # Swap current row with pivot row
                # Note: col = the pivot row as the "default pivot" is the diagonal 
                A[[i, col]] = A[[col, i]]
                row_swapped = True
                verbose_print(f"Swap rows {i} {col}")

        if row_swapped:
            verbose_print(A)
            verbose_print()

        # Loop through the rows after the pivot column row
        for i in range(col + 1, m):
            factor = np.float64(A[i, col]) / np.float64(A[col, col])
            FLOPs_Counter += 1

            verbose_print(f"Factor = {np.float16(A[i, col])} / {np.float16(A[col, col])}")
            verbose_print(f"E{i} -> E{i} - ({A[i, col]}/{A[col, col]})*E{col}")
            # Loop through the column values starting at the current column
            for j in range(col, n):
                # verbose_print(f"{A[i, j]} -= {factor * A[col, j]}")
                A[i, j] -= factor * A[col, j]
                FLOPs_Counter += 2

            verbose_print(A)
            verbose_print()

        verbose_print("=" * LINE_WIDTH)

    print("FINAL RESULT:")
    print(f"{m} x {n} matrix | DATA_TYPE = {DATA_TYPE}")
    print(f"FLOPs = {FLOPs_Counter}")
    print(A)


I = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
]

A = [
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 9],
]

A_2 = [
    [0.0001, 1, 1],
    [1, 1, 2],
]

A_3 = default_rng(340).random((11, 10))

gaussian_elimination_partial_pivoting(A_3)
print("-" * LINE_WIDTH)

gaussian_elimination_partial_pivoting(I)
print("-" * LINE_WIDTH)

gaussian_elimination_partial_pivoting(A)
print("-" * LINE_WIDTH)

gaussian_elimination_partial_pivoting(A_2)
print("-" * LINE_WIDTH)


# Compare to SymPy's Row Echelon Form
# print("COMPARE TO SYMPY'S ROW ECHELON FORM:")

# Don't display pivots in the output
# pprint(sp.Matrix(np.asarray(A_3, dtype=np.float64)).echelon_form(with_pivots=False))
# pprint(sp.Matrix(np.asarray(A_3, dtype=DATA_TYPE)).echelon_form(with_pivots=True))