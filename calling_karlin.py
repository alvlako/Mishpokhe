import os
import sys
import subprocess

import ctypes as ct
from ctypes import *
import numpy as np

subprocess.call(['cc', '-fPIC', '-shared', '-o', 'karlin_c.so', 'karlin_c.c'])

print('karlin')
so_file = "karlin_c.so"
my_functions = CDLL(so_file)

print('lambda')
arr = np.asarray([0.7, 0.0, 0.1, 0.0, 0.0, 0.2])

n_mem = len(arr)
arr_c = (c_double * n_mem)(*arr)

import sys
import io
from contextlib import redirect_stdout


ct.cast(arr_c, ct.POINTER(c_double))
for i in arr_c:
  print(i)
print('lambda bis')
pointer_to_middle = arr[2:].ctypes.data_as(ct.POINTER(ct.c_double))
print(pointer_to_middle.contents)
my_functions.BlastKarlinLambdaBis.argtypes = [POINTER(c_double),c_int32, c_int32]
my_functions.BlastKarlinLambdaBis.restype = c_double
print(my_functions.BlastKarlinLambdaBis(pointer_to_middle, -2, 3))


my_functions.BlastKarlinLambdaNR.argtypes = [POINTER(c_double),c_int32, c_int32]
my_functions.BlastKarlinLambdaNR.restype = c_double
BlastKarlin_lambda = my_functions.BlastKarlinLambdaNR(pointer_to_middle, -2, 3)
print("lambda", BlastKarlin_lambda)
BlastKarlin_lambda = ct.c_double(BlastKarlin_lambda)
print(type(BlastKarlin_lambda))

my_functions.BlastKarlinLtoH.restype = c_double
BlastKarlin_H = my_functions.BlastKarlinLtoH(pointer_to_middle, -2, 3, BlastKarlin_lambda)
print('H', BlastKarlin_H)
BlastKarlin_H = ct.c_double(BlastKarlin_H)

my_functions.BlastKarlinLHtoK.restype = c_double
BlastKarlin_K = my_functions.BlastKarlinLHtoK(pointer_to_middle, -2, 3, BlastKarlin_lambda, BlastKarlin_H)
print('K', BlastKarlin_K)