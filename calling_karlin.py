from ctypes import *
import ctypes as ct
import numpy as np





#so_file = "karlin_c.so"
so_file = "ex.so"
my_functions = CDLL(so_file)
 
print(type(my_functions))
print(my_functions.square(10))
print(my_functions.square(8))

print('karlin')
so_file = "karlin_c.so"
my_functions = CDLL(so_file)

print('lambda')
arr = [0.7, 0.0, 0.1, 0.0, 0.0, 0.2]

n_mem = len(arr)
arr_c = (c_double * n_mem)(*arr)

import sys
import io
from contextlib import redirect_stdout


ct.cast(arr_c, ct.POINTER(c_double))
for i in arr_c:
  print(i)
print('lambda bis')
my_functions.BlastKarlinLambdaBis.argtypes = [POINTER(c_double),c_int32, c_int32]
my_functions.BlastKarlinLambdaBis.restype = c_double
print(my_functions.BlastKarlinLambdaBis(arr_c, -2, 3))





#f = io.StringIO()
#with redirect_stdout(f):
    #so_file = "karlin_c.so"
    #my_functions = CDLL(so_file)
    #BlastKarlin_lambda = ct.c_double(BlastKarlin_lambda)
    #print(type(BlastKarlin_lambda))
    #my_functions.BlastKarlinLambdaNR.argtypes = [POINTER(c_double),c_int32, c_int32]
    #my_functions.BlastKarlinLambdaNR.restype = c_double
    
    
    #BlastKarlin_lambda = my_functions.BlastKarlinLambdaNR(ct.cast(arr_c, ct.POINTER(c_double)), -2, 3)
    #print("lambda", BlastKarlin_lambda)
    #print(type(BlastKarlin_lambda))
#s = f.getvalue()
#print('stdout')
#print('stdout', s)
#print("stdout lambda", BlastKarlin_lambda)

#INTP = ct.POINTER(ct.c_double)
#addr = ct.addressof(BlastKarlin_lambda)
#print ('address:', addr, type(addr))
#BlastKarlin_lambda_ptr = ct.cast(addr, INTP)
#print(type(BlastKarlin_lambda))

#my_functions.BlastKarlinLambdaNR.argtypes = [POINTER(c_double),c_int32, c_int32]
#my_functions.BlastKarlinLambdaNR.restype = c_double
#BlastKarlin_lambda = my_functions.BlastKarlinLambdaNR(ct.cast(arr_c, ct.POINTER(c_double)), -2, 3)
#print('lambda', BlastKarlin_lambda)
#BlastKarlin_lambda = c_double(BlastKarlin_lambda)

#my_functions.BlastKarlinLtoH.restype = c_double
#BlastKarlin_H = my_functions.BlastKarlinLtoH(arr_c, -2, 3, BlastKarlin_lambda)
#print('H', BlastKarlin_H)
#BlastKarlin_H = ct.c_double(BlastKarlin_H)

#my_functions.BlastKarlinLHtoK.restype = c_double
#BlastKarlin_K = my_functions.BlastKarlinLHtoK(arr_c, -2, 3, BlastKarlin_lambda, BlastKarlin_H)
#print('K', BlastKarlin_K)