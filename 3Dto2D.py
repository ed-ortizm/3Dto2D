#!/usr/bin/python
import numpy as np

# arguments for the script
import sys
# sys.arg[0] is the name of the script itself (careful)
file_name = sys.argv[1]
filter_name= sys.argv[2]

n_arguments = len(sys.argv)
# raising an assetion error when the script is not provided with the right number of n_arguments
#https://stackoverflow.com/questions/1964126/what-exception-to-raise-if-wrong-number-of-arguments-passed-in-to-args
assert n_arguments == 3, "Please provide only two arguments to the script: arg1 = file name, arg2= filter name "
