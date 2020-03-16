#!/usr/bin/python
import numpy as np

# arguments for the script
import sys
# sys.arg[0] is the name of the script itself (careful)
################################ code starts ###################################
# raising an assetion error when the script is not provided with the right number of n_arguments or the right format
#https://stackoverflow.com/questions/1964126/what-exception-to-raise-if-wrong-number-of-arguments-passed-in-to-args

n_arguments = len(sys.argv)
if n_arguments == 3:
    file_name = sys.argv[1]
    assert (file_name[-5:] == ".fits"), "The format of the file is .fits: file_name.fits (lower case)"
    filter_name= sys.argv[2]
#assert (filter_name[-5:] == ".fits"), "Please remember that the format for the filter name is fits: fiter_name.fits, (lower case)"
# or wavelenght range
elif n_arguments == 4:
    file_name = sys.argv[1]
    assert (file_name[-5:] == ".fits"), "The format of the file is .fits: file_name.fits (lower case)"
    lambda1 = float(sys.argv[2])
    lambda2 = float(sys.argv[3])
    assert lambda1 < lambda2, "first wavelenght must be smaller than the second one"
else:
    assert False, "# of args must be 2 (file and filter name) OR 3 (file name, lower wavelenght and higher wavelenght)"
