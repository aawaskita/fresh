from math import pi, pow
import matplotlib.pyplot as plt
import numpy as np

"""
Development notes:
(2jan2019)
- cekricek KUDIF execution: kenapa GES nilainya NOL?
- perlu cekricek perubahan index dari fortran ke python !!!
"""

#from triac2_utils import *
#from InputData import *
import triac2_utils
import InputData
#import KUDIFed

from triac2_utils import anfang, temp, difpo, difko, kudif
from InputData import readdata

# Dummy properties
inpf = 'inptriac2'
data = InputData.readdata(inpf)
#
print(data['geodi'])
## ANFANG soubroutine 
anfang(data)


