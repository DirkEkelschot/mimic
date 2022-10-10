import numpy as np
import Cpp2Py
import operator
import sys
import matplotlib.pyplot as plt
from mpi4py import MPI
import time

def main():
    comm      = MPI.COMM_WORLD
    comm_info = MPI.Info.Create()
    rank      = comm.Get_rank()
    
    #======================================================================
    #======================================================================
    #======================================================================
    I = np.ones((10,1))
    B = np.ones((4,1))
    value    = Cpp2Py.CalculateEnergy(10,I,4,B);
    print(value)
    #======================================================================
    #======================================================================
    #======================================================================

if __name__ == "__main__": main()
