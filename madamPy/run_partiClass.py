import numpy as np
import PartiClass as parti
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

    # Initializing the paritition/Mesh object
    pyPart    = parti.Partition("conn.h5","grid.h5","data.h5",0,"Mach",comm,comm_info)
    vertices  = pyPart.getVertices()
    nv        = len(vertices)
    Ustate    = pyPart.getUState()             #   Ustate ->  {key=El_globID,   value = Ustate}
    Ustate    = pyPart.addAdjUState(Ustate)    #   Ustate ->  {key=El_globID,   value = Ustate}
    gUstate   = pyPart.computeGradU(Ustate)    #   gUstate -> {key=El_globID,   value = tuple of [dUdXi,dUdYi,dUdZi]}
    gather    = pyPart.gatherMeshOnRoot(Ustate)
#    for i in range(0,nv):
#        print(vertices[i][0],vertices[i][1],vertices[i][2])
 
    
    #======================================================================
    #======================================================================
    #======================================================================




if __name__ == "__main__": main()
