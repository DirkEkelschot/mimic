import numpy as np
import madam as madam
import operator
import sys
import matplotlib.pyplot as plt
from mpi4py import MPI
import time

def main():

  comm      = MPI.COMM_WORLD
  comm_info = MPI.Info.Create();
  rank      = comm.Get_rank()
  # mesh_io contains ien_p, ief_p, iee_p
  
  # This routine reads the the mesh and the data in parallel where the layout of the mesh_read is essentialy a struct of the required mapping arrays. mesh_read is a list with 9 entries where the entries themselfs are arrays that represent the following:
  #  Vertices:
  #  mesh_read[0] = xcn_read
  #  Element to entity maps:
  #  mesh_read[1] = iet_read, mesh_read[2] = ien_read, mesh_read[3] = ief_read, mesh_read[4] = iee_read;
  #  Face to entity maps:
  #  mesh_read[5],= ifn_read, mesh_read[6] = ife_read, mesh_read[7] = if_ref_read, mesh_read[8] = if_nv_read
  #mesh_read = madam.ReadDomainData("conn.h5","grid.h5","data.h5",0,comm,comm_info);

  mesh_part = madam.Partition("conn.h5","grid.h5","data.h5",0,"Mach",comm,comm_info);
  
  #Vrts_part     = mesh_part[0];#         -> List of local vertex coordinates at rank;
  #LocEl_part    = mesh_part[1];#         -> list of global element ID at rank;
  #Ustate_part   = mesh_part[2];#         -> {key=El_globID, value = U at cell center}
  #gV2lV_part    = mesh_part[3];#         -> {key=Vrt_globID, value = Vrt_locID}
  #LocEl2Nf_part = mesh_part[4];#         -> {key=El_globID, value = Number of Faces}
  #LocEl2Nv_part = mesh_part[5];#         -> {key=El_globID, value = Number of Vertices}
  #gE2lV_part    = mesh_part[6];#         -> {key=El_globID, value = list of Vrt_locIDs}
  #iee_part      = mesh_part[7];#         -> {key=El_globID, value = list of Adjacent El_globIDs}
  #ief_part      = mesh_part[8];#         -> {key=El_globID, value = list of Face_globIDs}
  #ifn_part      = mesh_part[9];#         -> {key=Face_globID, value = list of Vrt_globIDs}
  #if_Nv_part    = mesh_part[10];#        -> {key=Face_globID, value = number of Verts for that face}
  #nGlobElem     = mesh_part[11];#        -> integer idicating the global number of elements.
  #Ustate_ghost  = mesh_part[12];#        -> List of ghost values. This copy exists on each rank.

  #dUdXi is a dictionary/map where the key is the global element ID (gID) and the value is the corresponding [dUdx_gID,dUdy_gID,dUdz_gID]^t
  
  #dUdXi = madam.GradU(mesh_part,Ustate_part,comm,comm_info);#   -> {key=El_globID,   value = tuple of [dUdXi,dUdYi,dUdZi]}
  
#  nvert = len(Vrts_part);
#  xv = np.zeros((nvert,1));
#  yv = np.zeros((nvert,1));
#  zv = np.zeros((nvert,1));
#
#  for i in range(0,nvert):
#    xv[i] =  Vrts_part[i][0];
#    yv[i] =  Vrts_part[i][1];
#    zv[i] =  Vrts_part[i][2];
#
#  nLoc = len(LocEl_part);
#  for i in range(0,nLoc):
#
#
##  fig = plt.figure()
##  ax = fig.add_subplot(projection='3d')
##  ax.scatter(xv, yv, zv, '.')
#
#
#  plt.show()
if __name__ == "__main__": main()
