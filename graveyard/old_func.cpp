int ChkHexorient(double* P, int* Pid) {
  int     changed,k,i;
  double  volref,volhex;
  int* Pnew = new int[8];
  volref = 1;
  double* c1 = new double[3];
  double* c2 = new double[3];
  double* c3 = new double[3];
  double* c4 = new double[3];
  c1[0] = P[0*8+0];c1[1] = P[0*8+1];c1[2] = P[0*8+2];
  c2[0] = P[1*8+0];c1[1] = P[1*8+1];c1[2] = P[1*8+2];
  c3[0] = P[3*8+0];c1[1] = P[3*8+1];c1[2] = P[3*8+2];
  c4[0] = P[4*8+0];c1[1] = P[4*8+1];c1[2] = P[4*8+2];
  /** check the orientability of the hexahedra : vol of tet p0 p1 p3 p4 */
  volhex = ComputeQuickVol(c1,c2,c3,c4);
  changed = 0;
  if ( volref*volhex < 0 )
  {
      changed = 1;
      Pnew[0] = Pid[0];
      Pnew[1] = Pid[3];
      Pnew[2] = Pid[2];
      Pnew[3] = Pid[1];
      Pnew[4] = Pid[4];
      Pnew[5] = Pid[7];
      Pnew[6] = Pid[6];
      Pnew[7] = Pid[5];
      
      Pid[0] = Pnew[0];
      Pid[1] = Pnew[1];
      Pid[2] = Pnew[2];
      Pid[3] = Pnew[3];
      Pid[4] = Pnew[4];
      Pid[5] = Pnew[5];
      Pid[6] = Pnew[6];
      Pid[7] = Pnew[7];
  }

  return changed;
}



//ien = 6;
//ief = 6;
//ifn = 6*4;



double ComputeQuickVol(double *c1,double *c2,double *c3,double *c4) {
  double   ax,ay,az,bx,by,bz,vol;

  ax = c3[0] - c1[0];
  ay = c3[1] - c1[1];
  az = c3[2] - c1[2];

  bx = c4[0] - c1[0];
  by = c4[1] - c1[1];
  bz = c4[2] - c1[2];

  vol =   (c2[0]-c1[0]) * (ay*bz - az*by) \
        + (c2[1]-c1[1]) * (az*bx - ax*bz) \
        + (c2[2]-c1[2]) * (ax*by - ay*bx);

  return vol;
}

