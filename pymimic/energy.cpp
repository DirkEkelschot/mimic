#include "energy.h"
#include <math.h>
#include <iostream>
#include <fstream>

double ComputeEnergy(int nB, double * B, int nI, double *I){
  
  double E = 0.0;
  for(int i=0;i<nB;i++){
    E += B[i];
  }
  for(int i=0;i<nI;i++){
    E += I[i];
  }
  return E;
}
