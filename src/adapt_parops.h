#include "adapt.h"
#include "adapt_partition.h"
#include "adapt_datatype.h"

#ifndef ADAPT_PAROPS_H
#define ADAPT_PAROPS_H

using namespace std;


Array<double>* GetOptimizedMMG3DMeshOnRoot(Partition* P, US3D* us3d, std::map<int,Array<double>*> mv_map, MPI_Comm comm);

#endif
