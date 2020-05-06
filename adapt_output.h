#include "adapt_datastruct.h"
#include "adapt_partition.h"
#include "adapt_io.h"
#include "adapt.h"
#include "adapt_array.h"

#ifndef ADAPT_OUTPUT_H
#define ADAPT_OUTPUT_H

using namespace std;

void OutputQuantityPartition(Partition* pa, Array<double>* Quan, MPI_Comm comm);

void OutputPartionVolumes(ParArray<int>* ien, Array<double>* xcn_on_root, MPI_Comm comm);

void OutputPartitionFaces();

void WriteBoundaryDataInSerial3(Array<double>* xcn);

void WriteBoundaryDataInSerial(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* detJ_verts, double* vol_verts, double* Jnorm_verts);

void WriteBoundaryDataInSerial2(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* vol_verts);


#endif
