#include "adapt_datatype.h"

#ifndef ADAPT_GEOMETRY_H
#define ADAPT_GEOMETRY_H

Vec3D* ComputeSurfaceNormal(Vec3D* a, Vec3D* b);
double ComputeQuadSurfaceArea(double *P);
double ComputeTriSurfaceArea(double* P);
#endif
