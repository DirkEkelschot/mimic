#include "adapt_datastruct.h"

#ifndef ADAPT_GEOMETRY_H
#define ADAPT_GEOMETRY_H

double CheckFaceOrientation(Vert* VcF, std::vector<Vert*> Vfaces, Vert* Vijk);
Vec3D* ComputeSurfaceNormal(Vec3D* a, Vec3D* b);
double ComputeQuadSurfaceArea(double *P);
double ComputeTriSurfaceArea(double* P);
#endif
