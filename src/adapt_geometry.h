#include "adapt_datastruct.h"

#ifndef ADAPT_GEOMETRY_H
#define ADAPT_GEOMETRY_H

double CheckFaceOrientation(std::vector<double> VcF,
                            std::vector<std::vector<double> > Vfaces,
                            std::vector<double> Vijk);

std::vector<double> ComputeSurfaceNormal(std::vector<double> a,
                                         std::vector<double> b);

double ComputeQuadSurfaceArea(double *P);
double ComputeTriSurfaceArea(double* P);
#endif
