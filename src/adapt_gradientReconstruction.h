#include "adapt_partobject.h"
#include "adapt_math.h"
#include "adapt_compute.h"
#include "adapt_meshtopology_lite.h"

#ifndef ADAPT_GRADIENTRECONSTRUCTION_H
#define ADAPT_GRADIENTRECONSTRUCTION_H


std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Lite(PartObject* RePa,
                                                             std::map<int,std::vector<double> > ghosts,
                                                             int Nel,
                                                             int variable,
                                                             int nvariables,
                                                             MPI_Comm comm,
                                                             int extrap);

std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Vrt_Lite(PartObject* RePa,
                                                             std::map<int,std::vector<double> > Uval,
                                                             std::map<int,std::vector<double> > ghosts,
                                                             int Nel,
                                                             int variable,
                                                             int nvariable,
                                                             MPI_Comm comm,
                                                             int extrap);

std::map<int,std::vector<double> > ComputedUdx_LSQ_LS_US3D(PartObject* RePa,
                                                           std::map<int,std::vector<double> > ghosts,
                                                           int Nel,
                                                           int variable,
                                                           int nvariables,
                                                           MPI_Comm comm,
                                                           int extrap);
                                                           
#endif
