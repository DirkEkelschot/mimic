#include "adapt_repartition.h"
#include "adapt_math.h"
#include "adapt_compute.h"
#include "adapt_meshtopology_lite.h"

#ifndef ADAPT_RECONGRAD_LITE_H
#define ADAPT_RECONGRAD_LITE_H

std::map<int,std::vector<double> > ComputedUdx_LSQ_LS_US3D_Lite(RepartitionObject* RePa,
                                                           PrismTetraTrace* trace,
                                                           std::map<int,std::vector<double> > ghosts,
                                                           int Nel,
                                                           int variable,
                                                           int approxOrder,
                                                           MPI_Comm comm);

std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Lite(RepartitionObject* RePa,
                                                            std::map<int,std::vector<double> > U,
                                                             PrismTetraTrace* trace,
                                                             std::map<int,std::vector<double> > ghosts,
                                                             int Nel,
                                                             int variable,
                                                             int approxOrder,
                                                             MPI_Comm comm,
                                                             int print);
#endif
