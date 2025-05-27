#include "adapt.h"
#include "adapt_partition.h"
#include "adapt_io.h"

#ifndef ADAPT_PAROPS_H
#define ADAPT_PAROPS_H

using namespace std;




inline MPI_Datatype get_mpi_datatype(const int &)
{
    return MPI_INT;
}

inline MPI_Datatype get_mpi_datatype(const double &)
{
    return MPI_DOUBLE;
}

template <typename T>
inline MPI_Datatype get_mpi_datatype() {
    return get_mpi_datatype(T());
}

void RedistributeMeshtThroughRoot(std::map<int,std::vector<int> > elements, int nvpelement, MPI_Comm comm);




#endif
