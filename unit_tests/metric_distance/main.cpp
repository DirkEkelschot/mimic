#include <vector>
#include <math.h>
#include "src/adapt_metricdistance.h"

int main(int argc, char** argv)
{

    clock_t start_total = clock();
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);


    std::vector<double> p0(3,0.0);
    std::vector<double> p1(3,1.0);
  

    std::vector<double> metric0(6,0.0);
    metric0[0] = 1.0;
    metric0[1] = 0.0;
    metric0[2] = 0.0;
    metric0[3] = 1.0;
    metric0[4] = 0.0;
    metric0[5] = 1.0;
    std::vector<double> metric1(6,0.0);
    metric1[0] = 1.0;
    metric1[1] = 0.0;
    metric1[2] = 0.0;
    metric1[3] = 1.0;
    metric1[4] = 0.0;
    metric1[5] = 1.0;

    std::vector<std::vector<double> > edge(2);
    edge[0] = p0;
    edge[1] = p1;

    std::vector<std::vector<double> > metrics(2);
    metrics[0] = metric0;
    metrics[1] = metric1;

    MetricDistance(edge,metrics);

}