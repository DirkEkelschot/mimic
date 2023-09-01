#include "adapt_prismtetratrace.h"
#include "adapt_distri_parstate.h"
PrismTetraTrace::PrismTetraTrace(MPI_Comm comm,
                                 std::map<int,std::vector<int> > ife, 
                                 std::map<int,int> iet,
                                 int nElem,
                                 int nFace,
                                 int nVert)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::map<int,int>::iterator itmii;
    std::map<int,std::vector<int> >::iterator itmiv;

    std::vector<int> iet_loc(nElem,0);
    std::vector<int> iet_glob(nElem,0);
    for(itmii=iet.begin();itmii!=iet.end();itmii++)
    {
        int el_id       = itmii->first;
        int type        = itmii->second; 
        iet_loc[el_id]  = type;
    }

    MPI_Allreduce(&iet_loc.data()[0], 
                  &iet_glob.data()[0], 
                  nElem, MPI_INT, MPI_SUM, comm);

    //Change this so that it reduces to root and it scatters the result of trace
    
    int trace_loc = 0;

    std::vector<int> elem_l;
    std::vector<int> elem_r;
    std::vector<int> trace_f;

    for(itmiv=ife.begin();itmiv!=ife.end();itmiv++)
    {
        int face_id = itmiv->first;

        int elem_0  = itmiv->second[0];
        int elem_1  = itmiv->second[1];
        int type_0  = iet_glob[elem_0];
        int type_1  = iet_glob[elem_1];

        if(elem_0<nElem && elem_1<nElem)
        {
            if(type_0==2 && type_1==6)
            {
                elem_l.push_back(elem_0);
                elem_r.push_back(elem_1);
                trace_f.push_back(face_id);
                trace_loc++;
            }
            if(type_0==6 && type_1==2)
            {
                elem_l.push_back(elem_1);
                elem_r.push_back(elem_0);
                trace_f.push_back(face_id);
                trace_loc++;
            }
        }   
    }

    int nTraceF = trace_f.size();

    DistributedParallelState* trace_comm = new DistributedParallelState(nTraceF,comm);
    int nTraceF_glob = trace_comm->getNel();
    std::vector<int> elem_l_glob(nTraceF_glob,0);
    std::vector<int> elem_r_glob(nTraceF_glob,0);
    std::vector<int> trace_f_glob(nTraceF_glob,0);

    MPI_Allgatherv(elem_l.data(),
                   nTraceF,
                   MPI_INT,
                   elem_l_glob.data(),
                   trace_comm->getNlocs(),
                   trace_comm->getOffsets(),
                   MPI_INT,comm);
    
    MPI_Allgatherv(elem_r.data(),
                   nTraceF,
                   MPI_INT,
                   elem_r_glob.data(),
                   trace_comm->getNlocs(),
                   trace_comm->getOffsets(),
                   MPI_INT,comm);

    MPI_Allgatherv(trace_f.data(),
                   nTraceF,
                   MPI_INT,
                   trace_f_glob.data(),
                   trace_comm->getNlocs(),
                   trace_comm->getOffsets(),
                   MPI_INT,comm);


    for(int i=0;i<nTraceF_glob;i++)
    {
        int trace_fid = trace_f_glob[i];
        
        if(trace.find(trace_fid)==trace.end())
        {
            std::vector<int> traceElem(2);
            traceElem[0]                = elem_l_glob[i];
            traceElem[1]                = elem_r_glob[i];
            trace[trace_fid]            = traceElem;
        }
    }

    // Test
    // if(rank==0)
    // {
    //     for(itmiv=trace.begin();itmiv!=trace.end();itmiv++)
    //     {
    //         std::cout << itmiv->first << " (" << itmiv->second[0] << " " << itmiv->second[1] << ")" << std::endl;
    //     }
    // }
}

std::map<int,std::vector<int> > PrismTetraTrace::GetTrace()
{
    return trace;
}
