#include "adapt_prismtetratrace.h"
#include "adapt_distri_parstate.h"


PrismTetraTrace::PrismTetraTrace(MPI_Comm comm,
                                 std::vector<int> element2rank,
                                 std::map<int,std::vector<int> > ife,
                                 std::map<int,std::vector<int> > ifn,
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
    std::vector<int> face_vf;
    std::vector<int> trace_f;

    
    for(itmiv=ife.begin();itmiv!=ife.end();itmiv++)
    {
        int face_id = itmiv->first;

        int elem_0  = itmiv->second[0];
        int elem_1  = itmiv->second[1];
        if(elem_0<nElem && elem_1<nElem)
        {
            int type_0  = iet_glob[elem_0];
            int type_1  = iet_glob[elem_1];
            int rank_0  = element2rank[elem_0];
            int rank_1  = element2rank[elem_1];
            std::vector<int> trave_v(3,0);
            if(elem_0<nElem && elem_1<nElem)
            {
                if(type_0==2 && type_1!=2)
                {
                    trace_f.push_back(face_id);

                    for(int q=0;q<3;q++)
                    {
                        face_vf.push_back(ifn[face_id][q]);
                    }
                    trace_loc++;
                }
                if(type_0!=2 && type_1==2)
                {
                    trace_f.push_back(face_id);

                    for(int q=0;q<3;q++)
                    {
                        face_vf.push_back(ifn[face_id][q]);
                    }
                    trace_loc++;
                }
            }  
         } 
    }
    
    int nTraceF  = trace_f.size();
    int nTraceFV = face_vf.size();
    DistributedParallelState* trace_comm = new DistributedParallelState(nTraceF,comm);
    DistributedParallelState* traceV_comm = new DistributedParallelState(nTraceFV,comm);

    int nTraceF_glob = trace_comm->getNel();
    int nTraceFV_glob = traceV_comm->getNel();
    std::vector<int> trace_f_glob(nTraceF_glob,0);
    std::vector<int> trace_fv_glob(nTraceFV_glob,0);

    MPI_Allgatherv(trace_f.data(),
                   nTraceF,
                   MPI_INT,
                   trace_f_glob.data(),
                   trace_comm->getNlocs(),
                   trace_comm->getOffsets(),
                   MPI_INT,comm);

    MPI_Allgatherv(face_vf.data(),
                   nTraceFV,
                   MPI_INT,
                   trace_fv_glob.data(),
                   traceV_comm->getNlocs(),
                   traceV_comm->getOffsets(),
                   MPI_INT,comm);

    int ref  = 100;
    int vref = 100;

    std::set<int> trace_elems;
    
    for(int i=0;i<nTraceF_glob;i++)
    {
        int trace_fid = trace_f_glob[i];
        
        if(trace_elems.find(trace_fid)==trace_elems.end())
        {
            trace_elems.insert(trace_fid);
            
            for(int k=0;k<3;k++)
            {
                int vid = trace_fv_glob[i*3+k];

                if(unique_trace_verts.find(vid)==unique_trace_verts.end())
                {  
                    unique_trace_verts[vid] = vref;

                    vref = vref + 1;
                }
                else
                {
                    int vreff = unique_trace_verts[vid];
                }
            }
        }
        
    }

    trace_elems.clear();
    trace_fv_glob.clear();
    trace_f_glob.clear();

}


std::map<int,int> PrismTetraTrace::GetUniqueTraceVerts2RefMap()
{
    return unique_trace_verts;
}


PrismTetraTrace::~PrismTetraTrace()
{
    unique_trace_verts.clear();
}
