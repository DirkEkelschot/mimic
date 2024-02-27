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
    std::vector<int> rank_l;
    std::vector<int> rank_r;
    std::vector<int> elem_l;
    std::vector<int> elem_r;
    std::vector<int> face_vf;
    std::vector<int> trace_f;

    for(itmiv=ife.begin();itmiv!=ife.end();itmiv++)
    {
        int face_id = itmiv->first;

        int elem_0  = itmiv->second[0];
        int elem_1  = itmiv->second[1];

        int type_0  = iet_glob[elem_0];
        int type_1  = iet_glob[elem_1];

        int rank_0  = element2rank[elem_0];
        int rank_1  = element2rank[elem_1];
        std::vector<int> trave_v(3,0);
        if(elem_0<nElem && elem_1<nElem)
        {
            if(type_0==2 && type_1!=2)
            {
                rank_l.push_back(rank_0);
                rank_r.push_back(rank_1);

                elem_l.push_back(elem_0);
                elem_r.push_back(elem_1);

                trace_f.push_back(face_id);

                for(int q=0;q<3;q++)
                {
                    face_vf.push_back(ifn[face_id][q]);
                }
                trace_loc++;
            }
            if(type_0!=2 && type_1==2)
            {
                rank_l.push_back(rank_1);
                rank_r.push_back(rank_0);

                elem_l.push_back(elem_1);
                elem_r.push_back(elem_0);

                trace_f.push_back(face_id);

                for(int q=0;q<3;q++)
                {
                    face_vf.push_back(ifn[face_id][q]);
                }

                trace_loc++;
            }
        }   
    }

    int nTraceF  = trace_f.size();
    int nTraceFV = face_vf.size();
    DistributedParallelState* trace_comm = new DistributedParallelState(nTraceF,comm);
    DistributedParallelState* traceV_comm = new DistributedParallelState(nTraceFV,comm);

    int nTraceF_glob = trace_comm->getNel();
    int nTraceFV_glob = traceV_comm->getNel();

    std::vector<int> rank_l_glob(nTraceF_glob,0);
    std::vector<int> rank_r_glob(nTraceF_glob,0);
    std::vector<int> elem_l_glob(nTraceF_glob,0);
    std::vector<int> elem_r_glob(nTraceF_glob,0);
    std::vector<int> trace_f_glob(nTraceF_glob,0);
    std::vector<int> trace_fv_glob(nTraceFV_glob,0);

    MPI_Allgatherv(rank_l.data(),
                   nTraceF,
                   MPI_INT,
                   rank_l_glob.data(),
                   trace_comm->getNlocs(),
                   trace_comm->getOffsets(),
                   MPI_INT,comm);

    MPI_Allgatherv(rank_r.data(),
                   nTraceF,
                   MPI_INT,
                   rank_r_glob.data(),
                   trace_comm->getNlocs(),
                   trace_comm->getOffsets(),
                   MPI_INT,comm);

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

    MPI_Allgatherv(face_vf.data(),
                   nTraceFV,
                   MPI_INT,
                   trace_fv_glob.data(),
                   traceV_comm->getNlocs(),
                   traceV_comm->getOffsets(),
                   MPI_INT,comm);

    int ref  = 100;
    int vref = 100;

   
    
    for(int i=0;i<nTraceF_glob;i++)
    {
        int trace_fid = trace_f_glob[i];
        
        if(trace_elems.find(trace_fid)==trace_elems.end())
        {
            std::map<int,int> tr;
            tr[elem_l_glob[i]] = rank_l_glob[i];
            tr[elem_r_glob[i]] = rank_r_glob[i];
            std::vector<int> LeftRight(2,0);
            LeftRight[0] = elem_l_glob[i];
            LeftRight[1] = elem_r_glob[i];
            std::vector<int> fvs(3,0);
            std::vector<int> refs(3,0);
            for(int k=0;k<3;k++)
            {
                fvs[k] = trace_fv_glob[i*3+k];

                if(unique_trace_verts.find(fvs[k])==unique_trace_verts.end())
                {  
                    unique_trace_verts[fvs[k]] = vref;
                    refs[k]=vref;

                    vref = vref + 1;
                }
                else
                {
                    int vreff = unique_trace_verts[fvs[k]];
                    refs[k] = vreff;
                }


            }

            FaceSharedPtr RefTraceFacePointer = std::shared_ptr<NekFace>(new NekFace(refs));
            RefTraceFacePointer->SetFaceID(trace_fid);
            RefTraceFacePointer->SetFaceLeftElement(LeftRight[0]);
            RefTraceFacePointer->SetFaceRightElement(LeftRight[1]);
            m_RefTraceFaces.insert(RefTraceFacePointer);
            

            trace_LR_elem[trace_fid]    = LeftRight;
            trace_elems[trace_fid]      = tr;
            trace_verts[trace_fid]      = fvs;
            trace_ref[trace_fid]        = ref;
            ref++;
        }
    }
    //std::cout << "trace_ref " << trace_verts.size() << std::endl;

}

// void PrismTetraTrace::GetRequiredPrisms()
// {
// }
FaceSetPointer PrismTetraTrace::GetRefTraceFaceSet()
{
    return m_RefTraceFaces;
}
std::map<int,int> PrismTetraTrace::GetUniqueTraceVerts2RefMap()
{
    return unique_trace_verts;
}

std::map<int,std::map<int,int> > PrismTetraTrace::GetTrace()
{
    return trace_elems;
}

std::map<int,std::vector<int> > PrismTetraTrace::GetTraceVerts()
{
    return trace_verts;
}

std::map<int,std::vector<int> > PrismTetraTrace::GetLeftRightElements()
{
    return trace_LR_elem;
}

std::map<int,int > PrismTetraTrace::GetTraceRef()
{
    return trace_ref;
}
