#include "adapt_parops.h"
#include "adapt_operations.h"
#include "adapt_parmetisstate_lite.h"


void RedistributeMeshtThroughRoot(std::map<int,std::vector<int> > elements, int nvpelement, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    int nelements = elements.size();
    //============================================================================================
    //============================================================================================
    //============================================================================================
    //============================================================================================
    std::vector<int> ien_tetra(nelements*nvpelement,0);
    std::map<int,std::vector<int> >::iterator itmiv;
    int cnt=0;
    for(itmiv=elements.begin();itmiv!=elements.end();itmiv++)
    {
        int nv = itmiv->second.size();
        for(int i=0;i<nv;i++)
        {
            ien_tetra[cnt+i] = itmiv->second[i];
        }
        cnt=cnt+nv;
    }

    DistributedParallelState* elements_dist = new DistributedParallelState(nelements*nvpelement,comm);
    int nelements_total_tmp                      = elements_dist->getNel();
    int nelements_total = (int) nelements_total_tmp/nvpelement;
    ParallelState* elements_ideal_dist = new ParallelState(nelements_total,comm);
    int nelements_ideal = elements_ideal_dist->getNloc(world_rank);

    std::vector<int> ien_tetra_glob;

    std::cout << "nelements_ideal " << nelements << " " << nelements_ideal << " " << nelements_total_tmp << " "   << nelements_total << std::endl; 
    if(world_rank == 0)
    {
        ien_tetra_glob = std::vector<int>(nelements_total*nvpelement);
    }
    else
    {
        ien_tetra_glob = std::vector<int>(1);
    }

    MPI_Gatherv(ien_tetra.data(),
                nelements*nvpelement,
                MPI_INT,
                ien_tetra_glob.data(),
                elements_dist->getNlocs(),
                elements_dist->getOffsets(),
                MPI_INT, 0, comm);


    std::vector<int> ien_tetra_ideal(nelements_ideal*nvpelement);
    std::vector<int> nlocs_new(world_size);
    std::vector<int> offsets_new(world_size);

    for(int i=0;i<world_size;i++)
    {
        nlocs_new[i]   = elements_ideal_dist->getNlocs()[i]*nvpelement;
        offsets_new[i] = elements_ideal_dist->getOffsets()[i]*nvpelement;
    }

    MPI_Scatterv(&ien_tetra_glob.data()[0],
                nlocs_new.data(),
                offsets_new.data(),
                MPI_INT,
                &ien_tetra_ideal.data()[0],
                nelements_ideal*nvpelement,
                MPI_INT, 0, comm);

    std::map<int,std::vector<int> > elements_ideal;

    int offset_el_ideal = elements_ideal_dist->getOffsets()[world_rank];
    for(int i=0;i<nelements_ideal;i++)
    {
        int gelid = i+offset_el_ideal;
        //std::cout << world_rank << " " << gelid << std::endl;
        std::vector<int> row(nvpelement,0);
        for(int j=0;j<nvpelement;j++)
        {
            row[j] = ien_tetra_ideal[i*nvpelement+j];
        }
        elements_ideal[gelid] = row;
    }

    std::vector<int> elTypes(3,0);
    elTypes[0] = 1;

    ParallelState_Parmetis_Lite* pstate_parmetis = new ParallelState_Parmetis_Lite(elements_ideal, 
                                elTypes,  
                                comm);

    //ParallelState_Parmetis* pstate_parmetis2 = new ParallelState_Parmetis(ien,comm,8);
//
    int nloc = elements_ideal.size();
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {pstate_parmetis->getNcommonNodes()};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj_par      = NULL;
    idx_t *adjncy_par    = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {2};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.1};
    real_t *ubvec    = ubvec_;

    std::vector<int> elmwgt = pstate_parmetis->getElmWgt();
    
    int np           = world_size;
    idx_t ncon_[]    = {1};
    idx_t *ncon      = ncon_;
    real_t *tpwgts   = new real_t[np*ncon[0]];

    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }

    idx_t nparts_[] = {np};
    idx_t *nparts = nparts_;
    int* part_arr = new int[nloc];
    real_t itr_[]    = {1.05};
    real_t *itr = itr_;

    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;

    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist().data(),
                          pstate_parmetis->getEptr().data(),
                          pstate_parmetis->getEind().data(),
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&comm);

    for(int i=0;i<nloc;i++)
    {
        std::cout << "RANK " << world_rank  << " iloc = " << i << " :: (" << xadj_par[i] << " " << xadj_par[i+1] << ") ";
        for(int j=xadj_par[i];j<xadj_par[i+1];j++)
        {
            std::cout << adjncy_par[j] << " ";
        }
        std::cout << std::endl;
    }

    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist().data(),
                         xadj_par,
                         adjncy_par,
                         elmwgt.data(), NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);

    // std::cout << tetra_id << " prims vs. " << prism_id << std::endl;

    //============================================================================================
    //============================================================================================
    //============================================================================================
    //============================================================================================
}


