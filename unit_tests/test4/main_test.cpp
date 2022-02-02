#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include <iomanip>

int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
   
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    ParArray<int>*   ien_test = ReadDataSetFromFileInParallel<int>("../test_mesh/test_mesh.h5","ien",comm,info);
    ParallelState_Parmetis* parm_test = new ParallelState_Parmetis(ien_test,comm,8);
    ParallelState* ienp_test = new ParallelState(ien_test->getNglob(),comm);
    int nrow = ien_test->getNrow();
    int nloc = nrow;
    
    idx_t numflag_[] = {0};
    idx_t *numflag = numflag_;
    idx_t ncommonnodes_[] = {4};
    idx_t *ncommonnodes = ncommonnodes_;
    int edgecut      = 0;
    idx_t *xadj_par      = NULL;
    idx_t *adjncy_par    = NULL;
    idx_t options_[] = {0, 0, 0};
    idx_t *options   = options_;
    idx_t wgtflag_[] = {0};
    idx_t *wgtflag   = wgtflag_;
    real_t ubvec_[]  = {1.1};
    real_t *ubvec    = ubvec_;

    idx_t *elmwgt = NULL;

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
    real_t *itr      = itr_;
    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;
    
    int status = ParMETIS_V3_Mesh2Dual(parm_test->getElmdist(), parm_test->getEptr(), parm_test->getEind(), numflag, ncommonnodes, &xadj_par, &adjncy_par, &comm);

//    ParMETIS_V3_AdaptiveRepart(parm_test->getElmdist(),
//                           xadj_par, adjncy_par,
//                               elmwgt, adjwgt,
//                       vsize, wgtflag,
//                   numflag, ncon, nparts,
//                   tpwgts, ubvec, itr, options,
//                   &edgecut, part_arr, &comm);
    
    
    ParMETIS_V3_PartKway(parm_test->getElmdist(),
                         xadj_par,
                         adjncy_par,
                         elmwgt, NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);
    
    ParArray<int>* part = new ParArray<int>(ien_test->getNglob(),1,comm);
    Array<int>* part_global = new Array<int>(ien_test->getNglob(),1);

    part->data = part_arr;
    int* xadj = xadj_par;
    int* adjcny = adjncy_par;
    
    MPI_Allgatherv(&part->data[0],
                   nloc, MPI_INT,
                   &part_global->data[0],
                   ienp_test->getNlocs(),
                   ienp_test->getOffsets(),
                   MPI_INT,comm);
    
    std::map<int,std::vector<int> > elms_to_send_to_ranks;
    int not_on_rank = 0;
    for(int i=0;i<part->getNrow();i++)
    {
        int p_id  = part->getVal(i,0);
    
        
        int el_id = part->getOffset(world_rank)+i;
        
        if(p_id!=world_rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id);
            
            not_on_rank++;
        }
    }
        
    ScheduleObj* sObj = DoScheduling(elms_to_send_to_ranks,comm);
    
    std::map<int, std::set<int> >::iterator its;
    
    if(world_rank == 0)
    {
        std::map<int,std::vector<int> >::iterator itm;
        for(itm=elms_to_send_to_ranks.begin();itm!=elms_to_send_to_ranks.end();itm++)
        {
            std::cout << "Rank  " << world_rank << " sends the following elements: ";
            for(int q=0;q<itm->second.size();q++)
            {
                std::cout <<itm->second[q] << " ";
            }
            std::cout << " to rank " << itm->first << std::endl;
        }
        std::cout << "Send Schedule " << sObj->SendFromRank2Rank.size() <<std::endl;

        for(its=sObj->SendFromRank2Rank.begin();its!=sObj->SendFromRank2Rank.end();its++)
        {
            std::cout << world_rank << " -> " << its->first << " ellies ";
            std::set<int>::iterator it;
            for(it=its->second.begin();it!=its->second.end();it++)
            {
                std::cout << *it << " ";
            }
            std::cout << std::endl;
        }
//
//
        std::cout << "Recv Schedule " << sObj->SendFromRank2Rank.size() <<std::endl;

        for(its=sObj->RecvRankFromRank.begin();its!=sObj->RecvRankFromRank.end();its++)
        {
            std::cout << world_rank << " -> " << its->first << " ellies ";
            std::set<int>::iterator it;
            for(it=its->second.begin();it!=its->second.end();it++)
            {
                std::cout << *it << " ";
            }
            std::cout << std::endl;
        }
    }
    
    MPI_Finalize();
}
