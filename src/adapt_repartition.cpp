#include "adapt_repartition.h"

RepartitionObject::RepartitionObject(mesh* meshInput,
                        std::map<int,std::vector<int> > elements,
                        std::map<int,std::vector<int> > trace,
                        std::map<int,std::vector<double> > data,
                        MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    //GetSharedTraces(elements,meshInput->ife,meshInput->if_ref,meshInput->iet,meshInput->element2rank,comm);
    
    std::map<int,std::vector<int> > elements_update = GetOptimalDistributionSchedule(elements,comm);

    DeterminePartitionLayout(elements_update,comm);
}

void RepartitionObject::GetSharedTraces(std::map<int,std::vector<int> > elements,
                                        std::map<int,std::vector<int> > ife,
                                        std::map<int,int > if_ref,
                                        std::map<int,int > iet,
                                        std::vector<int> element2rank,
                                        MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::map<int,int> sharedFaces;
    std::map<int,std::vector<int> >::iterator itmiv;
    int shf = 0;
    for(itmiv=ife.begin();itmiv!=ife.end();itmiv++)
    {
        int gfid        = itmiv->first;
        int fref        = if_ref[gfid];
        int el0         = itmiv->second[0];
        int el1         = itmiv->second[1];
        int eltype0     = iet[el0];
        int eltype1     = iet[el1];
        int r0          = element2rank[el0];
        int r1          = element2rank[el1];

        if(fref==2)
        {
            if(r0==rank && r1!=rank)
            {
                sharedFaces[gfid] = r0;
                shf++;
            }
            if(r0!=rank && r1==rank)
            {
                sharedFaces[gfid] = r1;
                shf++;
            }
        }
    }
}

std::map<int,std::vector<int> > RepartitionObject::GetOptimalDistributionSchedule(std::map<int,std::vector<int> > elements, MPI_Comm comm)
{

    int world_size;
    MPI_Comm_size(comm, &world_size);
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<int> iniDist(world_size,0);
    std::vector<int> red_iniDist(world_size,0);
    std::vector<int> ini_offsetDist(world_size,0);
    
    int nElem = elements.size();
    int nVpEl = elements[0].size();
    for(int i=0;i<world_size;i++)
    {        
        if(i==world_rank)
        {
            iniDist[i] = elements.size();
        }
    }
    MPI_Allreduce(iniDist.data(), red_iniDist.data(), world_size, MPI_INT, MPI_SUM, comm);

    int offsetEl = 0;
    
    for(int i=0;i<world_size;i++)
    {
        ini_offsetDist[i]   = offsetEl;
        offsetEl            = offsetEl + red_iniDist[i];
    }

    int nElemTotal          = offsetEl;
    
    int size                = world_size;
    int optimalSize         = int(nElemTotal/size) + ( world_rank < nElemTotal%size );
    int NtoRecv             = 0;
    int NtoSend             = 0;
    
    if(nElem>optimalSize)
    {
        NtoSend = nElem-optimalSize;
    }
    if(nElem<optimalSize)
    {
        NtoRecv = optimalSize-nElem;
    }

    int* toR_red_update         = new int[world_size];
    int* toS_red_update         = new int[world_size];
    int* toS_red                = new int[world_size];
    int* toR_red                = new int[world_size];
    int* optiSize               = new int[world_size];
    int* optiSize_red           = new int[world_size];
    int* old_nelem              = new int[world_size];
    int* old_nelem_red          = new int[world_size];
    int* toS                    = new int[world_size];
    int* toR                    = new int[world_size];
    

    for(int i=0;i<world_size;i++)
    {
        toS[i]              = 0;
        toR[i]              = 0;
        optiSize[i]         = 0;
        optiSize_red[i]     = 0;
        toS_red[i]          = 0;
        toR_red[i]          = 0;
        toS_red_update[i]   = 0;
        toR_red_update[i]   = 0;
        old_nelem[i]        = 0;
        old_nelem_red[i]    = 0;

        if(i==world_rank)
        {
            optiSize[i]     = optimalSize;
            toS[i]          = NtoSend;
            toR[i]          = NtoRecv;
            old_nelem[i]    = nElem;

        }
    }
    MPI_Allreduce(optiSize,            optiSize_red, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(toS,                      toS_red, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(toR,                      toR_red, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&old_nelem[0],  &old_nelem_red[0], world_size, MPI_INT, MPI_SUM, comm);
    int sent;
    int sendUpdate;
    
    int opti_Rank,nelem_Rank;
    int ns = 0;
    int nr = 0;
    int nc = 0;
    for(int i=0;i<world_size;i++)
    {
        opti_Rank  = optiSize_red[i];
        nelem_Rank = old_nelem_red[i];
        
        if(nelem_Rank>opti_Rank)
        {
            toS_red_update[i] = old_nelem_red[i]-optiSize_red[i];
            toR_red_update[i] = 0;
            
            ns++;
        }
        if(opti_Rank>nelem_Rank)
        {
            toS_red_update[i] = 0;
            toR_red_update[i] = optiSize_red[i]-old_nelem_red[i];
            nr++;
        }
        if(opti_Rank==nelem_Rank)
        {
            toS_red_update[i] = 0;
            toR_red_update[i] = 0;
            nc++;
        }
    }
    
    int* Psending   = new int[ns];
    int* Nsending   = new int[ns];
    int* Preceiving = new int[nr];
    int* Nreceiving = new int[nr];

    std::set<int> DoNotUpdatePartition;
    
    std::map<int,std::vector<int> > recvRa;
    std::map<int,std::vector<int> > recvNe;
    
    std::map<int,std::vector<int> > sendRa;
    std::map<int,std::vector<int> > sendNe;
    
    int r = 0;
    int s = 0;
    
    for(int i=0;i<world_size;i++)
    {
        if(toS_red_update[i]>0)
        {
            Psending[s] = i;
            Nsending[s] = toS_red_update[i];
            s++;
        }
        if(toR_red_update[i]>0)
        {
            Preceiving[r] = i;
            Nreceiving[r] = toR_red_update[i];
            r++;
        }
        if(toS_red_update[i]==0 && toR_red_update[i]==0)
        {
            DoNotUpdatePartition.insert(i);
        }
    }
    
    int adv         = 0;
    int st          = 0;
    int residual    = 0;
    int Psend;
    while(adv<ns)
    {
        int dist  = Nsending[adv];
        
        if(dist!=0)
        {
            Psend = Psending[adv];
        }

        std::vector<int> toRank;
        std::vector<int> NtoRank;

        int itte=0;
        while(dist!=0)
        {

            int PtoS = Preceiving[st];

            toRank.push_back(PtoS);
            
            if(residual != 0)
            {
                
                if(dist>residual)
                {
                    dist = dist - residual;
                    NtoRank.push_back(residual);
                    residual = 0;
                    st++;
                }
                else if(dist<residual)
                {
                    NtoRank.push_back(dist);
                    residual = residual - dist;
                    dist     = 0;
                }
                else if(dist==residual)
                {
                    NtoRank.push_back(dist);
                    residual = residual - dist;
                    dist     = 0;
                    st++;
                }
            }
            else if(dist>Nreceiving[st])
            {
                dist = dist - Nreceiving[st];
                NtoRank.push_back(Nreceiving[st]);
                st++;
            }
            else if(dist<Nreceiving[st])
            {
                NtoRank.push_back(dist);
                residual   = Nreceiving[st]-dist;
                dist       = 0;
            }
            else if(dist==Nreceiving[st])
            {
                NtoRank.push_back(dist);
                residual = Nreceiving[st]-dist;
                st++;
                dist = 0;
            }
            
            itte++;
        }
        
        //std::cout << std::endl;
//
        sendRa[Psend]=toRank;
        sendNe[Psend]=NtoRank;
        
        adv++;
    }
    

    
    std::map<int,std::vector<int> >::iterator its;

    for(its=sendRa.begin();its!=sendRa.end();its++)
    {
        for(int q=0;q<its->second.size();q++)
        {
            recvRa[its->second[q]].push_back(its->first);
            recvNe[its->second[q]].push_back(sendNe[its->first][q]);
        }
    }

    //==============================================================
    //==============================================================
    //==============================================================

    std::map<int,std::vector<int> > elements_update;
    std::map<int,std::vector<int> >::iterator itmiv;
    if(DoNotUpdatePartition.find(world_rank)!=DoNotUpdatePartition.end())
    {
        std::map<int,std::vector<int> >::iterator itmiv;
        int eloc = 0;
        for(itmiv=elements.begin();itmiv!=elements.end();itmiv++)
        {
            int gElId = itmiv->first;

            for(int j=0;j<elements.size();j++)
            {
                elements_update[gElId] = itmiv->second;
            }
        }
    }
    if(sendRa.find(world_rank)!=sendRa.end())
    {
        std::vector<int> toRanks    = sendRa[world_rank];
        std::vector<int> NeltoRanks = sendNe[world_rank];

        std::vector<std::vector<int> > elIDs;
        std::vector<std::vector<int> > elNvs;
        std::vector<std::vector<int> > ien_to_send;
        
        for(int i=0;i<toRanks.size();i++)
        {
            int Nel = NeltoRanks[i];

            std::vector<int> rowElID(Nel);
            elIDs.push_back(rowElID);
            std::vector<int> rowNvEl(Nel);
            elNvs.push_back(rowNvEl);

            std::vector<int> rows_ien(Nel*nVpEl);
            ien_to_send.push_back(rows_ien);
        }
                
        int offPrank = 0;
        int cntv     = 0;
        int t        = 0;
        int nuloc    = 0;
        int uloc     = 0;
        int u        = 0;
        int cc       = 0;

        for(itmiv=elements.begin();itmiv!=elements.end();itmiv++)
        {
            int nelPrank    = NeltoRanks[cc];
            int gElId       = itmiv->first;
            int nv_el       = itmiv->second.size();
            std::vector<int> ien_row(nv_el);

            if(u<toS_red_update[world_rank])
            {
                if(t<(nelPrank))
                {

                    elIDs[cc][t] = gElId;
                    elNvs[cc][t] = nv_el;

                    for(int j=0;j<nv_el;j++)
                    {
                        ien_to_send[cc][nVpEl*t+j] = itmiv->second[j];
                    }

                    if(t==nelPrank-1)
                    {
                        t = 0;
                        cc=cc+1;
                    }
                    else
                    {
                        t=t+1;
                    }
                }
            }
            else
            {
                std::vector<int> update_row(nv_el,0);
                for(int j=0;j<nv_el;j++)
                {
                    update_row[j] = itmiv->second[j];
                }

                elements_update[gElId] = update_row;

            }

            u++;
        }

        int acull = 0;
        for(int i=0;i<toRanks.size();i++)
        {
            int dest     = toRanks[i];
            int n_Elem   = NeltoRanks[i];
            int n_Vrt    = ien_to_send[i].size();

            std::vector<int> ElIDs = elIDs[i];
            std::vector<int> NvEl  = elNvs[i];
            std::vector<int> Elvec = ien_to_send[i];

            std::cout << "elIDs " << elNvs[i].size() << " " << n_Elem << " " << toS_red_update[world_rank] << " " << n_Vrt << std::endl;

            MPI_Send(&n_Elem         ,      1,     MPI_INT, dest, dest,          comm);
            MPI_Send(&ElIDs[0]       , n_Elem,     MPI_INT, dest, dest*123000,   comm);
            MPI_Send(&NvEl[0]        , n_Elem,     MPI_INT, dest, dest*492000,   comm);

            MPI_Send(&n_Vrt          ,      1,     MPI_INT, dest, dest*984000,   comm);
            MPI_Send(&Elvec[0]       ,  n_Vrt,     MPI_INT, dest, dest*1968000,  comm);

            acull = acull + n_Vrt;
        }

    }
    if(recvRa.find(world_rank)!=recvRa.end())
    {
        std::vector<int > expFromRank = recvRa[world_rank];
        
        std::map<int,std::vector<int> > recvd_elem_ids;
        std::map<int,std::vector<int> > recvd_elem_nvs;
        std::map<int,std::vector<int> > recvd_elem_vrt_ids;

        for(int i=0;i<expFromRank.size();i++)
        {
            int origin = expFromRank[i];

            int n_Elem;
            MPI_Recv(&n_Elem,   1, MPI_INT, origin, world_rank, comm, MPI_STATUS_IGNORE);   
            std::vector<int> rcvd_el_ids(n_Elem);
            MPI_Recv(&rcvd_el_ids[0], n_Elem, MPI_INT, origin, world_rank*123000, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_nvs(n_Elem);
            MPI_Recv(&rcvd_el_nvs[0], n_Elem, MPI_INT, origin, world_rank*492000, comm, MPI_STATUS_IGNORE);


            int n_Vrt;
            MPI_Recv(&n_Vrt,    1, MPI_INT, origin, world_rank*984000, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_vrt_ids(n_Vrt);
            MPI_Recv(&rcvd_el_vrt_ids[0], n_Vrt, MPI_INT, origin, world_rank*1968000, comm, MPI_STATUS_IGNORE);
            
            recvd_elem_ids[origin]   = rcvd_el_ids;
            recvd_elem_nvs[origin]   = rcvd_el_nvs;
            recvd_elem_vrt_ids[origin]  = rcvd_el_vrt_ids;
        }

        
        //============================================================
        for(itmiv=elements.begin();itmiv!=elements.end();itmiv++)
        {
            int gElId               = itmiv->first;
            elements_update[gElId]  = itmiv->second;
        }

        int ntot = elements.size();
        for(itmiv=recvd_elem_ids.begin();itmiv!=recvd_elem_ids.end();itmiv++)
        {
            int nel = itmiv->second.size();

            //int fromRank = collit->first;

            int accul = 0;
            for(int q=0;q<nel;q++)
            {
                int gElId = itmiv->second[q];
                int nvpel = recvd_elem_nvs[itmiv->first][q];
                std::vector<int> elements_update_row(nvpel,0);

                for(int s=0;s<nvpel;s++)
                {
                    elements_update_row[s] = recvd_elem_vrt_ids[itmiv->first][accul+s];
                }

                elements_update[gElId] = elements_update_row;

                accul = accul + nvpel;

            }
            
            ntot    = ntot + nel;
        }
    }

    return elements_update;
}




void RepartitionObject::DeterminePartitionLayout(std::map<int,std::vector<int> > elements,
                                                MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nrow = elements.size();
    int nvpEL = elements.begin()->second.size();
    int nloc = nrow;

    std::vector<int> elTypes(3,0);
    if(nvpEL == 4)
    {
        elTypes[0] = 1;
    }
    if(nvpEL == 6)
    {
        elTypes[1] = 1;
    }

    std::cout << " elTypes " << elTypes[0] << " " << elTypes[1] << " " << elTypes[2] << " " << nvpEL << std::endl;
    
    ParallelState_Parmetis_Lite* pstate_parmetis = new ParallelState_Parmetis_Lite(elements,  elTypes, comm);

    //=================================================================
    //=================================================================
    //=================================================================
    
    //ParallelState_Parmetis* pstate_parmetis2 = new ParallelState_Parmetis(ien,comm,8);
//
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
    
    int np           = size;
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
    
    // for(int i=0;i<pstate_parmetis->getElmWgt().size();i++)
    // {
    //     std::cout << "pstate_parmetis->getElmWgt().data() " << pstate_parmetis->getElmWgt().data()[i] << std::endl;
    // }

    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist().data(),
                          pstate_parmetis->getEptr().data(),
                          pstate_parmetis->getEind().data(),
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&comm);


    // for(int i=0;i<nloc;i++)
    // {
    //     std::cout << "RANK " << rank  << " iloc = " << i << " :: (" << xadj_par[i] << " " << xadj_par[i+1] << ") ";
    //     for(int j=xadj_par[i];j<xadj_par[i+1];j++)
    //     {
    //         std::cout << adjncy_par[j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist().data(),
                         xadj_par,
                         adjncy_par,
                         pstate_parmetis->getElmWgt().data(), NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr, &comm);
    
    part = std::vector<int>(elements.size(),0);
    int nElemTotal = pstate_parmetis->getNtotalElem();
    part_global = std::vector<int>(nElemTotal,0);

    for(int i=0;i<elements.size();i++)
    {
        part[i] = part_arr[i];
    }
    
    MPI_Allgatherv(&part.data()[0],
                   nloc, MPI_INT,
                   &part_global.data()[0],
                   pstate_parmetis->getNlocs().data(),
                   pstate_parmetis->getElmdist().data(),
                   MPI_INT,comm);
    
   if(rank == 0)
   {
        std::cout << "Succesfully found a redistribution of the elements." << std::endl;
   }
   delete[] xadj_par;
   delete[] adjncy_par;
}

// destructor
RepartitionObject::~RepartitionObject()
{
    
}

