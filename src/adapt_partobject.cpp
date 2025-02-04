#include "adapt_partobject.h"
#include "adapt_compute.h"
#include "adapt_geometry.h"


PartObject::PartObject(mesh* meshInput,
                    std::map<int,std::vector<int> > Elem2Vert_i,
                    std::map<int,std::vector<int> > Elem2Face_i,
                    std::map<int,std::vector<int> > Elem2Elem_i,
                    std::map<int,std::vector<double> > Elem2Data_i,
                    std::map<int,int> Elem2Type_i,
                    int nAdjLayer,
                    bool reconstruct_ifn,
                    MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    int rank;
    MPI_Comm_rank(comm, &rank);

    Ne_glob = meshInput->nElem;
    Nf_glob = meshInput->nFace;
    Nv_glob = meshInput->nVert;
   
    std::map<int,std::vector<int> > Elem2Vert_uniform;
    std::map<int,std::vector<int> > Elem2Face_uniform;
    std::map<int,std::vector<int> > Elem2Elem_uniform;
    std::map<int,std::vector<double> > Elem2Data_uniform;
    std::map<int,int> Elem2Type_uniform;

    clock_t start, end;
    double dur_max,time_taken;
    
    //===============================================================================================================
    start = clock();
    GetOptimalDistributionSchedule(Elem2Vert_i, 
                                    Elem2Face_i, 
                                    Elem2Elem_i, 
                                    Elem2Data_i,
                                    Elem2Type_i,
                                    Elem2Vert_uniform, 
                                    Elem2Face_uniform, 
                                    Elem2Elem_uniform, 
                                    Elem2Data_uniform,
                                    Elem2Type_uniform,
                                    comm);

    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << std::setprecision(3) << "Finished uniformly distributing Element2Entity data in:                              " << std::fixed
        << dur_max;
        std::cout << " sec " << std::endl;
    }
    //===============================================================================================================
    // DeterminePartitionLayout calls Parmetis and determines a part_map for each rank. 
    // part_map is a std::map<key,val> where key is the global element ID and val is the 
    // rank ID that this element should live on.
    // Next to that, DeterminePartitionLayout communicates the partitioning layout to root (rank = 0).
    // Root holds a global vector that is of length number of elements where for each element an rank id is allocated.
    start = clock();
    DeterminePartitionLayout(Elem2Vert_uniform,meshInput->element2rank,meshInput->elTypes,comm);
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << std::setprecision(3) << "Finished running Parmetis on the uniformly distributed Element2Entity data in:       " << std::fixed
        << dur_max;
        std::cout << " sec " << std::endl;
    }
    //===============================================================================================================
    // DetermineElement2ProcMap communicates the uniformly distributed Element2Entity maps to
    // apppropriate ranks based on part_map.
    start = clock();
    DetermineElement2ProcMap(Elem2Vert_uniform, 
                             Elem2Face_uniform, 
                             Elem2Elem_uniform, 
                             Elem2Data_uniform,
                             Elem2Type_uniform,
                             meshInput->xcn,
                             Nf_glob,
                             Nv_glob, 
                             comm);
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << std::setprecision(3) << "Based on Parmetis, communicating the Element2Entity data to appropriate ranks in:    " << std::fixed
        << dur_max;
        std::cout << " sec " << std::endl;
    }
    //===============================================================================================================
    start = clock();
    // std::map<int,std::vector<int> > Rank2Elem = CommunicateAdjacencyInfoLocalPartition(comm);
    std::map<int,std::vector<int> > Rank2Elem = CommunicateAdjacencyInfoLocalPartition(comm);
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << std::setprecision(3) << "For each rank query root where the adjacent data lives:                              " << std::fixed
        << dur_max;
        std::cout << " sec " << std::endl;
    }
    //===============================================================================================================
    start = clock();
    GenerateFace2ElementMap(meshInput->ife, comm);
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << std::setprecision(3) << "The Face2Element map is generated based on the Element2Entity data:                  " << std::fixed
        << dur_max;
        std::cout << " sec " << std::endl;
    }
    //===============================================================================================================
    start = clock();
    if(reconstruct_ifn==true)
    {
        // generate m_Face2Vert without communication since the mesh does not contain hexes.
        GenerateFace2VertexMap(comm);

    }
    else
    {
        // generate m_Face2Vert with communication because the mesh contains hexes.
         getFace2VertexPerPartitionMap(meshInput->ifn,
                                       meshInput->if_Nv,
                                       Nf_glob,
                                       comm);
    }
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << std::setprecision(3) << "The Face2Vertex map is generated based on the Element2Entity data:                   " << std::fixed
        << dur_max;
        std::cout << " sec " << std::endl;
    }
    //===============================================================================================================
    start = clock();
    GenerateTraceMap();
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << std::setprecision(3) << "Generating the trace map:                                                            " << std::fixed
        << dur_max;
        std::cout << " sec " << std::endl;
    }


    //===============================================================================================================
    for(int i=0;i<nAdjLayer;i++)
    {
        // std::cout << "m_Elem2Face.size() "<< rank << " -- " << m_Elem2Face.size() << std::endl;
        start = clock();
        std::map<int,std::vector<int> > adjElem2Face = getAdjacentElementLayer(meshInput->xcn,comm);
        end = clock();
        time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
        MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(rank == 0)
        {
            std::cout << std::setprecision(3) << "Extending the scheme by incorporating adjacent data layer : " << i << "                        " << std::fixed
            << dur_max;
            std::cout << " sec " << std::endl;
        }


        // this updates m_Face2Vert
        start = clock();
        updateFace2EntityPerPartition(adjElem2Face, 
                                      meshInput->ifn,
                                      meshInput->if_Nv,  
                                      comm);
        end = clock();
        time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
        MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(rank == 0)
        {
            std::cout << std::setprecision(3) << "After communicating additional adjacent data, update m_Face2Vert : " << i << "                 " << std::fixed
            << dur_max;
            std::cout << " sec " << std::endl;
        }
    }
    //===============================================================================================================

}



PartObject::~PartObject()
{

}



/**

* @brief This function takes in the the Element2Entity data structures that are read in through the 
* specified CFD solver IO routines available in MIMIC.

* This function then computes a uniform distribution of the Element2Entity data structures
* and takes care of communicating the data to each rank.
* Ultimately, this function fills m_Elem2Vert, m_Elem2Face, m_Elem2Elem, m_Elem2Data and m_Elem2Type
* with a uniformly distributed set of data.
 */

void PartObject::GetOptimalDistributionSchedule(std::map<int,std::vector<int> > Elem2Vert_i,
                                                std::map<int,std::vector<int> > Elem2Face_i,
                                                std::map<int,std::vector<int> > Elem2Elem_i,
                                                std::map<int,std::vector<double> > Elem2Data_i,
                                                std::map<int,int> Elem2Type_i,
				                                std::map<int,std::vector<int> >& Elem2Vert_uniform,
                                                std::map<int,std::vector<int> >& Elem2Face_uniform,
                                                std::map<int,std::vector<int> >& Elem2Elem_uniform,
                                                std::map<int,std::vector<double> >& Elem2Data_uniform,
                                                std::map<int,int>& Elem2Type_uniform,
                                                MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<int> iniDist(world_size,0);
    std::vector<int> red_iniDist(world_size,0);
    std::vector<int> ini_offsetDist(world_size,0);
    
    int nElem = Elem2Vert_i.size();
    int ndata = 0;

    Elem2Vert_uniform.clear();
    std::map<int,std::vector<int> > Elem2Vert_uniform_v2;

  
    // Check if this rank has at least an element allocated here.

    if(nElem != 0)
    {
        ndata = Elem2Data_i.begin()->second.size();
    }

    for(int i=0;i<world_size;i++)
    {        
        if(i==world_rank)
        {
            iniDist[i] = Elem2Vert_i.size();
        }
    }

    // Communicate the current distribution to all ranks and figuring out the uniform distribution.
    // "optimal" here means uniform.

    // Until the next comment, all this code aims to determine the communication pattern that is
    // necessary to uniformly lay out the Element2Entity data.


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

    // The actual communication of the data takes place below.
    // This is based on the communication patter that is established in the code above.
    //==============================================================
    //==============================================================
    //==============================================================

    std::map<int,std::vector<int> >::iterator itmiv;
    if(DoNotUpdatePartition.find(world_rank)!=DoNotUpdatePartition.end())
    {
        std::map<int,std::vector<int> >::iterator itmiv;
        int eloc = 0;
        for(itmiv=Elem2Vert_i.begin();itmiv!=Elem2Vert_i.end();itmiv++)
        {
            int gElId = itmiv->first;
            int gtype = Elem2Type_i[gElId];

            Elem2Vert_uniform[gElId]      = itmiv->second;
            Elem2Vert_uniform_v2[gElId]   = itmiv->second;
            Elem2Type_uniform[gElId]      = gtype;
            

            Elem2Face_uniform[gElId]     = Elem2Face_i[gElId];
            Elem2Elem_uniform[gElId]     = Elem2Elem_i[gElId];
            Elem2Data_uniform[gElId]     = Elem2Data_i[gElId];

        }
    }
    if(sendRa.find(world_rank)!=sendRa.end())
    {
        std::vector<int> toRanks    = sendRa[world_rank];
        std::vector<int> NeltoRanks = sendNe[world_rank];

        std::vector<std::vector<int> > elIDs;
        std::vector<std::vector<int> > elTypes;
        std::vector<std::vector<int> > elNvs;
        std::vector<std::vector<int> > elNfs;
        std::vector<std::vector<int> > elNdata;
        std::vector<std::vector<int> > ien_to_send(toRanks.size());
        std::vector<std::vector<int> > ief_to_send(toRanks.size());
        std::vector<std::vector<int> > iee_to_send(toRanks.size());
        std::vector<std::vector<double> > data_to_send;

        for(int i=0;i<toRanks.size();i++)
        {
            int Nel = NeltoRanks[i];

            std::vector<int> rowElID(Nel);
            elIDs.push_back(rowElID);
            std::vector<int> rowTypes(Nel);
            elTypes.push_back(rowTypes);
            std::vector<int> rowNvEl(Nel);
            elNvs.push_back(rowNvEl);
            std::vector<int> rowNfEl(Nel);
            elNfs.push_back(rowNfEl);
            std::vector<int> rowNdataEl(Nel);
            elNdata.push_back(rowNdataEl);

            std::vector<double> rows_data(Nel*ndata);
            data_to_send.push_back(rows_data);
        }
                
        int offPrank = 0;
        int cntv     = 0;
        int t        = 0;
        int nuloc    = 0;
        int uloc     = 0;
        int u        = 0;
        int cc       = 0;

        for(itmiv=Elem2Vert_i.begin();itmiv!=Elem2Vert_i.end();itmiv++)
        {
            int nelPrank    = NeltoRanks[cc];

            int gElId               = itmiv->first;
            int update_type         = Elem2Type_i[gElId];
            int nv_el               = itmiv->second.size();
            int nf_el               = Elem2Face_i[gElId].size();
            int ndata               = Elem2Data_i[gElId].size();

            std::vector<int> ien_row(nv_el);

            if(u<toS_red_update[world_rank])
            {
                if(t<(nelPrank))
                {
                    elTypes[cc][t]  = update_type;
                    elIDs[cc][t]    = gElId;
                    elNvs[cc][t]    = nv_el;
                    elNfs[cc][t]    = nf_el;
                    elNdata[cc][t]  = ndata;

                    for(int j=0;j<nv_el;j++)
                    {
                        ien_to_send[cc].push_back(itmiv->second[j]);
                    }

                    for(int j=0;j<nf_el;j++)
                    {

                        ief_to_send[cc].push_back(Elem2Face_i[gElId][j]);
                        iee_to_send[cc].push_back(Elem2Elem_i[gElId][j]);
                    
                    }


                    for(int j=0;j<ndata;j++)
                    {
                        data_to_send[cc][ndata*t+j] = Elem2Data_i[gElId][j];
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
                std::vector<int> update_e2v_row(nv_el,0);

                for(int j=0;j<nv_el;j++)
                {
                    update_e2v_row[j] = itmiv->second[j];
                }


                std::vector<int> update_e2f_row(nf_el,0);
                std::vector<int> update_e2e_row(nf_el,0);

                for(int j=0;j<nf_el;j++)
                {
                    update_e2f_row[j] = Elem2Face_i[gElId][j];
                    update_e2e_row[j] = Elem2Elem_i[gElId][j];
                }

                std::vector<double> update_data_row(ndata,0);

                for(int j=0;j<ndata;j++)
                {
                    update_data_row[j] = Elem2Data_i[gElId][j];
                }
                Elem2Type_uniform[gElId]      = update_type;
                Elem2Vert_uniform[gElId]      = update_e2v_row;
                Elem2Face_uniform[gElId]      = update_e2f_row;
                Elem2Elem_uniform[gElId]      = update_e2e_row;
                Elem2Data_uniform[gElId]      = update_data_row;
            }
          
            u++;
        }
        
        int acull = 0;
        for(int i=0;i<toRanks.size();i++)
        {
            int dest     = toRanks[i];
            int n_Elem   = NeltoRanks[i];
            int n_Vrt    = ien_to_send[i].size();
            int n_Fce    = ief_to_send[i].size();
            int n_Data   = data_to_send[i].size();

            std::vector<int> ElIDs   = elIDs[i];
            std::vector<int> ElTypes = elTypes[i];
            std::vector<int> NvEl    = elNvs[i];
            std::vector<int> NfEl    = elNfs[i];
            std::vector<int> NdataEl = elNdata[i];
            std::vector<int> El2v    = ien_to_send[i];
            std::vector<int> El2f    = ief_to_send[i];
            std::vector<int> El2e    = iee_to_send[i];
            std::vector<double> El2d = data_to_send[i];

            MPI_Send(&n_Elem         ,      1,     MPI_INT, dest, dest,             comm);
            MPI_Send(&ElIDs[0]       , n_Elem,     MPI_INT, dest, dest*123000,      comm);
            MPI_Send(&ElTypes[0]     , n_Elem,     MPI_INT, dest, dest*128000,      comm);
            MPI_Send(&NvEl[0]        , n_Elem,     MPI_INT, dest, dest*492000,      comm);
            MPI_Send(&NfEl[0]        , n_Elem,     MPI_INT, dest, dest*984000,      comm);

            MPI_Send(&n_Vrt          ,      1,     MPI_INT, dest, dest*1968000,     comm);
            MPI_Send(&El2v[0]       ,   n_Vrt,     MPI_INT, dest, dest*3936000,     comm);

            MPI_Send(&n_Fce          ,      1,     MPI_INT, dest, dest*7872000,     comm);
            MPI_Send(&El2f[0]       ,   n_Fce,     MPI_INT, dest, dest*15744000,    comm);
            MPI_Send(&El2e[0]       ,   n_Fce,     MPI_INT, dest, dest*31488000,    comm);

            MPI_Send(&n_Data        ,      1,      MPI_INT, dest, dest*62976000,    comm);
            MPI_Send(&El2d[0]       ,   n_Data,    MPI_DOUBLE, dest, dest*125952,   comm);
            MPI_Send(&NdataEl[0]    ,   n_Elem,    MPI_INT, dest, dest*251904,      comm);
            acull = acull + n_Vrt;
        }

    }
    if(recvRa.find(world_rank)!=recvRa.end())
    {
        std::vector<int > expFromRank = recvRa[world_rank];
        
        std::map<int,std::vector<int> > recvd_elem_ids;
        std::map<int,std::vector<int> > recvd_elem_type;
        std::map<int,std::vector<int> > recvd_elem_nvs;
        std::map<int,std::vector<int> > recvd_elem_nfs;
        std::map<int,std::vector<int> > recvd_elem_ndata;
        std::map<int,std::vector<int> > recvd_elem_vert_ids;
        std::map<int,std::vector<int> > recvd_elem_face_ids;
        std::map<int,std::vector<int> > recvd_elem_elem_ids;
        std::map<int,std::vector<double> > recvd_elem_data_ids;
        for(int i=0;i<expFromRank.size();i++)
        {
            int origin = expFromRank[i];

            int n_Elem;
            MPI_Recv(&n_Elem,   1, MPI_INT, origin, world_rank, comm, MPI_STATUS_IGNORE);   
            std::vector<int> rcvd_el_ids(n_Elem,0);
            MPI_Recv(&rcvd_el_ids[0], n_Elem, MPI_INT, origin, world_rank*123000, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_type(n_Elem,0);
            MPI_Recv(&rcvd_el_type[0], n_Elem, MPI_INT, origin, world_rank*128000, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_nvs(n_Elem,0);
            MPI_Recv(&rcvd_el_nvs[0], n_Elem, MPI_INT, origin, world_rank*492000, comm, MPI_STATUS_IGNORE);  
            std::vector<int> rcvd_el_nfs(n_Elem,0);
            MPI_Recv(&rcvd_el_nfs[0], n_Elem, MPI_INT, origin, world_rank*984000, comm, MPI_STATUS_IGNORE);


            int n_Vrt;
            MPI_Recv(&n_Vrt,    1, MPI_INT, origin, world_rank*1968000, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_vrt_ids(n_Vrt,0);
            MPI_Recv(&rcvd_el_vrt_ids[0], n_Vrt, MPI_INT, origin, world_rank*3936000, comm, MPI_STATUS_IGNORE);


            int n_Fce;
            MPI_Recv(&n_Fce,    1, MPI_INT, origin, world_rank*7872000, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_face_ids(n_Fce,0);
            MPI_Recv(&rcvd_el_face_ids[0], n_Fce, MPI_INT, origin, world_rank*15744000, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_elem_ids(n_Fce,0);
            MPI_Recv(&rcvd_el_elem_ids[0], n_Fce, MPI_INT, origin, world_rank*31488000, comm, MPI_STATUS_IGNORE);


            int n_Data;
            MPI_Recv(&n_Data,    1, MPI_INT, origin, world_rank*62976000, comm, MPI_STATUS_IGNORE);
            std::vector<double> rcvd_el_data_ids(n_Data,0.0);
            MPI_Recv(&rcvd_el_data_ids[0], n_Data, MPI_DOUBLE, origin, world_rank*125952, comm, MPI_STATUS_IGNORE);
            std::vector<int> rcvd_el_ndatas(n_Elem,0);
            MPI_Recv(&rcvd_el_ndatas[0], n_Elem, MPI_INT, origin, world_rank*251904, comm, MPI_STATUS_IGNORE);
            
            recvd_elem_ids[origin]       = rcvd_el_ids;
            recvd_elem_type[origin]      = rcvd_el_type;
            recvd_elem_nvs[origin]       = rcvd_el_nvs;
            recvd_elem_vert_ids[origin]  = rcvd_el_vrt_ids;
            recvd_elem_nfs[origin]       = rcvd_el_nfs;
            recvd_elem_face_ids[origin]  = rcvd_el_face_ids;
            recvd_elem_elem_ids[origin]  = rcvd_el_elem_ids;
            recvd_elem_ndata[origin]     = rcvd_el_ndatas;
            recvd_elem_data_ids[origin]  = rcvd_el_data_ids;

        }

        
        //============================================================
        for(itmiv=Elem2Vert_i.begin();itmiv!=Elem2Vert_i.end();itmiv++)
        {
            int gElId                       = itmiv->first;
            int type_update                 = Elem2Type_i[gElId]; 
            Elem2Vert_uniform[gElId]              = itmiv->second;
            Elem2Face_uniform[gElId]              = Elem2Face_i[gElId];
            Elem2Elem_uniform[gElId]              = Elem2Elem_i[gElId];
            Elem2Data_uniform[gElId]              = Elem2Data_i[gElId];
            Elem2Type_uniform[gElId]              = type_update; 
        }

        int ntot = Elem2Vert_i.size();
        for(itmiv=recvd_elem_ids.begin();itmiv!=recvd_elem_ids.end();itmiv++)
        {
            int nel  = itmiv->second.size();
            
            //int fromRank = collit->first;

            int accul_v = 0;
            int accul_f = 0;
            int accul_d = 0;
            for(int q=0;q<nel;q++)
            {
                int gElId = itmiv->second[q];
                int type  = recvd_elem_type[itmiv->first][q];
                int nvpel = recvd_elem_nvs[itmiv->first][q];
                int nfpel = recvd_elem_nfs[itmiv->first][q];
                int ndata = recvd_elem_ndata[itmiv->first][q];
                std::vector<int> Elem2Vert_uniform_row(nvpel,0);
                std::vector<int> Elem2Face_uniform_row(nfpel,0);
                std::vector<int> Elem2Elem_uniform_row(nfpel,0);
                std::vector<double> Elem2Data_uniform_row(ndata,0);
                for(int s=0;s<nvpel;s++)
                {
                    Elem2Vert_uniform_row[s] = recvd_elem_vert_ids[itmiv->first][accul_v+s];
                }
                for(int s=0;s<nfpel;s++)
                {
                    Elem2Face_uniform_row[s]    = recvd_elem_face_ids[itmiv->first][accul_f+s];
                    Elem2Elem_uniform_row[s] = recvd_elem_elem_ids[itmiv->first][accul_f+s];
                }

                for(int s=0;s<ndata;s++)
                {
                    Elem2Data_uniform_row[s]     = recvd_elem_data_ids[itmiv->first][accul_d+s];
                }


                Elem2Type_uniform[gElId]          = type;
                Elem2Vert_uniform[gElId]          = Elem2Vert_uniform_row;
                Elem2Face_uniform[gElId]          = Elem2Face_uniform_row;
                Elem2Elem_uniform[gElId]          = Elem2Elem_uniform_row;
                Elem2Data_uniform[gElId]          = Elem2Data_uniform_row;

                accul_v = accul_v + nvpel;
                accul_f = accul_f + nfpel;
                accul_d = accul_d + ndata;
            }

            ntot    = ntot + nel;
        }
    }
}





void PartObject::DeterminePartitionLayout(std::map<int,std::vector<int> > Elem2Vert_uniform,
                                         std::vector<int> element2rank,
                                         std::vector<int> elTypes,
                                         MPI_Comm comm)
{
    int root = 0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nrow    = Elem2Vert_uniform.size();
    int nvpEL   = Elem2Vert_uniform.begin()->second.size();
    int nloc    = nrow;

    //std::cout << " elTypes " << elTypes[0] << " " << elTypes[1] << " " << elTypes[2] << " " << nvpEL << std::endl;
    
    ParallelState_Parmetis_Lite* pstate_parmetis = new ParallelState_Parmetis_Lite(Elem2Vert_uniform,  elTypes, comm);

    //=================================================================
    //=================================================================
    //=================================================================
    
    //ParallelState_Parmetis* pstate_parmetis2 = new ParallelState_Parmetis(ien,comm,8);
//
    idx_t numflag_[]        = {0};
    idx_t *numflag          = numflag_;
    idx_t ncommonnodes_[]   = {pstate_parmetis->getNcommonNodes()};
    idx_t *ncommonnodes     = ncommonnodes_;
    int edgecut             = 0;
    idx_t *xadj_par         = NULL;
    idx_t *adjncy_par       = NULL;
    idx_t options_[]        = {0, 0, 0};
    idx_t *options          = options_;
    idx_t wgtflag_[]        = {2};
    idx_t *wgtflag          = wgtflag_;
    real_t ubvec_[]         = {1.1};
    real_t *ubvec           = ubvec_;

    std::vector<int> elmwgt = pstate_parmetis->getElmWgt();
    
    int np                  = size;
    idx_t ncon_[]           = {1};
    idx_t *ncon             = ncon_;
    real_t *tpwgts          = new real_t[np*ncon[0]];

    for(int i=0; i<np*ncon[0]; i++)
    {
        tpwgts[i] = 1.0/np;
    }

    idx_t nparts_[]         = {np};
    idx_t *nparts           = nparts_;
    std::vector<int> part_arr(nloc,0);
    real_t itr_[]           = {1.05};
    real_t *itr             = itr_;

    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;

    ParMETIS_V3_Mesh2Dual(pstate_parmetis->getElmdist().data(),
                          pstate_parmetis->getEptr().data(),
                          pstate_parmetis->getEind().data(),
                          numflag,ncommonnodes,
                          &xadj_par,&adjncy_par,&comm);

    ParMETIS_V3_PartKway(pstate_parmetis->getElmdist().data(),
                         xadj_par,
                         adjncy_par,
                         pstate_parmetis->getElmWgt().data(), NULL, wgtflag, numflag,
                         ncon, nparts,
                         tpwgts, ubvec, options,
                         &edgecut, part_arr.data(), &comm);

    m_part = std::vector<int>(Elem2Vert_uniform.size(),0);
    int nElemTotal = pstate_parmetis->getNtotalElem();
    

    std::vector<int> part_arr_ell_id(nloc,0);
    std::vector<int> part_arr_ell_Rid(nloc,0);

    std::map<int,std::vector<int> >::iterator itmiv;
    int i = 0;

    for(itmiv=Elem2Vert_uniform.begin();itmiv!=Elem2Vert_uniform.end();itmiv++)
    {
        int gid                 = itmiv->first;
        part_arr_ell_id[i]      = gid;
        part_arr_ell_Rid[i]     = part_arr[i];
        m_partMap[gid]          = part_arr[i];
        i++;
    }

    std::vector<int> m_partGlobalRoot_vec;
    std::vector<int> m_partGlobalRoot_Rid_vec;
    if ( rank == root) 
    {
        m_partGlobalRoot_vec.resize(nElemTotal,0);
        m_partGlobalRoot_Rid_vec.resize(nElemTotal,0);
    } 

    MPI_Gatherv(&part_arr_ell_id.data()[0],
                    nloc, MPI_INT,
                    &m_partGlobalRoot_vec.data()[0],
                    pstate_parmetis->getNlocs().data(),
                    pstate_parmetis->getElmdist().data(),
                    MPI_INT,root,comm);

    MPI_Gatherv(&part_arr_ell_Rid.data()[0],
                    nloc, MPI_INT,
                    &m_partGlobalRoot_Rid_vec.data()[0],
                    pstate_parmetis->getNlocs().data(),
                    pstate_parmetis->getElmdist().data(),
                    MPI_INT,root,comm);

    for(int i=0;i<m_partGlobalRoot_vec.size();i++)
    {
        int eid = m_partGlobalRoot_vec[i];
        int rid = m_partGlobalRoot_Rid_vec[i];

        m_partGlobalRoot[eid] = rid;
    }

    nElemGlobalPart = m_partGlobalRoot.size();
    // std::cout << "root " << rank << " nElemGlobalPart " << nElemGlobalPart << " " << m_part.size() << std::endl;
    MPI_Bcast(&nElemGlobalPart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] xadj_par;
    delete[] adjncy_par;
}






void PartObject::DetermineElement2ProcMap(std::map<int,std::vector<int> >   Elem2Vert_uniform, 
                                        std::map<int,std::vector<int> >     Elem2Face_uniform,
                                        std::map<int,std::vector<int> >     Elem2Elem_uniform,
                                        std::map<int,std::vector<double> >  Elem2Data_uniform,
                                        std::map<int,int>                   Elem2Type_uniform, 
                                        std::map<int,std::vector<double> >  vertices_i,
                                        int Nf_glob,
                                        int Nv_glob, 
                                        MPI_Comm comm)
{

    std::map<int,std::vector<int> > elms_to_send_to_ranks;
    std::map<int,std::vector<int> > nvPerElms_to_send_to_ranks;
    std::map<int,std::vector<int> > nfPerElms_to_send_to_ranks;
    std::map<int,std::vector<int> > ndataPerElms_to_send_to_ranks;
    int floc_tmp=0;
    int vloc_tmp=0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    int el_id;
    int p_id;
    int v_id;
    Vert V;
    std::vector<Vert> part_verts;
    std::vector<std::vector<int> > part_elem2verts;

    std::map<int,std::vector<int> > vertIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > faceIDs_to_send_to_ranks;
    std::map<int,std::vector<int> > elemIDs_to_send_to_ranks;
    std::map<int,std::vector<double> > data_to_send_to_ranks;
    std::map<int,std::vector<int> > elmtype_to_send_to_ranks;
    ParallelState* ife_pstate = new ParallelState(Nf_glob,comm);
    ParallelState* vert_pstate = new ParallelState(Nv_glob,comm);

    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_face;
    std::vector<int> faceIDs_on_rank;
    
    std::vector<int> vertIDs_on_rank;
    std::vector<int> part_v;
    
    int r     = 0;
    int lv_id = 0;
    int lf_id = 0;
    int f_id  = 0;
    int ea_id = 0;

    int not_on_rank=0;
    int on_rank = 0;
    int* new_V_offsets = new int[size];
    int* new_F_offsets = new int[size];
    for(i=0;i<size;i++)
    {
        new_V_offsets[i] = vert_pstate->getOffsets()[i]-1;
        new_F_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
    
    int nvPerEl;
    int nfPerEl;
    int ndataPerEl;
    int tett = 0;
    std::map<int,int>::iterator itmii;
    
    std::vector<double> loc_r_data;
    std::vector<int> loc_r_Ndata;
    std::vector<int> loc_r_elem;
    std::vector<int> loc_r_nf_elem;
    std::vector<int> loc_r_nv_elem;
    for(itmii=m_partMap.begin();itmii!=m_partMap.end();itmii++)
    {
        el_id   = itmii->first;
        p_id    = itmii->second;
        int el_type = Elem2Type_uniform[el_id];
    
        nvPerEl = Elem2Vert_uniform[el_id].size();
        nfPerEl = Elem2Face_uniform[el_id].size();
        ndataPerEl = Elem2Data_uniform[el_id].size();

        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            nvPerElms_to_send_to_ranks[p_id].push_back(nvPerEl);
            nfPerElms_to_send_to_ranks[p_id].push_back(nfPerEl);
            ndataPerElms_to_send_to_ranks[p_id].push_back(ndataPerEl);
            elmtype_to_send_to_ranks[p_id].push_back(el_type);

            // if(el_type == 0)
            // {
            //     std::cout << "Already wrong " << std::endl;
            // }

            //====================Hybrid=======================
            for(int k=0;k<nvPerEl;k++)//This works for hexes.
            {
                v_id = Elem2Vert_uniform[el_id][k];

                vertIDs_to_send_to_ranks[p_id].push_back(v_id);

            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            //std::cout << "Check Elem2Elem_uniform :: ";
            for(int k=0;k<nfPerEl;k++)//This works for hexes.
            {
                f_id = Elem2Face_uniform[el_id][k];
                ea_id = Elem2Elem_uniform[el_id][k];

                //std::cout << ea_id << " ";
                faceIDs_to_send_to_ranks[p_id].push_back(f_id);
                elemIDs_to_send_to_ranks[p_id].push_back(ea_id);
            }

            //std::cout << std::endl;
            for(int k=0;k<ndataPerEl;k++)//This works for hexes.
            {
                double dataV = Elem2Data_uniform[el_id][k];
                data_to_send_to_ranks[p_id].push_back(dataV);
            }
            //====================Hybrid=======================
            not_on_rank++;
        }
        else // Here we are storing the actual vertices/elements that are required by the current rank.
        {
            std::vector<int> elem;

            for(int k=0;k<nvPerEl;k++)// looping over the vertices for element "i".
            {
                v_id = Elem2Vert_uniform[el_id][k];
                
                
                //elem.push_back(v_id);
                if(m_vertIDs_on_rank.find( v_id ) == m_vertIDs_on_rank.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    m_vertIDs_on_rank.insert(v_id);
                    
                    r = FindRank(new_V_offsets,size,v_id);

                    
                    if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                    {
                        rank2req_vert[r].push_back(v_id); // add the vertex id that needs to be requested from rank r.
                    }
                    else
                    {
                        vertIDs_on_rank.push_back(v_id);  // add the vertex to list that is already available on rank.
                        vloc_tmp++;
                    }
                    lv_id++;
                }
            }

            


            // std::cout << "Check Elem2Elem_uniform :: ";
            // for(int k=0;k<nfPerEl;k++)//This works for hexes.
            // {
            //     ea_id = Elem2Elem_uniform[el_id][k];

            //     std::cout << ea_id << " ";
            // }

            // std::cout << std::endl;
            
            loc_r_elem.push_back(el_id);
            loc_r_Ndata.push_back(ndataPerEl);
            for(int k=0;k<ndataPerEl;k++)//This works for hexes.
            {
                double dataV = Elem2Data_uniform[el_id][k];
                //std::cout << "dataV " << dataV << std::endl; 
                loc_r_data.push_back(dataV);
            }

            loc_r_nv_elem.push_back(nvPerEl);
            loc_r_nf_elem.push_back(nfPerEl);
            m_elem2type_on_rank[el_id] = el_type;
            
            on_rank++;
        }
    }
    
    ScheduleObj* part_schedule_elem = DoScheduling(elms_to_send_to_ranks,comm);
    
    std::map<int,std::vector<int> > part_tot_recv_elIDs_map;
    std::map<int,std::vector<int> > part_tot_recv_elNVs_map;
    std::map<int,std::vector<int> > part_tot_recv_elNFs_map;
    
    std::map<int,std::vector<int> > part_tot_recv_tett_map;

    std::map<int,std::vector<double> > part_tot_recv_data_map;
    std::map<int,std::vector<int> > TotRecvElement_IDs_v_map;
    std::map<int,std::vector<int> > TotRecvElement_IDs_f_map;
    std::map<int,std::vector<int> > TotRecvElement_IDs_e_map;
    std::map<int,std::vector<int> > part_tot_recv_elNdata_map;
    std::map<int,std::vector<int> > part_tot_recv_elTypes_map;
    std::map<int,std::vector<int> >::iterator it;
    
    int n_req_recv;
    int n_req_recv_v;
    int n_req_recv_f;
    int n_req_recv_d;
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = elms_to_send_to_ranks.begin(); it != elms_to_send_to_ranks.end(); it++)
            {
                int n_req           = it->second.size();
                int n_req_v         = vertIDs_to_send_to_ranks[it->first].size();
                int n_req_f         = faceIDs_to_send_to_ranks[it->first].size();

                int n_data2send     = data_to_send_to_ranks[it->first].size();
                int dest            = it->first;
                                
                MPI_Send(&n_req  , 1, MPI_INT, dest, dest, comm);
                MPI_Send(&n_req_v, 1, MPI_INT, dest, dest*111, comm);
                MPI_Send(&n_req_f, 1, MPI_INT, dest, dest*222, comm);
                MPI_Send(&n_data2send, 1, MPI_INT, dest, dest*333, comm);

                //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 100+dest*2, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, dest*66666+5555, comm);
                MPI_Send(&nvPerElms_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*33333+7777, comm);
                MPI_Send(&nfPerElms_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*44444+8888, comm);
                MPI_Send(&ndataPerElms_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*55555+8888, comm);
                MPI_Send(&vertIDs_to_send_to_ranks[it->first][0], n_req_v, MPI_INT, dest, 9000+100+dest*2, comm);
                MPI_Send(&faceIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 229000+100+dest*2, comm);
                MPI_Send(&elemIDs_to_send_to_ranks[it->first][0], n_req_f, MPI_INT, dest, 449000+100+dest*2, comm);
                MPI_Send(&data_to_send_to_ranks[it->first][0], n_data2send, MPI_DOUBLE, dest, 339000+100+dest*2, comm);
                MPI_Send(&elmtype_to_send_to_ranks[it->first][0], n_req, MPI_INT, dest, dest*55555+7777, comm);

                //MPI_Send(&vrt_coords_to_send_to_ranks[it->first][0], n_req_v*3, MPI_DOUBLE, dest, 678000+100+dest*2, comm);

                i++;
            }
        }
        else if (part_schedule_elem->SendFromRank2Rank[q].find( rank ) != part_schedule_elem->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_req_recv,   1, MPI_INT, q,     rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_v, 1, MPI_INT, q, rank*111, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_f, 1, MPI_INT, q, rank*222, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&n_req_recv_d, 1, MPI_INT, q, rank*333, comm, MPI_STATUS_IGNORE);

            std::vector<int>    part_recv_el_id(n_req_recv,0);
            std::vector<int>    part_recv_el_type(n_req_recv,0);
            std::vector<int>    part_recv_el_nv(n_req_recv,0);
            std::vector<int>    part_recv_el_nf(n_req_recv,0);
            std::vector<int>    part_recv_el_ndata(n_req_recv,0);

            std::vector<int>    part_recv_vrt_id(n_req_recv_v,0);
            std::vector<int>    part_recv_face_id(n_req_recv_f,0);
            std::vector<int>    part_recv_elem_id(n_req_recv_f,0);

            std::vector<double>    part_recv_vrt_coords(n_req_recv_v*3,0);
            std::vector<double>    part_recv_data(n_req_recv_d,0.0);
            
            MPI_Recv(&part_recv_el_id[0], n_req_recv, MPI_INT, q, rank*66666+5555, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_nv[0], n_req_recv, MPI_INT, q, rank*33333+7777, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_nf[0], n_req_recv, MPI_INT, q, rank*44444+8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_ndata[0], n_req_recv, MPI_INT, q, rank*55555+8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_vrt_id[0],  n_req_recv_v, MPI_INT, q, 9000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_face_id[0], n_req_recv_f, MPI_INT, q, 229000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_elem_id[0], n_req_recv_f, MPI_INT, q, 449000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_data[0], n_req_recv_d, MPI_DOUBLE, q, 339000+100+rank*2, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&part_recv_el_type[0], n_req_recv, MPI_INT, q, rank*55555+7777, comm, MPI_STATUS_IGNORE);

            //MPI_Recv(&part_recv_vrt_coords[0], n_req_recv_d*3, MPI_DOUBLE, q, 678000+100+rank*2, comm, MPI_STATUS_IGNORE);
            
            TotRecvElement_IDs_v_map[q]     = part_recv_vrt_id;
            TotRecvElement_IDs_f_map[q]     = part_recv_face_id;
            TotRecvElement_IDs_e_map[q]     = part_recv_elem_id;

            part_tot_recv_elIDs_map[q]      = part_recv_el_id;
            part_tot_recv_elTypes_map[q]    = part_recv_el_type;

            part_tot_recv_elNVs_map[q]      = part_recv_el_nv;
            part_tot_recv_elNFs_map[q]      = part_recv_el_nf;
            part_tot_recv_elNdata_map[q]    = part_recv_el_ndata;
            part_tot_recv_data_map[q]       = part_recv_data;
        }
    }
    
    std::vector<int> TotRecvElement_IDs;
    std::vector<int> TotRecvElement_NVs;
    std::vector<int> TotRecvElement_NFs;
    std::vector<int> TotRecvElement_Types;
    std::vector<int> TotRecvVerts_IDs;
    std::vector<int> TotRecvFaces_IDs;
    std::vector<int> TotRecvElem_IDs;
    std::vector<double> TotRecvDatas;
    std::vector<int> TotRecvElement_Ndatas;
    std::map<int,std::vector<int> >::iterator totrecv;
    std::map<int,std::vector<double> >::iterator totrecv_double;
    //unpack the element IDs and their corresponding variable values.
    int TotNelem_recv = 0;
    for(totrecv=part_tot_recv_elIDs_map.begin();totrecv!=part_tot_recv_elIDs_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvElement_IDs.push_back(part_tot_recv_elIDs_map[totrecv->first][r]);
            TotRecvElement_NVs.push_back(part_tot_recv_elNVs_map[totrecv->first][r]);
            TotRecvElement_NFs.push_back(part_tot_recv_elNFs_map[totrecv->first][r]);
            TotRecvElement_Types.push_back(part_tot_recv_elTypes_map[totrecv->first][r]);

            TotRecvElement_Ndatas.push_back(part_tot_recv_elNdata_map[totrecv->first][r]);
        }
        TotNelem_recv = TotNelem_recv + totrecv->second.size();
    }

    //unpack the vertex IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_v_map.begin();totrecv!=TotRecvElement_IDs_v_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvVerts_IDs.push_back(TotRecvElement_IDs_v_map[totrecv->first][r]);
        }
    }
    //unpack the face IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_f_map.begin();totrecv!=TotRecvElement_IDs_f_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvFaces_IDs.push_back(TotRecvElement_IDs_f_map[totrecv->first][r]);
        }
    }

     //unpack the face IDs and their corresponding variable values.
    for(totrecv=TotRecvElement_IDs_e_map.begin();totrecv!=TotRecvElement_IDs_e_map.end();totrecv++)
    {
        for(int r=0;r<totrecv->second.size();r++)
        {
            TotRecvElem_IDs.push_back(TotRecvElement_IDs_e_map[totrecv->first][r]);
        }
    }

    for(totrecv_double=part_tot_recv_data_map.begin();totrecv_double!=part_tot_recv_data_map.end();totrecv_double++)
    {
        for(int r=0;r<totrecv_double->second.size();r++)
        {
            TotRecvDatas.push_back(totrecv_double->second[r]);
            //std::cout << "totrecv_double->second[r] " << totrecv_double->second[r] << std::endl;
        }
    }

    int Nel_extra = TotNelem_recv;
    int cnt_v = 0;
    int cnt_f = 0;
    int cnt_data = 0;
    for(int i=0;i<TotNelem_recv;i++)
    {
        std::vector<int> elem;

        int nvPerEl = TotRecvElement_NVs[i];
        int nfPerEl = TotRecvElement_NFs[i];
        int ndatas  = TotRecvElement_Ndatas[i];  

        for(int k=0;k<nvPerEl;k++)
        {
            int v_id_n = TotRecvVerts_IDs[cnt_v+k];
            //elem.push_back(v_id_n);
            
            if(m_vertIDs_on_rank.find( v_id_n ) == m_vertIDs_on_rank.end()) // add the required unique vertex for current rank.
            {
                m_vertIDs_on_rank.insert(v_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                
                r = FindRank(new_V_offsets,size,v_id_n);


                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    //std::cout << "cehck location " << r << " " << v_id_n  << " " << rank<< std::endl;
                    rank2req_vert[r].push_back(v_id_n); // add the vertex id that needs to be requested from rank r.
                }
                else
                {
                    vertIDs_on_rank.push_back(v_id_n);  // add the vertex to list that is already available on rank.
                    vloc_tmp++;
                }
                lv_id++;
            }
        }
        for(int k=0;k<nfPerEl;k++)// looping over the vertices for element "i".
        {
            int f_id_n = TotRecvFaces_IDs[cnt_v+k];
            //int e_id_n = TotRecvElem_IDs[cnt_v+k];
            if(m_faceIDs_on_rank.find( f_id_n ) == m_faceIDs_on_rank.end()) // add the required unique vertex for current rank.
            {
                m_faceIDs_on_rank.insert(f_id_n);
                //unique_verts_on_rank_vec.push_back(v_id);
                
                r = FindRank(new_F_offsets,size,f_id_n);

                if (r!=rank)// if vertex is present on other rank, add it to vertIDs_on_rank map..
                {
                    rank2req_face[r].push_back(f_id_n); // add the vertex id that needs to be requested from rank r.
                }
                else
                {
                    faceIDs_on_rank.push_back(f_id_n);  // add the vertex to list that is already available on rank.
                    floc_tmp++;
                }
                lf_id++;
            }
        }

        // std::vector<double> rowdata(ndatas,0.0);
        // for(int k=0;k<ndatas;k++)
        // {
        //     rowdata[k] = TotRecvDatas[cnt_data+k];
        // }
        
        cnt_v=cnt_v+nvPerEl;
        cnt_f=cnt_f+nfPerEl;
        cnt_data = cnt_data+ndatas;
        //part_elem2verts.push_back(elem);
        on_rank++;
        
    }
    
    // Loop over all received vertex IDs in order to determine the remaining required unique vertices on the current rank.
    
    
    
    // ==========================================================================================
    // ==========================================================================================
    // ==========================================================================================
    
    // At this point we have all the elements that are required on current rank and the vertex ids as well
    // However we are still missing the vertex coordinate data which is spread out equally over the available procs.
    // This rank2req_vert map essentially holds this information by mapping the rank_id from which we need to request a list/vector of vertex ids (hence the name "rank2req_vert" name.
    
    // At this point the perspective changes. When we were figuring out the layout of the elements, we knew the partition ID for each element on the current rank. This means that from the current rank, we needed to send a certain element to another rank since it is more logical to reside there. For the vertices this changes since we just figured out which vertices are required on the current rank. The logic here is first to send for each the current rank a list/vector<int> of vertex IDs that is requested from another rank. The other rank assembles the list of the required coordinates and sends it back.
    
    // ==========================================================================================
    // ==========================================================================================
    // ==========================================================================================
    int m = 0;
    int n_reqstd_ids;
    int n_req_recv_v2;
    
    // This thing needs to revised because for the verts it doesnt work.
    // The current rank does not have the verts_to_send_rank. Instead it has an request list.
    
    ScheduleObj* part_schedule = DoScheduling(rank2req_vert,comm);
    
    std::map<int,std::vector<int> >  reqstd_ids_per_rank;

    
    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_vert.begin(); it != rank2req_vert.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, comm);
                //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                
                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);
            //std::cout << rank << " n_reqstd_ids " << n_reqstd_ids << std::endl;
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2+rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_ids_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nverts;
    std::map<int,std::vector<double> > recv_back_verts;
    std::map<int,std::vector<int> > recv_back_verts_ids;
    int n_recv_back;
        
    int nfound=0;
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
            {
                int nv_send = it->second.size();
                int dest = it->first;
                //double* vert_send = new double[nv_send*3];
                std::vector<double> vert_send(nv_send*3);
                offset_xcn        = vert_pstate->getOffset(rank);
                for(int u=0;u<it->second.size();u++)
                {
                    if(vertices_i.find(it->second[u])==vertices_i.end())
                    {
                        std::cout << "NOT FOUND" << std::endl;
                        std::cout << " it->second[u] " << it->second[u] << " " << dest << " " << rank  << " " << Nv_glob << std::endl;
                        nfound++;
                    }
                    else
                    {
                        vert_send[u*3+0]=vertices_i[it->second[u]][0];
                        vert_send[u*3+1]=vertices_i[it->second[u]][1];
                        vert_send[u*3+2]=vertices_i[it->second[u]][2];
                    }
                    // vert_send[u*3+0]=xcn->getVal(it->second[u]-offset_xcn,0);
                    // vert_send[u*3+1]=xcn->getVal(it->second[u]-offset_xcn,1);
                    // vert_send[u*3+2]=xcn->getVal(it->second[u]-offset_xcn,2);
                }
                
                MPI_Send(&nv_send, 1, MPI_INT, dest, 9876+1000*dest, comm);
                // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);
            
                MPI_Send(&vert_send.data()[0], nv_send*3, MPI_DOUBLE, dest, 9876+dest*8888, comm);
                MPI_Send(&it->second.data()[0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                
                //delete[] vert_send;
            }
        }
        if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876+1000*rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<double> recv_back_arr(n_recv_back*3);
            std::vector<int> recv_back_arr_ids(n_recv_back);
            MPI_Recv(&recv_back_arr.data()[0], n_recv_back*3, MPI_DOUBLE, q, 9876+rank*8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_arr_ids.data()[0], n_recv_back, MPI_INT, q, 8888*9876+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nverts[q]     = n_recv_back;
            recv_back_verts[q]      = recv_back_arr;
            recv_back_verts_ids[q]  = recv_back_arr_ids;
        
         }
    }

    int vfor = 0;
    std::map<int,std::vector<double> >::iterator it_f;
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int c  = 0;
        vfor=vfor+recv_back_Nverts[it_f->first];

    }

    int gvid=0;
    int lvid=0;


    for(m=0;m<vloc_tmp;m++)
    {
        gvid = vertIDs_on_rank[m];

       
        std::vector<double> V(3,0.0);
        V[0] = vertices_i[gvid][0];
        V[1] = vertices_i[gvid][1];
        V[2] = vertices_i[gvid][2];
        m_LocalVertsMap[gvid]  = V;

        m_LocalV2GlobalV[lvid] = gvid;
        m_GlovalV2LocalV[gvid] = lvid;

        lvid++;
    }
    
    m = 0;
    
    for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
    {
        int Nv = recv_back_Nverts[it_f->first];
       
        for(int u=0;u<Nv;u++)
        {
            gvid = rank2req_vert[it_f->first][u];
            
            std::vector<double> V(3);
            V[0] = it_f->second[u*3+0];
            V[1] = it_f->second[u*3+1];
            V[2] = it_f->second[u*3+2];
            m_LocalVertsMap[gvid] = V;

            m_LocalV2GlobalV[lvid] = gvid;
            m_GlovalV2LocalV[gvid] = lvid;
           
            m++;
            lvid++;
        }
    }

    int nLoc_Verts = m_LocalVertsMap.size();

    int glob_v      = 0;
    int loc_v       = 0;
    int glob_f      = 0;
    int glob_e      = 0;
    double varia_v  = 0.0;
    int ndatas      = 0;

    std::map<int,int> LocElem2Ndata;

    for(m=0;m<loc_r_elem.size();m++)
    {
        el_id   = loc_r_elem[m];
        m_ElemSet.insert(el_id);
        nvPerEl = loc_r_nv_elem[m];
        nfPerEl = loc_r_nf_elem[m];
        ndatas  = loc_r_Ndata[m];

        // LocElem2Nv[el_id]    = nvPerEl;
        // LocElem2Nf[el_id]    = nfPerEl;
        LocElem2Ndata[el_id] = ndatas;

        int ndatas = Elem2Data_uniform[el_id].size();
        std::vector<double> rowdata(ndatas,0.0);
        for(int k=0;k<ndatas;k++)
        {
            rowdata[k] = Elem2Data_uniform[el_id][k];
        }

        std::vector<int> tmp_globv;
        std::vector<int> tmp_locv;
        std::vector<int> tmp_globf;
        std::vector<int> tmp_globe;
        std::vector<int> tmp_loce;
        std::vector<double> Pijk(nvPerEl*3);

        for(int p=0;p<nvPerEl;p++)
        {

            glob_v  = Elem2Vert_uniform[el_id][p];
            loc_v   = m_GlovalV2LocalV[glob_v];

            // if(glob_v>vert_pstate->getNel())
            // {
            //     std::cout << "Nel On error " << glob_v << std::endl;
            // }
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            Pijk[p*3+0]     = m_LocalVertsMap[glob_v][0];
            Pijk[p*3+1]     = m_LocalVertsMap[glob_v][1];
            Pijk[p*3+2]     = m_LocalVertsMap[glob_v][2];


            m_GlobalVert2Elem[glob_v].push_back(el_id); //globVerts2globElem[glob_v].push_back(el_id);
            m_Elem2LocalVert[el_id].push_back(loc_v); //globElem2locVerts[el_id].push_back(loc_v);
        }

        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,nvPerEl);

        m_Elem2Centroid[el_id] = Vijk;

        for(int p=0;p<nfPerEl;p++)
        {
            glob_f = Elem2Face_uniform[el_id][p];
            glob_e = Elem2Elem_uniform[el_id][p];

            m_globElem2globFaces[el_id].push_back(glob_f);
            m_globFace2GlobalElements[glob_f].push_back(el_id);

            tmp_globf.push_back(glob_f);
            tmp_globe.push_back(glob_e);
        }
    

        //=======================================
        m_Elem2Vert[el_id]      = tmp_globv;
        m_Elem2Face[el_id]      = tmp_globf;
        m_Elem2Elem[el_id]      = tmp_globe;
        m_Elem2Data[el_id]      = rowdata;

        //=======================================

        tmp_globv.clear();
        tmp_locv.clear();
    }

    int cnv = 0;
    int cnf = 0;
    cnt_data = 0;
    for(m=0;m<Nel_extra;m++)
    {
        
        el_id                             = TotRecvElement_IDs[m];
        m_ElemSet.insert(el_id);
        int el_type                       = TotRecvElement_Types[m];
        nvPerEl                           = TotRecvElement_NVs[m];
        nfPerEl                           = TotRecvElement_NFs[m];
        int ndatas                        = TotRecvElement_Ndatas[m];
        m_elem2type_on_rank[el_id]        = el_type;
        LocElem2Ndata[el_id]              = ndatas;

        // if(el_type == 0)
        // {
        //     std::cout << "wrong " << std::endl;
        // }

        std::vector<int> tmp_globv;
        std::vector<int> tmp_locv;
        std::vector<int> tmp_globf;
        std::vector<int> tmp_globe;
        std::vector<int> tmp_loce;
        std::vector<double> rowdata(ndatas,0.0);

        for(int k=0;k<ndatas;k++)
        {
            rowdata[k] = TotRecvDatas[cnt_data+k];
        }
        cnt_data = cnt_data + ndatas;
        std::vector<double> Pijk(nvPerEl*3);
        for(int p=0;p<nvPerEl;p++)
        {
            glob_v  = TotRecvVerts_IDs[cnv+p];
            loc_v   = m_GlovalV2LocalV[glob_v];

            // if(glob_v>vert_pstate->getNel())
            // {
            //     std::cout << "Nel Extra error " << glob_v << std::endl;
            // }

            m_Elem2LocalVert[el_id].push_back(loc_v);
            m_GlobalVert2Elem[glob_v].push_back(el_id);

            Pijk[p*3+0] = m_LocalVertsMap[glob_v][0];
            Pijk[p*3+1] = m_LocalVertsMap[glob_v][1];
            Pijk[p*3+2] = m_LocalVertsMap[glob_v][2];

            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            
        }

        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,nvPerEl);

        m_Elem2Centroid[el_id] = Vijk;

        for(int p=0;p<nfPerEl;p++)
        {
            glob_f = TotRecvFaces_IDs[cnf+p];
            glob_e = TotRecvElem_IDs[cnf+p];

            m_globElem2globFaces[el_id].push_back(glob_f);
            m_globFace2GlobalElements[glob_f].push_back(el_id);

            tmp_globf.push_back(glob_f);
            tmp_globe.push_back(glob_e);
        }

        //=======================================
        m_Elem2Vert[el_id]    = tmp_globv;
        m_Elem2Face[el_id]    = tmp_globf;
        m_Elem2Elem[el_id]    = tmp_globe;
        m_Elem2Data[el_id]    = rowdata;
        //=======================================
       

        cnv=cnv+nvPerEl;
        cnf=cnf+nfPerEl;

        tmp_globv.clear();
        tmp_locv.clear();
    }
    

    int nLoc_Elem   = m_Elem2Vert.size();
    vloc            = m_LocalVertsMap.size();
    floc            = cnf;


}













void PartObject::GenerateFace2ElementMap(std::map<int,std::vector<int> > Face2Elem_i,
                                            MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::map<int,std::vector<int> >::iterator itr;

    for(itr=m_Elem2Face.begin();itr!=m_Elem2Face.end();itr++)
    {
        int elid = itr->first;
        int nf   = m_Elem2Face[elid].size();
        int nv   = m_Elem2Vert[elid].size(); 
        int type = m_elem2type_on_rank[elid];

        for(int q = 0; q < m_Elem2Face[elid].size(); q++)
        {
            int fid = m_Elem2Face[elid][q];
            m_Face2Elem[fid].push_back(elid);
        }
    }

    std::map<int,std::vector<int> >::iterator itmiv;
    int tel = 0;

    ParallelState* ife_pstate = new ParallelState(Nf_glob,comm);    
    std::map<int,std::vector<int> > rank2req_Faces;
    std::vector<int> new_offsets(size,0);
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }


    for(itmiv=m_Face2Elem.begin();itmiv!=m_Face2Elem.end();itmiv++)
    {
        int faceid   = itmiv->first;
        int numel    = itmiv->second.size();

        if(numel == 1)
        {
            m_Face2Elem[faceid].clear();

            int r = FindRank(new_offsets.data(),size,faceid);

            if(r != rank)
            {
                rank2req_Faces[r].push_back(faceid);
            }
            else
            {
                if(Face2Elem_i.find(faceid)!=Face2Elem_i.end())
                {
                    std::vector<int> newrow(2,0);
                    newrow[0] = Face2Elem_i[faceid][0];
                    newrow[1] = Face2Elem_i[faceid][1];
                    m_Face2Elem[faceid] = newrow;
                }
            }
        }
    }


    ScheduleObj* ife_schedule = DoScheduling(rank2req_Faces,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_F_IDs_per_rank;

    for(int q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Faces.begin(); it != rank2req_Faces.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);

                i++;
            }
        }
        else if (ife_schedule->SendFromRank2Rank[q].find( rank ) != ife_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            
            reqstd_F_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IFE_Face_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IFE_recv   = 0;
    int eIFE_id             = 0;


    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nife;
    std::map<int,std::vector<int> > recv_back_face_ids;
    std::map<int,std::vector<int> > recv_back_ife;
    std::map<int,std::vector<int> > recv_back_face_Ne;
    int n_recv_back;

    for(int q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send = it->second.size();
                std::vector<int> ife_send;
                std::vector<int> fncol(it->second.size(),0);
                for(int u=0;u<it->second.size();u++)
                {
                    int ncol = 2;
                    ife_send.push_back(Face2Elem_i[it->second[u]][0]);
                    ife_send.push_back(Face2Elem_i[it->second[u]][1]);
                }

                int nfe_send = ife_send.size();

                int dest = it->first;
                MPI_Send(&nfe_send, 1, MPI_INT, dest, 223*6666+1000*dest, comm);
                MPI_Send(&ife_send[0], nfe_send, MPI_INT, dest, 9876*6666+dest*8888, comm);

            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            int nfe_recv;
            MPI_Recv(&nfe_recv, 1, MPI_INT, q, 223*6666+1000*rank, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_back_ife_arr(nfe_recv);
            MPI_Recv(&recv_back_ife_arr.data()[0], nfe_recv, MPI_INT, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_ife[q]        = recv_back_ife_arr;

         }
    }



    std::map<int,std::vector<int> >::iterator iter;
    for(iter=rank2req_Faces.begin();iter!=rank2req_Faces.end();iter++)
    {
        int recvdrank = iter->first;

        int L = iter->second.size();
        int offset = 0;
        int ncol = 0;
        for(int s=0;s<L;s++)
        {
            int face_id  = rank2req_Faces[recvdrank][s];
            ncol     = 2;

            std::vector<int> ife_loc_row(2,0);
            std::vector<int> ifref_loc_row(1,0);
            ife_loc_row[0]   = recv_back_ife[iter->first][offset+0];
            ife_loc_row[1]   = recv_back_ife[iter->first][offset+1];

            m_Face2Elem[face_id] = ife_loc_row;

            offset = offset + ncol;
        }
    }
}





void PartObject::GenerateFace2VertexMap(MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<std::vector<int> > tetra_faces   = getTetraFaceMap(); 
    std::vector<std::vector<int> > prism_faces   = getPrismFaceMap(); 
    std::vector<std::vector<int> > pyramid_faces = getPyramidFaceMap(); 
    // std::vector<std::vector<int> > hex_faces  = getHexFaceMap(); 

    //std::cout << "tetra_facestetra_faces " << tetra_faces.size() << std::endl;

    //FaceOnRankPointer m_FaceOnRankPointer;
    //start = clock();
    std::vector<std::vector<int> > map_faces;
    std::map<int,std::vector<int> >::iterator itr;
    for(itr=m_Elem2Face.begin();itr!=m_Elem2Face.end();itr++)
    {
        int elid = itr->first;
        int nf   = m_Elem2Face[elid].size();
        int nv   = m_Elem2Vert[elid].size(); 
        int type = m_elem2type_on_rank[elid];
        // if(type != 2 && type != 5 && type != 6)
        // {
        //     std::cout << "Different type " << type << std::endl;
        // }
        switch (type) {
        case 2:
            map_faces = tetra_faces;
            break;
        case 5:
            map_faces = pyramid_faces;
            break;
        case 6:
            map_faces = prism_faces;
            break;
        // case 4:
        //     map_faces = hex_faces;
        //     break;
        }
        //std::cout << map_faces.size() << " " << type << std::endl;
        for(int u=0;u<map_faces.size();u++)
        {
            int nv  = map_faces[u].size();
            int fid = m_Elem2Face[elid][u];

            std::vector<int> face(nv,0);
            //std::cout << "nv " << nv << std::endl;
            //std::cout << "m_Elem2Vert[elid] " << m_Elem2Vert[elid].size() << std::endl;
            for(int w=0;w<nv;w++)
            {
                face[w] = m_Elem2Vert[elid][map_faces[u][w]];
            }
            // if(rank == 0)
            // {
            //     std::cout << std::endl;
            // }
            
            m_Face2Vert[fid] = face;
        }
        // if(rank == 0)
        // {
        //     std::cout << std::endl;
        // }
    }
}








void PartObject::getFace2VertexPerPartitionMap(std::map<int,std::vector<int> > ifn_read, 
                                               std::map<int,std::vector<int> > if_Nv_read, 
                                               int Nf_glob, 
                                               MPI_Comm comm)
{
    
    ParallelState* ife_pstate = new ParallelState(Nf_glob,comm);    
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int ife_o = ife_pstate->getOffset(rank);
    int face_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Faces;
    int* new_offsets = new int[size];
    
    
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
    
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    int Nel = nElemGlobalPart;
    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;

    for(itefmap=m_Elem2Face.begin();itefmap!=m_Elem2Face.end();itefmap++)
    {
        for(int q=0;q<itefmap->second.size();q++)
        {
            int face_req = itefmap->second[q];
            
            r = FindRank(new_offsets,size,face_req);
            
            if(r != rank)
            {
                rank2req_Faces[r].push_back(face_req);
            }
            else
            {
                if(m_Face2Vert.find(face_req)==m_Face2Vert.end())
                {
                    int Nv = if_Nv_read[face_req][0];
                    std::vector<int> row(Nv,0);
                    for(int i=0;i<Nv;i++)
                    {
                        row[i] = ifn_read[face_req][i];
                    }
                    m_Face2Vert[face_req] = row;
                }
            }
        }
    }
    
    int own = m_Face2Vert.size();
    
    ScheduleObj* ife_schedule = DoScheduling(rank2req_Faces,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_F_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Faces.begin(); it != rank2req_Faces.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);

                i++;
            }
        }
        else if (ife_schedule->SendFromRank2Rank[q].find( rank ) != ife_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            
            reqstd_F_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IFE_Face_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IFE_recv   = 0;
    int eIFE_id             = 0;

    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nife;
    std::map<int,std::vector<int> > recv_back_face_ids;
    std::map<int,std::vector<int> > recv_back_ife;
    std::map<int,std::vector<int> > recv_back_face_Ne;
    int n_recv_back;

    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send = it->second.size();
                std::vector<int> ife_send;
                std::vector<int> fncol(it->second.size(),0);
                for(int u=0;u<it->second.size();u++)
                {
                    int ncol = if_Nv_read[it->second[u]][0];

                    for(int s=0;s<ncol;s++)
                    {
                        ife_send.push_back(ifn_read[it->second[u]][s]);
                    }

                    fncol[u] = ncol;
                }

                int nfe_send = ife_send.size();

                int dest = it->first;
                MPI_Send(&nf_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                MPI_Send(&nfe_send, 1, MPI_INT, dest, 223*6666+1000*dest, comm);

                MPI_Send(&fncol[0], nf_send, MPI_INT, dest, 12*6666+dest*8888, comm);
                MPI_Send(&ife_send[0], nfe_send, MPI_INT, dest, 9876*6666+dest*8888, comm);

            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            int nfe_recv;
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nfe_recv, 1, MPI_INT, q, 223*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            std::vector<int> recv_back_ife_arr(nfe_recv);
            std::vector<int> recv_back_ids_arr(n_recv_back);
            std::vector<int> fncol_rcv(n_recv_back);

            MPI_Recv(&fncol_rcv.data()[0], n_recv_back, MPI_INT, q, 12*6666+rank*8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_ife_arr.data()[0], nfe_recv, MPI_INT, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_face_Ne[q]    = fncol_rcv;
            recv_back_ife[q]        = recv_back_ife_arr;

         }
    }



    std::map<int,std::vector<int> >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=rank2req_Faces.begin();iter!=rank2req_Faces.end();iter++)
    {
        int recvdrank = iter->first;

        int L = iter->second.size();
        int offset = 0;
        int ncol = 0;
        for(int s=0;s<L;s++)
        {
            face_id  = rank2req_Faces[recvdrank][s];
            ncol     = recv_back_face_Ne[iter->first][s];

            std::vector<int> ifref_loc_row(ncol,0);

            for(int v=0;v<ncol;v++)
            {
                ifref_loc_row[v] = recv_back_ife[iter->first][offset+v];
            }

            m_Face2Vert[face_id]  = ifref_loc_row;

            offset = offset + ncol;
        }

        ntotal=ntotal+L;
    }
}


void PartObject::CommunicateAdjacencyInfo(MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::map<int,std::vector<int> >::iterator itm;
    int nreq;

    std::set<int> toquery_elids;
    for(itm=m_Elem2Elem.begin();itm!=m_Elem2Elem.end();itm++)
    {
        int elid = itm->first;

        for(int q=0;q<itm->second.size();q++)
        {
            int adjeid = itm->second[q];

            if(m_ElemSet.find(adjeid)==m_ElemSet.end() && adjeid<Ne_glob)
            {
                if(toquery_elids.find(adjeid)==toquery_elids.end())
                {
                    toquery_elids.insert(adjeid);
                }
            }
            else
            {
                m_Elem2Rank[adjeid] = rank;
            }
        }
    }

    std::vector<int> toquery_elids_vec(toquery_elids.size(),0);

    std::copy(toquery_elids.begin(), 
                toquery_elids.end(), 
                toquery_elids_vec.begin());

    if(rank != 0)
    {
        nreq = toquery_elids_vec.size();
        //std::cout << "nreq " << nreq << std::endl;
        MPI_Send(&nreq, 1, MPI_INT, 0, rank*100000, comm);
        MPI_Send(&toquery_elids_vec.data()[0], nreq, MPI_INT, 0, rank*1000000, comm);
        
        std::vector<int> pid_2recvback(nreq,0);
        MPI_Recv(&pid_2recvback[0], nreq, MPI_INT, 0, rank*20000, comm, MPI_STATUS_IGNORE);

        for(int i=0;i<toquery_elids_vec.size();i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id    = pid_2recvback[i];

            m_Elem2Rank[el_id] = p_id;
        }

        //std::cout << "adjacent2pid " << adjacent2pid.size() << " " << rank << std::endl;
    }
    else if(rank == 0)
    {
        std::vector<int> nrecv_toquery_elids(size-1,0);
        int accul = 0;
        for(int i=1;i<size;i++)
        {
            int dest = i*100000;
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, i, dest, comm, MPI_STATUS_IGNORE);
            nrecv_toquery_elids[i-1] = n_reqstd_ids;
            accul = accul + n_reqstd_ids;
        }

        for(int i=1;i<size;i++)
        {
            int n_reqstd_ids = nrecv_toquery_elids[i-1];
            std::vector<int> QueryOnRoot(n_reqstd_ids,0);
            std::vector<int> pid_2sendback(n_reqstd_ids,0);
            MPI_Recv(&QueryOnRoot[0], n_reqstd_ids, MPI_INT, i, i*1000000, comm, MPI_STATUS_IGNORE);
            for(int j=0;j<n_reqstd_ids;j++)
            {
                
                if(m_partGlobalRoot.find(QueryOnRoot[j])!=m_partGlobalRoot.end())
                {
                    pid_2sendback[j] = m_partGlobalRoot[QueryOnRoot[j]];
                }
                else
                {                    
                    pid_2sendback[j] = -1;
                }
            }

            MPI_Send(&pid_2sendback.data()[0], n_reqstd_ids, MPI_INT, i, i*20000, comm);
        }


        for(int i=0;i<toquery_elids_vec.size();i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id = -1;
            if(m_partGlobalRoot.find(el_id)!=m_partGlobalRoot.end())
            {
                p_id = m_partGlobalRoot[el_id];
            }

            m_Elem2Rank[el_id] = p_id;
        }
        
    }
}



std::map<int,std::vector<int> >  PartObject::CommunicateAdjacencyInfoLocalPartition(MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::map<int,std::vector<int> > Rank2Elem;
    std::map<int,std::vector<int> >::iterator itm;
    int nreq;

    std::set<int> toquery_elids;

    for(itm=m_Elem2Elem.begin();itm!=m_Elem2Elem.end();itm++)
    {
        int elid = itm->first;

        for(int q=0;q<itm->second.size();q++)
        {
            int adjeid = itm->second[q];

            if(m_ElemSet.find(adjeid)==m_ElemSet.end() && adjeid<Ne_glob)
            {
                if(toquery_elids.find(adjeid)==toquery_elids.end())
                {
                    toquery_elids.insert(adjeid);
                }
            }
            else
            {
                m_Elem2Rank[adjeid] = rank;
            }
        }
    }

    
    std::vector<int> toquery_elids_vec(toquery_elids.size(),0);
    std::copy(toquery_elids.begin(), 
                toquery_elids.end(), 
                toquery_elids_vec.begin());


    if(rank != 0)
    {
        nreq = toquery_elids_vec.size();
        MPI_Send(&nreq, 1, MPI_INT, 0, rank*100000, comm);
        MPI_Send(&toquery_elids_vec.data()[0], nreq, MPI_INT, 0, rank*1000000, comm);
        
        int nrec_b;
        MPI_Recv(&nrec_b, 1, MPI_INT, 0, rank*250000, comm, MPI_STATUS_IGNORE);
        std::vector<int> pid_2recvback(nrec_b,0);
        std::vector<int> el_2recvback(nrec_b,0);
        MPI_Recv(&el_2recvback[0], nrec_b,  MPI_INT, 0, rank*225000, comm, MPI_STATUS_IGNORE);
        MPI_Recv(&pid_2recvback[0], nrec_b, MPI_INT, 0, rank*20000, comm, MPI_STATUS_IGNORE);

        for(int i=0;i<nrec_b;i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id    = pid_2recvback[i];

            m_Elem2Rank[el_id] = p_id;

            if(p_id != -1)
            {
                Rank2Elem[p_id].push_back(el_id);
            }
        }
    }
    else if(rank == 0)
    {
        std::vector<int> nrecv_toquery_elids(size-1,0);
        int accul = 0;
        for(int i=1;i<size;i++)
        {
            int dest = i*100000;
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, i, dest, comm, MPI_STATUS_IGNORE);
            nrecv_toquery_elids[i-1] = n_reqstd_ids;
            accul = accul + n_reqstd_ids;
        }

        for(int i=1;i<size;i++)
        {
            int n_reqstd_ids = nrecv_toquery_elids[i-1];
            std::vector<int> QueryOnRoot(n_reqstd_ids,0);
            std::vector<int> pid_2sendback;
            std::vector<int> el_2sendback;
            MPI_Recv(&QueryOnRoot[0], n_reqstd_ids, MPI_INT, i, i*1000000, comm, MPI_STATUS_IGNORE);
            for(int j=0;j<n_reqstd_ids;j++)
            {
                
                if(m_partGlobalRoot.find(QueryOnRoot[j])!=m_partGlobalRoot.end())
                {
                    pid_2sendback.push_back(m_partGlobalRoot[QueryOnRoot[j]]);
                    el_2sendback.push_back(QueryOnRoot[j]);
                }
                else
                {                    
                    pid_2sendback.push_back(-1);
                    el_2sendback.push_back(QueryOnRoot[j]);
                }
            }
            int send_b = pid_2sendback.size();
            MPI_Send(&send_b, 1, MPI_INT, i, i*250000, comm);
            MPI_Send(&el_2sendback.data()[0], el_2sendback.size(), MPI_INT, i, i*225000, comm);
            MPI_Send(&pid_2sendback.data()[0], pid_2sendback.size(), MPI_INT, i, i*20000, comm);
        }


        for(int i=0;i<toquery_elids_vec.size();i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id    = -1;
            
            if(m_partGlobalRoot.find(el_id)!=m_partGlobalRoot.end())
            {
                p_id = m_partGlobalRoot[el_id];
                Rank2Elem[p_id].push_back(el_id);
            }
            m_Elem2Rank[el_id] = p_id;
      
        }
    }
    return Rank2Elem;  
}






std::map<int,std::vector<int> >  PartObject::CommunicateAdjacencyInfoExtendedPartition(MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::map<int,std::vector<int> > Rank2Elem;
    std::map<int,std::vector<int> >::iterator itm;
    int nreq;

    std::set<int> toquery_elids;
    int here = 0;
    for(itm=m_Elem2Elem.begin();itm!=m_Elem2Elem.end();itm++)
    {
        int elid = itm->first;

        for(int q=0;q<itm->second.size();q++)
        {
            int adjeid = itm->second[q];

            if(m_Elem2Elem.find(adjeid)==m_Elem2Elem.end() && adjeid<Ne_glob)
            {
                if(toquery_elids.find(adjeid)==toquery_elids.end())
                {
                    toquery_elids.insert(adjeid);
                }
            }
            // else
            // {
            //     m_Elem2Rank[adjeid] = rank;
            // }
        }
    }

    
    std::vector<int> toquery_elids_vec(toquery_elids.size(),0);
    std::copy(toquery_elids.begin(), 
                toquery_elids.end(), 
                toquery_elids_vec.begin());


    if(rank != 0)
    {
        nreq = toquery_elids_vec.size();
        MPI_Send(&nreq, 1, MPI_INT, 0, rank*100000, comm);
        MPI_Send(&toquery_elids_vec.data()[0], nreq, MPI_INT, 0, rank*1000000, comm);
        
        int nrec_b;
        MPI_Recv(&nrec_b, 1, MPI_INT, 0, rank*250000, comm, MPI_STATUS_IGNORE);
        std::vector<int> pid_2recvback(nrec_b,0);
        std::vector<int> el_2recvback(nrec_b,0);
        MPI_Recv(&el_2recvback[0], nrec_b,  MPI_INT, 0, rank*225000, comm, MPI_STATUS_IGNORE);
        MPI_Recv(&pid_2recvback[0], nrec_b, MPI_INT, 0, rank*20000, comm, MPI_STATUS_IGNORE);

        for(int i=0;i<nrec_b;i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id    = pid_2recvback[i];

            if(p_id != -1)
            {
                Rank2Elem[p_id].push_back(el_id);
            }
            //Rank2Elem[p_id].push_back(el_id);
        }
    }
    else if(rank == 0)
    {
        std::vector<int> nrecv_toquery_elids(size-1,0);
        int accul = 0;
        for(int i=1;i<size;i++)
        {
            int dest = i*100000;
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, i, dest, comm, MPI_STATUS_IGNORE);
            nrecv_toquery_elids[i-1] = n_reqstd_ids;
            accul = accul + n_reqstd_ids;
        }

        for(int i=1;i<size;i++)
        {
            int n_reqstd_ids = nrecv_toquery_elids[i-1];
            std::vector<int> QueryOnRoot(n_reqstd_ids,0);
            std::vector<int> pid_2sendback;
            std::vector<int> el_2sendback;
            MPI_Recv(&QueryOnRoot[0], n_reqstd_ids, MPI_INT, i, i*1000000, comm, MPI_STATUS_IGNORE);
            for(int j=0;j<n_reqstd_ids;j++)
            {
                
                if(m_partGlobalRoot.find(QueryOnRoot[j])!=m_partGlobalRoot.end())
                {
                    pid_2sendback.push_back(m_partGlobalRoot[QueryOnRoot[j]]);
                    el_2sendback.push_back(QueryOnRoot[j]);
                }
                else
                {                    
                    pid_2sendback.push_back(-1);
                    el_2sendback.push_back(QueryOnRoot[j]);
                }
            }
            int send_b = pid_2sendback.size();
            MPI_Send(&send_b, 1, MPI_INT, i, i*250000, comm);
            MPI_Send(&el_2sendback.data()[0], el_2sendback.size(), MPI_INT, i, i*225000, comm);
            MPI_Send(&pid_2sendback.data()[0], pid_2sendback.size(), MPI_INT, i, i*20000, comm);
        }


        for(int i=0;i<toquery_elids_vec.size();i++)
        {
            int el_id   = toquery_elids_vec[i];
            int p_id    = -1;
            
            if(m_partGlobalRoot.find(el_id)!=m_partGlobalRoot.end())
            {
                p_id = m_partGlobalRoot[el_id];
                Rank2Elem[p_id].push_back(el_id);
            }      
        }
    }

    return Rank2Elem;  
}





void PartObject::GenerateTraceMap()
{
    // std::cout << "m_Elem2Rank " << m_Elem2Rank.size() << " " << m_Elem2Elem.size() << std::endl; 
    std::map<int,std::vector<int> >::iterator itm;
    int cnt = 0;
    for(itm=m_Elem2Elem.begin();itm!=m_Elem2Elem.end();itm++)
    {
        int elid = itm->first;
        int ne = itm->second.size();
        for(int q=0;q<ne;q++)
        {
            int elemid = itm->second[q];
            int faceid = m_Elem2Face[elid][q];
            std::vector<int> leftright(2,0);
            // Change this so that it bases this of the element type.
            if(m_Elem2Rank.find(elemid)!=m_Elem2Rank.end())
            {
                // std::cout << "m_Elem2Rank[elemid] = " << m_Elem2Rank[elemid] << " " << elemid << std::endl;
                if(m_Elem2Rank[elemid]==-1)
                {
                    m_TraceFacesOnRank.insert(faceid);
                    // std::cout << "m_Elem2Rank[elemid] = 1" << std::endl;
                    int nv = m_Face2Vert[faceid].size();
                    for(int n=0;n<nv;n++)
                    {
                        int vid = m_Face2Vert[faceid][n];

                        if(m_TraceVertsOnRank.find(vid)==m_TraceVertsOnRank.end())
                        {
                            m_TraceVertsOnRank.insert(vid);
                            cnt++;
                        }
                        m_TraceFace2Vert[faceid] = m_Face2Vert[faceid];
                        leftright[0] = elid;
                        leftright[1] = elemid;
                        m_TraceFace2Elem[faceid]  = leftright;
                    }
                }
            }
        }
    }  
}















std::map<int, std::vector<int> > PartObject::getAdjacentElementLayer(std::map<int,std::vector<double> > xcn, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    ParallelState* ife_pstate = new ParallelState(Nf_glob,comm);
    ParallelState* xcn_pstate = new ParallelState(Nv_glob,comm);
    int* new_V_offsets = new int[size];
    for(int i=0;i<size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
    }
    std::map<int, std::vector<int> > adjEl2Face;

    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int e_id;
    int r;
        int itel = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    int ff   = 0;
    std::vector<int> faceIDs_on_rank;
        
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_face;
    std::map<int, std::vector<int> > adj_el;

    std::map<int,std::vector<int> > Rank2Elem = CommunicateAdjacencyInfoExtendedPartition(comm);
    std::map<int,std::vector<int> >::iterator iter;


    for(iter = Rank2Elem.begin();iter != Rank2Elem.end(); iter++)
    {
        int ra = iter->first;
        
        if(m_Rank2ReqElem.find(ra)==m_Rank2ReqElem.end())
        {
            std::set<int> el_list;
            for(int q=0;q<iter->second.size();q++)
            {
                el_list.insert(iter->second[q]);
            }
            m_Rank2ReqElem[ra] = el_list;
        }
        if(m_Rank2ReqElem.find(ra)!=m_Rank2ReqElem.end())
        {
            for(int q=0;q<iter->second.size();q++)
            {
                m_Rank2ReqElem[ra].insert(iter->second[q]);
            }
        }
    }

    //std::cout << "Funished " << adj_elements.size() << std::endl;
    // adj_elements = adjid2rank;
    // std::cout << "rank " << rank << " adjid2rank " << adjid2rank.size() << std::endl;

    // if(rank == 1)
    // {
    //     std::map<int,int>::iterator itc;
    //     std::cout << "Should ivert right " << std::endl;
    //     for(itc=adjid2rank.begin();itc!=adjid2rank.end();itc++)
    //     {
    //         std::cout << itc->first << " " << itc->second << std::endl;
    //     }
    // }

    ScheduleObj* adj_schedule = DoScheduling(Rank2Elem, comm);
    std::map<int,std::vector<int> > reqstd_adj_ids_per_rank;
    std::map<int,std::vector<int> >::iterator it;
    int n_reqstd_adj_ids;


    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = Rank2Elem.begin(); it != Rank2Elem.end(); it++)
            {
                int n_req_adj_el           = it->second.size();
                int dest                   = it->first;
                //std::cout << dest << " " << adj_elements.size() << std::endl;
                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req_adj_el, 1, MPI_INT, dest, 9876000+10*dest, comm);
                //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second.data()[0], n_req_adj_el, MPI_INT, dest, 9876000*2+dest*2, comm);
               i++;
            }
        }
        else if (adj_schedule->SendFromRank2Rank[q].find( rank ) != adj_schedule->SendFromRank2Rank[q].end())
        {
            MPI_Recv(&n_reqstd_adj_ids, 1, MPI_INT, q, 9876000+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_adj_ids(n_reqstd_adj_ids);
            MPI_Recv(&recv_reqstd_adj_ids[0], n_reqstd_adj_ids, MPI_INT, q, 9876000*2+rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_adj_ids_per_rank[q] = recv_reqstd_adj_ids;
        }
    }


    //std::cout << reqstd_adj_ids_per_rank.size() << " " << rank << std::endl;

   
    std::map<int,std::vector<int> >::iterator itv;
    std::map<int,std::vector<int> > send_adj_verts_IDs;
    std::map<int,std::vector<int> > send_adj_faces_IDs;
    std::map<int,std::vector<int> > send_adj_element_IDs;

    std::map<int,std::vector<int> > send_adj_NvertsPel;
    std::map<int,std::vector<int> > send_adj_NfacesPel;
    std::map<int,std::vector<int> > send_adj_NdatasPel;
    std::map<int,std::vector<double> > send_adj_data;
    int TotNelem_adj_recv = 0;
    std::vector<int> TotAdj_El_IDs;
    int adj_id;

    int lelem = 0;

    for(itv=reqstd_adj_ids_per_rank.begin();itv!=reqstd_adj_ids_per_rank.end();itv++)
    {
        int dest = itv->first;
        for(int j=0;j<itv->second.size();j++)
        {
            adj_id = itv->second[j];
            TotAdj_El_IDs.push_back(adj_id);

            int nvPerEl = m_Elem2Vert[adj_id].size();
            int nfPerEl = m_Elem2Face[adj_id].size();
            int ndata   = m_Elem2Data[adj_id].size();
            //Array<double>* Var  = LocAndAdjElemVaria[adj_id];
            
            send_adj_NvertsPel[dest].push_back(nvPerEl);
            send_adj_NfacesPel[dest].push_back(nfPerEl);
            send_adj_NdatasPel[dest].push_back(ndata);

            for(int k=0;k<ndata;k++)
            {
                double dataV = m_Elem2Data[adj_id][k];
                send_adj_data[dest].push_back(dataV);
            }

            //send_adj_VarPel[dest].push_back();
            
            for(int k=0;k<nvPerEl;k++)
            {
                v_id = m_Elem2Vert[adj_id][k];
                send_adj_verts_IDs[dest].push_back(v_id);
            }
            for(int k=0;k<nfPerEl;k++)
            {
                f_id = m_Elem2Face[adj_id][k];
                send_adj_faces_IDs[dest].push_back(f_id);
                e_id = m_Elem2Elem[adj_id][k];
                send_adj_element_IDs[dest].push_back(e_id);
            }
        }

        TotNelem_adj_recv = TotNelem_adj_recv + itv->second.size();
    }
    
    int offset_adj_xcn = 0;
    int nloc_adj_xcn   = 0;
    std::map<int,int  > recv_adj_back_Nverts;
    std::map<int,std::vector<int> > recv_adj_back_verts_ids;
    std::map<int,int  > recv_adj_back_Nfaces;
    std::map<int,std::vector<int> > recv_adj_back_faces_ids;
    std::map<int,std::vector<int> > recv_adj_back_element_ids;

    //std::map<int,int > recv_adj_back_Nrhos;
    std::map<int,std::vector<int>  > recv_adj_NvPel;
    std::map<int,std::vector<int>  > recv_adj_NfPel;
    std::map<int,std::vector<int>  > recv_adj_NdataPel;
    std::map<int,std::vector<double>  >recv_adj_NVarPel;
    int n_adj_vert_recv_back;
    int n_adj_face_recv_back;
     
    // This sends the right vertices of the requested elements to correct processor.
    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = send_adj_verts_IDs.begin(); it != send_adj_verts_IDs.end(); it++)
            {
                int nv_adj_send       = it->second.size();
                int dest = it->first;
                MPI_Send(&nv_adj_send, 1, MPI_INT, dest, 98760000+1000*dest, comm);
                MPI_Send(&it->second.data()[0], it->second.size(), MPI_INT, dest, 19999*9876+dest*8888,comm);
                
                int NnvPel    = send_adj_NvertsPel[it->first].size();
                int NnfPel    = send_adj_NfacesPel[it->first].size();
                int NndataPel = send_adj_NdatasPel[it->first].size();
                int NVarPel   = send_adj_data[it->first].size();
                MPI_Send(&NnvPel,  1, MPI_INT, dest, 98764444+5000*dest, comm);
                MPI_Send(&NnfPel,  1, MPI_INT, dest, 98764444-5000*dest, comm);
                MPI_Send(&NndataPel,  1, MPI_INT, dest, 98765555-5000*dest, comm);

                MPI_Send(&NVarPel, 1, MPI_INT, dest, 98766666-5000*dest, comm);

                MPI_Send(&send_adj_NvertsPel[it->first][0], NnvPel, MPI_INT, dest, 98364444+15000*dest, comm);
                MPI_Send(&send_adj_NfacesPel[it->first][0], NnfPel, MPI_INT, dest, 98364444-15000*dest, comm);
                MPI_Send(&send_adj_NdatasPel[it->first][0], NndataPel, MPI_INT, dest, 98364444-15000*dest, comm);

                MPI_Send(&send_adj_data[it->first][0], NVarPel, MPI_DOUBLE, dest, 98366666-15000*dest, comm);

                int nf_adj_send = send_adj_faces_IDs[it->first].size();
                MPI_Send(&nf_adj_send, 1, MPI_INT, dest, 3333*9876+dest*8888,comm);
                MPI_Send(&send_adj_faces_IDs[it->first][0], nf_adj_send, MPI_INT, dest, 2222*9876+dest*8888,comm);
                MPI_Send(&send_adj_element_IDs[it->first][0], nf_adj_send, MPI_INT, dest, 5555*9876+dest*8888,comm);

                //int n_adj_rhos = send_adj_rhos[it->first].size();
                //MPI_Send(&n_adj_rhos, 1, MPI_INT, dest, 4444*9876+dest*8888,comm);
                //MPI_Send(&send_adj_rhos[it->first][0], n_adj_rhos, MPI_DOUBLE, dest, 5555*9876+dest*8888,comm);


            }
        }
        if(adj_schedule->RecvRankFromRank[q].find( rank ) != adj_schedule->RecvRankFromRank[q].end())
        {
            MPI_Recv(&n_adj_vert_recv_back, 1, MPI_INT, q, 98760000+1000*rank, comm, MPI_STATUS_IGNORE);
            //int* recv_adj_back_arr_ids = new int[n_adj_vert_recv_back];
            std::vector<int> recv_adj_back_arr_ids(n_adj_vert_recv_back);
            MPI_Recv(&recv_adj_back_arr_ids.data()[0], n_adj_vert_recv_back, MPI_INT, q, 19999*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            
            int NnvPel_recv_back,NnfPel_recv_back,NndataPel_recv_back,NVarPel_recv_back;
            MPI_Recv(&NnvPel_recv_back, 1, MPI_INT, q, 98764444+5000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&NnfPel_recv_back, 1, MPI_INT, q, 98764444-5000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&NndataPel_recv_back, 1, MPI_INT, q, 98765555-5000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&NVarPel_recv_back, 1, MPI_INT, q, 98766666-5000*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> Nnv_RB(NnvPel_recv_back);
            std::vector<int> Nnf_RB(NnfPel_recv_back);
            std::vector<int> Nndata_RB(NndataPel_recv_back);
            std::vector<double> NnVar_RB(NVarPel_recv_back);

            MPI_Recv(&Nnv_RB[0], NnvPel_recv_back, MPI_INT, q, 98364444+15000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&Nnf_RB[0], NnfPel_recv_back, MPI_INT, q, 98364444-15000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&Nndata_RB[0], NndataPel_recv_back, MPI_INT, q, 98364444-15000*rank, comm, MPI_STATUS_IGNORE);

            MPI_Recv(&NnVar_RB[0], NVarPel_recv_back, MPI_DOUBLE, q, 98366666-15000*rank, comm, MPI_STATUS_IGNORE);

            int n_adj_face_recv_back;
            MPI_Recv(&n_adj_face_recv_back, 1, MPI_INT, q, 3333*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            //int* recv_adj_back_arr_face_ids = new int[n_adj_face_recv_back];
            std::vector<int> recv_adj_back_arr_face_ids(n_adj_face_recv_back);
            std::vector<int> recv_adj_back_arr_element_ids(n_adj_face_recv_back);

            MPI_Recv(&recv_adj_back_arr_face_ids.data()[0], n_adj_face_recv_back, MPI_INT, q, 2222*9876+rank*8888, comm,   MPI_STATUS_IGNORE);
            MPI_Recv(&recv_adj_back_arr_element_ids.data()[0], n_adj_face_recv_back, MPI_INT, q, 5555*9876+rank*8888, comm,   MPI_STATUS_IGNORE);

            //int n_adj_rho_recv_back;
            //MPI_Recv(&n_adj_rho_recv_back, 1, MPI_INT, q, 4444*9876+rank*8888, comm, MPI_STATUS_IGNORE);
            //double* recv_adj_back_arr_rho = new double[n_adj_rho_recv_back];
            //MPI_Recv(&recv_adj_back_arr_rho[0], n_adj_rho_recv_back, MPI_DOUBLE, q, 5555*9876+rank*8888, comm,   MPI_STATUS_IGNORE);


            recv_adj_back_Nverts[q]     = n_adj_vert_recv_back;
            recv_adj_back_verts_ids[q]  = recv_adj_back_arr_ids;

            recv_adj_back_Nfaces[q]     = n_adj_face_recv_back;
            recv_adj_back_faces_ids[q]  = recv_adj_back_arr_face_ids;
            recv_adj_back_element_ids[q]  = recv_adj_back_arr_element_ids;

            recv_adj_NvPel[q] = Nnv_RB;
            recv_adj_NfPel[q] = Nnf_RB;
            recv_adj_NdataPel[q] = Nndata_RB;
            recv_adj_NVarPel[q] = NnVar_RB;
        }
    }
     
    int TotNvert_adj_recv = 0;
    int TotNface_adj_recv = 0;
    int TotNrho_adj_recv  = 0;
    
    
    std::map<int,int >::iterator itm;
    std::vector<int> adj_verts;
    for(itm=recv_adj_back_Nverts.begin();itm!=recv_adj_back_Nverts.end();itm++)
    {
        TotNvert_adj_recv = TotNvert_adj_recv+itm->second;
        for(int i=0;i<itm->second;i++)
        {
            adj_verts.push_back(recv_adj_back_verts_ids[itm->first][i]);
        }
    }

    std::vector<int> adj_faces;
    std::vector<int> adj_element;
    for(itm=recv_adj_back_Nfaces.begin();itm!=recv_adj_back_Nfaces.end();itm++)
    {
        TotNface_adj_recv = TotNface_adj_recv+itm->second;
        for(int i=0;i<itm->second;i++)
        {
            adj_faces.push_back(recv_adj_back_faces_ids[itm->first][i]);
            adj_element.push_back(recv_adj_back_element_ids[itm->first][i]);
        }
    }

    std::vector<int> Rank2Elem_vec;
    std::vector<int> NvPEl_rb;
    std::vector<int> NfPEl_rb;
    std::vector<int> NdataPEl_rb;
    std::vector<double> NVarPEl_rb;
    std::map<int,std::vector<int> >::iterator itm_el;
    int offvvv = 0;
   
    for(itm_el=Rank2Elem.begin();itm_el!=Rank2Elem.end();itm_el++)
    {
        
        for(int i=0;i<itm_el->second.size();i++)
        {
            Rank2Elem_vec.push_back(Rank2Elem[itm_el->first][i]);
            int Nv    = recv_adj_NvPel[itm_el->first][i];
            int Nf    = recv_adj_NfPel[itm_el->first][i];
            int ndata = recv_adj_NdataPel[itm_el->first][i];
            for(int q=0;q<ndata;q++)
            {
                NVarPEl_rb.push_back(recv_adj_NVarPel[itm_el->first][i*ndata+q]);
            }

            NvPEl_rb.push_back(Nv);
            NfPEl_rb.push_back(Nf);
            NdataPEl_rb.push_back(ndata);
            
            offvvv=offvvv+Nv;

        }
    }
    
    for(int i=0;i<adj_verts.size();i++)
    {
        int v_id_n = adj_verts[i];
        r = FindRank(new_V_offsets,size,v_id_n);

        if(m_vertIDs_on_rank.find( v_id_n ) == m_vertIDs_on_rank.end()) // add the required unique vertex for current rank.
        {
            m_vertIDs_on_rank.insert(v_id_n);
            //unique_verts_on_rank_vec.push_back(v_id_n);
            //part_v.push_back(r);

            if (r!=rank)// check whether this vertex is already present on current rank. if vertex is present on other rank, add it to vertIDs_on_rank map.
            {
                rank2req_vert[r].push_back(v_id_n); // add vertex to rank2req_vert map.
            }
            else
            {
                vertIDs_on_rank.push_back(v_id_n); // add the vertex to list that is already available on rank.
                vloc_tmp++;
            }
        }
    }
    
    for(int i=0;i<adj_faces.size();i++)
    {
        int f_id_n = adj_faces[i];

        if(m_faceIDs_on_rank.find( f_id_n ) == m_faceIDs_on_rank.end()) // add the required unique vertex for current rank.
        {
            m_faceIDs_on_rank.insert(f_id_n);
            faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
            floc_tmp++;

            if (r!=rank)// check whether this vertex is already present on current rank. if vertex is present on other rank, add it to vertIDs_on_rank map.
            {
                rank2req_face[r].push_back(f_id_n); // add vertex to rank2req_vert map.
            }
            else
            {
                faceIDs_on_rank.push_back(f_id_n); // add the vertex to list that is already available on rank.
                floc_tmp++;
            }
        }
        
    }
    
    

    // ==========================================================================================
       // ==========================================================================================
       // ==========================================================================================

       // At this point we have all the elements that are required on current rank and the vertex ids as well
       // However we are still missing the vertex coordinate data which is spread out equally over the available procs.
       // This rank2req_vert map essentially holds this information by mapping the rank_id from which we need to request a list/vector of vertex ids (hence the name "rank2req_vert" name.

       // At this point the perspective changes. When we were figuring out the layout of the elements, we knew the partition ID for each element on the current rank. This means that from the current rank, we needed to send a certain element to another rank since it is more logical to reside there. For the vertices this changes since we just figured out which vertices are required on the current rank. The logic here is first to send for each the current rank a list/vector<int> of vertex IDs that is requested from another rank. The other rank assembles the list of the required coordinates and sends it back.

       // ==========================================================================================
       // ==========================================================================================
       // ==========================================================================================

       int n_reqstd_ids;
       int n_req_recv_v2;

       // This thing needs to revised because for the verts it doesnt work.
       // The current rank does not have the verts_to_send_rank. Instead it has an request list.

       ScheduleObj* part_schedule = DoScheduling(rank2req_vert,comm);

       std::map<int,std::vector<int> >  reqstd_ids_per_rank;

       for(q=0;q<size;q++)
       {
           if(rank==q)
           {
               int i=0;
               for (it = rank2req_vert.begin(); it != rank2req_vert.end(); it++)
               {
                   int n_req           = it->second.size();
                   int dest            = it->first;

                   //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                   MPI_Send(&n_req, 1, MPI_INT, dest, 6547+10*dest, comm);
                   //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                   MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 6547*2+dest*2, comm);

                   i++;
               }
           }
           else if (part_schedule->SendFromRank2Rank[q].find( rank ) != part_schedule->SendFromRank2Rank[q].end())
           {
               MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6547+10*rank, comm, MPI_STATUS_IGNORE);
               //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

               std::vector<int> recv_reqstd_ids(n_reqstd_ids);
               MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 6547*2+rank*2, comm, MPI_STATUS_IGNORE);
               reqstd_ids_per_rank[q] = recv_reqstd_ids;
           }
       }

       int offset_xcn = 0;
       int nloc_xcn = 0;
       std::map<int,int > recv_back_Nverts;
       std::map<int,std::vector<double> > recv_back_verts;
       std::map<int,std::vector<int> > recv_back_verts_ids;
       int n_recv_back;
    
    
    

        
    
       for(q=0;q<size;q++)
       {
           if(rank == q)
           {
               for (it = reqstd_ids_per_rank.begin(); it != reqstd_ids_per_rank.end(); it++)
               {
                   int nv_send = it->second.size();
                   //double* vert_send = new double[nv_send*3];
                   std::vector<double> vert_send(nv_send*3);
                   offset_xcn        = xcn_pstate->getOffset(rank);
                   for(int u=0;u<it->second.size();u++)
                   {
                       // vert_send[u*3+0]=xcn->getVal(it->second[u]-offset_xcn,0);
                       // vert_send[u*3+1]=xcn->getVal(it->second[u]-offset_xcn,1);
                       // vert_send[u*3+2]=xcn->getVal(it->second[u]-offset_xcn,2);

                        vert_send[u*3+0]=xcn[it->second[u]][0];
                        vert_send[u*3+1]=xcn[it->second[u]][1];
                        vert_send[u*3+2]=xcn[it->second[u]][2];
                   }

                   int dest = it->first;
                   MPI_Send(&nv_send, 1, MPI_INT, dest, 6547+1000*dest, comm);
                   // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);

                   MPI_Send(&vert_send.data()[0], nv_send*3, MPI_DOUBLE, dest, 6547+dest*8888, comm);
                   MPI_Send(&it->second.data()[0], it->second.size(), MPI_INT, dest, 8888*6547+dest*8888,comm);

                   //delete[] vert_send;
               }
           }
           if(part_schedule->RecvRankFromRank[q].find( rank ) != part_schedule->RecvRankFromRank[q].end())
            {
               MPI_Recv(&n_recv_back, 1, MPI_INT, q, 6547+1000*rank, comm, MPI_STATUS_IGNORE);

               //double* recv_back_arr = new double[n_recv_back*3];
               std::vector<double> recv_back_arr(n_recv_back*3);
               //int* recv_back_arr_ids = new int[n_recv_back];
               std::vector<int> recv_back_arr_ids(n_recv_back);
               //MPI_Recv(&recv_back_vec[0], n_recv_back, MPI_DOUBLE, q, 9876+rank*888, comm, MPI_STATUS_IGNORE);
               MPI_Recv(&recv_back_arr.data()[0], n_recv_back*3, MPI_DOUBLE, q, 6547+rank*8888, comm, MPI_STATUS_IGNORE);
               MPI_Recv(&recv_back_arr_ids.data()[0], n_recv_back, MPI_INT, q, 8888*6547+rank*8888, comm, MPI_STATUS_IGNORE);

               recv_back_Nverts[q]     = n_recv_back;
               recv_back_verts[q]      = recv_back_arr;
               recv_back_verts_ids[q]  = recv_back_arr_ids;

               }
       }


       
       int vfor = 0;
       std::map<int,std::vector<double> >::iterator it_f;
       for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
       {
           int c  = 0;
           vfor=vfor+recv_back_Nverts[it_f->first];
       }


       int gvid=0;
       int lvid=vloc;

       int m=0;
       for(m=0;m<vloc_tmp;m++)
       {
           gvid = vertIDs_on_rank[m];

           if(m_GlovalV2LocalV.find(gvid)==m_GlovalV2LocalV.end())
           {
               
               std::vector<double> V(3);

               V[0] = xcn[gvid][0];
               V[1] = xcn[gvid][1];
               V[2] = xcn[gvid][2];
               
               m_LocalVertsMap[gvid] = V;

               lvid++;
           }
       }

       //int o = 3*vloc_tmp;
       int u = 0;
       for(it_f=recv_back_verts.begin();it_f!=recv_back_verts.end();it_f++)
       {
           int Nv = recv_back_Nverts[it_f->first];

           for(u=0;u<Nv;u++)
           {
               gvid = recv_back_verts_ids[it_f->first][u];
               
               if(m_GlovalV2LocalV.find(gvid)==m_GlovalV2LocalV.end())
               {
                   
                   std::vector<double> V(3);
                   V[0] = it_f->second[u*3+0];
                   V[1] = it_f->second[u*3+1];
                   V[2] = it_f->second[u*3+2];

                   //LocalVerts.push_back(V);
                   m_LocalVertsMap[gvid] = V;
                //    o_lvertex2gvertex[lvid]=gvid;
                //    m_GlovalV2LocalV[gvid]=lvid;
                   
                   lvid++;
               }
            }
       }

    int nLoc_Verts = m_LocalVertsMap.size();

    // ================================== Faces on Rank =========================================


    int cnv = 0;
    int cnf = 0;
    int idsave = 0;
    //double rho_v;
    int loc_v;
    int glob_f;
    int loc_f;
    int glob_v;
    int glob_e=0;
    int before = m_Elem2Elem.size();

    double varia_v = 0.0;
    std::set<int> LocAdjElemSet;
    std::vector<int> adjElLayer(3*itel);
    int offvv = 0;
    for(int m=0;m<Rank2Elem_vec.size();m++)
    {
        el_id         = Rank2Elem_vec[m];
        
        int Nv        = NvPEl_rb[m];
        int Nf        = NfPEl_rb[m];
        int ndatas    = NdataPEl_rb[m];

        std::vector<int> tmp_globv(Nv,0);
        std::vector<int> tmp_locv(Nv,0);
        std::vector<int> tmp_globf(Nf,0);
        std::vector<int> tmp_globe(Nf,0);

        std::vector<double> datarow(ndatas,0.0);
        for(int q=0;q<ndatas;q++)
        {
            datarow[q] = NVarPEl_rb[m*ndatas+q];
        }
    
        
        std::vector<double> Padj(Nv*3);
        for(int p=0;p<Nv;p++)
        {
            glob_v = adj_verts[offvv+p];
            loc_v  = m_GlovalV2LocalV[glob_v];

            Padj[p*3+0] = m_LocalVertsMap[glob_v][0];
            Padj[p*3+1] = m_LocalVertsMap[glob_v][1];
            Padj[p*3+2] = m_LocalVertsMap[glob_v][2];

            tmp_globv[p] = glob_v;
            tmp_locv[p] = loc_v;

            m_Elem2LocalVert[el_id].push_back(loc_v);//globElem2locVerts[el_id].push_back(loc_v);
            cnv++;
            
        }
        std::vector<double> Vadj  = ComputeCentroidCoord(Padj,Nv);
        m_Elem2Centroid[el_id] = Vadj;

        for(int p=0;p<Nf;p++)
        {
            glob_f = adj_faces[cnf];
            glob_e = adj_element[cnf];
            tmp_globf[p] = glob_f;
            tmp_globe[p] = glob_e;
            m_globFace2GlobalElements[glob_f].push_back(el_id);
            cnf++;
        }
    
        //=======================================
        m_Elem2Vert[el_id]      = tmp_globv;
        m_Elem2Face[el_id]      = tmp_globf;
        m_Elem2Elem[el_id]      = tmp_globe;
        m_Elem2Data[el_id]      = datarow;
        //=======================================
        adjEl2Face[el_id]       = tmp_globf;

        offvv = offvv+Nv;
        tmp_locv.clear();
    }

    delete[] new_V_offsets;

    reqstd_ids_per_rank.clear();
    recv_back_Nverts.clear();
    recv_back_verts.clear();
    recv_back_verts_ids.clear();
    
    NvPEl_rb.clear();
    NfPEl_rb.clear();
    
    rank2req_vert.clear();
    recv_back_Nverts.clear();
    recv_back_verts.clear();
    recv_back_verts_ids.clear();


    return adjEl2Face;

}








void PartObject::updateFace2EntityPerPartition(std::map<int,std::vector<int> > adjacent_ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > if_Nv_read,
                                                   MPI_Comm comm)
{
    
    ParallelState* ife_pstate = new ParallelState(Nf_glob,comm);    
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int ife_o = ife_pstate->getOffset(rank);
    int face_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Faces;
    int* new_offsets = new int[size];
    
    
    for(int i=0;i<size;i++)
    {
        new_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
    
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    //int Nel = part_global.size();

    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;
    
    for(itefmap=adjacent_ief.begin();itefmap!=adjacent_ief.end();itefmap++)
    {
        for(int q=0;q<itefmap->second.size();q++)
        {
            int face_req = itefmap->second[q];
            
            r = FindRank(new_offsets,size,face_req);
            
            if(r != rank && m_Face2Vert.find(face_req)==m_Face2Vert.end())
            {
                rank2req_Faces[r].push_back(face_req);
            }
            else
            {
                if(m_Face2Vert.find(face_req)==m_Face2Vert.end())
                {
                    int ncol = ife_read[face_req].size();
                    int ncol_d = if_Nv_read[face_req][0];
                    for(int j=0;j<ncol_d;j++)
                    {
                        int vrtid = ife_read[face_req][j];
                        m_Face2Vert[face_req].push_back(vrtid);
                        //std::cout << vrtid << " ";
                    }
                    //std::cout << std::endl;
                }
            }
        }
    }
    
    int own = m_Face2Vert.size();
    
    ScheduleObj* ife_schedule = DoScheduling(rank2req_Faces,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_F_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Faces.begin(); it != rank2req_Faces.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                //MPI_Send(&dest, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+10*dest, comm);
                //MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876+dest*2, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876*2*7654+dest*2, comm);

                i++;
            }
        }
        else if (ife_schedule->SendFromRank2Rank[q].find( rank ) != ife_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+10*rank, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&TotRecvVert_IDs[RecvAlloc_offset_map_v[q]], n_reqstd_ids, MPI_INT, q, 9876+rank*2, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*2, comm, MPI_STATUS_IGNORE);
            
            reqstd_F_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IFE_Face_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IFE_recv   = 0;
    int eIFE_id             = 0;


    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Nife;
    std::map<int,std::vector<int> > recv_back_face_ids;
    std::map<int,std::vector<int> > recv_back_ife;
    std::map<int,std::vector<int> > recv_back_face_Ne;
    int n_recv_back;

    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_F_IDs_per_rank.begin(); it != reqstd_F_IDs_per_rank.end(); it++)
            {
                int nf_send = it->second.size();
                std::vector<int> ife_send;
                std::vector<int> fncol(it->second.size(),0);
                for(int u=0;u<it->second.size();u++)
                {
                    //int ncol = ife_read[it->second[u]].size();
                    int ncol = if_Nv_read[it->second[u]][0];
                    for(int q=0;q<ncol;q++)
                    {
                        ife_send.push_back(ife_read[it->second[u]][q]);
                    }

                    fncol[u] = ncol;

                }

                int nfe_send = ife_send.size();


                int dest = it->first;
                MPI_Send(&nf_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                MPI_Send(&nfe_send, 1, MPI_INT, dest, 223*6666+1000*dest, comm);

                MPI_Send(&it->second.data()[0], nf_send, MPI_INT, dest, 9876*7777+dest*888, comm);
                MPI_Send(&fncol[0], nf_send, MPI_INT, dest, 12*6666+dest*8888, comm);
                MPI_Send(&ife_send[0], nfe_send, MPI_INT, dest, 9876*6666+dest*8888, comm);

            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            int nfe_recv;
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nfe_recv, 1, MPI_INT, q, 223*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            std::vector<int> recv_back_ife_arr(nfe_recv);
            std::vector<int> recv_back_ids_arr(n_recv_back);
            std::vector<int> fncol_rcv(n_recv_back);

            MPI_Recv(&recv_back_ids_arr.data()[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&fncol_rcv.data()[0], n_recv_back, MPI_INT, q, 12*6666+rank*8888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_ife_arr.data()[0], nfe_recv, MPI_INT, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nife[q]       = n_recv_back;
            recv_back_face_ids[q]   = recv_back_ids_arr;
            recv_back_face_Ne[q]    = fncol_rcv;
            recv_back_ife[q]        = recv_back_ife_arr;

         }
    }



    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=recv_back_Nife.begin();iter!=recv_back_Nife.end();iter++)
    {
        int L = iter->second;
        int offset = 0;
        for(int s=0;s<L;s++)
        {
            face_id  = recv_back_face_ids[iter->first][s];
            int ncol = recv_back_face_Ne[iter->first][s];
            
            std::vector<int> ife_row_loc(ncol,0);
            for(int r=0;r<ncol;r++)
            {
                ife_row_loc[r]          = recv_back_ife[iter->first][offset+r];
                //std::cout << "communicated  " << ife_row_loc[r] << " ";
                m_Face2Vert[face_id]    = ife_row_loc;
            }
            //std::cout << std::endl;
                
            
            offset = offset+ncol;
            
        }
        ntotal=ntotal+L;
    }
}






void PartObject::buildExtendedAdjacencyData(MPI_Comm comm, 
                                            std::map<int,std::vector<double> > ghosts)
{

    //int faceidtrace = 2353748;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<std::vector<double> > face;
    int tel = 0;
    std::set<int>::iterator its;

    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {
        int elid        = *its;
        int NvPEl       = m_Elem2Vert[elid].size();
        int nadj        = m_Elem2Elem[elid].size();

        std::set<int> AdjElem;

        std::vector<double> Vijk = m_Elem2Centroid[elid];
        for(int j=0;j<nadj;j++)
        {
            int adjid  = m_Elem2Elem[elid][j];
            int fadj   = m_Elem2Face[elid][j];

            if(m_Elem2Vert.find(adjid)!=m_Elem2Vert.end()
                && ghosts.find(adjid)==ghosts.end())
            {
                AdjElem.insert(adjid);
            }
            else if(ghosts.find(adjid)!=ghosts.end())
            {   
                if(m_GhostFaceVert.find(adjid)==m_GhostFaceVert.end())
                {    
                    std::vector<double> Vc = ComputeGhostCentroid(m_Face2Vert[fadj],m_LocalVertsMap,Vijk);
                    m_GhostFaceVert[adjid] = Vc;            
                }
                tel++;
            }
            if(m_Elem2Elem.find(adjid)!=m_Elem2Elem.end()
                && ghosts.find(adjid)==ghosts.end())
            {
                int n_adjid         = m_Elem2Elem[adjid].size();
                int NvPElnew        = m_Elem2Vert[adjid].size();

                for(int k=0;k<n_adjid;k++)
                {
                    int adjadj      = m_Elem2Elem[adjid][k];
                    int fadjadj     = m_Elem2Face[adjid][k];

                    if(m_Elem2Vert.find(adjadj)!=m_Elem2Vert.end() 
                    && ghosts.find(adjadj)==ghosts.end() && adjadj!=elid)
                    {
                        AdjElem.insert(adjadj);
                    }
                    else if(ghosts.find(adjadj)!=ghosts.end())
                    {
                        std::vector<double> Vadjadj = m_Elem2Centroid[adjid];

                        if(AdjElem.find(adjadj)==AdjElem.end())
                        {
                            AdjElem.insert(adjadj);
                        }
                        
                        if(m_GhostFaceVert.find(adjadj)==m_GhostFaceVert.end())
                        {
                            std::vector<double> Vc = ComputeGhostCentroid(m_Face2Vert[fadjadj],m_LocalVertsMap,Vadjadj);
                            m_GhostFaceVert[adjadj] = Vc;
                        }

                        tel++;

                    }
                    
                    if(m_Elem2Elem.find(adjadj)!=m_Elem2Elem.end()
                    && ghosts.find(adjadj)==ghosts.end())
                    {
                        int n_adjadj        = m_Elem2Elem[adjadj].size();
                        int NvPElnewnew     = m_Elem2Vert[adjadj].size();
                        for(int k=0;k<n_adjadj;k++)
                        {
                            int adjadjadj  = m_Elem2Elem[adjadj][k];
                            int fadjadjadj = m_Elem2Face[adjadj][k];

                            if(m_Elem2Vert.find(adjadjadj)!=m_Elem2Vert.end() 
                                && ghosts.find(adjadjadj)==ghosts.end() && adjadjadj!=elid)
                            {
                                AdjElem.insert(adjadjadj);
                            }
                            else if(ghosts.find(adjadjadj)!=ghosts.end())
                            {
                                std::vector<double> Vadjadjadj = m_Elem2Centroid[adjadj];
                                if(AdjElem.find(adjadjadj)==AdjElem.end())
                                {
                                    AdjElem.insert(adjadjadj);
                                }
                                

                                if(m_GhostFaceVert.find(adjadjadj)==m_GhostFaceVert.end())
                                {
                                    std::vector<double> Vc = ComputeGhostCentroid(m_Face2Vert[fadjadjadj],m_LocalVertsMap,Vadjadjadj);
                                    m_GhostFaceVert[adjadjadj] = Vc;
                                }

                                tel++;
                            }
                        }   
                    }
                }
            }
        }

        m_Elem2AdjElem[elid] = AdjElem;
        
    }
}







void PartObject::AddStateVecForAdjacentElements(std::map<int,std::vector<double> > &U, int nvar, MPI_Comm comm)
{
    
    int nanhere =0;
    int nothere = 0;
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);


    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > Rank2ReqElem_vec;
    std::map<int,std::vector<double> > U_loc;
    
    
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    
    std::vector<int> ee;

    std::map<int,std::set<int> >::iterator its;
    for(its=m_Rank2ReqElem.begin();its!=m_Rank2ReqElem.end();its++)
    {
        std::vector<int> elids(its->second.size(),0);
        std::copy(its->second.begin(), 
            its->second.end(), 
            elids.begin());

        Rank2ReqElem_vec[its->first] = elids;
    }
    
    ScheduleObj* iee_schedule = DoScheduling(Rank2ReqElem_vec,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_E_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = Rank2ReqElem_vec.begin(); it != Rank2ReqElem_vec.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876*7654+20*dest, comm);
                MPI_Send(&it->second.data()[0], n_req, MPI_INT, dest, 9876*2*7654+dest*40, comm);

                i++;
            }
        }
        else if (iee_schedule->SendFromRank2Rank[q].find( rank ) != iee_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876*7654+20*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2*7654+rank*40, comm, MPI_STATUS_IGNORE);
            
            reqstd_E_IDs_per_rank[q] = recv_reqstd_ids;
        }
    }
        
    std::map<int,std::vector<int> >::iterator ite;
    std::map<int,std::vector<int> > send_IEE_Elem_IDs;
    std::vector<int> TotIEE_El_IDs;

    int TotNelem_IEE_recv   = 0;
    int eIEE_id             = 0;
    
    int offset_xcn = 0;
    int nloc_xcn = 0;
    std::map<int,int > recv_back_Niee;
    std::map<int,std::vector<int> > recv_back_el_ids;
    std::map<int,std::vector<double> > recv_back_iee;
    int n_recv_back;

    for(q=0;q<size;q++)
    {
        if(rank == q)
        {
            for (it = reqstd_E_IDs_per_rank.begin(); it != reqstd_E_IDs_per_rank.end(); it++)
            {
                int ne_send             = it->second.size();
                std::vector<double> iee_send(ne_send*nvar);
                
                for(int u=0;u<ne_send;u++)
                {
                    
                    for(int s=0;s<nvar;s++)
                    {
                        if(U.find(it->second[u])==U.end())
                        {
                            nothere++;
                        }
                        
                        iee_send[u*nvar+s]=U[it->second[u]][s];
                        //std::cout << "it->second.size(); " << iee_send[u*nvar+s] << std::endl;
                        if(std::isnan(U[it->second[u]][s]))
                        {
                            nanhere++;
                        }
                        
                        
                    }
                }

                int dest = it->first;
                MPI_Send(&ne_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                
                MPI_Send(&it->second.data()[0], ne_send, MPI_INT, dest, 9876*7777+dest*888, comm);

                MPI_Send(&iee_send.data()[0], ne_send*nvar, MPI_DOUBLE, dest, 9876*6666+dest*8888, comm);

                //delete[] iee_send;
            }
        }
        if(iee_schedule->RecvRankFromRank[q].find( rank ) != iee_schedule->RecvRankFromRank[q].end())
         {
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
//          int*    recv_back_ids_arr   = new int[n_recv_back];
//          double* recv_back_iee_arr   = new double[n_recv_back*nvar];
             
            std::vector<int> recv_back_ids_arr(n_recv_back,0.0);
            std::vector<double> recv_back_iee_arr(n_recv_back*nvar,0.0);
             
            MPI_Recv(&recv_back_ids_arr.data()[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&recv_back_iee_arr.data()[0], n_recv_back*nvar, MPI_DOUBLE, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Niee[q]     = n_recv_back;
            recv_back_el_ids[q]   = recv_back_ids_arr;
            recv_back_iee[q]      = recv_back_iee_arr;

         }
    }
    
    
//
    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=recv_back_Niee.begin();iter!=recv_back_Niee.end();iter++)
    {
        int L = iter->second;
        
        for(int s=0;s<L;s++)
        {
            el_id = recv_back_el_ids[iter->first][s];
            std::vector<double> StateVec(nvar,0.0);
            for(int p=0;p<nvar;p++)
            {
                StateVec[p]=recv_back_iee[iter->first][s*nvar+p];
                //std::cout << "receive " << recv_back_iee[iter->first][s*nvar+p] << " ";
            }
            //std::cout << std::endl;
            
            U[el_id] = StateVec;
        }
    }
    /**/
}


std::map<int,std::map<int,double> > PartObject::GetNode2ElementMap()
{
    std::map<int,std::map<int,double> > Vert2ElemMap;
    std::map<int,std::vector<int> >::iterator itm;
    for(itm=m_Elem2Elem.begin();itm!=m_Elem2Elem.end();itm++)
    {
        int elid = itm->first;
        int Nv   = m_Elem2Vert[elid].size();
        int nadj = m_Elem2Elem[elid].size();

        std::vector<double> Pijk(Nv*3);
        for(int k=0;k<Nv;k++)
        {
            int gvid     = m_Elem2Vert[elid][k];
            Pijk[k*3+0]  = m_LocalVertsMap[gvid][0];
            Pijk[k*3+1]  = m_LocalVertsMap[gvid][1];
            Pijk[k*3+2]  = m_LocalVertsMap[gvid][2];
        }

        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,Nv);

        for(int j=0;j<Nv;j++)
        {
            int gvid = m_Elem2Vert[elid][j];

            if(m_LocalVertsMap.find(gvid)!=m_LocalVertsMap.end()
                    && Vert2ElemMap[gvid].find(elid)==Vert2ElemMap[gvid].end())
            {
                double dx = Vijk[0]-m_LocalVertsMap[gvid][0];
                double dy = Vijk[1]-m_LocalVertsMap[gvid][1];
                double dz = Vijk[2]-m_LocalVertsMap[gvid][2];
                double d  = sqrt(dx*dx+dy*dy+dz*dz);
                Vert2ElemMap[gvid].insert(std::pair<int,double>(elid,d));
            }
        } 
    }

    return Vert2ElemMap;
}



std::map<int,std::vector<double> > PartObject::ReduceStateVecToVertices(std::map<int,std::map<int,double> > Vert2ElemMap,
                                                                                std::map<int,std::vector<double> > Umap,
                                                                                int nvar)
{
    std::map<int,std::map<int,double> >::iterator n2eit;
    std::map<int,std::vector<double> > nodevals;
    
    for(n2eit=Vert2ElemMap.begin();n2eit!=Vert2ElemMap.end();n2eit++)
    {
        int gvid                        = n2eit->first;
        std::map<int,double> elidsdist  = n2eit->second;
        double distsum                  = 0.0;

        std::map<int,double>::iterator its;
        std::vector<double> uval(nvar,0.0);
        
        for(its=elidsdist.begin();its!=elidsdist.end();its++)
        {
            int elid        = its->first;
            double dist     = its->second;

            for(int q=0;q<nvar;q++)
            {
                uval[q] = uval[q] + Umap[elid][q]*fabs(dist);
            }
            
            distsum  = distsum+fabs(dist);
        }
        for(int q=0;q<nvar;q++)
        {
            uval[q] = uval[q]/distsum;
        }

        nodevals[gvid] = uval;
        
    }
    return nodevals;
}




std::map<int,std::map<int, double> > PartObject::getElem2ConnectedVertMap()
{
    std::map<int,std::map<int, double> > elem2nodemap;
    std::set<int>::iterator its;
    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {
        int elid = *its;
        int Ne   = m_Elem2Elem[elid].size();
        int Nv   = m_Elem2Vert[elid].size();
        int nadj = m_Elem2Elem[elid].size();

        std::vector<double> Pijk(Nv*3);

        for(int k=0;k<Nv;k++)
        {
            int gvid     = m_Elem2Vert[elid][k];
            Pijk[k*3+0]  = m_LocalVertsMap[gvid][0];
            Pijk[k*3+1]  = m_LocalVertsMap[gvid][1];
            Pijk[k*3+2]  = m_LocalVertsMap[gvid][2];
        }
        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,Nv);

        for(int k=0;k<Nv;k++)
        {
            int gvid = m_Elem2Vert[elid][k];
            if(m_LocalVertsMap.find(gvid)!=m_LocalVertsMap.end()
                    && elem2nodemap[elid].find(gvid)==elem2nodemap[elid].end())
            {
                double dx = Vijk[0]-m_LocalVertsMap[gvid][0];
                double dy = Vijk[1]-m_LocalVertsMap[gvid][1];
                double dz = Vijk[2]-m_LocalVertsMap[gvid][2];
                double d  = sqrt(dx*dx+dy*dy+dz*dz);
                elem2nodemap[elid].insert(std::pair<int,double>(gvid,d));
            }
        }

        for(int s=0;s<Ne;s++)
        {
            int adjid = m_Elem2Elem[elid][s];
            int Nvadj = m_Elem2Vert[adjid].size();

            std::map<int, double> vrt2dist;

            for(int j=0;j<Nvadj;j++)
            {
                int gvid = m_Elem2Vert[adjid][j];
                if(m_LocalVertsMap.find(gvid)!=m_LocalVertsMap.end()
                    && elem2nodemap[elid].find(gvid)==elem2nodemap[elid].end())
                {
                    double dx = Vijk[0]-m_LocalVertsMap[gvid][0];
                    double dy = Vijk[1]-m_LocalVertsMap[gvid][1];
                    double dz = Vijk[2]-m_LocalVertsMap[gvid][2];
                    double d  = sqrt(dx*dx+dy*dy+dz*dz);
                    elem2nodemap[elid].insert(std::pair<int,double>(gvid,d));
                }
            }
        }
    }

    return elem2nodemap;

}




std::map<int,int> PartObject::getPartMap()
{
    return m_partMap;
}

std::map<int,std::vector<int> > PartObject::getElem2LocalVertMap()
{
    return m_Elem2LocalVert;
}

std::map<int,std::vector<int> > PartObject::getElem2VertMap()
{
    return m_Elem2Vert;
}

std::map<int,std::vector<int> > PartObject::getElem2FaceMap()
{
   return m_Elem2Face;
}

std::map<int,std::vector<int> > PartObject::getFace2VertMap()
{
    return m_Face2Vert;
}

std::map<int,std::vector<int> > PartObject::getElem2ElemMap()
{
    return m_Elem2Elem;
}

std::map<int,std::vector<int> > PartObject::getFace2ElemMap()
{
    return m_Face2Elem;
}

std::map<int,std::vector<double> > PartObject::getElem2DataMap()
{
    return m_Elem2Data;
}

std::map<int,std::set<int> > PartObject::getExtendedAdjacencyData()
{
    return m_Elem2AdjElem;
}

std::map<int, std::vector<double> > PartObject::getElem2CentroidData()
{
    return m_Elem2Centroid;
}

std::map<int,std::vector<double> > PartObject::getGhostFaceVert()
{
    return m_GhostFaceVert;
}

std::map<int,int> PartObject::getLocalVert2GlobalVert()
{
    return m_LocalV2GlobalV;
}

std::map<int,int > PartObject::getElem2RankMap()
{
    return m_Elem2Rank;
}

std::set<int> PartObject::getTraceVertsOnRankMap()
{
    return m_TraceVertsOnRank;
}

std::set<int> PartObject::getTraceFacesOnRankMap()
{
    return m_TraceFacesOnRank;
}

std::set<int> PartObject::getLocalElemSet()
{
    return m_ElemSet;
}

void PartObject::SetStateVec(std::map<int,std::vector<double> > U, int nvar)
{
    m_Elem2Data.clear();
    m_Elem2Data = U;
}

std::map<int,std::vector<double> > PartObject::getLocalVertsMap()
{
    return m_LocalVertsMap;
}

std::set<int> PartObject::GetLocalTraceFacesSet()
{
    return m_TraceFacesOnRank;
}

std::set<int> PartObject::GetLocalTraceVertsSet()
{
    return m_TraceVertsOnRank;
}


std::map<int,std::vector<int> > PartObject::GetLocalTraceFace2LeftRight()
{
    return m_TraceFace2Elem;
}
std::map<int,int> PartObject::getGlobalElement2Rank()
{
    return m_partGlobalRoot;
}
std::map<int,int> PartObject::GetElement2TypeOnRankMap()
{
    return m_elem2type_on_rank;
}

