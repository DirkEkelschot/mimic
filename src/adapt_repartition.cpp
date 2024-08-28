#include "adapt_repartition.h"
#include "adapt_compute.h"
RepartitionObject::RepartitionObject(mesh* meshInput,
                        std::map<int,std::vector<int> > elements2verts,
                        std::map<int,std::vector<int> > elements2faces,
                        std::map<int,std::vector<int> > elements2elements,
                        std::map<int,int> element2type,
                        std::map<int,std::vector<double> > data,
                        int nAdjLayer,
                        bool reconstruct_ifn,
                        MPI_Comm comm)
{
    //int faceidtrace = 2353748;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    eloc = 0;
    vloc = 0;
    floc = 0;

    Ne_glob = meshInput->nElem;
    Nf_glob = meshInput->nFace;
    Nv_glob = meshInput->nVert;

    std::vector<int> elTypes = meshInput->elTypes;

    // elType[0] = tet;
    // elType[1] = prism;
    // elType[2] = pyramid;
    // elType[3] = hex;

    std::map<int, std::vector<int> > elements2verts_optimal;
    std::map<int, std::vector<int> > elements2faces_optimal;
    std::map<int, std::vector<int> > elements2elements_optimal;
    std::map<int, std::vector<double> > elements2data_optimal;
    std::map<int, int> elements2type_optimal;
    
    clock_t start1, end1;
    double dur_max1,time_taken1;

    start1 = clock();

    GetOptimalDistributionSchedule(elements2verts, 
                                   elements2faces, 
                                   elements2elements, 
                                   element2type,
                                   data,
                                   comm,
                                   elements2type_optimal, 
                                   elements2verts_optimal, 
                                   elements2faces_optimal, 
                                   elements2elements_optimal,
                                   elements2data_optimal);

    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << "Time taken to execute GetOptimalDistributionSchedule : " << std::fixed
        << dur_max1 << std::setprecision(16);
        std::cout << " sec " << std::endl;
    }
    
    

    start1 = clock();


    DeterminePartitionLayout(elements2verts_optimal,meshInput->element2rank,elTypes,comm);
    

    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << "Time taken to execute DeterminePartitionLayout : " << std::fixed
        << dur_max1 << std::setprecision(16);
        std::cout << " sec " << std::endl;
    }

    start1 = clock();

    DetermineElement2ProcMap(elements2verts_optimal, 
                             elements2faces_optimal,
                             elements2elements_optimal, 
                             elements2data_optimal, 
                             elements2type_optimal,
                             meshInput->xcn,
                             Nf_glob,
                             Nv_glob, 
                             comm,
                             elements2verts_update, 
                             elements2faces_update, 
                             elements2elements_update,
                             elements2data_update);
    
    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << "Time taken to execute DetermineElement2ProcMap : " << std ::fixed
        << dur_max1 << std::setprecision(16);
        std::cout << " sec " << std::endl;
    }
    


     CommunicateAdjacencyInfo(comm, elements2elements_update, Ne_glob);





    // The partitioning happens over the elements. 
    // The faces still live equally distributed over the available ranks.
    // They need to be distributed over the ranks.

    // std::map<int,std::vector<int> > face2elements_update;
    start1 = clock();

    // std::map<int,std::vector<int> > face2reference_update2;
    // std::map<int,std::vector<int> > face2Nverts_update2;
    
    // std::map<int,std::vector<int> > face2elements_update2;


    // getFace2EntityPerPartitionMap(elements2faces_update,
    //                               meshInput->ife,
    //                               meshInput->if_ref,
    //                               meshInput->if_Nv,
    //                               meshInput->ifn,
    //                               Nf_glob,
    //                               face2elements_update2,
    //                               face2reference_update2,
    //                               face2verts_update2,
    //                               comm);
    
    // end1 = clock();
    // time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    // MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    // if(rank == 0)
    // {
    //     std::cout << "Time taken to execute getFace2EntityPerPartitionMap : " << std::fixed
    //     << dur_max1 << std::setprecision(16);
    //     std::cout << " sec " << std::endl;
    // }

    start1 = clock();
    ReconstructFace2ElementMap(elements2faces_update,
                               meshInput->ife,
                               face2elements_update,
                               comm);

    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << "Time taken to execute ReconstructFace2ElementMap : " << std::fixed
        << dur_max1 << std::setprecision(16);
        std::cout << " sec " << std::endl;
    }

    start1 = clock();
    std::map<int,std::vector<int> > face2verts_update2;

    if(reconstruct_ifn==true)
    {
        // getFace2VertexPerPartitionMap(elements2faces_update,
        //                           meshInput->ifn,
        //                           meshInput->if_Nv,
        //                           Nf_glob,
        //                           face2verts_update,
        //                           comm);

        ReconstructFace2VertexMap(elements2faces_update,
                                  face2verts_update,
                                  comm);

        //std::cout << "face2verts_update2 " << face2verts_update.size() << " " << face2verts_update2.size() << std::endl;
    }
    else
    {
         getFace2VertexPerPartitionMap(elements2faces_update,
                                  meshInput->ifn,
                                  meshInput->if_Nv,
                                  Nf_glob,
                                  face2verts_update,
                                  comm);
    }
   
    
    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(rank == 0)
    {
        std::cout << "Time taken to execute getFace2EntityPerPartitionMap : " << std::fixed
        << dur_max1 << std::setprecision(16);
        std::cout << " sec " << std::endl;
    }


    std::map<int,std::vector<int> >::iterator itm;
    int cnt = 0;
    for(itm=elements2elements_update.begin();itm!=elements2elements_update.end();itm++)
    {
        int elid = itm->first;
        int ne = itm->second.size();

        for(int q=0;q<ne;q++)
        {
            int elemid = itm->second[q];
            int faceid = elements2faces_update[elid][q];
            std::vector<int> leftright(2,0);
            // Change this so that it bases this of the element type.
            if(adjacent2pid.find(elemid)!=adjacent2pid.end())
            {
                if(adjacent2pid[elemid]==-1)
                {
                    m_loc_trace_faces.insert(faceid);

                    int nv = face2verts_update[faceid].size();
                    for(int n=0;n<nv;n++)
                    {
                        int vid = face2verts_update[faceid][n];

                        if(m_loc_trace_verts.find(vid)==m_loc_trace_verts.end())
                        {
                            m_loc_trace_verts.insert(vid);
                            cnt++;
                        }
                        m_loc_trace_face2vertmap[faceid]    = face2verts_update[faceid];
                        leftright[0] = elid;
                        leftright[1] = elemid;
                        m_loc_trace_face2leftright[faceid]  = leftright;
                    }
                }
            }

            /*
            if(rank != 0)
            {
                if(adjacent2pid.find(elemid)!=adjacent2pid.end())
                {
                    if(adjacent2pid[elemid]==-1)
                    {
                        m_loc_trace_faces.insert(faceid);

                        int nv = face2verts_update[faceid].size();
                        for(int n=0;n<nv;n++)
                        {
                            int vid = face2verts_update[faceid][n];

                            if(m_loc_trace_verts.find(vid)==m_loc_trace_verts.end())
                            {
                                m_loc_trace_verts.insert(vid);
                                cnt++;
                            }
                            m_loc_trace_face2vertmap[faceid]    = face2verts_update[faceid];
                            leftright[0] = elid;
                            leftright[1] = elemid;
                            m_loc_trace_face2leftright[faceid]  = leftright;
                        }
                    }
                }
            }
            else
            {
                if(part_global.find(elemid)==part_global.end() && elemid < Ne_glob)
                {
                    m_loc_trace_faces.insert(faceid);

                    int nv = face2verts_update[faceid].size();
                    for(int n=0;n<nv;n++)
                    {
                        int vid = face2verts_update[faceid][n];

                        if(m_loc_trace_verts.find(vid)==m_loc_trace_verts.end())
                        {
                            m_loc_trace_verts.insert(vid);
                            cnt++;
                        }
                        m_loc_trace_face2vertmap[faceid]    = face2verts_update[faceid];
                        leftright[0] = elid;
                        leftright[1] = elemid;
                        m_loc_trace_face2leftright[faceid]  = leftright;
                    }
                }

            }

            */
            

            // if(part_global.find(elemid)==part_global.end() && elemid < Ne_glob)
            // {
            //     m_loc_trace_faces.insert(faceid);

            //     int nv = face2verts_update[faceid].size();
            //     for(int n=0;n<nv;n++)
            //     {
            //         int vid = face2verts_update[faceid][n];

            //         if(m_loc_trace_verts.find(vid)==m_loc_trace_verts.end())
            //         {
            //             m_loc_trace_verts.insert(vid);
            //             cnt++;
            //         }
            //         m_loc_trace_face2vertmap[faceid]    = face2verts_update[faceid];
            //         leftright[0] = elid;
            //         leftright[1] = elemid;
            //         m_loc_trace_face2leftright[faceid]  = leftright;
            //     }
            // }
        }
    }


    // if(m_loc_trace_faces.find(faceidtrace)!=m_loc_trace_faces.end())
    // {
    //     std::cout << "faceidtrace is here " << rank << std::endl;
    // }
    // if(face2verts_update.find(faceidtrace)!=face2verts_update.end())
    // {
    //     std::cout << "face2verts_update is here " << rank << std::endl;
    // }

    // if(m_loc_trace_verts.find(807212)!=m_loc_trace_verts.end())
    // {
    //     std::cout << "807212 found in m_loc_trace_verts " << rank << std::endl;
    // }
    // if(unique_tracevert2ref.find(807212)!=unique_tracevert2ref.end())
    // {
    //     std::cout << "807212 found in unique_tracevert2ref " << rank << std::endl;
    // }
    // m_loc_trace_verts.clear();
    // std::map<int,int>::iterator itim;
    // for(itim=unique_tracevert2ref.begin();itim!=unique_tracevert2ref.end();itim++)
    // {
    //     m_loc_trace_verts.insert(itim->first);
    // }


    //m_loc_trace_verts_map_glob  = AllGatherMap_T(m_loc_trace_verts_map,comm);

    // m_loc_trace_verts.clear();
    // std::map<int,int>::iterator itim;
    // for(itim=m_loc_trace_verts_map_glob.begin();itim!=m_loc_trace_verts_map_glob.end();itim++)
    // {
    //     m_loc_trace_verts.insert(itim->first);
    // }
        
    int tracefoundRed = 0;
    int nfacetrace      = m_loc_trace_faces.size();
    //MPI_Allreduce(&nfacetrace, &tracefoundRed, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nfacetrace, &tracefoundRed, 1, MPI_INT, MPI_SUM, comm);
    //if(rank == 0)
    //{
    //    std::cout << "tracefoundRed  " << tracefoundRed << std::endl;
    //}


    // getFace2EntityPerPartition(elements2faces_update, 
    //                            meshInput->ife, 
    //                            Nf_glob, 
    //                            face2elements_update, 
    //                            comm);


    // std::map<int,std::vector<int> > face2reference_update;
    // getFace2EntityPerPartition(elements2faces_update, 
    //                            meshInput->if_ref, 
    //                            Nf_glob, 
    //                            face2reference_update, 
    //                            comm);

    // //std::map<int,std::vector<int> > face2verts_update;
    // getFace2EntityPerPartition(elements2faces_update, 
    //                            meshInput->if_Nv, 
    //                            Nf_glob, 
    //                            face2Nverts_update, 
    //                            comm);

    // //std::map<int,std::vector<int> > face2verts_update;
    // getFace2EntityPerPartition(elements2faces_update, 
    //                            meshInput->ifn, 
    //                            Nf_glob, 
    //                            face2verts_update, 
    //                            comm);




    // vloc = LocalVertsMap.size();
    
    for(int i=0;i<nAdjLayer;i++)
    {
        
        //std::cout << rank << " BEFORE elements2verts_update " << i << " " << elements2verts_update.size() << std::endl;
        start1 = clock();
        std::map<int,std::vector<int> > adjEl2Face = getAdjacentElementLayer(elements2verts_update,
                                                                        elements2faces_update,
                                                                        elements2elements_update,
                                                                        meshInput->xcn, 
                                                                        elements2data_update,
                                                                        Ne_glob,
                                                                        Nf_glob,
                                                                        Nv_glob, 
                                                                        comm,
                                                                        elements2verts_update, 
                                                                        elements2faces_update, 
                                                                        elements2elements_update,
                                                                        elements2data_update);
        
        //std::cout << rank << " AFTER elements2verts_update " << i << " " << elements2verts_update.size() << std::endl;

        end1 = clock();
        time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
        MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(rank == 0)
        {
            std::cout << "Time taken to execute getAdjacentElementLayer : " << std::fixed
            << dur_max1 << std::setprecision(16);
            std::cout << " sec " << std::endl;
        }

        start1 = clock();
        

        updateFace2EntityPerPartition(adjEl2Face, 
                                    meshInput->ifn,
                                    meshInput->if_Nv,  
                                    Nf_glob, 
                                    face2verts_update, 
                                    comm);
        
        end1 = clock();
        time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
        MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(rank == 0)
        {
            std::cout << "Time taken to execute updateFace2EntityPerPartition : " << std::fixed
            << dur_max1 << std::setprecision(16);
            std::cout << " sec " << std::endl;
        }
    }

    
}



void RepartitionObject::GetOptimalDistributionSchedule(std::map<int,std::vector<int> > elements2verts,
                                                       std::map<int,std::vector<int> > elements2faces,
                                                       std::map<int,std::vector<int> > elements2elements,
                                                       std::map<int,int> element2type,
                                                       std::map<int,std::vector<double> > data,
                                                       MPI_Comm comm,
                                                       std::map<int,int> &element2type_update,
                                                       std::map<int,std::vector<int> >& elements2verts_update,
                                                       std::map<int,std::vector<int> >& elements2faces_update,
                                                       std::map<int,std::vector<int> >& elements2elements_update,
                                                       std::map<int,std::vector<double> >& elements2data_update)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::vector<int> iniDist(world_size,0);
    std::vector<int> red_iniDist(world_size,0);
    std::vector<int> ini_offsetDist(world_size,0);
    
    int nElem = elements2verts.size();
    // int nVpEl = 0;
    // int nFpEl = 0;
    // int nEpEl = 0;
    int ndata = 0;

    elements2verts_update.clear();
    std::map<int,std::vector<int> > elements2verts_update_v2;

  

    if(nElem != 0)
    {
        // std::map<int,std::vector<int> >::iterator itm;
        
        // // for(itm=elements2verts.begin();itm!=elements2verts.end();itm++)
        // // {
        // //     int elid = itm->first;
        // // }
        // nVpEl = elements2verts.begin()->second.size();
        // nFpEl = elements2faces.begin()->second.size();
        ndata = data.begin()->second.size();
    }

    for(int i=0;i<world_size;i++)
    {        
        if(i==world_rank)
        {
            iniDist[i] = elements2verts.size();
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

    std::map<int,std::vector<int> >::iterator itmiv;
    if(DoNotUpdatePartition.find(world_rank)!=DoNotUpdatePartition.end())
    {
        std::map<int,std::vector<int> >::iterator itmiv;
        int eloc = 0;
        for(itmiv=elements2verts.begin();itmiv!=elements2verts.end();itmiv++)
        {
            int gElId = itmiv->first;
            int gtype = element2type[gElId];

            elements2verts_update[gElId]      = itmiv->second;
            elements2verts_update_v2[gElId]   = itmiv->second;
            element2type_update[gElId]        = gtype;
            // for(int j=0;j<elements2verts.size();j++)
            // {
                

            //     // std::vector<int> row(itmiv->second.size(),0);
            //     // for(int p=0;p<itmiv->second.size();p++)
            //     // {
            //     //     row[p] = itmiv->second[p];
            //     // }
            //     // elements2verts_update[gElId] = row;

            // }
            // for(int j=0;j<elements2faces.size();j++)
            // {
                
            // }

            elements2faces_update[gElId]    = elements2faces[gElId];
            elements2elements_update[gElId] = elements2elements[gElId];
            elements2data_update[gElId]     = data[gElId];

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

            //std::vector<int> rows_ien(Nel*nVpEl);
            //ien_to_send.push_back(rows_ien);
            //std::vector<int> rows_ief(Nel*nFpEl);
            //ief_to_send.push_back(rows_ief);
            //std::vector<int> rows_iee(Nel*nFpEl);
            //iee_to_send.push_back(rows_iee);

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

        for(itmiv=elements2verts.begin();itmiv!=elements2verts.end();itmiv++)
        {
            int nelPrank    = NeltoRanks[cc];

            int gElId               = itmiv->first;
            int update_type         = element2type[gElId];
            int nv_el               = itmiv->second.size();
            int nf_el               = elements2faces[gElId].size();
            int ndata               = data[gElId].size();

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
                        //ien_to_send[cc][nVpEl*t+j] = itmiv->second[j];
                        ien_to_send[cc].push_back(itmiv->second[j]);
                        // ien_to_send[cc][nv_el*t+j] = itmiv->second[j];
                    }

                    for(int j=0;j<nf_el;j++)
                    {
                        //ief_to_send[cc][nFpEl*t+j] = elements2faces[gElId][j];
                        //iee_to_send[cc][nFpEl*t+j] = elements2elements[gElId][j];

                        ief_to_send[cc].push_back(elements2faces[gElId][j]);
                        iee_to_send[cc].push_back(elements2elements[gElId][j]);
                        
                        
                        // ief_to_send[cc][nf_el*t+j] = elements2faces[gElId][j];
                        // iee_to_send[cc][nf_el*t+j] = elements2elements[gElId][j];
                    }

                    // if(world_rank == 1)
                    // {
                    //     std::cout << "SEND 1 ::";
                    //     for(int s=0;s<nf_el;s++)
                    //     {
                    //         std::cout << elements2elements[gElId][s] << " ";
                    //     }
                    //     std::cout << std::endl;
                    // }

                    for(int j=0;j<ndata;j++)
                    {
                        data_to_send[cc][ndata*t+j] = data[gElId][j];
                        //std::cout << "data[gElId][j] " << data[gElId][j] << std::endl;
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
                    update_e2f_row[j] = elements2faces[gElId][j];
                    update_e2e_row[j] = elements2elements[gElId][j];
                }

                // if(world_rank == 1)
                // {
                //     std::cout << "SEND 2 :: ";
                //     for(int s=0;s<nf_el;s++)
                //     {
                //         std::cout << elements2elements[gElId][s] << " ";
                //     }
                //     std::cout << std::endl;
                // }

                std::vector<double> update_data_row(ndata,0);

                for(int j=0;j<ndata;j++)
                {
                    update_data_row[j] = data[gElId][j];
                    //std::cout << "update_data_row[j] " << update_data_row[j] << std::endl;
                }
                element2type_update[gElId]      = update_type;
                elements2verts_update[gElId]    = update_e2v_row;
                elements2faces_update[gElId]    = update_e2f_row;
                elements2elements_update[gElId] = update_e2e_row;
                elements2data_update[gElId]     = update_data_row;
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

            MPI_Send(&n_Elem         ,      1,     MPI_INT, dest, dest,           comm);
            MPI_Send(&ElIDs[0]       , n_Elem,     MPI_INT, dest, dest*123000,    comm);
            MPI_Send(&ElTypes[0]     , n_Elem,     MPI_INT, dest, dest*128000,    comm);
            MPI_Send(&NvEl[0]        , n_Elem,     MPI_INT, dest, dest*492000,    comm);
            MPI_Send(&NfEl[0]        , n_Elem,     MPI_INT, dest, dest*984000,    comm);

            MPI_Send(&n_Vrt          ,      1,     MPI_INT, dest, dest*1968000,   comm);
            MPI_Send(&El2v[0]       ,   n_Vrt,     MPI_INT, dest, dest*3936000,   comm);

            MPI_Send(&n_Fce          ,      1,     MPI_INT, dest, dest*7872000,   comm);
            MPI_Send(&El2f[0]       ,   n_Fce,     MPI_INT, dest, dest*15744000,  comm);
            MPI_Send(&El2e[0]       ,   n_Fce,     MPI_INT, dest, dest*31488000,  comm);

            MPI_Send(&n_Data        ,      1,      MPI_INT, dest, dest*62976000,   comm);
            MPI_Send(&El2d[0]       ,   n_Data,    MPI_DOUBLE, dest, dest*125952,  comm);
            MPI_Send(&NdataEl[0]    ,   n_Elem,    MPI_INT, dest, dest*251904,  comm);
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
        for(itmiv=elements2verts.begin();itmiv!=elements2verts.end();itmiv++)
        {
            int gElId                        = itmiv->first;
            int type_update                  = element2type[gElId]; 
            elements2verts_update[gElId]     = itmiv->second;
            elements2faces_update[gElId]     = elements2faces[gElId];
            elements2elements_update[gElId]  = elements2elements[gElId];
            elements2data_update[gElId]      = data[gElId];
            element2type_update[gElId]       = type_update; 
        }

        int ntot = elements2verts.size();
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
                int type = recvd_elem_type[itmiv->first][q];
                int nvpel = recvd_elem_nvs[itmiv->first][q];
                int nfpel = recvd_elem_nfs[itmiv->first][q];
                int ndata = recvd_elem_ndata[itmiv->first][q];
                std::vector<int> elements2verts_update_row(nvpel,0);
                std::vector<int> elements2faces_update_row(nfpel,0);
                std::vector<int> elements2elements_update_row(nfpel,0);
                std::vector<double> elements2data_update_row(ndata,0);
                for(int s=0;s<nvpel;s++)
                {
                    elements2verts_update_row[s] = recvd_elem_vert_ids[itmiv->first][accul_v+s];
                }
                for(int s=0;s<nfpel;s++)
                {
                    elements2faces_update_row[s]    = recvd_elem_face_ids[itmiv->first][accul_f+s];
                    elements2elements_update_row[s] = recvd_elem_elem_ids[itmiv->first][accul_f+s];
                }

                for(int s=0;s<ndata;s++)
                {
                    elements2data_update_row[s]     = recvd_elem_data_ids[itmiv->first][accul_d+s];
                    //std::cout << "recvd_elem_data_ids[itmiv->first][accul+s] " << recvd_elem_data_ids[itmiv->first][accul+s] << std::endl;
                }


                element2type_update[gElId]          = type;
                elements2verts_update[gElId]        = elements2verts_update_row;
                elements2faces_update[gElId]        = elements2faces_update_row;
                elements2elements_update[gElId]     = elements2elements_update_row;
                elements2data_update[gElId]         = elements2data_update_row;

                accul_v = accul_v + nvpel;
                accul_f = accul_f + nfpel;
                accul_d = accul_d + ndata;
            }

            ntot    = ntot + nel;
        }
    }

}




void RepartitionObject::DeterminePartitionLayout(std::map<int,std::vector<int> > elements2verts,
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
    
    int nrow = elements2verts.size();
    int nvpEL = elements2verts.begin()->second.size();
    int nloc = nrow;

    //std::cout << " elTypes " << elTypes[0] << " " << elTypes[1] << " " << elTypes[2] << " " << nvpEL << std::endl;
    
    ParallelState_Parmetis_Lite* pstate_parmetis = new ParallelState_Parmetis_Lite(elements2verts,  elTypes, comm);

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

    part = std::vector<int>(elements2verts.size(),0);
    int nElemTotal = pstate_parmetis->getNtotalElem();
    //std::vector<int> part_global_vec (nElemTotal,0);
    

    std::vector<int> part_arr_ell_id(nloc,0);
    
    std::vector<int> part_arr_ell_Rid(nloc,0);

    std::map<int,std::vector<int> >::iterator itmiv;
    int i = 0;

    for(itmiv=elements2verts.begin();itmiv!=elements2verts.end();itmiv++)
    {
        int gid                 = itmiv->first;
        part_arr_ell_id[i]      = gid;
        part_arr_ell_Rid[i]     = part_arr[i];
        part_map[gid]           = part_arr[i];
        i++;
    }

    std::vector<int> part_global_vec;
    std::vector<int> part_global_Rid_vec;
    if ( rank == root) 
    {
        part_global_vec.resize(nElemTotal,0);
        part_global_Rid_vec.resize(nElemTotal,0);
    } 

    
    MPI_Gatherv(&part_arr_ell_id.data()[0],
                    nloc, MPI_INT,
                    &part_global_vec.data()[0],
                    pstate_parmetis->getNlocs().data(),
                    pstate_parmetis->getElmdist().data(),
                    MPI_INT,root,comm);

    MPI_Gatherv(&part_arr_ell_Rid.data()[0],
                    nloc, MPI_INT,
                    &part_global_Rid_vec.data()[0],
                    pstate_parmetis->getNlocs().data(),
                    pstate_parmetis->getElmdist().data(),
                    MPI_INT,root,comm);

    // MPI_Allgatherv(&part_arr_ell_id.data()[0],
    //                 nloc, MPI_INT,
    //                 &part_global_vec.data()[0],
    //                 pstate_parmetis->getNlocs().data(),
    //                 pstate_parmetis->getElmdist().data(),
    //                 MPI_INT,comm);

    // MPI_Allgatherv(&part_arr_ell_Rid.data()[0],
    //                 nloc, MPI_INT,
    //                 &part_global_Rid_vec.data()[0],
    //                 pstate_parmetis->getNlocs().data(),
    //                 pstate_parmetis->getElmdist().data(),
    //                 MPI_INT,comm);

    for(int i=0;i<part_global_vec.size();i++)
    {
        int eid = part_global_vec[i];
        int rid = part_global_Rid_vec[i];

        part_global[eid] = rid;
    }

    glob_elem_size = part_global.size();
    MPI_Bcast(&glob_elem_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] xadj_par;
    delete[] adjncy_par;
}



void RepartitionObject::DetermineElement2ProcMap(std::map<int,std::vector<int> >     ien, 
                                                 std::map<int,std::vector<int> >     ief,
                                                 std::map<int,std::vector<int> >     iee,
                                                 std::map<int,std::vector<double> >  data,
                                                 std::map<int,int> element2type, 
                                                 std::map<int,std::vector<double> >  xcn,
                                                 int Nf_glob,
                                                 int Nv_glob, 
                                                 MPI_Comm comm,
                                                 std::map<int,std::vector<int> >& elements2verts_update,
                                                 std::map<int,std::vector<int> >& elements2faces_update,
                                                 std::map<int,std::vector<int> >& elements2elements_update,
                                                 std::map<int,std::vector<double> >& elements2data_update)
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
    ParallelState* xcn_pstate = new ParallelState(Nv_glob,comm);

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
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
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

    for(itmii=part_map.begin();itmii!=part_map.end();itmii++)
    {
        el_id   = itmii->first;
        p_id    = itmii->second;
        int el_type = element2type[el_id];
    
        nvPerEl = ien[el_id].size();
        nfPerEl = ief[el_id].size();
        ndataPerEl = data[el_id].size();

        if(p_id!=rank) // If element is not on this rank and needs to be send to other rank (p_id), add it to rank to element map.
        {
            elms_to_send_to_ranks[p_id].push_back(el_id); // rank to element map.
            nvPerElms_to_send_to_ranks[p_id].push_back(nvPerEl);
            nfPerElms_to_send_to_ranks[p_id].push_back(nfPerEl);
            ndataPerElms_to_send_to_ranks[p_id].push_back(ndataPerEl);
            elmtype_to_send_to_ranks[p_id].push_back(el_type);
            //====================Hybrid=======================
            for(int k=0;k<nvPerEl;k++)//This works for hexes.
            {
                v_id = ien[el_id][k];

                vertIDs_to_send_to_ranks[p_id].push_back(v_id);

            }// We do not care about the vertices for these elements since they are needed on other ranks anyways.
            //std::cout << "Check iee :: ";
            for(int k=0;k<nfPerEl;k++)//This works for hexes.
            {
                f_id = ief[el_id][k];
                ea_id = iee[el_id][k];

                //std::cout << ea_id << " ";
                faceIDs_to_send_to_ranks[p_id].push_back(f_id);
                elemIDs_to_send_to_ranks[p_id].push_back(ea_id);
            }

            //std::cout << std::endl;
            for(int k=0;k<ndataPerEl;k++)//This works for hexes.
            {
                double dataV = data[el_id][k];
                //std::cout << "dataV " << dataV << std::endl; 
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
                v_id = ien[el_id][k];
                
                
                //elem.push_back(v_id);
                if(unique_vertIDs_on_rank_set.find( v_id ) == unique_vertIDs_on_rank_set.end() && v_id != -1)// find the unique vertices that need to be send to other partitions.
                {
                    unique_vertIDs_on_rank_set.insert(v_id);
                    //unique_verts_on_rank_vec.push_back(v_id);
                    
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

            


            // std::cout << "Check iee :: ";
            // for(int k=0;k<nfPerEl;k++)//This works for hexes.
            // {
            //     ea_id = iee[el_id][k];

            //     std::cout << ea_id << " ";
            // }

            // std::cout << std::endl;
            
            loc_r_elem.push_back(el_id);
            loc_r_Ndata.push_back(ndataPerEl);
            for(int k=0;k<ndataPerEl;k++)//This works for hexes.
            {
                double dataV = data[el_id][k];
                //std::cout << "dataV " << dataV << std::endl; 
                loc_r_data.push_back(dataV);
            }
            

            loc_r_nv_elem.push_back(nvPerEl);
            loc_r_nf_elem.push_back(nfPerEl);
            elem_set.insert(el_id);
            loc_r_elem_set.insert(el_id);
            elem_map[el_id] = on_rank;
            elem2type_on_rank[el_id] = el_type;
            
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
            
            if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_vertIDs_on_rank_set.insert(v_id_n);
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
            if(unique_faceIDs_on_rank_set.find( f_id_n ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
            {
                unique_faceIDs_on_rank_set.insert(f_id_n);
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
        elem_set.insert(TotRecvElement_IDs[i]);
        loc_r_elem_set.insert(TotRecvElement_IDs[i]);
        elem_map[el_id] = on_rank;
        //loc_data[el_id] = rowdata;
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
                offset_xcn        = xcn_pstate->getOffset(rank);
                for(int u=0;u<it->second.size();u++)
                {
                    if(xcn.find(it->second[u])==xcn.end())
                    {
                        std::cout << "NOT FOUND" << std::endl;
                        std::cout << " it->second[u] " << it->second[u] << " " << dest << " " << rank  << " " << Nv_glob << std::endl;
                        nfound++;
                    }
                    else
                    {
                        vert_send[u*3+0]=xcn[it->second[u]][0];
                        vert_send[u*3+1]=xcn[it->second[u]][1];
                        vert_send[u*3+2]=xcn[it->second[u]][2];
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
        V[0] = xcn[gvid][0];
        V[1] = xcn[gvid][1];
        V[2] = xcn[gvid][2];
        LocalVertsMap[gvid] = V;

        o_lvertex2gvertex_part[lvid] = gvid;
        //o_gvertex2lvertex_part[gvid] = lvid;

        // o_lvertex2gvertex[lvid] = gvid;
        // o_gvertex2lvertex[gvid] = lvid;
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
            LocalVertsMap[gvid] = V;

            o_lvertex2gvertex_part[lvid] = gvid;
            //o_gvertex2lvertex_part[gvid] = lvid;

            // o_lvertex2gvertex[lvid]=gvid;
            // o_gvertex2lvertex[gvid]=lvid;
           
            m++;
            lvid++;
        }
    }

    nLoc_Verts = LocalVertsMap.size();
    // ================================== Faces on Rank =========================================
    
    int lfid = 0;
    int gfid = 0;
    for(m=0;m<floc_tmp;m++)
    {
        gfid = faceIDs_on_rank[m];
    
        lface2gface[lfid] = gfid;
        gface2lface[gfid] = lfid;
        lfid++;
    }
    
    // ================================== Faces on Rank =========================================
    //NlocElem             = loc_r_elem.size()+Nel_extra+itel;
    //LocalElem2GlobalVert = new Array<int>(NlocElem,8);
    //lelement2lvertex  = new Array<int>(NlocElem,8);
    //std::vector<double> U0vert;
    //U0Elem               = new Array<double>(NlocElem,1);
    //U0Vert               = new Array<double>(LocalVerts.size(),1);
    //ElemPart             = new Array<int>(NlocElem,1);

    int glob_v = 0;
    int loc_v  = 0;
    int glob_f = 0;
    int glob_e = 0;
    int loc_f  = 0;
    double varia_v = 0.0;
    int ndatas = 0;
    for(m=0;m<loc_r_elem.size();m++)
    {
        el_id   = loc_r_elem[m];
        nvPerEl = loc_r_nv_elem[m];
        nfPerEl = loc_r_nf_elem[m];
        ndatas  = loc_r_Ndata[m];
        Loc_Elem.push_back(el_id);
        Loc_Elem_Set.insert(el_id);

        LocElem2Nv[el_id]    = nvPerEl;
        LocElem2Nf[el_id]    = nfPerEl;
        LocElem2Ndata[el_id] = ndatas;

        Loc_Elem_Nv.push_back(nvPerEl);
        Loc_Elem_Nf.push_back(nfPerEl);
        LocAndAdj_Elem.push_back(el_id);
        LocAndAdj_Elem_Nv.push_back(nvPerEl);
        LocAndAdj_Elem_Nf.push_back(nfPerEl);
        loc_data[el_id] = data[el_id];
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        int ndatas = data[el_id].size();
        std::vector<double> rowdata(ndatas,0.0);
        for(int k=0;k<ndatas;k++)
        {
            rowdata[k] = data[el_id][k];
        }

        eloc++;
        std::vector<int> tmp_globv;
        std::vector<int> tmp_locv;
        std::vector<int> tmp_globf;
        std::vector<int> tmp_globe;
        std::vector<int> tmp_loce;

        for(int p=0;p<nvPerEl;p++)
        {
            //gface2lface
            //o_lvertex2gvertex
            glob_v = ien[el_id][p];
            //glob_v = ien->getVal(el_id-ien_o,p);
            if(glob_v>xcn_pstate->getNel())
            {
                std::cout << "Nel On error " << glob_v << std::endl;
            }
            
            loc_v  = o_gvertex2lvertex[glob_v];
            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);

            //LocalElem2GlobalVert->setVal(m,p,glob_v);
            //lelement2lvertex->setVal(m,p,loc_v);
            //collect_var[loc_v].push_back(rho_v);

            //globElem2globVerts[el_id].push_back(glob_v);
            globVerts2globElem[glob_v].push_back(el_id);
            globElem2locVerts[el_id].push_back(loc_v);
        }

        // if(rank==0)
        // {
        //     std::cout << "Check iee 22 : ";
        // }

        for(int p=0;p<nfPerEl;p++)
        {
            glob_f = ief[el_id][p];
            glob_e = iee[el_id][p];
            loc_f  = gface2lface[glob_f];
            globElem2localFaces[el_id].push_back(loc_f);
            globElem2globFaces[el_id].push_back(glob_f);
            globFace2GlobalElements[glob_f].push_back(el_id);
            tmp_globf.push_back(glob_f);
            tmp_globe.push_back(glob_e);
            int ea_id = iee[el_id][p];

        }
        // if(rank==0)
        // {
        //     std::cout << std::endl;
        // }
        //elem2type_on_rank[el_id] = el_type;


        //=======================================
        elements2verts_update[el_id]    = tmp_globv;
        elements2faces_update[el_id]    = tmp_globf;
        elements2elements_update[el_id] = tmp_globe;
        elements2data_update[el_id]     = rowdata;
        //=======================================

        //LocalElem2GlobalVert.push_back(tmp_globv);
        lelement2lvertex.push_back(tmp_locv);
        tmp_globv.clear();
        tmp_locv.clear();
    }

    int cnv = 0;
    int cnf = 0;
    cnt_data = 0;
    for(m=0;m<Nel_extra;m++)
    {
        el_id                             = TotRecvElement_IDs[m];
        int el_type                       = TotRecvElement_Types[m];
        nvPerEl                           = TotRecvElement_NVs[m];
        nfPerEl                           = TotRecvElement_NFs[m];
        int ndatas                        = TotRecvElement_Ndatas[m];
        elem2type_on_rank[el_id]          = el_type;
        Loc_Elem_Set.insert(el_id);
        LocElem2Nv[el_id]                 = nvPerEl;
        LocElem2Nf[el_id]                 = nfPerEl;
        LocElem2Ndata[el_id]              = ndatas;
        LocalElement2GlobalElement[eloc]  = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;
        Loc_Elem.push_back(el_id);
        Loc_Elem_Nv.push_back(nvPerEl);
        Loc_Elem_Nf.push_back(nfPerEl);
        LocAndAdj_Elem.push_back(el_id);
        LocAndAdj_Elem_Nv.push_back(nvPerEl);
        LocAndAdj_Elem_Nf.push_back(nfPerEl);

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
        loc_data[el_id] = rowdata;
        cnt_data = cnt_data + ndatas;
        for(int p=0;p<nvPerEl;p++)
        {
            glob_v = TotRecvVerts_IDs[cnv+p];
            if(glob_v>xcn_pstate->getNel())
            {
                std::cout << "Nel Extra error " << glob_v << std::endl;
            }
            loc_v = o_gvertex2lvertex[glob_v];
            //globElem2globVerts[el_id].push_back(glob_v);
            globElem2locVerts[el_id].push_back(loc_v);
            globVerts2globElem[glob_v].push_back(el_id);

            tmp_globv.push_back(glob_v);
            tmp_locv.push_back(loc_v);
            
        }
        for(int p=0;p<nfPerEl;p++)
        {
            glob_f = TotRecvFaces_IDs[cnf+p];
            glob_e = TotRecvElem_IDs[cnf+p];
            tmp_globf.push_back(glob_f);
            tmp_globe.push_back(glob_e);
            loc_f  = gface2lface[glob_f];
            globElem2localFaces[el_id].push_back(loc_f);
            globElem2globFaces[el_id].push_back(glob_f);
            globFace2GlobalElements[glob_f].push_back(el_id);
            

        }

        //=======================================
        elements2verts_update[el_id]    = tmp_globv;
        elements2faces_update[el_id]    = tmp_globf;
        elements2elements_update[el_id] = tmp_globe;
        elements2data_update[el_id]     = rowdata;
        //=======================================


        cnv=cnv+nvPerEl;
        cnf=cnf+nfPerEl;
        //LocalElem2GlobalVert.push_back(tmp_globv);
        lelement2lvertex.push_back(tmp_locv);
        tmp_globv.clear();
        tmp_locv.clear();
    }
    
    // std::map<int,std::vector<int> >::iterator itmvvv;

    // for(itmvvv=globFace2GlobalElements.begin();itmvvv!=globFace2GlobalElements.end();itmvvv++)
    // {
    //     std::cout << " size " << itmvvv->first << " " << itmvvv->second.size() << std::endl;
    // }


    nLoc_Elem = Loc_Elem.size();
    vloc = LocalVertsMap.size();
    floc = cnf;

    /**/

}


void RepartitionObject::CommunicateAdjacencyInfo(MPI_Comm comm, std::map<int,std::vector<int> > element2elements, int Nel_glob)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::map<int,std::vector<int> >::iterator itm;
    int nreq;

    std::set<int> toquery_elids;
        
    for(itm=element2elements.begin();itm!=element2elements.end();itm++)
    {
        int elid = itm->first;

        for(int q=0;q<itm->second.size();q++)
        {
            int adjeid = itm->second[q];

            if(part_map.find(adjeid)==part_map.end() && adjeid<Nel_glob)
            {
                if(toquery_elids.find(adjeid)==toquery_elids.end())
                {
                    toquery_elids.insert(adjeid);
                }
            }
        }
    }

    std::vector<int> toquery_elids_vec(toquery_elids.size(),0);
    //std::cout << "collect_prism_ids_toberequested " << collect_prism_ids_toberequested.size() << " " << world_rank << std::endl;
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

            adjacent2pid[el_id] = p_id;
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
                
                if(part_global.find(QueryOnRoot[j])!=part_global.end())
                {
                    pid_2sendback[j] = part_global[QueryOnRoot[j]];
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
            if(part_global.find(el_id)!=part_global.end())
            {
                p_id = part_global[el_id];
            }

            adjacent2pid[el_id] = p_id;
        }
        
    }

}


void RepartitionObject::updateFace2EntityPerPartition(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > if_Nv_read, 
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ife_loc,
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
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    //int Nel = part_global.size();
    int Nel = glob_elem_size;

    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;
    
    for(itefmap=ief.begin();itefmap!=ief.end();itefmap++)
    {
        for(int q=0;q<itefmap->second.size();q++)
        {
            int face_req = itefmap->second[q];
            
            r = FindRank(new_offsets,size,face_req);
            
            if(r != rank && ife_loc.find(face_req)==ife_loc.end())
            {
                rank2req_Faces[r].push_back(face_req);
            }
            else
            {
                if(ife_loc.find(face_req)==ife_loc.end())
                {
                    int ncol = ife_read[face_req].size();
                    int ncol_d = if_Nv_read[face_req][0];
                    //std::cout << "ncol first " << ncol << " " << ncol_d << std::endl;
                    for(int j=0;j<ncol_d;j++)
                    {
                        int vrtid = ife_read[face_req][j];
                        ife_loc[face_req].push_back(vrtid);
                    }
                }
            }
        }
    }
    
    int own = ife_loc.size();
    
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
                ife_row_loc[r]      = recv_back_ife[iter->first][offset+r];
                ife_loc[face_id]    = ife_row_loc;
            }
                
            
            offset = offset+ncol;
            
        }
        ntotal=ntotal+L;
    }
}










void RepartitionObject::getFace2EntityPerPartition(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ife_read, 
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ife_loc,
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
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    //int Nel = part_global.size();
    int Nel = glob_elem_size;
    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;
    
    for(itefmap=ief.begin();itefmap!=ief.end();itefmap++)
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
                if(ife_loc.find(face_req)==ife_loc.end())
                {
                    int ncol = ife_read[face_req].size();
                    //std::cout << "ncol first " << ncol << std::endl;
                    for(int j=0;j<ncol;j++)
                    {
                        int vrtid = ife_read[face_req][j];
                        ife_loc[face_req].push_back(vrtid);
                        //std::cout << "vrtid " << vrtid << std::endl;
                    }
                }
            }
        }
    }
    
    int own = ife_loc.size();
    
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
                    int ncol = ife_read[it->second[u]].size();
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
        
        for(int s=0;s<L;s++)
        {
            face_id = recv_back_face_ids[iter->first][s];
            int ncol = recv_back_face_Ne[iter->first][s];
            if(ife_loc.find(face_id)==ife_loc.end())
            {

                //std::cout << "recv_back_ife[iter->first][s*ncol+r] ";
                for(int r=0;r<ncol;r++)
                {
                    int vrtid_n = recv_back_ife[iter->first][s*ncol+r];
                    ife_loc[face_id].push_back(recv_back_ife[iter->first][s*ncol+r]);
                    //std::cout << recv_back_ife[iter->first][s*ncol+r] << " ";
                    
                }
                //std::cout << std::endl;
            }
            
        }
        ntotal=ntotal+L;
    }
}


void RepartitionObject::getFace2VertexPerPartitionMap(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ifn_read, 
                                                   std::map<int,std::vector<int> > if_Nv_read, 
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ifn_loc,
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
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    //int Nel = part_global.size();
    int Nel = glob_elem_size;
    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;

    for(itefmap=ief.begin();itefmap!=ief.end();itefmap++)
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
                if(ifn_loc.find(face_req)==ifn_loc.end())
                {
                    int Nv = if_Nv_read[face_req][0];
                    std::vector<int> row(Nv,0);
                    for(int i=0;i<Nv;i++)
                    {
                        row[i] = ifn_read[face_req][i];
                    }
                    ifn_loc[face_req] = row;
                }
            }
        }
    }
    
    int own = ifn_loc.size();
    
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

            ifn_loc[face_id]  = ifref_loc_row;

            offset = offset + ncol;
        }

        ntotal=ntotal+L;
    }
}


void RepartitionObject::ReconstructFace2VertexMap(std::map<int,std::vector<int> > elements2faces_update,
                                  std::map<int,std::vector<int> > &face2verts_update_new,
                                  MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    std::vector<std::vector<int> > tetra_faces   = getTetraFaceMap(); 
    std::vector<std::vector<int> > prism_faces   = getPrismFaceMap(); 
    std::vector<std::vector<int> > pyramid_faces = getPyramidFaceMap(); 
    // std::vector<std::vector<int> > hex_faces     = getHexFaceMap(); 

    //FaceSetPointer m_FaceSetPointer;
    //start = clock();
    std::vector<std::vector<int> > map_faces;
    std::map<int,std::vector<int> > face2vertexmap_a;
    std::map<int,std::vector<int> >::iterator itr;
    for(itr=elements2faces_update.begin();itr!=elements2faces_update.end();itr++)
    {
        int elid = itr->first;
        int nf   = elements2faces_update[elid].size();
        int nv   = elements2verts_update[elid].size(); 
        int type = elem2type_on_rank[elid];

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

        // if(nf!=map_faces.size())
        // {
        //     std::cout << "Error: facemap for element type does not correspond." << std::endl; 
        // }

        for(int u=0;u<map_faces.size();u++)
        {
            int nv  = map_faces[u].size();
            int fid = elements2faces_update[elid][u];

            std::vector<int> face(nv,0);
            for(int w=0;w<nv;w++)
            {
                face[w] = elements2verts_update[elid][map_faces[u][w]];
            }
            face2verts_update_new[fid] = face;
        }
    }
}


void RepartitionObject::ReconstructFace2ElementMap(std::map<int,std::vector<int> > ief,
                                                std::map<int,std::vector<int> > ife_read,
                                                std::map<int,std::vector<int> > &face2elements_update_new,
                                                MPI_Comm comm)
{

    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    // std::vector<std::vector<int> > tetra_faces   = getTetraFaceMap(); 
    // std::vector<std::vector<int> > prism_faces   = getPrismFaceMap(); 
    // std::vector<std::vector<int> > pyramid_faces = getPyramidFaceMap(); 
    // std::vector<std::vector<int> > hex_faces     = getHexFaceMap(); 

    //FaceSetPointer m_FaceSetPointer;
    //start = clock();
    std::vector<std::vector<int> > map_faces;
    std::map<int,std::vector<int> > face2vertexmap_a;
    
    std::map<int,std::vector<int> >::iterator itr;
    for(itr=elements2faces_update.begin();itr!=elements2faces_update.end();itr++)
    {
        int elid = itr->first;
        int nf   = elements2faces_update[elid].size();
        int nv   = elements2verts_update[elid].size(); 
        int type = elem2type_on_rank[elid];

        // switch (type) {
        // case 2:
        //     map_faces = tetra_faces;
        //     break;
        // case 4:
        //     map_faces = hex_faces;
        //     break;
        // case 5:
        //     map_faces = pyramid_faces;
        //     break;
        // case 6:
        //     map_faces = prism_faces;
        //     break;
        // }

        // if(nf!=map_faces.size())
        // {
        //     std::cout << "Error: facemap for element type does not correspond." << std::endl; 
        // }

        for(int q = 0;q < elements2faces_update[elid].size(); q++)
        {
            int fid = elements2faces_update[elid][q];
            face2elements_update_new[fid].push_back(elid);

            // FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face));
            // std::pair<FaceSetPointer::iterator, bool> testInsPointer;
            // testInsPointer = m_FaceSetPointer.insert(facePointer);
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

    for(itmiv=face2elements_update_new.begin();itmiv!=face2elements_update_new.end();itmiv++)
    {
        int faceid   = itmiv->first;
        int numel    = itmiv->second.size();

        if(numel == 1)
        {
            face2elements_update_new[faceid].clear();

            int r = FindRank(new_offsets.data(),size,faceid);

            if(r != rank)
            {
                rank2req_Faces[r].push_back(faceid);
            }
            else
            {
                if(ife_read.find(faceid)!=ife_read.end())
                {
                    std::vector<int> newrow(2,0);
                    newrow[0] = ife_read[faceid][0];
                    newrow[1] = ife_read[faceid][1];
                    face2elements_update_new[faceid] = newrow;
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
                    ife_send.push_back(ife_read[it->second[u]][0]);
                    ife_send.push_back(ife_read[it->second[u]][1]);
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

            face2elements_update_new[face_id] = ife_loc_row;

            offset = offset + ncol;
        }
    }

}


void RepartitionObject::getFace2EntityPerPartitionRef(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ifref_read, 
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ifref_loc,
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
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    //int Nel = part_global.size();
    int Nel = glob_elem_size;
    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;

    for(itefmap=ief.begin();itefmap!=ief.end();itefmap++)
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
                if(ifref_loc.find(face_req)==ifref_loc.end())
                {
                    ifref_loc[face_req] = ifref_read[face_req];
                }
            }
        }
    }
    
    int own = ifref_loc.size();
    
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
                    int ncol = ifref_read[it->second[u]].size();

                    for(int s=0;s<ncol;s++)
                    {
                        ife_send.push_back(ifref_read[it->second[u]][s]);
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

            ifref_loc[face_id]  = ifref_loc_row;

            offset = offset + ncol;
        }

        ntotal=ntotal+L;
    }
}




                                //     elements2faces_update,
                                //   meshInput->ife,
                                //   meshInput->if_ref,
                                //   meshInput->if_Nv,
                                //   meshInput->ifn,
                                //   Nf_glob,
                                //   face2elements_update,
                                //   face2reference_update,
                                //   face2Nverts_update,
                                //   face2verts_update,
                                //   comm);


void RepartitionObject::getFace2EntityPerPartitionMap(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > ifref_read, 
                                                   std::map<int,std::vector<int> > ifNv_read, 
                                                   std::map<int,std::vector<int> > ifn_read,  
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ife_loc,
                                                   std::map<int,std::vector<int> > &ifref_loc,
                                                   std::map<int,std::vector<int> > &ifn_loc,
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
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    //int Nel = part_global.size();
    int Nel = glob_elem_size;
    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;

    for(itefmap=ief.begin();itefmap!=ief.end();itefmap++)
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
                if(ife_loc.find(face_req)==ife_loc.end())
                {
                    ife_loc[face_req]   = ife_read[face_req];
                    ifref_loc[face_req] = ifref_read[face_req];
                    //ifNv_loc[face_req]  = ifNv_read[face_req];

                    int Nv              = ifNv_read[face_req][0];
                    std::vector<int> ifn_row(Nv,0);
                    for(int j=0;j<Nv;j++)
                    {
                        ifn_row[j] = ifn_read[face_req][j];
                    }

                    ifn_loc[face_req]  = ifn_row;
                }
            }
        }
    }
    
    int own = ife_loc.size();
    
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
                    int Nv   = ifNv_read[it->second[u]][0];
                    int ncol = 2+1+Nv;
                    // ife_send.push_back(ncol);
                    ife_send.push_back(ife_read[it->second[u]][0]);
                    ife_send.push_back(ife_read[it->second[u]][1]);
                    ife_send.push_back(ifref_read[it->second[u]][0]);
                    //ife_send.push_back(ifNv_read[it->second[u]][0]);
                    for(int s=0;s<Nv;s++)
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

            std::vector<int> ife_loc_row(2,0);
            std::vector<int> ifref_loc_row(1,0);
            ife_loc_row[0]   = recv_back_ife[iter->first][offset+0];
            ife_loc_row[1]   = recv_back_ife[iter->first][offset+1];
            ifref_loc_row[0] = recv_back_ife[iter->first][offset+2];
            int Nv           = ncol-3;//recv_back_ife[iter->first][s*ncol+3];
            std::vector<int> ifn_loc_row(Nv,0);
            for(int v=0;v<Nv;v++)
            {
                ifn_loc_row[v] = recv_back_ife[iter->first][offset+3+v];
            }

            ife_loc[face_id]    = ife_loc_row;
            ifref_loc[face_id]  = ifref_loc_row;
            ifn_loc[face_id]    = ifn_loc_row;

            offset = offset + ncol;
        }

        ntotal=ntotal+L;
    }
}






void RepartitionObject::getFace2EntityPerPartitionVector(std::map<int,std::vector<int> > ief, 
                                                   std::map<int,std::vector<int> > ife_read,
                                                   std::map<int,std::vector<int> > ifref_read, 
                                                   std::map<int,std::vector<int> > ifNv_read, 
                                                   std::map<int,std::vector<int> > ifn_read,  
                                                   int Nf_glob, 
                                                   std::map<int,std::vector<int> > &ife_loc,
                                                   std::map<int,std::vector<int> > &ifref_loc,
                                                   std::map<int,std::vector<int> > &ifNv_loc,
                                                   std::map<int,std::vector<int> > &ifn_loc,
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
    
    //std::cout << " " << rank << " LocalVerts.size() before " << LocalVerts.size() << std::endl;
    std::map<int,std::vector<int> > req_face;
    int itel = 0;
    
    //int Nel = part_global.size();
    int Nel = glob_elem_size;
    std::vector<int> ee;
    std::map<int,std::vector<int> >::iterator itefmap;
    std::set<int> collectFids;

    std::vector<int> FIDs;
    std::vector<int> ifeLoc;
    std::vector<int> ifrefLoc;
    std::vector<int> ifNvLoc;
    std::vector<int> ifnLoc;
    for(itefmap=ief.begin();itefmap!=ief.end();itefmap++)
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
                if(collectFids.find(face_req)==collectFids.end())
                {

                    FIDs.push_back(face_req);

                    ifeLoc.push_back(ife_read[face_req][0]);
                    ifeLoc.push_back(ife_read[face_req][1]);

                    ifrefLoc.push_back(ifref_read[face_req][0]);

                    int Nv = ifNv_read[face_req][0];
                    ifNvLoc.push_back(Nv);
                    
                    for(int j=0;j<Nv;j++)
                    {
                        // ifn_row[j] = ifn_read[face_req][j];
                        ifnLoc.push_back(ifn_read[face_req][j]);
                    }
                    // ife_loc[face_req]   = ife_read[face_req];
                    // ifref_loc[face_req] = ifref_read[face_req];
                    // ifNv_loc[face_req]  = ifNv_read[face_req];
                    // int Nv              = ifNv_read[face_req][0];
                    // std::vector<int> ifn_row(Nv,0);
                    
                    // ifn_loc[face_req] = ifn_row;
                }
            }
        }
    }
    
    int own = ife_loc.size();
    
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
                    int Nv   = ifNv_read[it->second[u]][0];
                    int ncol = 2+1+1+Nv;
                    ife_send.push_back(ncol);
                    ife_send.push_back(ife_read[it->second[u]][0]);
                    ife_send.push_back(ife_read[it->second[u]][1]);
                    ife_send.push_back(ifref_read[it->second[u]][0]);
                    ife_send.push_back(ifNv_read[it->second[u]][0]);
                    
                    for(int s=0;s<Nv;s++)
                    {
                        ife_send.push_back(ifn_read[it->second[u]][s]);
                    }

                }

                int nfe_send = ife_send.size();

                int dest = it->first;
                MPI_Send(&nf_send, 1, MPI_INT, dest, 9876*6666+1000*dest, comm);
                MPI_Send(&nfe_send, 1, MPI_INT, dest, 223*6666+1000*dest, comm);

                //MPI_Send(&it->second.data()[0], nf_send, MPI_INT, dest, 9876*7777+dest*888, comm);
                //MPI_Send(&fncol[0], nf_send, MPI_INT, dest, 12*6666+dest*8888, comm);
                //MPI_Send(&ife_send[0], nfe_send, MPI_INT, dest, 9876*6666+dest*8888, comm);

            }
        }
        if(ife_schedule->RecvRankFromRank[q].find( rank ) != ife_schedule->RecvRankFromRank[q].end())
         {
            int nfe_recv;
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876*6666+1000*rank, comm, MPI_STATUS_IGNORE);
            MPI_Recv(&nfe_recv, 1, MPI_INT, q, 223*6666+1000*rank, comm, MPI_STATUS_IGNORE);
             
            std::vector<int> recv_back_ife_arr(nfe_recv);
            //std::vector<int> recv_back_ids_arr(n_recv_back);
            //std::vector<int> fncol_rcv(n_recv_back);

            //MPI_Recv(&recv_back_ids_arr.data()[0], n_recv_back, MPI_INT, q, 9876*7777+rank*888, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&fncol_rcv.data()[0], n_recv_back, MPI_INT, q, 12*6666+rank*8888, comm, MPI_STATUS_IGNORE);
            //MPI_Recv(&recv_back_ife_arr.data()[0], nfe_recv, MPI_INT, q, 9876*6666+rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_Nife[q]       = n_recv_back;
            //recv_back_face_ids[q]   = recv_back_ids_arr;
            //recv_back_face_Ne[q]    = fncol_rcv;
            recv_back_ife[q]        = recv_back_ife_arr;

         }
    }
/*


    std::map<int,int >::iterator iter;
    int ntotal=0;
    ee.clear();
    for(iter=recv_back_Nife.begin();iter!=recv_back_Nife.end();iter++)
    {
        int L = iter->second;
        
        for(int s=0;s<L;s++)
        {
            face_id  = recv_back_face_ids[iter->first][s];
            int ncol = recv_back_face_Ne[iter->first][s];
            // if(ife_loc.find(face_id)==ife_loc.end())
            if(collectFids.find(face_id)==collectFids.end())
            {

                FIDs.push_back(face_id);

                ifeLoc.push_back(recv_back_ife[iter->first][s*ncol+0]);
                ifeLoc.push_back(recv_back_ife[iter->first][s*ncol+1]);

                ifrefLoc.push_back(recv_back_ife[iter->first][s*ncol+2]);

                int Nv = recv_back_ife[iter->first][s*ncol+3];
                ifNvLoc.push_back(Nv);
                
                for(int j=0;j<Nv;j++)
                {
                    // ifn_row[j] = ifn_read[face_req][j];
                    ifnLoc.push_back(recv_back_ife[iter->first][s*ncol+4+j]);
                }


                // std::vector<int> ife_loc_row(2,0);
                // std::vector<int> ifref_loc_row(1,0);
                // std::vector<int> ifNv_loc_row(1,0);
                // ife_loc_row[0]   = recv_back_ife[iter->first][s*ncol+0];
                // ife_loc_row[1]   = recv_back_ife[iter->first][s*ncol+1];
                // ifref_loc_row[0] = recv_back_ife[iter->first][s*ncol+2];
                // ifNv_loc_row[0]  = recv_back_ife[iter->first][s*ncol+3];
                // int Nv = recv_back_ife[iter->first][s*ncol+3];
                // std::vector<int> ifn_loc_row(Nv,0);
                // for(int v=0;v<Nv;v++)
                // {
                //     ifn_loc_row[v] = recv_back_ife[iter->first][s*ncol+4+v];
                // }

                // ife_loc[face_id] = ife_loc_row;
                // ifref_loc[face_id] = ifref_loc_row;
                // ifNv_loc[face_id] = ifNv_loc_row;
                // ifn_loc[face_id] = ifn_loc_row;
            }   
        }

        ntotal=ntotal+L;
    }
    */

}


std::map<int, std::vector<int> > RepartitionObject::getAdjacentElementLayer(std::map<int,std::vector<int> > element2verts,
                                                std::map<int,std::vector<int> > element2faces,
                                                std::map<int,std::vector<int> > element2element,
                                                std::map<int,std::vector<double> > xcn, 
                                                std::map<int,std::vector<double> > data,
                                                int Ne_glob,
                                                int Nf_glob,
                                                int Nv_glob, 
                                                MPI_Comm comm,
                                                std::map<int,std::vector<int> >& elements2verts_update_output,
                                                std::map<int,std::vector<int> >& elements2faces_update_output,
                                                std::map<int,std::vector<int> >& elements2elements_update_output,
                                                std::map<int,std::vector<double> >& elements2data_update_output)
{
    ParallelState* ife_pstate = new ParallelState(Nf_glob,comm);
    ParallelState* xcn_pstate = new ParallelState(Nv_glob,comm);

    std::map<int, std::vector<int> > adjEl2Face;

    std::map<int,std::vector<int> > adj_elements;

    //std::map<int,std::map<int,int> > trace2elements = trace->GetTrace();

    //adj_elements.clear();
    int floc_tmp = 0;
    int vloc_tmp = 0;
    int q=0;
    int i=0;
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    int xcn_o = xcn_pstate->getOffset(rank);

    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int e_id;
    int r;
    std::vector<int> faceIDs_on_rank;
        
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_vert;
    std::map<int,std::vector<int> > rank2req_face;
    int* new_V_offsets = new int[size];
    //int* new_F_offsets = new int[size];

    for(int i=0;i<size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
        //new_F_offsets[i] = ife_pstate->getOffsets()[i]-1;
    }
        
    std::map<int,std::vector<int> > req_elem;
    std::map<int,std::vector<int> > req_trace_elem;
    std::map<int,std::vector<int> > on_trace_elem;

    int itel = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    int ff   = 0;


    std::set<int> treated_set;
    std::set<int> trace_set;


    int ti=0;
    int ti2=0;
    int ti3=0;
    int t = 0;
    for(itmiv=element2faces.begin();itmiv!=element2faces.end();itmiv++)
    {
        int elId     = itmiv->first;
        int nfPerEl  = itmiv->second.size();
        //double varia = LocElemVaria[elId];
        int k        = 0;
        
        for(int j=0;j<nfPerEl;j++)
        {
            int faceid   = itmiv->second[j];
            int adjEl_id = element2element[elId][j];

            if(treated_set.find(faceid)==treated_set.end() && adjEl_id < Ne_glob)
            {
                treated_set.insert(faceid);

                if(part_map.find(adjEl_id)!=part_map.end())
                {
                    p_id = part_map[adjEl_id];
                    //std::cout << "p_id partmap " << p_id << " " << rank << std::endl;

                    if((elem_set.find(adjEl_id)==elem_set.end()))
                    {
                        elem_set.insert(adjEl_id);
                        
                        if(p_id != rank)
                        {
                            adj_elements[p_id].push_back(adjEl_id);
                            req_elem[p_id].push_back(adjEl_id);
                            itel++;
                        }
                    }
                }
                //std::cout << "adjacent2pid " << adjacent2pid.size() << " " << rank << std::endl;
                if(adjacent2pid.find(adjEl_id)!=adjacent2pid.end())
                {

                    int p_id = adjacent2pid[adjEl_id];
                    //std::cout << "p_id " << p_id << " " << rank << std::endl;
                    if(p_id != -1)
                    {
                        if((elem_set.find(adjEl_id)==elem_set.end()))
                        {
                            elem_set.insert(adjEl_id);
                            
                            if(p_id != rank)
                            {
                                adj_elements[p_id].push_back(adjEl_id);
                                req_elem[p_id].push_back(adjEl_id);
                                itel++;
                            }
                            t++;
                        }
                    }
                    
                }

                // if(trace2elements.find(faceid)!=trace2elements.end())
                // {
                //     p_id = trace2elements[faceid][adjEl_id];
                //     //std::cout << "trace2element pid = " << p_id << " rank " << rank  << std::endl;
                //     if(p_id != rank)
                //     {
                //         req_trace_elem[p_id].push_back(adjEl_id);
                //     }
                //     else
                //     {
                //         on_trace_elem[p_id].push_back(adjEl_id);
                //     }
                     
                //     ff++;
                // }
                          
            }
                
            // This condition essentially means a prism is required adjEl_id<Nel && adjEl_id<part_global.size()
            
            
        }
    }


    //std::cout << "rank = " << rank << " " << t << " " <<  adjacent2pid.size() << std::endl;
    

    ScheduleObj* adj_schedule = DoScheduling(req_elem, comm);
    std::map<int,std::vector<int> > reqstd_adj_ids_per_rank;
    std::map<int,std::vector<int> >::iterator it;
    int n_reqstd_adj_ids;


    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = req_elem.begin(); it != req_elem.end(); it++)
            {
                int n_req_adj_el           = it->second.size();
                int dest                   = it->first;
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

            int nvPerEl = LocElem2Nv[adj_id];
            int nfPerEl = LocElem2Nf[adj_id];
            int ndata   = LocElem2Ndata[adj_id];
            //Array<double>* Var  = LocAndAdjElemVaria[adj_id];
            
            send_adj_NvertsPel[dest].push_back(nvPerEl);
            send_adj_NfacesPel[dest].push_back(nfPerEl);
            send_adj_NdatasPel[dest].push_back(ndata);

            for(int k=0;k<ndata;k++)
            {
                double dataV = data[adj_id][k];
                send_adj_data[dest].push_back(dataV);
            }

            //send_adj_VarPel[dest].push_back();
            
            for(int k=0;k<nvPerEl;k++)
            {
                //v_id = globElem2globVerts[adj_id][k];
                v_id = element2verts[adj_id][k];
                send_adj_verts_IDs[dest].push_back(v_id);
            }
        
            for(int k=0;k<nfPerEl;k++)
            {
                f_id = globElem2globFaces[adj_id][k];
                send_adj_faces_IDs[dest].push_back(f_id);
                e_id = globFace2GlobalElements[f_id][k];
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

    std::vector<int> adj_elements_vec;
    std::vector<int> NvPEl_rb;
    std::vector<int> NfPEl_rb;
    std::vector<int> NdataPEl_rb;
    std::vector<double> NVarPEl_rb;
    std::map<int,std::vector<int> >::iterator itm_el;
    int offvvv = 0;

    for(itm_el=adj_elements.begin();itm_el!=adj_elements.end();itm_el++)
    {
        //std::cout << "recv_adj_NvPel " << recv_adj_NvPel[itm_el->first].size() << " " <<  recv_adj_NvPel[itm_el->first].size() << " " << recv_adj_NfPel[itm_el->first].size() << " " << recv_adj_NVarPel[itm_el->first].size() << std::endl;
        //TotNrho_adj_recv = TotNrho_adj_recv+itm->second;
        for(int i=0;i<itm_el->second.size();i++)
        {
            adj_elements_vec.push_back(adj_elements[itm_el->first][i]);
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
    
    //std::cout << "TotNelem_adj_recv " << TotNelem_adj_recv << " should be equal to " << TotNrho_adj_recv << std::endl;
    //std::cout << " Compare " <<  TotNelem_adj_recv << " " << itel << " " << adj_rhos.size() << std::endl;

    for(int i=0;i<adj_verts.size();i++)
    {
        int v_id_n = adj_verts[i];
        r = FindRank(new_V_offsets,size,v_id_n);

        if(unique_vertIDs_on_rank_set.find( v_id_n ) == unique_vertIDs_on_rank_set.end()) // add the required unique vertex for current rank.
        {
            unique_vertIDs_on_rank_set.insert(v_id_n);
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

        if(unique_faceIDs_on_rank_set.find( f_id_n ) == unique_faceIDs_on_rank_set.end()) // add the required unique vertex for current rank.
        {
            unique_faceIDs_on_rank_set.insert(f_id_n);
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

           if(o_gvertex2lvertex.find(gvid)==o_gvertex2lvertex.end())
           {
               
               std::vector<double> V(3);

               // V[0] = xcn->getVal(gvid-xcn_o,0);
               // V[1] = xcn->getVal(gvid-xcn_o,1);
               // V[2] = xcn->getVal(gvid-xcn_o,2);

               V[0] = xcn[gvid][0];
               V[1] = xcn[gvid][1];
               V[2] = xcn[gvid][2];
               
            //    o_lvertex2gvertex[lvid] = gvid;
            //    o_gvertex2lvertex[gvid] = lvid;
               //LocalVerts.push_back(V);
               LocalVertsMap[gvid] = V;

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
               
               if(o_gvertex2lvertex.find(gvid)==o_gvertex2lvertex.end())
               {
                   
                   std::vector<double> V(3);
                   V[0] = it_f->second[u*3+0];
                   V[1] = it_f->second[u*3+1];
                   V[2] = it_f->second[u*3+2];

                   //LocalVerts.push_back(V);
                   LocalVertsMap[gvid] = V;
                //    o_lvertex2gvertex[lvid]=gvid;
                //    o_gvertex2lvertex[gvid]=lvid;
                   
                   lvid++;
               }
            }
       }

       nLoc_Verts = LocalVertsMap.size();
    //std::cout << " " << rank << " LocalVerts.size() after " << LocalVerts.size() << std::endl;

    //NlocElem = eloc;
    //std::cout << "eloc " << eloc << std::endl;
    // ================================== Faces on Rank =========================================

    int lfid = floc;
    int gfid = 0;
    for(int m=0;m<floc_tmp;m++)
    {
        gfid = faceIDs_on_rank[m];
        lface2gface[lfid] = gfid;
        gface2lface[gfid] = lfid;
        lfid++;
    }

    int cnv = 0;
    int cnf = 0;
    int idsave = 0;
    //double rho_v;
    int loc_v;
    int glob_f;
    int loc_f;
    int glob_v;
    int glob_e=0;
    //std::cout << adj_verts.size() << " " << Nel_extra2 <<std::endl;


    

    
    
    double varia_v = 0.0;
    std::set<int> LocAdjElemSet;
    std::vector<int> adjElLayer(3*itel);
    int offvv = 0;
    for(int m=0;m<adj_elements_vec.size();m++)
    {
        el_id         = adj_elements_vec[m];
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
        
        LocalElement2GlobalElement[eloc] = el_id;
        GlobalElement2LocalElement[el_id] = eloc;
        eloc++;
        //rho_v = adj_rhos[m];

        //U0Elem.push_back(rho_v);
        // Array<double>* VariaV_arr = new Array<double>(1,1);
        // VariaV_arr->setVal(0,0,VariaV);
        //LocAndAdjElemVaria[el_id] = VariaV_arr;
        LocAdjElemSet.insert(el_id);
        LocAndAdj_Elem.push_back(el_id);
        LocAndAdj_Elem_Nv.push_back(Nv);
        LocAndAdj_Elem_Nf.push_back(Nf);
        adjElLayer[3*m+0] = el_id;
        adjElLayer[3*m+1] = Nv;
        adjElLayer[3*m+2] = Nf;
        
        //LocAndAdj_Elem_Varia.push_back(varia_v);
//      U0Elem->setVal(m+o,0,rho_v);
//      ElemPart->setVal(m+o,0,el_id);

        for(int p=0;p<Nv;p++)
        {
            glob_v = adj_verts[offvv+p];
            loc_v  = o_gvertex2lvertex[glob_v];
            //LocalElem2GlobalVert->setVal(m+o,p,glob_v);
            //lelement2lvertex->setVal(m+o,p,loc_v);
            tmp_globv[p] = glob_v;
            tmp_locv[p] = loc_v;
            globElem2globVerts[el_id].push_back(glob_v);
            globVerts2globElem[glob_v].push_back(el_id);

            globElem2locVerts[el_id].push_back(loc_v);
            //collect_var[loc_v].push_back(rho_v);
            cnv++;
            
        }
        for(int p=0;p<Nf;p++)
        {
            glob_f = adj_faces[cnf];
            glob_e = adj_element[cnf];
            tmp_globf[p] = glob_f;
            tmp_globe[p] = glob_e;
            loc_f  = gface2lface[glob_f];
            globElem2localFaces[el_id].push_back(loc_f);
            globElem2globFaces[el_id].push_back(glob_f);
            globFace2GlobalElements[glob_f].push_back(el_id);
            cnf++;
        }
        
        //=======================================
        elements2verts_update_output[el_id]    = tmp_globv;
        elements2faces_update_output[el_id]    = tmp_globf;
        elements2elements_update_output[el_id] = tmp_globe;
        elements2data_update_output[el_id]     = datarow;
        //=======================================
        adjEl2Face[el_id] = tmp_globf;

        offvv = offvv+Nv;
        //LocalElem2GlobalVert.push_back(tmp_globv);
        lelement2lvertex.push_back(tmp_locv);
        //tmp_globv.clear();
        tmp_locv.clear();
    }

    delete[] new_V_offsets;
    //delete[] new_F_offsets;

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







void RepartitionObject::AddStateVecForAdjacentElements(std::map<int,std::vector<double> > &U, int nvar, MPI_Comm comm)
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
    //std::cout << xcn->getOffset(rank) << " " << xcn_pstate->getOffset(rank) << std::endl;
    
    //int ien_o = ien_pstate->getOffset(rank);

    int el_id;
    int p_id;
    int v_id;
    int f_id;
    int r;
    std::vector<int> faceIDs_on_rank;
    std::vector<int> vertIDs_on_rank;
    std::map<int,std::vector<int> > rank2req_Elems;
    std::map<int,std::vector<double> > U_loc;
    
    
    std::map<int,std::vector<int> > req_elem;
    int itel = 0;
    
    std::vector<int> ee;

    for(int i=0;i<LocAndAdj_Elem.size();i++)
    {
        int el_req = LocAndAdj_Elem[i];

        if(part_map.find(el_req)!=part_map.end())
        {
            r = part_map[el_req];
        }

        if(adjacent2pid.find(el_req)!=adjacent2pid.end())
        {
            if(adjacent2pid[el_req]!=-1)
            {
                r = adjacent2pid[el_req];
            }
        }

        if(r != rank)
        {
            rank2req_Elems[r].push_back(el_req);
        }
//        else
//        {
//            Array<double>* StateVec = new Array<double>(nvar,1);
//            for(int q=0;q<nvar;q++)
//            {
//                StateVec->setVal(q,0,U[el_req]->getVal(q,0));
//            }
//            U_loc[el_req] = StateVec;
//        }
    }
    
    ScheduleObj* iee_schedule = DoScheduling(rank2req_Elems,comm);

    std::map<int,std::vector<int> >::iterator it;
    std::map<int,std::vector<int> >  reqstd_E_IDs_per_rank;

    for(q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (it = rank2req_Elems.begin(); it != rank2req_Elems.end(); it++)
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
                //double* iee_send      = new double[ne_send*nvar];
                std::vector<double> iee_send(ne_send*nvar);
                //int offset_iee        = ien_pstate->getOffset(rank);
                
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







void RepartitionObject::buildCommunicationMaps(MPI_Comm comm)
{

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    int element, fid;
    int Nel_loc = Loc_Elem.size();
    DistributedParallelState* newElemDist = new DistributedParallelState(Nel_loc,comm);
    
    std::map<int,std::vector<int> >::iterator itf;
    std::vector<int> sharedFonRank;
    std::vector<int> interiorFonRank;
    //std::vector<std::vector<double> > locVs = LocalVerts;

    for(itf=o_face2elements_global.begin();itf!=o_face2elements_global.end();itf++)
    {
        if(itf->second.size()==1)
        {
            sharedFonRank.push_back(itf->first);
        }
        if(itf->second.size()==2)
        {
            interiorFonRank.push_back(itf->first);
        }
    }
    
    int nSharedFonRank = sharedFonRank.size();
    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFonRank,comm);
    // std::cout << "sharedFonRank " << sharedFonRank.size() << std::endl;
    // std::cout << "interiorFonRank " << interiorFonRank.size() << std::endl;

    int Nt_shFaces               = distSharedFaces->getNel();
    int* shFace_offsets          = distSharedFaces->getOffsets();
    int* shFace_nlocs            = distSharedFaces->getNlocs();
    std::vector<int> shFacesIDs(nSharedFonRank,0);
    std::vector<int> shFaces_RankIDs(nSharedFonRank,0);

    for(int i=0;i<nSharedFonRank;i++)
    {
        shFacesIDs[i]      = sharedFonRank[i];
        shFaces_RankIDs[i] = world_rank;
    }
    
    std::vector<int> TotalSharedFaces(Nt_shFaces,0);
    std::vector<int> TotalSharedFaces_RankID(Nt_shFaces,0);
    
    // Communicate face map to all ranks.
    MPI_Allgatherv(&shFacesIDs.data()[0],
                   nSharedFonRank,
                   MPI_INT,
                   &TotalSharedFaces.data()[0],
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, comm);
    
    MPI_Allgatherv(&shFaces_RankIDs.data()[0],
                   nSharedFonRank,
                   MPI_INT,
                   &TotalSharedFaces_RankID.data()[0],
                   shFace_nlocs,
                   shFace_offsets,
                   MPI_INT, comm);
    
    std::map<int,std::vector<int> > face2rank;
    
    for(int i=0;i<Nt_shFaces;i++)
    {
        int key = TotalSharedFaces[i];
        int val = TotalSharedFaces_RankID[i];
        face2rank[key].push_back(val);
    }
    
    //delete[] TotalSharedFaces;
    //delete[] TotalSharedFaces_RankID;
    delete[] shFace_offsets;
    delete[] shFace_nlocs;
    //delete[] shFacesIDs;
    //delete[] shFaces_RankIDs;
    //delete distSharedFaces;
    
    std::map<int,std::vector<int> >::iterator itff;
    int shf = 0;
    int bf  = 0;
    std::map<int,std::vector<int> > Boundary_Ref2Face;
    std::set<int> uSharedVert;
    std::set<int> uBoundVert;
    int f = 0;

    // std::map<int,int> o_globShF2locShF;
    
    std::vector<int> m_faces4parmmg;

    for(itff=face2rank.begin();itff!=face2rank.end();itff++)
    {
        int faceID = itff->first;
        
        if(itff->second.size()==2)
        {
            if(itff->second[0]==world_rank)
            {
                o_ColorsFaces[itff->second[1]].push_back(faceID);
            }
            else if(itff->second[1]==world_rank)
            {
                o_ColorsFaces[itff->second[0]].push_back(faceID);
            }

            if(o_face2verts_global.find(faceID)!=o_face2verts_global.end())
            {
                o_faces4parmmg.push_back(faceID);
                o_globShF2locShF[faceID] = f;
                o_locShF2globShF[f] = faceID;
                f++;
                
                shf++;
            }
            
            
        }
        
        if(itff->second.size()==1)
        {         
             if(o_face2verts_global.find(faceID)!=o_face2verts_global.end())
            {
                o_faces4parmmg.push_back(faceID);
                o_globShF2locShF[faceID] = f;
                o_locShF2globShF[f] = faceID;

                f++;
            }   
            
        }
    }
    
    o_ncomm           = o_ColorsFaces.size();
    o_color_face = (int *) malloc(o_ncomm*sizeof(int));
    o_ntifc = (int *) malloc(o_ncomm*sizeof(int));
    
    o_ifc_tria_loc = (int **)malloc(o_ncomm*sizeof(int *));
    o_ifc_tria_glo = (int **)malloc(o_ncomm*sizeof(int *));

    int icomm=0;
    std::map<int,std::vector<int> >::iterator itc;
    
    
    for(itc=o_ColorsFaces.begin();itc!=o_ColorsFaces.end();itc++)
    {
        o_color_face[icomm]     = itc->first;
        o_ntifc[icomm]          = itc->second.size();
        o_ifc_tria_loc[icomm]   = (int *) malloc(itc->second.size()*sizeof(int));
        o_ifc_tria_glo[icomm]  = (int *) malloc(itc->second.size()*sizeof(int));

        for(int q=0;q<itc->second.size();q++)
        {
            o_ifc_tria_glo[icomm][q] = itc->second[q]+1;
            o_ifc_tria_loc[icomm][q]  = o_globShF2locShF[itc->second[q]]+1;
        }

        //std::cout << std::endl;


        icomm++;
    }
}






void RepartitionObject::buildUpdatedVertexAndFaceNumbering(MPI_Comm comm, 
                                                           std::map<int,std::vector<int> > ranges_id,
                                                           std::map<int,std::vector<int> > ranges_ref)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);

    //std::map<int,std::map<int,int> > trace_elem     = trace->GetTrace();
    //std::map<int,std::vector<int> > trace_verts     = trace->GetTraceVerts();
    //std::map<int,int> uniqure_trace_verts2ref       = trace->GetUniqueTraceVerts2RefMap();

    std::set<int> m_loc_trace_verts_glob = AllGatherSet(m_loc_trace_verts,comm);
    int Nel_loc  = Loc_Elem.size();
    int ref;
    
    std::set<int> gvid_set;
    int inf = 0;
    int shf = 0;
    int tracef = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    std::map<int,int> sharedFace2Rank;
    std::set<int> unique_vertex_set_on_rank;
    int bface = 0;

    std::set<int> refsc;

    std::map<int,int> sharedFaces;
    std::map<int,int> interiorFaces;
    std::map<int,int> boundaryFaces;

    std::set<int>::iterator its;
    int ntrace = 0;

    for(its=Loc_Elem_Set.begin();its!=Loc_Elem_Set.end();its++)
    {
        int elid = *its;
        int Nf = elements2faces_update[elid].size();

        for(int j=0;j<Nf;j++)
        {
            int gfid = elements2faces_update[elid][j];

            if(m_loc_trace_faces.find(gfid)!=m_loc_trace_faces.end())
            {
                ref = 13;
                ntrace++;
            }
            else
            {
                if(face2verts_update.find(gfid)!=face2verts_update.end())
                {
                    // ref      = face2reference_update[gfid][0];
                    ref      = ProvideBoundaryRef(gfid,ranges_ref);
                    int e0   = face2elements_update[gfid][0];
                    int e1   = face2elements_update[gfid][1];
                    int Nv   = face2verts_update[gfid].size();

                    // copy the correct amount of nodes into a new face vector.
                    std::vector<int> face(Nv,0);
                    for(int k=0;k<Nv;k++)
                    {
                        int vid = face2verts_update[gfid][k];
                        face[k] = vid;

                        if(unique_vertex_set_on_rank.find(vid)==unique_vertex_set_on_rank.end() &&
                        m_loc_trace_verts_glob.find(vid)==m_loc_trace_verts_glob.end())
                        {
                            unique_vertex_set_on_rank.insert(vid);
                        }
                    }
                    if(e0<Ne_glob && e1<Ne_glob)
                    {

                        int r0;// = part_global[e0]; // rank of first adjacent element.
                        int r1;// = part_global[e1]; // rank of second adjacent element.

                        if(part_map.find(e0)!=part_map.end())
                        {
                            r0 = part_map[e0];
                        }
                        if(adjacent2pid.find(e0)!=adjacent2pid.end())
                        {
                            r0 = adjacent2pid[e0];
                        }

                        if(part_map.find(e1)!=part_map.end())
                        {
                            r1 = part_map[e1];
                        }
                        if(adjacent2pid.find(e1)!=adjacent2pid.end())
                        {
                            r1 = adjacent2pid[e1];
                        }
                        // int r0 = part_global[e0]; // rank of first adjacent element.
                        // int r1 = part_global[e1]; // rank of second adjacent element.


                        if(ref==2)// Internal and Shared faces are here.
                        {
                            if(r0==rank && r1!=rank)
                            {
                                sharedFace2Rank[gfid] = r0;
                            }
                            if(r0!=rank && r1==rank)
                            {
                                sharedFace2Rank[gfid] = r1;
                            }
                        }
                    }
                }
            }
        }
    }

    


    // Take care of shared faces and vertices
    std::set<int> UniqueSharedVertsOnRank_set;
    std::vector<int> UniqueSharedVertsOnRank_RankID;

    std::set<int> UniqueSharedVertsOnRank_set_ref;
    std::vector<int> UniqueSharedVertsOnRank_RankID_ref;

    std::map<int,int>::iterator itmii;
    for(itmii=sharedFace2Rank.begin();itmii!=sharedFace2Rank.end();itmii++)
    {
        int gfid = itmii->first;
        // int Nv = face2verts_update[gfid].size();
        int Nv   = face2verts_update[gfid].size();


        if(gfid==2353748)
        {
            std::cout << "2353748 on rank queries" << std::endl;
        }

        for(int q=0;q<Nv;q++)
        {
            int vid = face2verts_update[gfid][q];
            //std::cout << vid << " ";

            if(vid==807212)
            {
                std::cout << "found here " << gfid << " " << rank << " " << Nv << std::endl;
            }
            if(UniqueSharedVertsOnRank_set.find(vid)==UniqueSharedVertsOnRank_set.end() &&
                   m_loc_trace_verts_glob.find(vid)==m_loc_trace_verts_glob.end()
            )
            {
                UniqueSharedVertsOnRank_set.insert(vid);
                UniqueSharedVertsOnRank_RankID.push_back(rank);
            }

            if(UniqueSharedVertsOnRank_set_ref.find(vid)==UniqueSharedVertsOnRank_set_ref.end() &&
                   m_loc_trace_verts_glob.find(vid)==m_loc_trace_verts_glob.end()
            )
            {
                UniqueSharedVertsOnRank_set_ref.insert(vid);
                UniqueSharedVertsOnRank_RankID_ref.push_back(rank);
            }
        }
        //std::cout << std::endl;
    }


    // if(UniqueSharedVertsOnRank_set_ref.size()!=UniqueSharedVertsOnRank_set.size())
    // {   
    //     int notf    = 0;
    //     int fo      = 0;
    //     std::set<int>::iterator its;
    //     for(its=UniqueSharedVertsOnRank_set.begin();its!=UniqueSharedVertsOnRank_set.end();its++)
    //     {
    //         int vid = *its;
    //         if(UniqueSharedVertsOnRank_set_ref.find(vid)==UniqueSharedVertsOnRank_set_ref.end())
    //         {
    //             std::cout << "rank " << rank << " vid " << vid << std::endl; 
    //             notf;
    //         }
    //         else
    //         {
    //             fo++;
    //         }
    //     }
    //     std::cout << "UniqueSharedVertsOnRank_set_ref " << UniqueSharedVertsOnRank_set_ref.size() << " " << UniqueSharedVertsOnRank_set.size() << " " << notf<< " " << fo << std::endl;
    // }


    int nLocalVertsTot = unique_vertex_set_on_rank.size();
    std::vector<int> nLocalVertsTot_loc_test(size,0);
    nLocalVertsTot_loc_test[rank] = nLocalVertsTot;
    std::vector<int> nLocalVertsTot_red_test(size,0);

    MPI_Allreduce(nLocalVertsTot_loc_test.data(), 
                  nLocalVertsTot_red_test.data(), 
                  size, MPI_INT, MPI_SUM, comm);

    int N_un_shareV = UniqueSharedVertsOnRank_set.size();
    int summednLocalVertsTot = 0;

    for(int i=0;i<size;i++)
    {
        summednLocalVertsTot = summednLocalVertsTot+nLocalVertsTot_red_test[i];
    }
    int testNv = summednLocalVertsTot - N_un_shareV;
    
    //=====================================================================================
    //================================SHARED FACES=========================================
    //=====================================================================================
    int nSharedFaces = sharedFace2Rank.size();
    int N_localFaces = face2elements_update.size();
    DistributedParallelState* distSharedFaces = new DistributedParallelState(nSharedFaces,comm);
    DistributedParallelState* distLocalFaces  = new DistributedParallelState(N_localFaces,comm);

    int Nt_shFaces   = distSharedFaces->getNel();
    std::vector<int> sharedFaceIDs_tmp(sharedFace2Rank.size(),0);
    std::vector<int> sharedFaceRankIDs_tmp(sharedFace2Rank.size(),0);
    
    int c = 0;

    // Copy the face ID in the sharedface2mode map into temporary vector in order to communicate 
    
    for(itmii=sharedFace2Rank.begin();itmii!=sharedFace2Rank.end();itmii++)
    {
        sharedFaceIDs_tmp[c]     = itmii->first;
        sharedFaceRankIDs_tmp[c] = itmii->second;
        c++;
    }

    std::vector<int> TotalSharedFaces(Nt_shFaces,0);
    std::vector<int> TotalSharedFaces_RankID(Nt_shFaces,0);
    // Communicate face map to all ranks.
    MPI_Allgatherv(&sharedFaceIDs_tmp.data()[0],
                   nSharedFaces,
                   MPI_INT,
                   &TotalSharedFaces.data()[0],
                   distSharedFaces->getNlocs(),
                   distSharedFaces->getOffsets(),
                   MPI_INT, comm);
    
    MPI_Allgatherv(&sharedFaceRankIDs_tmp.data()[0],
                   nSharedFaces,
                   MPI_INT,
                  &TotalSharedFaces_RankID.data()[0],
                   distSharedFaces->getNlocs(),
                   distSharedFaces->getOffsets(),
                   MPI_INT, comm);

    int tmp;
    
    
    for(int i=0;i<Nt_shFaces;i++)
    {
        int key = TotalSharedFaces[i];
        int val = TotalSharedFaces_RankID[i];

        o_SharedFace2Rank[key].push_back(val);

        if(o_sharedFace2RankMap.find(key)==o_sharedFace2RankMap.end())
        {
            o_sharedFace2RankMap[key]=val;
        }
        else
        {
            tmp = o_sharedFace2RankMap[key];
            if(val<tmp)
            {
                o_sharedFace2RankMap[key]=val;
            }
            if(val>tmp)
            {
                o_sharedFace2RankMap[key]=tmp;
            }
        }
    }

    std::map<int,int> sharedFmap;
    int iFshared = distLocalFaces->getNel()-o_sharedFace2RankMap.size();
    std::map<int,int>::iterator itii;
    // for(itii=o_sharedFace2RankMap.begin();
    //         itii!=o_sharedFace2RankMap.end();itii++)
    // {
    //     sharedFmap[itii->first] = iFshared;
    //     iFshared++;
    // }

    
    //=====================================================================================
    //===============================SHARED VERTEX=========================================
    //=====================================================================================
    int nSharedVerts = UniqueSharedVertsOnRank_set.size();
    DistributedParallelState* distSharedVerts = new DistributedParallelState(nSharedVerts,comm);
    
    int Nt_shVerts = distSharedVerts->getNel();
    std::vector<int> TotalSharedVerts(Nt_shVerts,0);
    std::vector<int> TotalSharedVerts_RankID(Nt_shVerts,0);

    std::vector<int> UniqueSharedVertsOnRank_vec(UniqueSharedVertsOnRank_set.begin(),
                                                 UniqueSharedVertsOnRank_set.end());


    MPI_Allgatherv(&UniqueSharedVertsOnRank_vec[0],
                    nSharedVerts,
                    MPI_INT,
                   &TotalSharedVerts[0],
                    distSharedVerts->getNlocs(),
                    distSharedVerts->getOffsets(),
                    MPI_INT, comm);

    MPI_Allgatherv(&UniqueSharedVertsOnRank_RankID[0],
                    nSharedVerts,
                    MPI_INT,
                   &TotalSharedVerts_RankID[0],
                    distSharedVerts->getNlocs(),
                    distSharedVerts->getOffsets(),
                    MPI_INT, comm);


     
    for(int i=0;i<Nt_shVerts;i++)
    {
        int key = TotalSharedVerts[i];
        int val = TotalSharedVerts_RankID[i];
        
        if(o_sharedVertex2RankMap.find(key)==o_sharedVertex2RankMap.end())
        {
            o_sharedVertex2RankMap[key]=val;
        }
        else
        {
            tmp = o_sharedVertex2RankMap[key];
            
            if(val<tmp)
            {
                o_sharedVertex2RankMap[key]=val;
            }
            if(val>tmp)
            {
                o_sharedVertex2RankMap[key]=tmp;
            }
        }
    }

    std::map<int,int>::iterator itmm;
    int nOwnedSharedVerts = 0;
    for(itmm=o_sharedVertex2RankMap.begin();itmm!=o_sharedVertex2RankMap.end();itmm++)
    {
        if(itmm->second==rank)
        {
            nOwnedSharedVerts++;
        }
    }
    
    DistributedParallelState* ownedSharedVrtsDist = new DistributedParallelState(nOwnedSharedVerts,comm);

    int nNonSharedVerts          = nLocalVertsTot-nSharedVerts;
    DistributedParallelState* distLocalVerts = new DistributedParallelState(nLocalVertsTot,comm);

    //int nNonSharedVerts          = nLocalVerts-nSharedVerts;
    //int nNonSharedFaces          = nLocalFaces-nSharedFaces;
    int nNonSharedFaces          = N_localFaces-nSharedFaces;

    DistributedParallelState* nonSharedVertDistr = new DistributedParallelState(nNonSharedVerts,comm);
    DistributedParallelState* nonSharedFaceDistr = new DistributedParallelState(nNonSharedFaces,comm);


    int iVshared = distLocalVerts->getNel()-o_sharedVertex2RankMap.size();

    std::map<int,int >::iterator itvv;
    std::map<int,int> sharedVmap;
    for(itvv=o_sharedVertex2RankMap.begin();itvv!=o_sharedVertex2RankMap.end();itvv++)
    {
        // std::cout << "iVshared "  << itvv->first << "  --  " <<  iVshared << " " << distLocalVerts->getNel() << " " << o_sharedVertex2RankMap.size()<< std::endl;
        sharedVmap[itvv->first] = iVshared;
        iVshared++;
    }


    std::vector<int> ownedVs(size,0);

    for(int i=0;i<size;i++)
    {
        ownedVs[i] = 0;
    }

    int nNonSharedVertsTot = nonSharedVertDistr->getNel();

    
    for(itmm=o_sharedVertex2RankMap.begin();itmm!=o_sharedVertex2RankMap.end();itmm++)
    {
        int globid = nNonSharedVertsTot+ownedSharedVrtsDist->getOffsets()[itmm->second]+ownedVs[itmm->second]+1;
        
        o_sharedVertexMapUpdatedGlobalID[itmm->first] = globid;
        
        if(itmm->second==rank &&
           m_loc_trace_verts_glob.find(itmm->first)==m_loc_trace_verts_glob.end())
        {
            o_SharedVertsOwned[itmm->first] = globid;
        }
        
        ownedVs[itmm->second] = ownedVs[itmm->second]+1;
    }

    
    int Fid_shared = distLocalFaces->getNel()-o_sharedFace2RankMap.size();

    for(itii=o_sharedFace2RankMap.begin();itii!=o_sharedFace2RankMap.end();itii++)
    {
        o_sharedFaceMapUpdatedGlobalID[itii->first] = Fid_shared;
        Fid_shared++;
    }    


    DistributedParallelState* ElementDistr = new DistributedParallelState(part.size(),comm);
    int u = 0;
    int gvidd = 0;

    std::map<int,int> tag2ElementID;
    std::set<int> faceset;
    // int N_localFaces = face2elements_update.size();
    // int nNonSharedFaces = N_localFaces-nSharedFaces;
    std::vector<int> nNonSharedFacesArray(size,0);
    std::vector<int> nNonSharedFacesArrayRed(size,0);
    std::vector<int> nNonSharedArrayRed(size,0);
    std::vector<int> nNonSharedVertsArrayOff(size,0);
    for(int i=0;i<size;i++)
    {
        nNonSharedFacesArray[i] = 0;
        if(i==rank)
        {
            nNonSharedFacesArray[i] = nNonSharedFaces;
        }
    }

    MPI_Allreduce(nNonSharedFacesArray.data(),
                nNonSharedFacesArrayRed.data(),
                size,
                MPI_INT, MPI_SUM, comm);

    std::vector<int> nNonSharedFacesArrayOff(size,0);
    int nonFacesSharedOff     = 0;
    int nonSharedOff = 0;
    for(int i=0;i<size;i++)
    {
        nNonSharedVertsArrayOff[i] = nonSharedOff;
        nNonSharedFacesArrayOff[i] = nonFacesSharedOff;
        
        nonSharedOff = nonSharedOff + nNonSharedArrayRed[i];
        nonFacesSharedOff = nonFacesSharedOff + nNonSharedFacesArrayRed[i];
    }

    int lfid = nNonSharedFacesArrayOff[rank];

    // int Nel_loc = Loc_Elem.size();
    
    int gloVid = nNonSharedFacesArrayOff[rank];
    int gloFid = nNonSharedFacesArrayOff[rank];
    int locVid = 0;
    int locFid = 0;
    
    for(int i = 0;i < Nel_loc;i++)
    {
        int gelid   = Loc_Elem[i];       
        int gEl     = ElementDistr->getOffsets()[rank]+u+1;
        int Nf      = elements2faces_update[gelid].size();
        int Nv      = elements2verts_update[gelid].size();

        // std::vector<int> vertices(Nv,0);

        for(int p=0;p<Nv;p++)
        {
            int TagVid = elements2verts_update[gelid][p];

            if(o_sharedVertexMapUpdatedGlobalID.find(TagVid)!=o_sharedVertexMapUpdatedGlobalID.end())
            {
                int GlobVID             = o_sharedVertexMapUpdatedGlobalID[TagVid];

                if(o_tag2globV.find(TagVid)==o_tag2globV.end())
                {
                    o_tag2globV[TagVid] = GlobVID;
                    //o_glob2tagV[GlobVID] = TagVid;
                    //o_loc2globV[locVid] = GlobVID;
                    o_glob2locV[GlobVID] = locVid;
                    locVid = locVid + 1;
                }
                             
            }
            else
            {
                if(o_tag2globV.find(TagVid)==o_tag2globV.end())
                {
                    o_tag2globV[TagVid] = gloVid;
                    //o_glob2tagV[gloVid] = TagVid;
                    //o_loc2globV[locVid] = gloVid;
                    o_glob2locV[gloVid] = locVid;
                    gloVid = gloVid + 1;
                    locVid = locVid + 1;
                }
                
            }         
        }

        // o_element2verts_global[gEl] = vertices;

        std::vector<int> faces(Nf,0);
        for(int q=0;q<Nf;q++)
        {
            int TagFid = elements2faces_update[gelid][q];

            if(o_sharedFaceMapUpdatedGlobalID.find(TagFid)!=o_sharedFaceMapUpdatedGlobalID.end())
            {
                int GlobFID = o_sharedFaceMapUpdatedGlobalID[TagFid]; 
                faces[q]    = GlobFID; 
                o_face2elements_global[GlobFID].push_back(gEl);

                if(o_tag2globF.find(TagFid)==o_tag2globF.end())
                {
                    o_tag2globF[TagFid] = GlobFID;
                    o_glob2tagF[GlobFID] = TagFid;
                    // o_loc2globF[locFid] = GlobFID;
                    // o_glob2locF[GlobFID] = locFid;
                }

            }
            else
            {
                if(o_tag2globF.find(TagFid)==o_tag2globF.end())
                {
                    o_tag2globF[TagFid] = gloFid;
                    o_glob2tagF[gloFid] = TagFid;
                    // o_loc2globF[locFid] = gloFid;
                    // o_glob2locF[gloFid] = locFid;

                    faces[q]     = gloFid;
                    o_face2elements_global[gloFid].push_back(gEl);
                    gloFid       = gloFid + 1;
                }
                else
                {
                    int gloFid_tmp = o_tag2globF[TagFid];
                    faces[q]       = gloFid_tmp;

                    o_face2elements_global[gloFid_tmp].push_back(gEl);
                }

                lfid++;
            }
        }

        o_element2faces_global[gEl] = faces;

        u = u + 1;

    }
    
    for(itmiv=o_element2faces_global.begin();itmiv!=o_element2faces_global.end();itmiv++)
    {
        int eid = itmiv->first;
        int nf = itmiv->second.size();
        
        for(int i=0;i<nf;i++)
        {
            int fid    = o_element2faces_global[eid][i];
            int tagfid = o_glob2tagF[fid];
            // int nv     = face2verts_update[tagfid].size();


            if(face2verts_update.find(tagfid)!=face2verts_update.end())
            {
                int nv     = face2verts_update[tagfid].size();
                std::vector<int> verts(nv,0);
                for(int j=0;j<nv;j++)
                {
                    int tagvid  =  face2verts_update[tagfid][j];
                    verts[j]    =  o_tag2globV[tagvid];
                }
                if(o_face2verts_global.find(fid)==o_face2verts_global.end())
                {
                    o_face2verts_global[fid] = verts;
                }
            }
        }
    }

    m_loc_trace_verts_glob.clear();
    /**/
}









void RepartitionObject::buildInteriorSharedAndBoundaryFaceMaps(MPI_Comm comm, 
                                                               std::map<int,std::vector<int> > ranges_id,
                                                               std::map<int,std::vector<int> > ranges_ref)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    std::map<int,std::vector<int> >::iterator itmiv;
    //==============================
    // BEGIN building o_tagE2gE
    // BEGIN building o_gE2tagE
    // BEGIN building LeftHandRightHandElementVec lhp/rhp
    //==============================

    // buildInteriorSharedAndBoundaryFaceMaps()
    // std::map<int,std::map<int,int> > trace_elem     = trace->GetTrace();
    // //std::map<int,std::vector<int> > trace_verts     = trace->GetTraceVerts();
    // std::map<int,int> uniqure_trace_verts2ref       = trace->GetUniqueTraceVerts2RefMap();
    
    DistributedParallelState* ElementDistr = new DistributedParallelState(Loc_Elem_Set.size(),comm);

    std::set<int> m_loc_trace_verts_glob = AllGatherSet(m_loc_trace_verts, comm);

    std::map<int,int> tag2ElementID;
    std::set<int> faceset;
    int nSharedFaces = o_sharedFace2Nodes.size();
    //std::cout << "nSharedFaces from long " << nSharedFaces << std::endl;
    int N_localFaces = face2elements_update.size();
    int nNonSharedFaces = N_localFaces-nSharedFaces;
    std::vector<int> nNonSharedFacesArray(size,0);
    std::vector<int> nNonSharedFacesArrayRed(size,0);
    std::vector<int> nNonSharedArrayRed(size,0);
    std::vector<int> nNonSharedVertsArrayOff(size,0);
    for(int i=0;i<size;i++)
    {
        nNonSharedFacesArray[i] = 0;
        if(i==rank)
        {
            nNonSharedFacesArray[i] = nNonSharedFaces;
        }
    }

    MPI_Allreduce(nNonSharedFacesArray.data(),
                nNonSharedFacesArrayRed.data(),
                size,
                MPI_INT, MPI_SUM, comm);

    std::vector<int> nNonSharedFacesArrayOff(size,0);
    int nonFacesSharedOff     = 0;
    int nonSharedOff = 0;
    for(int i=0;i<size;i++)
    {
        nNonSharedVertsArrayOff[i] = nonSharedOff;
        nNonSharedFacesArrayOff[i] = nonFacesSharedOff;
        
        nonSharedOff = nonSharedOff + nNonSharedArrayRed[i];
        nonFacesSharedOff = nonFacesSharedOff + nNonSharedFacesArrayRed[i];
    }

    int lvid  = nNonSharedFacesArrayOff[rank];
    int lfid  = nNonSharedFacesArrayOff[rank];
    // lfid = 0;
    std::map<int,int> tagF2locFID;
    int fownedInmap = 0;
    std::set<int>::iterator its;
    int ref = -1;

    int u = 0;
    int gvidd = 0;
    for(its=Loc_Elem_Set.begin();its!=Loc_Elem_Set.end();its++)
    {
        int gelid           = *its;

        int lEl             = ElementDistr->getOffsets()[rank]+u+1;
        int Nf              = elements2faces_update[gelid].size();
        int Nv              = elements2verts_update[gelid].size();

        o_tagE2gE[gelid]    = lEl;
        o_gE2tagE[lEl]      = gelid;
        //o_tagE2lE[gelid]    = u;
        //o_lE2tagE[u]        = gelid;

        for(int q=0;q<Nf;q++)
        {
            int gfid = elements2faces_update[gelid][q];

            if(face2verts_update.find(gfid)!=face2verts_update.end())
            {
                int Nv   = face2verts_update[gfid].size();

                if(m_loc_trace_faces.find(gfid)!=m_loc_trace_faces.end())
                {
                    ref  = 13;
                }
                else
                {
                    // ref  = face2reference_update[gfid][0];   
                    ref      = ProvideBoundaryRef(gfid,ranges_ref);         
                }
                
                if(faceset.find(gfid)==faceset.end())
                {
                    faceset.insert(gfid);
                
                    if(o_sharedFace2RankMap.find(gfid)!=o_sharedFace2RankMap.end() && ref!=13)
                    {

                        if(o_sharedFace2RankMap[gfid] == rank)
                        {                                
                            if(o_sharedFace2Nodes.find(gfid)==o_sharedFace2Nodes.end())
                            {
                                int e0   = face2elements_update[gfid][0];
                                int e1   = face2elements_update[gfid][1];
                                
                                int r0;// = part_global[e0]; // rank of first adjacent element.
                                int r1;// = part_global[e1]; // rank of second adjacent element.
                                // int r0 = part_global[e0];
                                // int r1 = part_global[e1];
                                if(part_map.find(e0)!=part_map.end())
                                {
                                    r0 = part_map[e0];

                                }
                                if(adjacent2pid.find(e0)!=adjacent2pid.end())
                                {
                                    r0 = adjacent2pid[e0];
                                }

                                if(part_map.find(e1)!=part_map.end())
                                {
                                    r1 = part_map[e1];
                                }
                                if(adjacent2pid.find(e1)!=adjacent2pid.end())
                                {
                                    r1 = adjacent2pid[e1];
                                }


                                // int r0 = part_global[e0]; // rank of first adjacent element.
                                // int r1 = part_global[e1]; // rank of second adjacent element.
                                
                                if(r0==rank && r1!=rank)
                                {
                                    o_colorRh[r1].push_back(e1);
                                    o_colorFh[r1].push_back(gfid);
                                }
                                else if(r1==rank && r0!=rank)
                                {
                                    o_colorRh[r0].push_back(e0);
                                    o_colorFh[r0].push_back(gfid);
                                }
                                
                                std::vector<int> fn_tag(Nv,0);
                                for(int n=0;n<Nv;n++)
                                {
                                    fn_tag[n] = face2verts_update[gfid][n];

                                    if(o_SharedVertsOwned.find(fn_tag[n])==o_SharedVertsOwned.end() 
                                            && m_loc_trace_verts_glob.find(fn_tag[n])==m_loc_trace_verts_glob.end())
                                    {
                                        if(o_SharedVertsNotOwned.find(fn_tag[n])==o_SharedVertsNotOwned.end())
                                        {
                                            o_SharedVertsNotOwned[fn_tag[n]] = o_sharedVertexMapUpdatedGlobalID[fn_tag[n]];
                                        }
                                    }
                                }

                                o_sharedFace2Nodes[gfid] = fn_tag;
                                
                                o_lhp[gfid] = lEl;

                            }
                        }
                    }
                    else
                    {
                        if(ref == 2)
                        {
                            if(o_interiorFace2Nodes.find(gfid)==o_interiorFace2Nodes.end())
                            {
                                std::vector<int> fn_tag(Nv,0);
                                for(int n=0;n<Nv;n++)
                                {
                                    fn_tag[n] = face2verts_update[gfid][n];
                                    
                                    if(o_sharedVertex2RankMap.find(fn_tag[n])!=o_sharedVertex2RankMap.end())
                                    {
                                        if(o_SharedVertsOwned.find(fn_tag[n])==o_SharedVertsOwned.end() 
                                                && m_loc_trace_verts_glob.find(fn_tag[n])==m_loc_trace_verts_glob.end())
                                        {
                                            if(o_SharedVertsNotOwned.find(fn_tag[n])==o_SharedVertsNotOwned.end())
                                            {
                                                o_SharedVertsNotOwned[fn_tag[n]] = o_sharedVertexMapUpdatedGlobalID[fn_tag[n]];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if(o_NonSharedVertsOwned.find(fn_tag[n])==o_NonSharedVertsOwned.end() 
                                        && m_loc_trace_verts_glob.find(fn_tag[n])==m_loc_trace_verts_glob.end())
                                        {
                                            o_NonSharedVertsOwned[fn_tag[n]]  = gvidd;
                                            gvidd++;
                                        }
                                    }
                                }
                                o_lhp[gfid] = lEl;
                                o_interiorFace2Nodes[gfid] = fn_tag;
                            }
                        }
                        if(ref != 2 && ref!=13)
                        {

                            int fzone = ProvideBoundaryID(gfid,ranges_id);
                            o_zone2bcface[fzone].push_back(gfid);

                            if(o_boundaryFace2Nodes.find(gfid)==o_boundaryFace2Nodes.end())
                            {
                                std::vector<int> fn_tag(Nv,0);
                                for(int n=0;n<Nv;n++)
                                {
                                    fn_tag[n] = face2verts_update[gfid][n];

                                    if(o_sharedVertex2RankMap.find(fn_tag[n])!=o_sharedVertex2RankMap.end())
                                    {
                                        if(o_SharedVertsOwned.find(fn_tag[n])==o_SharedVertsOwned.end() 
                                                && m_loc_trace_verts_glob.find(fn_tag[n])==m_loc_trace_verts_glob.end())
                                        {
                                            if(o_SharedVertsNotOwned.find(fn_tag[n])==o_SharedVertsNotOwned.end())
                                            {
                                                o_SharedVertsNotOwned[fn_tag[n]] = o_sharedVertexMapUpdatedGlobalID[fn_tag[n]];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if(o_NonSharedVertsOwned.find(fn_tag[n])==o_NonSharedVertsOwned.end() 
                                        && m_loc_trace_verts_glob.find(fn_tag[n])==m_loc_trace_verts_glob.end())
                                        {
                                            o_NonSharedVertsOwned[fn_tag[n]]  = gvidd;
                                            gvidd++;
                                        }
                                    }
                                    
                                }

                                o_lhp[gfid] = lEl;
                                o_boundaryFace2Nodes[gfid] = fn_tag;
                                
                            }
                        }
                    }
                }
                else
                {
                    tag2ElementID[gelid] = lEl;

                    o_rhp[gfid]          = lEl;
                }

            }

            
        }
        u++;    
    }


    std::map<int,std::vector<int> >::iterator itv;
    ScheduleObj* rh_schedule = DoScheduling(o_colorRh,comm);
    std::map<int,std::vector<int> > recv_rhElIDs;
    for(int q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            
            for (itv = o_colorRh.begin(); itv != o_colorRh.end(); itv++)
            {
                int n_req           = itv->second.size();
                int dest            = itv->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
                MPI_Send(&itv->second[0], n_req, MPI_INT, dest, 14876+dest, comm);
                i++;
            }
        }
        else if (rh_schedule->SendFromRank2Rank[q].find( rank ) != rh_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 
                    14876+rank, comm, MPI_STATUS_IGNORE);

            recv_rhElIDs[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> > sendEl;
    std::map<int,std::vector<int> >::iterator rcvit;
    for(rcvit=recv_rhElIDs.begin();rcvit!=recv_rhElIDs.end();rcvit++)
    {
        int frank = rcvit->first;
        int nE    = rcvit->second.size();
        
        for(int j=0;j<nE;j++)
        {
            int gEl = o_tagE2gE[rcvit->second[j]];
            sendEl[frank].push_back(gEl);
        }
    }
    
    ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);

    std::map<int,std::vector<int> > adj_ids;
    for(int q=0;q<size;q++)
    {
        if(rank==q)
        {
            int i=0;
            for (itv = sendEl.begin(); itv != sendEl.end(); itv++)
            {
                int n_req           = itv->second.size();
                int dest            = itv->first;

                MPI_Send(&n_req, 1,
                        MPI_INT, dest,
                        6798+78000*dest, comm);
                
                MPI_Send(&itv->second[0],
                        n_req, MPI_INT,
                        dest, 14876000+dest, comm);

                i++;
            }
        }
        else if (ishBack_schedule->SendFromRank2Rank[q].find( rank ) != ishBack_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            
            MPI_Recv(&n_reqstd_ids,
                    1, MPI_INT, q,
                    6798+78000*rank,
                    comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0],
                    n_reqstd_ids,
                    MPI_INT, q,
                    14876000+rank,
                    comm, MPI_STATUS_IGNORE);
            
            adj_ids[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,int> shFid2el_rh;
    int adde = 0;
    int adde2 = 0;
    for(rcvit=adj_ids.begin();rcvit!=adj_ids.end();rcvit++)
    {
        
        int rrank = rcvit->first;
        int nelem = rcvit->second.size();
        for(int q=0;q<nelem;q++)
        {
            int fid = o_colorFh[rrank][q];
            
            if(shFid2el_rh.find(fid)==shFid2el_rh.end())
            {
                shFid2el_rh[fid] = rcvit->second[q];
                adde2++;
                if(o_rhp.find(fid)==o_rhp.end())
                {
                    adde++;
                    o_rhp[fid] = rcvit->second[q];
                }
            }
        }
    }


    //std::cout << "nLocShVrts " << nLocShVrts << std::endl;
    m_loc_trace_verts_glob.clear();
}



std::map<int,std::set<int> > RepartitionObject::GetNode2ElementMap()
{
    std::map<int,std::set<int> > node2elem_map;

    //std::map<int,std::vector<int> > loc_prisms = prisms->getElement2VertexMap();

    for(int i=0;i<Loc_Elem.size();i++)
    {
        int elid = Loc_Elem[i];
        int Nv   = elements2verts_update[elid].size();
        int nadj = elements2elements_update[elid].size();

        for(int j=0;j<Nv;j++)
        {
            int gvid = elements2verts_update[elid][j];
            node2elem_map[gvid].insert(elid);
        }

        /*
        for(int j=0;j<nadj;j++)
        {
            int adjID = elements2verts_update[elid][j];
        
            if(adjID < Ne_glob && adjID!=elid)
            {
                if(elements2verts_update.find(adjID)!=elements2verts_update.end())
                {
                    int Nvadj = elements2verts_update[adjID].size();

                    for(int j=0;j<Nvadj;j++)
                    {
                        int gvid = elements2verts_update[adjID][j];
                        node2elem_map[gvid].insert(adjID);
                    }

                    int nadjadj = elements2elements_update[adjID].size();

                    for(int j=0;j<nadjadj;j++)
                    {
                        int adjadjID    = elements2verts_update[adjID][j];

                        if(adjadjID < Ne_glob && adjadjID!=elid)
                        {
                            if(elements2verts_update.find(adjadjID)!=elements2verts_update.end())
                            {
                                int Nvadjadj = elements2verts_update[adjadjID].size();

                                for(int j=0;j<Nvadjadj;j++)
                                {
                                    int gvid = elements2verts_update[adjadjID][j];
                                    node2elem_map[gvid].insert(adjadjID);
                                }
                            }
                        }
                        if(loc_prisms.find(adjadjID)!=loc_prisms.end())
                        {
                            //std::cout << "adjadjID in prisms " << std::endl;  
                            int Nvadjadj = loc_prisms[adjadjID].size();

                            for(int j=0;j<Nvadjadj;j++)
                            {
                                int gvid = loc_prisms[adjadjID][j];
                                node2elem_map[gvid].insert(adjadjID);
                            }
                        }       
                    }
                }
            }
            if(loc_prisms.find(adjID)!=loc_prisms.end())
            {
                int Nvadj = loc_prisms[adjID].size();

                for(int j=0;j<Nvadj;j++)
                {
                    int gvid = loc_prisms[adjID][j];
                    node2elem_map[gvid].insert(adjID);
                }
            }
        }*/
        
    }

    return node2elem_map;

}



std::map<int,std::map<int,double> > RepartitionObject::GetNode2ElementMapV2()
{
    std::map<int,std::map<int,double> > node2elem_map;

    for(int i=0;i<LocAndAdj_Elem.size();i++)
    {
        int elid = LocAndAdj_Elem[i];
        int Nv   = elements2verts_update[elid].size();
        int nadj = elements2elements_update[elid].size();

        std::vector<double> Pijk(Nv*3);
        for(int k=0;k<Nv;k++)
        {
            int gvid     = elements2verts_update[elid][k];
            Pijk[k*3+0]  = LocalVertsMap[gvid][0];
            Pijk[k*3+1]  = LocalVertsMap[gvid][1];
            Pijk[k*3+2]  = LocalVertsMap[gvid][2];
        }

        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,Nv);

        for(int j=0;j<Nv;j++)
        {
            int gvid = elements2verts_update[elid][j];

            if(LocalVertsMap.find(gvid)!=LocalVertsMap.end()
                    && node2elem_map[gvid].find(elid)==node2elem_map[gvid].end())
            {
                double dx = Vijk[0]-LocalVertsMap[gvid][0];
                double dy = Vijk[1]-LocalVertsMap[gvid][1];
                double dz = Vijk[2]-LocalVertsMap[gvid][2];
                double d  = sqrt(dx*dx+dy*dy+dz*dz);
                node2elem_map[gvid].insert(std::pair<int,double>(elid,d));
            }
        } 
    }

    return node2elem_map;
}







std::map<int,std::map<int, double> > RepartitionObject::GetElement2NodeMap()
{
    std::map<int,std::map<int, double> > elem2nodemap;

    for(int i=0;i<Loc_Elem.size();i++)
    {
        int elid = Loc_Elem[i];
        int Ne   = elements2elements_update[elid].size();
        int Nv   = elements2verts_update[elid].size();
        int nadj = elements2elements_update[elid].size();

        std::vector<double> Pijk(Nv*3);

        for(int k=0;k<Nv;k++)
        {
            int gvid     = elements2verts_update[elid][k];
            Pijk[k*3+0]  = LocalVertsMap[gvid][0];
            Pijk[k*3+1]  = LocalVertsMap[gvid][1];
            Pijk[k*3+2]  = LocalVertsMap[gvid][2];
        }
        std::vector<double> Vijk     = ComputeCentroidCoord(Pijk,Nv);

        for(int k=0;k<Nv;k++)
        {
            int gvid = elements2verts_update[elid][k];
            if(LocalVertsMap.find(gvid)!=LocalVertsMap.end()
                    && elem2nodemap[elid].find(gvid)==elem2nodemap[elid].end())
            {
                double dx = Vijk[0]-LocalVertsMap[gvid][0];
                double dy = Vijk[1]-LocalVertsMap[gvid][1];
                double dz = Vijk[2]-LocalVertsMap[gvid][2];
                double d  = sqrt(dx*dx+dy*dy+dz*dz);
                elem2nodemap[elid].insert(std::pair<int,double>(gvid,d));
            }
        }

        for(int s=0;s<Ne;s++)
        {
            int adjid = elements2elements_update[elid][s];
            int Nvadj = elements2verts_update[adjid].size();

            std::map<int, double> vrt2dist;

            for(int j=0;j<Nvadj;j++)
            {
                int gvid = elements2verts_update[adjid][j];
                if(LocalVertsMap.find(gvid)!=LocalVertsMap.end()
                    && elem2nodemap[elid].find(gvid)==elem2nodemap[elid].end())
                {
                    double dx = Vijk[0]-LocalVertsMap[gvid][0];
                    double dy = Vijk[1]-LocalVertsMap[gvid][1];
                    double dz = Vijk[2]-LocalVertsMap[gvid][2];
                    double d  = sqrt(dx*dx+dy*dy+dz*dz);
                    elem2nodemap[elid].insert(std::pair<int,double>(gvid,d));
                }
            }
        }
    }

    return elem2nodemap;

}







std::map<int,std::vector<double> > RepartitionObject::ReduceStateVecToVertices(std::map<int,std::map<int,double> > node2elem_map,
                                                                                std::map<int,std::vector<double> > Umap,
                                                                                int nvar)
{
    std::map<int,std::map<int,double> >::iterator n2eit;
    std::map<int,std::vector<double> > nodevals;
    
    for(n2eit=node2elem_map.begin();n2eit!=node2elem_map.end();n2eit++)
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
                uval[q] = uval[q] + Umap[elid][q];//*dist;
            }
            
            distsum  = distsum+dist;
        }
        for(int q=0;q<nvar;q++)
        {
            uval[q] = uval[q]/elidsdist.size();// / distsum;
            //std::cout << "uval[q] " << uval[q] << " " << q << std::endl;
        }

        nodevals[gvid] = uval;
        
    }

    
    return nodevals;
}







std::map<int,std::vector<int> > RepartitionObject::GetLocalTraceFace2LeftRight()
{
    return m_loc_trace_face2leftright;
}
std::set<int> RepartitionObject::GetLocalTraceFacesSet()
{
    return m_loc_trace_faces;
}
std::set<int> RepartitionObject::GetLocalTraceVertSet()
{
    return m_loc_trace_verts;
}


void RepartitionObject::SetStateVec(std::map<int,std::vector<double> > U, int nvar)
{
    elements2data_update.clear();
    elements2data_update = U;
}

std::map<int,std::vector<int> > RepartitionObject::getFace2ElementMap()
{
    return face2elements_update;
}

std::map<int,int> RepartitionObject::GetElement2TypeOnRankMap()
{
    return elem2type_on_rank;
}

std::map<int,std::vector<int> > RepartitionObject::getZone2boundaryFaceID()
{
    return o_zone2bcface;
}
std::map<int,int> RepartitionObject::GetLeftHandFaceElementMap()
{
    return o_lhp;
}
std::map<int,int> RepartitionObject::GetRightHandFaceElementMap()
{
    return o_rhp;
}
std::map<int,int> RepartitionObject::GetLocalSharedFace2GlobalSharedFace()
{
    return o_locShF2globShF;
}
std::map<int,int> RepartitionObject::getSharedVertsNotOwned()
{
    return o_SharedVertsNotOwned;
}
std::map<int,int> RepartitionObject::getNonSharedVertsOwned()
{
    return o_NonSharedVertsOwned;
}
std::map<int,int> RepartitionObject::getSharedVertsOwned()
{
    return o_SharedVertsOwned;
}
std::map<int,int> RepartitionObject::getUpdatedTag2GlobalVMap()
{
    return o_tag2globV;
}
// std::map<int,int> RepartitionObject::getUpdatedGlobal2TagVMap()
// {
//     return o_glob2tagV;
// }
// std::map<int,int> RepartitionObject::getUpdatedGlobal2LocalFMap()
// {
//     return o_glob2locF;
// }
// std::map<int,int> RepartitionObject::getUpdatedLocal2GlobalFMap()
// {
//     return o_loc2globF;
// }
std::map<int,int> RepartitionObject::getUpdatedGlobal2LocalVMap()
{
    return o_glob2locV;
}
// std::map<int,int> RepartitionObject::getUpdatedLocal2GlobalVMap()
// {
//     return o_loc2globV;
// }
std::map<int,std::vector<int> > RepartitionObject::getFaceTag2VertTagMap()
{
    return o_face2verts_global;
}
std::map<int,int> RepartitionObject::getGlobal2TagFMap()
{
    return o_glob2tagF;
}
std::vector<int> RepartitionObject::getFace4ParMMG()
{
    return o_faces4parmmg;
}
int** RepartitionObject::getParMMGCommFace2GlobalVertMap()
{
    return o_ifc_tria_glo;
}

int** RepartitionObject::getParMMGCommFace2LocalVertMap()
{
    return o_ifc_tria_loc;
}

int* RepartitionObject::getParMMGCommColorFace()
{
    return o_color_face;
}

int* RepartitionObject::getParMMGCommNFacesPerColor()
{
    return o_ntifc;
}
int RepartitionObject::getParMMGNComm()
{
    return o_ncomm;
}
std::map<int,int> RepartitionObject::getLocalVert2VertTag()
{
    return o_lvertex2gvertex_part;
}
// std::map<int,int> RepartitionObject::getVertTag2LocalVert()
// {
//     return o_gvertex2lvertex_part;
// }

std::map<int,int> RepartitionObject::getGlobalElement2ElementTag()
{
    return o_gE2tagE;
}
std::map<int,int> RepartitionObject::getElementTag2GlobalElement()
{
    return o_tagE2gE;
}

// std::map<int,int> RepartitionObject::getElementTag2LocalElement()
// {
//     return o_tagE2lE;
// }
// std::map<int,int> RepartitionObject::getLocalElement2ElementTag()
// {
//     return o_lE2tagE;
// }

// std::map<int,int> RepartitionObject::getGlob2TagElementID()
// {
//     return o_tagE2gE;
// }
std::map<int, std::vector<int> > RepartitionObject::getSharedFaceMap()
{
    return o_sharedFace2Nodes;
}
std::map<int, std::vector<int> > RepartitionObject::getInteriorFaceMap()
{
    return o_interiorFace2Nodes;
}
std::map<int, std::vector<int> > RepartitionObject::getBoundaryFaceMap()
{
    return o_boundaryFace2Nodes;
}
std::map<int, std::vector<int> > RepartitionObject::getFace2VertexMap()
{
    return face2verts_update;
}

// std::map<int, std::vector<int> > RepartitionObject::getFace2NVertexMap()
// {
//     return face2Nverts_update;
// }

std::map<int,std::vector<int> > RepartitionObject::getElement2VertexMap()
{
    return elements2verts_update;
}

std::map<int,std::vector<int> > RepartitionObject::getElement2ElementMap()
{
    return elements2elements_update;
}

std::map<int,std::vector<int> > RepartitionObject::getElement2FacesMap()
{
    return elements2faces_update;
}


std::map<int,std::vector<double> > RepartitionObject::getElement2DataMap()
{
    return elements2data_update;
}


std::map<int,std::vector<int> > RepartitionObject::getGlobalElement2LocalVertMap()
{
    return globElem2locVerts;
}


std::map<int, std::vector<double> >  RepartitionObject::getLocalVertsMap()
{
    return LocalVertsMap;
}


std::map<int,int> RepartitionObject::getGlobalElement2Rank()
{
    return part_global;
}

std::vector<int> RepartitionObject::getLocElem()
{
    return Loc_Elem;
}

// std::map<int,int> RepartitionObject::getTag2ElementTrace()
// {
//     return tag2element_trace;
// }

// destructor
RepartitionObject::~RepartitionObject()
{
    // LocalVertsMap.clear();
    // Loc_Elem.clear();
    // part_global.clear();
    // globElem2locVerts.clear();
    // elements2data_update.clear();
    // elements2faces_update.clear();
    // elements2elements_update.clear();
    // elements2verts_update.clear();
    // face2reference_update.clear();
    // face2verts_update.clear();

    // o_boundaryFace2Nodes.clear();
    // o_sharedFace2Nodes.clear();
    // o_interiorFace2Nodes.clear();
    // o_tagE2gE.clear();
    // o_gE2tagE.clear();
    // o_lvertex2gvertex_part.clear();
    // o_face2verts_global.clear();
    // o_glob2locV.clear();
    // o_tag2globV.clear();
    // o_SharedVertsOwned.clear();
    // o_NonSharedVertsOwned.clear();
    // o_SharedVertsNotOwned.clear();
    // o_locShF2globShF.clear();
    // o_rhp.clear();
    // o_lhp.clear();
    // o_zone2bcface.clear();
    // elem2type_on_rank.clear();
    // face2elements_update.clear();
    // m_loc_trace_face2leftright.clear();
    // m_loc_trace_faces.clear();
    // m_loc_trace_verts.clear();
}

