#include "adapt_runparmmg.h"



#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))




PMMG_pParMesh InitializeParMMGmesh(MPI_Comm comm, 
                                   std::set<int> loc_trace_faces,
                                   std::set<int> loc_trace_verts,
                                   std::map<int,std::vector<int> > loc_trace_face2leftright,
                                   PrepareAdaption* prepare,
                                   std::map<int,std::vector<int> > ranges_id,
                                   int bndIDmax,
                                   std::map<int, std::vector<double> > metric_vmap)
{
    
    int ier;

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    PMMG_pParMesh   parmesh;
    int k;
    
    // This test assumes that the following two build requests on the tetra_repart object have been called before entering this routine
    // tetra_repart->buildInteriorSharedAndBoundaryFacesMaps(comm,pttrace, ranges_id);
    // tetra_repart->buildTag2GlobalElementMapsAndFace2LeftRightElementMap(comm,pttrace, ranges_id);

    // tetra_repart->buildCommunicationMaps(comm);

    // std::map<int,std::vector<int> > face2vertsMap               = tetra_repart->getFace2VertMap();
    // std::set<int> loc_trace_faces                               = tetra_repart->GetLocalTraceFacesSet();
    // std::map<int,std::vector<int> > loc_trace_face2leftright    = tetra_repart->GetLocalTraceFace2LeftRight();
    // std::set<int> loc_trace_verts                               = tetra_repart->GetLocalTraceVertsSet();

    std::map<int,std::vector<int> > gE2gV_t                     = prepare->getElem2VertMap();
    std::map<int,std::vector<int> > face2vertsMap               = prepare->getFace2VertMap();
    std::map<int,int> locv2tagvID                               = prepare->getLocalVert2GlobalVert();
    std::set<int> Owned_Elem_t                                  = prepare->getElemSet();
    std::map<int,int> globalv2localvID                          = prepare->getNewGlobalVert2LocalVertMap();
    std::map<int,int> tag2globalV                               = prepare->getTagVert2GlobalVertMap();
    
    std::vector<int> face4parmmg                                = prepare->getFace4ParMMG(); // checked
    std::map<int, std::vector<double> > LocalVertsMap_t         = prepare->getLocalVertsMap();
    std::map<int,int> global2tagF                               = prepare->getGlobal2TagFMap();
    
    int** ifc_tria_glob                                         = prepare->getParMMGCommFace2GlobalVertMap();
    int** ifc_tria_loc                                          = prepare->getParMMGCommFace2LocalVertMap();
    int* color_face                                             = prepare->getParMMGCommColorFace();
    int *ntifc                                                  = prepare->getParMMGCommNFacesPerColor();
    int ncomm                                                   = prepare->getParMMGNComm();


    std::set<int> glob_trace_verts = AllGatherSet(loc_trace_verts,comm);
    
    int nVertices   = locv2tagvID.size();
    int nTetrahedra = Owned_Elem_t.size();
    int nEdges      = 0;
    int nTriangles  = face4parmmg.size();

    std::cout << nTriangles << " = nTriangles   " << nVertices << " = nVertices " << nTetrahedra << " = nTetrahedra" << world_rank << std::endl;
    
    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);

    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,0,nTriangles,0,nEdges) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    
    
    // std::cout << "nVertices " << nVertices << " " << world_rank << std::endl;
    //PMMG_Set_metSize(PMMG_pParMesh parmesh,int typEntity,int np,int typSol)
    if ( PMMG_Set_metSize(parmesh,MMG5_Vertex,nVertices,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);

    for (k=0; k<nVertices; ++k )
    {
        int tagvid  = locv2tagvID[k];
        int glovid  = tag2globalV[tagvid];
        int locvid  = globalv2localvID[glovid];
        //std::cout << "locvid = " << locvid << " --- " << glovid << std::endl;
        double vx = LocalVertsMap_t[tagvid][0];
        double vy = LocalVertsMap_t[tagvid][1];
        double vz = LocalVertsMap_t[tagvid][2];

        //std::cout <<"coords = "<< vx << " " << vy << " " << vz << std::endl;
        
        int vref = 86;

        if ( PMMG_Set_vertex(parmesh,vx,vy,vz, vref, locvid+1) != 1 )
        {
        MPI_Finalize();
        exit(EXIT_FAILURE);
        }

        std::vector<double> tensor(6,0);
        tensor[0] = metric_vmap[tagvid][0];
        tensor[1] = metric_vmap[tagvid][1];
        tensor[2] = metric_vmap[tagvid][2];
        tensor[3] = metric_vmap[tagvid][3];
        tensor[4] = metric_vmap[tagvid][4];
        tensor[5] = metric_vmap[tagvid][5];

        // std::cout << metric_vmap[tagvid][0] << " " << metric_vmap[tagvid][1] << " " << metric_vmap[tagvid][2] << " " << metric_vmap[tagvid][3] << " " << metric_vmap[tagvid][4] << " " << metric_vmap[tagvid][5] << std::endl;  

        if(PMMG_Set_tensorMet(parmesh,tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],locvid+1)!=1)
        {
         MPI_Finalize();
         exit(EXIT_FAILURE);
        }
    }
    

    int refer  = 0;
    //double* c0 = new double[3];
    int iref;
    int vertref = 86;
    int vertref2 = 86;
    int locs = 0;
    
    std::map<int,int> shell_g2l;
    std::vector<std::vector<double> > unshellVin;
    std::vector<std::vector<int> > unshellTin;
    int outflowbc  = 0;
    int inflowbc   = 0;
    int tracebc    = 0;
    int symmetrybc = 0;
    int wallbc = 0;
    int nothere = 0;

    
    //std::map<int,std::vector<int> > f2refmap = tetra_repart->getFace2RefMap();
    //std::map<int,std::vector<int> > leftright_trace     = pttrace->GetLeftRightElements();
    int fff = 0;


    std::vector<int> bndid_vec;
    std::map<int,std::vector<int> >::iterator itr;
    for(itr=ranges_id.begin();itr!=ranges_id.end();itr++)
    {
        bndid_vec.push_back(itr->first);
    }

    //bndIDmax = *std::max_element(bndid_vec.begin(), bndid_vec.end());
    //std::cout << "bndIDmax " << bndIDmax << " -- " << bndid_vec.size() << std::endl;
    int found1 = 0;
    int found2 = 0;
    std::cout << world_rank << " nTriangles. " << nTriangles << std::endl;
    for ( k=0; k<nTriangles; ++k )
    {
        
        int faceID      = face4parmmg[k];
        int tagFaceID   = global2tagF[faceID];
        // int ref      = f2refmap[tagFaceID][0];
        int ref = -1;
        
        if(ref==3)
        {
            wallbc++;
        }
        if(ref==7)
        {
            symmetrybc++;
        }
        if(ref==10)
        {
            inflowbc++;
        }
        if(ref==36)
        {
            outflowbc++;
        }

        if(loc_trace_faces.find(tagFaceID)!=loc_trace_faces.end())
        {
            
            // ref = 13;//-(leftright_trace[tagFaceID][1]+1);
            //ref = 1.0e09+(leftright_trace[tagFaceID][1]);
            ref = (bndIDmax+1)+(loc_trace_face2leftright[tagFaceID][1]);
            // if(leftright_trace.find(tagFaceID)!=leftright_trace.end())
            // {
            //     ref = -(leftright_trace[tagFaceID][1]+1);
            //     fff++;
            // }
            // else
            // {
            //     ref = -13;
            // }
            for(int s=0;s<3;s++)
            {
                //int ref   = 13;

                int tagvid      = face2vertsMap[tagFaceID][s];
                double vx       = LocalVertsMap_t[tagvid][0];
                double vy       = LocalVertsMap_t[tagvid][1];
                double vz       = LocalVertsMap_t[tagvid][2];
                int glovid      = tag2globalV[tagvid];
                int locvid      = globalv2localvID[glovid];

                if(glob_trace_verts.find(tagvid)!=glob_trace_verts.end())
                {
                    //vertref = loc_trace2ref[tagvid];
                    vertref = tagvid+100;
                    //std::cout << "vertref " << vertref << " " << tagvid << std::endl;
                    found1++;
                }
                else
                {
                    std::cout << "NOT HERE " << tagvid << std::endl;
                }
                // if(unique_trace_verts2refmap.find(tagvid)!=unique_trace_verts2refmap.end())
                // {
                //     vertref2 = unique_trace_verts2refmap[tagvid];
                //     found2++;
                // }
                // else
                // {
                //     std::cout << "Warning:: reference value on trace vertID " << tagvid << " is wrong!" << std::endl;
                // }

                //std::cout << "vertref " << vertref << " " << vertref2 << std::endl;

                if ( PMMG_Set_vertex(parmesh, vx, vy, vz, vertref, locvid+1) != 1 )
                {
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
            }
            
            // int tagvid0   = trace_verts[tagFaceID][0];
            // int tagvid1   = trace_verts[tagFaceID][1];
            // int tagvid2   = trace_verts[tagFaceID][2];

            int tagvid0   = face2vertsMap[tagFaceID][0];
            int tagvid1   = face2vertsMap[tagFaceID][1];
            int tagvid2   = face2vertsMap[tagFaceID][2];

            int glovid0   = tag2globalV[tagvid0];
            int glovid1   = tag2globalV[tagvid1];
            int glovid2   = tag2globalV[tagvid2];

            int locvid0   = globalv2localvID[glovid0];
            int locvid1   = globalv2localvID[glovid1];
            int locvid2   = globalv2localvID[glovid2];
        
            if ( PMMG_Set_triangle(parmesh,locvid0+1,locvid1+1,locvid2+1,ref,k+1) != 1 )
            {
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }

            PMMG_Set_requiredTriangle( parmesh, k+1 );

            tracebc++;
        }
        else
        {
            ref 		  = ProvideBoundaryID(tagFaceID,ranges_id);
            

            int tagvid0   = face2vertsMap[tagFaceID][0];
            int tagvid1   = face2vertsMap[tagFaceID][1];
            int tagvid2   = face2vertsMap[tagFaceID][2];

            int glovid0   = tag2globalV[tagvid0];
            int glovid1   = tag2globalV[tagvid1];
            int glovid2   = tag2globalV[tagvid2];

            int locvid0   = globalv2localvID[glovid0];
            int locvid1   = globalv2localvID[glovid1];
            int locvid2   = globalv2localvID[glovid2];
            
            if ( PMMG_Set_triangle(parmesh,locvid0+1,locvid1+1,locvid2+1,ref,k+1) != 1 )
            {
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }
        // std::cout << "face ref " << ref << std::endl;        
    }

    std::cout << "tracebc " << tracebc << " " << world_rank << std::endl;
    glob_trace_verts.clear();
    //std::cout << "found1 found2 " << found1 << " " << found2 << std::endl;
    //std::cout << "fff " << fff << std::endl; 
    
    int outflowbc_sum   = 0;
    int inflowbc_sum    = 0;
    int tracebc_sum     = 0;
    int wallbc_sum      = 0;
    int symmetrybc_sum  = 0;

    MPI_Allreduce(&outflowbc, &outflowbc_sum, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&inflowbc, &inflowbc_sum, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&tracebc, &tracebc_sum, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&wallbc, &wallbc_sum, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&symmetrybc, &symmetrybc_sum, 1, MPI_INT, MPI_SUM, comm);

    // if(world_rank == 0)
    // {
    //     std::cout << "=============Boundary Face Stats for Tetrahedra ============" << std::endl;
    //     std::cout << "                  Outflow: " << outflowbc_sum << std::endl;
    //     std::cout << "                  Inflow: " << inflowbc_sum << std::endl;
    //     std::cout << "                  Trace: " << tracebc_sum << std::endl;
    //     std::cout << "                  Wall: " << wallbc_sum << std::endl;
    //     std::cout << "                  Symmetry: " << symmetrybc_sum << std::endl;
    //     std::cout << "============================================================" << std::endl;
    // }

    int API_mode = 0;
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };


    // tetra_repart->buildParMMGCommunicationMaps(comm);


        
    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);

   
    for(int icomm=0; icomm<ncomm; icomm++ )
    {
      //std::cout << "Initialization " << icomm << " " << ncomm << " " << ntifc[icomm] << " " << world_rank << std::endl;
      // Set nb. of entities on interface and rank of the outward proc
      ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                             color_face[icomm],
                                             ntifc[icomm]);

        //Set local and global index for each entity on the interface
      ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                               ifc_tria_loc[icomm],
                                               ifc_tria_glob[icomm], 1);
    }

    
   
    // //==========================Clear from here==============================
    
    std::map<int,std::vector<int> >::iterator ittet;
    k = 0;
    int vloc0,vloc1,vloc2,vloc3;
    int vglo0,vglo1,vglo2,vglo3;
    int vtag0,vtag1,vtag2,vtag3;
    // for ( int t = 0;t < nTetrahedra; t++  )
    std::set<int>::iterator its;
    int t = 0;
    
    for(its=Owned_Elem_t.begin();its!=Owned_Elem_t.end();its++)
    {
        //From local element ID get the tag ID (initial global ID).
        int gEid = *its; //Owned_Elem_t[t];
        //Using the tag Element ID, get the tag Vertex ID.
        vtag0 = gE2gV_t[gEid][0];
        vtag1 = gE2gV_t[gEid][1];
        vtag2 = gE2gV_t[gEid][2];
        vtag3 = gE2gV_t[gEid][3];

        //From the tag vertex ID, get the local vertex ID.

        vglo0 = tag2globalV[vtag0];
        vglo1 = tag2globalV[vtag1];
        vglo2 = tag2globalV[vtag2];
        vglo3 = tag2globalV[vtag3];

        vloc0 = globalv2localvID[vglo0];
        vloc1 = globalv2localvID[vglo1];
        vloc2 = globalv2localvID[vglo2];
        vloc3 = globalv2localvID[vglo3];

        std::vector<double> P(4*3);
        P[0*3+0]=LocalVertsMap_t[vtag0][0];   P[0*3+1]=LocalVertsMap_t[vtag0][1];    P[0*3+2]=LocalVertsMap_t[vtag0][2];
        P[1*3+0]=LocalVertsMap_t[vtag1][0];   P[1*3+1]=LocalVertsMap_t[vtag1][1];    P[1*3+2]=LocalVertsMap_t[vtag1][2];
        P[2*3+0]=LocalVertsMap_t[vtag2][0];   P[2*3+1]=LocalVertsMap_t[vtag2][1];    P[2*3+2]=LocalVertsMap_t[vtag2][2];
        P[3*3+0]=LocalVertsMap_t[vtag3][0];   P[3*3+1]=LocalVertsMap_t[vtag3][1];    P[3*3+2]=LocalVertsMap_t[vtag3][2];

        double Vtet = GetQualityTetrahedra(P);

        //std::cout << "Vtet " << Vtet << "  " << vloc1+1 << " " << vloc2+1 << " " << vloc3+1 << std::endl;
        
        if(Vtet<0.0)
        {
            std::cout << " negative volume in Element " << t << " on rank " << world_rank  <<std::endl;
        }
        if ( PMMG_Set_tetrahedron(parmesh,vloc0+1,vloc1+1,vloc2+1,vloc3+1,1.0,t+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        t++;
    }
    
    

    return parmesh;
}



void RunParMMGandWriteTetraUS3Dformat(MPI_Comm comm, 
                PMMG_pParMesh &parmesh,
                std::map<int,int> part_global_bl,
                int bndIDmax,
                Inputs* inputs, 
                int nElemsGlob_P, int nVertsGlob_P, 
                std::map<int,int>& tracerefV2globalV, 
                std::vector<int> &ifn_T,
                std::map<int,std::vector<int> > &bcmap,
                std::map<int,std::vector<int> > &BoundaryFaces_T,
                std::map<int,int> &glob2locVid,
                std::map<int,int> &LocationSharedVert_update,
                std::vector<double> &vertices_output,
                std::vector<std::vector<int> > &new_tetrahedra,
                std::map<int,int> &oldglob2newglob,
                std::map<int,int> &lh_T_bc,
                std::map<int,int> &new_globE2locE,
                int &nLocIntVrts,
                int &nLocShVrts,
                std::map<int,int> tagE2gE_P)
{
    int i,k;
    std::map<int,std::vector<int> >::iterator itmiv;
    FILE            *inm;

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int ier;

    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, inputs->niter ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hausd, inputs->hausd) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hgrad, inputs->hgrad) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if( !PMMG_Set_dparameter( parmesh,  PMMG_DPARAM_hgradreq , -1.0 ) ){
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
 
    int nVerticesIN   = 0;
    int nTetrahedraIN = 0;
    int nTrianglesIN  = 0;
    int nEdgesIN      = 0;
    if ( PMMG_Get_meshSize(parmesh,&nVerticesIN,&nTetrahedraIN,NULL,&nTrianglesIN,NULL,
                           &nEdgesIN) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }

    // std::cout << "nVertices input " << nVerticesIN << std::endl;
    // std::cout << "nTetrahedra input " << nTetrahedraIN << std::endl;
    // std::cout << "nTriangles input " << nTrianglesIN << std::endl;
    // std::cout << "nEdges input " << nEdgesIN << std::endl;
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_globalNum, 1 ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };

    int ierlib = PMMG_parmmglib_distributed( parmesh );
    
    if(ierlib==0 && world_rank == 0)
    {
        std::cout << "SUCCESFULLY adapted the mesh in parallel!" << std::endl;
    }
    else if(ierlib==1 && world_rank == 0)
    {
        std::cout << "FAILED to adapt the mesh in parallel! "<< std::endl;
    }
    FaceSetPointer::iterator ftit;
    // std::map<int,std::vector<int> > leftright_trace     = pttrace->GetLeftRightElements();
    // FaceSetPointer FaceTraceRefs                        = pttrace->GetRefTraceFaceSet();

    // 
    // std::vector<int> ownedTracePrisms;
    // int lr=0;
    // int found = 0;
    // std::set<int> unique_prisms;
    // std::map<int,int> glob2tag_prisms_on_trace;
    // std::map<int,int> tag2glob_prisms_on_trace;

    // for(ftit=FaceTraceRefs.begin();ftit!=FaceTraceRefs.end();ftit++)
    // {
    //     int prismID = (*ftit)->GetFaceRightElement();

    //     unique_prisms.insert(prismID);

    //     if(tagE2gE_P.find(prismID)!=tagE2gE_P.end())
    //     {
    //         // (*ftit)->SetFaceRightElement(tagE2gE_P[prismID]);
    //         glob2tag_prisms_on_trace[tagE2gE_P[prismID]] = prismID;
    //         tag2glob_prisms_on_trace[prismID]            = tagE2gE_P[prismID];
    //         // ownedTracePrisms.push_back(prismID);
    //         // // if(tagE2gE_P[prismID]>279070)
    //         // {
    //         //     lr++;
    //         // }
    //         // found++;
    //     }
    // }

    // //std::cout << "ownedTracePrisms " << lr << " " << unique_prisms.size() << " " << ownedTracePrisms.size() << " " << world_rank << " " << FaceTraceRefs.size() << std::endl;
    // std::map<int,int> glob2tag_prisms_on_trace_glob = AllGatherMap_T(glob2tag_prisms_on_trace,comm);
    // std::map<int,int> tag2glob_prisms_on_trace_glob = AllGatherMap_T(tag2glob_prisms_on_trace,comm);

    //std::cout << "LR = " << lr << " unique_prisms " << unique_prisms.size() << std::endl; 

    // std::map<int,int>::iterator itmii;
    // k = 0;
    // if(world_rank == 1)
    // {
    //     for(itmii=prisms_on_trace_glob.begin();itmii!=prisms_on_trace_glob.end();itmii++)
    //     {
    //         std::cout << k << " :=> " << itmii->first << " " << itmii->second << std::endl;
    //         k++;
    //     }
    // }
    //  if(world_rank == 1)
    // {
    //     for(itmii=prisms_on_trace_glob.begin();itmii!=prisms_on_trace_glob.end();itmii++)
    //     {
    //         std::cout << world_rank <<  " " << itmii->first << " " << itmii->second << std::endl;
    //     }
    // }
    
    //std::cout << "prisms_on_trace_glob " << prisms_on_trace_glob.size() << std::endl;

    int nVerticesOUT   = 0;
    int nTetrahedraOUT = 0;
    int nTrianglesOUT  = 0;
    int nEdgesOUT      = 0;
    
    if ( PMMG_Get_meshSize(parmesh,&nVerticesOUT,&nTetrahedraOUT,NULL,&nTrianglesOUT,NULL,
                        &nEdgesOUT) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }
    
    int             nodeGloNumber,nodeOwner;
    std::map<int,int> loc2globVid;


    //std::cout << "NVERt ice " << nVerticesOUT << std::endl;

    for( k = 1; k <= nVerticesOUT; k++ )
    {
        if( !PMMG_Get_vertexGloNum( parmesh, &nodeGloNumber, &nodeOwner ) )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    
        loc2globVid[k]=nodeGloNumber;
        glob2locVid[nodeGloNumber]=k;
    }
    
    int *required = (int*)calloc(MAX4(nVerticesOUT,
                                        nTetrahedraOUT,
                                        nTrianglesOUT,
                                        nEdgesOUT),sizeof(int));
    
    int *refOUT = (int*)calloc(MAX4(nVerticesOUT,
                                    nTetrahedraOUT,
                                    nTrianglesOUT,
                                    nEdgesOUT),sizeof(int));
    
    int *corner = (int*)calloc(nVerticesOUT,sizeof(int));

    //int *refOUT = new int[nVerticesOUT];

    std::vector<std::vector<int> > outT;
    //===================================================================
    // Collect Vertices
    //===================================================================
    // std::vector<std::vector<double> > new_vertices(nVerticesOUT);
    // new_vertices.size(nVerticesOUT);
    // std::fill(new_vertices.begin(), new_vertices.end(), 0);
    std::vector<std::vector<double> > new_vertices(nVerticesOUT);



    
    for ( k=0; k<nVerticesOUT; k++ )
    {
        std::vector<double> crds(3,0.0);

        if ( PMMG_Get_vertex(parmesh,&(crds[0]),&(crds[1]),&(crds[2]),
                            &(refOUT[k]),&(corner[k]),&(required[k])) != 1 ) {
            fprintf(inm,"Unable to get mesh vertex %d \n",k);
            ier = PMMG_STRONGFAILURE;
        }

        new_vertices[k]=crds;
    }

    //===================================================================
    // Collect Faces    ->      (triangles)
    //===================================================================
    
    std::vector<std::vector<int> > newFaces(nTrianglesOUT);
    std::vector<int> faceRefs(nTrianglesOUT);
    std::vector<int> requiredFaces(nTrianglesOUT,0);
    FaceSetPointer m_PMMG_Face2RefPointer;
    FaceSetPointer m_PMMG_TraceFacePointer;
    // std::cout << "leftright_trace " << leftright_trace.size() << std::endl;
    int tracef  = 0;
    int nreq2   = 0;
    int *required2  = (int*)calloc(nTrianglesOUT,sizeof(int));
    FaceSetPointer m_PMMG_TraceFace2PrismPointer;
    FaceSetPointer m_PMMG_TraceFace2PrismPointer2Add;
    int pf      = 0;
    int pftot   = 0;
    int pfn     = 0;
    int fount   = 0;
    int found1  = 0;
    std::map<int,std::vector<int> > rank2req_prism;
    std::set<int> collect_prism_ids_toberequested;
    for ( k=0; k<nTrianglesOUT; k++ )
    {
        std::vector<int> face(3,0);
        if ( PMMG_Get_triangle(parmesh,&(face[0]),
                                        &(face[1]),
                                        &(face[2]),
                                        &(faceRefs[k]),
                                        &(required2[k])) != 1 )
        {
            fprintf(inm,"Unable to get mesh triangle %d \n",k);
            ier = PMMG_STRONGFAILURE;
        }

        face[0] = loc2globVid[face[0]];
        face[1] = loc2globVid[face[1]];
        face[2] = loc2globVid[face[2]];
        if(faceRefs[k]!=0 && faceRefs[k]<(bndIDmax+1))
        {
            //std::cout << "faceRefs[k] " << faceRefs[k] << std::endl;

            FaceSharedPtr Face2RefPointer = std::shared_ptr<NekFace>(new NekFace(face));
            std::pair<FaceSetPointer::iterator, bool> testFace2RefPointer;
            testFace2RefPointer = m_PMMG_Face2RefPointer.insert(Face2RefPointer);
            if(testFace2RefPointer.second)
            {
                (*testFace2RefPointer.first)->SetFaceRef(faceRefs[k]);
            }
            found1++;
        }
        if(faceRefs[k] >= (bndIDmax+1))
        {
            FaceSharedPtr TraceFacePointer = std::shared_ptr<NekFace>(new NekFace(face));
            std::pair<FaceSetPointer::iterator, bool> testTraceFacePointer;
            testTraceFacePointer = m_PMMG_TraceFacePointer.insert(TraceFacePointer);
            
            if(testTraceFacePointer.second)
            {
                (*testTraceFacePointer.first)->SetFaceID(faceRefs[k]);
            }

            std::vector<int> refs(3,0);
            for(int v=0;v<3;v++)
            {
                refs[v] = refOUT[glob2locVid[face[v]]-1];
                if(refs[v]==0)
                {
                    std::cout << "WarNING " << std::endl;
                }
            }

            int PrismID_tmp = faceRefs[k]-(bndIDmax+1);
            int PrismID     = -1;
            if(tagE2gE_P.find(PrismID_tmp)!=tagE2gE_P.end())
            {
                int PrismID = tagE2gE_P[PrismID_tmp];

                std::pair<FaceSetPointer::iterator, bool> Trace2PrismPointer;
                Trace2PrismPointer = m_PMMG_TraceFace2PrismPointer.insert(TraceFacePointer);

                if(Trace2PrismPointer.second)
                {
                    (*Trace2PrismPointer.first)->SetFaceRightElement(PrismID);
                }
            }
            else
            {
                collect_prism_ids_toberequested.insert(PrismID_tmp);
                if(world_rank == 0)
                {
                    int pid = part_global_bl[PrismID_tmp];
                    rank2req_prism[pid].push_back(PrismID_tmp);
                }
                
                TraceFacePointer->SetFaceRightElement(PrismID_tmp);

                std::pair<FaceSetPointer::iterator, bool> Trace2PrismPointer;
                Trace2PrismPointer = m_PMMG_TraceFace2PrismPointer2Add.insert(TraceFacePointer);
            }

            // FaceSharedPtr TraceFaceRefPointer = std::shared_ptr<NekFace>(new NekFace(refs));
            // FaceSetPointer::iterator RefsOnTracePointer = FaceTraceRefs.find(TraceFaceRefPointer);
            // // Mechanism to match up the shell faces from the adapted tetrahedra to the fixed prisms.

            // if(RefsOnTracePointer!=FaceTraceRefs.end())
            // {
            //     // int PrismID = (*RefsOnTracePointer)->GetFaceRightElement();

            //     int PrismID_tmp = (*RefsOnTracePointer)->GetFaceRightElement();
            //     int PrismID     = tag2glob_prisms_on_trace_glob[PrismID_tmp];
        
            //     std::pair<FaceSetPointer::iterator, bool> Trace2PrismPointer;
            //     Trace2PrismPointer = m_PMMG_TraceFace2PrismPointer.insert(TraceFacePointer);

            //     if(Trace2PrismPointer.second)
            //     {
            //         (*Trace2PrismPointer.first)->SetFaceRightElement(PrismID);
            //     }
            // }
            //found++;
        }

        newFaces[k] = face;

        if ( required2 && required2[k] )  nreq2++;
    }

    int nreq;

    if(world_rank != 0)
    {
        std::vector<int> collect_prism_ids_toberequested_vec(collect_prism_ids_toberequested.size(),0);
        //std::cout << "collect_prism_ids_toberequested " << collect_prism_ids_toberequested.size() << " " << world_rank << std::endl;
        std::copy(collect_prism_ids_toberequested.begin(), 
                    collect_prism_ids_toberequested.end(), 
                    collect_prism_ids_toberequested_vec.begin());

        nreq = collect_prism_ids_toberequested_vec.size();
        MPI_Send(&nreq, 1, MPI_INT, 0, world_rank*1000, comm);
        MPI_Send(&collect_prism_ids_toberequested_vec.data()[0], nreq, MPI_INT, 0, world_rank*10000, comm);

        std::vector<int> pid_2recvback(nreq,0);
        MPI_Recv(&pid_2recvback[0], nreq, MPI_INT, 0, world_rank*20000, comm, MPI_STATUS_IGNORE);

        for(int j=0;j<nreq;j++)
        {
            rank2req_prism[pid_2recvback[j]].push_back(collect_prism_ids_toberequested_vec[j]);
        }

        // std::map<int,std::vector<int> >::iterator itmmm;
        // for(itmmm=rank2req_prism.begin();itmmm!=rank2req_prism.end();itmmm++)
        // {
        //     std::cout << world_rank << " " <<  itmmm->first << " " << itmmm->second.size() << std::endl;
        // }

        
    }
    if(world_rank == 0)
    {
        
        std::vector<int> nrecv(world_size-1,0);
        int accul = 0;
        for(int i=1;i<world_size;i++)
        {
            int dest = i*1000;
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, i, dest, comm, MPI_STATUS_IGNORE);
            nrecv[i-1] = n_reqstd_ids;
            accul = accul + n_reqstd_ids;
        }

        for(int i=1;i<world_size;i++)
        {
            int n_reqstd_ids = nrecv[i-1];
            std::vector<int> nreq_prismids(n_reqstd_ids,0);
            std::vector<int> pid_2sendback(n_reqstd_ids,0);
            MPI_Recv(&nreq_prismids[0], n_reqstd_ids, MPI_INT, i, i*10000, comm, MPI_STATUS_IGNORE);

            for(int j=0;j<n_reqstd_ids;j++)
            {
                pid_2sendback[j] = part_global_bl[nreq_prismids[j]];
            }

            MPI_Send(&pid_2sendback.data()[0], n_reqstd_ids, MPI_INT, i, i*20000, comm);

        }        
    }
    
    
    


    std::map<int,std::vector<int> >::iterator itje;
    ScheduleObj* part_schedule = DoScheduling(rank2req_prism,comm);
    std::map<int,std::vector<int> >  reqstd_prisms_per_rank;
    
    for(int q=0;q<world_size;q++)
    {
        if(world_rank==q)
        {
            int i=0;
            for (itje = rank2req_prism.begin(); itje != rank2req_prism.end(); itje++)
            {
                int n_req           = itje->second.size();
                int dest            = itje->first;
                MPI_Send(&n_req, 1, MPI_INT, dest, 9876+10*dest, comm);
                MPI_Send(&itje->second.data()[0], n_req, MPI_INT, dest, 9876*2+dest*2, comm);
                
                i++;
            }
        }
        else if (part_schedule->SendFromRank2Rank[q].find( world_rank ) != part_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 9876+10*world_rank, comm, MPI_STATUS_IGNORE);
            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 9876*2+world_rank*2, comm, MPI_STATUS_IGNORE);
            reqstd_prisms_per_rank[q] = recv_reqstd_ids;
        }
    }

    std::map<int,std::vector<int> > recv_back_newprismids;
    for(int q=0;q<world_size;q++)
    {
        if(world_rank == q)
        {
            for (itje = reqstd_prisms_per_rank.begin(); itje != reqstd_prisms_per_rank.end(); itje++)
            {
                int nv_send = itje->second.size();
                int dest    = itje->first;
                std::vector<int> newprismids_send(nv_send,0);
                for(int u=0;u<itje->second.size();u++)
                {
                    newprismids_send[u] = tagE2gE_P[itje->second[u]];
                }
                
                MPI_Send(&nv_send, 1, MPI_INT, dest, 9876+1000*dest, comm);
                // MPI_Send(&vert_send[0], nv_send, MPI_DOUBLE, dest, 9876+dest*888, comm);
            
                MPI_Send(&newprismids_send.data()[0], nv_send, MPI_INT, dest, 9876+dest*8888, comm);
                //MPI_Send(&it->second.data()[0], it->second.size(), MPI_INT, dest, 8888*9876+dest*8888,comm);
                
                //delete[] vert_send;
            }
        }
        if(part_schedule->RecvRankFromRank[q].find( world_rank ) != part_schedule->RecvRankFromRank[q].end())
         {
            int n_recv_back;
            MPI_Recv(&n_recv_back, 1, MPI_INT, q, 9876+1000*world_rank, comm, MPI_STATUS_IGNORE);
            
            std::vector<int> recv_back_arr(n_recv_back);
            MPI_Recv(&recv_back_arr.data()[0], n_recv_back, MPI_INT, q, 9876+world_rank*8888, comm, MPI_STATUS_IGNORE);

            recv_back_newprismids[q] = recv_back_arr;
         }
    }
    std::map<int,int> old_prism_new_prism;
    for(itje=recv_back_newprismids.begin();itje!=recv_back_newprismids.end();itje++)
    {
        int pid     =  itje->first;
        int nrecv   = itje->second.size();

        for(int q=0;q<nrecv;q++)
        {
            int prismID    = rank2req_prism[pid][q];
            int prismIDnew = itje->second[q];
            old_prism_new_prism[prismID] = prismIDnew;
        }
    }
   

    for(ftit=m_PMMG_TraceFace2PrismPointer2Add.begin();ftit!=m_PMMG_TraceFace2PrismPointer2Add.end();ftit++)
    {

        FaceSharedPtr facenew = *ftit;

        int prismID     = (*ftit)->GetFaceRightElement();
        int prismID_new = old_prism_new_prism[prismID];
        (*ftit)->SetFaceRightElement(prismID_new);

        std::pair<FaceSetPointer::iterator, bool> Trace2PrismPointer;
        Trace2PrismPointer = m_PMMG_TraceFace2PrismPointer.insert(facenew);        
    }
    
    int **out_tria_loc, **out_node_loc;
    int *nitem_face_comm,*nitem_node_comm;
    int next_face_comm, next_node_comm;
    int *color_node_out,*color_face_out;
    
    ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&next_node_comm);
    
    color_node_out  = (int *) malloc(next_node_comm*sizeof(int));
    nitem_node_comm = (int *) malloc(next_node_comm*sizeof(int));
    for( int icomm=0; icomm<next_node_comm; icomm++ )
        ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                                &color_node_out[icomm],
                                                &nitem_node_comm[icomm]);
    
    // Get IDs of nodes on each interface //
    out_node_loc = (int **) malloc(next_node_comm*sizeof(int *));
    for( int icomm=0; icomm<next_node_comm; icomm++ )
        out_node_loc[icomm] = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
    ier = PMMG_Get_NodeCommunicator_nodes(parmesh, out_node_loc);
    
    //===================================================================
    ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
    color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
    nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
    for( int icomm=0; icomm<next_face_comm; icomm++ )
        ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                &color_face_out[icomm],
                                                &nitem_face_comm[icomm]);
    
    
    
    out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
    for( int icomm=0; icomm<next_face_comm; icomm++ )
    {
        out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
    //   std::cout << "nitem_face_comm[icomm] " << icomm << " " << nitem_face_comm[icomm] << " " << world_rank << std::endl;  
    }
    ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
    
    int nPartFace       = 0;
    int nTshared_owned  = 0;
    std::map<int,int> rank2icomm;
    std::set<int> PMMG_SharedVerticesOwned_set;
    std::vector<int> PMMG_SharedVerticesOwned;
    std::vector<int> PMMG_SharedVerticesOwnedRank;
    std::set<int> PMMG_SharedVertices; // Overal number of shared vertices.
    std::map<int,std::vector<int> > Color_SharedOwned;
    FaceSetPointer m_PMMG_SharedFacePointer;
    FaceSetPointer m_PMMG_OwnedSharedFacePointer;


    // This logic essentially determines th at the 
    for(int icomm = 0; icomm < next_face_comm; icomm++ )
    {
        nPartFace = nPartFace + nitem_face_comm[icomm];
        rank2icomm[color_face_out[icomm]] = icomm;
        
        if(color_face_out[icomm] > world_rank)
        {
            nTshared_owned = nTshared_owned + nitem_face_comm[icomm];
            
            for( i = 0; i < nitem_face_comm[icomm]; i++ )
            {
                int ft       = out_tria_loc[icomm][i];
                int face_ref = faceRefs[ft-1];
                //std::cout << "face_refs " << face_ref << std::endl;
                //Color_SharedOwned[icomm].push_back(i);
                std::vector<int> sharedFace(3,0);
                for(int k=0;k<3;k++)
                {
                    int gvt  = newFaces[ft-1][k];
                    // int gvt = loc2globVid[vt];
                    sharedFace[k] = gvt;
                    
                    if(PMMG_SharedVertices.find(gvt)==PMMG_SharedVertices.end())
                    {
                        PMMG_SharedVertices.insert(gvt);
                    }
                    
                    if(PMMG_SharedVerticesOwned_set.find(gvt)==PMMG_SharedVerticesOwned_set.end())
                    {
                        PMMG_SharedVerticesOwned_set.insert(gvt);
                        PMMG_SharedVerticesOwned.push_back(gvt);
                        PMMG_SharedVerticesOwnedRank.push_back(world_rank);
                    }
                }
                
                Color_SharedOwned[color_face_out[icomm]].push_back(i);
                
                FaceSharedPtr sharedFacePointer = std::shared_ptr<NekFace>(new NekFace(sharedFace));
                std::pair<FaceSetPointer::iterator, bool> SharedFPointer;
                SharedFPointer      = m_PMMG_SharedFacePointer.insert(sharedFacePointer);
                std::pair<FaceSetPointer::iterator, bool> OwnedSharedFPointer;
                OwnedSharedFPointer = m_PMMG_OwnedSharedFacePointer.insert(sharedFacePointer);
                
                if(SharedFPointer.second)
                {
                    (*SharedFPointer.first)->SetFaceID(ft);
                }
                
                if(OwnedSharedFPointer.second)
                {
                    (*OwnedSharedFPointer.first)->SetFaceID(ft);
                }                    
            }
        }
        else
        {
            for( i = 0; i < nitem_face_comm[icomm]; i++ )
            {
                int ft       = out_tria_loc[icomm][i];
                
                std::vector<int> sharedFace(3,0);
                for(int k=0;k<3;k++)
                {
                    int gvt  = newFaces[ft-1][k];
                    // int gvt = loc2globVid[vt];
                    sharedFace[k]=gvt;
                    if(PMMG_SharedVertices.find(gvt)==PMMG_SharedVertices.end())
                    {
                        PMMG_SharedVertices.insert(gvt);
                    }
                }
                
                FaceSharedPtr sharedFacePointer = std::shared_ptr<NekFace>(new NekFace(sharedFace));
                std::pair<FaceSetPointer::iterator, bool> testInsPointer;
                testInsPointer = m_PMMG_SharedFacePointer.insert(sharedFacePointer);
                if(testInsPointer.second)
                {
                    (*testInsPointer.first)->SetFaceID(ft);
                }
            }
        }
    }

    int nLocallyOwnedVerts = PMMG_SharedVerticesOwned.size();

    DistributedParallelState* locallyOwnedVerts = new DistributedParallelState(nLocallyOwnedVerts,comm);
    std::vector<int> OwnedSharedVertices_distri(locallyOwnedVerts->getNel(),0);
    std::vector<int> OwnedSharedVertices_distriRank(locallyOwnedVerts->getNel(),0);


    MPI_Allgatherv(&PMMG_SharedVerticesOwned.data()[0],
                    nLocallyOwnedVerts,
                    MPI_INT,
                    &OwnedSharedVertices_distri.data()[0],
                    locallyOwnedVerts->getNlocs(),
                    locallyOwnedVerts->getOffsets(),
                    MPI_INT, comm);
    
    MPI_Allgatherv(&PMMG_SharedVerticesOwnedRank.data()[0],
                    nLocallyOwnedVerts,
                    MPI_INT,
                    &OwnedSharedVertices_distriRank.data()[0],
                    locallyOwnedVerts->getNlocs(),
                    locallyOwnedVerts->getOffsets(),
                    MPI_INT, comm);


    std::map<int,int> ActualOwnedVertDistr;
    std::map<int,std::vector<int> > ActualOwnedVertDistr_map;
    
    for(int i=0;i<locallyOwnedVerts->getNel();i++)
    {
        int gvd = OwnedSharedVertices_distri[i];
        
        if(ActualOwnedVertDistr.find(gvd)==ActualOwnedVertDistr.end())
        {
            ActualOwnedVertDistr_map[gvd].push_back(OwnedSharedVertices_distriRank[i]);

            ActualOwnedVertDistr[gvd] = OwnedSharedVertices_distriRank[i];
        }
        else
        {
            if(OwnedSharedVertices_distriRank[i]<ActualOwnedVertDistr[gvd])
            {
                ActualOwnedVertDistr[gvd] = OwnedSharedVertices_distriRank[i];
            }
            else
            {
                ActualOwnedVertDistr[gvd] = OwnedSharedVertices_distriRank[i];
            }
        }
    }

    std::map<int,std::vector<int> > ActualOwnedVertDistr_map_update;
    std::vector<int> tells(world_size,0);
    std::map<int,std::vector<int> >::iterator itm;
    std::map<int,int>::iterator iitm;
    int tel = 0;
    std::map<int,int> LocationSharedVert;
    for(itm=ActualOwnedVertDistr_map.begin();itm!=ActualOwnedVertDistr_map.end();itm++)
    {
        int gv      = itm->first;
        int ra      = *min_element(itm->second.begin(), itm->second.end());
        tells[ra]   = tells[ra]+1;
        ActualOwnedVertDistr_map_update[ra].push_back(gv);
        LocationSharedVert[gv]=ra;
        tel++;
    }


    //===================================================================
    //Collect Elements  ->      (tetrahedra)
    //===================================================================
    // Read all the adapted tetrahedra.
    // std::vector<std::vector<int> > new_tetrahedra(nTetrahedraOUT);
    // new_tetrahedra.resize(nTetrahedraOUT);
    std::vector<int> refTET(nTetrahedraOUT,0);
    // std::vector<int> required(nTetrahedraOUT,0);
    int tetra_faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};

    // Initialize counters and data structures that will be filled in the coming loop.
    int fid = 0;
    FaceSetPointer m_FaceSetPointer;
    int fv0,fv1,fv2;        
    std::map<int,int> new_locSharedF2globSharedF;
    // Determine the global element numbers of the adapted tetrahedra on this rank.
    DistributedParallelState* distTetraOut = new DistributedParallelState(nTetrahedraOUT,comm);
    int* TetraOUT_offsets = distTetraOut->getOffsets();

    int nPrismsTot = nElemsGlob_P;
    int curElID = TetraOUT_offsets[world_rank]+1+nPrismsTot;

    std::map<int,int> lhsh;
    std::map<int,std::vector<int> > InternalFaces_T;
    std::map<int,std::vector<int> > SharedFaces_T;
    std::map<int,std::vector<int> > TraceFaces_T;
    std::map<int,int> lh_T;
    std::map<int,int> rh_T;
    std::map<int,int> lhsh_owned;
    std::map<int,int> rhsh_owned;
    std::map<int,int> lhtrace;
    std::map<int,int> rhtrace;
    std::map<int,int> NonSharedVrts;
    
    int ngv = 0;
    int cnt = 0;
    // std::map<int,int> new_globE2locE;
    std::map<int,int> tagV2traceVref;
    std::map<int,int> traceref2newtag;
    std::vector<double> Vface(3);
    std::vector<double> r0(3);
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    int negit  = 0;
    int inhere = 0;

    int larg = 0;
    int largf = 0;
    int larg2 = 0;
    int larg2f = 0;
    int largtot = 0;
    int larg3 = 0;
    int largf3 = 0;
    for ( k=0; k<nTetrahedraOUT; k++ )  
    {
        std::vector<int> vertids(4,0);

        if ( PMMG_Get_tetrahedron(parmesh,
                                &(vertids[0]),&(vertids[1]),
                                &(vertids[2]),&(vertids[3]),
                                &(refTET[k]),&(required[k])) != 1 )
        {
            fprintf(inm,"Unable to get mesh tetra %d \n",k);
            ier = PMMG_STRONGFAILURE;
        }

        new_tetrahedra.push_back(vertids);
        new_globE2locE[curElID] = k;

        std::vector<double> P(4*3);
        P[0*3+0]=new_vertices[(vertids[0]-1)][0];
        P[0*3+1]=new_vertices[(vertids[0]-1)][1];
        P[0*3+2]=new_vertices[(vertids[0]-1)][2];
        
        P[1*3+0]=new_vertices[(vertids[1]-1)][0];
        P[1*3+1]=new_vertices[(vertids[1]-1)][1];
        P[1*3+2]=new_vertices[(vertids[1]-1)][2];
        
        P[2*3+0]=new_vertices[(vertids[2]-1)][0];
        P[2*3+1]=new_vertices[(vertids[2]-1)][1];
        P[2*3+2]=new_vertices[(vertids[2]-1)][2];
        
        P[3*3+0]=new_vertices[(vertids[3]-1)][0];
        P[3*3+1]=new_vertices[(vertids[3]-1)][1];
        P[3*3+2]=new_vertices[(vertids[3]-1)][2];
        
        double Vtet = GetQualityTetrahedra(P);
        
        if(Vtet<0.0)
        {
            std::cout << " negative volume in Element " << k << " on rank " << world_rank  <<std::endl;
        }
        std::vector<double> vCenter = ComputeCentroidCoord(P, 4);

        
        for(int u=0;u<4;u++)
        {

            std::vector<int> face(3,0);


            face[0] = loc2globVid[vertids[tetra_faces[u][0]]];
            face[1] = loc2globVid[vertids[tetra_faces[u][1]]];
            face[2] = loc2globVid[vertids[tetra_faces[u][2]]];

            FaceSharedPtr facePointer = std::shared_ptr<NekFace>(new NekFace(face));
            std::pair<FaceSetPointer::iterator, bool> testInsPointer;
            testInsPointer = m_FaceSetPointer.insert(facePointer);

            if(testInsPointer.second)
            {


                double v0x = new_vertices[vertids[tetra_faces[u][0]]-1][0];
                double v0y = new_vertices[vertids[tetra_faces[u][0]]-1][1];
                double v0z = new_vertices[vertids[tetra_faces[u][0]]-1][2];

                double v1x = new_vertices[vertids[tetra_faces[u][1]]-1][0];
                double v1y = new_vertices[vertids[tetra_faces[u][1]]-1][1];
                double v1z = new_vertices[vertids[tetra_faces[u][1]]-1][2];

                double v2x = new_vertices[vertids[tetra_faces[u][2]]-1][0];
                double v2y = new_vertices[vertids[tetra_faces[u][2]]-1][1];
                double v2z = new_vertices[vertids[tetra_faces[u][2]]-1][2];

                Vface[0] = (v0x+v1x+v2x)/3.0;
                Vface[1] = (v0y+v1y+v2y)/3.0;
                Vface[2] = (v0z+v1z+v2z)/3.0;
                
                r0[0] = (Vface[0]-vCenter[0]);///Lr;
                r0[1] = (Vface[1]-vCenter[1]);///Lr;
                r0[2] = (Vface[2]-vCenter[2]);///Lr;
                
                v0[0] = v1x-v0x;
                v0[1] = v1y-v0y;
                v0[2] = v1z-v0z;

                v1[0] = v2x-v0x;
                v1[1] = v2y-v0y;
                v1[2] = v2z-v0z;
                
                std::vector<double> n0        = ComputeSurfaceNormal(v0,v1);
                double orient0   = DotVec3D(r0,n0);
                if(orient0<0.0)
                {
                    negit++;
                }
                
                (*testInsPointer.first)->SetFaceID(fid);

                FaceSetPointer::iterator SharedFPointer         = m_PMMG_SharedFacePointer.find(facePointer);
                FaceSetPointer::iterator OwnedSharedFPointer    = m_PMMG_OwnedSharedFacePointer.find(facePointer);
                FaceSetPointer::iterator testTracePointer       = m_PMMG_TraceFacePointer.find(facePointer);
                FaceSetPointer::iterator testFace2RefPointer    = m_PMMG_Face2RefPointer.find(facePointer);
                FaceSetPointer::iterator testTracePrismPointer  = m_PMMG_TraceFace2PrismPointer.find(facePointer);


                if(SharedFPointer != m_PMMG_SharedFacePointer.end())
                {
                    int lshf                             = (*SharedFPointer)->GetFaceID();
                    new_locSharedF2globSharedF[lshf]     = fid;
                    lhsh[fid]                            = curElID;
                }




                if(SharedFPointer == m_PMMG_SharedFacePointer.end()
                    && testFace2RefPointer == m_PMMG_Face2RefPointer.end()
                    && testTracePointer == m_PMMG_TraceFacePointer.end())
                {
                    InternalFaces_T[fid] = face;
                    //fmInt[fid]  = fce;
                    lh_T[fid]           = curElID;
                    
                    for(int s=0;s<3;s++)
                    {
                        int gvm2 = face[s];
                        if(NonSharedVrts.find(gvm2)==NonSharedVrts.end() &&
                                PMMG_SharedVertices.find(gvm2)==PMMG_SharedVertices.end())
                        {
                            NonSharedVrts[gvm2] = ngv;
                            ngv++;
                        }
                    }
                }
                


                if(OwnedSharedFPointer != m_PMMG_OwnedSharedFacePointer.end())
                {
                    SharedFaces_T[fid]            = face;
                    lhsh_owned[fid]             = curElID;
                }
                
                
                if(testTracePointer != m_PMMG_TraceFacePointer.end())
                {
                    int FaceRef          = (*testTracePointer)->GetFaceRef();
                    TraceFaces_T[fid]    = face;
                    lhtrace[fid]         = curElID;

                    if(testTracePrismPointer != m_PMMG_TraceFace2PrismPointer.end())
                    {
                        rhtrace[fid] = (*testTracePrismPointer)->GetFaceRightElement();

                        for(int y=0;y<3;y++)
                        {
                            int traceVertRef = refOUT[vertids[tetra_faces[u][y]]-1];
                            int vtag         = face[y];

                            //std::cout << "traceVertRef " << traceVertRef << std::endl;
                            
                            if(tagV2traceVref.find(vtag)==tagV2traceVref.end())
                            {
                                tagV2traceVref[vtag]=traceVertRef;
                                traceref2newtag[traceVertRef]=vtag;
                            }
                        }

                        inhere++;
                    }
                }
                
                
                if(testFace2RefPointer != m_PMMG_Face2RefPointer.end())
                {
                    int FaceRef          = (*testFace2RefPointer)->GetFaceRef();
                    BoundaryFaces_T[fid] = face;
                    lh_T_bc[fid]         = curElID;
                    bcmap[FaceRef].push_back(fid);
                }
            }
            else
            {
                int fid_n       = (*testInsPointer.first)->GetFaceID();
                rh_T[fid_n]     = curElID;
            }


            fid++;
        }

        curElID++;
    }


    
    std::map<int,int> tagV2traceVref_glob  = AllGatherMap_T(tagV2traceVref,comm);
    //std::map<int,int> traceref2newtag_glob = AllGatherMap_T(traceref2newtag,comm);

    // std::cout << "===============================================" << std::endl;
    // std::cout << "BoundaryFaces_T " << BoundaryFaces_T.size() << " " << world_rank << std::endl;
    // std::cout << "TraceFaces_T " << TraceFaces_T.size() << " " << world_rank << std::endl;
    // std::cout << "SharedFaces_T " << SharedFaces_T.size() << " " << world_rank << std::endl;
    // std::cout << "InternalFaces_T " << InternalFaces_T.size() << " " << world_rank << std::endl;
    // std::cout << "===============================================" << std::endl;

    PMMG_Free_all(PMMG_ARG_start,
            PMMG_ARG_ppParMesh,&parmesh,
            PMMG_ARG_end);



    nLocIntVrts             = NonSharedVrts.size();
    nLocShVrts              = ActualOwnedVertDistr_map_update[world_rank].size();


    int nLocTotVrts         = nLocIntVrts+nLocShVrts;
    DistributedParallelState* distnLocTotVrts = new DistributedParallelState(nLocTotVrts,comm);

    int vert                = 0;//distnLocTotVrts->getOffsets()[world_rank];
    int ToTVrts             =  distnLocTotVrts->getNel();
    int* ToTVrts_offsets    =  distnLocTotVrts->getOffsets();
    int ToTVrts_offset      =  ToTVrts_offsets[world_rank];

    std::map<int,int> newglob2oldglob;
    
    int nPrismVerts_tmp = nVertsGlob_P;//distPrismIntVerts->getNel()+distPrismShaVerts->getNel();

    vertices_output.resize(nLocTotVrts*3);

    std::fill(vertices_output.begin(), vertices_output.end(), 0.0);
    for(iitm=NonSharedVrts.begin();iitm!=NonSharedVrts.end();iitm++)
    {
        int globvid = iitm->first;

        vertices_output[vert*3+0] = new_vertices[glob2locVid[globvid] - 1][0];
        vertices_output[vert*3+1] = new_vertices[glob2locVid[globvid] - 1][1];
        vertices_output[vert*3+2] = new_vertices[glob2locVid[globvid] - 1][2];

        oldglob2newglob[globvid] = vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp;
        newglob2oldglob[vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp] = globvid;
        
        vert++;
    }
        
    for(int i=0;i<nLocShVrts;i++)
    {
        int globvid = ActualOwnedVertDistr_map_update[world_rank][i];

        vertices_output[vert*3+0] = new_vertices[glob2locVid[globvid] - 1][0];
        vertices_output[vert*3+1] = new_vertices[glob2locVid[globvid] - 1][1];
        vertices_output[vert*3+2] = new_vertices[glob2locVid[globvid] - 1][2];

        oldglob2newglob[globvid] = vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp;
        newglob2oldglob[vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp] = globvid;

        vert++;
    }


    // The vertex IDs have just been updated so that they correspond to a global mesh definition.
    // The new global vertex indexing, that has been done above, includes the prisms etc.
    // However, the vertices that lie on the shared interfaces are counted double when defining the new vertex index number.
    // This is why, using the same logic as before, we need to go through all shared vertices, register their rank,
    // and determine that the lowest rank will own that particular shared vertex ID.
    // Perhaps revisit this logic and see if it can be done in an easier/more direct way.

    std::vector<int> updateGlobSharedVrtID(LocationSharedVert.size(),0);
    std::vector<int> updateGlobSharedVrtID_red(LocationSharedVert.size(),0);
    std::vector<int> originalGlobSharedVrtID(LocationSharedVert.size(),0);
    std::vector<int> originalGlobSharedVrtID_red(LocationSharedVert.size(),0);
    
    int ig = 0;

    for(iitm=LocationSharedVert.begin();iitm!=LocationSharedVert.end();iitm++)
    {
        int old_globid = iitm->first;
        int ra  = iitm->second;
        
        originalGlobSharedVrtID_red[ig] = 0;
        updateGlobSharedVrtID_red[ig]   = 0;
        
        if(ra == world_rank)
        {
            originalGlobSharedVrtID[ig] = old_globid;
            updateGlobSharedVrtID[ig]   = oldglob2newglob[old_globid];
        }
        else
        {
            originalGlobSharedVrtID[ig] = 0;
            updateGlobSharedVrtID[ig]   = 0;
        }
        ig++;
    }

    MPI_Allreduce(&originalGlobSharedVrtID.data()[0],
                    &originalGlobSharedVrtID_red.data()[0],
                    LocationSharedVert.size(),
                    MPI_INT, MPI_SUM, comm);
    
    MPI_Allreduce(&updateGlobSharedVrtID.data()[0],
                    &updateGlobSharedVrtID_red.data()[0],
                    LocationSharedVert.size(),
                    MPI_INT, MPI_SUM, comm);
    
    
    for(int i=0;i<LocationSharedVert.size();i++)
    {
        LocationSharedVert_update[originalGlobSharedVrtID_red[i]] = updateGlobSharedVrtID_red[i];
    }
    
    // Generate datastructures to output
    int nLocFace = InternalFaces_T.size()+SharedFaces_T.size()+TraceFaces_T.size();
    // std::map<int,std::vector<int> >::iterator itmiv;
    // int* ifn_T = new int[nLocFace*8]; 
    ifn_T.resize(nLocFace*8);
    std::fill(ifn_T.begin(), ifn_T.end(), 0);

    std::vector<double> VcF(3,0.0);
    std::vector<double> Vijk(3,0.0);
    std::vector<int> fq(3,0);
    
    int ftot = 0;
    int nfound = 0;
    for(itmiv=InternalFaces_T.begin();itmiv!=InternalFaces_T.end();itmiv++)
    {
        int fid = itmiv->first;

        // set the type of face. In case of a triagle type=3;
        int type = 3;
        ifn_T[ftot*8+0] = type;
        // The face orientation needs to be determined in order to make sure that
        // the face is stored in the global array in the appropriate way.
        VcF[0] = 0.0;
        VcF[1] = 0.0;
        VcF[2] = 0.0;
        
        std::vector<std::vector<double> > Vfaces;
        for(int i=0;i<itmiv->second.size();i++)
        {
            int gvert = itmiv->second[i];
            int lvert = glob2locVid[gvert];

            std::vector<double> Vf(3,0.0);

            Vf[0]  = new_vertices[lvert-1][0];
            Vf[1]  = new_vertices[lvert-1][1];
            Vf[2]  = new_vertices[lvert-1][2];
            
            VcF[0] = VcF[0] + Vf[0];
            VcF[1] = VcF[1] + Vf[1];
            VcF[2] = VcF[2] + Vf[2];

            Vfaces.push_back(Vf);

            if(LocationSharedVert_update.find(gvert)!=LocationSharedVert_update.end())
            {
                fq[i] = LocationSharedVert_update[gvert];
                if(fq[i] == 0)
                {
                    std::cout << "tracetag2glob Internal " << world_rank << " " << i << " " << std::endl;
                }
            }
            else
            {
                fq[i] = oldglob2newglob[gvert];
                if(fq[i] == 0)
                {
                    std::cout << "tracetag2glob Internal " << world_rank << " " << i << " " << std::endl;
                }
            }
        }

        VcF[0] = VcF[0]/3.0;
        VcF[1] = VcF[1]/3.0;
        VcF[2] = VcF[2]/3.0;

        int gv0 = fq[0];
        int gv1 = fq[1];
        int gv2 = fq[2];
        
        Vijk[0] = 0.0;
        Vijk[1] = 0.0;
        Vijk[2] = 0.0;

        int curElID = lh_T[fid];
        int locElID = new_globE2locE[curElID];

        for(int u=0;u<new_tetrahedra[locElID].size();u++)
        {
            Vijk[0] = Vijk[0] + new_vertices[new_tetrahedra[locElID][u]-1][0];
            Vijk[1] = Vijk[1] + new_vertices[new_tetrahedra[locElID][u]-1][1];
            Vijk[2] = Vijk[2] + new_vertices[new_tetrahedra[locElID][u]-1][2];  
        }
        
        Vijk[0] = Vijk[0]/new_tetrahedra[locElID].size();
        Vijk[1] = Vijk[1]/new_tetrahedra[locElID].size();
        Vijk[2] = Vijk[2]/new_tetrahedra[locElID].size();

        double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);

        // Storing the face appropriately depending on its original orientation.
        if(orient0 < 0.0)
        {
            ifn_T[ftot*8+1]=gv0;
            ifn_T[ftot*8+2]=gv2;
            ifn_T[ftot*8+3]=gv1;
            ifn_T[ftot*8+4]=0;
        }
        else
        {
            ifn_T[ftot*8+1]=gv0;
            ifn_T[ftot*8+2]=gv1;
            ifn_T[ftot*8+3]=gv2;
            ifn_T[ftot*8+4]=0;
        }
        
        ifn_T[ftot*8+5]=rh_T[fid];
        ifn_T[ftot*8+6]=lh_T[fid];
        ifn_T[ftot*8+7]=2;

        ftot++;
    }
    //===================================================================
    
    // std::map<int,int> locShF2globShF = tetra_repart->GetLocalSharedFace2GlobalSharedFace();
    // Preparation to store the SharedFaces to the general array ifn_T
    ScheduleObj* ish_schedule = DoScheduling(Color_SharedOwned,comm);
    std::map<int,std::vector<int> > recv_ids;
    std::map<int,std::vector<int> >::iterator it;
        
    for(int q=0;q<world_size;q++)
    {
        if(world_rank==q)
        {
            int i=0;
            for (it = Color_SharedOwned.begin(); it != Color_SharedOwned.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
                MPI_Send(&it->second[0], n_req, MPI_INT, dest, 14876+dest, comm);
                i++;
            }
        }
        else if (ish_schedule->SendFromRank2Rank[q].find( world_rank ) != ish_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*world_rank, comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 14876+world_rank, comm, MPI_STATUS_IGNORE);
            recv_ids[q] = recv_reqstd_ids;
        }
    }
    
    std::map<int,std::vector<int> > sendEl;
    std::map<int,std::vector<int> >::iterator rcvit;
    
    for(rcvit=recv_ids.begin();rcvit!=recv_ids.end();rcvit++)
    {
        int frank = rcvit->first;
        
        int icomm = rank2icomm[frank];
        int nF    = rcvit->second.size();
        
        for(int j=0;j<nF;j++)
        {
            int fidInt = rcvit->second[j];
            int ofid   = out_tria_loc[icomm][fidInt];
            int nfid   = new_locSharedF2globSharedF[ofid];
            sendEl[frank].push_back(lhsh[nfid]);
            
            if(lhsh[nfid] == 0)
            {
                std::cout << "yep commi not correct " << std::endl;
            }
        }
    }
    
    ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);

    std::map<int,std::vector<int> > adj_ids;
    for(int q=0;q<world_size;q++)
    {
        if(world_rank==q)
        {
            int i=0;
            for (it = sendEl.begin(); it != sendEl.end(); it++)
            {
                int n_req           = it->second.size();
                int dest            = it->first;

                MPI_Send(&n_req, 1,
                        MPI_INT, dest,
                        6798+78000*dest, comm);
                
                MPI_Send(&it->second[0],
                        n_req, MPI_INT,
                        dest, 14876000+dest, comm);

                i++;
            }
        }
        else if (ishBack_schedule->SendFromRank2Rank[q].find( world_rank ) != ishBack_schedule->SendFromRank2Rank[q].end())
        {
            int n_reqstd_ids;
            
            MPI_Recv(&n_reqstd_ids,
                    1, MPI_INT, q,
                    6798+78000*world_rank,
                    comm, MPI_STATUS_IGNORE);

            std::vector<int> recv_reqstd_ids(n_reqstd_ids);
            
            MPI_Recv(&recv_reqstd_ids[0],
                    n_reqstd_ids,
                    MPI_INT, q,
                    14876000+world_rank,
                    comm, MPI_STATUS_IGNORE);
            
            adj_ids[q] = recv_reqstd_ids;

        }
    }
    
    DistributedParallelState* rhbefore = new DistributedParallelState(rh_T.size(),comm);
    DistributedParallelState* lhbefore = new DistributedParallelState(lh_T.size(),comm);

    std::map<int,int> adjElements;
    int fid_loc,fid_glo;
    int telli = 0;
    std::set<int> frh;
    int outside = 0;
    for(itm=Color_SharedOwned.begin();itm!=Color_SharedOwned.end();itm++)
    {
        int rrank = itm->first;
        
        for(int j=0;j<itm->second.size();j++)
        {
            fid_loc              = itm->second[j];
            int icomm            = rank2icomm[rrank];
            int ft               = out_tria_loc[icomm][fid_loc];
            fid_glo              = new_locSharedF2globSharedF[ft];
            rh_T[fid_glo]        = adj_ids[itm->first][j];
            
            if(frh.find(fid_glo)==frh.end())
            {
                adjElements[fid_glo] = adj_ids[itm->first][j];
                
                frh.insert(fid_glo);
                outside++;
            }
            telli++;
        }
    }
    

    int elLh    = -1;
    int elRh    = -1;
    for(itmiv=SharedFaces_T.begin();itmiv!=SharedFaces_T.end();itmiv++)
    {
        fid             = itmiv->first;
        ifn_T[ftot*8+0] = 3;


        std::vector<std::vector<double> > Vfaces;
        VcF[0] = 0.0;
        VcF[1] = 0.0;
        VcF[2] = 0.0;
        for(int q=0;q<3;q++)
        {
            // int lvert = glob2locVid[itmiv->second->GetEdgeIDs()[q]];
            int gvert = itmiv->second[q];
            int lvert = glob2locVid[gvert];

            std::vector<double> Vf(3,0.0);

            Vf[0]  = new_vertices[lvert-1][0];
            Vf[1]  = new_vertices[lvert-1][1];
            Vf[2]  = new_vertices[lvert-1][2];
            
            VcF[0] = VcF[0] + Vf[0];
            VcF[1] = VcF[1] + Vf[1];
            VcF[2] = VcF[2] + Vf[2];

            Vfaces.push_back(Vf);
            
            if(LocationSharedVert_update.find(gvert)!=LocationSharedVert_update.end())
            {
                fq[q] = LocationSharedVert_update[gvert];

                if(fq[q] == 0)
                {
                    std::cout << "tracetag2glob Shared " << world_rank << " " << q << " " << std::endl;
                }
            }
            else
            {
                fq[q] = oldglob2newglob[gvert];
                if(fq[q] == 0)
                {
                    std::cout << "tracetag2glob Shared " << world_rank << " " << q << " " << std::endl;
                }
            }
        }
        
        VcF[0] = VcF[0]/3.0;
        VcF[1] = VcF[1]/3.0;
        VcF[2] = VcF[2]/3.0;
        
        int gv0 = fq[0];
        int gv1 = fq[1];
        int gv2 = fq[2];
        
        if(lhsh_owned.find(fid)!=lhsh_owned.end())
        {
            elLh = lhsh_owned[fid];
            elRh = adjElements[fid];
        }
        
        Vijk[0] = 0.0;
        Vijk[1] = 0.0;
        Vijk[2] = 0.0;
        // compute element center;
        //ienOUT[curElID] = Elvrts

        int locElID = new_globE2locE[elLh];
        for(int u=0;u<new_tetrahedra[locElID].size();u++)
        {
            Vijk[0] = Vijk[0] + new_vertices[new_tetrahedra[locElID][u]-1][0];
            Vijk[1] = Vijk[1] + new_vertices[new_tetrahedra[locElID][u]-1][1];
            Vijk[2] = Vijk[2] + new_vertices[new_tetrahedra[locElID][u]-1][2];  
        }
        
        
        Vijk[0] = Vijk[0]/new_tetrahedra[locElID].size();
        Vijk[1] = Vijk[1]/new_tetrahedra[locElID].size();
        Vijk[2] = Vijk[2]/new_tetrahedra[locElID].size();

        double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
        
        if(orient0 < 0.0)
        {
            ifn_T[ftot*8+1] = gv0;
            ifn_T[ftot*8+2] = gv2;
            ifn_T[ftot*8+3] = gv1;
            ifn_T[ftot*8+4] = 0;
        }
        else
        {
            ifn_T[ftot*8+1] = gv0;
            ifn_T[ftot*8+2] = gv1;
            ifn_T[ftot*8+3] = gv2;
            ifn_T[ftot*8+4] = 0;
        }

        ifn_T[ftot*8+5] = elRh;
        ifn_T[ftot*8+6] = elLh;
        ifn_T[ftot*8+7] = 2;

        
        if(elRh == 0)
        {    
            std::cout <<"Found the face in shared -> " << fid << " " << world_rank << " " << elLh << std::endl;
        } 

        ftot++;
    }    



    std::map<int,int> tracetag2glob;
    int ff = 0;
    int fff = 0;
    for(itmiv=TraceFaces_T.begin();itmiv!=TraceFaces_T.end();itmiv++)
    {
        fid             = itmiv->first;
        ifn_T[ftot*8+0] = 3; 
        
        VcF[0] = 0.0;
        VcF[1] = 0.0;
        VcF[2] = 0.0;
        
        std::vector<std::vector<double> > Vfaces;
        std::vector<int> reference(3);

        for(int q=0;q<3;q++)
        {
            int gvert = itmiv->second[q];
            int lvert = glob2locVid[gvert];

            std::vector<double> Vf(3,0.0);

            Vf[0]           = new_vertices[lvert-1][0];
            Vf[1]           = new_vertices[lvert-1][1];
            Vf[2]           = new_vertices[lvert-1][2];

            reference[q]    = refOUT[(lvert-1)];
            
            VcF[0] = VcF[0] + Vf[0];
            VcF[1] = VcF[1] + Vf[1];
            VcF[2] = VcF[2] + Vf[2];
            
            Vfaces.push_back(Vf);
            
            if(LocationSharedVert_update.find(itmiv->second[q])!=LocationSharedVert_update.end())
            {
                fq[q]           = LocationSharedVert_update[itmiv->second[q]];
                
                int traceRef    = tagV2traceVref_glob[itmiv->second[q]];
                
                if(tracetag2glob.find(traceRef)==tracetag2glob.end())
                {
                    tracetag2glob[traceRef] = fq[q];
                    
                    if(fq[q] == 0)
                    {
                        std::cout << "tracetag2glob Trace " << world_rank << " " << q << " " << traceRef << " " << itmiv->second[q] << std::endl;
                    }
                }
            }
            else
            {
                if(oldglob2newglob.find(itmiv->second[q])==oldglob2newglob.end())
                {
                    std::cout << "itmiv->second[q] " << q << " " << itmiv->second[q] << " " << std::endl;  
                }

                fq[q]           = oldglob2newglob[itmiv->second[q]];

                int traceRef    = tagV2traceVref_glob[itmiv->second[q]];
                
                if(tracetag2glob.find(traceRef)==tracetag2glob.end())
                {
                    tracetag2glob[traceRef] = fq[q];

                    if(fq[q] == 0)
                    {
                        std::cout << "Noot in shared tracetag2glob Trace " << world_rank << " " << q << " " << traceRef << " " << itmiv->second[q] << std::endl;
                    }
                }
            }
        }
        
        VcF[0] = VcF[0]/3.0;
        VcF[1] = VcF[1]/3.0;
        VcF[2] = VcF[2]/3.0;
        
        
        int gv0 = fq[0];
        int gv1 = fq[1];
        int gv2 = fq[2];
        
        if(rhtrace.find(fid)!=rhtrace.end())
        {
            fff++;
            elRh = rhtrace[fid];
            elLh = lhtrace[fid];
        }
        else
        {
            //std::cout << "fid not here " << fid << std::endl;
            elRh = rh_T[fid];
            elLh = lhtrace[fid];
            ff++;
        }
        
        std::vector<double> Vijk(3);
        
        Vijk[0] = 0.0;
        Vijk[1] = 0.0;
        Vijk[2] = 0.0;
        
        // compute element center;
        //ienOUT[curElID] = Elvrts

        int locElID = new_globE2locE[elLh];
        for(int u=0;u<new_tetrahedra[locElID].size();u++)
        {
            Vijk[0] = Vijk[0] + new_vertices[new_tetrahedra[locElID][u]-1][0];
            Vijk[1] = Vijk[1] + new_vertices[new_tetrahedra[locElID][u]-1][1];
            Vijk[2] = Vijk[2] + new_vertices[new_tetrahedra[locElID][u]-1][2];  
        }
        
        Vijk[0] = Vijk[0]/new_tetrahedra[locElID].size();
        Vijk[1] = Vijk[1]/new_tetrahedra[locElID].size();
        Vijk[2] = Vijk[2]/new_tetrahedra[locElID].size();

        double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
        
        if(orient0 < 0.0)
        {
            ifn_T[ftot*8+1] = gv0; 
            ifn_T[ftot*8+2] = gv2; 
            ifn_T[ftot*8+3] = gv1; 
            ifn_T[ftot*8+4] = 0; 
        }
        else
        {
            ifn_T[ftot*8+1] = gv0; 
            ifn_T[ftot*8+2] = gv1; 
            ifn_T[ftot*8+3] = gv2; 
            ifn_T[ftot*8+4] = 0; 
        }
        
        ifn_T[ftot*8+5] = elRh; 
        ifn_T[ftot*8+6] = elLh; 
        ifn_T[ftot*8+7] = 2;

        Vfaces.clear();

        if(elRh == 0)
        {     
            std::cout <<"Found the face in shell -> " << fid << " " << world_rank << " " << elLh << std::endl;
        } 
        
        ftot++;
    }  
    
    tracerefV2globalV = AllGatherMap_T(tracetag2glob,comm);
    
}

void RunParMMGAndTestPartitioning(MPI_Comm comm,
                                  PMMG_pParMesh parmesh, 
                                  PrepareAdaption* prepare,
                                  std::map<int,std::vector<int> > ranges_id,
                                  Inputs* inputs)
{

    FILE *inm;

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    int ier, opt;

    


    std::set<int> Owned_Elem_t                          = prepare->getElemSet();

    // This test assumes that the following two build requests on the tetra_repart object have been called before entering this routine
    // tetra_repart->buildInteriorSharedAndBoundaryFacesMaps(comm,pttrace, ranges_id);
    // tetra_repart->buildTag2GlobalElementMapsAndFace2LeftRightElementMap(comm,pttrace, ranges_id);

    //std::map<int,int> locv2tagvID                       = tetra_repart->getLocalVert2VertTag();
    std::map<int,int> locv2tagvID                       = prepare->getLocalVert2GlobalVert();
    //std::map<int,int> globalv2localvID                  = prepare->getUpdatedGlobal2LocalVMap();
    std::map<int,int> globalv2localvID                  = prepare->getNewGlobalVert2LocalVertMap();
    //std::map<int,int> tag2globalV                     = tetra_repart->getUpdatedTag2GlobalVMap();

    // This test assumes that the following two build requests on the tetra_repart object have been called before entering this routine
    // tetra_repart->buildCommunicationMaps(comm);
    std::vector<int> face4parmmg                        = prepare->getFace4ParMMG(); // checked
    std::map<int,int> global2tagF                       = prepare->getGlobal2TagFMap();
    int** ifc_tria_glob                                 = prepare->getParMMGCommFace2GlobalVertMap();
    int** ifc_tria_loc                                  = prepare->getParMMGCommFace2LocalVertMap();
    int* color_face                                     = prepare->getParMMGCommColorFace();
    int *ntifc                                          = prepare->getParMMGCommNFacesPerColor();
    int ncomm                                           = prepare->getParMMGNComm();

    std::map<int,int> locShF2globShF                    = prepare->getLocalSharedFace2GlobalSharedFace();
    std::map<int,std::vector<int> > face2node           = prepare->getFace2VertNewMap();

    int nVertices   = locv2tagvID.size();
    int nTetrahedra = Owned_Elem_t.size();
    int nEdges      = 0;
    int nTriangles  = face4parmmg.size();

    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, inputs->niter ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hausd, inputs->hausd) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hgrad, inputs->hgrad) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if( !PMMG_Set_dparameter( parmesh,  PMMG_DPARAM_hgradreq , -1.0 ) ){
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
 
    int nVerticesIN   = 0;
    int nTetrahedraIN = 0;
    int nTrianglesIN  = 0;
    int nEdgesIN      = 0;
    
    if ( PMMG_Get_meshSize(parmesh,&nVerticesIN,&nTetrahedraIN,NULL,&nTrianglesIN,NULL,
                           &nEdgesIN) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }
    
    // std::cout << "nVertices input " << nVerticesIN << std::endl;
    // std::cout << "nTetrahedra input " << nTetrahedraIN << std::endl;
    // std::cout << "nTriangles input " << nTrianglesIN << std::endl;
    // std::cout << "nEdges input " << nEdgesIN << std::endl;
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_globalNum, 1 ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };

    int ierlib = PMMG_parmmglib_distributed( parmesh );
    
    if(ierlib==0 && world_rank == 0)
    {
        std::cout << "SUCCESFULLY adapted the mesh in parallel!" << std::endl;
    }
    else if(ierlib==1 && world_rank == 0)
    {
        std::cout << "FAILED to adapt the mesh in parallel! "<< std::endl;
    }

    int nVerticesOUT   = 0;
    int nTetrahedraOUT = 0;
    int nTrianglesOUT  = 0;
    int nEdgesOUT      = 0;


    if ( PMMG_Get_meshSize(parmesh,&nVerticesOUT,&nTetrahedraOUT,NULL,&nTrianglesOUT,NULL,
                           &nEdgesOUT) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }
    
    int             nodeGloNumber,nodeOwner;
    std::map<int,int> loc2globVid;
    std::map<int,int> glob2locVid;
    std::vector<int> globIDs;
    for(int k = 1; k <= nVerticesOUT; k++ )
    {
        if( !PMMG_Get_vertexGloNum( parmesh, &nodeGloNumber, &nodeOwner ) )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    
        loc2globVid[k]=nodeGloNumber;
        glob2locVid[nodeGloNumber]=k;
        globIDs.push_back(nodeGloNumber);
    }
    
    int *required = (int*)calloc(MAX4(nVerticesOUT,
                                      nTetrahedraOUT,
                                      nTrianglesOUT,
                                      nEdgesOUT),sizeof(int));
    
    int *refOUT = (int*)calloc(MAX4(nVerticesOUT,
                                    nTetrahedraOUT,
                                    nTrianglesOUT,
                                    nEdgesOUT),sizeof(int));
    
    int *corner = (int*)calloc(nVerticesOUT,sizeof(int));
    int pos;
    //int *refOUT = new int[nVerticesOUT];

    std::vector<std::vector<int> > outT;
    
    double *vertOUT = (double*)calloc((nVerticesOUT)*3,sizeof(double));
    //std::map<int,std::vector<double> > ref2coordinatesOUT;
    for ( int k=0; k<nVerticesOUT; k++ )
    {
          pos = 3*k;
          if ( PMMG_Get_vertex(parmesh,&(vertOUT[pos]),&(vertOUT[pos+1]),&(vertOUT[pos+2]),
                               &(refOUT[k]),&(corner[k]),&(required[k])) != 1 ) {
            fprintf(inm,"Unable to get mesh vertex %d \n",k);
            ier = PMMG_STRONGFAILURE;
          }
    }

    // RUN TEST TO SEE WHETHER ADAPTATION WAS SUCCESFULL



    if(inputs->niter==0)
    {
        std::cout << "Check the input and outputted shared faces." << std::endl;
        
        int **out_tria_loc;
        int *nitem_face_comm;
        int next_face_comm;
        ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
        int *color_node_out,*color_face_out;
        color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
        nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                 &color_face_out[icomm],
                                                 &nitem_face_comm[icomm]);
        
        out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
        
        // Check matching of input interface nodes with the set ones

        // FaceSharedPtr TestFacePointer = std::shared_ptr<NekFace>(new NekFace(vidstest));
        // FaceSetPointer m_TestFaceSet;

        // Get input triangle nodes
        int** faceNodes2 = (int **) malloc(ncomm*sizeof(int *));
        for( int icomm = 0; icomm < ncomm; icomm++ ) {
          faceNodes2[icomm] = (int *) malloc(3*ntifc[icomm]*sizeof(int));
          for(int i = 0; i < ntifc[icomm]; i++ ) {
              
            int faceID = ifc_tria_loc[icomm][i]-1;
            int faceIDg = ifc_tria_glob[icomm][i]-1;
            int faceID2 = locShF2globShF[faceID];

            int v0 = face2node[faceID2][0];
            int v1 = face2node[faceID2][1];
            int v2 = face2node[faceID2][2];
    
            int v0l = globalv2localvID[v0];
            int v1l = globalv2localvID[v1];
            int v2l = globalv2localvID[v2];

            // FaceSharedPtr FacePointer = std::shared_ptr<NekFace>(new NekFace(vids));
            // std::pair<FaceSetPointer::iterator, bool> testInsPointer;
            // testInsPointer = m_TestFaceSet.insert(FacePointer);

            faceNodes2[icomm][3*i]     = v0l+1; // tria_vert[3*(pos-1)];
            faceNodes2[icomm][3*i+1]   = v1l+1; // tria_vert[3*(pos-1)+1];
            faceNodes2[icomm][3*i+2]   = v2l+1; // tria_vert[3*(pos-1)+2];
          }
        }


        // Check matching of input interface triangles with the set ones
        if( !PMMG_Check_Set_FaceCommunicators(parmesh,ncomm,ntifc,
                                              color_face,faceNodes2) ) {
          printf("### FAILED:: Wrong set face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: Set the correct face communicators\n");
        }
        
        int *ref2       = (int*)calloc(nTriangles,sizeof(int));
        int *required2  = (int*)calloc(nTriangles,sizeof(int));
        int *triaNodes2 = (int*)calloc(3*nTriangles,sizeof(int));

        if ( PMMG_Get_triangles(parmesh,triaNodes2,ref2,required2) != 1 ) {
          fprintf(stderr,"FAILED:: Unable to get mesh triangles\n");
          ier = PMMG_STRONGFAILURE;
        }
        else
        {
            printf("### SUCCES:: retrieved all mesh triangles\n");
        }

        int** faceNodes_out = (int **) malloc(next_face_comm*sizeof(int *));
        for( int icomm = 0; icomm < next_face_comm; icomm++ )
        {
              faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
              for(int i = 0; i < nitem_face_comm[icomm]; i++ )
              {
                  int pos = out_tria_loc[icomm][i];
                  faceNodes_out[icomm][3*i]   = triaNodes2[3*(pos-1)];
                  faceNodes_out[icomm][3*i+1] = triaNodes2[3*(pos-1)+1];
                  faceNodes_out[icomm][3*i+2] = triaNodes2[3*(pos-1)+2];
              }
        }

        
        // Check matching of input interface triangles with the output ones
        if( !PMMG_Check_Get_FaceCommunicators(parmesh,ncomm,ntifc,
                                              color_face,faceNodes2,
                                              next_face_comm,nitem_face_comm,
                                              color_face_out,faceNodes_out) )
        {
          printf("### FAILED:: Wrong retrieved face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: retrieved the correct face communicators\n");
        }
        
        int lvt = 0;
        int lvt_o = 0;
        std::map<int,int> sharedVrts_Owned;
        std::map<int,int> sharedVert;
        for( int icomm = 0; icomm < next_face_comm; icomm++ )
        {
            faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
            int nPartFace                 = nPartFace + nitem_face_comm[icomm];
            //rank2icomm[color_face_out[icomm]] = icomm;

            if(world_rank < color_face_out[icomm])
            {
                int nTshared_owned = nTshared_owned + nitem_face_comm[icomm];

                for(int i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft   = out_tria_loc[icomm][i];
                    int reff = ref2[ft-1];

                    for(int k=0;k<3;k++)
                    {
                        int vt = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];

                        if(sharedVrts_Owned.find(vt)==sharedVrts_Owned.end())
                        {
                            sharedVrts_Owned[vt] = lvt_o;
                            sharedVert[vt]       = lvt;

                            lvt++;
                            lvt_o++;
                        }
                    }
                }
            }
            else
            {
                for(int i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft = out_tria_loc[icomm][i];

                    for(int k=0;k<3;k++)
                    {
                        int vt3 = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];
                        if(sharedVert.find(vt3)==sharedVert.end())
                        {
                            sharedVert[vt3] = lvt;
                            lvt++;
                        }
                    }
                }
            }
        }
    }
}









PMMG_pParMesh InitializeParMMGmeshFromHyperSolveInputs(MPI_Comm comm, 
                                                        PartObjectLite* partition,
                                                        FaceSetPointer ActualSharedFaces,
                                                        FaceSetPointer InterFaceFaces,
                                                        FaceSetPointer OwnedBoundaryFaces,
                                                        int nInteriorSharedFaces,
                                                        std::vector<std::vector<double> > t_metric)
{
    int ier;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    PMMG_pParMesh   parmesh;
    int k = 0;

    std::map<int,std::vector<int> > gE2gV_t                     = partition->getElem2VertMap();
    std::map<int,int> locv2tagvID                               = partition->getLocalVert2GlobalVert();
    std::map<int,int> tagv2locvID                               = partition->getGlobalVert2LocalVert();
    std::set<int> Owned_Elem_t                                  = partition->getLocalElemSet();
    std::map<int,std::vector<double> > LocalVertsMap_t          = partition->getLocalVertsMap();
    FaceSetPointer ExternalFacesForRank                         = partition->getExternalFacesForRankFaceMap();
    FaceSetPointer OwnedBoundaryFaceFaceMap                     = partition->getOwnedBoundaryFaceFaceMap();
    FaceSetPointer ExternalInterFaceFaceMap                     = partition->getExternalInterFaceFaceMap();
    FaceSetPointer SharedFacesForRank                           = partition->getSharedFacesForRankMap();
    std::map<int,std::vector<int> > sharedFace2Verts            = partition->getSharedFaceMap();

    //std::cout << "numbers check " << world_rank << " " << ExternalFacesForRank.size() << " " << ExternalInterFaceFaceMap.size() << std::endl;

    //int nTriangles  = OverallFaceMapOnRank.size()+OwnedBoundaryFaces.size();//ExternalFacesForRank.size()+OwnedBoundaryFaceFaceMap.size();
    int nTriangles  = SharedFacesForRank.size()+OwnedBoundaryFaceFaceMap.size()+ExternalInterFaceFaceMap.size();//+ExternalInterFaceFaceMap.size();
    int nTetrahedra = Owned_Elem_t.size();
    int nVertices   = locv2tagvID.size();
    int nEdges      = 0;
    
    //std::cout << "SharedFacesForRank " << SharedFacesForRank.size() << " " << OwnedBoundaryFaceFaceMap.size() << std::endl;


    //std::cout << "SharedFacesForRank " << SharedFacesForRank.size() << " " << OwnedBoundaryFaceFaceMap.size() << " " << world_rank << std::endl;
    
    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);

    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,0,nTriangles,0,nEdges) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    //std::cout << nTriangles << " = nTriangles   " << nVertices << " = nVertices " << nTetrahedra << " = nTetrahedra" << world_rank << std::endl;

    //std::cout << "LocalVertsMap_t " << LocalVertsMap_t.size() << " " <<   tagv2locvID.size() << " " << locv2tagvID.size() << std::endl;
    // std::cout << "nVertices " << nVertices << " " << world_rank << std::endl;
    //PMMG_Set_metSize(PMMG_pParMesh parmesh,int typEntity,int np,int typSol)
    if ( PMMG_Set_metSize(parmesh,MMG5_Vertex,nVertices,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);
    std::map<int,std::vector<double> >::iterator vrtit;
    
    //for (k=0; k<nVertices; ++k )
    for(vrtit=LocalVertsMap_t.begin();vrtit!=LocalVertsMap_t.end();vrtit++)
    {
        int tagvid  = vrtit->first;
        //int glovid  = tag2globalV[tagvid];
        //int locvid  = globalv2localvID[glovid];
        //std::cout << "locvid = " << locvid << " --- " << glovid << std::endl;
        double vx  = LocalVertsMap_t[tagvid][0];
        double vy  = LocalVertsMap_t[tagvid][1];
        double vz  = LocalVertsMap_t[tagvid][2];
        int locvid = tagv2locvID[tagvid];
        // //std::cout <<"coords = "<< vx << " " << vy << " " << vz << std::endl;
        
        int vref = 86;

        if ( PMMG_Set_vertex(parmesh,vx,vy,vz, vref, locvid+1) != 1 )
        {
        MPI_Finalize();
        exit(EXIT_FAILURE);
        }

        std::vector<double> tensor(6,0);

        // tensor[0] = t_metric[tagvid][0];
        // tensor[1] = t_metric[tagvid][1];
        // tensor[2] = t_metric[tagvid][2];
        // tensor[3] = t_metric[tagvid][3];
        // tensor[4] = t_metric[tagvid][4];
        // tensor[5] = t_metric[tagvid][5];

        tensor[0] = 1.0;
        tensor[1] = 0.0;
        tensor[2] = 0.0;
        tensor[3] = 1.0;
        tensor[4] = 0.0;
        tensor[5] = 1.0;

        if(PMMG_Set_tensorMet(parmesh,tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],locvid+1)!=1)
        {
         MPI_Finalize();
         exit(EXIT_FAILURE);
        }
    }
    
    FaceSetPointer::iterator ftit;
    
    int faceID          = 0;
    int bcf             = 0;
    int SharedF         = 0;
    int InternalF       = 0;
    int InternalFreal   = 0;
    int here            = 0;
    int foundaswell     = 0;
    int dd              = 0;
    int ee              = 0;

    int nVrts = locv2tagvID.size();
    k = 1;

    for(ftit=ExternalFacesForRank.begin();ftit!=ExternalFacesForRank.end();ftit++)
    {
        FaceSetPointer::iterator ConsideredSharedFace = ActualSharedFaces.find((*ftit));
        // Mechanism to match up the shell faces from the adapted tetrahedra to the fixed prisms.
        if(ConsideredSharedFace!=ActualSharedFaces.end())
        {
            int ref 		          = (*ConsideredSharedFace)->GetFaceRef();
            std::vector<int> faceVrts = (*ConsideredSharedFace)->GetEdgeIDs();

            // std::cout << "REF " << ref << std::endl;
        
            int tagvid0   = faceVrts[0];
            int tagvid1   = faceVrts[1];
            int tagvid2   = faceVrts[2];

            int locvid0   = tagv2locvID[tagvid0];
            int locvid1   = tagv2locvID[tagvid1];
            int locvid2   = tagv2locvID[tagvid2];

            if ( PMMG_Set_triangle(parmesh,locvid0+1,locvid1+1,locvid2+1,ref,k) != 1 )
            {
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
            k = k + 1;
        }

        FaceSetPointer::iterator ConsideredIntFace = ExternalInterFaceFaceMap.find((*ftit));

        if(ConsideredIntFace!=ExternalInterFaceFaceMap.end())
        {
            int ref 		          = 10;//(*InterFaceFPointer)->GetFaceRef();
            std::vector<int> faceVrts = (*ConsideredIntFace)->GetEdgeIDs();

            int tagvid0   = faceVrts[0];
            int tagvid1   = faceVrts[1];
            int tagvid2   = faceVrts[2];

            int locvid0   = tagv2locvID[tagvid0];
            int locvid1   = tagv2locvID[tagvid1];
            int locvid2   = tagv2locvID[tagvid2];

            // if ( PMMG_Set_triangle(parmesh,locvid0+1,locvid1+1,locvid2+1,ref,k) != 1 )
            // {
            //     MPI_Finalize();
            //     exit(EXIT_FAILURE);
            // }

            // PMMG_Set_requiredTriangle( parmesh, k );
            // k = k + 1;  
            //InternalFreal++;  
        }
    }

    // for(ftit=SharedFacesForRank.begin();ftit!=SharedFacesForRank.end();ftit++)
    // {
    //     int ref 		          = (*ftit)->GetFaceRef();
    //     std::vector<int> faceVrts = (*ftit)->GetEdgeIDs();

    //     // std::cout << "REF " << ref << std::endl;
    
    //     int tagvid0   = faceVrts[0];
    //     int tagvid1   = faceVrts[1];
    //     int tagvid2   = faceVrts[2];

    //     int locvid0   = tagv2locvID[tagvid0];
    //     int locvid1   = tagv2locvID[tagvid1];
    //     int locvid2   = tagv2locvID[tagvid2];

    //     if ( PMMG_Set_triangle(parmesh,locvid0+1,locvid1+1,locvid2+1,ref,k) != 1 )
    //     {
    //         MPI_Finalize();
    //         exit(EXIT_FAILURE);
    //     }
    //     k = k + 1;
    // }


    for(ftit=ExternalInterFaceFaceMap.begin();ftit!=ExternalInterFaceFaceMap.end();ftit++)
    {
        int ref 		          = 10;//(*InterFaceFPointer)->GetFaceRef();
        std::vector<int> faceVrts = (*ftit)->GetEdgeIDs();

        int tagvid0   = faceVrts[0];
        int tagvid1   = faceVrts[1];
        int tagvid2   = faceVrts[2];

        int locvid0   = tagv2locvID[tagvid0];
        int locvid1   = tagv2locvID[tagvid1];
        int locvid2   = tagv2locvID[tagvid2];
        if ( PMMG_Set_triangle(parmesh,locvid0+1,locvid1+1,locvid2+1,ref,k) != 1 )
        {

            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        PMMG_Set_requiredTriangle( parmesh, k );
        k = k + 1;  
        InternalFreal++;  
    }
    //std::cout << "InternalFreal " << InternalFreal << std::endl;
    //std::cout << "ExternalFacesForRank " << ExternalFacesForRank.size() << " " << k<< " " << SharedFacesForRank.size() << " " <<ExternalInterFaceFaceMap.size() << std::endl;  
    //std::cout << "InterInterFaceFPointer " << foundaswell << " " << InternalFreal << " " << dd << " " << ee << std::endl;
    //std::cout << " here " << here << " " << nTriangles  << " " << world_rank << " InternalFreal " << InternalFreal << " " << ExternalFacesForRank.size() << " " << ExternalInterFaceFaceMap.size() << " " << bcf << std::endl;
        
    // if(ExternalFacesForRank.size()==k)
    // {
    //     std::cout << "all Shared interfaces are set correctly" << std::endl;
    // }

    int ibc2 = 0;
    int ibc3 = 0;
    int ibcother = 0;
    for(ftit=OwnedBoundaryFaces.begin();ftit!=OwnedBoundaryFaces.end();ftit++)
    {
        int ref 		          = (*ftit)->GetFaceRef();
        std::vector<int> faceVrts = (*ftit)->GetEdgeIDs();
    
        int tagvid0   = faceVrts[0];
        int tagvid1   = faceVrts[1];
        int tagvid2   = faceVrts[2];

        int locvid0   = tagv2locvID[tagvid0];
        int locvid1   = tagv2locvID[tagvid1];
        int locvid2   = tagv2locvID[tagvid2];


        if ( PMMG_Set_triangle(parmesh,locvid0+1,locvid1+1,locvid2+1,ref,k) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }

        k = k + 1;

        if(ref == 2)
        {
            ibc2++;
        }
        if(ref == 3)
        {
            ibc3++;
        }
    }

    //std::cout << "ibcs = " << ibc2 << " " << ibc3 << " " << ibcother << " " << world_rank << std::endl;
    //std::cout << "nTriangles " << world_rank << " " << OwnedBoundaryFaces.size() << " " << k << " " << OverallFaceMapOnRank.size() << " " << nTriangles<< std::endl;


    int** m_ifc_tria_glo                                  = partition->getParMMGCommFace2GlobalVertMap();
    int** m_ifc_tria_loc                                  = partition->getParMMGCommFace2LocalVertMap();
    int* m_color_face                                     = partition->getParMMGCommColorFace();
    int *m_ntifc                                          = partition->getParMMGCommNFacesPerColor();
    int m_ncomm                                           = partition->getParMMGNComm();

        
    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, m_ncomm);

   
    for(int icomm=0; icomm<m_ncomm; icomm++ )
    {
      //std::cout << "Initialization " << icomm << " " << ncomm << " " << ntifc[icomm] << " " << world_rank << std::endl;
      // Set nb. of entities on interface and rank of the outward proc
      ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                             m_color_face[icomm],
                                             m_ntifc[icomm]);

        //Set local and global index for each entity on the interface
      ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                               m_ifc_tria_loc[icomm],
                                               m_ifc_tria_glo[icomm], 1);
    }

    
    
    

    std::map<int,std::vector<int> >::iterator ittet;
    k = 0;
    int vloc0,vloc1,vloc2,vloc3;
    int vglo0,vglo1,vglo2,vglo3;
    int vtag0,vtag1,vtag2,vtag3;
    // for ( int t = 0;t < nTetrahedra; t++  )
    std::set<int>::iterator its;
    int t = 0;
    
    for(its=Owned_Elem_t.begin();its!=Owned_Elem_t.end();its++)
    {
        //From local element ID get the tag ID (initial global ID).
        int gEid = *its; //Owned_Elem_t[t];
        //Using the tag Element ID, get the tag Vertex ID.
        vtag0 = gE2gV_t[gEid][0];
        vtag1 = gE2gV_t[gEid][1];
        vtag2 = gE2gV_t[gEid][2];
        vtag3 = gE2gV_t[gEid][3];

        //From the tag vertex ID, get the local vertex ID.

        int vloc0 = tagv2locvID[vtag0];
        int vloc1 = tagv2locvID[vtag1];
        int vloc2 = tagv2locvID[vtag2];
        int vloc3 = tagv2locvID[vtag3];

        std::vector<double> P(4*3);
        P[0*3+0]=LocalVertsMap_t[vtag0][0];   P[0*3+1]=LocalVertsMap_t[vtag0][1];    P[0*3+2]=LocalVertsMap_t[vtag0][2];
        P[1*3+0]=LocalVertsMap_t[vtag1][0];   P[1*3+1]=LocalVertsMap_t[vtag1][1];    P[1*3+2]=LocalVertsMap_t[vtag1][2];
        P[2*3+0]=LocalVertsMap_t[vtag2][0];   P[2*3+1]=LocalVertsMap_t[vtag2][1];    P[2*3+2]=LocalVertsMap_t[vtag2][2];
        P[3*3+0]=LocalVertsMap_t[vtag3][0];   P[3*3+1]=LocalVertsMap_t[vtag3][1];    P[3*3+2]=LocalVertsMap_t[vtag3][2];

        double Vtet = GetQualityTetrahedra(P);

        //std::cout << "Vtet " << Vtet << "  " << vloc1+1 << " " << vloc2+1 << " " << vloc3+1 << std::endl;
        
        if(Vtet<0.0)
        {
            std::cout << " negative volume in Element " << t << " on rank " << world_rank  <<std::endl;
        }
        if ( PMMG_Set_tetrahedron(parmesh,vloc0+1,vloc1+1,vloc2+1,vloc3+1,1.0,t+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        t++;
    }
    

    return parmesh;
    
}





void RunParMMGAndTestPartitioningFromHyperSolveInputs(MPI_Comm comm,
                                  PMMG_pParMesh parmesh, 
                                  PartObjectLite* partition,
                                  FaceSetPointer ActualSharedFaces,
                                  Inputs* inputs)
{

    FILE *inm;

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    int ier, opt;

    


    std::set<int> Owned_Elem_t                          = partition->getLocalElemSet();
    std::map<int,int> locv2tagvID                       = partition->getLocalVert2GlobalVert();
    std::map<int,int> tagv2locvID                       = partition->getGlobalVert2LocalVert();
    std::vector<int> face4parmmg                        = partition->getFace4ParMMG(); // checked
    int** ifc_tria_glob                                 = partition->getParMMGCommFace2GlobalVertMap();
    int** ifc_tria_loc                                  = partition->getParMMGCommFace2LocalVertMap();
    int* color_face                                     = partition->getParMMGCommColorFace();
    int *ntifc                                          = partition->getParMMGCommNFacesPerColor();
    int ncomm                                           = partition->getParMMGNComm();

    std::map<int,int> locShF2globShF                    = partition->getLocalSharedFace2GlobalSharedFace();
    std::map<int,std::vector<int> > sharedface2vert     = partition->getSharedFaceMap();
    //std::map<int,std::vector<int> > face2node         = partition->getFace2VertNewMap();
    //FaceSetPointer OverallFaceMapOnRank                         = partition->getOverallFaceMapOnRank();
    FaceSetPointer ExternalFacesForRank                         = partition->getExternalFacesForRankFaceMap();
    FaceSetPointer OwnedBoundaryFaceFaceMap                     = partition->getOwnedBoundaryFaceFaceMap();
    FaceSetPointer ExternalInterFaceFaceMap                     = partition->getExternalInterFaceFaceMap();
    std::map<int,std::vector<int> > sharedFace2Verts            = partition->getSharedFaceMap();
    FaceSetPointer SharedFacesForRank                           = partition->getSharedFacesForRankMap();

    int nTriangles  = SharedFacesForRank.size()+OwnedBoundaryFaceFaceMap.size()+ExternalInterFaceFaceMap.size();//+ExternalInterFaceFaceMap.size();
    // int nTriangles = OverallFaceMapOnRank.size();
    int nVertices   = locv2tagvID.size();
    int nTetrahedra = Owned_Elem_t.size();
    int nEdges      = 0;
    
 
    int nVerticesIN   = 0;
    int nTetrahedraIN = 0;
    int nTrianglesIN  = 0;
    int nEdgesIN      = 0;
    
    if ( PMMG_Get_meshSize(parmesh,&nVerticesIN,&nTetrahedraIN,NULL,&nTrianglesIN,NULL,
                           &nEdgesIN) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }


    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, 0 ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hausd, inputs->hausd) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hgrad, inputs->hgrad) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if( !PMMG_Set_dparameter( parmesh,  PMMG_DPARAM_hgradreq , -1.0 ) ){
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // std::cout << "nVertices input " << nVerticesIN << std::endl;
    // std::cout << "nTetrahedra input " << nTetrahedraIN << std::endl;
    // std::cout << "nTriangles input " << nTrianglesIN << std::endl;
    // std::cout << "nEdges input " << nEdgesIN << std::endl;
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_globalNum, 1 ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };

    
    int ierlib = PMMG_parmmglib_distributed( parmesh );
    
    
    if(ierlib==0 && world_rank == 0)
    {
        std::cout << "SUCCESFULLY adapted the mesh in parallel!" << std::endl;
    }
    else if(ierlib==1 && world_rank == 0)
    {
        std::cout << "FAILED to adapt the mesh in parallel! "<< std::endl;
    }

    int nVerticesOUT   = 0;
    int nTetrahedraOUT = 0;
    int nTrianglesOUT  = 0;
    int nEdgesOUT      = 0;


    if ( PMMG_Get_meshSize(parmesh,&nVerticesOUT,&nTetrahedraOUT,NULL,&nTrianglesOUT,NULL,
                           &nEdgesOUT) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }

    /*
    int             nodeGloNumber,nodeOwner;
    std::map<int,int> loc2globVid;
    std::map<int,int> glob2locVid;
    std::vector<int> globIDs;
    for(int k = 1; k <= nVerticesOUT; k++ )
    {
        if( !PMMG_Get_vertexGloNum( parmesh, &nodeGloNumber, &nodeOwner ) )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    
        loc2globVid[k]=nodeGloNumber;
        glob2locVid[nodeGloNumber]=k;
        globIDs.push_back(nodeGloNumber);
    }
    
    int *required = (int*)calloc(MAX4(nVerticesOUT,
                                      nTetrahedraOUT,
                                      nTrianglesOUT,
                                      nEdgesOUT),sizeof(int));
    
    int *refOUT = (int*)calloc(MAX4(nVerticesOUT,
                                    nTetrahedraOUT,
                                    nTrianglesOUT,
                                    nEdgesOUT),sizeof(int));
    
    int *corner = (int*)calloc(nVerticesOUT,sizeof(int));
    int pos;
    //int *refOUT = new int[nVerticesOUT];

    std::vector<std::vector<int> > outT;
    
    double *vertOUT = (double*)calloc((nVerticesOUT)*3,sizeof(double));
    //std::map<int,std::vector<double> > ref2coordinatesOUT;
    for ( int k=0; k<nVerticesOUT; k++ )
    {
          pos = 3*k;
          if ( PMMG_Get_vertex(parmesh,&(vertOUT[pos]),&(vertOUT[pos+1]),&(vertOUT[pos+2]),
                               &(refOUT[k]),&(corner[k]),&(required[k])) != 1 ) {
            fprintf(inm,"Unable to get mesh vertex %d \n",k);
            ier = PMMG_STRONGFAILURE;
          }
    }

    // RUN TEST TO SEE WHETHER ADAPTATION WAS SUCCESFULL



    if(inputs->niter==0)
    {
        std::cout << "Check the input and outputted shared faces." << std::endl;
        
        int **out_tria_loc;
        int *nitem_face_comm;
        int next_face_comm;
        ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
        int *color_node_out,*color_face_out;
        color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
        nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                 &color_face_out[icomm],
                                                 &nitem_face_comm[icomm]);
        
        out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
        
        // Check matching of input interface nodes with the set ones

        // FaceSharedPtr TestFacePointer = std::shared_ptr<NekFace>(new NekFace(vidstest));
        // FaceSetPointer m_TestFaceSet;

        // Get input triangle nodes
        int** faceNodes2 = (int **) malloc(ncomm*sizeof(int *));
        for( int icomm = 0; icomm < ncomm; icomm++ ) {
          faceNodes2[icomm] = (int *) malloc(3*ntifc[icomm]*sizeof(int));
          for(int i = 0; i < ntifc[icomm]; i++ ) {
              
            int faceID = ifc_tria_loc[icomm][i]-1;
            int faceIDg = ifc_tria_glob[icomm][i]-1;
            int faceID2 = locShF2globShF[faceID];

            int v0 = sharedface2vert[faceID2][0];
            int v1 = sharedface2vert[faceID2][1];
            int v2 = sharedface2vert[faceID2][2];

            int v0l   = tagv2locvID[v0];
            int v1l   = tagv2locvID[v1];
            int v2l   = tagv2locvID[v2];
    
            // int v0l = globalv2localvID[v0];
            // int v1l = globalv2localvID[v1];
            // int v2l = globalv2localvID[v2];

            // FaceSharedPtr FacePointer = std::shared_ptr<NekFace>(new NekFace(vids));
            // std::pair<FaceSetPointer::iterator, bool> testInsPointer;
            // testInsPointer = m_TestFaceSet.insert(FacePointer);

            faceNodes2[icomm][3*i]     = v0l+1; // tria_vert[3*(pos-1)];
            faceNodes2[icomm][3*i+1]   = v1l+1; // tria_vert[3*(pos-1)+1];
            faceNodes2[icomm][3*i+2]   = v2l+1; // tria_vert[3*(pos-1)+2];
          }
        }


        // Check matching of input interface triangles with the set ones
        if( !PMMG_Check_Set_FaceCommunicators(parmesh,ncomm,ntifc,
                                              color_face,faceNodes2) ) {
          printf("### FAILED:: Wrong set face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: Set the correct face communicators\n");
        }
        
        int *ref2       = (int*)calloc(nTrianglesIN,sizeof(int));
        int *required2  = (int*)calloc(nTrianglesIN,sizeof(int));
        int *triaNodes2 = (int*)calloc(3*nTrianglesIN,sizeof(int));
        
        if ( PMMG_Get_triangles(parmesh,triaNodes2,ref2,required2) != 1 ) {
          fprintf(stderr,"FAILED:: Unable to get mesh triangles\n");
          ier = PMMG_STRONGFAILURE;
        }
        else
        {
            printf("### SUCCES:: retrieved all mesh triangles\n");
        }
        
        int** faceNodes_out = (int **) malloc(next_face_comm*sizeof(int *));
        for( int icomm = 0; icomm < next_face_comm; icomm++ )
        {
              faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
              for(int i = 0; i < nitem_face_comm[icomm]; i++ )
              {
                  int pos = out_tria_loc[icomm][i];
                  faceNodes_out[icomm][3*i]   = triaNodes2[3*(pos-1)];
                  faceNodes_out[icomm][3*i+1] = triaNodes2[3*(pos-1)+1];
                  faceNodes_out[icomm][3*i+2] = triaNodes2[3*(pos-1)+2];
              }
        }

        
        // Check matching of input interface triangles with the output ones
        if( !PMMG_Check_Get_FaceCommunicators(parmesh,ncomm,ntifc,
                                              color_face,faceNodes2,
                                              next_face_comm,nitem_face_comm,
                                              color_face_out,faceNodes_out) )
        {
          printf("### FAILED:: Wrong retrieved face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: retrieved the correct face communicators\n");
        }
        
        int lvt = 0;
        int lvt_o = 0;
        std::map<int,int> sharedVrts_Owned;
        std::map<int,int> sharedVert;
        for( int icomm = 0; icomm < next_face_comm; icomm++ )
        {
            faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
            int nPartFace                 = nPartFace + nitem_face_comm[icomm];
            //rank2icomm[color_face_out[icomm]] = icomm;

            if(world_rank < color_face_out[icomm])
            {
                int nTshared_owned = nTshared_owned + nitem_face_comm[icomm];

                for(int i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft   = out_tria_loc[icomm][i];
                    int reff = ref2[ft-1];

                    for(int k=0;k<3;k++)
                    {
                        int vt = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];

                        if(sharedVrts_Owned.find(vt)==sharedVrts_Owned.end())
                        {
                            sharedVrts_Owned[vt] = lvt_o;
                            sharedVert[vt]       = lvt;

                            lvt++;
                            lvt_o++;
                        }
                    }
                }
            }
            else
            {
                for(int i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft = out_tria_loc[icomm][i];

                    for(int k=0;k<3;k++)
                    {
                        int vt3 = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];
                        if(sharedVert.find(vt3)==sharedVert.end())
                        {
                            sharedVert[vt3] = lvt;
                            lvt++;
                        }
                    }
                }
            }
        }
    }
    */
}



