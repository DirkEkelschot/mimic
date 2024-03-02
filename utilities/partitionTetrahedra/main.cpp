
#include <chrono>
#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
#include "../../src/adapt_prismtetratrace.h"
#include "../../src/adapt_repartition.h"
#include "../../src/adapt_output_vtk.h"
#include "../../src/adapt_meshtopology_lite.h"
#include "../../src/adapt_gradreconstruct_lite.h"
#include "../../src/adapt_runparmmg.h"
#include "../../src/adapt_inputs.h"
#include "../../src/adapt_writeus3ddata.h"


#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;






int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j,k;
    clock_t t0_met = clock();

    int ier,opt;
    int debug = 1;
    
    const char* fn_grid="inputs/grid.h5";
    const char* fn_conn="inputs/conn.h5";
    const char* fn_data="inputs/data.h5";
    // Read in the inputs from metric.xml
    Inputs* inputs = ReadXmlFile(comm, "inputs/metric.xml");

    //===========================================================================
    // Read in the data from grid.h5/conn.h5/data.h5 in parallel using HDF5.
    // the outputted data in meshRead contains uniformly distributed data structures that will have to be partitioned.
    mesh* meshRead = ReadUS3DMeshData(fn_conn,fn_grid,fn_data,
                                      inputs->ReadFromStats,
                                      inputs->StateVar,
                                      comm,info);

    int Nel_loc = meshRead->ien.size();
    //Start building the trace object that contains the information regarding the prism-tetra interfaces.
    // It contains the unique vertex and face information.
    //====================================================================================
    

    clock_t start1, end1, start, end;
    double dur_max1,time_taken1;

    start1 = clock();

    PrismTetraTrace* pttrace = new PrismTetraTrace(comm, 
                                                   meshRead->element2rank, 
                                                   meshRead->ife,
                                                   meshRead->ifn,
                                                   meshRead->iet, 
                                                   meshRead->nElem, 
                                                   meshRead->nFace, 
                                                   meshRead->nVert);

    
    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        cout << "Time taken to broadcast boundary layer/tetra trace data is : " << fixed
        << dur_max1 << setprecision(16);
        cout << " sec " << endl;
    }


    std::map<int,std::map<int,int> > trace_elem     = pttrace->GetTrace();
    std::map<int,std::vector<int> > trace_verts     = pttrace->GetTraceVerts();
    std::map<int,int> unique_trace_verts2refmap     = pttrace->GetUniqueTraceVerts2RefMap();
    std::map<int,std::vector<int> > leftright_trace = pttrace->GetLeftRightElements();
    std::map<int,int> trace_ref                     = pttrace->GetTraceRef();
    FaceSetPointer FaceTraceRefs                    = pttrace->GetRefTraceFaceSet();

    //Filter out the tetrahedra and prisms into seperate maps from the IO data structures (meshRead).
    //=====================================================================================
    std::map<int,std::vector<int> > tetras_e2v,tetras_e2f,tetras_e2e;
    std::map<int,std::vector<double> > tetras_data;
    std::map<int,std::vector<int> > prisms_e2v,prisms_e2f,prisms_e2e;
    std::map<int,std::vector<double> > prisms_data;

    int ntetra      = meshRead->ntetra;
    int nprism      = meshRead->nprism;
    std::map<int,std::vector<int> >::iterator itmiv;
    int foundte     = 0;
    int foundpr     = 0;

    int ndata    = meshRead->interior.begin()->second.size();
    int nTrcFace = trace_verts.size();

    std::map<int,double> tracePrismData;
    std::map<int,double> traceTetraData;
    std::map<int,int> tetra2type;
    std::map<int,int> prism2type;   

    if(world_rank == 0)
    {
        std::cout << "Start filtering the element types..." << std::endl; 
    }

    for(itmiv=meshRead->ien.begin();itmiv!=meshRead->ien.end();itmiv++)
    {
        int elid   = itmiv->first;
        int eltype = meshRead->iet[elid];

        if(eltype == 2)
        {
            tetras_e2v[elid]  = itmiv->second;
            tetras_e2f[elid]  = meshRead->ief[elid];
            tetras_e2e[elid]  = meshRead->iee[elid];
            tetras_data[elid] = meshRead->interior[elid]; 
            tetra2type[elid]  = eltype;
        }
        else
        {
            prisms_e2v[elid]  = itmiv->second;
            prisms_e2f[elid]  = meshRead->ief[elid];
            prisms_e2e[elid]  = meshRead->iee[elid];
            prisms_data[elid] = meshRead->interior[elid];
            prism2type[elid] = eltype;
            int nf = 5;
            for(int j=0;j<nf;j++)
            {
                int faceID = meshRead->ief[elid][j];
                if(trace_verts.find(faceID)!=trace_verts.end())
                {
                    if(tracePrismData.find(elid)==tracePrismData.end())
                    {
                        tracePrismData[elid] = meshRead->interior[elid][1];
                    }
                }
            }
        }   
    }

    if(world_rank == 0)
    {
        std::cout << "Done filtering the element types..." << std::endl; 
    }
    
    //I am adding the prism elements and their data to the ghost map so that that data is in the boundaries data structures.
    std::map<int,double> tracePrismData_glob = AllGatherMap_T(tracePrismData,comm);
    
    std::map<int,double>::iterator itr;
    for(itr=tracePrismData_glob.begin();itr!=tracePrismData_glob.end();itr++)
    {
        int elid    = itr->first;
        double data = itr->second;

        std::vector<double> rowghost(2,0.0);
        rowghost[0] = 0.0;
        rowghost[1] = data;

        if(meshRead->ghost.find(elid)==meshRead->ghost.end())
        {
            meshRead->ghost[elid] = rowghost;
        }
    }
    int nLocTetra  = tetras_e2v.size();
    int nLocPrism  = prisms_e2v.size();
    int nElemsGlob_T = 0;
    int nElemsGlob_P = 0;
    MPI_Allreduce(&nLocTetra, &nElemsGlob_T, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nLocPrism, &nElemsGlob_P, 1, MPI_INT, MPI_SUM, comm);

    

    if(world_rank == 0)
    {
        std::cout << "Done trace operation..." << std::endl; 
    }

    //=========END FILTERING OUT TETRA AND PRISMS FROM IO DATA STRUCTURES===============

    // we need to pass the number of verts per element in case the partition has no elements of this type.

 
    //  You can call it like this : start = time(NULL); 
    // in both the way start contain total time in seconds 
    // since the Epoch. 


    double dur_max,time_taken;

    start = clock();
    RepartitionObject* tetra_repart = new RepartitionObject(meshRead, 
                                                        tetras_e2v, 
                                                        tetras_e2f,
                                                        tetras_e2e,
                                                        tetra2type,
                                                        pttrace, 
                                                        tetras_data,
                                                        2,
                                                        comm);

    tetras_e2v.clear();
    tetras_e2f.clear();
    tetras_e2e.clear();
    tetra2type.clear();
    tetras_data.clear();
                                                        
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;

    
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        cout << "Time taken to execute repartioning tetrahedera is : " << fixed 
        << dur_max << setprecision(16); 
        cout << " sec " << endl;
    }

    
    
    tetras_e2v.clear();tetras_e2f.clear();tetras_e2e.clear();

    tetra_repart->buildUpdatedVertexAndFaceNumbering(comm, pttrace, meshRead->ranges_id);
    tetra_repart->buildInteriorSharedAndBoundaryFaceMaps(comm, pttrace, meshRead->ranges_id);

    // std::vector<int> face4parmmg                        = tetra_repart->getFace4ParMMG(); // checked
    // std::map<int,int> global2tagF                       = tetra_repart->getGlobal2TagFMap();
    // int** ifc_tria_glob                                 = tetra_repart->getParMMGCommFace2GlobalVertMap();
    // int** ifc_tria_loc                                  = tetra_repart->getParMMGCommFace2LocalVertMap();
    // int* color_face                                     = tetra_repart->getParMMGCommColorFace();
    // int *ntifc                                          = tetra_repart->getParMMGCommNFacesPerColor();
    // int ncomm                                           = tetra_repart->getParMMGNComm();

    // std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();
    // std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
    // std::map<int, std::vector<double> > LocalVertsMap_t = tetra_repart->getLocalVertsMap();

    //=======================================================================================
    //=============================COMPUTE GRADIENTS/HESSIAN=================================
    //=======================================================================================
    // std::map<int,std::vector<double> > U                = tetra_repart->getElement2DataMap();
    start = clock();
    std::map<int,std::vector<double> > tetra_grad = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
                                                                              pttrace,
                                                                              meshRead->ghost,
                                                                              meshRead->nElem,
                                                                              1, 
                                                                              1,
                                                                              comm);
    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;

    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        cout << "Time taken to execute calculating first and second order gradients : " << fixed 
        << dur_max << setprecision(16); 
        cout << " sec " << endl;
    }

    //==============================================

    std::map<int,std::set<int> >::iterator itmis;

    if (inputs->niter == -1)
    {
        std::map<int,std::set<int> > node2element_map = tetra_repart->GetNode2ElementMap();
        std::map<int,std::vector<std::vector<double> > > metric_vmap;
        // ====================== Generating dummy metric ==========================
        for(itmis=node2element_map.begin();itmis!=node2element_map.end();itmis++)
        {
            int gvid = itmis->first;
            std::vector<std::vector<double> > metric(3);
            for(int i=0;i<3;i++)
            {
                std::vector<double> row(3,0);
                row[i] = 1.0;

                metric[i] = row;
            }

            metric_vmap[gvid] = metric;
        }
        // ====================== Done generating dummy metric =====================
        PMMG_pParMesh parmesh = InitializeParMMGmesh(comm, tetra_repart, pttrace, meshRead->ranges_id, metric_vmap);

        RunParMMGAndTestPartitioning(comm, parmesh, tetra_repart, pttrace,  meshRead->ranges_id, inputs);
        delete tetra_repart;
    }
    else
    {   
        
        std::map<int,std::vector<std::vector<double> > > metric_vmap = ComputeMetric_Lite(comm, 
                                                                            tetra_repart,
                                                                            tetra_grad, 
                                                                            inputs);

        tetra_grad.clear();

        

        PMMG_pParMesh parmesh = InitializeParMMGmesh(comm, tetra_repart, pttrace, meshRead->ranges_id, metric_vmap);
        delete tetra_repart;

        metric_vmap.clear();

        start = clock();
        RepartitionObject* prism_repart = new RepartitionObject(meshRead, 
                                                        prisms_e2v, 
                                                        prisms_e2f,
                                                        prisms_e2e,
                                                        prism2type,
                                                        pttrace, 
                                                        prisms_data,
                                                        0,
                                                        comm);
        prisms_e2v.clear();
        prisms_e2f.clear();
        prisms_e2e.clear();
        prism2type.clear();
        prisms_data.clear();

        end = clock();
        time_taken = ( end - start) / (double) CLOCKS_PER_SEC;
        dur_max;
        MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
            cout << "Time taken to execute repartioning the boundary layer mesh is : " << fixed 
            << dur_max << setprecision(16); 
            cout << " sec " << endl;
        }
        
    
        meshRead->ien.clear();
        meshRead->ife.clear();
        meshRead->ifn.clear();
        meshRead->iet.clear();
        meshRead->iee.clear();
        meshRead->if_ref.clear();
        meshRead->element2rank.clear();
        meshRead->ie_Nv.clear();
        meshRead->ie_Nf.clear();
        meshRead->if_Nv.clear();
        meshRead->interior.clear();
        meshRead->ghost.clear();

        prism_repart->buildUpdatedVertexAndFaceNumbering(comm,pttrace, meshRead->ranges_id);
        prism_repart->buildInteriorSharedAndBoundaryFaceMaps(comm,pttrace, meshRead->ranges_id);

        std::map<int,int> NonSharedVertsOwned_P       = prism_repart->getNonSharedVertsOwned();
        std::map<int,int> SharedVertsOwned_P          = prism_repart->getSharedVertsOwned();
        int nNonSharedVertsOwned_P                    = NonSharedVertsOwned_P.size();
        int nSharedVertsOwned_P                       = SharedVertsOwned_P.size();
        DistributedParallelState* distPrismIntVerts   = new DistributedParallelState(nNonSharedVertsOwned_P,comm);
        DistributedParallelState* distPrismShaVerts   = new DistributedParallelState(nSharedVertsOwned_P,comm);
        int nVertsGlob_P                              = distPrismIntVerts->getNel()+distPrismShaVerts->getNel();
        
        // Next to all the output data structures for us3d like xcn and ifn we also need to build the tracerefV2globalVidInTotalMesh map 
        // since this is defined by the tetrahedra ordering.
        std::map<int,int> tagE2gE_P = prism_repart->getElementTag2GlobalElement();

        std::vector<int> ifn_adaptedTetra;
        std::map<int,std::vector<int> > bcref2bcface_adaptedTetra;
        std::map<int,std::vector<int> > BoundaryFaces_adaptedTetra;
        std::map<int,int> glob2locVid_adaptedTetra;
        std::map<int,int> LocationSharedVert_adaptedTetra;
        std::vector<double> vertices_AdaptedTetra;
        std::vector<std::vector<int> > elements_AdaptedTetra;
        std::map<int,int> lh_bc_adaptedTetra;
        std::map<int,int> new_globE2locE_adaptedTetra;

        std::map<int,int> adaptedTetraVert2TotalAdaptedMeshVertMap;
        std::map<int,int> tracerefV2globalVidInTotalMesh;

        

        if(world_rank == 0)
        {
            std::cout << "running and writing tetra data..." << std::endl;
        }

        int nLocSharedVerts_adaptedTetra = 0;
        int nLocInteriorVerts_adaptedTetra = 0;

        RunParMMGandWriteTetraUS3Dformat(comm, 
                                        parmesh, 
                                        pttrace, 
                                        inputs, 
                                        nElemsGlob_P, 
                                        nVertsGlob_P, 
                                        tracerefV2globalVidInTotalMesh, 
                                        ifn_adaptedTetra,
                                        bcref2bcface_adaptedTetra,
                                        BoundaryFaces_adaptedTetra,
                                        glob2locVid_adaptedTetra,
                                        LocationSharedVert_adaptedTetra,
                                        vertices_AdaptedTetra,
                                        elements_AdaptedTetra,
                                        adaptedTetraVert2TotalAdaptedMeshVertMap,
                                        lh_bc_adaptedTetra,
                                        new_globE2locE_adaptedTetra,
                                        nLocInteriorVerts_adaptedTetra,
                                        nLocSharedVerts_adaptedTetra,
                                        tagE2gE_P);
        
        int nLocTotVrts_T                           = nLocInteriorVerts_adaptedTetra+nLocSharedVerts_adaptedTetra;
        DistributedParallelState* distnLocTotVrts_T = new DistributedParallelState(nLocTotVrts_T,comm);
        int* TotVrts_offsets_T                      = distnLocTotVrts_T->getOffsets();
        int nTotVertsTets                           = distnLocTotVrts_T->getNel();
        int TotVrts_offset_T                        = TotVrts_offsets_T[world_rank];
        std::vector<int> ifn_P;

        if(world_rank == 0)
        {
            std::cout << "writing prism data..." << std::endl;
        }

        WritePrismsUS3DFormat(comm,prism_repart,pttrace,tracerefV2globalVidInTotalMesh,ifn_P,meshRead->ranges_id);
                
        delete pttrace;

        int ftott = (int)ifn_adaptedTetra.size()/8;
        int ftotp = (int)ifn_P.size()/8;

        if(world_rank == 0)
        {
            std::cout << "writing boundary data..." << std::endl;
        }

        std::vector<std::vector<int> > bcArrays;
        std::map<int,int> bcsizing;
        std::vector<int> bciTot_offsets;
        int nTotBCFaces = 0;
        std::set<int> bcsToT;
        std::vector<int> nlbc;
        std::vector<int> bcid;
        std::vector<int> bci_offsets;
        std::vector<int> zdefs_new;
        WriteBoundaryDataUS3DFormat(comm,
                                    prism_repart,
                                    ifn_adaptedTetra,
                                    ifn_P,
                                    bcref2bcface_adaptedTetra,
                                    unique_trace_verts2refmap,
                                    tracerefV2globalVidInTotalMesh,
                                    BoundaryFaces_adaptedTetra,
                                    glob2locVid_adaptedTetra,
                                    elements_AdaptedTetra,
                                    vertices_AdaptedTetra,
                                    LocationSharedVert_adaptedTetra,
                                    adaptedTetraVert2TotalAdaptedMeshVertMap,
                                    lh_bc_adaptedTetra,
                                    new_globE2locE_adaptedTetra,
                                    bcArrays,
                                    bcsizing,
                                    meshRead->zdefs,
                                    meshRead->zone2bcref,
                                    meshRead->zone2name,
                                    meshRead->znames,
                                    bciTot_offsets,
                                    nTotBCFaces,
                                    bcsToT,
                                    nlbc,
                                    bcid,
                                    bci_offsets,
                                    zdefs_new);
        
        
        std::vector<double> interiorverts_prism(NonSharedVertsOwned_P.size()*3,0);
        std::vector<double> sharedverts_prism(SharedVertsOwned_P.size()*3,0);
        int ivi = 0;
        std::map<int, std::vector<double> > LocalVertsMap_P = prism_repart->getLocalVertsMap();
        std::map<int,int>::iterator itmii;
        for(itmii=NonSharedVertsOwned_P.begin();itmii!=NonSharedVertsOwned_P.end();itmii++)
        {
            int tag                      = itmii->first;

            interiorverts_prism[ivi*3+0] = LocalVertsMap_P[tag][0];
            interiorverts_prism[ivi*3+1] = LocalVertsMap_P[tag][1];
            interiorverts_prism[ivi*3+2] = LocalVertsMap_P[tag][2];

            ivi++;            
        }


        int ivs = 0;
        for(itmii=SharedVertsOwned_P.begin();itmii!=SharedVertsOwned_P.end();itmii++)
        {
            int tag                    = itmii->first;

            sharedverts_prism[ivs*3+0] =  LocalVertsMap_P[tag][0];
            sharedverts_prism[ivs*3+1] =  LocalVertsMap_P[tag][1];
            sharedverts_prism[ivs*3+2] =  LocalVertsMap_P[tag][2];

            ivs++;
        }
        
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        
        int nTetraLoc = elements_AdaptedTetra.size();
        int nPrismLoc = prism_repart->getLocElem().size();

        std::vector<int> prisms_type(nPrismLoc,0);
        DistributedParallelState* PrismElementDistr = new DistributedParallelState(nPrismLoc,comm);
        std::map<int,int> element2type_bl = prism_repart->GetElement2TypeOnRankMap();
        std::map<int,int> new2tag = prism_repart->getGlobalElement2ElementTag();

        for(int i=0;i<nPrismLoc;i++)
        {
            int newid = PrismElementDistr->getOffsets()[world_rank] + i + 1;
            int tag_el_id = new2tag[newid];
            prisms_type[i] = element2type_bl[tag_el_id];
        }

        std::vector<int> tetra_type2(nTetraLoc,0);
        for(int i=0;i<nTetraLoc;i++)
        {
            tetra_type2[i] = 2;
        }
        

        DistributedParallelState* distTetra         = new DistributedParallelState(elements_AdaptedTetra.size(),comm);
        int ToTElements_offset_tetra                = distTetra->getOffsets()[world_rank];

        DistributedParallelState* distPrism         = new DistributedParallelState(prism_repart->getLocElem().size(),comm);
        int ToTElements_offset_prism                = distPrism->getOffsets()[world_rank];

        DistributedParallelState* distFaceTetra     = new DistributedParallelState(ftott,comm);
        DistributedParallelState* distFacePrisms    = new DistributedParallelState(ftotp,comm);

        int ToTElements_prism                       = distPrism->getNel();  
        int ToTElements_tetra                       = distTetra->getNel();
        int nTotElements                            = ToTElements_prism + ToTElements_tetra;
        int nTotInteriorFaces_prism                 = distFacePrisms->getNel();
        int nTotInteriorFaces_tetra                 = distFaceTetra->getNel();
        int* TotIntFaces_offsets_prism              = distFacePrisms->getOffsets();
        int* TotIntFaces_offsets_tetra              = distFaceTetra->getOffsets();
        int nvertices                               = (int)vertices_AdaptedTetra.size()/3;

        DistributedParallelState* distTetraVerts    = new DistributedParallelState(nvertices,comm);
        DistributedParallelState* distPrismVerts    = new DistributedParallelState(nNonSharedVertsOwned_P+nSharedVertsOwned_P,comm);
        int nTotFaces                               = nTotInteriorFaces_prism + nTotInteriorFaces_tetra + nTotBCFaces;
        int nTotTetraVerts                          = distTetraVerts->getNel();
        int nTotPrismVerts                          = distPrismVerts->getNel();
        int nTotVertsPrismTetra                     = nTotPrismVerts+nTotTetraVerts;
        int nTotPrismIntVerts                       = distPrismIntVerts->getNel();
        int nTotPrismShaVerts                       = distPrismShaVerts->getNel();
        int TotPrismVerts_offset_int                = distPrismIntVerts->getOffsets()[world_rank];
        int TotPrismVerts_offset_sha                = distPrismShaVerts->getOffsets()[world_rank];



        

        hid_t ret;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);
        hid_t file_id = H5Fcreate("grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

        hid_t    dset_id;
        hid_t    filespace;
        hid_t    memspace;
        hid_t    status;
        hsize_t     dimsf[2];
        hsize_t     countH5[2];
        hsize_t     offsetH5[2];
        
        hsize_t dimsf_att = 1;
        hid_t att_space = H5Screate_simple(1, &dimsf_att, NULL);
        hid_t type =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type, 14);
        ret = H5Tset_strpad(type,H5T_STR_SPACEPAD);
        hid_t attr_id   = H5Acreate (file_id, "filetype", type, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri[] = "US3D Grid File";
        status = H5Awrite(attr_id, type, &stri);
        H5Aclose(attr_id);
        
        hsize_t dimsf_att2 = 1;
         att_space = H5Screate_simple(1, &dimsf_att2, NULL);
        hid_t type2 =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type2, 5);
        ret = H5Tset_strpad(type2,H5T_STR_SPACEPAD);
        attr_id   = H5Acreate (file_id, "filevers", type2, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri2[] = "1.1.8";
        status = H5Awrite(attr_id, type2, &stri2);
        H5Aclose(attr_id);
        
        hid_t group_info_id  = H5Gcreate(file_id, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

        hsize_t dimsf_att3 = 1;
        att_space = H5Screate_simple(1, &dimsf_att3, NULL);
        hid_t type3 =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type3, 10);
        ret = H5Tset_strpad(type3,H5T_STR_SPACEPAD);
        attr_id   = H5Acreate (group_info_id, "date", type3, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri3[] = "27-05-1987";
        status = H5Awrite(attr_id, type3, &stri3);
        H5Aclose(attr_id);
        
        hid_t group_grid_id  = H5Gcreate(group_info_id, "grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        dimsf_att = 1;
        att_space = H5Screate_simple(1, &dimsf_att, NULL);
        attr_id   = H5Acreate (group_grid_id, "nc", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        int value = nTotElements;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "nf", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotFaces;// nTotInteriorFaces_prism+nTotInteriorFaces+nTotBCFaces;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "ng", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotBCFaces;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "nn", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotPrismVerts+nTotTetraVerts;//ToTVrts;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        //std::cout << " NN " << nTotVertsPrismTetra << std::endl;
        //====================================================================================
        // Add iet map to the grid.h5 file
        //====================================================================================

        
        dimsf[0] = nTotElements;
        dimsf[1] = 1;
        filespace = H5Screate_simple(2, dimsf, NULL);

        dset_id = H5Dcreate(file_id, "iet",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = prisms_type.size();
        countH5[1]  = 1;
        
        offsetH5[0] = ToTElements_offset_prism;
        offsetH5[1] = 0;

        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, prisms_type.data());
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        countH5[0]  = tetra_type2.size();
        countH5[1]  = 1;

        offsetH5[0] = ToTElements_prism+ToTElements_offset_tetra;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        //filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, tetra_type2.data());
        //====================================================================================
        // Add xcn map to the grid.h5 file
        //====================================================================================
        
        
        dimsf[0] = nTotPrismIntVerts+nTotPrismShaVerts+nTotTetraVerts;
        dimsf[1] = 3;
        filespace = H5Screate_simple(2, dimsf, NULL);
        
        dset_id = H5Dcreate(file_id, "xcn",
                            H5T_NATIVE_DOUBLE, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = NonSharedVertsOwned_P.size();
        countH5[1]  = 3;
        
        offsetH5[0] = TotPrismVerts_offset_int;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, interiorverts_prism.data());

        
        // delete xcn_prisms_int;
    
//        dimsf[0] = nTotPrismShaVerts_v2;//+nTotTetraVerts_v2;
//        dimsf[1] = xcn_prisms_shared->getNcol();
//        filespace = H5Screate_simple(2, dimsf, NULL);
        
//        dset_id = H5Dcreate(file_id, "xcn",
//                            H5T_NATIVE_DOUBLE, filespace,
//                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        countH5[0]  = SharedVertsOwned_P.size();
        countH5[1]  = 3;
        
        
        offsetH5[0] = nTotPrismIntVerts+TotPrismVerts_offset_sha;
        offsetH5[1] = 0;
        // countH5[0]  = nTotPrismShaVerts;
        // countH5[1]  = 3;
        
        
        // offsetH5[0] = nTotPrismIntVerts+TotPrismVerts_offset_sha;
        // offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, sharedverts_prism.data());
        
        
        // delete xcn_prisms_shared;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
        countH5[0]  = nvertices;
        countH5[1]  = 3;

        offsetH5[0] = nTotPrismIntVerts+nTotPrismShaVerts+TotVrts_offset_T;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);
        //filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &vertices_AdaptedTetra.data()[0]);
        // delete xcn_parmmg;
        
       
        //===================================================================================
//        int nTotInteriorFaces          = distftot->getNel();
//        int* TotIntFaces_offsets       = distftot->getOffsets();
//
//        int nTotInteriorFaces_prism    = distfptot->getNel();
//        int* TotIntFaces_offsets_prism = distfptot->getOffsets();
        
        dimsf[0]  = nTotInteriorFaces_prism+nTotInteriorFaces_tetra+nTotBCFaces;
        dimsf[1]  = 8;
        
        filespace = H5Screate_simple(2, dimsf, NULL);
        dset_id   = H5Dcreate(file_id, "ifn",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//
        countH5[0]  = (int)ifn_P.size()/8;
        countH5[1]  = dimsf[1];
        
        offsetH5[0] = TotIntFaces_offsets_prism[world_rank];
        offsetH5[1] = 0;
        
        memspace      = H5Screate_simple(2, countH5, NULL);
        filespace     = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id      = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
//
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, &ifn_P.data()[0]);
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        countH5[0]  = (int)ifn_adaptedTetra.size()/8;
        countH5[1]  = dimsf[1];

        offsetH5[0] = nTotInteriorFaces_prism+TotIntFaces_offsets_tetra[world_rank];
        offsetH5[1] = 0;

        memspace     = H5Screate_simple(2, countH5, NULL);
        filespace     = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
//
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, &ifn_adaptedTetra.data()[0]);
        
        //===================================================================================
        int nbo = bcArrays.size();
        for(int i=0;i<bcsToT.size();i++)
        {
            int bc_id      = bcid[i];
            DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);

            int NelLoc_bci = nlbc[i];
            int NelTot_bci = distBCi->getNel();

            std::vector<int> ifn_bc_i = bcArrays[i];

            countH5[0]    = (int)ifn_bc_i.size()/8;
            countH5[1]    = 8;
            
            offsetH5[0]   = nTotInteriorFaces_prism+nTotInteriorFaces_tetra+bciTot_offsets[i]+bci_offsets[i];
            offsetH5[1]   = 0;
            memspace      = H5Screate_simple(2, countH5, NULL);
            filespace     = H5Dget_space(dset_id);

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);

            plist_id     = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        //
            status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                          plist_id, &ifn_bc_i.data()[0]);

        }
        
        
        // Create group;
        //====================================================================================
        hid_t group_zones_id  = H5Gcreate(file_id, "zones", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        // Add attribute to group:
        //====================================================================================
        dimsf_att = 1;
        att_space = H5Screate_simple(1, &dimsf_att, NULL);
        attr_id   = H5Acreate (group_zones_id, "nz", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = 3+nbo;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        //====================================================================================
        dimsf[0] = 3+nbo;
        dimsf[1] = 7;
        filespace = H5Screate_simple(2, dimsf, NULL);
        hid_t dset_zdefs_id = H5Dcreate(group_zones_id, "zdefs", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);
        
        countH5[0]  = dimsf[0];
        countH5[1]  = dimsf[1];
        offsetH5[0] = 0;
        offsetH5[1] = 0;
        memspace  = H5Screate_simple(2, countH5, NULL);
        filespace = H5Dget_space(dset_zdefs_id);
        
        //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = H5Dwrite(dset_zdefs_id, H5T_NATIVE_INT, memspace, filespace, plist_id, zdefs_new.data());

        //====================================================================================
        
        dimsf_att 			 = meshRead->znames.size();
        filespace 			 = H5Screate_simple(1, &dimsf_att, NULL);
        type      			 = H5Tcopy (H5T_C_S1);
        ret       	         = H5Tset_size (type, 20);
        ret       			 = H5Tset_strpad(type, H5T_STR_SPACEPAD);
        hid_t dset_znames_id = H5Dcreate(group_zones_id, "znames", type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Sclose(filespace);

        hsize_t cnt2 = meshRead->znames.size();
        memspace  = H5Screate_simple(1, &cnt2, NULL);
        filespace = H5Dget_space(dset_znames_id);

        std::vector<char> znames_new;

        for(int i=0;i<meshRead->znames.size();i++)
        {
            for(int j=0;j<meshRead->znames[0].size();j++)
            {
                znames_new.push_back(meshRead->znames[i][j]);
            }
        }

        status = H5Dwrite(dset_znames_id, type, memspace, filespace, plist_id, znames_new.data());

    
        
        //===================================================================================
        //===================================================================================
        //===================================================================================
    
        // clock_t t1_ijk = clock();
        // double duration = ( t1_ijk - t0_ijk) / (double) CLOCKS_PER_SEC;
        // double dur_max;
        // MPI_Allreduce(&duration, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        
        // if(world_rank==0)
        // {
        //     //std::cout << std::setprecision(16) << "Computing the metric takes " << dur_max_met << " seconds using " << world_size << " procs. " << std::endl;
        //     std::cout << std::setprecision(16) << "Redistributing the tetrahedra takes " << dur_max_redis << " seconds using " << world_size << " procs. " << std::endl;
        //     std::cout << std::setprecision(16) << "Adapting the tetrahedra takes " << dur_max_adapt << " seconds using " << world_size << " procs. " << std::endl;
        //     std::cout << std::setprecision(16) << "Writing out the grid takes " << dur_max << " seconds using " << world_size << "procs. " << std::endl;
        //     std::cout << "Finalizing process" << std::endl;     
        // }
      
    }
    
    MPI_Finalize();
        
}

