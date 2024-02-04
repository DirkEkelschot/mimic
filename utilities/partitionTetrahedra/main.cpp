
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
// #include "../../src/NekFace.h"
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
    PrismTetraTrace* pttrace = new PrismTetraTrace(comm, 
                                                   meshRead->element2rank, 
                                                   meshRead->ife,
                                                   meshRead->ifn, 
                                                   meshRead->iet, 
                                                   meshRead->nElem, 
                                                   meshRead->nFace, 
                                                   meshRead->nVert);

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
    int ntetra     = meshRead->ntetra;
    int nprism     = meshRead->nprism;
    std::map<int,std::vector<int> >::iterator itmiv;
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
            
        }
        if(eltype == 6)
        {
            prisms_e2v[elid]  = itmiv->second;
            prisms_e2f[elid]  = meshRead->ief[elid];
            prisms_e2e[elid]  = meshRead->iee[elid];
            prisms_data[elid] = meshRead->interior[elid];
        }
    }

    int nLocTetra  = tetras_e2v.size();
    int nLocPrism  = prisms_e2v.size();
    int nElemsGlob_T = 0;
    int nElemsGlob_P = 0;
    MPI_Allreduce(&nLocTetra, &nElemsGlob_T, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&nLocPrism, &nElemsGlob_P, 1, MPI_INT, MPI_SUM, comm);

    //=========END FILTERING OUT TETRA AND PRISMS FROM IO DATA STRUCTURES===============

    // we need to pass the number of verts per element in case the partition has no elements of this type.

    time_t start, end; 
 
    /* You can call it like this : start = time(NULL); 
    in both the way start contain total time in seconds 
    since the Epoch. */

    time(&start); 
    RepartitionObject* tetra_repart = new RepartitionObject(meshRead, 
                                                        tetras_e2v, 
                                                        tetras_e2f,
                                                        tetras_e2e,
                                                        pttrace, 
                                                        tetras_data,
                                                        comm);
    time(&end); 
    double time_taken = double(end - start); 
    cout << "Time taken to repartion tetrahedera is : " << fixed 
        << time_taken << setprecision(16); 
    cout << " sec " << endl;

    tetras_e2v.clear();tetras_e2f.clear();tetras_e2e.clear();



    std::map<int,std::vector<double> > loc_data_t       = tetra_repart->getElement2DataMap();
    std::map<int,std::vector<int> > gE2lV_t             = tetra_repart->getGlobalElement2LocalVertMap();
    std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
    std::map<int, std::vector<double> > LocalVertsMap_t = tetra_repart->getLocalVertsMap();
    std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();

    tetra_repart->buildInteriorSharedAndBoundaryFacesMaps(comm,pttrace, meshRead->ranges_id);
    tetra_repart->buildTag2GlobalElementMapsAndFace2LeftRightElementMap(comm,pttrace, meshRead->ranges_id);

    std::map<int,int> locv2tagvID                       = tetra_repart->getLocalVert2VertTag();
    std::map<int,int> tagv2locvID                       = tetra_repart->getVertTag2LocalVert();
    std::map<int,int> le2tagID                          = tetra_repart->getLocalElement2ElementTag();
    std::map<int,int> globalv2localvID                  = tetra_repart->getUpdatedGlobal2LocalVMap();
    std::map<int,int> oldglob2newglob_TalV              = tetra_repart->getUpdatedLocal2GlobalVMap();

    tetra_repart->buildCommunicationMaps(comm);

    std::vector<int> face4parmmg                        = tetra_repart->getFace4ParMMG(); // checked
    std::map<int,int> global2tagF                       = tetra_repart->getGlobal2TagFMap();
    int** ifc_tria_glob                                 = tetra_repart->getParMMGCommFace2GlobalVertMap();
    int** ifc_tria_loc                                  = tetra_repart->getParMMGCommFace2LocalVertMap();
    int* color_face                                     = tetra_repart->getParMMGCommColorFace();
    int *ntifc                                          = tetra_repart->getParMMGCommNFacesPerColor();
    int ncomm                                           = tetra_repart->getParMMGNComm();

    //=======================================================================================
    //=======================================================================================
    //=======================================================================================

    // std::map<int,std::vector<double> > tetra_grad = ComputedUdx_LSQ_US3D_Lite(tetra_repart, 
    //                                                                           pttrace,
    //                                                                           meshRead->ghost,
    //                                                                           meshRead->nElem,
    //                                                                           1, 
    //                                                                           1,
    //                                                                           comm);
    // time(&end); 
    // time_taken = double(end - start); 
    // cout << "Time taken to calculate gradients for tetrahedra is : " << fixed 
    //     << time_taken << setprecision(16); 
    // cout << " sec " << endl;
    // string filename_tg = "tetraGrad_" + std::to_string(world_rank) + ".vtu";

    // OutputTetraMeshPartitionVTK(filename_tg, 
    //                             Owned_Elem_t, 
    //                             gE2gV_t, 
    //                             tetra_grad, 
    //                             varnamesGrad, 
    //                             LocalVertsMap_t);


    //=======================================================================================
    //=======================================================================================
    //=======================================================================================


    //=======================================================================================
    //==================================OUTPUT TETRA PARTITIONS==============================
    //=======================================================================================
    // string filename_pgNew = "tetraData2_" + std::to_string(world_rank) + ".vtu";
    // std::map<int,std::string > varnames;
    // varnames[0]         = "TKE";
    // varnames[1]         = "Temperature";

    // OutputTetraMeshPartitionVTK(comm,
    //                             filename_pgNew, 
    //                             Owned_Elem_t_new, 
    //                             gE2gV_t, 
    //                             loc_data_t, 
    //                             varnames, 
    //                             LocalVertsMap_t);
    //=======================================================================================
    //==================================OUTPUT TETRA PARTITIONS==============================
    //=======================================================================================
    std::cout << "Initializing the data structure for ParMMG..." << std::endl;
    PMMG_pParMesh parmesh = InitializeParMMGmesh(comm, tetra_repart, pttrace, meshRead->ranges_id);

    if (inputs->niter == 0)
    {
        RunParMMGAndTestPartitioning(comm, parmesh, tetra_repart, pttrace,  meshRead->ranges_id, inputs);
    }
    else
    {   
        
        RepartitionObject* prism_repart = new RepartitionObject(meshRead, 
                                                        prisms_e2v, 
                                                        prisms_e2f,
                                                        prisms_e2e,
                                                        pttrace, 
                                                        prisms_data,
                                                        comm);

        
        prisms_e2v.clear();prisms_e2f.clear();prisms_e2e.clear();                                                
        
        prism_repart->buildInteriorSharedAndBoundaryFacesMaps(comm,pttrace, meshRead->ranges_id);
        prism_repart->buildTag2GlobalElementMapsAndFace2LeftRightElementMap(comm,pttrace, meshRead->ranges_id);

        int nNonSharedVertsOwned_P                          = prism_repart->getNonSharedVertsOwned().size();
        int nSharedVertsOwned_P                             = prism_repart->getSharedVertsOwned().size();
        DistributedParallelState* distPrismIntVerts         = new DistributedParallelState(nNonSharedVertsOwned_P,comm);
        DistributedParallelState* distPrismShaVerts         = new DistributedParallelState(nSharedVertsOwned_P,comm);
        int nVertsGlob_P                                    = distPrismIntVerts->getNel()+distPrismShaVerts->getNel();;
        
        
        // Next to all the output data structures for us3d like xcn and ifn we also need to build the tracetagV2globalV map 
        // since this is defined by the tetrahedra ordering.
        
        std::vector<int> ifn_T;
        std::map<int,int> tracetagV2globalV;
        std::map<int,std::vector<int> > bcref2bcface_T;
        std::map<int,std::vector<int> > BoundaryFaces_T;
        std::map<int,int> glob2locVid_T;
        std::map<int,int> LocationSharedVert_T;
        std::vector<std::vector<double> > new_vertices_T;
        std::vector<std::vector<int> > new_tetrahedra_T;
        std::map<int,int> oldglob2newglob_T;
        std::map<int,int> lh_T_bc;
        std::map<int,int> new_globE2locE_T;
        std::cout << "running and writing tetra data..." << std::endl;
        RunParMMGandWriteTetraUS3Dformat(comm, 
                                        parmesh, 
                                        pttrace, 
                                        inputs, 
                                        nElemsGlob_P, 
                                        nVertsGlob_P, 
                                        tracetagV2globalV, 
                                        ifn_T,
                                        bcref2bcface_T,
                                        BoundaryFaces_T,
                                        glob2locVid_T,
                                        LocationSharedVert_T,
                                        new_vertices_T,
                                        new_tetrahedra_T,
                                        oldglob2newglob_T,
                                        lh_T_bc,
                                        new_globE2locE_T);

         
        std::vector<int> ifn_P;
        std::cout << "writing prism data..." << std::endl;
        WritePrismsUS3DFormat(comm,prism_repart,pttrace,tracetagV2globalV,ifn_P,meshRead->ranges_id);
        
        int ftott = (int)ifn_T.size()/8;
        int ftotp = (int)ifn_P.size()/8;

        std::vector<std::vector<int> > bcArrays;
        std::map<int,int> bcsizing;
        std::cout << "writing boundary data..." << std::endl;
        WriteBoundaryDataUS3DFormat(comm,
                                    prism_repart,
                                    pttrace,
                                    ifn_T,
                                    ifn_P,
                                    bcref2bcface_T,
                                    unique_trace_verts2refmap,
                                    tracetagV2globalV,
                                    BoundaryFaces_T,
                                    glob2locVid_T,
                                    new_tetrahedra_T,
                                    new_vertices_T,
                                    LocationSharedVert_T,
                                    oldglob2newglob_T,
                                    lh_T_bc,
                                    new_globE2locE_T,
                                    bcArrays,
                                    bcsizing,
                                    meshRead->zdefs,
                                    meshRead->zone2bcref,
                                    meshRead->zone2name,
                                    meshRead->znames);

        
    }

    MPI_Finalize();    
}

