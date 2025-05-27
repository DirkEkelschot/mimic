
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
#include "../../src/adapt_inputs.h"
#include "../../src/adapt_gradientReconstruction.h"
#include "../../src/adapt_runparmmg.h"

#include "../../src/adapt_writeus3ddata.h"
#include "../../src/adapt_operations.h"
#include "../../src/adapt_partobject.h"
#include "../../src/adapt_prepareadaption.h"

//#include "/Users/dekelsch/mimic_libmesh/utilities/partitionTetrahedra/build/ThirdParty/dist/include/libmeshb7.h"


#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;

int main(int argc, char** argv)
{
    clock_t start_total = clock();
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

    //int64_t LibIdx;
    //int ver, dim, NmbVer, NmbTri, (*Nodes)[4], *Domains;
    //float (*Coords)[3];

    // Open the mesh file for reading
    //LibIdx = GmfOpenMesh( "triangles.meshb", GmfRead, &ver, &dim );

    // Read the egads tree stored as a raw byte flow
    // cad = GmfReadByteFlow(InpMsh, &NmbBytes);


    const char* fn_grid="inputs/grid.h5";
    const char* fn_conn="inputs/conn.h5";
    const char* fn_data="inputs/data.h5";
    // Read in the inputs from metric.xml
    Inputs* inputs = ReadXmlFile(comm, "inputs/metric.xml");
    mesh* meshRead;
    std::map<int,std::vector<int> > tetras_e2v,tetras_e2f,tetras_e2e;
    std::map<int,std::vector<double> > tetras_data;
    std::map<int,std::vector<int> > prisms_e2v,prisms_e2f,prisms_e2e;
    std::map<int,std::vector<double> > prisms_data;
    PartObject* tetra_repart;
    PartObject* prism_repart;
    PrepareAdaption* adaptionObject;
    PrepareAdaption* PrismPrepare;
    std::map<int,std::vector<double> > tetra_grad;
    clock_t start1, end1, start, end;
    double dur_max1,time_taken1;
    int nElemsGlob_P = 0;
    std::map<int,int> tetra2type;
    std::map<int,int> prism2type;   
    double dur_max,time_taken;
    bool prism_ifn;

    //===========================================================================
    // Read in the data from grid.h5/conn.h5/data.h5 in parallel using HDF5.
    // the outputted data in meshRead contains uniformly distributed data structures that will have to be partitioned.
    meshRead = ReadUS3DMeshData(fn_conn,fn_grid,fn_data,
                                    inputs->ReadFromStats,
                                    inputs->StateVar,
                                    inputs->RunNumber,
                                    comm,info);

    int Nel_loc = meshRead->ien.size();
    //Start building the trace object that contains the information regarding the prism-tetra interfaces.
    // It contains the unique vertex and face information.
    //====================================================================================

    start1 = clock();

    // PrismTetraTrace* pttrace = new PrismTetraTrace(comm, 
    //                                                meshRead->element2rank, 
    //                                                meshRead->ife,
    //                                                meshRead->ifn,
    //                                                meshRead->iet, 
    //                                                meshRead->nElem, 
    //                                                meshRead->nFace, 
    //                                                meshRead->nVert);

    // std::map<int,int> un_tracevert2ref = pttrace->GetUniqueTraceVerts2RefMap();
    
    end1 = clock();
    time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        cout << setprecision(3) << "Time taken to broadcast boundary layer/tetra trace data is :          " << fixed
        << dur_max1 << setprecision(3);
        cout << " sec " << endl;
    }


    //std::map<int,std::vector<int> > trace_verts = pttrace->GetTraceVerts();

    //Filter out the tetrahedra and prisms into seperate maps from the IO data structures (meshRead).
    //=====================================================================================


    int ntetra      = meshRead->ntetra;
    int nprism      = meshRead->nprism;
    std::map<int,std::vector<int> >::iterator itmiv;
    int foundte     = 0;
    int foundpr     = 0;

    int ndata    = meshRead->interior.begin()->second.size();

    std::map<int,double> tracePrismData;
    std::map<int,double> traceTetraData;


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
            prism2type[elid]  = eltype;

            // int nf = meshRead->ief[elid].size();

            // for(int j=0;j<nf;j++)
            // {
            //     int faceID = meshRead->ief[elid][j];

            //     if(trace_verts.find(faceID)!=trace_verts.end())
            //     {
            //         if(tracePrismData.find(elid)==tracePrismData.end())
            //         {
            //             tracePrismData[elid] = meshRead->interior[elid][1];
            //         }
            //     }
            // }
        }   
    }




    bool tetra_ifn = true;
    // bool prism_ifn = true;
    if(meshRead->elTypes[3]>0)
    {
        prism_ifn = false;
    }

    if(world_rank == 0)
    {
        std::cout << "Done filtering the element types..." << std::endl; 
    }
    
    //I am adding the prism elements and their data to the ghost map so that that data is in the boundaries data structures.
    // std::map<int,double> tracePrismData_glob = AllGatherMap_T(tracePrismData,comm);
    
    // std::map<int,double>::iterator itr;
    // for(itr=tracePrismData_glob.begin();itr!=tracePrismData_glob.end();itr++)
    // {
    //     int elid    = itr->first;
    //     double data = itr->second;

    //     std::vector<double> rowghost(2,0.0);
    //     rowghost[0] = 0.0;
    //     rowghost[1] = data;

    //     if(meshRead->ghost.find(elid)==meshRead->ghost.end())
    //     {
    //         meshRead->ghost[elid] = rowghost;
    //     }
    // }
    int nLocTetra  = tetras_e2v.size();
    int nLocPrism  = prisms_e2v.size();
    int nElemsGlob_T = 0;
    // int nElemsGlob_P = 0;
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


    




    
    std::map<int,int> loc_trace2ref;

    int* new_V_offsets = new int[world_size];

    ParallelState* xcn_pstate = new ParallelState(meshRead->nVert,comm);
    for(int i=0;i<world_size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
    }

    std::map<int,int> vertrefmap_pack;

    

    //=========================================================================================
    start = clock();
    tetra_repart = new PartObject(meshRead, 
                                tetras_e2v, 
                                tetras_e2f,
                                tetras_e2e,
                                tetras_data,
                                tetra2type,
                                3,
                                tetra_ifn,
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
        << dur_max << setprecision(3); 
        cout << " sec " << endl;
    }


    std::map<int,std::vector<int> > m_Elem2Vert         = tetra_repart->getElem2VertMap();
    std::map<int,std::vector<int> > m_Elem2Face         = tetra_repart->getElem2FaceMap();
    std::map<int,std::vector<int> > m_Face2Vert         = tetra_repart->getFace2VertMap();
    std::map<int,std::vector<int> > m_Face2Elem         = tetra_repart->getFace2ElemMap();
    std::map<int,int> m_partMap                         = tetra_repart->getPartMap();
    std::map<int,int> m_Elem2Rank                       = tetra_repart->getElem2RankMap();
    std::set<int> m_TraceVertsOnRank                    = tetra_repart->getTraceVertsOnRankMap();
    std::set<int> m_TraceFacesOnRank                    = tetra_repart->getTraceFacesOnRankMap();
    std::set<int> m_ElemSet                             = tetra_repart->getLocalElemSet();
    std::map<int,int> localV2globalV                    = tetra_repart->getLocalVert2GlobalVert();
    std::map<int,std::vector<double> > LocalVertsMap    = tetra_repart->getLocalVertsMap();

    start = clock();
    adaptionObject = new PrepareAdaption(tetra_repart, 
                                        comm,
                                        meshRead->nElem,
                                        std::move(LocalVertsMap),
                                        std::move(localV2globalV),
                                        std::move(m_Elem2Face),
                                        std::move(m_Elem2Vert),
                                        std::move(m_Face2Vert),
                                        std::move(m_Face2Elem),
                                        std::move(m_Elem2Rank),
                                        m_TraceVertsOnRank,
                                        m_TraceFacesOnRank,
                                        std::move(m_ElemSet),
                                        meshRead->ranges_id, 
                                        meshRead->ranges_ref);

    end = clock();
    time_taken = ( end - start) / (double) CLOCKS_PER_SEC;

    
    MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        cout << "Time taken to prepare for adaptation : " << fixed 
        << dur_max << setprecision(3); 
        cout << " sec " << endl;
    }

    tetra_repart->buildExtendedAdjacencyData(comm,meshRead->ghost);

    std::map<int,std::vector<double> > Ue = tetra_repart->getElem2DataMap();
    tetra_repart->AddStateVecForAdjacentElements(Ue,2,comm);
    tetra_repart->SetStateVec(Ue,2);
    start = clock();

    QRdata* qrd = Collect_QR_Data(tetra_repart, 
                                meshRead->ghost, 
                                meshRead->nElem,
                                1,2,
                                comm,
                                1);

    //std::cout << world_rank << " qrd->Nentries " << qrd->Nentries << std::endl;

    std::map<int,std::vector<double> > Amat_data = GatherJaggedGlobalMapOnRoot_T(qrd->Amat,comm);
    std::map<int,std::vector<double> > bvec_data = GatherJaggedGlobalMapOnRoot_T(qrd->bvec,comm);
    std::map<int,int> Am_data = AllGatherMap_T(qrd->Am,comm);
    std::map<int,int> An_data = AllGatherMap_T(qrd->An,comm);
    if(world_rank == 0)  
    {
        ofstream myfile;
        ofstream myfile3;
        ofstream myfile4;
        myfile.open("Amats.dat");
        myfile3.open("elID.dat");
        myfile4.open("bvecs.dat");
        std::cout << "glob_data  " << Amat_data.size() << " ntetra " << nElemsGlob_T << std::endl;
        std::map<int,std::vector<double> >::iterator itd;

        for(itd=Amat_data.begin();itd!=Amat_data.end();itd++)
        {
            int elid = itd->first;
            myfile3 << itd->first;
            for(int q=0;q<itd->second.size();q++)
            {
                myfile << itd->second[q] << " ";
            }
            for(int q=0;q<bvec_data[itd->first].size();q++)
            {
                myfile4 << bvec_data[itd->first][q] << " ";
            }
            myfile  << (double) Am_data[elid] << " ";
            myfile  << (double) An_data[elid] << " ";
            myfile  << std::endl;
            myfile3 << std::endl;
            myfile4 << std::endl;

        }
        myfile.close();
        myfile3.close();
        myfile4.close();

    }


    MPI_Finalize();
        
}

