
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
#include "../../src/adapt_operations.h"
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;




void GetUniqueTraceData(MPI_Comm comm, 
                        RepartitionObject* partition,
                        std::vector<int> &unique_loc_trace_verts_vec,
                        std::vector<int> &unique_duplicate_loc_trace_verts_vec)
{

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::map<int,int> NonSharedVertsOwned_T                 = partition->getNonSharedVertsOwned();
    std::map<int,int> SharedVertsOwned_T                    = partition->getSharedVertsOwned();
    std::set<int> loc_trace_verts                           = partition->GetLocalTraceVertSet();
    std::map<int,std::vector<int> > shF2vert                = partition->getSharedFaceMap();
    std::map<int,std::vector<int> >::iterator its;
    std::vector<int> shared_trace_verts;
    std::vector<int> shared_trace_verts_rank;

    std::map<int,int> unique_trace_vert_map;
    //std::cout << "loc_trace_verts " << loc_trace_verts.size() << std::endl;
    for(its=shF2vert.begin();its!=shF2vert.end();its++)
    {
        int fid = its->first;
        int nv  = its->second.size();

        for(int i=0;i<nv;i++)
        {
            int vid = its->second[i];

            if(loc_trace_verts.find(vid)!=loc_trace_verts.end())
            {
                shared_trace_verts.push_back(vid);
                shared_trace_verts_rank.push_back(world_rank);
            }
        }
    }

    DistributedParallelState* distSharedTraceVerts = new DistributedParallelState(shared_trace_verts.size(),comm);
    int Ntotal_shared_trace_verts = distSharedTraceVerts->getNel();
    std::vector<int> tot_shared_trace_vert(Ntotal_shared_trace_verts,0);
    // Communicate face map to all ranks.

    MPI_Allgatherv(&shared_trace_verts.data()[0],
                    shared_trace_verts.size(),
                    MPI_INT,
                    &tot_shared_trace_vert.data()[0],
                    distSharedTraceVerts->getNlocs(),
                    distSharedTraceVerts->getOffsets(),
                    MPI_INT, comm);


    std::set<int> o_tracevert;
    for(int i=0;i<Ntotal_shared_trace_verts;i++)
    {
        int key = tot_shared_trace_vert[i];
        if(o_tracevert.find(key)==o_tracevert.end())
        {
            o_tracevert.insert(key);
            unique_duplicate_loc_trace_verts_vec.push_back(key);
        }
    }


    std::vector<int> loc_trace_verts_vec;
    int found = 0;
    std::set<int>::iterator itss;
    for(itss=loc_trace_verts.begin();itss!=loc_trace_verts.end();itss++)
    {
        int vid = *itss;

        if(o_tracevert.find(vid)==o_tracevert.end())
        {
            unique_loc_trace_verts_vec.push_back(vid);
        }
        loc_trace_verts_vec.push_back(vid);
    }

    std::set<int> global_set = AllGatherSet(loc_trace_verts,comm);

}

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
        cout << "Time taken to broadcast boundary layer/tetra trace data is : " << fixed
        << dur_max1 << setprecision(16);
        cout << " sec " << endl;
    }


    //std::map<int,std::vector<int> > trace_verts = pttrace->GetTraceVerts();

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
    bool prism_ifn = true;
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

    RepartitionObject* tetra_repart;
    std::map<int,std::vector<double> > tetra_grad;
    std::map<int,int> loc_trace2ref;

    int* new_V_offsets = new int[world_size];

    ParallelState* xcn_pstate = new ParallelState(meshRead->nVert,comm);
    for(int i=0;i<world_size;i++)
    {
        new_V_offsets[i] = xcn_pstate->getOffsets()[i]-1;
    }

    std::map<int,int> vertrefmap_pack;




    if(inputs->recursive == 0)
    {
        start = clock();
        tetra_repart = new RepartitionObject(meshRead, 
                                            tetras_e2v, 
                                            tetras_e2f,
                                            tetras_e2e,
                                            tetra2type,
                                            tetras_data,
                                            2,
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
            << dur_max << setprecision(16); 
            cout << " sec " << endl;
        }

        tetra_repart->buildUpdatedVertexAndFaceNumbering(comm, 
                                                        meshRead->ranges_id, 
                                                        meshRead->ranges_ref);

        tetra_repart->buildInteriorSharedAndBoundaryFaceMaps(comm, 
                                                            meshRead->ranges_id,
                                                            meshRead->ranges_ref);
        start = clock();

        std::map<int,std::vector<double> > Ue = tetra_repart->getElement2DataMap();
        tetra_repart->AddStateVecForAdjacentElements(Ue,2,comm);
        tetra_repart->SetStateVec(Ue,2);

        tetra_grad = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
                                                meshRead->ghost,
                                                meshRead->nElem,
                                                1,
                                                2,
                                                comm,
                                                0);
        end = clock();



        time_taken = ( end - start) / (double) CLOCKS_PER_SEC;

        MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
            cout << "Time taken to execute calculating first and second order gradients using quadratic reconstruction: " << fixed 
            << dur_max << setprecision(16); 
            cout << " sec " << endl;
        }

        std::map<int,std::string > varnamesGrad;

        varnamesGrad[0]     = "dUdx";
        varnamesGrad[1]     = "dUdy";
        varnamesGrad[2]     = "dUdz";
        varnamesGrad[3]     = "dU2dx2";
        varnamesGrad[4]     = "dU2dxy";
        varnamesGrad[5]     = "dU2dxz";
        varnamesGrad[6]     = "dU2dy2";
        varnamesGrad[7]     = "dU2dyz";
        varnamesGrad[8]     = "dU2dz2";

        string filename_tg = "mimic_hessian.vtu";
        std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();
        std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
        std::map<int, std::vector<double> > LocalVertsMap_t = tetra_repart->getLocalVertsMap();
        std::cout << "Starting to write mimic_hessian " << std::endl;
        OutputTetraMeshOnRootVTK(comm,
                                filename_tg, 
                                Owned_Elem_t, 
                                gE2gV_t, 
                                tetra_grad, 
                                varnamesGrad, 
                                LocalVertsMap_t);

        std::map<int,std::vector<double> > eigvalues;

        std::map<int,std::vector<std::vector<double> > > metric_vmap = ComputeElementMetric_Lite(comm, 
                                                                    tetra_repart,
                                                                    tetra_grad,
                                                                    eigvalues, 
                                                                    inputs);
        
        std::map<int,std::vector<std::vector<double> > >::iterator itmm;

        std::map<int,std::vector<double> > metric_vmap_new;

        for(itmm=metric_vmap.begin();itmm!=metric_vmap.end();itmm++)
        {
            int tagvid = itmm->first;

            std::vector<double> tensor(6,0);
            tensor[0] = metric_vmap[tagvid][0][0];
            tensor[1] = metric_vmap[tagvid][0][1];
            tensor[2] = metric_vmap[tagvid][0][2];
            tensor[3] = metric_vmap[tagvid][1][1];
            tensor[4] = metric_vmap[tagvid][1][2];
            tensor[5] = metric_vmap[tagvid][2][2];

            metric_vmap_new[tagvid] = tensor;

        }
        std::map<int,std::string > varnamesGrad_met;

        varnamesGrad_met[0]     = "M00";
        varnamesGrad_met[1]     = "M01";
        varnamesGrad_met[2]     = "M02";
        varnamesGrad_met[3]     = "M11";
        varnamesGrad_met[4]     = "M12";
        varnamesGrad_met[5]     = "M22";
    
        string filename_met = "mimic_metric.vtu";

        OutputTetraMeshOnRootVTK(comm,
                        filename_met, 
                        Owned_Elem_t, 
                        gE2gV_t, 
                        metric_vmap_new, 
                        varnamesGrad_met, 
                        LocalVertsMap_t);

        std::map<int,std::string > varnamesGrad_diag;

        varnamesGrad_diag[0]     = "D0";
        varnamesGrad_diag[1]     = "D1";
        varnamesGrad_diag[2]     = "D2";

    
        string filename_eig = "mimic_eigvalues.vtu";

        OutputTetraMeshOnRootVTK(comm,
                        filename_eig, 
                        Owned_Elem_t, 
                        gE2gV_t, 
                        eigvalues, 
                        varnamesGrad_diag, 
                        LocalVertsMap_t); 



    }
    else if(inputs->recursive == 1)
    {   
        //=========================================================================================
        start = clock();
        tetra_repart = new RepartitionObject(meshRead, 
                                            tetras_e2v, 
                                            tetras_e2f,
                                            tetras_e2e,
                                            tetra2type,
                                            tetras_data,
                                            1,
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
            << dur_max << setprecision(16); 
            cout << " sec " << endl;
        }

        tetra_repart->buildUpdatedVertexAndFaceNumbering(comm, 
                                                        meshRead->ranges_id, 
                                                        meshRead->ranges_ref);

        tetra_repart->buildInteriorSharedAndBoundaryFaceMaps(comm, 
                                                            meshRead->ranges_id,
                                                            meshRead->ranges_ref);


        std::map<int,std::vector<double> > Ue = tetra_repart->getElement2DataMap();
        tetra_repart->AddStateVecForAdjacentElements(Ue,2,comm);
        tetra_repart->SetStateVec(Ue,2);
        

        start = clock();
        // std::map<int,std::vector<double> > tetra_grad_v2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
        //                                                                             meshRead->ghost,
        //                                                                             meshRead->nElem,
        //                                                                             1,
        //                                                                             2,
        //                                                                             comm);

        std::map<int,std::vector<double> > tetra_grad_v2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
                                                                                        meshRead->ghost,
                                                                                        meshRead->nElem,
                                                                                        1,
                                                                                        2,
                                                                                        comm,
                                                                                        0);                                                                           

        std::map<int,std::vector<double> >::iterator iti;
        std::map<int,std::vector<double> > dudx_map;
        std::map<int,std::vector<double> > dudy_map;
        std::map<int,std::vector<double> > dudz_map;
        for(iti=tetra_grad_v2.begin();iti!=tetra_grad_v2.end();iti++)
        {
            int elid = iti->first;
            dudx_map[elid].push_back(iti->second[0]);
            dudy_map[elid].push_back(iti->second[1]);
            dudz_map[elid].push_back(iti->second[2]);
        }

        // tetra_repart->AddStateVecForAdjacentElements(dudx_map,1,comm);
        // tetra_repart->AddStateVecForAdjacentElements(dudy_map,1,comm);
        // tetra_repart->AddStateVecForAdjacentElements(dudz_map,1,comm);

        tetra_repart->AddStateVecForAdjacentElements(dudx_map,1,comm);
        tetra_repart->SetStateVec(dudx_map,1);

        // std::map<int,std::vector<double> > dU2dx2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
        //                                                                             meshRead->ghost,
        //                                                                             meshRead->nElem,
        //                                                                             0,
        //                                                                             1,
        //                                                                             comm);
        std::map<int,std::vector<double> > dU2dx2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart, 
                                                                        meshRead->ghost, 
                                                                        meshRead->nElem,
                                                                        0,1,
                                                                        comm,
                                                                        0);

        tetra_repart->AddStateVecForAdjacentElements(dudy_map,1,comm);
        tetra_repart->SetStateVec(dudy_map,1);
        // std::map<int,std::vector<double> > dU2dy2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
        //                                                                             meshRead->ghost,
        //                                                                             meshRead->nElem,
        //                                                                             0,
        //                                                                             1,
        //                                                                             comm);

        std::map<int,std::vector<double> > dU2dy2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart, 
                                                                        meshRead->ghost, 
                                                                        meshRead->nElem,
                                                                        0,1,
                                                                        comm,
                                                                        0);

        tetra_repart->AddStateVecForAdjacentElements(dudz_map,1,comm);
        tetra_repart->SetStateVec(dudz_map,1);
        // std::map<int,std::vector<double> > dU2dz2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
        //                                                                             meshRead->ghost,
        //                                                                             meshRead->nElem,
        //                                                                             0,
        //                                                                             1,
        //                                                                             comm);
        std::map<int,std::vector<double> > dU2dz2 = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart, 
                                                                meshRead->ghost, 
                                                                meshRead->nElem,
                                                                0,1,
                                                                comm,
                                                                0);

        for(iti=dU2dx2.begin();iti!=dU2dx2.end();iti++)
        {
            int elid = iti->first;
            std::vector<double> row(10,0.0);
            row[0] = dudx_map[elid][0];
            row[1] = dudy_map[elid][0];
            row[2] = dudz_map[elid][0];
            row[3] = dU2dx2[elid][0];
            row[4] = dU2dx2[elid][1];
            row[5] = dU2dx2[elid][2];
            row[6] = dU2dy2[elid][1];
            row[7] = dU2dy2[elid][2];
            row[8] = dU2dz2[elid][2];
            row[9] = Ue[elid][1];


            tetra_grad[elid] = row;
        }

        tetra_grad_v2.clear();

        dudx_map.clear();
        dudy_map.clear();
        dudz_map.clear();

        dU2dx2.clear();
        dU2dy2.clear();
        dU2dz2.clear();
        Ue.clear();
        end = clock();

        time_taken = ( end - start) / (double) CLOCKS_PER_SEC;

        MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
            cout << "Time taken to execute calculating first and second order gradients using linear reconstruction iteratively: " << fixed 
            << dur_max << setprecision(16); 
            cout << " sec " << endl;
        }

        std::map<int,std::string > varnamesGrad;

        varnamesGrad[0]     = "dUdx";
        varnamesGrad[1]     = "dUdy";
        varnamesGrad[2]     = "dUdz";
        varnamesGrad[3]     = "dU2dx2";
        varnamesGrad[4]     = "dU2dxy";
        varnamesGrad[5]     = "dU2dxz";
        varnamesGrad[6]     = "dU2dy2";
        varnamesGrad[7]     = "dU2dyz";
        varnamesGrad[8]     = "dU2dz2";
        varnamesGrad[9]     = "U";

        string filename_tg = "mimic_hessian.vtu";
        std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();
        std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
        std::map<int, std::vector<double> > LocalVertsMap_t = tetra_repart->getLocalVertsMap();

        OutputTetraMeshOnRootVTK(comm,
                                filename_tg, 
                                Owned_Elem_t, 
                                gE2gV_t, 
                                tetra_grad, 
                                varnamesGrad, 
                                LocalVertsMap_t);

        std::map<int,std::vector<double> > eigvalues;

        std::map<int,std::vector<std::vector<double> > > metric_vmap = ComputeElementMetric_Lite(comm, 
                                                                    tetra_repart,
                                                                    tetra_grad,
                                                                    eigvalues, 
                                                                    inputs);
        
        std::map<int,std::vector<std::vector<double> > >::iterator itmm;

        std::map<int,std::vector<double> > metric_vmap_new;

        for(itmm=metric_vmap.begin();itmm!=metric_vmap.end();itmm++)
        {
            int tagvid = itmm->first;

            std::vector<double> tensor(6,0);
            tensor[0] = metric_vmap[tagvid][0][0];
            tensor[1] = metric_vmap[tagvid][0][1];
            tensor[2] = metric_vmap[tagvid][0][2];
            tensor[3] = metric_vmap[tagvid][1][1];
            tensor[4] = metric_vmap[tagvid][1][2];
            tensor[5] = metric_vmap[tagvid][2][2];

            metric_vmap_new[tagvid] = tensor;

        }
        std::map<int,std::string > varnamesGrad_met;

        varnamesGrad_met[0]     = "M00";
        varnamesGrad_met[1]     = "M01";
        varnamesGrad_met[2]     = "M02";
        varnamesGrad_met[3]     = "M11";
        varnamesGrad_met[4]     = "M12";
        varnamesGrad_met[5]     = "M22";
    
        string filename_met = "mimic_metric.vtu";

        OutputTetraMeshOnRootVTK(comm,
                        filename_met, 
                        Owned_Elem_t, 
                        gE2gV_t, 
                        metric_vmap_new, 
                        varnamesGrad_met, 
                        LocalVertsMap_t);


        std::map<int,std::string > varnamesGrad_diag;

        varnamesGrad_diag[0]     = "D0";
        varnamesGrad_diag[1]     = "D1";
        varnamesGrad_diag[2]     = "D2";

    
        string filename_eig = "mimic_eigvalues.vtu";

        OutputTetraMeshOnRootVTK(comm,
                        filename_eig, 
                        Owned_Elem_t, 
                        gE2gV_t, 
                        eigvalues, 
                        varnamesGrad_diag, 
                        LocalVertsMap_t);

    
    }
    
    clock_t end_total = clock();

    end1 = clock();
    time_taken1 = ( end_total - start_total) / (double) CLOCKS_PER_SEC;
    MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        cout << "Time taken to execute MIMIC : " << fixed
        << dur_max1 << setprecision(16);
        cout << " sec " << endl;
    }
    
    MPI_Finalize();
        
}

