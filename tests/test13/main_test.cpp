#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"

#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// This is basically textbook recursive merge sort using std::merge_inplace
// but it considers the offsets of segments that are already sorted



std::map<int,int> AllGatherMap(std::map<int,int> mappie, MPI_Comm mpi_comm)
{
    int mapSizeLoc = mappie.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot = distrimap->getNel();
    int* key_loc = new int[mapSizeLoc];
    int* val_loc = new int[mapSizeLoc];
    int* key_tot = new int[mapSizeTot];
    int* val_tot = new int[mapSizeTot];
    int i = 0;
    
    std::map<int,int>::iterator itred;
    for(itred=mappie.begin();itred!=mappie.end();itred++)
    {
        key_loc[i] = itred->first;
        val_loc[i] = itred->second;
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    
    MPI_Allgatherv(key_loc,
                   mapSizeLoc,
                   MPI_INT,
                   key_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    
    MPI_Allgatherv(val_loc,
                   mapSizeLoc,
                   MPI_INT,
                   val_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    int key,val;
    std::map<int,int> mappie_glob;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = key_tot[i];
        val = val_tot[i];
        
        if(mappie_glob.find(key)==mappie_glob.end())
        {
        	mappie_glob[key] = val;
        }
    }
    
    return mappie_glob;
}







//void OutputTetrahedralMeshOnPartition(TetrahedraMesh* tmesh, MPI_Comm comm)
//{
//
//    int world_size;
//    MPI_Comm_size(comm, &world_size);
//    // Get the rank of the process
//    int world_rank;
//    MPI_Comm_rank(comm, &world_rank);
//
//    std::vector<int> lverts;
//    std::map<int,int> lpartv2gv_v2;
//    std::map<int,int> gv2lpv2;
//
//    std::set<int> gv_set;
//    int lcv2 = 0;
//    Array<int>* ien_part_tetra     = tmesh->ien_part_tetra;
//    Array<int>* ien_part_hybrid    = tmesh->ien_part_hybrid;
//    std::vector<Vert*> locVs       = tmesh->LocalVerts;
//    int nElonRank = ien_part_tetra->getNrow();
//
//    Array<int>* locelem2locnode= new Array<int>(nElonRank,4);
//
//    std::vector<Vert*> printVs;
//
//    for(int i=0;i<ien_part_tetra->getNrow();i++)
//    {
//        for(int q=0;q<ien_part_tetra->getNcol();q++)
//        {
//            int gv = ien_part_tetra->getVal(i,q);
//            int lvv = tmesh->globV2locV[gv];
//
//            if(gv_set.find(gv)==gv_set.end())
//            {
//                gv_set.insert(gv);
//                lverts.push_back(lvv);
//                lpartv2gv_v2[lvv]=gv;
//                gv2lpv2[gv]=lcv2;
//                locelem2locnode->setVal(i,q,lcv2);
//
//                printVs.push_back(locVs[lvv]);
//
//                lcv2=lcv2+1;
//            }
//            else
//            {
//                int lcv_u = gv2lpv2[gv];
//                locelem2locnode->setVal(i,q,lcv_u);
//            }
//        }
//    }
//
//    std::vector<Vert*> lv = tmesh->LocalVerts;
//    std::string filename = "checkPart_" + std::to_string(world_rank) + ".dat";
//    std::ofstream myfile;
//    myfile.open(filename);
//    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
//    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
//    myfile <<"ZONE N = " << printVs.size() << ", E = " << nElonRank << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
//
//    for(int i=0;i<printVs.size();i++)
//    {
//        myfile << printVs[i]->x << " " << printVs[i]->y << " " << printVs[i]->z << std::endl;
//    }
//    int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
//    int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
//    for(int i=0;i<ien_part_hybrid->getNrow();i++)
//    {
//        myfile <<   locelem2locnode->getVal(i,0)+1 << "  " <<
//        locelem2locnode->getVal(i,1)+1 << "  " <<
//        locelem2locnode->getVal(i,2)+1 << "  " <<
//        locelem2locnode->getVal(i,3)+1 << "  " << std::endl;
//    }
//
//
//    myfile.close();
//}












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
    
    int ier,opt;
    int debug = 1;
    const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
    const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
    const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    const char* fn_metric = "metric.inp";
    
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
    
    int ReadFromStats = 0;
    if(metric_inputs.size()==6)
    {
        ReadFromStats=metric_inputs[5];
    }
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve = us3d->xcn->getNglob();
    
    int Nel_part = us3d->ien->getNrow();
    
    //Array<double>* xcn_ref = ReadDataSetFromFile<double>(fn_grid,"xcn");
    //Array<int>* ien_ref    = ReadDataSetFromFile<int>(fn_conn,"ien");

    Array<double>* Ui = new Array<double>(Nel_part,1);
    int varia = 4;
    double rhoState,uState,vState,wState,TState,VtotState,aState,MState;
    for(int i=0;i<Nel_part;i++)
    {
        rhoState  = us3d->interior->getVal(i,0);
        uState    = us3d->interior->getVal(i,1);
        vState    = us3d->interior->getVal(i,2);
        wState    = us3d->interior->getVal(i,3);
        TState    = us3d->interior->getVal(i,4);
        VtotState = sqrt(uState*uState+vState*vState+wState*wState);
        aState    = sqrt(1.4*287.05*TState);
        MState    = VtotState/aState;
        Ui->setVal(i,0,MState);
    }
    
    delete us3d->interior;

    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,comm);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
      
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, us3d->ie_tetCnt, comm);

    std::vector<int> LocElem = P->getLocElem();
    std::vector<double> Uvaria  = P->getLocElemVaria();
    std::map<int,Array<double>*> Uvaria_map;
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid = LocElem[i];
        UvariaV   = Uvaria[i];
        Array<double>* Uarr = new Array<double>(1,1);
        Uarr->setVal(0,0,UvariaV);
        Uvaria_map[gid] = Uarr;
    }
    
    Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);
    
    std::map<int,double> Volumes = meshTopo->getVol();
    P->AddStateVecForAdjacentElements(Uvaria_map,1,comm);
    
    std::map<int,Array<double>* > var_vmap = P->ReduceStateVecToAllVertices(Uvaria_map,1);
    std::map<int,Array<double>* >::iterator vm;
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }
    
    delete us3d->ghost;
    
    
    //std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,Uvaria_map,gB,comm);
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_Vrt_US3D(P,Uvaria_map,var_vmap,meshTopo,gB,comm);
    
    P->AddStateVecForAdjacentElements(dUdXi,3,comm);
    
    std::map<int,Array<double>* >::iterator grit;
    std::map<int,Array<double>* > dUidxi_map;
    std::map<int,Array<double>* > dUidyi_map;
    std::map<int,Array<double>* > dUidzi_map;
    std::map<int,Array<double>* > dUidXi_map;
    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
        
        Array<double>* dUdx_E = new Array<double>(1,1);
        dUdx_E->setVal(0,0,grit->second->getVal(0,0));
        Array<double>* dUdy_E = new Array<double>(1,1);
        dUdy_E->setVal(0,0,grit->second->getVal(1,0));
        Array<double>* dUdz_E = new Array<double>(1,1);
        dUdz_E->setVal(0,0,grit->second->getVal(2,0));
        
        dUidxi_map[grit->first]=dUdx_E;
        dUidyi_map[grit->first]=dUdy_E;
        dUidzi_map[grit->first]=dUdz_E;
        
        delete grit->second;
    }

    std::map<int,Array<double>* > dudx_vmap = P->ReduceStateVecToAllVertices(dUidxi_map,1);
    std::map<int,Array<double>* > dudy_vmap = P->ReduceStateVecToAllVertices(dUidyi_map,1);
    std::map<int,Array<double>* > dudz_vmap = P->ReduceStateVecToAllVertices(dUidzi_map,1);


    //std::cout << "second gradient "<<std::endl;
//    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_US3D(P,dUidxi_map,gB,comm);
//    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_US3D(P,dUidyi_map,gB,comm);
//    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_US3D(P,dUidzi_map,gB,comm);

    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidxi_map,dudx_vmap,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidyi_map,dudy_vmap,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidzi_map,dudz_vmap,meshTopo,gB,comm);

//      Array<double>* dU2dXi2 = ComputedUdx_MGG(P,dUdxauxNew,meshTopo,gB,comm);
//      Array<double>* dU2dYi2 = ComputedUdx_MGG(P,dUdyauxNew,meshTopo,gB,comm);
//      Array<double>* dU2dZi2 = ComputedUdx_MGG(P,dUdzauxNew,meshTopo,gB,comm);
            
    std::map<int,Array<double>*> Hess_map;
    //std::cout << "second gradient 2"<<std::endl;

    std::map<int,Array<double>* >::iterator itgg;
    int te = 0;
    
    for(itgg=dU2dXi2.begin();itgg!=dU2dXi2.end();itgg++)
    {
        int gid = itgg->first;
        
        Array<double>* Hess = new Array<double>(6,1);
        
        Hess->setVal(0,0,dU2dXi2[gid]->getVal(0,0));
        Hess->setVal(1,0,dU2dXi2[gid]->getVal(1,0));
        Hess->setVal(2,0,dU2dXi2[gid]->getVal(2,0));

        Hess->setVal(3,0,dU2dYi2[gid]->getVal(1,0));
        Hess->setVal(4,0,dU2dYi2[gid]->getVal(2,0));
        Hess->setVal(5,0,dU2dZi2[gid]->getVal(2,0));
        
        Hess_map[gid] = Hess;
        
        delete dU2dXi2[gid];
        delete dU2dYi2[gid];
        delete dU2dZi2[gid];
        
        t++;
    }
    
    dU2dXi2.clear();
    dU2dYi2.clear();
    dU2dZi2.clear();
    
    
    double* Hessie = new double[9];
    double * WRn = new double[3];
    Array<double>* DR  = new Array<double>(3,3);
    Array<double>* UR  = new Array<double>(3,3);
    //+++++++++++++++++++++++++++++++++++++++++++
    //++++  Scaling eigenvalues/eigenvectors ++++
    double hmin         = metric_inputs[1];
    double hmax         = metric_inputs[2];
    double f            = metric_inputs[3];
    //+++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++
    double cmplxty = 0.0;
    double cmplxty2 = 0.0;
    double Volu=0.0,cmplxty_red=0.0,cmplxty_tmp=0.0,cmplxty_tmp2=0.0;
    double po = 6.0;
    for(int j=0;j<3;j++)
    {
        for(int k=0;k<3;k++)
        {
            DR->setVal(j,k,0.0);
        }
    }
    for(itgg=Hess_map.begin();itgg!=Hess_map.end();itgg++)
    {
        Hessie[0] = itgg->second->getVal(0,0);
        Hessie[1] = itgg->second->getVal(1,0);
        Hessie[2] = itgg->second->getVal(2,0);
        
        Hessie[3] = itgg->second->getVal(1,0);
        Hessie[4] = itgg->second->getVal(3,0);
        Hessie[5] = itgg->second->getVal(4,0);
        
        Hessie[6] = itgg->second->getVal(2,0);
        Hessie[7] = itgg->second->getVal(4,0);
        Hessie[8] = itgg->second->getVal(5,0);
        
        Eig* eig = ComputeEigenDecomp(3, Hessie);
        
        WRn[0] = std::min(std::max(f*fabs(eig->Dre[0]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[1] = std::min(std::max(f*fabs(eig->Dre[1]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        WRn[2] = std::min(std::max(f*fabs(eig->Dre[2]),1.0/(hmax*hmax)),1.0/(hmin*hmin));
        
        DR->setVal(0,0,WRn[0]);DR->setVal(0,1,0.0);DR->setVal(0,1,0.0);
        DR->setVal(1,0,0.0);DR->setVal(1,1,WRn[1]);DR->setVal(1,2,0.0);
        DR->setVal(2,0,0.0);DR->setVal(2,1,0.0);DR->setVal(2,2,WRn[2]);
        
        UR->setVal(0,0,eig->V[0]);UR->setVal(0,1,eig->V[1]);UR->setVal(0,2,eig->V[2]);
        UR->setVal(1,0,eig->V[3]);UR->setVal(1,1,eig->V[4]);UR->setVal(1,2,eig->V[5]);
        UR->setVal(2,0,eig->V[6]);UR->setVal(2,1,eig->V[7]);UR->setVal(2,2,eig->V[8]);

        Array<double>* iVR = MatInv(UR);
        Array<double>* Rs = MatMul(UR,DR);
        Array<double>* Rf = MatMul(Rs,iVR);
        double detRf = sqrt( Rf->getVal(0,0)*(Rf->getVal(1,1)*Rf->getVal(2,2)-Rf->getVal(2,1)*Rf->getVal(1,2))
                            -Rf->getVal(0,1)*(Rf->getVal(1,0)*Rf->getVal(2,2)-Rf->getVal(2,0)*Rf->getVal(1,2))
                            +Rf->getVal(0,2)*(Rf->getVal(1,0)*Rf->getVal(2,1)-Rf->getVal(2,0)*Rf->getVal(1,1)));
        
        Volu = Volumes[itgg->first];
        cmplxty_tmp = detRf;
        cmplxty_tmp2 = Volu;
        //std::pow(cmplxty_tmp,(po+1.0)/(2.0*po+3.0));
        cmplxty = cmplxty + cmplxty_tmp;
        cmplxty2 = cmplxty2 + cmplxty_tmp2;
        delete iVR;
        delete Rs;
        delete Rf;
    }
    

    MPI_Allreduce(&cmplxty, &cmplxty_red, 1, MPI_DOUBLE, MPI_SUM, comm);
    //cmplxty_red = std::pow(cmplxty_red,(po+1)/(2.0*po+3.0));

    P->AddStateVecForAdjacentElements(Hess_map,6,comm);

    std::map<int,Array<double>* > hess_vmap = P->ReduceStateVecToAllVertices(Hess_map,6);

    cmplxty_red = Nve/cmplxty_red;
    
    ComputeMetric(P, metric_inputs, comm, hess_vmap, 1.0, po);
    
    //std::cout << "understand bitshift " << world_rank << "  " << (world_rank & (~(1 << 3)))<< std::endl;
    //std::vector<int> LocElem     = P->getLocElem();

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Domain* pDom = P->getPartitionDomain();
    std::map<int,std::vector<int> > tetras     = pDom->GTetras;
    std::map<int,std::vector<int> > prisms     = pDom->GPrisms;
    std::map<int,std::vector<int> > ushell     = pDom->ushell;
    i_part_map* if_Nv_part_map                 = P->getIF_Nvpartmap();
    i_part_map* ifn_part_map                   = P->getIFNpartmap();
    i_part_map* ife_part_map                   = P->getIFEpartmap();
    i_part_map* ief_part_map                   = P->getIEFpartmap();
    i_part_map* ien_part_map                   = P->getIENpartmap();
    i_part_map* if_ref_part_map                = P->getIFREFpartmap();
    Array<int>* part_global                    = P->getGlobalPartition();
    //std::map<int,Array<double>* > M_vmap;
    // Based on a partitioned hybrid mesh, we are extracting the tetra and redistribute them
    // uniformly. The hybrid mesh is partitioned without any weights which means that it might happen that the tetrahedra are initial distributed non-uniformly where a significant number of processors end up with any tetrahedra. This function distributed them equally.
    //std::cout << "Starting extracting..." << tetras.size() <<std::endl;
    //std::cout << "n_prism = " << prisms.size() << " " << tetras.size() << std::endl;
    // local face2vert_map for a prism {0,1,2},{3,5,4},{0,2,4,3},{1,5,4,2},{0,3,5,1} };
    // rowmap {3,3,4,4,4}
    
    
//    std::vector<std::vector<int> > prism_faces(4);
//    prism_faces[0].push_back(0);prism_faces[0].push_back(1);prism_faces[0].push_back(2);
//    prism_faces[1].push_back(3);prism_faces[1].push_back(5);prism_faces[1].push_back(4);
//
//    prism_faces[2].push_back(0);prism_faces[2].push_back(2);prism_faces[2].push_back(4);prism_faces[2].push_back(3);
//    prism_faces[3].push_back(1);prism_faces[3].push_back(5);prism_faces[3].push_back(4);prism_faces[3].push_back(2);
//    prism_faces[4].push_back(0);prism_faces[4].push_back(3);prism_faces[4].push_back(5);prism_faces[4].push_back(1);
    
//    std::map<int,std::vector<int> >::iterator itmm;
//    for(itmm=prisms.begin();itmm!=prism.end();itmm++)
//    {
//        pElid                   = itmm->first;
//        std::vector<int> prism  = itmm->second;
//
//        for(int p=0;p<prism_faces.size();p++)
//        {
//            int nvf = prism_faces[p].size();
//            std::vector<int> fce(nvf);
//            std::set<int> fceset;
//            for(int s=0;s<nvf;s++)
//            {
//                fce[s] = prism[prism_faces[p][s]];
//                fceset.insert(prism[prism_faces[p][s]]);
//            }
//            fceset.clear();
//
//        }
//    }
    
    
//    RedistributePartitionObject* prism_distri = new RedistributePartitionObject(us3d,part_global,
//                                                                                prisms,
//                                                                                ief_part_map->i_map,
//                                                                                ifn_part_map->i_map,
//                                                                                ife_part_map->i_map,
//                                                                                if_ref_part_map->i_map,
//                                                                                ushell,
//                                                                                hess_vmap, comm);
    
    
    RedistributePartitionObject* tetra_distri = new RedistributePartitionObject(us3d,part_global,
                                                                                tetras,
                                                                                ief_part_map->i_map,
                                                                                ifn_part_map->i_map,
                                                                                ife_part_map->i_map,
                                                                                if_ref_part_map->i_map,
                                                                                ushell,
                                                                                hess_vmap, comm);
    
    Array<int>* element2node                        = tetra_distri->GetElement2NodeMap();
    std::map<int,Array<double>* > metric            = tetra_distri->GetVert2MetricMap();
    int** ifc_tria_glob                             = tetra_distri->GetFace2GlobalNode();
    int** ifc_tria_loc                              = tetra_distri->GetFace2LocalNode();
    int nFaces                                      = tetra_distri->GetNBoundaryFaces();
    std::vector<Vert*> locVs                        = tetra_distri->GetLocalVertices();
    std::vector<int> faces4parmmg                   = tetra_distri->GetFaces4ParMMG();
    std::map<int,int*> face2node                    = tetra_distri->GetFace2NodeMap();
    std::map<int,std::vector<int> > face2element    = tetra_distri->GetFace2ElementMap();
    std::map<int,int> globV2locV                    = tetra_distri->GetGlobalVert2LocalVertMap();
    std::map<int,int> locV2globV                    = tetra_distri->GetLocalVert2GlobalVertMap();
    int ncomm                                       = tetra_distri->GetNcomm();
    int* color_face                                 = tetra_distri->GetColorFace();
    //int** face2globnode                           = tetra_distri->GetFace2GlobalNode();
    int *ntifc                                      = tetra_distri->GetNFacesPerColor();
    std::map<int,int> locShF2globShF                = tetra_distri->GetLocalSharedFace2GlobalSharedFace();
    std::map<int,int> face2ref                      = tetra_distri->GetFace2RefMap();
    std::map<int,int> shell_tet2hybF                = tetra_distri->GetShellTet2HybFaceMap();
    std::map<int,int> shellvert2ref                 = tetra_distri->GetShellVert2RefMap_Global();
    std::map<int,std::set<int> > shellface2vertref  = tetra_distri->GetShellFace2VertRefMap();

    std::map<int,int> tetF2hybF                     = tetra_distri->GetTetF2HybFMap();
    std::map<int,int> tetV2tagV						= tetra_distri->GetTet2TagVertMap();
    
    std::map<int,int> shellvertTag2ref;
    
    std::map<int,int>::iterator itt;
    for(itt=shellvert2ref.begin();itt!=shellvert2ref.end();itt++)
    {
    	int vtet = itt->first;
    	int vhyb = tetV2tagV[vtet];
    	int ref  = itt->second;
    	shellvertTag2ref[vhyb] = ref;
    }
    
    
    std::map<int,int> shellvertOriginalTag2ref_Glob = AllGatherMap(shellvert2ref,comm);
    
    
    //std::cout << " shell_Vert2Ref " << world_rank << " " << shellvertOriginalTag2ref_Glob.size() << " " << shellvert2ref.size() << std::endl;
    
    //int icomm;
    //Based on the new local tetrahedra mesh, we output a tecplot file per processor that has the geometry of the computational domain that is owned by world_rank.
    
    //OutputTetrahedralMeshOnPartition(tmesh,comm);
    
    Array<int>* ien_part_tetra   = element2node;
    int nTetrahedra              = element2node->getNrow();
    int nVertices                = locVs.size();
    int nTriangles               = nFaces;
    int nEdges                   = 0;
    int nPrisms                  = 0;
    int nQuadrilaterals          = 0;
    
    //std::cout << " shell_tet2hybF.size --> " << shell_tet2hybF.size() << std::endl;

    std::map<int,int> tag2locV_map = P->getGlobalVert2LocalVert();
    std::vector<Vert*> localVsPartition = P->getLocalVerts();
    
    newNumberingNodesFaces* nnf = DetermineNewNumberingOfElementSubset_Test(part_global,
                                                                       prisms,
                                                                       ief_part_map->i_map,
                                                                       ifn_part_map->i_map,
                                                                       ife_part_map->i_map,
                                                                       if_ref_part_map->i_map,
                                                                       if_Nv_part_map->i_map,
                                                                       ushell,
																	   tag2locV_map,
																	   localVsPartition,
                                                                       shellvertOriginalTag2ref_Glob,
                                                                       comm);
    
    
    nPrisms                 = nnf->ien.size();
    int nPrismFaces         = nnf->face2node.size();
    std::map<int,int> ifref = nnf->ifref;
    std::map<int,int> tagE2gE = nnf->tagE2gE;
    std::map<int,int> gE2tagE = nnf->gE2tagE;
    DistributedParallelState* distPrismIN      = new DistributedParallelState(nPrisms,comm);
    int nPrismsTot = distPrismIN->getNel();
    std::vector<std::vector<int> > prism_faces;
    std::map<int,std::vector<int> > face2node_prisms        = nnf->face2node;
    std::map<int,int > locface2globface                     = nnf->loc2globF;
    std::map<int,int > globface2locface                     = nnf->glob2locF;
    std::map<int,int > locface2globface_sha                 = nnf->loc2globF_sha;
    std::map<int,int > globface2locface_sha                 = nnf->glob2locF_sha;
    std::map<int,int> rhp                                   = nnf->rh;
    std::map<int,int> lhp                                   = nnf->lh;
    std::map<int,int> rhp_sha                               = nnf->rh_sha;
    std::map<int,int> lhp_sha                               = nnf->lh_sha;
    std::map<int,std::vector<int> > pbcmap                  = nnf->pbcmap;
    std::map<int,int> tag2element                           = nnf->tag2ElementID;
    std::map<int,std::vector<int> > shell_face2node_prism   = nnf->shellFace2Node;
    std::map<int,std::vector<int> > shared_face2node_prism  = nnf->sharedFace2Node;
    std::map<int,std::vector<int> > int_face2node_prism     = nnf->intFace2Node;
    std::map<int,int> tag2glob_prism = nnf->tag2glob_prism;
    //std::cout << "int_face2node_prism size just after being created " << int_face2node_prism.size() << std::endl;
    
    std::map<int,std::vector<int> > bc_face2node_prism      = nnf->bcFace2Node;
    int nTotUniqueNonShellVerts						        = nnf->nTotUniqueNonShellVerts;
    std::map<int,int> tagV2localV_prism                     = nnf->tagV2localV;
    std::map<int,int> localV2tagV_prism                     = nnf->localV2tagV;
    std::map<int,int> sharedvert2rank_prism                 = nnf->sharedvert2rank;
    Array<int>* ifnOUT_prism = new Array<int>(int_face2node_prism.size()+shared_face2node_prism.size(),8);
    Array<int>* parmmg_iet_prisms                           = nnf->iet;
    // std::cout << bcFace2Node_prism.size() << " " << nnf->rh.size() << " " << nnf->lh.size() << std::endl;
    Array<double>* xcn_prisms_int                               = nnf->xcn_int;
    Array<double>* xcn_prisms_shared                               = nnf->xcn_shared;


    
    std::cout << "int_face2node_prism.size() and shared_face2node_prism.size()  " << int_face2node_prism.size() << " " <<shared_face2node_prism.size() << " " << world_rank << std::endl;
    std::map<int,std::vector<int> >::iterator itmm;
    int foundU = 0;
    
    std::map<std::set<int>, int > vertref2shell_prism;
    int nshell = 0;
    for(itmm=ushell.begin();itmm!=ushell.end();itmm++)
    {
        int fhyb      = itmm->first;
        int Etettag   = itmm->second[0];
        int Eprismtag = itmm->second[1];
        int EprismNew = tag2element[Eprismtag];
        
        if(shellface2vertref.find(fhyb)!=shellface2vertref.end())
        {
            vertref2shell_prism[shellface2vertref[fhyb]] = EprismNew;
            foundU++;
        }
        nshell++;
    }

    int nShellVertsTot = shellvertOriginalTag2ref_Glob.size();
    int nLocVertPrism  = tagV2localV_prism.size();
    
	DistributedParallelState* dist_Vprism = new DistributedParallelState(nLocVertPrism,comm);
	int nTotPrismVerts          = dist_Vprism->getNel();
        
    int nTetVertsOffset_prism   = nTotUniqueNonShellVerts;
    
    //std::cout << "FoundU =============>>>>>>>>>>>  " << shellface2vertref.size() << " " << vertref2shell_prism.size() << " " << nshell << " " << ushell.size() << " " << foundU << " " << nShellVertsTot << " " << dist_Vprism->getNel() << std::endl;
    //std::cout << "wfafkadslka " << world_rank << " " << nTotPrismVerts << " " << tagV2localV_prism.size() << " " << localV2tagV_prism.size() << " " << nTotUniqueNonShellVerts << " " << nTetVertsOffset_prism << std::endl;
    
    
        
        std::map<int,std::vector<int> >::iterator prit;
        int pid   = 0;
        int pfid  = 0;
        int f13   = 0;
        int fptot = 0; //   IntFprims_offsets[world_rank];
        int llvid = 0;
        std::map<int,int> fftell;
        int ftell = 0;
        int gftel = 0;
        int lbb = 1;
        std::map<int,int> Tag2globPrimsVid;
        std::map<int,int> allPnodes;
        int ggvid;
        std::map<int,int> local2globalVertMap = nnf->local2globalVertMap;
        
        DistributedParallelState* distriPrismVertmap = new DistributedParallelState(local2globalVertMap.size(),comm);
        int cc       = 0;
        int notFound = 0;
        int tagnew   = 0;
        std::map<int,int> t2g_tmp;
        std::map<int,int> testmap;
        int notfo  = 0;
        int fo     = 0;
        int notfo2 = 0;
        
        std::map<int,int> sharedVmap = nnf->sharedVmap;
        int redflag = 0;
        int inbitch = 0;
        int swbitch = 0;
        int inbitchint = 0;
    
    
        std::map<int,int> SharedVertsOwned     =  nnf->SharedVertsOwned;
        std::map<int,int> NonSharedVertsOwned  =  nnf->NonSharedVertsOwned;
        std::map<int,int> SharedVertsNotOwned  =  nnf->SharedVertsNotOwned;
    
        for( prit=int_face2node_prism.begin();prit!=int_face2node_prism.end();prit++)
        {
            int gfid    = prit->first;
            int npf     = prit->second.size();
            std::vector<int> fce(npf);
            ifnOUT_prism->setVal(fptot,0,npf);
            
            std::vector<Vert*> Vfaces;
            Vert* VcF = new Vert;
            
            for(int g=0;g<npf;g++)
            {
                int oldtag = prit->second[g];
                
                if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = localVsPartition[lvp]->x;
                    Vf->y = localVsPartition[lvp]->y;
                    Vf->z = localVsPartition[lvp]->z;
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = tag2glob_prism[oldtag];
                    fce[g]     = globid;
                
                }
                else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = localVsPartition[lvp]->x;
                    Vf->y = localVsPartition[lvp]->y;
                    Vf->z = localVsPartition[lvp]->z;
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = SharedVertsNotOwned[oldtag];
                    fce[g]     = globid;
                    
                }
            }
            
            VcF->x = VcF->x/npf;
            VcF->y = VcF->y/npf;
            VcF->z = VcF->z/npf;
            
//            if(npf == 3)
//            {
//                ifnOUT_prism->setVal(fptot,4,0);
//            }
            
            if(rhp[gfid] == 0 || lhp[gfid] == 0)
            {
                std::cout << world_rank << " rhp[gfid] and lhp[gfid] are zero " << std::endl;
            }
            
            int leftEl  = lhp[gfid];
            int leftTag = gE2tagE[leftEl];
            Vert* Vijk = new Vert;
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            int nvp = prisms[leftTag].size();

            for(int q=0;q<nvp;q++)
            {
                int tag  = prisms[leftTag][q];
                int lvp  = tag2locV_map[tag];

                Vijk->x = Vijk->x + localVsPartition[lvp]->x;
                Vijk->y = Vijk->y + localVsPartition[lvp]->y;
                Vijk->z = Vijk->z + localVsPartition[lvp]->z;
            }

            Vijk->x = Vijk->x/nvp;
            Vijk->y = Vijk->y/nvp;
            Vijk->z = Vijk->z/nvp;

            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            
            if(orient0 < 0.0)
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[2]);
                    ifnOUT_prism->setVal(fptot,3,fce[1]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[3]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[1]);
                }
                
            }
            else
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[3]);
                }
            }
//            if(npf == 4 && ifnOUT_prism->getVal(fptot,4)==0)
//            {
//                std::cout << "INT wrong" <<std::endl;
//            }
            if(lhp[gfid]==4 || rhp[gfid] == 4)
            {
                std::cout << "Int wrong " << lhp[gfid] << " " << rhp[gfid] << " " << orient0 << " Vijk " << Vijk->x << ", " << Vijk->y << ", " << Vijk->z << " VcF " << VcF->x << ", " << VcF->y << ", " << VcF->z << std::endl;
                
                
                std::cout << " = [";
                for(int q=0;q<Vfaces.size();q++)
                {
                    std::cout << "[" << Vfaces[q]->x << ", " << Vfaces[q]->y << ", " << Vfaces[q]->z << "], "<<std::endl;
                }
                std::cout << "]"<<std::endl;
            }

            ifnOUT_prism->setVal(fptot,5,rhp[gfid]);
            ifnOUT_prism->setVal(fptot,6,lhp[gfid]);
            ifnOUT_prism->setVal(fptot,7,2);

            fptot++;
            pid++;
        }
    
        //std::cout << " notfo " << notfo2 << " " << notfo <<  " " << fo << " " << fptot << " " << int_face2node_prism.size() << " " << world_rank << std::endl;
            
        int cnttt     = 0;
        int inbitchsh = 0;
        int notany = 0;
        for( prit=shared_face2node_prism.begin();prit!=shared_face2node_prism.end();prit++)
        {
            int gfid  = prit->first;
            int npf   = prit->second.size();
            
            ifnOUT_prism->setVal(fptot,0,npf);
            std::vector<int> fce(npf);
            std::vector<Vert*> Vfaces;
            Vert* VcF = new Vert;
            for(int g=0;g<npf;g++)
            {
                int oldtag = prit->second[g];
                if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = localVsPartition[lvp]->x;
                    Vf->y = localVsPartition[lvp]->y;
                    Vf->z = localVsPartition[lvp]->z;
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = tag2glob_prism[oldtag];
                    fce[g]     = globid;
                    
                }
                else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = localVsPartition[lvp]->x;
                    Vf->y = localVsPartition[lvp]->y;
                    Vf->z = localVsPartition[lvp]->z;
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = SharedVertsNotOwned[oldtag];
                    fce[g]     = globid;
                    
                }
            }
            
            VcF->x = VcF->x/npf;
            VcF->y = VcF->y/npf;
            VcF->z = VcF->z/npf;
            if(npf == 3)
            {
                ifnOUT_prism->setVal(fptot,4,0);
            }
            
            int leftEl  = lhp[gfid];
            int leftTag = gE2tagE[leftEl];
            
            Vert* Vijk = new Vert;
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            int nvp = prisms[leftTag].size();

            for(int q=0;q<nvp;q++)
            {
                int tag  = prisms[leftTag][q];
                int lvp  = tag2locV_map[tag];

                Vijk->x = Vijk->x + localVsPartition[lvp]->x;
                Vijk->y = Vijk->y + localVsPartition[lvp]->y;
                Vijk->z = Vijk->z + localVsPartition[lvp]->z;
            }

            Vijk->x = Vijk->x/nvp;
            Vijk->y = Vijk->y/nvp;
            Vijk->z = Vijk->z/nvp;

            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
//
            if(orient0 < 0.0)
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[2]);
                    ifnOUT_prism->setVal(fptot,3,fce[1]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[3]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[1]);
                }
            }
            else
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[3]);
                }
            }
            
            
            ifnOUT_prism->setVal(fptot,5,rhp[gfid]);
            ifnOUT_prism->setVal(fptot,6,lhp[gfid]);
            ifnOUT_prism->setVal(fptot,7,2);
            
                if(lhp[gfid]==4 || rhp[gfid] == 4)
                {
                    std::cout << "Shared wrong " << lhp[gfid] << " " << rhp[gfid] << " " << orient0 << " Vijk " << Vijk->x << ", " << Vijk->y << ", " << Vijk->z << " VcF " << VcF->x << ", " << VcF->y << ", " << VcF->z << std::endl;
                    
                    
                    std::cout << " = [";
                    for(int q=0;q<Vfaces.size();q++)
                    {
                        std::cout << "[" << Vfaces[q]->x << ", " << Vfaces[q]->y << ", " << Vfaces[q]->z << "], "<<std::endl;
                    }
                    std::cout << "]"<<std::endl;
                }
                
                
            fptot++;
            pid++;
        }
    
    //std::cout << "notany " << world_rank << " " << notany << std::endl;
        //std::cout << "ftoto " << world_rank << " " << ftot << " "  << ifnOUT->getNrow() << "  fptot " << cnttt << std::endl;

        
        //std::cout << local2globalVertMap.size() << " " << distriPrismVertmap->getNel() << " " << world_rank << " " << fptot << " " << int_face2node_prism.size()+shared_face2node_prism.size() <<  std::endl;
        
        DistributedParallelState* distriUniquePrismVertmap = new DistributedParallelState(local2globalVertMap.size(),comm);
        
        
        //std::cout << "cntty  " << " " << local2globalVertMap.size() << " " << distriUniquePrismVertmap->getNel() << " wr " << world_rank << std::endl;
        
        
        // End reducing the map.
        
        
        DistributedParallelState* distfptot  = new DistributedParallelState(fptot,comm);

//        int nTotInteriorFaces          = distftot->getNel();
//        int* TotIntFaces_offsets       = distftot->getOffsets();
        
        int nTotInteriorFaces_prism    = distfptot->getNel();
        int* TotIntFaces_offsets_prism = distfptot->getOffsets();
        
        std::map<int,std::vector<int> >::iterator bit;
        std::set<int> sorted_BCid;
        std::map<int,int> sorted_BCid_map;
        std::map<int,int> sorted_NBCid_map;
        int ii = 0;
        int Nbf;
        int nloc_bcs_p = pbcmap.size();
        std::set<int> bcids_tot;
        
        std::vector<int> Lbcs;
        
        for(bit=pbcmap.begin();bit!=pbcmap.end();bit++)
		{
            //std::cout << "bitp " << world_rank << " " << bit->first << " " << bit->second.size() << std::endl;

			if(bcids_tot.find(bit->first)==bcids_tot.end())
			{
				bcids_tot.insert(bit->first);
				Lbcs.push_back(bit->first);
			}
		}
        
        int nloc_bcs   = Lbcs.size();

        DistributedParallelState* distLocBCs = new DistributedParallelState(nloc_bcs,comm);
        
        int Nt_BCs                 = distLocBCs->getNel();
        int* BCs_offsets           = distLocBCs->getOffsets();
        int* BCs_nlocs             = distLocBCs->getNlocs();
        int* BCs_arr               = new int[Nt_BCs];
        
        MPI_Allgatherv(&Lbcs[0],
                       nloc_bcs,
                       MPI_INT,
                       BCs_arr,
                       BCs_nlocs,
                       BCs_offsets,
                       MPI_INT, comm);

        std::set<int> bcsToT;
        std::map<int,int> bcentry;
        for(int i=0;i<Nt_BCs;i++)
        {
            if(bcsToT.find(BCs_arr[i])==bcsToT.end())
            {
                bcsToT.insert(BCs_arr[i]);
            }
        }
        
        int* bcid = new int[bcsToT.size()];
        int* nlbc = new int[bcsToT.size()];
        std::set<int>::iterator entry;
        int cnt = 0;
        
        int lac;
        int q = 0;
        for(entry=bcsToT.begin();entry!=bcsToT.end();entry++)
        {
            int ee = *entry;
            int Nbft = 0;
            int Nbfp = 0;
            
            if(pbcmap.find(ee)!=pbcmap.end())
            {
            	Nbfp = pbcmap[ee].size();
            }

            bcid[q]  = ee;
            nlbc[q]  = Nbfp;
            
            //std::cout << " asff " << world_rank << " " << ee << " -> " << Nbft+Nbfp  << " " << Nbft << " " << Nbfp << std::endl;
            
            q++;
        }
        
        
            
        
        std::vector<Array<int>* > bcArrays;
        std::map<int,int> bcsizing;
        std::vector<int> bci_offsets;
        std::vector<int> bciTot_offsets;
        int nTotBCFaces = 0;
        int nTotBCFaces_offset = 0;
        
//
//        std::map<int,int> tagV2localV_prism                     = nnf->tagV2localV;
//        std::map<int,int> localV2tagV_prism                     = nnf->localV2tagV;
//
//
        int failbc = 0;
        int globalVid;
        for(int i=0;i<bcsToT.size();i++)
        {
            int bc_id = bcid[i];
            DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);
            int NelLoc_bci = nlbc[i];
            int NelTot_bci = distBCi->getNel();
            
            Array<int>* ifn_bc_i = new Array<int>(NelLoc_bci,8);
            int offsetbci        = distBCi->getOffsets()[world_rank];
            int fbc  = 0;
            int Nbft = 0;
			int Nbfp = 0;

			if(pbcmap.find(bc_id)!=pbcmap.end())
			{
				Nbfp = pbcmap[bc_id].size();
			}
            
            if(Nbfp!=0)
            {
                int sk = 0;

                for(int q=0;q<Nbfp;q++)
				{
                    
					int bcface = pbcmap[bc_id][q];
					int flag = -1;
					
					int nppf = bc_face2node_prism[bcface].size();
					
					ifn_bc_i->setVal(fbc,0,nppf);
                    std::vector<int> face_tmp(nppf);
                    std::vector<int> fce(nppf);
                    std::vector<Vert*> Vfaces;
                    Vert* VcF = new Vert;
                    VcF->x = 0.0;
                    VcF->y = 0.0;
                    VcF->z = 0.0;
                    
                    for(int g=0;g<nppf;g++)
                    {
                        int oldtag = bc_face2node_prism[bcface][g];
                        
                        if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end())
                        {
                            int lvp  = tag2locV_map[oldtag];
                            Vert* Vf = new Vert;
                            Vf->x = localVsPartition[lvp]->x;
                            Vf->y = localVsPartition[lvp]->y;
                            Vf->z = localVsPartition[lvp]->z;
                            VcF->x = VcF->x + Vf->x;
                            VcF->y = VcF->y + Vf->y;
                            VcF->z = VcF->z + Vf->z;
                            
                            Vfaces.push_back(Vf);
                            int globid = tag2glob_prism[oldtag];
                            fce[g]     = globid;
                        }
                        else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end())
                        {
                            int lvp  = tag2locV_map[oldtag];
                            Vert* Vf = new Vert;
                            Vf->x = localVsPartition[lvp]->x;
                            Vf->y = localVsPartition[lvp]->y;
                            Vf->z = localVsPartition[lvp]->z;
                            VcF->x = VcF->x + Vf->x;
                            VcF->y = VcF->y + Vf->y;
                            VcF->z = VcF->z + Vf->z;
                            
                            Vfaces.push_back(Vf);
                            int globid = SharedVertsNotOwned[oldtag];
                            fce[g]     = globid;
                        }
                    }
                    
                    VcF->x = VcF->x/nppf;
                    VcF->y = VcF->y/nppf;
                    VcF->z = VcF->z/nppf;
                    
                    int leftEl  = lhp[bcface];
                    int leftTag = gE2tagE[leftEl];
                    
                    Vert* Vijk = new Vert;
                    Vijk->x = 0.0;
                    Vijk->y = 0.0;
                    Vijk->z = 0.0;
                    // compute element center;
                    int nvp = prisms[leftTag].size();

                    for(int q=0;q<nvp;q++)
                    {
                        int tag  = prisms[leftTag][q];
                        int lvp  = tag2locV_map[tag];

                        Vijk->x = Vijk->x + localVsPartition[lvp]->x;
                        Vijk->y = Vijk->y + localVsPartition[lvp]->y;
                        Vijk->z = Vijk->z + localVsPartition[lvp]->z;
                    }

                    Vijk->x = Vijk->x/nvp;
                    Vijk->y = Vijk->y/nvp;
                    Vijk->z = Vijk->z/nvp;

                    double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
                    
                    if(orient0 < 0.0)
                    {
                        if(nppf == 3)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[2]);
                            ifn_bc_i->setVal(fbc,3,fce[1]);
                            ifn_bc_i->setVal(fbc,4,0);
                        }
                        if(nppf == 4)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[3]);
                            ifn_bc_i->setVal(fbc,3,fce[2]);
                            ifn_bc_i->setVal(fbc,4,fce[1]);
                        }
                        
                    }
                    else
                    {
                        if(nppf == 3)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[1]);
                            ifn_bc_i->setVal(fbc,3,fce[2]);
                            ifn_bc_i->setVal(fbc,4,0);
                        }
                        if(nppf == 4)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[1]);
                            ifn_bc_i->setVal(fbc,3,fce[2]);
                            ifn_bc_i->setVal(fbc,4,fce[3]);
                        }
                    }
                    
                    if(lhp[bcface] ==4 )
                    {
                        std::cout << "bc wrong " << lhp[bcface] << " " << orient0 << " Vijk " << Vijk->x << ", " << Vijk->y << ", " << Vijk->z << " VcF " << VcF->x << ", " << VcF->y << ", " << VcF->z << std::endl;
                        
                        std::cout << " = [";
                        for(int q=0;q<Vfaces.size();q++)
                        {
                            std::cout << "[" << Vfaces[q]->x << ", " << Vfaces[q]->y << ", " << Vfaces[q]->z << "], "<<std::endl;
                        }
                        std::cout << "]"<<std::endl;
                    }
					ifn_bc_i->setVal(fbc,5,0);
					ifn_bc_i->setVal(fbc,6,lhp[bcface]);
					ifn_bc_i->setVal(fbc,7,bc_id);
					
					fbc++;
				}
                
                
            }
//            
                         
            
            
            int nbt = Nbfp;
            
//            if(nbt!=(fbc-1))
//            {
//                std::cout << "Check lengths of bcfaces " << bc_id << std::endl;
//            }
//
            bcsizing[bc_id] = NelTot_bci;
            bci_offsets.push_back(offsetbci);
            bciTot_offsets.push_back(nTotBCFaces_offset);
            bcArrays.push_back(ifn_bc_i);
            
            //std::cout << " asff " << world_rank << " " << bc_id << " -> " << NelLoc_bci << " " << NelTot_bci << " " << Nbft << " " << Nbfp << " " << nTotBCFaces_offset << std::endl;

            //nLocBCFaces_offset = nLocBCFaces_offset + NelLoc_bci;
            nTotBCFaces_offset = nTotBCFaces_offset + NelTot_bci;
            nTotBCFaces        = nTotBCFaces + NelTot_bci;
            
            
        }
        
        //std::cout << world_rank << " nTotBCFaces " << " " << nTotBCFaces << " on -> " << world_rank  << std::endl;
        
        std::cout << "failbc  " << failbc << std::endl;
            
    
        int nPrismOUT = parmmg_iet_prisms->getNrow();
        //std::cout << "nPrismOUTnPrismOUTnPrismOUTnPrismOUT " << nPrismOUT << " " << world_rank << std::endl;
//        DistributedParallelState* distTetra      = new DistributedParallelState(nTetrahedraOUT,comm);
        DistributedParallelState* distPrism      = new DistributedParallelState(nPrismOUT,comm);
        //DistributedParallelState* distTetraVerts = new DistributedParallelState(xcn_parmmg->getNrow(),comm);
        DistributedParallelState* distPrismVerts = new DistributedParallelState(xcn_prisms_int->getNrow()+xcn_prisms_shared->getNrow(),comm);
        DistributedParallelState* distPrismIntVerts = new DistributedParallelState(xcn_prisms_int->getNrow(),comm);
        DistributedParallelState* distPrismShaVerts = new DistributedParallelState(xcn_prisms_shared->getNrow(),comm);
        
        int ToTElements_prism           = distPrism->getNel();
        int ToTElements_offset_prism    = distPrism->getOffsets()[world_rank];
        
        //int ToTElements                 = distTetra->getNel();
        //int ToTElements_offset          = distTetra->getOffsets()[world_rank];
        
        int nTotElements                = ToTElements_prism;
        int nTotFaces                   = nTotInteriorFaces_prism + nTotBCFaces;
        int nTotIntFaces                = nTotInteriorFaces_prism;
        
        int nTotPrismVerts_v2           = distPrismVerts->getNel();
        int nTotPrismIntVerts_v2        = distPrismIntVerts->getNel();
        int nTotPrismShaVerts_v2        = distPrismShaVerts->getNel();
    
        int TotPrismVerts_offset_int    = distPrismIntVerts->getOffsets()[world_rank];
        int TotPrismVerts_offset_sha    = distPrismShaVerts->getOffsets()[world_rank];

        int TotPrismVerts_offset        = distPrismVerts->getOffsets()[world_rank];
        
        int nTotVertsPrismTetra = nTotPrismVerts_v2;
        
        int nbo = bcArrays.size();
        //std::cout << "-- Constructing the zdefs array..."<<std::endl;
        Array<int>* adapt_zdefs = new Array<int>(3+nbo,7);
        // Collect node data (10) . Starting index-ending index Nodes
        adapt_zdefs->setVal(0,0,10);
        adapt_zdefs->setVal(0,1,-1);
        adapt_zdefs->setVal(0,2,1);
        adapt_zdefs->setVal(0,3,1);
        adapt_zdefs->setVal(0,4,nTotVertsPrismTetra);
        adapt_zdefs->setVal(0,5,us3d->zdefs->getVal(0,5));
        adapt_zdefs->setVal(0,6,us3d->zdefs->getVal(0,6));
        // Collect element data (12) . Starting index-ending index Element
        adapt_zdefs->setVal(1,0,12);
        adapt_zdefs->setVal(1,1,-1);
        adapt_zdefs->setVal(1,2,2);
        adapt_zdefs->setVal(1,3,1);
        adapt_zdefs->setVal(1,4,nTotElements);
        adapt_zdefs->setVal(1,5,us3d->zdefs->getVal(1,5));
        adapt_zdefs->setVal(1,6,2);
        // Collect internal face data (13) . Starting index-ending index internal face.
        adapt_zdefs->setVal(2,0,13);
        adapt_zdefs->setVal(2,1,-1);
        adapt_zdefs->setVal(2,2, 3);
        adapt_zdefs->setVal(2,3, 1);
        adapt_zdefs->setVal(2,4,nTotIntFaces);
        adapt_zdefs->setVal(2,5,us3d->zdefs->getVal(2,5));
        adapt_zdefs->setVal(2,6,2);
        
        int qq  = 1;
        int nb = 0;
        int face_start = nTotIntFaces+1;
        int face_end;
        std::map<int,int>::iterator itr;
        for(itr=bcsizing.begin();itr!=bcsizing.end();itr++)
        {
            int bnd_ref  = itr->first;
            int bnd_size = itr->second;
            face_end = face_start+bnd_size-1;
            adapt_zdefs->setVal(3+nb,0,13);
            adapt_zdefs->setVal(3+nb,1,-1);
            adapt_zdefs->setVal(3+nb,2,3+qq);
            adapt_zdefs->setVal(3+nb,3,face_start);
            adapt_zdefs->setVal(3+nb,4,face_end);
            adapt_zdefs->setVal(3+nb,5,bnd_ref);
            adapt_zdefs->setVal(3+nb,6,2);
            face_start = face_end+1;

            nb++;
            qq++;
        }
        
        
        
        if(world_rank == 0)
        {
            
            std::cout << "Total faces = " << nTotFaces << " " << nTotInteriorFaces_prism  << std::endl;
            
            PlotBoundaryData(us3d->znames,adapt_zdefs);
        }
        


        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        
        hid_t ret;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);
        hid_t file_id = H5Fcreate("par_grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

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
    value = nTotPrismVerts_v2;//+nTotTetraVerts_v2;//ToTVrts;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        //std::cout << " NN " << nTotVertsPrismTetra << std::endl;
        //====================================================================================
        // Add iet map to the grid.h5 file
        //====================================================================================
        dimsf[0] = nTotElements;
        dimsf[1] = parmmg_iet_prisms->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);

        dset_id = H5Dcreate(file_id, "iet",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = parmmg_iet_prisms->getNrow();
        countH5[1]  = parmmg_iet_prisms->getNcol();
        
        offsetH5[0] = ToTElements_offset_prism;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, parmmg_iet_prisms->data);
        delete parmmg_iet_prisms;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//        countH5[0]  = parmmg_iet->getNrow();
//        countH5[1]  = parmmg_iet->getNcol();
//
//        offsetH5[0] = ToTElements_prism+ToTElements_offset;
//        offsetH5[1] = 0;
//        memspace = H5Screate_simple(2, countH5, NULL);
//
//        //filespace = H5Dget_space(dset_id);
//        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
//        plist_id     = H5Pcreate(H5P_DATASET_XFER);
//        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//
//        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, parmmg_iet->data);
//        delete parmmg_iet;
        //====================================================================================
        // Add xcn map to the grid.h5 file
        //====================================================================================
        

        
        
        dimsf[0] = nTotPrismIntVerts_v2+nTotPrismShaVerts_v2;//+nTotTetraVerts_v2;
        dimsf[1] = xcn_prisms_int->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);
        
        dset_id = H5Dcreate(file_id, "xcn",
                            H5T_NATIVE_DOUBLE, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = xcn_prisms_int->getNrow();
        countH5[1]  = xcn_prisms_int->getNcol();
        
        offsetH5[0] = TotPrismVerts_offset_int;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_prisms_int->data);
        delete xcn_prisms_int;
    
//        dimsf[0] = nTotPrismShaVerts_v2;//+nTotTetraVerts_v2;
//        dimsf[1] = xcn_prisms_shared->getNcol();
//        filespace = H5Screate_simple(2, dimsf, NULL);
        
//        dset_id = H5Dcreate(file_id, "xcn",
//                            H5T_NATIVE_DOUBLE, filespace,
//                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = xcn_prisms_shared->getNrow();
        countH5[1]  = xcn_prisms_shared->getNcol();
        
            std::cout << "nTotPrismIntVerts_v2+TotPrismVerts_offset_sha " << nTotPrismIntVerts_v2+TotPrismVerts_offset_sha << " " << world_rank << std::endl;
        offsetH5[0] = nTotPrismIntVerts_v2+TotPrismVerts_offset_sha;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_prisms_shared->data);
    std::cout << "world " << nTotPrismIntVerts_v2 << " " << TotPrismVerts_offset_sha << " " << xcn_prisms_shared->getNrow() << std::endl;
        delete xcn_prisms_shared;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
//        countH5[0]  = xcn_parmmg->getNrow();
//        countH5[1]  = xcn_parmmg->getNcol();
//
//        offsetH5[0] = nTotPrismVerts_v2+ToTVrts_offset;
//        offsetH5[1] = 0;
//        memspace = H5Screate_simple(2, countH5, NULL);
//        //filespace = H5Dget_space(dset_id);
//        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
//        plist_id     = H5Pcreate(H5P_DATASET_XFER);
//        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//
//        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_parmmg->data);
//        delete xcn_parmmg;
        //===================================================================================
//        int nTotInteriorFaces          = distftot->getNel();
//        int* TotIntFaces_offsets       = distftot->getOffsets();
//
//        int nTotInteriorFaces_prism    = distfptot->getNel();
//        int* TotIntFaces_offsets_prism = distfptot->getOffsets();
        
        dimsf[0]  = nTotInteriorFaces_prism+nTotBCFaces;
        dimsf[1]  = ifnOUT_prism->getNcol();
        
        filespace = H5Screate_simple(2, dimsf, NULL);
        dset_id   = H5Dcreate(file_id, "ifn",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//
        countH5[0]  = ifnOUT_prism->getNrow();
        countH5[1]  = dimsf[1];
        
        offsetH5[0] = TotIntFaces_offsets_prism[world_rank];
        offsetH5[1] = 0;
        
        memspace      = H5Screate_simple(2, countH5, NULL);
        filespace     = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id      = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, ifnOUT_prism->data);
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
//        countH5[0]  = ifnOUT->getNrow();
//        countH5[1]  = dimsf[1];
//
//        offsetH5[0] = nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank];
//        offsetH5[1] = 0;
//        //std::cout << "nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank ] " << TotIntFaces_offsets_prism[world_rank] << " " << ifnOUT_prism->getNrow() << " " << nTotInteriorFaces_prism << " " << TotIntFaces_offsets[world_rank ] << " " << ifnOUT->getNrow() << world_rank << " -> " <<nTotInteriorFaces_prism << " + " << nTotInteriorFaces << " " << ftot << " " << fptot  << " " << world_rank << " " << nTotInteriorFaces_prism+nTotInteriorFaces << std::endl;
//
//        memspace     = H5Screate_simple(2, countH5, NULL);
//        filespace     = H5Dget_space(dset_id);
//        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
//        plist_id     = H5Pcreate(H5P_DATASET_XFER);
//        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
////
//        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
//                      plist_id, ifnOUT->data);
//
        //std::cout << "world " << nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank ]+ifnOUT->getNrow() << " --> " << nTotInteriorFaces_prism+nTotInteriorFaces << " +++ "<< world_rank << std::endl;
        //===================================================================================

        for(int i=0;i<bcsToT.size();i++)
        {
            
            
            int bc_id = bcid[i];
            DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);

            int NelLoc_bci = nlbc[i];
            int NelTot_bci = distBCi->getNel();

            Array<int>* ifn_bc_i = bcArrays[i];

            countH5[0]    = ifn_bc_i->getNrow();
            countH5[1]    = dimsf[1];
            
//            if(i==0)
//            {
//                std::cout << " first bnd offsets = " << world_rank << " -> " << i << " " << NelLoc_bci << " " << NelTot_bci<< " " << ifn_bc_i->getNrow() << " offeis " << bciTot_offsets[i]<<" "<<bci_offsets[i] << std::endl;
//
//            }
//            if(i==0)
//            {
//                std::cout << " zero bnd offsets = " << world_rank << " -> " << i << " " << NelLoc_bci << " " << NelTot_bci << " " << ifn_bc_i->getNrow() << " offeis " << bciTot_offsets[i] <<" "<< bci_offsets[i] << " for id " << bcid[i] << " Final " << nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i] << std::endl;
//            }
//            if(i==1)
//            {
//                std::cout << " first bnd offsets = " << world_rank << " -> " << i << " " << NelLoc_bci << " " << NelTot_bci << " " << ifn_bc_i->getNrow() << " offeis " << bciTot_offsets[i] <<" "<< bci_offsets[i] << " for id " << bcid[i] << " Final " << nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i]<< std::endl;
//            }
//            if(i==2)
//            {
//                std::cout << " secon bnd offsets = " << world_rank << " -> " << i << " " << NelLoc_bci << " " << NelTot_bci << " " << ifn_bc_i->getNrow() << " offeis " << bciTot_offsets[i] <<" "<< bci_offsets[i] << " for id " << bcid[i] << " Final " << nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i]<< std::endl;
//            }
//            if(i==3)
//            {
//                std::cout << " third bnd offsets = " << world_rank << " -> " << i << " " << NelLoc_bci << " " << NelTot_bci << " " << ifn_bc_i->getNrow() << " offeis " << bciTot_offsets[i] <<" "<< bci_offsets[i] << " for id " << bcid[i] << " Final " << nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i]<< std::endl;
//            }

            //std::cout << "nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i] " << world_rank << " " << nTotInteriorFaces_prism+nTotInteriorFaces << " " << nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i] << std::endl;
            offsetH5[0]   = nTotInteriorFaces_prism+bciTot_offsets[i]+bci_offsets[i];
            offsetH5[1]   = 0;
            memspace      = H5Screate_simple(2, countH5, NULL);
            filespace     = H5Dget_space(dset_id);

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);

            plist_id     = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    //
            status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                          plist_id, ifn_bc_i->data);

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
        dimsf[0] = adapt_zdefs->getNrow();
        dimsf[1] = adapt_zdefs->getNcol();
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

        status = H5Dwrite(dset_zdefs_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_zdefs->data);
        //====================================================================================
        
        dimsf_att = us3d->znames->getNrow();
        filespace = H5Screate_simple(1, &dimsf_att, NULL);
        type =  H5Tcopy (H5T_C_S1);
        ret  = H5Tset_size (type, 20);
        ret  = H5Tset_strpad(type, H5T_STR_SPACEPAD);
        hid_t dset_znames_id = H5Dcreate(group_zones_id, "znames", type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Sclose(filespace);

        hsize_t cnt2 = us3d->znames->getNrow();
        memspace  = H5Screate_simple(1, &cnt2, NULL);
        filespace = H5Dget_space(dset_znames_id);


        status = H5Dwrite(dset_znames_id, type, memspace, filespace, plist_id, us3d->znames->data);
    

        //===================================================================================
        //===================================================================================
        //===================================================================================
        
    /* */
       
    
    
     
    
    MPI_Finalize();
    
}

