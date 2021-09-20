#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"

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
    std::map<int,int> shellvert2ref_local           = tetra_distri->GetShellVert2RefMap_Local();
    std::map<int,int> tetF2hybF                     = tetra_distri->GetTetF2HybFMap();
    std::map<int,int> tetV2tagV						= tetra_distri->GetTet2TagVertMap();
    std::map<int,int> shellvertOriginalTag2ref_Glob = tetra_distri->GetShellVertTag2RefMap_Global();

    std::map<int,int> shellvertTag2ref;
    std::map<int,int> shellvertTag2ref2;

    std::map<int,int>::iterator itt;
    for(itt=shellvert2ref.begin();itt!=shellvert2ref.end();itt++)
    {
    	int vtet = itt->first;
    	int vhyb = tetV2tagV[vtet];
    	int ref  = itt->second;
    	shellvertTag2ref[vhyb] = ref;
    }
    
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
    

    std::map<int,int> tag2locV_map = P->getGlobalVert2LocalVert();
    std::vector<Vert*> localVsPartition = P->getLocalVerts();
    
    
    PrismaticLayer* prsmLyr = new PrismaticLayer(part_global, prisms,
                                                 ief_part_map->i_map, ifn_part_map->i_map,
                                                 ife_part_map->i_map, if_ref_part_map->i_map,
                                                 if_Nv_part_map->i_map, ushell, tag2locV_map,
                                                 localVsPartition, shellvertOriginalTag2ref_Glob, comm);
    
    nPrisms                                                 = prisms.size();
    DistributedParallelState* distPrismIN                   = new DistributedParallelState(nPrisms,comm);
    int nPrismsTot                                          = distPrismIN->getNel();
    
    std::map<int,int> tagE2gE                               = prsmLyr->getTag2GlobalElementMap();
    std::map<int,int> gE2tagE                               = prsmLyr->getGlobal2TagElementMap();
    std::map<int,int> rhp                                   = prsmLyr->getRightElementGlobalIDForFaces();
    std::map<int,int> lhp                                   = prsmLyr->getLeftElementGlobalIDForFaces();
    std::map<int,std::vector<int> > pbcmap                  = prsmLyr->getBoundaryCondition2FaceID();
    std::map<int,int> tag2element_shell                     = prsmLyr->getTag2Element4TetPrismInterface();
    std::map<int,std::vector<int> > shared_face2node_prism  = prsmLyr->getOwnedSharedFace2NodeMap();
    std::map<int,std::vector<int> > int_face2node_prism     = prsmLyr->getInternalFace2NodeMap();
    std::map<int,std::vector<int> > bc_face2node_prism      = prsmLyr->getBoundaryFace2NodeMap();
    std::map<int,int> tag2glob_prism                        = prsmLyr->getVertexTag2GlobalMap();
    Array<double>* xcn_prisms_int                           = prsmLyr->getInternalCoordinates();
    Array<double>* xcn_prisms_shared                        = prsmLyr->getSharedCoordinates();
    std::map<int,int> SharedVertsNotOwned                   = prsmLyr->getNotOwnedSharedVerticesMap();
    Array<int>* parmmg_iet_prisms                           = prsmLyr->getElementType();
    std::map<int,int> sharedVmap                            = prsmLyr->getSharedVertexMap();

    
    Array<int>* ifnOUT_prism = new Array<int>(int_face2node_prism.size()+shared_face2node_prism.size(),8);

    DistributedParallelState* distPrismIntVerts = new DistributedParallelState(xcn_prisms_int->getNrow(),comm);
    DistributedParallelState* distPrismShaVerts = new DistributedParallelState(xcn_prisms_shared->getNrow(),comm);

    std::map<int,std::vector<int> >::iterator itmm;
    int foundU = 0;
    
    std::map<std::set<int>, int > vertref2shell_prism;
    int nshell = 0;
    for(itmm=ushell.begin();itmm!=ushell.end();itmm++)
    {
        int fhyb      = itmm->first;
        int Etettag   = itmm->second[0];
        int Eprismtag = itmm->second[1];
        int EprismNew = tag2element_shell[Eprismtag];

        if(shellface2vertref.find(fhyb)!=shellface2vertref.end())
        {
            vertref2shell_prism[shellface2vertref[fhyb]] = EprismNew;
            foundU++;
        }
        nshell++;
    }
    
    
    //===============================================================================================
    //===============================================================================================
    //===============================================================================================
    
    
    PMMG_pParMesh   parmesh;
    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);

    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                              nQuadrilaterals,nEdges) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    
    //PMMG_Set_metSize(PMMG_pParMesh parmesh,int typEntity,int np,int typSol)
    if ( PMMG_Set_metSize(parmesh,MMG5_Vertex,nVertices,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);

    
    int vrefmax = -1;
    for ( k=0; k<nVertices; ++k )
    {
        double vx = locVs[k]->x;
        double vy = locVs[k]->y;
        double vz = locVs[k]->z;
        
        int vert = locV2globV[k];

        int vref = 1;
        


        //std::cout << world_rank << " vref = " << vref << std::endl;
        if ( PMMG_Set_vertex(parmesh,vx,vy,vz, vref, k+1) != 1 )
        {
        MPI_Finalize();
        exit(EXIT_FAILURE);
        }

        double* tensor = new double[6];


        tensor[0] = metric[vert]->getVal(0,0);
        tensor[1] = metric[vert]->getVal(1,0);
        tensor[2] = metric[vert]->getVal(2,0);
        tensor[3] = metric[vert]->getVal(3,0);
        tensor[4] = metric[vert]->getVal(4,0);
        tensor[5] = metric[vert]->getVal(5,0);

        if(PMMG_Set_tensorMet(parmesh,tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],k+1)!=1)
        {
         MPI_Finalize();
         exit(EXIT_FAILURE);
        }
    }
    
   
    int v0,v1,v2,v3;
    int v0l,v1l,v2l,v3l;
    int teller = 0;
    int refer  = 0;
    int cref36 = 0;
    double* c0 = new double[3];
    std::map<int,int> vidConsidered;
    std::map<std::set<int>, int> faceUvids;
    int iref;
    int c13 = 0;
    int suc = 0;
    int buggi = 0;
    int shelfound = 0;
    int shelfound2 = 0;
    for ( k=0; k<nTriangles; ++k )
    {
        int faceID  = faces4parmmg[k];
        int facetag = tetF2hybF[faceID];
    
        v0      = face2node[faceID][0];
        v1      = face2node[faceID][1];
        v2      = face2node[faceID][2];
        
        v0l     = globV2locV[v0];
        v1l     = globV2locV[v1];
        v2l     = globV2locV[v2];
        
        refer   = face2ref[faceID];
        
        if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,refer,k+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        
        if(refer == 13)
        {
            if(shell_tet2hybF.find(faceID)!=shell_tet2hybF.end())
            {
                if(shell_tet2hybF[faceID]!=facetag)
                {
                    std::cout << "Error:: Tet2Hyb maps are incorrect." << std::endl;
                }
                else
                {
                    std::set<int>::iterator its;
                    std::set<int> sf    = shellface2vertref[facetag];
                    
                    std::set<int> sft_t;
                    for(int u=0;u<3;u++)
                    {
                        int vid = face2node[faceID][u];
                        int vidg = globV2locV[vid];
                        double* c0 = new double[3];
                        c0[0] = locVs[vidg]->x;
                        c0[1] = locVs[vidg]->y;
                        c0[2] = locVs[vidg]->z;
                        int vertref = shellvert2ref[vid];
                        
                        PMMG_Set_vertex( parmesh, c0[0], c0[1], c0[2], vertref, vidg+1 );

                        sft_t.insert(vertref);
                        
                    }

//                    if(vertref2shell_prism.find(sf)!=vertref2shell_prism.end())
//                    {
//                        shelfound++;
//                    }
                    
                    sft_t.clear();
                }
                suc++;
            }
            PMMG_Set_requiredTriangle( parmesh, k+1 );
        }
    }
    
    
    
    int API_mode = 0;
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);
    
    for(int icomm=0; icomm<ncomm; icomm++ )
    {
      // Set nb. of entities on interface and rank of the outward proc
     
      ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                             color_face[icomm],
                                             ntifc[icomm]);

    //Set local and global index for each entity on the interface
      ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                               ifc_tria_loc[icomm],
                                               ifc_tria_glob[icomm], 1);
    }
    
    
    std::map<int,std::vector<int> >::iterator ittet;
    k = 0;
    for ( int t = 0;t < nTetrahedra; t++  )
    {
        v0 = ien_part_tetra->getVal(t,0);
        v1 = ien_part_tetra->getVal(t,1);
        v2 = ien_part_tetra->getVal(t,2);
        v3 = ien_part_tetra->getVal(t,3);
        
        v0l = globV2locV[v0];
        v1l = globV2locV[v1];
        v2l = globV2locV[v2];
        v3l = globV2locV[v3];
        
        double* P = new double[4*3];
        
        P[0*3+0]=locVs[v0l]->x;   P[0*3+1]=locVs[v0l]->y;    P[0*3+2]=locVs[v0l]->z;
        P[1*3+0]=locVs[v1l]->x;   P[1*3+1]=locVs[v1l]->y;    P[1*3+2]=locVs[v1l]->z;
        P[2*3+0]=locVs[v2l]->x;   P[2*3+1]=locVs[v2l]->y;    P[2*3+2]=locVs[v2l]->z;
        P[3*3+0]=locVs[v3l]->x;   P[3*3+1]=locVs[v3l]->y;    P[3*3+2]=locVs[v3l]->z;

        double Vtet = GetQualityTetrahedra(P);
        
        if(Vtet<0.0)
        {
            std::cout << " negative volume in Element " << t << " on rank " << world_rank  <<std::endl;
        }
        if ( PMMG_Set_tetrahedron(parmesh,v0l+1,v1l+1,v2l+1,v3l+1,1.0,t+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }
    
    int niter = 3;
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, niter ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hgrad, 2.0) != 1 )
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
    for( k = 1; k <= nVerticesOUT; k++ )
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
    
    int *required = (int*)calloc(MAX4(nVerticesOUT,nTetrahedraOUT,nTrianglesOUT,nEdgesOUT),sizeof(int));
    //int *refOUT = (int*)calloc(MAX4(nVerticesOUT,nTetrahedraOUT,nTrianglesOUT,nEdgesOUT),sizeof(int));
    int *corner = (int*)calloc(nVerticesOUT,sizeof(int));
    int pos;
    int *refOUT = new int[nVerticesOUT];

    std::vector<std::vector<int> > outT;
    
    double *vertOUT = (double*)calloc((nVerticesOUT)*3,sizeof(double));
    
    for ( k=0; k<nVerticesOUT; k++ ) {
          pos = 3*k;
          if ( PMMG_Get_vertex(parmesh,&(vertOUT[pos]),&(vertOUT[pos+1]),&(vertOUT[pos+2]),
                               &(refOUT[k]),&(corner[k]),&(required[k])) != 1 ) {
            fprintf(inm,"Unable to get mesh vertex %d \n",k);
            ier = PMMG_STRONGFAILURE;
          }
    }
    
    std::vector<int> rf;
    int kn = 0;
    int foundRref = 0;
    std::map<int,int> frf;
    for ( k=0; k<nVerticesOUT; k++ )
    {
        if(refOUT[k]!=0 && refOUT[k]!=1)
        {
            if(frf.find(refOUT[k]) == frf.end())
            {
                frf[refOUT[k]] = k;
                rf.push_back(refOUT[k]);
                foundRref++;
            }
        }
    }
    
    std::map<int,int> outShellVerts = AllGatherMap(frf,comm);
    
    double minrf = *min_element(rf.begin(), rf.end());
    double maxrf = *max_element(rf.begin(), rf.end());

    //std::cout << "found Ref " << foundRref << " " << minrf << " " << maxrf << " " << world_rank << " " << outShellVerts.size() << std::endl;
    
    if(niter==0)
    {
        
        std::cout << "Check the input and outputted shard faces." << std::endl;
        
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
        
        // Get input triangle nodes
        int** faceNodes2 = (int **) malloc(ncomm*sizeof(int *));
        for( int icomm = 0; icomm < ncomm; icomm++ ) {
          faceNodes2[icomm] = (int *) malloc(3*ntifc[icomm]*sizeof(int));
          for( i = 0; i < ntifc[icomm]; i++ ) {
              
            int faceID = ifc_tria_loc[icomm][i]-1;
            int faceID2 = locShF2globShF[faceID];

            v0 = face2node[faceID2][0];
            v1 = face2node[faceID2][1];
            v2 = face2node[faceID2][2];
    
            v0l = globV2locV[v0];
            v1l = globV2locV[v1];
            v2l = globV2locV[v2];

            //pos = ifc_tria_loc[icomm][i];
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
              for( i = 0; i < nitem_face_comm[icomm]; i++ )
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
                                              color_face_out,faceNodes_out) ) {
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

                for( i = 0; i < nitem_face_comm[icomm]; i++ )
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
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft = out_tria_loc[icomm][i];

                    std::set<int> faceSh;

                    for(int k=0;k<3;k++)
                    {
                        int vt3 = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];
                        faceSh.insert(vt3);
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
    else
    {
        //std::cout << "Outputting the new triangles..." << std::endl;
        //===================================================================
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
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
        
//
//        PMMG_pGrp      listgrp,grpI;
//        listgrp         = parmesh->listgrp;
//        grpI            = &listgrp[0];
//        int            poi_id_int,poi_id_glo,imsh,k;
//
//        std::map<int,int> locSh2globSh;
//        std::map<int,int> globSh2locSh;
//
//        std::map<int,int> globSh2Rank;
//        for ( k = 0; k < grpI->nitem_int_node_comm; ++k )
//        {
//            poi_id_int = grpI->node2int_node_comm_index1[k];
//            poi_id_glo = grpI->node2int_node_comm_index2[k];
//            locSh2globSh[poi_id_int]=poi_id_glo;
//            globSh2locSh[poi_id_int]=poi_id_int;
//            globSh2Rank[poi_id_glo]=world_rank;
//        }
//
    
        // Check matching of input interface nodes with the set ones
        int *ref2       = (int*)calloc(nTrianglesOUT,sizeof(int));
        int *required2  = (int*)calloc(nTrianglesOUT,sizeof(int));
        int *triaNodes2 = (int*)calloc(3*nTrianglesOUT,sizeof(int));
        
//        if ( PMMG_Get_triangles(parmesh,triaNodes2,ref2,required2) != 1 ) {
//          fprintf(stderr,"FAILED:: Unable to get mesh triangles\n");
//          ier = PMMG_STRONGFAILURE;
//        }
        
        int nreq2   = 0;
        int pos     = 0;
        int cnt_bnd = 0;
        int cnt_th  = 0;
        int cnt_int = 0;
        std::map<std::set<int>, int> PMMG_Face2Ref;
        std::set<std::set<int> > TotalFaces;
        std::set<std::set<int> > InteriorFaces;
        std::map<std::set<int>, int> PMMG_Shell2Prism;
        std::map<int,std::vector<int> > PMMG_ShellFace;

        int c36 = 0;
        int gv0,gv1,gv2;
        int prismFound = 0;
        std::map<std::set<int>, int> PMMG_ShellFaces;
        int faceCovered = 0;
        std::vector<int> vref1;
        for ( k=0; k<nTrianglesOUT; k++ )
        {
            pos = 3*k;
            if ( PMMG_Get_triangle(parmesh,&(triaNodes2[pos]),
                                         &(triaNodes2[pos+1]),
                                         &(triaNodes2[pos+2]),
                                         &(ref2[k]),
                                         &(required2[k])) != 1 )
            {
            fprintf(inm,"Unable to get mesh triangle %d \n",k);
            ier = PMMG_STRONGFAILURE;
            }

            std::set<int> faceSh;
            gv0 = loc2globVid[triaNodes2[pos]];
            gv1 = loc2globVid[triaNodes2[pos+1]];
            gv2 = loc2globVid[triaNodes2[pos+2]];
        
            faceSh.insert(gv0);
            faceSh.insert(gv1);
            faceSh.insert(gv2);
            
            TotalFaces.insert(faceSh);
            
            if(ref2[k]==0)
            {
                InteriorFaces.insert(faceSh);
            }
            if(ref2[k]!=0 && ref2[k]!=13)
            {
                PMMG_Face2Ref[faceSh] = ref2[k];
            }
            if(ref2[k]==13)
            {
                PMMG_ShellFaces[faceSh] = ref2[k];
                
                std::set<int> test_set;
                test_set.insert(refOUT[triaNodes2[pos]-1]);
                test_set.insert(refOUT[triaNodes2[pos+1]-1]);
                test_set.insert(refOUT[triaNodes2[pos+2]-1]);
//                std::vector<int> test_set_order;
//                test_set_order.push_back(refOUT[triaNodes2[pos]-1]);
//                test_set_order.push_back(refOUT[triaNodes2[pos+1]-1]);
//                test_set_order.push_back(refOUT[triaNodes2[pos+2]-1]);
                std::vector<int> shFace(3);
                shFace[0] = triaNodes2[pos];
                shFace[1] = triaNodes2[pos+1];
                shFace[2] = triaNodes2[pos+2];
                
                if(vertref2shell_prism.find(test_set)!=vertref2shell_prism.end())
                {
                    PMMG_Shell2Prism[faceSh] = vertref2shell_prism[test_set];
                    //std::cout << "vertref2shell_prism[test_set] " << vertref2shell_prism[test_set] << std::endl;
                    PMMG_ShellFace[vertref2shell_prism[test_set]] = shFace;
                    prismFound++;
                }
//
                test_set.clear();
                
                faceCovered++;
                // Need to catch the references here //
            }
         
            faceSh.clear();
    
            if ( required2 && required2[k] )  nreq2++;
        }
        
        //std::cout << world_rank << " PMMG_ShellFaces size  " << PMMG_ShellFaces.size() << " prismFound-> " << prismFound << " " << faceCovered << " " << PMMG_Shell2Prism.size() <<  std::endl;
        
        int itt2            = 0;
        int nTshared_owned  = 0;

        int vt,ft,gvt,gft;
        
        std::map<std::set<int>, int> PMMG_SharedFaces;
        std::map<std::set<int>, int> PMMG_SharedFacesOwned;
        std::set<int> PMMG_SharedVertsOwned;
        int nPartFace = 0;
        std::map<int,int> rank2icomm;
        std::map<int,std::vector<int> > LocateOppositeFace;
        std::map<int,int> vertonrank;
        int locID_NotShVrt = 0;
        std::vector<int> PMMG_SharedVertsOwned_vec;
        std::vector<int> PMMG_SharedVertsOwnedRank_vec;
        std::set<int> PMMG_SharedVertices;
        std::map<int,std::vector<int> > Color_SharedOwned;

        for(int icomm = 0; icomm < next_face_comm; icomm++ )
        {
            nPartFace                           = nPartFace + nitem_face_comm[icomm];
            rank2icomm[color_face_out[icomm]] = icomm;
                    
            if(world_rank < color_face_out[icomm])
            {
                nTshared_owned = nTshared_owned + nitem_face_comm[icomm];
                
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft       = out_tria_loc[icomm][i];
                    int face_ref = ref2[ft-1];
                    //Color_SharedOwned[icomm].push_back(i);
                    std::set<int> faceSh;
                    for(int k=0;k<3;k++)
                    {
                        int vt  = triaNodes2[3*(ft-1)+k];
                        int gvt = loc2globVid[vt];
                        faceSh.insert(gvt);
                        
                        if(PMMG_SharedVertices.find(gvt)==PMMG_SharedVertices.end())
                        {
                            PMMG_SharedVertices.insert(gvt);
                        }
                        
                        if(PMMG_SharedVertsOwned.find(gvt)==PMMG_SharedVertsOwned.end())
                        {
                            PMMG_SharedVertsOwned.insert(gvt);
                            PMMG_SharedVertsOwned_vec.push_back(gvt);
                            PMMG_SharedVertsOwnedRank_vec.push_back(world_rank);
                        }
                    }
                    
                    
                    Color_SharedOwned[color_face_out[icomm]].push_back(i);
                    PMMG_SharedFacesOwned[faceSh]=ft;
                    PMMG_SharedFaces[faceSh]=ft;
                    faceSh.clear();
                }
            }
            else
            {
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft       = out_tria_loc[icomm][i];
                    
                    std::set<int> faceSh;
                    for(int k=0;k<3;k++)
                    {
                        int vt  = triaNodes2[3*(ft-1)+k];
                        int gvt = loc2globVid[vt];
                        faceSh.insert(gvt);
                        if(PMMG_SharedVertices.find(gvt)==PMMG_SharedVertices.end())
                        {
                            PMMG_SharedVertices.insert(gvt);
                        }
                    }
                    if(PMMG_SharedFaces.find(faceSh)==PMMG_SharedFaces.end())
                    {
                        PMMG_SharedFaces[faceSh]=ft;
                    }
                    
                    faceSh.clear();
                }
            }
        }
       
        int nLocallyOwnedVerts                      = PMMG_SharedVertsOwned.size();
        DistributedParallelState* locallyOwnedVerts = new DistributedParallelState(nLocallyOwnedVerts,comm);
        int* OwnedVertsDistri                       = new int[locallyOwnedVerts->getNel()];
        int* OwnedVertsDistriRank                   = new int[locallyOwnedVerts->getNel()];
        int* PMMG_SharedVertsOwned_arr              = new int[nLocallyOwnedVerts];
        int* PMMG_SharedVertsOwnedRank_arr          = new int[nLocallyOwnedVerts];
        
        for(int i=0;i<nLocallyOwnedVerts;i++)
        {
            PMMG_SharedVertsOwned_arr[i]            = PMMG_SharedVertsOwned_vec[i];
            PMMG_SharedVertsOwnedRank_arr[i]        = PMMG_SharedVertsOwnedRank_vec[i];
        }
        
        MPI_Allgatherv(&PMMG_SharedVertsOwned_arr[0],
                       nLocallyOwnedVerts,
                       MPI_INT,
                       OwnedVertsDistri,
                       locallyOwnedVerts->getNlocs(),
                       locallyOwnedVerts->getOffsets(),
                       MPI_INT, comm);
        
        MPI_Allgatherv(&PMMG_SharedVertsOwnedRank_arr[0],
                       nLocallyOwnedVerts,
                       MPI_INT,
                       OwnedVertsDistriRank,
                       locallyOwnedVerts->getNlocs(),
                       locallyOwnedVerts->getOffsets(),
                       MPI_INT, comm);
        
        std::map<int,int> ActualOwnedVertDistr;
        std::map<int,std::vector<int> > ActualOwnedVertDistr_map;
        int ownedID;
        
        for(int i=0;i<locallyOwnedVerts->getNel();i++)
        {
            int gvd = OwnedVertsDistri[i];
            
            if(ActualOwnedVertDistr.find(gvd)==ActualOwnedVertDistr.end())
            {
                ActualOwnedVertDistr_map[gvd].push_back(OwnedVertsDistriRank[i]);

                ActualOwnedVertDistr[gvd] = OwnedVertsDistriRank[i];
                ownedID = OwnedVertsDistriRank[i];
            }
            else
            {
                if(OwnedVertsDistriRank[i]<ActualOwnedVertDistr[gvd])
                {
                    ActualOwnedVertDistr[gvd] = OwnedVertsDistriRank[i];
                    ownedID = OwnedVertsDistriRank[i];
                }
                else
                {
                    ActualOwnedVertDistr[gvd] = OwnedVertsDistriRank[i];
                    ownedID = ActualOwnedVertDistr[gvd];
                }
            }
        }
        
        std::map<int,std::vector<int> > ActualOwnedVertDistr_map_update;

        std::map<int,std::vector<int> >::iterator itm;
        std::map<int,int>::iterator iitm;
        int tel = 0;
        int* tells = new int[world_size];
        for(int u=0;u<world_size;u++)
        {
            tells[u] = 0;
        }
        
        std::map<int,int> LocationSharedVert;
        for(itm=ActualOwnedVertDistr_map.begin();itm!=ActualOwnedVertDistr_map.end();itm++)
        {
            int gv = itm->first;
            int ra = *min_element(itm->second.begin(), itm->second.end());
            tells[ra] = tells[ra]+1;
            ActualOwnedVertDistr_map_update[ra].push_back(gv);
            LocationSharedVert[gv]=ra;
            tel++;
        }

        
        
        
        int* tetraOUT = new int[nTetrahedraOUT*4];
        std::vector<std::vector<int> > tetrasOUT;
        std::vector<int> NonSharedVrts_vec;

        std::map<int,std::vector<double> > coords_int;
        std::map<int,int> lhshown;
        std::map<int,int> lhbnd;
        std::map<int,int> lhshell;
        std::map<int,int> rhshell;
        std::map<int,int> lhsh;
        std::map<int,int> lh;
        std::map<int,int> rh;
        
        DistributedParallelState* distTetraOut = new DistributedParallelState(nTetrahedraOUT,comm);
        int* TetraOUT_offsets = distTetraOut->getOffsets();
        
        std::map<int,std::vector<int> > ienOUT;
        
        std::map<int,std::vector<int> > fmBnd;
        std::map<int,std::vector<int> > fmInt;
        std::map<int,std::vector<int> > fmSha;
        std::map<int,std::vector<int> > fmShell;
        std::map<int,std::vector<int> > fm;
        
        std::map<std::set<int>, int> facemap;
        std::map<int,std::set<int> > facemap_inv;
        std::map<int,std::vector<int> > bcmap;
        
        std::map<int,std::vector<int> > pfmInt;
        std::map<int,std::vector<int> > pfmSha;
        std::map<int,std::vector<int> > pfm;
        
        std::map<std::set<int>, int> pfacemap;
        std::map<int,std::set<int> > pfacemap_inv;

        Array<int>* parmmg_iet = new Array<int>(nTetrahedraOUT,1);
        int q;
        
        std::map<int,int> gvid2shid;
        int fid  = 0;
        int lshf = 0;
        int tetra_faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
        
        int NoNShFaces = 0;
        
        std::map<int,int> locShF2globShF;
        std::map<int,int> locBF2globBF;
        std::map<int,int> NonSharedVrts;
        int ngv = 0;
        
        std::vector<int> Elvrts(4);
        int fv0,fv1,fv2;
        Vec3D* v0 = new Vec3D;
        Vec3D* v1 = new Vec3D;
        Vec3D* r0 = new Vec3D;

        Vert* Vface = new Vert;
        int negit = 0;
        
        
        int curElID = TetraOUT_offsets[world_rank]+1+nPrismsTot;

        int* refTET = new int[nTetrahedraOUT];
        
        std::map<int,int> tag2shelltag;
        std::map<int,int> shelltag2tag;
        int hellofound = 0;
        
        std::cout << "Start ID for tets on rank = " << curElID << std::endl;
        
        for ( k=0; k<nTetrahedraOUT; k++ )
        {
            parmmg_iet->setVal(k,0,2);
            pos = 4*k;
            if ( PMMG_Get_tetrahedron(parmesh,
                                        &(tetraOUT[pos  ]),&(tetraOUT[pos+1]),
                                        &(tetraOUT[pos+2]),&(tetraOUT[pos+3]),
                                        &(refTET[k]),&(required[k])) != 1 )
            {
                fprintf(inm,"Unable to get mesh tetra %d \n",k);
                ier = PMMG_STRONGFAILURE;
            }
            
            Elvrts[0] = tetraOUT[pos];
            Elvrts[1] = tetraOUT[pos+1];
            Elvrts[2] = tetraOUT[pos+2];
            Elvrts[3] = tetraOUT[pos+3];
            
            double* P = new double[4*3];
            
            P[0*3+0]=vertOUT[(Elvrts[0]-1)*3];
            P[0*3+1]=vertOUT[(Elvrts[0]-1)*3+1];
            P[0*3+2]=vertOUT[(Elvrts[0]-1)*3+2];
            
            P[1*3+0]=vertOUT[(Elvrts[1]-1)*3];
            P[1*3+1]=vertOUT[(Elvrts[1]-1)*3+1];
            P[1*3+2]=vertOUT[(Elvrts[1]-1)*3+2];
            
            P[2*3+0]=vertOUT[(Elvrts[2]-1)*3];
            P[2*3+1]=vertOUT[(Elvrts[2]-1)*3+1];
            P[2*3+2]=vertOUT[(Elvrts[2]-1)*3+2];
            
            P[3*3+0]=vertOUT[(Elvrts[3]-1)*3];
            P[3*3+1]=vertOUT[(Elvrts[3]-1)*3+1];
            P[3*3+2]=vertOUT[(Elvrts[3]-1)*3+2];
            
            ienOUT[curElID] = Elvrts;
            tetrasOUT.push_back(Elvrts);
            double Vtet = GetQualityTetrahedra(P);
            Vert* vCenter = ComputeCentroidCoord(P, 4);
            
            if(Vtet<0.0)
            {
                std::cout << "Error " << Vtet << std::endl;
            }
            
            
            for(int u=0;u<4;u++)
            {
                fv0 = loc2globVid[Elvrts[tetra_faces[u][0]]];
                fv1 = loc2globVid[Elvrts[tetra_faces[u][1]]];
                fv2 = loc2globVid[Elvrts[tetra_faces[u][2]]];
                
                std::set<int> Face;
                Face.insert(fv0);
                Face.insert(fv1);
                Face.insert(fv2);
                
                if(facemap.find(Face)==facemap.end())
                {
                    facemap[Face] = fid;
                    std::vector<int> fce(3);
                    fce[0]  = fv0;
                    fce[1]  = fv1;
                    fce[2]  = fv2;
                    
                    double v0x = vertOUT[(Elvrts[tetra_faces[u][0]]-1)*3];
                    double v0y = vertOUT[(Elvrts[tetra_faces[u][0]]-1)*3+1];
                    double v0z = vertOUT[(Elvrts[tetra_faces[u][0]]-1)*3+2];
                    
                    double v1x = vertOUT[(Elvrts[tetra_faces[u][1]]-1)*3];
                    double v1y = vertOUT[(Elvrts[tetra_faces[u][1]]-1)*3+1];
                    double v1z = vertOUT[(Elvrts[tetra_faces[u][1]]-1)*3+2];
                    
                    double v2x = vertOUT[(Elvrts[tetra_faces[u][2]]-1)*3];
                    double v2y = vertOUT[(Elvrts[tetra_faces[u][2]]-1)*3+1];
                    double v2z = vertOUT[(Elvrts[tetra_faces[u][2]]-1)*3+2];
                    
                    Vface->x = (v0x+v1x+v2x)/3.0;
                    Vface->y = (v0y+v1y+v2y)/3.0;
                    Vface->z = (v0z+v1z+v2z)/3.0;
                    
                    r0->c0 = (Vface->x-vCenter->x);///Lr;
                    r0->c1 = (Vface->y-vCenter->y);///Lr;
                    r0->c2 = (Vface->z-vCenter->z);///Lr;
                    
                    v0->c0 = v1x-v0x;
                    v0->c1 = v1y-v0y;
                    v0->c2 = v1z-v0z;

                    v1->c0 = v2x-v0x;
                    v1->c1 = v2y-v0y;
                    v1->c2 = v2z-v0z;
                    
                    Vec3D* n0        = ComputeSurfaceNormal(v0,v1);
                    double orient0   = DotVec3D(r0,n0);
                    if(orient0<0.0)
                    {
                        negit++;
                    }
                    
                    fm[fid] = fce;
                    if(PMMG_SharedFaces.find(Face) != PMMG_SharedFaces.end())
                    {
                        lshf                     = PMMG_SharedFaces[Face];
                        locShF2globShF[lshf]     = fid;
                        lhsh[fid]                = curElID;
                    }
                    
                    
                    if(PMMG_SharedFaces.find(Face)   == PMMG_SharedFaces.end()
                       && PMMG_Face2Ref.find(Face)   == PMMG_Face2Ref.end()
                       && PMMG_ShellFaces.find(Face) == PMMG_ShellFaces.end())
                    {
                        fmInt[fid]  = fce;
                        lh[fid]     = curElID;
                        
                        for(int s=0;s<3;s++)
                        {
                            int gvm2 = fce[s];
                            if(NonSharedVrts.find(gvm2)==NonSharedVrts.end() &&
                                    PMMG_SharedVertices.find(gvm2)==PMMG_SharedVertices.end())
                            {
                                NonSharedVrts[gvm2] = ngv;
                                ngv++;
                            }
                        }
                    }
                    
                    
                    
                    
                    if(PMMG_SharedFacesOwned.find(Face) != PMMG_SharedFacesOwned.end())
                    {
                        fmSha[fid]               = fce;
                        lhshown[fid]             = curElID;
                    }
                    
                    
                    
                    if(PMMG_ShellFaces.find(Face) != PMMG_ShellFaces.end())
                    {
                        fmShell[fid]        = fce;
                        lhshell[fid]        = curElID;
                        
                        if(PMMG_Shell2Prism.find(Face)!=PMMG_Shell2Prism.end())
                        {
                            rhshell[fid] = PMMG_Shell2Prism[Face];
                            //std::cout << "Do we find it ? " << PMMG_Shell2Prism[Face] << std::endl;

                            // std::cout << "rhshell[fid] "  << rhshell[fid] << std::endl;
                            
                            for(int y=0;y<3;y++)
                            {
                                int shelltag = refOUT[Elvrts[tetra_faces[u][y]]-1];
                                int vtag     = fce[y];
                                
                                if(tag2shelltag.find(vtag)==tag2shelltag.end())
                                {
                                    tag2shelltag[vtag]=shelltag;
                                    shelltag2tag[shelltag]=vtag;
                                }
                            }
                            
                            hellofound++;
                            
                        }
                        else
                        {
                            
                        }
                    }
                    
                    
                    if(PMMG_Face2Ref.find(Face) != PMMG_Face2Ref.end())
                    {
                        int FaceRef         = PMMG_Face2Ref[Face];
                        fmBnd[fid]          = fce;
                        lhbnd[fid]          = curElID;
                        bcmap[FaceRef].push_back(fid);
                    }
                }
                else
                {
                    int fid_n     = facemap[Face];
                    rh[fid_n]     = curElID;
                }
                
                    
                fid++;
                
                Face.clear();
            }
            
            curElID++;

            delete[] P;
        }
        
        
        std::map<int,int> tag2shelltag_glob = AllGatherMap(tag2shelltag,comm);
        std::map<int,int> shelltag2tag_glob = AllGatherMap(shelltag2tag,comm);
        //std::cout << negit << " normals are pointing inwards. " << tag2shelltag.size() << " " << tag2shelltag_glob.size() << " vs " << shelltag2tag_glob.size() << " " << hellofound << " on rnk " << world_rank << std::endl;
        int nLocIntVrts         = NonSharedVrts.size();
        int nLocShVrts          = ActualOwnedVertDistr_map_update[world_rank].size();
        int nLocTotVrts         = nLocIntVrts+nLocShVrts;
        DistributedParallelState* distnLocTotVrts = new DistributedParallelState(nLocTotVrts,comm);

        int vert = 0;//distnLocTotVrts->getOffsets()[world_rank];
        int ToTVrts =  distnLocTotVrts->getNel();
        int* ToTVrts_offsets = distnLocTotVrts->getOffsets();
        int ToTVrts_offset = ToTVrts_offsets[world_rank];
        Array<double>* xcn_parmmg = new Array<double>(nLocTotVrts,3);
        std::map<int,int> tag2glob;
        std::map<int,int> glob2tag;
        
        
        
        int nPrismVerts_tmp = distPrismIntVerts->getNel()+distPrismShaVerts->getNel();
        

        for(iitm=NonSharedVrts.begin();iitm!=NonSharedVrts.end();iitm++)
        {
            int tag = iitm->first;
            int lvert = glob2locVid[tag];

            double xc = vertOUT[(lvert-1)*3+0];
            double yc = vertOUT[(lvert-1)*3+1];
            double zc = vertOUT[(lvert-1)*3+2];

            
            xcn_parmmg->setVal(vert,0,xc);
            xcn_parmmg->setVal(vert,1,yc);
            xcn_parmmg->setVal(vert,2,zc);
            
            tag2glob[tag] = vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp;
            glob2tag[vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp] = tag;
            
            vert++;
        }
        
        
        for(int i=0;i<nLocShVrts;i++)
        {
            int tag = ActualOwnedVertDistr_map_update[world_rank][i];
            int lvert = glob2locVid[tag];

            double xc = vertOUT[(lvert-1)*3+0];
            double yc = vertOUT[(lvert-1)*3+1];
            double zc = vertOUT[(lvert-1)*3+2];

            xcn_parmmg->setVal(vert,0,xc);
            xcn_parmmg->setVal(vert,1,yc);
            xcn_parmmg->setVal(vert,2,zc);

            tag2glob[tag] = vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp;
            glob2tag[vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp] = tag;

            vert++;
        }
        
        int* updateGlobSharedVrtID       = new int[LocationSharedVert.size()];
        int* updateGlobSharedVrtID_red   = new int[LocationSharedVert.size()];
        
        int* originalGlobSharedVrtID     = new int[LocationSharedVert.size()];
        int* originalGlobSharedVrtID_red = new int[LocationSharedVert.size()];
        int ig = 0;

        for(iitm=LocationSharedVert.begin();iitm!=LocationSharedVert.end();iitm++)
        {
            int tag = iitm->first;
            int ra  = iitm->second;
            
            originalGlobSharedVrtID_red[ig] = 0;
            updateGlobSharedVrtID_red[ig]   = 0;
            
            if(ra == world_rank)
            {
                originalGlobSharedVrtID[ig] = tag;
                updateGlobSharedVrtID[ig]   = tag2glob[tag];
            }
            else
            {
                originalGlobSharedVrtID[ig] = 0;
                updateGlobSharedVrtID[ig]   = 0;
                
            }
            ig++;
        }
        
        //std::cout << "Performing an AllReduce of the global vertex IDs of the shared Vertices between partitions." << std::endl;
        
        MPI_Allreduce(originalGlobSharedVrtID,
                      originalGlobSharedVrtID_red,
                      LocationSharedVert.size(),
                      MPI_INT, MPI_SUM, comm);
        
        MPI_Allreduce(updateGlobSharedVrtID,
                      updateGlobSharedVrtID_red,
                      LocationSharedVert.size(),
                      MPI_INT, MPI_SUM, comm);
        
        std::map<int,int> LocationSharedVert_update;
        for(int i=0;i<LocationSharedVert.size();i++)
        {
            LocationSharedVert_update[originalGlobSharedVrtID_red[i]] = updateGlobSharedVrtID_red[i];
        }
    
        int nBndFaces       = PMMG_Face2Ref.size();
        int nShaFaces       = PMMG_SharedFacesOwned.size();
        int nIntFaces       = fmInt.size();
        int nLocFaceNBnd    = fmInt.size()+fmSha.size()+fmShell.size();
        Array<int>* ifnOUT  = new Array<int>(nLocFaceNBnd,8);
        std::map<int,int> updateID;
        int ftot            =  0;
        std::vector<int> fq(3);
        int flag            = -1;
        int flagged         =  0;
        
        //std::cout <<  "right-handed " << world_rank << " " << rh.size() << " " << nShaFaces << " VS " << fmInt.size() << " " << fmSha.size() << " " << fmShell.size() << std::endl;
        
        int nothere2 = 0;
        
        int testElem = 9;
        int testComp = 10;
        
        
        Vert* Vijk = new Vert;
        Vert* VcF = new Vert;

        for(itm=fmInt.begin();itm!=fmInt.end();itm++)
        {
            fid             = itm->first;
            updateID[fid]   = ftot;
            ifnOUT->setVal(ftot,0,3);
            flag            = -1;
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            for(int q=0;q<3;q++)
            {
                int lvert = glob2locVid[itm->second[q]];
                Vert* Vf = new Vert;
                Vf->x  = vertOUT[(lvert-1)*3+0];
                Vf->y  = vertOUT[(lvert-1)*3+1];
                Vf->z  = vertOUT[(lvert-1)*3+2];
                
                VcF->x = VcF->x + Vf->x;
                VcF->y = VcF->y + Vf->y;
                VcF->z = VcF->z + Vf->z;
                
                Vfaces.push_back(Vf);
                
                
                if(LocationSharedVert_update.find(itm->second[q])!=LocationSharedVert_update.end())
                {
                    fq[q] = LocationSharedVert_update[itm->second[q]];
                    flag  = q;
                }
                else
                {
                    fq[q] = tag2glob[itm->second[q]];
                }
            }
            
            VcF->x = VcF->x/3.0;
            VcF->y = VcF->y/3.0;
            VcF->z = VcF->z/3.0;
            
            int gv0 = fq[0];
            int gv1 = fq[1];
            int gv2 = fq[2];
            
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            //ienOUT[curElID] = Elvrts
            for(int u=0;u<ienOUT[lh[fid]].size();u++)
            {
                
                Vijk->x = Vijk->x + vertOUT[(ienOUT[lh[fid]][u]-1)*3];
                Vijk->y = Vijk->y + vertOUT[(ienOUT[lh[fid]][u]-1)*3+1];
                Vijk->z = Vijk->z + vertOUT[(ienOUT[lh[fid]][u]-1)*3+2];
                
            }
            
            Vijk->x = Vijk->x/ienOUT[lh[fid]].size();
            Vijk->y = Vijk->y/ienOUT[lh[fid]].size();
            Vijk->z = Vijk->z/ienOUT[lh[fid]].size();
            
            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            if(orient0 < 0.0)
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv2);
                ifnOUT->setVal(ftot,3,gv1);
                ifnOUT->setVal(ftot,4,0);
            }
            else
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv1);
                ifnOUT->setVal(ftot,3,gv2);
                ifnOUT->setVal(ftot,4,0);
            }
            
            ifnOUT->setVal(ftot,5,rh[fid]);
            ifnOUT->setVal(ftot,6,lh[fid]);
            ifnOUT->setVal(ftot,7,2);

            ftot++;
        }
        
        //Color_SharedOwned
        
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
                int nfid   = locShF2globShF[ofid];
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
        
        DistributedParallelState* rhbefore = new DistributedParallelState(rh.size(),comm);
        DistributedParallelState* lhbefore = new DistributedParallelState(lh.size(),comm);

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
                fid_glo              = locShF2globShF[ft];
                rh[fid_glo]          = adj_ids[itm->first][j];
                
                if(frh.find(fid_glo)==frh.end())
                {
                    adjElements[fid_glo] = adj_ids[itm->first][j];
                    
                    frh.insert(fid_glo);
                    outside++;
                }
                
                telli++;
            }
        }
        
        DistributedParallelState* rhafter  = new DistributedParallelState(rh.size(),comm);
        DistributedParallelState* lhafter  = new DistributedParallelState(lh.size(),comm);
        DistributedParallelState* adjafter = new DistributedParallelState(adjElements.size(),comm);
        
        int nothere =  0;
        int elLh    = -1;
        int elRh    = -1;
        int buthere =  0;
        
        int ftot_inter = ftot;
        for(itm=fmSha.begin();itm!=fmSha.end();itm++)
        {
            fid             = itm->first;
            updateID[fid]   = ftot;
            ifnOUT->setVal(ftot,0,3);
            flag = -1;
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            for(int q=0;q<3;q++)
            {
                int lvert = glob2locVid[itm->second[q]];
                Vert* Vf = new Vert;
                Vf->x  = vertOUT[(lvert-1)*3+0];
                Vf->y  = vertOUT[(lvert-1)*3+1];
                Vf->z  = vertOUT[(lvert-1)*3+2];
                
                VcF->x = VcF->x + Vf->x;
                VcF->y = VcF->y + Vf->y;
                VcF->z = VcF->z + Vf->z;
                
                Vfaces.push_back(Vf);
                
                
                if(LocationSharedVert_update.find(itm->second[q])!=LocationSharedVert_update.end())
                {
                    fq[q] = LocationSharedVert_update[itm->second[q]];
                    flag = q;
                }
                else
                {
                    fq[q] = tag2glob[itm->second[q]];
                }
            }
            
            VcF->x = VcF->x/3.0;
            VcF->y = VcF->y/3.0;
            VcF->z = VcF->z/3.0;
            
            int gv0 = fq[0];
            int gv1 = fq[1];
            int gv2 = fq[2];
            
            if(lhshown.find(fid)!=lhshown.end())
            {
                elLh = lhshown[fid];
                elRh = adjElements[fid];
            }
            
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            //ienOUT[curElID] = Elvrts
            for(int u=0;u<ienOUT[elLh].size();u++)
            {
                
                Vijk->x = Vijk->x + vertOUT[(ienOUT[elLh][u]-1)*3];
                Vijk->y = Vijk->y + vertOUT[(ienOUT[elLh][u]-1)*3+1];
                Vijk->z = Vijk->z + vertOUT[(ienOUT[elLh][u]-1)*3+2];
                
            }
            
            Vijk->x = Vijk->x/ienOUT[elLh].size();
            Vijk->y = Vijk->y/ienOUT[elLh].size();
            Vijk->z = Vijk->z/ienOUT[elLh].size();
            
            
            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            
            if(orient0 < 0.0)
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv2);
                ifnOUT->setVal(ftot,3,gv1);
                ifnOUT->setVal(ftot,4,0);
            }
            else
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv1);
                ifnOUT->setVal(ftot,3,gv2);
                ifnOUT->setVal(ftot,4,0);
            }
            
            ifnOUT->setVal(ftot,5,elRh);
            ifnOUT->setVal(ftot,6,elLh);
            ifnOUT->setVal(ftot,7,2);

            ftot++;
        }
        
        
        std::map<int,int> shelltag2glob;
        int alhere      = 0;
        int nothere3    = 0;
        int tetshe      = 0;
        int sna         = 0;
        int notsh       = 0;
        int tellie      = 0;
        int ormin       = 0;
        int tellie2     = 0;
        for(itm=fmShell.begin();itm!=fmShell.end();itm++)
        {
            fid             = itm->first;
            updateID[fid]   = ftot;
            ifnOUT->setVal(ftot,0,3);
            flag = -1;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            std::vector<Vert*> Vfaces;
            std::vector<int> reference(3);
            for(int q=0;q<3;q++)
            {
                int lvert = glob2locVid[itm->second[q]];
                Vert* Vf = new Vert;
                Vf->x  = vertOUT[(lvert-1)*3+0];
                Vf->y  = vertOUT[(lvert-1)*3+1];
                Vf->z  = vertOUT[(lvert-1)*3+2];
                
                reference[q] = refOUT[(lvert-1)];
                
                VcF->x = VcF->x + Vf->x;
                VcF->y = VcF->y + Vf->y;
                VcF->z = VcF->z + Vf->z;
                
                Vfaces.push_back(Vf);
                
                if(LocationSharedVert_update.find(itm->second[q])!=LocationSharedVert_update.end())
                {
                    fq[q] = LocationSharedVert_update[itm->second[q]];
                    
                    int shelltag = tag2shelltag_glob[itm->second[q]];
                    
                    if(shelltag2glob.find(shelltag)==shelltag2glob.end())
                    {
                        shelltag2glob[shelltag] = fq[q];
                    }
                                        
                    flag = q;
                }
                else
                {
                    fq[q] = tag2glob[itm->second[q]];
                    
                    int shelltag = tag2shelltag_glob[itm->second[q]];
                    
                    if(shelltag2glob.find(shelltag)==shelltag2glob.end())
                    {
                        shelltag2glob[shelltag] = fq[q];
                    }
                    else
                    {
                        alhere++;
                    }
                }
            }
            
            VcF->x = VcF->x/3.0;
            VcF->y = VcF->y/3.0;
            VcF->z = VcF->z/3.0;
            
            
            int gv0 = fq[0];
            int gv1 = fq[1];
            int gv2 = fq[2];
            
            
            if(rhshell.find(fid)==rhshell.end())
            {
                elRh = rh[fid];
                elLh = lhshell[fid];
                notsh++;
            }
            else
            {
                elRh = rhshell[fid];
                elLh = lhshell[fid];
            }
            
            
            Vert* Vijk   = new Vert;
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            
            // compute element center;
            //ienOUT[curElID] = Elvrts
            for(int u=0;u<ienOUT[elLh].size();u++)
            {
                Vijk->x = Vijk->x + vertOUT[(ienOUT[elLh][u]-1)*3];
                Vijk->y = Vijk->y + vertOUT[(ienOUT[elLh][u]-1)*3+1];
                Vijk->z = Vijk->z + vertOUT[(ienOUT[elLh][u]-1)*3+2];
            }
            
            Vijk->x = Vijk->x/ienOUT[elLh].size();
            Vijk->y = Vijk->y/ienOUT[elLh].size();
            Vijk->z = Vijk->z/ienOUT[elLh].size();
            
            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            
            
            if(orient0 < 0.0)
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv2);
                ifnOUT->setVal(ftot,3,gv1);
                ifnOUT->setVal(ftot,4,0);
            }
            else
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv1);
                ifnOUT->setVal(ftot,3,gv2);
                ifnOUT->setVal(ftot,4,0);
            }
            
            ifnOUT->setVal(ftot,5,elRh);
            ifnOUT->setVal(ftot,6,elLh);
            ifnOUT->setVal(ftot,7,2);
            
            Vfaces.clear();

            ftot++;
        }
        
        
        std::map<int,int> shelltag2glob_glob = AllGatherMap(shelltag2glob,comm);
        
        int mapSizeLoc = shelltag2glob.size();
        DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,comm);
        int mapSizeTot = distrimap->getNel();
        int* shelltag_loc = new int[mapSizeLoc];
        int* tag_loc = new int[mapSizeLoc];
        int* shelltag_tot = new int[mapSizeTot];
        int* tag_tot = new int[mapSizeTot];

        int i = 0;
        
        std::map<int,int>::iterator itred;
        
        for(itred=shelltag2glob.begin();itred!=shelltag2glob.end();itred++)
        {
            shelltag_loc[i] = itred->first;
            tag_loc[i] = itred->second;
            i++;
        }
        
        int* offsets = distrimap->getOffsets();
        int* nlocs   = distrimap->getNlocs();
        
        
        MPI_Allgatherv(shelltag_loc,
                       mapSizeLoc,
                       MPI_INT,
                       shelltag_tot,
                       nlocs,
                       offsets,
                       MPI_INT, comm);
        
        
        MPI_Allgatherv(tag_loc,
                       mapSizeLoc,
                       MPI_INT,
                       tag_tot,
                       nlocs,
                       offsets,
                       MPI_INT, comm);
        
        int key,val;
        std::map<int,int> shelltag2glob_global;
        for(int i=0;i<mapSizeTot;i++)
        {
            key = tag_tot[i];
            val = shelltag_tot[i];
            
            if(shelltag2glob_global.find(val)==shelltag2glob_global.end())
            {
                shelltag2glob_global[val] = key;
            }
        }
        
        
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
        
        int cc       = 0;
        int notFound = 0;
        int tagnew   = 0;
        std::map<int,int> t2g_tmp;
        std::map<int,int> testmap;
        int notfo  = 0;
        int fo     = 0;
        int notfo2 = 0;
        int redflag = 0;
    
        std::map<int,int> collect;
        std::map<int,int> collect2;
        int cn = 0;
        int fshell = 0;
        for( prit=int_face2node_prism.begin();prit!=int_face2node_prism.end();prit++)
        {
            int gfid    = prit->first;
            int npf     = prit->second.size();
            std::vector<int> fce(npf);
            ifnOUT_prism->setVal(fptot,0,npf);
            
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            fshell = 0;
            for(int g=0;g<npf;g++)
            {
                int oldtag = prit->second[g];
                
                if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end() &&
                   shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
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
                else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end() &&
                        shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
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
                else if(shellvertOriginalTag2ref_Glob.find(oldtag)!=shellvertOriginalTag2ref_Glob.end())
                {
                    int ref     = shellvertOriginalTag2ref_Glob[oldtag];
                    int globid  = shelltag2glob_global[ref];
                    fce[g]      = globid;
                    fshell = 1;
                    if(globid>nPrismVerts_tmp)
                    {
                        if(collect2.find(globid)==collect2.end())
                        {
                            collect2[globid] = cn;
                            cn++;
                        }
                    }
                    if(collect.find(globid)==collect.end())
                    {
                        collect[globid] = g;
                    }
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

            
            ifnOUT_prism->setVal(fptot,5,rhp[gfid]);
            ifnOUT_prism->setVal(fptot,6,lhp[gfid]);
            ifnOUT_prism->setVal(fptot,7,2);

            fptot++;
            pid++;
        }
    

        int cnttt     = 0;
        int notany    = 0;
        for( prit=shared_face2node_prism.begin();prit!=shared_face2node_prism.end();prit++)
        {
            int gfid  = prit->first;
            int npf   = prit->second.size();
            
            ifnOUT_prism->setVal(fptot,0,npf);
            std::vector<int> fce(npf);
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            for(int g=0;g<npf;g++)
            {
                int oldtag = prit->second[g];
                if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end() &&
                   shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
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
                else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end() &&
                        shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
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
                else if(shellvertOriginalTag2ref_Glob.find(oldtag)!=shellvertOriginalTag2ref_Glob.end())
                {
                    int ref     = shellvertOriginalTag2ref_Glob[oldtag];
                    int globid = shelltag2glob_global[ref];
                    fce[g]     = globid;
                    if(globid>nPrismVerts_tmp)
                    {
                        if(collect2.find(globid)==collect2.end())
                        {
                            collect2[globid] = cn;
                            cn++;
                        }
                    }
                    if(collect.find(globid)==collect.end())
                    {
                        collect[globid] = g;
                    }
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
            
            fptot++;
            pid++;
        }
        

        
        // End reducing the map.
        
        DistributedParallelState* distftot   = new DistributedParallelState(ftot,comm);
        DistributedParallelState* distfptot  = new DistributedParallelState(fptot,comm);

        int nTotInteriorFaces          = distftot->getNel();
        int* TotIntFaces_offsets       = distftot->getOffsets();
        
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
        
        for(bit=bcmap.begin();bit!=bcmap.end();bit++)
        {
            if(bcids_tot.find(bit->first)==bcids_tot.end())
            {
                bcids_tot.insert(bit->first);
                Lbcs.push_back(bit->first);
            }
        }
        
        for(bit=pbcmap.begin();bit!=pbcmap.end();bit++)
		{
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
        q = 0;
        for(entry=bcsToT.begin();entry!=bcsToT.end();entry++)
        {
            int ee = *entry;
            int Nbft = 0;
            int Nbfp = 0;
            
            if(bcmap.find(ee)!=bcmap.end())
            {
                Nbft = bcmap[ee].size();
            }
            if(pbcmap.find(ee)!=pbcmap.end())
            {
            	Nbfp = pbcmap[ee].size();
            }

            bcid[q]  = ee;
            nlbc[q]  = Nbft+Nbfp;
            q++;
        }
        
        
        
        std::vector<Array<int>* > bcArrays;
        std::map<int,int> bcsizing;
        std::vector<int> bci_offsets;
        std::vector<int> bciTot_offsets;
        int nTotBCFaces = 0;
        int nTotBCFaces_offset = 0;
    
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

            if(bcmap.find(bc_id)!=bcmap.end())
            {
                Nbft = bcmap[bc_id].size();
            }
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
                    VcF->x = 0.0;
                    VcF->y = 0.0;
                    VcF->z = 0.0;
                    
                    for(int g=0;g<nppf;g++)
                    {
                        int oldtag = bc_face2node_prism[bcface][g];
                        
                        if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end() &&
                           shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
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
                        else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end() &&
                                shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
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
                        else if(shellvertOriginalTag2ref_Glob.find(oldtag)!=shellvertOriginalTag2ref_Glob.end())
                        {
                            int ref     = shellvertOriginalTag2ref_Glob[oldtag];
                            int globid  = shelltag2glob_global[ref];
                            fce[g]      = globid;
                            if(globid>nPrismVerts_tmp)
                            {
                                if(collect2.find(globid)==collect2.end())
                                {
                                    collect2[globid] = cn;
                                    cn++;
                                }
                            }
                            if(collect.find(globid)==collect.end())
                            {
                                collect[globid] = g;
                            }
                        }
                    }
                    
                    VcF->x = VcF->x/nppf;
                    VcF->y = VcF->y/nppf;
                    VcF->z = VcF->z/nppf;
                    
                    int leftEl  = lhp[bcface];
                    int leftTag = gE2tagE[leftEl];
                    
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
                    
					ifn_bc_i->setVal(fbc,5,0);
					ifn_bc_i->setVal(fbc,6,lhp[bcface]);
					ifn_bc_i->setVal(fbc,7,bc_id);
					
					fbc++;
				}
                
                
            }
            
            if(Nbft!=0 )
            {
                for(int q=0;q<Nbft;q++)
                {
                    int bcface = bcmap[bc_id][q];
                    ifn_bc_i->setVal(fbc,0,3);

                    std::vector<Vert*> Vfaces;
                    
                    VcF->x = 0.0;
                    VcF->y = 0.0;
                    VcF->z = 0.0;
                    flag = -1;
                    for(int w=0;w<3;w++)
                    {
                        
                        int lvert = glob2locVid[fmBnd[bcface][w]];
                        Vert* Vf = new Vert;
                        Vf->x  = vertOUT[(lvert-1)*3+0];
                        Vf->y  = vertOUT[(lvert-1)*3+1];
                        Vf->z  = vertOUT[(lvert-1)*3+2];
                        
                        VcF->x = VcF->x + Vf->x;
                        VcF->y = VcF->y + Vf->y;
                        VcF->z = VcF->z + Vf->z;
                        
                        Vfaces.push_back(Vf);
                        
                        if(LocationSharedVert_update.find(fmBnd[bcface][w])!=LocationSharedVert_update.end())
                        {
                            fq[w] = LocationSharedVert_update[fmBnd[bcface][w]];
                            ifn_bc_i->setVal(fbc,w+1,LocationSharedVert_update[fmBnd[bcface][w]]);
                            if(ifn_bc_i->getVal(fbc,w+1)==0)
                            {
                                std::cout << " TET Boundary V zero in the tet local LocationSharedVert_update " << world_rank << " " << fbc << std::endl;
                            }
                            flag = w;
                        }
                        else
                        {
                            fq[w] = tag2glob[fmBnd[bcface][w]];
                            ifn_bc_i->setVal(fbc,w+1,tag2glob[fmBnd[bcface][w]]);
                            if(ifn_bc_i->getVal(fbc,w+1)==0)
                            {
                                std::cout << " TET Boundary V zero in the tet local tag2glob " << world_rank << " " << fbc << std::endl;
                            }
                        }
                    }
                    
                    int elLh = lhbnd[bcface];
                    
                    Vijk->x = 0.0;
                    Vijk->y = 0.0;
                    Vijk->z = 0.0;
                    // compute element center;
                    //ienOUT[curElID] = Elvrts
                    for(int u=0;u<ienOUT[elLh].size();u++)
                    {
                        
                        Vijk->x = Vijk->x + vertOUT[(ienOUT[elLh][u]-1)*3];
                        Vijk->y = Vijk->y + vertOUT[(ienOUT[elLh][u]-1)*3+1];
                        Vijk->z = Vijk->z + vertOUT[(ienOUT[elLh][u]-1)*3+2];
                        
                    }
                    
                    Vijk->x = Vijk->x/ienOUT[elLh].size();
                    Vijk->y = Vijk->y/ienOUT[elLh].size();
                    Vijk->z = Vijk->z/ienOUT[elLh].size();
                    
                    double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);

                    if(orient0 < 0.0)
                    {
                        std::cout << "Weve got negative faces " << orient0 << std::endl;
                    }
                    
                    Vfaces.clear();
                    
                    ifn_bc_i->setVal(fbc,4,0);
                    ifn_bc_i->setVal(fbc,5,0);
                    ifn_bc_i->setVal(fbc,6,lhbnd[bcface]);
                    ifn_bc_i->setVal(fbc,7,bc_id);
                    
                    fbc++;
                }
            }
            
            
            int nbt = Nbft+Nbfp;
            
            bcsizing[bc_id] = NelTot_bci;
            bci_offsets.push_back(offsetbci);
            bciTot_offsets.push_back(nTotBCFaces_offset);
            bcArrays.push_back(ifn_bc_i);
            
            nTotBCFaces_offset = nTotBCFaces_offset + NelTot_bci;
            nTotBCFaces        = nTotBCFaces + NelTot_bci;
            
            
        }

        int nPrismOUT = parmmg_iet_prisms->getNrow();
        DistributedParallelState* distTetra      = new DistributedParallelState(nTetrahedraOUT,comm);
        DistributedParallelState* distPrism      = new DistributedParallelState(nPrismOUT,comm);
        DistributedParallelState* distTetraVerts = new DistributedParallelState(xcn_parmmg->getNrow(),comm);
        DistributedParallelState* distPrismVerts = new DistributedParallelState(xcn_prisms_int->getNrow()+xcn_prisms_shared->getNrow(),comm);
        
        
        int ToTElements_prism           = distPrism->getNel();
        int ToTElements_offset_prism    = distPrism->getOffsets()[world_rank];
        
        int ToTElements                 = distTetra->getNel();
        int ToTElements_offset          = distTetra->getOffsets()[world_rank];
        int nTotTetraVerts_v2           = distTetraVerts->getNel();

        int nTotElements                = ToTElements_prism+ToTElements;
        int nTotFaces                   = nTotInteriorFaces_prism + nTotInteriorFaces + nTotBCFaces;
        int nTotIntFaces                = nTotInteriorFaces_prism + nTotInteriorFaces;
        
        int nTotPrismVerts_v2           = distPrismVerts->getNel();
        int nTotPrismIntVerts_v2        = distPrismIntVerts->getNel();
        int nTotPrismShaVerts_v2        = distPrismShaVerts->getNel();
    
        int TotPrismVerts_offset_int    = distPrismIntVerts->getOffsets()[world_rank];
        int TotPrismVerts_offset_sha    = distPrismShaVerts->getOffsets()[world_rank];

        int TotPrismVerts_offset        = distPrismVerts->getOffsets()[world_rank];
        
        int nTotVertsPrismTetra = nTotPrismVerts_v2+nTotTetraVerts_v2;
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
        value = nTotPrismVerts_v2+nTotTetraVerts_v2;//ToTVrts;
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
        countH5[0]  = parmmg_iet->getNrow();
        countH5[1]  = parmmg_iet->getNcol();

        offsetH5[0] = ToTElements_prism+ToTElements_offset;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        //filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, parmmg_iet->data);
        delete parmmg_iet;
        //====================================================================================
        // Add xcn map to the grid.h5 file
        //====================================================================================
        

        
        
        dimsf[0] = nTotPrismIntVerts_v2+nTotPrismShaVerts_v2+nTotTetraVerts_v2;
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
        
        
        offsetH5[0] = nTotPrismIntVerts_v2+TotPrismVerts_offset_sha;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_prisms_shared->data);
        //std::cout << "world " << nTotPrismIntVerts_v2 << " " << TotPrismVerts_offset_sha << " " << xcn_prisms_shared->getNrow() << std::endl;
        delete xcn_prisms_shared;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
        countH5[0]  = xcn_parmmg->getNrow();
        countH5[1]  = xcn_parmmg->getNcol();

        offsetH5[0] = nTotPrismIntVerts_v2+nTotPrismShaVerts_v2+ToTVrts_offset;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);
        //filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_parmmg->data);
        delete xcn_parmmg;
        
       
        //===================================================================================
//        int nTotInteriorFaces          = distftot->getNel();
//        int* TotIntFaces_offsets       = distftot->getOffsets();
//
//        int nTotInteriorFaces_prism    = distfptot->getNel();
//        int* TotIntFaces_offsets_prism = distfptot->getOffsets();
        
        dimsf[0]  = nTotInteriorFaces_prism+nTotInteriorFaces+nTotBCFaces;
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
        
        countH5[0]  = ifnOUT->getNrow();
        countH5[1]  = dimsf[1];

        offsetH5[0] = nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank];
        offsetH5[1] = 0;
        //std::cout << "nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank ] " << TotIntFaces_offsets_prism[world_rank] << " " << ifnOUT_prism->getNrow() << " " << nTotInteriorFaces_prism << " " << TotIntFaces_offsets[world_rank ] << " " << ifnOUT->getNrow() << world_rank << " -> " <<nTotInteriorFaces_prism << " + " << nTotInteriorFaces << " " << ftot << " " << fptot  << " " << world_rank << " " << nTotInteriorFaces_prism+nTotInteriorFaces << std::endl;

        memspace     = H5Screate_simple(2, countH5, NULL);
        filespace     = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
//
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, ifnOUT->data);

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
            
            offsetH5[0]   = nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i];
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
        
        /**/
    }
    
    
     
    
    MPI_Finalize();
    
}

