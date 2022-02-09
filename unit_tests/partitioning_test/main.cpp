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





std::map<int,Array<double>*> ReadReferenceData(int world_rank)
{
    std::map<int,Array<double>*> output;
    std::ifstream fin_v;
    fin_v.open("test_data/ref_data/UVariaValues_"+ std::to_string(world_rank) + ".txt");
    
    std::vector<double> row_v(2);
    while(fin_v >> row_v[0] >> row_v[1])
    {
        Array<double>* entry = new Array<double>(1,1);
        entry->setVal(0,0,row_v[1]);
        int key = (int) row_v[0];
        output[key] = entry;
        
    }
    fin_v.close();
    
    return output;
    
}






void OutputMesh_PMMG(int nV, double* VertOUT, int nE, int* tetraOUT, string fname)
{
    int pos;
    
    std::ofstream myfile;
    myfile.open(fname);
    myfile << "TITLE=\"new_volume.tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << nV << ", E = " << nE << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<nV;i++)
    {
        pos = 3*i;
        myfile << VertOUT[pos] << " " <<VertOUT[pos+1] << " " << VertOUT[pos+2] <<  std::endl;
    }
    for(int i=0;i<nE;i++)
    {
        pos=4*i;
        
        myfile << tetraOUT[pos] << " " << tetraOUT[pos+1]  << " " << tetraOUT[pos+2]  << " " << tetraOUT[pos+3]  << std::endl;
        
    }
    myfile.close();

}


void OutputMesh_PMMG_V2(int nV, double* VertOUT, int nE, int* tetraOUT, string fname)
{
    int pos;
    
    std::vector<std::vector<int> > plotnodes;
    for(int i=0;i<nE;i++)
    {
        pos=4*i;
        if(VertOUT[(tetraOUT[pos]-1)*3+2]<0.0 &&
           VertOUT[(tetraOUT[pos+1]-1)*3+2]<0.0 &&
           VertOUT[(tetraOUT[pos+2]-1)*3+2]<0.0 &&
           VertOUT[(tetraOUT[pos+3]-1)*3+2]<0.0)
        {
            std::vector<int> roww(4);
            roww[0] = tetraOUT[pos];
            roww[1] = tetraOUT[pos+1];
            roww[2] = tetraOUT[pos+2];
            roww[3] = tetraOUT[pos+3];
            
            plotnodes.push_back(roww);
        }
        
    }
        
    if(plotnodes.size() > 0 )
    {
        std::ofstream myfile;
        myfile.open(fname);
        myfile << "TITLE=\"new_volume.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
        myfile <<"ZONE N = " << nV << ", E = " << plotnodes.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

        for(int i=0;i<nV;i++)
        {
            pos = 3*i;
            myfile << VertOUT[pos] << " " <<VertOUT[pos+1] << " " << VertOUT[pos+2] <<  std::endl;
        }
        for(int i=0;i<plotnodes.size();i++)
        {
            pos=4*i;
    //        if(VertOUT[tetraOUT[pos]*3+2]>0.25 &&
    //           VertOUT[tetraOUT[pos+1]*3+2]>0.25 &&
    //           VertOUT[tetraOUT[pos+2]*3+2]>0.25 &&
    //           VertOUT[tetraOUT[pos+3]*3+2]>0.25)
    //        {
            
                //myfile << tetraOUT[pos] << " " << tetraOUT[pos+1]  << " " << tetraOUT[pos+2]  << " " << tetraOUT[pos+3]  << std::endl;
            myfile << plotnodes[i][0] << " " << plotnodes[i][1]  << " " << plotnodes[i][2]  << " " << plotnodes[i][3]  << std::endl;
    //        }
            
        }
        myfile.close();
    }
    
}


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



std::map<int,std::vector<double> > AllGatherMapDoubleVec(std::map<int,std::vector<double> > mappie, MPI_Comm mpi_comm)
{
    int mapSizeLoc = mappie.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    
    DistributedParallelState* distrimapVal = new DistributedParallelState(mapSizeLoc*3,mpi_comm);
    
    int mapSizeTot = distrimap->getNel();
    int* key_loc = new int[mapSizeLoc];
    double* val_loc = new double[mapSizeLoc*3];
    int* key_tot = new int[mapSizeTot];
    double* val_tot = new double[mapSizeTot*3];
    int i = 0;
    
    std::map<int,std::vector<double> >::iterator itred;
    for(itred=mappie.begin();itred!=mappie.end();itred++)
    {
        key_loc[i] = itred->first;
        int nrow   = itred->second.size();
        for(int q=0;q<nrow;q++)
        {
            val_loc[i*3+q] = itred->second[q];
            //std::cout << "itred->second[q] " << itred->second[q] << std::endl;
        }
        
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
    
    int* offsetsVal = distrimapVal->getOffsets();
    int* nlocsVal   = distrimapVal->getNlocs();
    
    MPI_Allgatherv(val_loc,
                   mapSizeLoc*3,
                   MPI_DOUBLE,
                   val_tot,
                   nlocsVal,
                   offsetsVal,
                   MPI_DOUBLE, mpi_comm);
    
    int key,val;
    std::map<int,std::vector<double> > mappie_glob;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = key_tot[i];
        
        std::vector<double> values(3);
        for(int q=0;q<3;q++)
        {
            values[q] = val_tot[i*3+q];
        }
        
        if(mappie_glob.find(key)==mappie_glob.end())
        {
            
            mappie_glob[key] = values;
            //std::cout << "itred->second[q] " << val[0] << " " << val[1] << " " << val[2] << std::endl;
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
    clock_t t0_met = clock();

    int ier,opt;
    int debug = 1;
//    const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
//    const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
//    const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    
    const char* fn_grid="test_data/inputs/grid.h5";
    const char* fn_conn="test_data/inputs/conn.h5";
    const char* fn_data="test_data/inputs/data.h5";
    const char* fn_metric = "test_data/inputs/metric.inp";
    // inputs 
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);

    
    //===========================================================================
    
    double hgrad         = metric_inputs[0];
    double hmin          = metric_inputs[1];
    double hmax          = metric_inputs[2];
    double MetScale      = metric_inputs[3];
    int ReadFromStats    = metric_inputs[4];
    int RunWakRefinement = metric_inputs[5];
    double hwake         = metric_inputs[6];
    int niter            = metric_inputs[7];
    int recursive	     = metric_inputs[8];
    int extended         = metric_inputs[9];
    
    if(world_rank == 0)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "============== Metric parameters ==================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "Nproc	    = " << world_size << std::endl;
        std::cout << "hgrad     = " << hgrad << std::endl;
        std::cout << "hmin      = " << hmin << std::endl;
        std::cout << "hmax      = " << hmax << std::endl;
        std::cout << "MetScale  = " << MetScale << std::endl;
        std::cout << "NiterPart = " << niter << std::endl;
        if(ReadFromStats == 0)
        {
            std::cout << "Reading statistics? -> NO (5th entry in the metric.inp file is set to 0.)" << std::endl;
            std::cout << "The metric is reconstructed based on instantaneous Mach number"<<std::endl;
        }
        if(ReadFromStats == 1)
        {
            std::cout << "Reading statistics? -> YES (5th entry in the metric.inp file is set to 1.)" << std::endl;
            std::cout << "The metric is reconstructed based on the mean of Mach number."<<std::endl;

        }
        if(RunWakRefinement==0)
        {
            std::cout << "Wake refinement is switch OFF. (6th entry in the metric.inp file is set to 0. hwake, the 7th entry defined in the metric.inp file, is being ignored)" << std::endl;
            
        }
        if(RunWakRefinement==1)
        {
            std::cout << "Wake refinement is switch ON with hwake = " << hwake << "(6th entry in the metric.inp file is set to 1 and hwake is set equal to the 7th entry defined in the metric.inp file.) " << std::endl;
        }
        
        std::cout << "===================================================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "  " << std::endl;
    }
    //===========================================================================
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve       = us3d->xcn->getNglob();
    
    int Nel_part  = us3d->ien->getNrow();
    
    Array<double>* Ui = new Array<double>(Nel_part,1);
    Array<double>* TKEi;
    int varia = 4;
    double TKE, MState;
    
    if(ReadFromStats==0)
    {
        for(int i=0;i<Nel_part;i++)
        {
            MState   = us3d->interior->getVal(i,0);
            Ui->setVal(i,0,MState);
        }
    }
    
    if(ReadFromStats==1)
    {
        TKEi = new Array<double>(Nel_part,1);

        for(int i=0;i<Nel_part;i++)
        {
            TKE      = us3d->interior->getVal(i,0);
            MState   = us3d->interior->getVal(i,1);
            Ui->setVal(i,0,MState);
            TKEi->setVal(i,0,TKE);
        }
    }
    
    
    delete us3d->interior;
 
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }
    int ngho = us3d->ghost->getNrow();
    int ngval;
    MPI_Allreduce(&ngho, &ngval, 1, MPI_INT, MPI_MAX, comm);
       
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
    
    std::vector<int> LocElem                    = P->getLocElem();
    std::vector<double> Uvaria                  = P->getLocElemVaria();
    std::map<int,Array<double>*> Uvaria_map     = P->getLocAndAdjElemVaria();
    
    std::map<int,Array<double>*> Uvaria_map2;
    double UvariaV              = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid   = LocElem[i];
        UvariaV   = Uvaria[i];

        Array<double>* Uarr = new Array<double>(1,1);
        Uarr->setVal(0,0,UvariaV);
        Uvaria_map2[gid] = Uarr;
    }
    
    P->AddStateVecForAdjacentElements(Uvaria_map2,1,comm);
    
    std::map<int,Array<double>*> Uvaria_ref = ReadReferenceData(world_rank);

    std::map<int,Array<double>* >::iterator pltV;
    int correctVal   = 0;
    int incorrectVal = 0;
    int notThere     = 0;
    for(pltV=Uvaria_map.begin();pltV!=Uvaria_map.end();pltV++)
    {
        if(Uvaria_ref.find(pltV->first)!=Uvaria_ref.end())
        {
            double err = fabs(pltV->second->getVal(0,0)-Uvaria_ref[pltV->first]->getVal(0,0));
            if(err > 1.0e-05)
            {
                incorrectVal++;
            }
            else
            {
                correctVal++;
            }
        }
    }

    
    int incorrectValGlobal;
    MPI_Allreduce(&incorrectVal, &incorrectValGlobal, 1, MPI_INT, MPI_SUM, comm);
    int notThereGlobal;
    MPI_Allreduce(&notThere, &notThereGlobal, 1, MPI_INT, MPI_SUM, comm);
    
    if(world_rank == 0)
    {
        if(incorrectValGlobal==0)
        {
            std::cout << "Partitioning test has PASSED." << std::endl;
            
            std::ofstream myfile;
            myfile.open("partition_test.PASSED");
            myfile << "SUCCES!" << std::endl;
            myfile.close();
            
        }
        else
        {
            std::cout << "Partitioning test has FAILED." << " " << incorrectValGlobal << std::endl;
            std::ofstream myfile;
            myfile.open("partition_test.FAILED");
            myfile << "SUCCES!" << std::endl;
            myfile.close();
        }
    }
    
    MPI_Finalize();
    
}

