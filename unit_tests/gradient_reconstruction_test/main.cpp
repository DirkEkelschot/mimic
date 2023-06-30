#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
#include <math.h>
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// This is basically textbook recursive merge sort using std::merge_inplace
// but it considers the offsets of segments that are already sorted


std::vector<double> ReadReferenceErrors(const char* fn_errors)
{
    std::ifstream fin;
    fin.open(fn_errors);
    if(!fin.is_open())
    {
        std::cout << "Error:: Make sure there is a errors.ref file in the directory where test.cpp resides. "<<std::endl;
        exit(0);
    }
    
    double v=0.0;
    std::vector<double> errors_inputs;
    int t=0;
    while(fin >> v)
    {
        errors_inputs.push_back(v);
       t++;
    }
    return errors_inputs;
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
    int debug = 0;
//    const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
//    const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
//    const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    
    
    
    const char* fn_grid   = "box_inputs/boxes_cas/grid_h0p01.h5";
    const char* fn_conn   = "box_inputs/boxes_cas/conn_h0p01.h5";
    const char* fn_metric = "inputs_h0p05/metric.inp";
    
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
    US3D* us3d    = ReadUS3DGrid(fn_conn,fn_grid,ReadFromStats,comm,info);

    //US3D* us3d  = ReadUS3DData(fn_conn,fn_grid,fn_data,ReadFromStats,comm,info);
    int Nve       = us3d->xcn->getNglob();
    
    int Nel_part  = us3d->ien->getNrow();
    int Nel_glob  = us3d->ien->getNglob();

    Array<double>* Ui = new Array<double>(Nel_part,1);
    
    for(int i=0;i<Nel_part;i++)
    {
        Ui->setVal(i,0,1.0);
    }

    delete us3d->interior;
 
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,0.0);
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
      
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief,
                                 us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife,
                                 us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate,
                                 ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, us3d->ie_tetCnt, comm);
    if(world_rank == 0)
    {
        std::cout << "Done Partitioning..." << std::endl;
    }

    
    std::vector<int> LocElem    = P->getLocElem();
    std::vector<int> LocAndAdjElem    = P->getLocAndAdjElem();

    i_part_map* ien_pmap        = P->getIENpartmap();
    i_part_map* iee_pmap        = P->getIEEpartmap();
    i_part_map* ief_pmap        = P->getIEFpartmap();
    i_part_map* ifn_pmap        = P->getIFNpartmap();
    
    std::map<int,int> gV2lV     = P->getGlobalVert2LocalVert();
    std::vector<Vert*> locVs    = P->getLocalVerts();
    
    std::map<int,Array<double>*> Uvaria_map;
    std::map<int,Array<double>*> dUdxAnalytical;
    std::map<int,Array<double>*> dUdyAnalytical;
    std::map<int,Array<double>*> dUdzAnalytical;
    std::map<int,Array<double>*> dU2dx2Analytical;
    std::map<int,Array<double>*> dU2dxyAnalytical;
    std::map<int,Array<double>*> dU2dxzAnalytical;
    
    std::map<int,double> gbMap;
    std::map<int,double> gbMap_dUdx;
    std::map<int,double> gbMap_dUdy;
    std::map<int,double> gbMap_dUdz;
    std::map<int,double> gbMap_dU2dx2;
    std::map<int,double> gbMap_dU2dxy;
    std::map<int,double> gbMap_dU2dxz;
    
    std::vector<Vert*> face;
    double rdotn;
    Vec3D* v0 = new Vec3D;
    Vec3D* v1 = new Vec3D;
    Vec3D* n0 = new Vec3D;
    int loc_vid;
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    double dUdx,dUdy,dUdz,dU2dx2,dU2dxy,dU2dxz;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid   = LocElem[i];
        int nvrts = ien_pmap->i_map[gid].size();
        double* Pv = new double[nvrts*3];
        for(int q=0;q<nvrts;q++)
        {
            int gvid = ien_pmap->i_map[gid][q];
            int lvid = gV2lV[gvid];
            Pv[q*3+0] = locVs[lvid]->x;
            Pv[q*3+1] = locVs[lvid]->y;
            Pv[q*3+2] = locVs[lvid]->z;
        }
        
        Vert* Vijk = ComputeCentroidCoord(Pv,nvrts);

//        double U    = 0.1*sin(50*Vm->x*Vm->z)+atan(0.1/((sin(5.0*Vm->y)-2.0*Vm->x*Vm->z)));
//        double nom  = (Vm->x*Vm->x*Vm->z*Vm->z-Vm->x*Vm->z*sin(5.0*Vm->y)+0.25*sin(5.0*Vm->y)*sin(5.0*Vm->y)+0.0025);
//        double dUdx = 0.05*Vm->z/nom+0.5*Vm->z*cos(5.0*Vm->x*Vm->z);
        //double dUdy = -0.125*cos(5.0*Vm->y)/nom;
        //double dUdz = 0.05*Vm->x/nom+0.5*Vm->x*cos(5.0*Vm->x*Vm->z);
        double vmx = Vijk->x;
        double vmy = Vijk->y;
        double vmz = Vijk->z;

//        double U = 2.0*vmx*vmx+0.5*vmy+0.25*vmz;
//        double dUdx = 4.0*vmx;
//        double dUdy = 0.5;
//        double dUdz = 0.25;
//        double dU2dx2 = 4.0;
        // sqrt(x*x+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5))
        double r    = sqrt(vmx*vmx+(vmy-0.5)*(vmy-0.5)+(vmz-0.5)*(vmz-0.5));
        double U    = 0.1*tanh(50*(r-0.5))+1.0;
        
        dUdx = 5*vmx/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
        dUdy = 5*(vmy-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
        dUdz = 5*(vmz-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
        dU2dx2 = (5*(-100*vmx*vmx*r*tanh(50*(r-0.5)))+(vmy-0.5)*(vmy-0.5)+(vmz-0.5)*(vmz-0.5))/(cosh(50*(r-0.5))*cosh(50*(r-0.5))*pow(r,3.0/2.0));
        dU2dxy =  5*vmx*(vmy-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));
        dU2dxz = 5*vmx*(vmz-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));

        Array<double>* Uarr = new Array<double>(1,1);
        Uarr->setVal(0,0,U);
        Uvaria_map[gid] = Uarr;
        Array<double>* dUdx_a = new Array<double>(1,1);
        dUdx_a->setVal(0,0,dUdx);
        dUdxAnalytical[gid] = dUdx_a;
        Array<double>* dUdy_a = new Array<double>(1,1);
        dUdy_a->setVal(0,0,dUdy);
        dUdyAnalytical[gid] = dUdy_a;
        Array<double>* dUdz_a = new Array<double>(1,1);
        dUdz_a->setVal(0,0,dUdz);
        Array<double>* dU2dx2_a = new Array<double>(1,1);
        dU2dx2_a->setVal(0,0,dU2dx2);
        dU2dx2Analytical[gid] = dU2dx2_a;
        Array<double>* dU2dxy_a = new Array<double>(1,1);
        dU2dxy_a->setVal(0,0,dU2dxy);
        dU2dxyAnalytical[gid] = dU2dxy_a;
        Array<double>* dU2dxz_a = new Array<double>(1,1);
        dU2dxz_a->setVal(0,0,dU2dxz);
        dU2dxzAnalytical[gid] = dU2dxz_a;
        delete[] Pv;
        
        int nadj = iee_pmap->i_map[gid].size();

        for(int j=0;j<nadj;j++)
        {
            int adjID = iee_pmap->i_map[gid][j];
            
            if(adjID>=Nel_glob)
            {
                int fid = ief_pmap->i_map[gid][j];
                int nvf = ifn_pmap->i_map[fid].size();
                double* Fa = new double[nvf*3];
                Vert* vface = new Vert;
                for(int k=0;k<nvf;k++)
                {
                    int gvid = ifn_pmap->i_map[fid][k];
                    int lvid = gV2lV[gvid];
                    vface->x = vface->x + locVs[lvid]->x;
                    vface->y = vface->y + locVs[lvid]->y;
                    vface->z = vface->z + locVs[lvid]->z;
                    
                    Vert* V = new Vert;
                    V->x    = locVs[lvid]->x;
                    V->y    = locVs[lvid]->y;
                    V->z    = locVs[lvid]->z;
                    face.push_back(V);
                }
                
                vface->x = vface->x/nvf;
                vface->y = vface->y/nvf;
                vface->z = vface->z/nvf;
                
                Vec3D* r0 = new Vec3D;
                r0->c0 = (vface->x-Vijk->x);
                r0->c1 = (vface->y-Vijk->y);
                r0->c2 = (vface->z-Vijk->z);
                
                v0->c0 = face[1]->x-face[0]->x;
                v0->c1 = face[1]->y-face[0]->y;
                v0->c2 = face[1]->z-face[0]->z;

                v1->c0 = face[2]->x-face[0]->x;
                v1->c1 = face[2]->y-face[0]->y;
                v1->c2 = face[2]->z-face[0]->z;
                
                n0 = ComputeSurfaceNormal(v0,v1);
                double orient0   = DotVec3D(r0,n0);
                
                if(orient0<0.0)
                {
                    NegateVec3D(n0);
                }
                
                rdotn = DotVec3D(r0,n0);
                
                Vec3D* reflect = new Vec3D;
                reflect->c0 = r0->c0-2.0*(rdotn)*n0->c0;
                reflect->c1 = r0->c1-2.0*(rdotn)*n0->c1;
                reflect->c2 = r0->c2-2.0*(rdotn)*n0->c2;
                
                double vgx = vface->x - reflect->c0;
                double vgy = vface->y - reflect->c1;
                double vgz = vface->z - reflect->c2;
                
                r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                U = 0.1*tanh(50*(r-0.5))+1.0;
                dUdx = 5*vgx/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                dUdy = 5*(vgy-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                dUdz = 5*(vgz-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                
                dU2dx2 = (5*(-100*vgx*vgx*r*tanh(50*(r-0.5)))+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5))/(cosh(50*(r-0.5))*cosh(50*(r-0.5))*pow(r,3.0/2.0));
                dU2dxy =  5*vgx*(vgy-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                dU2dxz = 5*vgx*(vmz-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                
                gbMap[adjID]=U;
                gbMap_dUdx[adjID]=dUdx;
                gbMap_dUdy[adjID]=dUdy;
                gbMap_dUdz[adjID]=dUdz;
                
                gbMap_dU2dx2[adjID]=dU2dx2;
                gbMap_dU2dxy[adjID]=dU2dxy;
                gbMap_dU2dxz[adjID]=dU2dxz;
                face.clear();
            }
            
            if(iee_pmap->i_map.find(adjID)!=iee_pmap->i_map.end())
            {
                int NadjadjID = iee_pmap->i_map[adjID].size();
                
                int NvPadjadjID = gE2lV[adjID].size();
                double* PadjadjID = new double[NvPadjadjID*3];
                
                for(int k=0;k<NvPadjadjID;k++)
                {
                    loc_vid     = gE2lV[adjID][k];
                    PadjadjID[k*3+0] = locVs[loc_vid]->x;
                    PadjadjID[k*3+1] = locVs[loc_vid]->y;
                    PadjadjID[k*3+2] = locVs[loc_vid]->z;
                }
                
                Vert* VadjadjID   = ComputeCentroidCoord(PadjadjID,NvPadjadjID);
                
                delete[] PadjadjID;
                
                for(int k=0;k<NadjadjID;k++)
                {
                    int adjadjID = iee_pmap->i_map[adjID][k];
                    
                    if(adjadjID>=Nel_glob)
                    {
                        int fid = ief_pmap->i_map[adjID][k];
                        int nvf = ifn_pmap->i_map[fid].size();
                        Vert* vface = new Vert;
                        for(int k=0;k<nvf;k++)
                        {
                            int gvid = ifn_pmap->i_map[fid][k];
                            int lvid = gV2lV[gvid];
                            vface->x = vface->x + locVs[lvid]->x;
                            vface->y = vface->y + locVs[lvid]->y;
                            vface->z = vface->z + locVs[lvid]->z;
                            
                            Vert* V = new Vert;
                            V->x    = locVs[lvid]->x;
                            V->y    = locVs[lvid]->y;
                            V->z    = locVs[lvid]->z;
                            face.push_back(V);
                        }
                        
                        vface->x = vface->x/nvf;
                        vface->y = vface->y/nvf;
                        vface->z = vface->z/nvf;
                        
                        Vec3D* r0 = new Vec3D;
                        r0->c0 = (vface->x-VadjadjID->x);
                        r0->c1 = (vface->y-VadjadjID->y);
                        r0->c2 = (vface->z-VadjadjID->z);
                        
                        v0->c0 = face[1]->x-face[0]->x;
                        v0->c1 = face[1]->y-face[0]->y;
                        v0->c2 = face[1]->z-face[0]->z;

                        v1->c0 = face[2]->x-face[0]->x;
                        v1->c1 = face[2]->y-face[0]->y;
                        v1->c2 = face[2]->z-face[0]->z;
                        
                        n0 = ComputeSurfaceNormal(v0,v1);
                        double orient0   = DotVec3D(r0,n0);
                        
                        if(orient0<0.0)
                        {
                            NegateVec3D(n0);
                        }
                        
                        rdotn = DotVec3D(r0,n0);
                        
                        Vec3D* reflect = new Vec3D;
                        reflect->c0 = r0->c0-2.0*(rdotn)*n0->c0;
                        reflect->c1 = r0->c1-2.0*(rdotn)*n0->c1;
                        reflect->c2 = r0->c2-2.0*(rdotn)*n0->c2;
                        
                        double vgx = vface->x - reflect->c0;
                        double vgy = vface->y - reflect->c1;
                        double vgz = vface->z - reflect->c2;
                        
                        r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                        U = 0.1*tanh(50*(r-0.5))+1.0;
                        dUdx = 5*vgx/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                        dUdy = 5*(vgy-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                        dUdz = 5*(vgz-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                        
                        dU2dx2 = (5*(-100*vgx*vgx*r*tanh(50*(r-0.5)))+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5))/(cosh(50*(r-0.5))*cosh(50*(r-0.5))*pow(r,3.0/2.0));
                        dU2dxy =  5*vgx*(vgy-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                        dU2dxz = 5*vgx*(vmz-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                        
                        gbMap[adjadjID]=U;
                        
                        gbMap_dUdx[adjadjID]=dUdx;
                        gbMap_dUdy[adjadjID]=dUdy;
                        gbMap_dUdz[adjadjID]=dUdz;
                        
                        gbMap_dU2dx2[adjadjID]=dU2dx2;
                        gbMap_dU2dxy[adjadjID]=dU2dxy;
                        gbMap_dU2dxz[adjadjID]=dU2dxz;
                        
                        face.clear();
                    }
                    
                    if(iee_pmap->i_map.find(adjadjID)!=iee_pmap->i_map.end())
                    {
                        int NadjadjadjID = iee_pmap->i_map[adjadjID].size();
//
                        int NvPadjadjadjID = gE2lV[adjadjID].size();
                        double* PadjadjadjID = new double[NvPadjadjadjID*3];

                        for(int k=0;k<NvPadjadjadjID;k++)
                        {
                            loc_vid     = gE2lV[adjadjID][k];
                            PadjadjadjID[k*3+0] = locVs[loc_vid]->x;
                            PadjadjadjID[k*3+1] = locVs[loc_vid]->y;
                            PadjadjadjID[k*3+2] = locVs[loc_vid]->z;
                        }

                        Vert* VadjadjadjID   = ComputeCentroidCoord(PadjadjadjID,NvPadjadjadjID);
//
                        delete[] PadjadjadjID;


                        for(int f=0;f<NadjadjadjID;f++)
                        {
                            int adjadjadjID = iee_pmap->i_map[adjadjID][f];

                            if(adjadjadjID>=Nel_glob)
                            {
                                int fid = ief_pmap->i_map[adjadjID][f];
                                int nvf = ifn_pmap->i_map[fid].size();
                                Vert* vface = new Vert;
                                for(int ve=0;ve<nvf;ve++)
                                {
                                    int gvid = ifn_pmap->i_map[fid][ve];
                                    int lvid = gV2lV[gvid];
                                    vface->x = vface->x + locVs[lvid]->x;
                                    vface->y = vface->y + locVs[lvid]->y;
                                    vface->z = vface->z + locVs[lvid]->z;

                                    Vert* V = new Vert;
                                    V->x    = locVs[lvid]->x;
                                    V->y    = locVs[lvid]->y;
                                    V->z    = locVs[lvid]->z;
                                    face.push_back(V);
                                }

                                vface->x = vface->x/nvf;
                                vface->y = vface->y/nvf;
                                vface->z = vface->z/nvf;

                                Vec3D* r0 = new Vec3D;
                                r0->c0 = (vface->x-VadjadjadjID->x);
                                r0->c1 = (vface->y-VadjadjadjID->y);
                                r0->c2 = (vface->z-VadjadjadjID->z);

                                v0->c0 = face[1]->x-face[0]->x;
                                v0->c1 = face[1]->y-face[0]->y;
                                v0->c2 = face[1]->z-face[0]->z;

                                v1->c0 = face[2]->x-face[0]->x;
                                v1->c1 = face[2]->y-face[0]->y;
                                v1->c2 = face[2]->z-face[0]->z;

                                n0 = ComputeSurfaceNormal(v0,v1);
                                double orient0   = DotVec3D(r0,n0);

                                if(orient0<0.0)
                                {
                                    NegateVec3D(n0);
                                }

                                rdotn = DotVec3D(r0,n0);

                                Vec3D* reflect = new Vec3D;
                                reflect->c0 = r0->c0-2.0*(rdotn)*n0->c0;
                                reflect->c1 = r0->c1-2.0*(rdotn)*n0->c1;
                                reflect->c2 = r0->c2-2.0*(rdotn)*n0->c2;

                                double vgx = vface->x - reflect->c0;
                                double vgy = vface->y - reflect->c1;
                                double vgz = vface->z - reflect->c2;

                                r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                                U = 0.1*tanh(50*(r-0.5))+1.0;
                                dUdx = 5*vgx/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                                dUdy = 5*(vgy-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                                dUdz = 5*(vgz-0.5)/(r*cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                                
                                dU2dx2 = (5*(-100*vgx*vgx*r*tanh(50*(r-0.5)))+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5))/(cosh(50*(r-0.5))*cosh(50*(r-0.5))*pow(r,3.0/2.0));
                                dU2dxy =  5*vgx*(vgy-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                                dU2dxz = 5*vgx*(vmz-0.5)*(-1.0/(pow(r*r,3.0/2.0))-100*tanh(50*(r-0.5))/(r*r))/(cosh(50*(r-0.5))*cosh(50*(r-0.5)));
                                
                                gbMap[adjadjadjID]=U;
                                
                                gbMap_dUdx[adjadjadjID]=dUdx;
                                gbMap_dUdy[adjadjadjID]=dUdy;
                                gbMap_dUdz[adjadjadjID]=dUdz;
                                
                                gbMap_dU2dx2[adjadjadjID]=dU2dx2;
                                gbMap_dU2dxy[adjadjadjID]=dU2dxy;
                                gbMap_dU2dxz[adjadjadjID]=dU2dxz;
                                face.clear();
                            }
                        }
                    }
                }
            }
        }
    }

    if(world_rank == 0)
    {
        std::cout << "Done Preparing Boundary Data..." << std::endl;
    }

    
    double sum      = 0.0;
    double sum_dist = 0.0;
    double di       = 0.0;
    double uval     = 0.0;
    double du       = 0.0;
    
    Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);

    P->AddStateVecForAdjacentElements(Uvaria_map,1,comm);
        
    std::map<int,Array<double>* > Mvar_vmap = P->ReduceStateVecToAllVertices(Uvaria_map,1);
    
    std::map<int,Array<double>* >::iterator vm;
    
    delete us3d->ghost;
    delete us3d->ien;
    delete us3d->iee;
    delete us3d->ief;
    delete us3d->ie_Nv;
    delete us3d->ie_Nf;
//  delete us3d->ifn;
    delete us3d->ife;
    //delete us3d->if_ref;
    delete us3d->if_Nv;
//  delete us3d->xcn;
    
    std::map<int,Array<double>* > dUdXi_HO  = ComputedUdx_LSQ_HO_US3D(P,Uvaria_map,meshTopo,gbMap,comm);
    //std::map<int,Array<double>* > dUdXi_LS  = ComputedUdx_LSQ_US3D_LargeStencil(P,Uvaria_map,meshTopo,gbMap,comm);
    std::map<int,Array<double>* > dUdXi     = ComputedUdx_LSQ_US3D(P,Uvaria_map,meshTopo,gbMap,comm);
    std::map<int,Array<double>* > dUdXi_mgg = ComputedUdx_MGG(P,Uvaria_map,meshTopo,gbMap,comm);
    std::map<int,Array<double>* > dUdx_lsqv1;
    std::map<int,Array<double>* >::iterator itsol;
    double L2_dudx_lsq_1   = 0.0;
    double L2_dudx_mgg     = 0.0;
    double L2_dudx_lsq_2   = 0.0;
    double L2_d2udx2_lsq_2 = 0.0;
    double L2_d2udxy_lsq_2 = 0.0;
    double L2_d2udxz_lsq_2 = 0.0;
    int Nell=0;
    for(itsol=dUdXi_HO.begin();itsol!=dUdXi_HO.end();itsol++)
    {
        int elId = itsol->first;
        
        Array<double>* entry=new Array<double>(1,1);
        entry->setVal(0,0,dUdXi[elId]->getVal(0,0));
        dUdx_lsqv1[elId] = entry;
        L2_dudx_lsq_1   = L2_dudx_lsq_1+(dUdXi[elId]->getVal(0,0)-dUdxAnalytical[elId]->getVal(0,0))*(dUdXi[elId]->getVal(0,0)-dUdxAnalytical[elId]->getVal(0,0));
        L2_dudx_mgg     = L2_dudx_mgg+(dUdXi_mgg[elId]->getVal(0,0)-dUdxAnalytical[elId]->getVal(0,0))*(dUdXi_mgg[elId]->getVal(0,0)-dUdxAnalytical[elId]->getVal(0,0));
        L2_dudx_lsq_2   = L2_dudx_lsq_2+(itsol->second->getVal(0,0)-dUdxAnalytical[elId]->getVal(0,0))*(itsol->second->getVal(0,0)-dUdxAnalytical[elId]->getVal(0,0));
        L2_d2udx2_lsq_2 = L2_d2udx2_lsq_2+(itsol->second->getVal(3,0)-dU2dx2Analytical[elId]->getVal(0,0))*(itsol->second->getVal(3,0)-dU2dx2Analytical[elId]->getVal(0,0));
        L2_d2udxy_lsq_2 = L2_d2udxy_lsq_2+(itsol->second->getVal(4,0)-dU2dxyAnalytical[elId]->getVal(0,0))*(itsol->second->getVal(4,0)-dU2dxyAnalytical[elId]->getVal(0,0));
        L2_d2udxz_lsq_2 = L2_d2udxz_lsq_2+(itsol->second->getVal(5,0)-dU2dxzAnalytical[elId]->getVal(0,0))*(itsol->second->getVal(5,0)-dU2dxzAnalytical[elId]->getVal(0,0));
        
        Nell++;
    }
    double L2_dudx_lsq_1_tot;
    MPI_Allreduce(&L2_dudx_lsq_1, &L2_dudx_lsq_1_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    double L2_dudx_mgg_tot;
    MPI_Allreduce(&L2_dudx_mgg, &L2_dudx_mgg_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    double L2_dudx_lsq_2_tot;
    MPI_Allreduce(&L2_dudx_lsq_2, &L2_dudx_lsq_2_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    double L2_d2udx2_lsq_2_tot;
    MPI_Allreduce(&L2_d2udx2_lsq_2, &L2_d2udx2_lsq_2_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    double L2_d2udxy_lsq_2_tot;
    MPI_Allreduce(&L2_d2udxy_lsq_2, &L2_d2udxy_lsq_2_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    double L2_d2udxz_lsq_2_tot;
    MPI_Allreduce(&L2_d2udxz_lsq_2, &L2_d2udxz_lsq_2_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    int Nell_tot;
    MPI_Allreduce(&Nell, &Nell_tot, 1, MPI_INT, MPI_SUM, comm);
    
    P->AddStateVecForAdjacentElements(dUdx_lsqv1,1,comm);
    std::map<int,Array<double>* > dU2dXi2   = ComputedUdx_LSQ_US3D(P,dUdx_lsqv1,meshTopo,gbMap_dUdx,comm);
    int Nell2 = 0;
    double L2_d2udx2_lsq_1 = 0.0;
    double L2_d2udxy_lsq_1 = 0.0;
    double L2_d2udxz_lsq_1 = 0.0;
    for(itsol=dU2dXi2.begin();itsol!=dU2dXi2.end();itsol++)
    {
        int elId = itsol->first;
        L2_d2udx2_lsq_1 = L2_d2udx2_lsq_1+(itsol->second->getVal(0,0)-dU2dx2Analytical[elId]->getVal(0,0))*(itsol->second->getVal(0,0)-dU2dx2Analytical[elId]->getVal(0,0));
        L2_d2udxy_lsq_1 = L2_d2udxy_lsq_1+(itsol->second->getVal(1,0)-dU2dxyAnalytical[elId]->getVal(0,0))*(itsol->second->getVal(1,0)-dU2dxyAnalytical[elId]->getVal(0,0));
        L2_d2udxz_lsq_1 = L2_d2udxz_lsq_1+(itsol->second->getVal(2,0)-dU2dxzAnalytical[elId]->getVal(0,0))*(itsol->second->getVal(2,0)-dU2dxzAnalytical[elId]->getVal(0,0));
        Nell2++;
    }
    double L2_d2udx2_lsq_1_tot;
    double L2_d2udxy_lsq_1_tot;
    double L2_d2udxz_lsq_1_tot;
    
    MPI_Allreduce(&L2_d2udx2_lsq_1, &L2_d2udx2_lsq_1_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&L2_d2udxy_lsq_1, &L2_d2udxy_lsq_1_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&L2_d2udxz_lsq_1, &L2_d2udxz_lsq_1_tot, 1, MPI_DOUBLE, MPI_SUM, comm);
    int Nell_tot2;
    MPI_Allreduce(&Nell2, &Nell_tot2, 1, MPI_INT, MPI_SUM, comm);
    
    std::vector<double> errors = ReadReferenceErrors("errors.ref");
    
    if(world_rank == 0)
    {
        double err0 = sqrt(L2_dudx_lsq_1_tot)/Nell_tot-errors[0];
        double err1 = sqrt(L2_dudx_lsq_2_tot)/Nell_tot-errors[1];
        double err2 = sqrt(L2_dudx_mgg_tot)/Nell_tot-errors[2];
        double err3 = sqrt(L2_d2udx2_lsq_1_tot)/Nell_tot2-errors[3];
        double err4 = sqrt(L2_d2udxy_lsq_1_tot)/Nell_tot2-errors[4];
        double err5 = sqrt(L2_d2udxz_lsq_1_tot)/Nell_tot2-errors[5];
        double err6 = sqrt(L2_d2udx2_lsq_2_tot)/Nell_tot-errors[6];
        double err7 = sqrt(L2_d2udxy_lsq_2_tot)/Nell_tot-errors[7];
        double err8 = sqrt(L2_d2udxz_lsq_2_tot)/Nell_tot-errors[8];
        
        std::cout << std::setprecision(16) << "L2norm dUdx   lsq v1   = " << sqrt(L2_dudx_lsq_1_tot)/Nell_tot << std::endl;
        std::cout << std::setprecision(16) << "L2norm dUdx   lsq v2   = " << sqrt(L2_dudx_lsq_2_tot)/Nell_tot << std::endl;
        std::cout << std::setprecision(16) << "L2norm dUdx   mgg      = " << sqrt(L2_dudx_mgg_tot)/Nell_tot << std::endl;
        std::cout << std::setprecision(16) << "L2norm dU2dx2 lsq v1   = " << sqrt(L2_d2udx2_lsq_1_tot)/Nell_tot2 << std::endl;
        std::cout << std::setprecision(16) << "L2norm dU2dxy lsq v1   = " << sqrt(L2_d2udxy_lsq_1_tot)/Nell_tot2 << std::endl;
        std::cout << std::setprecision(16) << "L2norm dU2dxz lsq v1   = " << sqrt(L2_d2udxz_lsq_1_tot)/Nell_tot2 << std::endl;
        std::cout << std::setprecision(16) << "L2norm dU2dx2 lsq v2   = " << sqrt(L2_d2udx2_lsq_2_tot)/Nell_tot << std::endl;
        std::cout << std::setprecision(16) << "L2norm dU2dxy lsq v2   = " << sqrt(L2_d2udxy_lsq_2_tot)/Nell_tot << std::endl;
        std::cout << std::setprecision(16) << "L2norm dU2dxz lsq v2   = " << sqrt(L2_d2udxz_lsq_2_tot)/Nell_tot << std::endl;
        
        double total_error = err0+err1+err2+err3+err4+err5+err6+err7+err8;
        std::cout << "Total error " << total_error << std::endl;
        
        if(fabs(total_error) < 1.0e-16)
        {
            std::cout << "The gradient reconstruction tests have PASSED!" << std::endl;
	    std::cout << "Partitioning test has PASSED." << std::endl;
            
            std::ofstream myfile;
            myfile.open("gradient_reconstruction_test.PASSED");
            myfile << "SUCCES!" << std::endl;
            myfile.close();
        }
	else
        {
            std::cout << "The gradient reconstruction tests have FAILED. " << total_error <<  std::endl;
            std::ofstream myfile;
            myfile.open("gradient_reconstruction.FAILED");
            myfile << "FAILED! error = " << total_error << std::endl;
            myfile.close();
        }
        
    }
    
    std::map<int,std::vector<int> >::iterator ienit;
    
    P->AddStateVecForAdjacentElements(dUdXi_HO,9,comm);
    std::map<int,Array<double>* > dudx_ho_vmap = P->ReduceStateVecToAllVertices(dUdXi_HO,9);
    
    P->AddStateVecForAdjacentElements(dUdXi_mgg,3,comm);
    std::map<int,Array<double>* > dudx_mgg_vmap = P->ReduceStateVecToAllVertices(dUdXi_mgg,3);
    
    P->AddStateVecForAdjacentElements(dUdXi,3,comm);
    std::map<int,Array<double>* > dudx_vmap = P->ReduceStateVecToAllVertices(dUdXi,3);
    
    P->AddStateVecForAdjacentElements(dU2dXi2,3,comm);
    std::map<int,Array<double>* > du2dx2_vmap = P->ReduceStateVecToAllVertices(dU2dXi2,3);

    P->AddStateVecForAdjacentElements(dUdxAnalytical,1,comm);
    std::map<int,Array<double>* > dudx_a_vmap = P->ReduceStateVecToAllVertices(dUdxAnalytical,1);
    
    P->AddStateVecForAdjacentElements(dU2dx2Analytical,1,comm);
    std::map<int,Array<double>* > du2dx2_a_vmap = P->ReduceStateVecToAllVertices(dU2dx2Analytical,1);
    
    std::set<int> gv_set;
    std::vector<int> lverts;
    std::map<int,int> lpartv2gv_v2;
    std::map<int,int> gv2lpv2;
    std::vector<std::vector<int> > locelem2locnode;
    std::vector<Vert*> printVs;
    int lcv2 = 0;
    std::vector<std::vector<double>> variable;
    
    for(int p=0;p<LocElem.size();p++)
    {
        int gEl = LocElem[p];
        int nv = ien_pmap->i_map[gEl].size();
        
        if(nv == 4) // tetrahedra
        {
            std::vector<int> nodes(4);
            for(int q=0;q<nv;q++)
            {
                int gv = ien_pmap->i_map[gEl][q];
                int lvv = gV2lV[gv];
                if(gv_set.find(gv)==gv_set.end())
                {
                    gv_set.insert(gv);
                    lverts.push_back(lvv);
                    lpartv2gv_v2[lvv]=gv;
                    gv2lpv2[gv]=lcv2;
                    nodes[q] = lcv2;
    
                    printVs.push_back(locVs[lvv]);
                    std::vector<double> comp(8);
                    comp[0] = Mvar_vmap[gv]->getVal(0,0);
                    comp[1] = dudx_vmap[gv]->getVal(0,0);
                    comp[2] = dudx_mgg_vmap[gv]->getVal(0,0);
                    comp[3] = dudx_a_vmap[gv]->getVal(0,0);
                    comp[4] = dudx_ho_vmap[gv]->getVal(0,0);
                    comp[5] = dudx_ho_vmap[gv]->getVal(3,0);
                    comp[6] = du2dx2_vmap[gv]->getVal(0,0);
                    comp[7] = du2dx2_a_vmap[gv]->getVal(0,0);
                    variable.push_back(comp);
                    
                    lcv2=lcv2+1;
                }
                else
                {
                    int lcv_u = gv2lpv2[gv];
                    nodes[q] = lcv_u;
                }
            }
            
            locelem2locnode.push_back(nodes);
        }
    }
    
    if(locelem2locnode.size()>0)
    {
        std::string filename = "checkPart_" + std::to_string(world_rank) + ".dat";
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"dUdx\", \"dUdx_mgg\", \"dUdxa\", \"dUdx_ho\", \"d2Udx2_ho\", \"dU2dx2_lo\", \"dU2dx2_a\"" << std::endl;
        myfile <<"ZONE N = " << printVs.size() << ", E = " << locelem2locnode.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

        for(int i=0;i<printVs.size();i++)
        {
            myfile << printVs[i]->x << " " << printVs[i]->y << " " << printVs[i]->z  << " " << variable[i][0] << " " << variable[i][1]<< " " << variable[i][2] << " " << variable[i][3] << " " << variable[i][4] << " " << variable[i][5] << " " << variable[i][6] << " " << variable[i][7] << std::endl;
        }
        int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
        int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
        for(int i=0;i<locelem2locnode.size();i++)
        {
            myfile <<   locelem2locnode[i][0]+1 << "  " <<
                        locelem2locnode[i][1]+1 << "  " <<
                        locelem2locnode[i][2]+1 << "  " <<
                        locelem2locnode[i][3]+1 << "  " << std::endl;
        }


        myfile.close();
    }
    
    
    
//    std::vector<int> lverts;
//    std::map<int,int> lpartv2gv_v2;
//    std::map<int,int> gv2lpv2;
//
//    std::set<int> gv_set;
//    int lcv2 = 0;
//    Array<int>* ien_part_tetra     = tmesh->ien_part_tetra;
//    Array<int>* ien_part_hybrid    = tmesh->ien_part_hybrid;
//    std::vector<Vert*> locVs       = tmesh->LocalVerts;
//    int nElonRank                  = ien_part_tetra->getNrow();
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
    
    
    MPI_Finalize();
    
}

