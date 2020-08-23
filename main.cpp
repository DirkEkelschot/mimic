#include "adapt_io.h"
#include "adapt_recongrad.h"
#include "adapt_recongrad2.h"
int mpi_size, mpi_rank;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int main(int argc, char** argv) {
    
    MPI_Init(NULL, NULL);
   
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    //======================================================================================
    //====================Parsing the inputs from the command line==========================
    //======================================================================================
    std::string fn_grid_t;
    std::string fn_conn_t;
    std::string fn_data_t;
    std::string destination;
    
    char buffer[20];
    
    const char* fn_grid;
    const char* fn_conn;
    const char* fn_data;
    for(int i = 1;i<argc;i++)
    {
        std::string str = argv[i];
        std::size_t l = str.copy(buffer,6,0);
        buffer[l]='\0';
        
        if (str.compare(0,6,"--grid") == 0)
        {
            char fn_grid_n[20];
            std::size_t ln = str.copy(fn_grid_n,str.size(),7);
            fn_grid_n[ln]='\0';
            fn_grid = fn_grid_n;
        }
        else if (str.compare(0,6,"--conn") == 0)
        {
            char fn_conn_n[20];
            std::size_t ln = str.copy(fn_conn_n,str.size(),7);
            fn_conn_n[ln]='\0';
            fn_conn = fn_conn_n;
        }
        else if (str.compare(0,6,"--data") == 0)
        {
            char fn_data_n[20];
            std::size_t ln = str.copy(fn_data_n,str.size(),7);
            fn_data_n[ln]='\0';
            fn_data = fn_data_n;
        }
        else
        {
            std::cout << "set correct path for " << buffer  << std::endl;
        }
    }
    //======================================================================================
    //======================================================================================
    //======================================================================================
    US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);
    
    int Nel_part = us3d->ien->getNrow();
    
    ParallelState* pstate = new ParallelState(us3d->ien->getNglob(),comm);

    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,comm,8);

    ParallelState* xcn_parstate = new ParallelState(us3d->xcn->getNglob(),comm);
    Array<double>* Uivar = new Array<double>(Nel_part,1);
    for(int i=0;i<Nel_part;i++)
    {
        Uivar->setVal(i,0,us3d->interior->getVal(i,0));
    }
    
    delete us3d->interior;
    
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, parmetis_pstate,pstate, us3d->xcn,xcn_parstate,Uivar,comm);

    //Array<double>* dUdXi_v2 = ComputedUdx(P, pstate, iee_copy, iee_loc, ief_loc, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm, ife_copy);
    
    std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Uivar,comm);
    
    Mesh_Topology* meshTopo = new Mesh_Topology(P,us3d->ifn,us3d->ghost,UauxNew,comm);
    
    Gradients* dudxObj = new Gradients(P,meshTopo,UauxNew,us3d->ghost,"us3d","MMG",comm);
    
    Array<double>* dUdXi = dudxObj->getdUdXi();
    
    
    Array<double>* dUidxi = new Array<double>(dUdXi->getNrow(),1);
    Array<double>* dUidyi = new Array<double>(dUdXi->getNrow(),1);
    Array<double>* dUidzi = new Array<double>(dUdXi->getNrow(),1);

    for(i=0;i<dUdXi->getNrow();i++)
    {
        dUidxi->setVal(i,0,dUdXi->getVal(i,0));
        dUidyi->setVal(i,0,dUdXi->getVal(i,1));
        dUidzi->setVal(i,0,dUdXi->getVal(i,2));
    }
    
    //std::map<int,std::vector<double> > dUdxiauxNew = P->CommunicateAdjacentDataUS3DNew(dUdXi,comm);
//    std::map<int,std::vector<double> > dUdxauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidxi,comm);
//    std::map<int,std::vector<double> > dUdyauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidyi,comm);
//    std::map<int,std::vector<double> > dUdzauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidzi,comm)
    ;
    std::map<int,double > dUdxauxNew  = P->CommunicateAdjacentDataUS3D(dUidxi,comm);
    std::map<int,double > dUdyauxNew  = P->CommunicateAdjacentDataUS3D(dUidyi,comm);
    std::map<int,double > dUdzauxNew  = P->CommunicateAdjacentDataUS3D(dUidzi,comm);
    
//    std::map<int,std::vector<double> >::iterator itmv;
//    int t=0;
//    for(itmv=dUdxiauxNew.begin();itmv!=dUdxiauxNew.end();itmv++)
//    {
//        std::cout << itmv->second[0]-dUdxauxNew[itmv->first][0] << " " <<                  itmv->second[1]-dUdyauxNew[itmv->first][0]  << " " << itmv->second[2]-dUdzauxNew[itmv->first][0]  << std::endl;
//        t++;
//    }
    delete dUdXi;
    
    Array<double>* dU2dXi2 = ComputedUdx_MGG(P,dUdxauxNew,meshTopo,us3d->ghost,comm);
    Array<double>* dU2dYi2 = ComputedUdx_MGG(P,dUdyauxNew,meshTopo,us3d->ghost,comm);
    Array<double>* dU2dZi2 = ComputedUdx_MGG(P,dUdzauxNew,meshTopo,us3d->ghost,comm);
    
    Array<double>* d2udx2 = new Array<double>(dU2dXi2->getNrow(),1);
    Array<double>* d2udxy = new Array<double>(dU2dXi2->getNrow(),1);
    Array<double>* d2udxz = new Array<double>(dU2dXi2->getNrow(),1);
    
    Array<double>* d2udyx = new Array<double>(dU2dYi2->getNrow(),1);
    Array<double>* d2udy2 = new Array<double>(dU2dYi2->getNrow(),1);
    Array<double>* d2udyz = new Array<double>(dU2dYi2->getNrow(),1);

    Array<double>* d2udzx = new Array<double>(dU2dZi2->getNrow(),1);
    Array<double>* d2udzy = new Array<double>(dU2dZi2->getNrow(),1);
    Array<double>* d2udz2 = new Array<double>(dU2dZi2->getNrow(),1);

    for(int i=0;i<dU2dXi2->getNrow();i++)
    {
        d2udx2->setVal(i,0,dU2dXi2->getVal(i,0));
        d2udxy->setVal(i,0,dU2dXi2->getVal(i,1));
        d2udxz->setVal(i,0,dU2dXi2->getVal(i,2));
        
        d2udyx->setVal(i,0,dU2dYi2->getVal(i,0));
        d2udy2->setVal(i,0,dU2dYi2->getVal(i,1));
        d2udyz->setVal(i,0,dU2dYi2->getVal(i,2));
        
        d2udzx->setVal(i,0,dU2dZi2->getVal(i,0));
        d2udzy->setVal(i,0,dU2dZi2->getVal(i,1));
        d2udz2->setVal(i,0,dU2dZi2->getVal(i,2));
    }
    
    std::vector<double> u_v = meshTopo->ReduceUToVertices(Uivar);
////
    std::vector<double> dudx_v = meshTopo->ReduceUToVertices(dUidxi);
    std::vector<double> dudy_v = meshTopo->ReduceUToVertices(dUidyi);
    std::vector<double> dudz_v = meshTopo->ReduceUToVertices(dUidzi);
    
    std::vector<double> d2udx2_v = meshTopo->ReduceUToVertices(d2udx2);
    std::vector<double> d2udxy_v = meshTopo->ReduceUToVertices(d2udxy);
    std::vector<double> d2udxz_v = meshTopo->ReduceUToVertices(d2udxz);
    
    std::vector<double> d2udyx_v = meshTopo->ReduceUToVertices(d2udyx);
    std::vector<double> d2udy2_v = meshTopo->ReduceUToVertices(d2udy2);
    std::vector<double> d2udyz_v = meshTopo->ReduceUToVertices(d2udyz);
    
    std::vector<double> d2udzx_v = meshTopo->ReduceUToVertices(d2udzx);
    std::vector<double> d2udzy_v = meshTopo->ReduceUToVertices(d2udzy);
    std::vector<double> d2udz2_v = meshTopo->ReduceUToVertices(d2udz2);

    std::vector<Vert> Verts =  P->getLocalVerts();
    int nVerts = Verts.size();
    Array<double>* hessian = new Array<double>(nVerts,9);
    Array<double>* grad    = new Array<double>(nVerts,3);
    for(i=0;i<nVerts;i++)
    {
        grad->setVal(i,0,dudx_v[i]);grad->setVal(i,1,dudy_v[i]);grad->setVal(i,2,dudz_v[i]);
        hessian->setVal(i,0,d2udx2_v[i]); hessian->setVal(i,1,d2udxy_v[i]); hessian->setVal(i,2,d2udxz_v[i]);
        hessian->setVal(i,3,d2udyx_v[i]); hessian->setVal(i,4,d2udy2_v[i]); hessian->setVal(i,5,d2udyz_v[i]);
        hessian->setVal(i,6,d2udzx_v[i]); hessian->setVal(i,7,d2udzy_v[i]); hessian->setVal(i,8,d2udz2_v[i]);
    }
    
    std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
        
    ComputeMetric(Verts,grad,hessian);
    
//    //=============================================================================
//    //=====================Output the data in Tecplot format=======================
//    //=============================================================================
    string filename = "gradients_hess_" + std::to_string(world_rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
    
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"u\", \"dx\", \"dy\", \"dz\", \"d2udx2\", \"d2udxy\", \"d2udxz\", \"d2udyx\", \"d2udy2\", \"d2udyz\", \"d2udzx\", \"d2udzy\", \"d2udz2\", \"M00\", \"M01\", \"M02\", \"M10\", \"M11\", \"M12\", \"M20\", \"M21\", \"M22\"" << std::endl;
    int nvert = Verts.size();
    myfile <<"ZONE N = " << nvert << ", E = " << us3d->ien->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;

    for(int i=0;i<nvert;i++)
    {
        myfile<< Verts[i].x << " " << Verts[i].y << " " << Verts[i].z <<
        " " << u_v[i] << " "<< dudx_v[i] << " " << dudy_v[i] << " " << dudz_v[i] <<
        " " << d2udx2_v[i] << " " << d2udxy_v[i] << " " << d2udxz_v[i] <<
        " " << d2udyx_v[i] << " " << d2udy2_v[i] << " " << d2udyz_v[i] <<
        " " << d2udzx_v[i] << " " << d2udzy_v[i] << " " << d2udz2_v[i] <<
        " " << hessian->getVal(i,0) << " " << hessian->getVal(i,1) << " " << hessian->getVal(i,2) <<
        " " << hessian->getVal(i,3) << " " << hessian->getVal(i,4) << " " << hessian->getVal(i,5) <<
        " " << hessian->getVal(i,6) << " " << hessian->getVal(i,7) << " " << hessian->getVal(i,8) << std::endl;
    }

    for(int i=0;i<us3d->ien->getNrow();i++)
    {
       myfile << loc_elem2verts_loc[i][0]+1 << "  " <<
                 loc_elem2verts_loc[i][1]+1 << "  " <<
                 loc_elem2verts_loc[i][2]+1 << "  " <<
                 loc_elem2verts_loc[i][3]+1 << "  " <<
                 loc_elem2verts_loc[i][4]+1 << "  " <<
                 loc_elem2verts_loc[i][5]+1 << "  " <<
                 loc_elem2verts_loc[i][6]+1 << "  " <<
                 loc_elem2verts_loc[i][7]+1 << std::endl;
    }
    myfile.close();
//

    delete d2udx2;
    delete d2udxy;
    delete d2udxz;

    delete d2udyx;
    delete d2udy2;
    delete d2udyz;

    delete d2udzx;
    delete d2udzy;
    delete d2udz2;


    d2udx2_v.erase(d2udx2_v.begin(),d2udx2_v.end());
    d2udxy_v.clear();
    d2udxz_v.clear();

    d2udyx_v.clear();
    d2udy2_v.clear();
    d2udyz_v.clear();

    d2udzx_v.clear();
    d2udzy_v.clear();
    d2udz2_v.clear();
    delete us3d->ifn;
    delete us3d->ien;
    delete us3d->iee;
    
    delete hessian;
    delete grad;
    delete parmetis_pstate;
    
    MPI_Finalize();
    
    return 0;
     
}
