#include "adapt_io.h"
#include "adapt_recongrad.h"
#include "adapt_recongrad2.h"
int mpi_size, mpi_rank;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

Array<double>* GetSolutionOnRoot(Partition* P, US3D* us3d, Array<double>* Hv, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j;
    
    ParallelState* xcn_pstate = P->getXcnParallelState();
    ParallelState* ien_pstate = P->getIenParallelState();
    
    int nlElem = us3d->ien->getNrow();
    int nElem = us3d->ien->getNglob();
    int nvg = us3d->xcn->getNglob();
    
    std::map<int,int> gV2lV = P->getGlobalVert2LocalVert();

    std::vector<double> Uvloc;
    std::vector<int> vids;
    std::vector<double> Uids;
    int gid,lid;
    int nval = 6;
    std::set<int> vdone;
    
    for(i=0;i<nlElem;i++)
    {
       for(j=0;j<8;j++)
       {
           gid = us3d->ien->getVal(i,j);
           lid = gV2lV[gid];
           int lrange = xcn_pstate->getOffsets()[world_rank];
           int hrange = xcn_pstate->getOffsets()[world_rank]+xcn_pstate->getNlocs()[world_rank];
           if(vdone.find(lid)==vdone.end() && (lrange<gid<hrange))
           {
               vdone.insert(lid);
               vids.push_back(gid);

               Uids.push_back(Hv->getVal(lid,0));
               Uids.push_back(Hv->getVal(lid,1));
               Uids.push_back(Hv->getVal(lid,2));
               Uids.push_back(Hv->getVal(lid,4));
               Uids.push_back(Hv->getVal(lid,5));
               Uids.push_back(Hv->getVal(lid,8));
           }
       }
    }
   
    int* vglob_nloc     = new int[world_size];
    int* vglob_nlocs    = new int[world_size];
    int* vglob_offsets  = new int[world_size];
   
    int* Uvglob_nloc    = new int[world_size];
    int* Uvglob_nlocs   = new int[world_size];
    int* Uvglob_offsets = new int[world_size];
   

    for(i=0;i<world_size;i++)
    {
       if(i == world_rank)
       {
          vglob_nloc[i] = vids.size();
       }
       else
       {
           vglob_nloc[i] = 0;
       }
    }
   
    MPI_Allreduce(vglob_nloc,  vglob_nlocs, world_size, MPI_INT, MPI_SUM, comm);
   
    int offset  = 0;
    int vglobt  = 0;
    int Uvglobt = 0;
    int Uoffset = 0;
    
    for(i=0;i<world_size;i++)
    {
        vglob_offsets[i]=offset;
        offset=offset+vglob_nlocs[i];
        vglobt = vglobt+vglob_nlocs[i];
       
        Uvglob_nlocs[i] = vglob_nlocs[i]*nval;
       
        Uvglob_offsets[i]=Uoffset;
        Uoffset=Uoffset+Uvglob_nlocs[i];
        Uvglobt = Uvglobt+Uvglob_nlocs[i];
       
    }
   
    int* vids_t = 0;
    double* Uvids_t = 0;
   
    if(world_rank == 0)
    {
       vids_t  = new int[vglobt];
       Uvids_t = new double[Uvglobt];
    }

    MPI_Gatherv(&vids[0],
                vids.size(),
                MPI_INT,
                &vids_t[0],
                vglob_nlocs,
                vglob_offsets,
                MPI_INT, 0, comm);
   
    MPI_Gatherv(&Uids[0],
                Uids.size(),
                MPI_DOUBLE,
                &Uvids_t[0],
                Uvglob_nlocs,
                Uvglob_offsets,
                MPI_DOUBLE, 0, comm);
   
    Array<double>* Ug;
    
    if(world_rank == 0)
    {
        Ug = new Array<double>(us3d->xcn->getNglob(),nval);
        int cid=0;
        for(i=0;i<world_size;i++)
        {
            for(j=0;j<vglob_nlocs[i];j++)
            {
                cid = vids_t[vglob_offsets[i]+j];
            
                for(int k=0;k<nval;k++)
                {
                    Ug->setVal(cid,k,Uvids_t[(vglob_offsets[i]+j)*nval+k]);
                }
            }
        }
    }
    else
    {
        Ug = new Array<double>(1,1);
    }
   
    Array<double>* xcn_g;
    Array<int>* ien_g;
   
    if(world_rank == 0)
    {
        xcn_g = new Array<double>(nvg,3);
        ien_g = new Array<int>(nElem,8);
    }
    else
    {
        xcn_g = new Array<double>(1,1);
        ien_g = new Array<int>(1,1);
    }
   
    int* ien_nlocs      = new int[world_size];
    int* ien_offsets    = new int[world_size];
    int* xcn_nlocs      = new int[world_size];
    int* xcn_offsets    = new int[world_size];
   
    for(int i=0;i<world_size;i++)
    {
        
        xcn_nlocs[i]   = xcn_pstate->getNlocs()[i]*3;
        xcn_offsets[i] = xcn_pstate->getOffsets()[i]*3;
       
        ien_nlocs[i]   = ien_pstate->getNlocs()[i]*8;
        ien_offsets[i] = ien_pstate->getOffsets()[i]*8;
    }
   
    MPI_Gatherv(&us3d->xcn->data[0],
                us3d->xcn->getNrow()*3,
                MPI_DOUBLE,
                &xcn_g->data[0],
                xcn_nlocs,
                xcn_offsets,
                MPI_DOUBLE, 0, comm);
   

    MPI_Gatherv(&us3d->ien->data[0],
                us3d->ien->getNrow()*8,
                MPI_INT,
                &ien_g->data[0],
                ien_nlocs,
                ien_offsets,
                MPI_INT, 0, comm);
    if(world_rank == 0)
    {
        string filename = "d2UdX2i_global.dat";
        ofstream myfile;
        myfile.open(filename);
        myfile << "TITLE=\"volume.tec\"" << std::endl;
       
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"gxx\", \"gxy\", \"gxz\", \"gyy\", \"gyz\", \"gzz\"" << std::endl;
        myfile <<"ZONE N = " << nvg << ", E = " << nElem << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
       
        
        for(int i=0;i<nvg;i++)
        {
            myfile << xcn_g->getVal(i,0) << " " << xcn_g->getVal(i,1) << " " << xcn_g->getVal(i,2)
            << " " << Ug->getVal(i,0) << " " << Ug->getVal(i,1) << " " << Ug->getVal(i,2)
            << " " << Ug->getVal(i,3) << " " << Ug->getVal(i,4) << " " << Ug->getVal(i,5) << std::endl;
        }

        for(int i=0;i<nElem;i++)
        {
            myfile << ien_g->getVal(i,0)+1 << " " <<
                      ien_g->getVal(i,1)+1 << " " <<
                      ien_g->getVal(i,2)+1 << " " <<
                      ien_g->getVal(i,3)+1 << " " <<
                      ien_g->getVal(i,4)+1 << " " <<
                      ien_g->getVal(i,5)+1 << " " <<
                      ien_g->getVal(i,6)+1 << " " <<
                      ien_g->getVal(i,7)+1 << std::endl;
        }
        myfile.close();
    }
    
    delete[] vglob_nloc;
    delete[] vglob_nlocs;
    delete[] vglob_offsets;
    
    delete[] Uvglob_nloc;
    delete[] Uvglob_nlocs;
    delete[] Uvglob_offsets;
    delete xcn_g;
    delete ien_g;
    
    return Ug;
}

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
    std::map<int,const char*> fnames;
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
            fnames[i]=fn_grid;
        }
        else if (str.compare(0,6,"--conn") == 0)
        {
            char fn_conn_n[20];
            std::size_t ln = str.copy(fn_conn_n,str.size(),7);
            fn_conn_n[ln]='\0';
            fn_conn = fn_conn_n;
            fnames[i]=fn_conn;
        }
        else if (str.compare(0,6,"--data") == 0)
        {
            char fn_data_n[20];
            std::size_t ln = str.copy(fn_data_n,str.size(),7);
            fn_data_n[ln]='\0';
            fn_data = fn_data_n;
            fnames[i]=fn_data;
        }
    }
    if(fnames.size()!=3)
    {
        if(world_rank == 0)
        {
           std::cout << "required command line inputs are : --grid=\"grid_name\".h5 --conn=\"conn_name\".h5 --data=\"data_name\".h5"  << std::endl;
        }
        
        MPI_Finalize();
        return 0;
    }
    else
    {
        //========================================================================
        //========================================================================
        //========================================================================
        US3D* us3d = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);
        
        int Nel_part = us3d->ien->getNrow();
        
        ParallelState* ien_pstate = new ParallelState(us3d->ien->getNglob(),comm);

        ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,comm,8);

        ParallelState* xcn_pstate = new ParallelState(us3d->xcn->getNglob(),comm);
        Array<double>* Uivar = new Array<double>(Nel_part,1);
        for(int i=0;i<Nel_part;i++)
        {
            Uivar->setVal(i,0,us3d->interior->getVal(i,0));
        }
        
        delete us3d->interior;
        
        Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, parmetis_pstate,ien_pstate, us3d->xcn,xcn_pstate,Uivar,comm);
        
        //Array<double>* dUdXi_v2 = ComputedUdx(P, pstate, iee_copy, iee_loc, ief_loc, ifn_copy, ief_copy, Nel, Uaux, ghost, bound, comm, ife_copy);
        
        std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Uivar,comm);
        
        Mesh_Topology* meshTopo = new Mesh_Topology(P,us3d->ifn,us3d->ghost,UauxNew,comm);
        
        //Gradients* dudxObj = new Gradients(P,meshTopo,UauxNew,us3d->ghost,"us3d","MMG",comm);
        
        //Array<double>* dUdXi = dudxObj->getdUdXi();
        //clock_t t;
        double tmax = 0.0;
        double tn = 0.0;
        //t = clock();
        Array<double>* dUdXi   = ComputedUdx_MGG(P,UauxNew,meshTopo,us3d->ghost,comm);
//        Gradients* dudxObj   = new Gradients(P,meshTopo,UauxNew,us3d->ghost,"us3d","LSQ",comm);
//        Array<double>* dUdXi = dudxObj->getdUdXi();
        
        //Array<double>* dUdXi = dudxObj->getdUdXi();

        Array<double>* dUidxi = new Array<double>(dUdXi->getNrow(),1);
        Array<double>* dUidyi = new Array<double>(dUdXi->getNrow(),1);
        Array<double>* dUidzi = new Array<double>(dUdXi->getNrow(),1);

        for(int i=0;i<dUdXi->getNrow();i++)
        {
            dUidxi->setVal(i,0,dUdXi->getVal(i,0));
            dUidyi->setVal(i,0,dUdXi->getVal(i,1));
            dUidzi->setVal(i,0,dUdXi->getVal(i,2));
        }
        
        //std::map<int,std::vector<double> > dUdxiauxNew = P->CommunicateAdjacentDataUS3DNew(dUdXi,comm);
    //    std::map<int,std::vector<double> > dUdxauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidxi,comm);
    //    std::map<int,std::vector<double> > dUdyauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidyi,comm);
    //    std::map<int,std::vector<double> > dUdzauxNew  = P->CommunicateAdjacentDataUS3DNew(dUidzi,comm);
        
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
        
        std::vector<double> u_v    = meshTopo->ReduceUToVertices(Uivar);
        
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
        for(int i=0;i<nVerts;i++)
        {
            grad->setVal(i,0,dudx_v[i]);grad->setVal(i,1,dudy_v[i]);grad->setVal(i,2,dudz_v[i]);
            hessian->setVal(i,0,d2udx2_v[i]); hessian->setVal(i,1,d2udxy_v[i]); hessian->setVal(i,2,d2udxz_v[i]);
            hessian->setVal(i,3,d2udyx_v[i]); hessian->setVal(i,4,d2udy2_v[i]); hessian->setVal(i,5,d2udyz_v[i]);
            hessian->setVal(i,6,d2udzx_v[i]); hessian->setVal(i,7,d2udzy_v[i]); hessian->setVal(i,8,d2udz2_v[i]);
        }
        
        Array<double>* UgRoot = GetSolutionOnRoot(P,us3d,hessian,comm);
        
        std::vector<std::vector<int> > loc_elem2verts_loc = P->getLocalElem2LocalVert();
            
        ComputeMetric(Verts,grad,hessian);
        
        //=================================================================
        //==================Output the data in Tecplot format==============
        //=================================================================
        
        string filename = "d2UdX2i_" + std::to_string(world_rank) + ".dat";
        ofstream myfile;
        myfile.open(filename);
        myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
        
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"u\", \"dx\", \"dy\", \"dz\", \"d2udx2\", \"d2udxy\", \"d2udxz\", \"d2udyx\", \"d2udy2\", \"d2udyz\", \"d2udzx\", \"d2udzy\", \"d2udz2\", \"M00\", \"M01\", \"M02\", \"M10\", \"M11\", \"M12\", \"M20\", \"M21\", \"M22\"" << std::endl;
        int nvert = Verts.size();
        myfile <<"ZONE N = " << nvert << ", E = " << us3d->ien->getNrow() << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
        
        std::map<int,int> gV2lV = P->getGlobalVert2LocalVert();
        
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
           myfile << loc_elem2verts_loc[i][0]+1 << " " <<
                     loc_elem2verts_loc[i][1]+1 << " " <<
                     loc_elem2verts_loc[i][2]+1 << " " <<
                     loc_elem2verts_loc[i][3]+1 << " " <<
                     loc_elem2verts_loc[i][4]+1 << " " <<
                     loc_elem2verts_loc[i][5]+1 << " " <<
                     loc_elem2verts_loc[i][6]+1 << " " <<
                     loc_elem2verts_loc[i][7]+1 << std::endl;
        }
        myfile.close();
         
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
        delete ien_pstate;
        delete parmetis_pstate;
        delete UgRoot;
        MPI_Finalize();
        
        return 0;
    }
}
