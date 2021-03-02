#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_parmetisstate2.h"
#include <iomanip>



std::vector<std::vector<double> > ReadRefMetricData(const char* filename)
{
    std::ifstream fin;
    fin.open(filename);

    // Read the file row by row
    std::vector<double> row(9);
    std::vector<std::vector<double> > m_ref;
    int t=0;
    while(fin >> row[0] >> row[1] >> row[2] >> row[3] >> row[4] >> row[5])
    {
       m_ref.push_back(row);
       t++;
    }
    
    return m_ref;
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
    
    const char* fn_grid="../test_mesh/cylinder_hex/grid.h5";
    const char* fn_conn="../test_mesh/cylinder_hex/conn.h5";
    const char* fn_data="../test_mesh/cylinder_hex/data.h5";
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,comm,info);
    const char* fn_metric = "metric.inp";
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
        
    int Nel_part = us3d->ien->getNrow();

    Array<double>* Ui = new Array<double>(Nel_part,1);
    int varia = 4;
    for(int i=0;i<Nel_part;i++)
    {
        Ui->setVal(i,0,us3d->interior->getVal(i,varia));
    }
    
    
    ParallelState* ien_pstate                   = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate                   = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState* xcn_pstate                   = new ParallelState(us3d->xcn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,comm);
//
    clock_t t;
    double tn = 0.0;
    t = clock();
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, comm);
//
    double duration = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Ptime = 0.0;
    MPI_Allreduce(&duration, &Ptime, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
    std::cout << "Timing partitioning: " << duration << std::endl;
    }

    std::vector<int> LocElem    = P->getLocElem();
    std::vector<double> Uvaria  = P->getLocElemVaria();

    std::map<int,double> Ui_map;
    double UvariaV = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid     = LocElem[i];
        UvariaV     = Uvaria[i];
        Ui_map[gid] = UvariaV;

    }

    Mesh_Topology* meshTopo = new Mesh_Topology(P,comm);

    Domain* pDom = P->getPartitionDomain();

    std::map<int,double> u_v     = P->ReduceFieldToVertices(Ui_map);

    std::vector<Vert> Verts      = P->getLocalVerts();

    std::map<int,double> UauxNew = P->CommunicateStateAdjacentElements(Ui_map,comm);

    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }

    t = clock();
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,UauxNew,gB,comm);
    double Gtiming = ( std::clock() - t) / (double) CLOCKS_PER_SEC;
    double Gmax_time = 0.0;
    MPI_Allreduce(&Gtiming, &Gmax_time, 1, MPI_DOUBLE, MPI_MAX, comm);
    if(world_rank == 0)
    {
        std::cout << "Timing gradient reconstruction... " << Gmax_time << std::endl;
    }

    Array<int>* lE2gE     = new Array<int>(dUdXi.size(),1);
    Array<double>* dUidxi = new Array<double>(dUdXi.size(),1);
    Array<double>* dUidyi = new Array<double>(dUdXi.size(),1);
    Array<double>* dUidzi = new Array<double>(dUdXi.size(),1);

    std::map<int,double> dUidxi_map;
    std::map<int,double> dUidyi_map;
    std::map<int,double> dUidzi_map;

    std::map<int,Array<double>* >::iterator grit;
    int ti=0;
    for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
    {
        lE2gE->setVal(ti,0,grit->first);
        dUidxi->setVal(ti,0,grit->second->getVal(0,0));
        dUidyi->setVal(ti,0,grit->second->getVal(1,0));
        dUidzi->setVal(ti,0,grit->second->getVal(2,0));
        dUidxi_map[grit->first]=grit->second->getVal(0,0);
        dUidyi_map[grit->first]=grit->second->getVal(1,0);
        dUidzi_map[grit->first]=grit->second->getVal(2,0);

        ti++;
    }
    std::map<int,double > dUdxauxNew  = P->CommunicateStateAdjacentElements(dUidxi_map,comm);
    std::map<int,double > dUdyauxNew  = P->CommunicateStateAdjacentElements(dUidyi_map,comm);
    std::map<int,double > dUdzauxNew  = P->CommunicateStateAdjacentElements(dUidzi_map,comm);

    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_US3D(P,dUdxauxNew,gB,comm);
    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_US3D(P,dUdyauxNew,gB,comm);
    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_US3D(P,dUdzauxNew,gB,comm);

    std::map<int,double> d2udx2_map,d2udxy_map,d2udxz_map,
                         d2udyx_map,d2udy2_map,d2udyz_map,
                         d2udzx_map,d2udzy_map,d2udz2_map;

    std::map<int,Array<double>*> Hess_map;

    std::map<int,Array<double>* >::iterator itgg;
    int te = 0;

    for(itgg=dU2dXi2.begin();itgg!=dU2dXi2.end();itgg++)
    {
        int gid = itgg->first;
        Array<double>* Hess = new Array<double>(3,3);

        Hess->setVal(0,0,dU2dXi2[gid]->getVal(0,0));
        Hess->setVal(0,1,dU2dXi2[gid]->getVal(1,0));
        Hess->setVal(0,2,dU2dXi2[gid]->getVal(2,0));

        Hess->setVal(1,0,dU2dYi2[gid]->getVal(0,0));
        Hess->setVal(1,1,dU2dYi2[gid]->getVal(1,0));
        Hess->setVal(1,2,dU2dYi2[gid]->getVal(2,0));

        Hess->setVal(2,0,dU2dZi2[gid]->getVal(0,0));
        Hess->setVal(2,1,dU2dZi2[gid]->getVal(1,0));
        Hess->setVal(2,2,dU2dZi2[gid]->getVal(2,0));

        Hess_map[gid] = Hess;
        t++;
    }

    std::map<int,Array<double>* > metric_e = ComputeMetric(P,metric_inputs, Hess_map, comm);

    Array<double>* lM2gM     = new Array<double>(metric_e.size(),6);

    for(grit=metric_e.begin();grit!=metric_e.end();grit++)
    {
        lM2gM->setVal(i,0,grit->second->getVal(0,0));
        lM2gM->setVal(i,1,grit->second->getVal(0,1));
        lM2gM->setVal(i,2,grit->second->getVal(0,2));
        lM2gM->setVal(i,3,grit->second->getVal(1,1));
        lM2gM->setVal(i,4,grit->second->getVal(1,2));
        lM2gM->setVal(i,5,grit->second->getVal(2,2));

        i++;
    }

    int nlElem = us3d->ien->getNrow();
    int nElem  = us3d->ien->getNglob();
    int nvg    = us3d->xcn->getNglob();


    Array<int>*     lE2gE_g;
    Array<double>*  lM2gM_g;
    Array<double>*  M_g;

    int nEl_glob = us3d->ien->getNglob();

    if(world_rank == 0)
    {
        lE2gE_g     = new Array<int>(nEl_glob,1);
        lM2gM_g     = new Array<double>(nEl_glob,6);
        M_g         = new Array<double>(nEl_glob,6);
    }
    else
    {
        lE2gE_g     = new Array<int>(1,1);
        lM2gM_g     = new Array<double>(1,1);
        M_g         = new Array<double>(1,1);
    }

    int* G_nlocs      = new int[world_size];
    int* red_G_nlocs  = new int[world_size];
    int* G_offsets    = new int[world_size];

    int* GM_nlocs      = new int[world_size];
    int* red_GM_nlocs  = new int[world_size];
    int* GM_offsets    = new int[world_size];

    for(i=0;i<world_size;i++)
    {
        G_nlocs[i] = 0;
        GM_nlocs[i] = 0;

        if(i==world_rank)
        {
            G_nlocs[i] = lM2gM->getNrow();
            GM_nlocs[i] = lM2gM->getNrow()*6;
        }
        else
        {
            G_nlocs[i] = 0;
            GM_nlocs[i] = 0;

        }
    }

    MPI_Allreduce(G_nlocs, red_G_nlocs, world_size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(GM_nlocs, red_GM_nlocs, world_size, MPI_INT, MPI_SUM, comm);

    int offset  = 0;
    int offsetM = 0;

    for(i=0;i<world_size;i++)
    {
        G_offsets[i]  = offset;
        offset        = offset+red_G_nlocs[i];

        GM_offsets[i] = offsetM;
        offsetM       = offsetM+red_GM_nlocs[i];
    }

    MPI_Gatherv(&lE2gE->data[0],
                lE2gE->getNrow(),
                MPI_INT,
                &lE2gE_g->data[0],
                red_G_nlocs,
                G_offsets,
                MPI_INT, 0, comm);

    MPI_Gatherv(&lM2gM->data[0],
                lM2gM->getNrow()*6,
                MPI_DOUBLE,
                &lM2gM_g->data[0],
                red_GM_nlocs,
                GM_offsets,
                MPI_DOUBLE, 0, comm);

    if(world_rank == 0)
    {
        std::vector<std::vector<double> > mv_ref = ReadRefMetricData("../test_mesh/metric_E_ref.dat");

        for(int i=0;i<nEl_glob;i++)
        {
            int gid = lE2gE_g->getVal(i,0);
            M_g->setVal(gid,0,lM2gM_g->getVal(i,0));
            M_g->setVal(gid,1,lM2gM_g->getVal(i,1));
            M_g->setVal(gid,2,lM2gM_g->getVal(i,2));
            M_g->setVal(gid,3,lM2gM_g->getVal(i,3));
            M_g->setVal(gid,4,lM2gM_g->getVal(i,4));
            M_g->setVal(gid,5,lM2gM_g->getVal(i,5));
        }
//        string filename = "metric_E.dat";
//        ofstream myfile;
//        myfile.open(filename);
//
//        for(int i=0;i<nEl_glob;i++)
//        {
//            myfile <<  std::setprecision(16) << M_g->getVal(i,0) << " " <<
//                        M_g->getVal(i,1) << " " <<
//                        M_g->getVal(i,2) << " " <<
//                        M_g->getVal(i,3) << " " <<
//                        M_g->getVal(i,4) << " " <<
//                        M_g->getVal(i,5) << std::endl;
//        }
//        myfile.close();

        int flip = 0;
        double err = 1.0e-08;
        for(int i=0;i<nEl_glob;i++)
        {
            for(int j=0;j<6;j++)
            {
                double diff = fabs(M_g->getVal(i,j)-mv_ref[i][j]);

                if(diff>err)
                {
                    std::cout << std::setprecision(16) << i << " " << diff << " " << M_g->getVal(i,j) << " " << mv_ref[i][j] << std::endl;
                    flip = 1;
                }
            }
        }

        if(flip == 1)
        {
            std::cout << " --::-- Parallel metric reconstruction test has FAILED. --::-- " << std::endl;
        }
        if(flip == 0)
        {
            std::cout << " --::-- Parallel metric reconstruction test has PASSED. --::-- " << std::endl;
        }
    }



    if(world_size == 4)
    {
        std::map<int,Array<double>* > Hess_vmap = P->ReduceMetricToVertices(Hess_map);
        std::map<int,Array<double>* > metric_v  = ComputeMetric(P,metric_inputs, Hess_vmap, comm);

        Array<double>* mv_g = GetOptimizedMMG3DMeshOnRoot(P, us3d, metric_v, comm);

        if(world_rank == 0)
        {
            std::vector<std::vector<double> > mv_ref = ReadRefMetricData("../test_mesh/metric_ref.dat");
            int flip = 0;
            int cnt = 0;
            double err = 1.0e-08;

            for(int i=0;i<mv_g->getNrow();i++)
            {
                for(int j=0;j<mv_g->getNcol();j++)
                {
                    double diff = fabs(mv_g->getVal(i,j)-mv_ref[i][j]);

                    if(diff>err)
                    {
                        std::cout << std::setprecision(16)<< "("<<i<<", "<<j<<") -> error = " << diff << " computed " << mv_g->getVal(i,j) << " ref " << mv_ref[i][j] << std::endl;
                        cnt++;
                        flip = 1;
                    }
                }
            }

            if(flip == 1)
            {
                std::cout << " --::-- Parallel metric computation+reduction to the vertices test has FAILED. --::-- " << cnt << " entries failed..." << std::endl;
            }
            if(flip == 0)
            {
                std::cout << " --::-- Parallel metric computation+reduction to the vertices test has PASSED. --::-- " << std::endl;
            }
        }
    }
//
//
//
//
//
    MPI_Finalize();
    
}
