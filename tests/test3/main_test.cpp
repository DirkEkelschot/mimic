#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include <iomanip>



std::vector<std::vector<double> > ReadRefMetricData()
{
    std::ifstream fin;
    fin.open("../test_mesh/metric_ref.dat");

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
    
    const char* fn_grid="../test_mesh/grid.h5";
    const char* fn_conn="../test_mesh/conn.h5";
    const char* fn_data="../test_mesh/data.h5";
    
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
    
    
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,comm,8);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief,
                                 us3d->ifn, us3d->ife, us3d->if_ref,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, comm);
    
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

    std::vector<double> u_v      = meshTopo->ReduceUToVertices(pDom,Ui_map);
    
    std::vector<Vert> Verts      = P->getLocalVerts();
    
    std::map<int,double> UauxNew = P->CommunicateAdjacentDataUS3D(Ui_map,comm);

    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    }
    
    t = clock();
    std::map<int,Array<double>* > dUdXi = ComputedUdx_LSQ_US3D(P,UauxNew,meshTopo,gB,comm);
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
        lE2gE->setVal(i,0,grit->first);
        dUidxi->setVal(i,0,grit->second->getVal(0,0));
        dUidyi->setVal(i,0,grit->second->getVal(1,0));
        dUidzi->setVal(i,0,grit->second->getVal(2,0));
        dUidxi_map[grit->first]=grit->second->getVal(0,0);
        dUidyi_map[grit->first]=grit->second->getVal(1,0);
        dUidzi_map[grit->first]=grit->second->getVal(2,0);
        
        ti++;
    }
    std::map<int,double > dUdxauxNew  = P->CommunicateAdjacentDataUS3D(dUidxi_map,comm);
    std::map<int,double > dUdyauxNew  = P->CommunicateAdjacentDataUS3D(dUidyi_map,comm);
    std::map<int,double > dUdzauxNew  = P->CommunicateAdjacentDataUS3D(dUidzi_map,comm);
    
    std::map<int,Array<double>* > dU2dXi2 = ComputedUdx_LSQ_US3D(P,dUdxauxNew,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dYi2 = ComputedUdx_LSQ_US3D(P,dUdyauxNew,meshTopo,gB,comm);
    std::map<int,Array<double>* > dU2dZi2 = ComputedUdx_LSQ_US3D(P,dUdzauxNew,meshTopo,gB,comm);
            
    std::map<int,double> d2udx2_map,d2udxy_map,d2udxz_map,
                         d2udyx_map,d2udy2_map,d2udyz_map,
                         d2udzx_map,d2udzy_map,d2udz2_map;
    
    std::map<int,Array<double>* >::iterator itgg;
    int te = 0;

    for(itgg=dU2dXi2.begin();itgg!=dU2dXi2.end();itgg++)
    {
        int gid = itgg->first;
        
        d2udx2_map[gid] = dU2dXi2[gid]->getVal(0,0);
        d2udxy_map[gid] = dU2dXi2[gid]->getVal(1,0);
        d2udxz_map[gid] = dU2dXi2[gid]->getVal(2,0);
        
        d2udyx_map[gid] = dU2dYi2[gid]->getVal(0,0);
        d2udy2_map[gid] = dU2dYi2[gid]->getVal(1,0);
        d2udyz_map[gid] = dU2dYi2[gid]->getVal(2,0);
        
        d2udzx_map[gid] = dU2dZi2[gid]->getVal(0,0);
        d2udzy_map[gid] = dU2dZi2[gid]->getVal(1,0);
        d2udz2_map[gid] = dU2dZi2[gid]->getVal(2,0);
        
        t++;
    }
    std::map<int,double> d2udx2_vm = meshTopo->ReduceFieldToVertices(pDom,d2udx2_map);
    std::map<int,double> d2udxy_vm = meshTopo->ReduceFieldToVertices(pDom,d2udxy_map);
    std::map<int,double> d2udxz_vm = meshTopo->ReduceFieldToVertices(pDom,d2udxz_map);

    std::map<int,double> d2udyx_vm = meshTopo->ReduceFieldToVertices(pDom,d2udyx_map);
    std::map<int,double> d2udy2_vm = meshTopo->ReduceFieldToVertices(pDom,d2udy2_map);
    std::map<int,double> d2udyz_vm = meshTopo->ReduceFieldToVertices(pDom,d2udyz_map);
    
    std::map<int,double> d2udzx_vm = meshTopo->ReduceFieldToVertices(pDom,d2udzx_map);
    std::map<int,double> d2udzy_vm = meshTopo->ReduceFieldToVertices(pDom,d2udzy_map);
    std::map<int,double> d2udz2_vm = meshTopo->ReduceFieldToVertices(pDom,d2udz2_map);
    
    std::map<int,Array<double>*> metric = ComputeMetric(P,metric_inputs,
                                                          d2udx2_vm,d2udxy_vm,d2udxz_vm,
                                                          d2udyx_vm,d2udy2_vm,d2udyz_vm,
                                                          d2udzx_vm,d2udzy_vm,d2udz2_vm,comm);

    Array<double>* mv_g = GetOptimizedMMG3DMeshOnRoot(P, us3d, metric, comm);
    
    if(world_rank == 0)
    {
        std::vector<std::vector<double> > mv_ref = ReadRefMetricData();
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
                    //std::cout << std::setprecision(16)<< "("<<i<<", "<<j<<") -> error = " << diff << " computed " << mv_g->getVal(i,j) << " ref " << mv_ref[i][j] << std::endl;
                    cnt++;
                    flip = 1;
                }
            }
        }

        if(flip == 1)
        {
            std::cout << " --::-- Parallel metric computation test has FAILED. --::-- " << cnt << " entries failed..." << std::endl;
        }
        if(flip == 0)
        {
            std::cout << " --::-- Parallel metric computation test has PASSED. --::-- " << std::endl;
        }
    }
    
    
    
    MPI_Finalize();
    
}
