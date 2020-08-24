#include "adapt_partition.h"
#include "adapt_math.h"
#include "adapt_compute.h"
#include "adapt_topology.h"

#ifndef ADAPT_RECONGRAD2_H
#define ADAPT_RECONGRAD2_H

class Gradients {
    public:
        Gradients(){};
        Gradients(Partition* Pa, Mesh_Topology* meshTopo, std::map<int,double> U, Array<double>* ghost, const char* solver_name, const char* recon_type, MPI_Comm comm);
        Array<double>* ComputedUdx_MGG(Partition* Pa, std::map<int,double> U,
                        Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm);
        Array<double>* ComputedUdx_LSQ_US3D_v2(Partition* Pa, std::map<int,double> U,
                        Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm);
        Array<double>* getdUdXi();
    private:
        //stored inputs
        std::vector<double> Un;
        Partition* P;
        MPI_Comm c;
    
        // outputs;
        Array<double>* dudx;
};

inline Gradients::Gradients(Partition* Pa, Mesh_Topology* meshTopo, std::map<int,double> U, Array<double>* ghost, const char* solver_name, const char* recon_type, MPI_Comm comm)
{
    if (strcmp(recon_type, "MGG") == 0)
    {
        dudx = ComputedUdx_MGG(Pa,U,meshTopo,ghost,comm);
    }
    if (strcmp(recon_type, "LSQ") == 0)
    {
        
        dudx = ComputedUdx_LSQ_US3D_v2(Pa,U,meshTopo,ghost,comm);
        //std::cout << "We still need to copy the correct routine for the LSQ grad recon." << std::endl;
    }
}



inline Array<double>* Gradients::ComputedUdx_MGG(Partition* Pa, std::map<int,double> U,
Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm)
{

        int lid, gEl, adjID, l_adjid, size, rank;
    double u_c, u_nb, gu_c_vx, gu_c_vy, gu_c_vz, gu_nb_vx, gu_nb_vy, gu_nb_vz,sum_phix,sum_phiy,sum_phiz,dphi_dn,Vol;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    int Nel = Pa->getLocalPartition()->getNglob();
    
    std::map<int,int> gE2lE                 = Pa->getGlobalElement2LocalElement();
    std::vector<int> Loc_Elem               = Pa->getLocElem();
    int nLoc_Elem                           = Loc_Elem.size();
    
    Array<double>* gu_c_x      = new Array<double>(nLoc_Elem,1);
    Array<double>* gu_c_y      = new Array<double>(nLoc_Elem,1);
    Array<double>* gu_c_z      = new Array<double>(nLoc_Elem,1);
    std::map<int,std::vector<int> > iee_vec = Pa->getIEEpartmap();
    dudx    = new Array<double>(nLoc_Elem,3);
    
    std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
    
    for(int i=0;i<nLoc_Elem;i++)
    {
        gu_c_x->setVal(i,0,0.0);
        gu_c_y->setVal(i,0,0.0);
        gu_c_z->setVal(i,0,0.0);
        
        dudx->setVal(i,0,0.0);
        dudx->setVal(i,1,0.0);
        dudx->setVal(i,2,0.0);
    }
    
    std::map<int,vector<Vec3D*> > normals   = meshTopo->getNormals();
    std::map<int,vector<Vec3D*> > rvector   = meshTopo->getRvectors();
    std::map<int,vector<Vec3D*> > dxfxc     = meshTopo->getdXfXc();
    std::map<int,vector<double> > dS        = meshTopo->getdS();
    std::map<int,vector<double> > dr        = meshTopo->getdr();
    std::map<int,double > vol               = meshTopo->getVol();

    int it = 0;
    double alpha   = 0.0;
    double L2normx = 0.0;
    double L2normy = 0.0;
    double L2normz = 0.0;
    double L2normx_max = 0.0;
    double L2normy_max = 0.0;
    double L2normz_max = 0.0;
    std::vector<Vec3D*> n_grads;
    clock_t t;
    
    for(int it=0;it<100;it++)
    {
        t = clock();
        //communicate grad phi!!!
        
        std::map<int,double> dUdx_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_x, comm);
        std::map<int,double> dUdy_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_y, comm);
        std::map<int,double> dUdz_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_z, comm);
        
        //std::map<int,std::vector<double> > dUdxi_p_bnd = Pa->CommunicateAdjacentDataUS3DNew(gu_c_old, comm);
                
        L2normx = 0.0;
        L2normy = 0.0;
        L2normz = 0.0;
        
        for(int i=0;i<nLoc_Elem;i++)
         {
             gEl = Loc_Elem[i];
             lid = gE2lE[gEl];
             u_c = U[gEl];
             
             gu_c_vx = gu_c_x->getVal(lid,0);
             gu_c_vy = gu_c_y->getVal(lid,0);
             gu_c_vz = gu_c_z->getVal(lid,0);

             sum_phix = 0.0;
             sum_phiy = 0.0;
             sum_phiz = 0.0;
             
             for(int j=0;j<6;j++)
             {
                 adjID   = iee_vec[gEl][j];
                 
                 if(adjID<Nel)
                 {
                     //l_adjid = gE2lE[adjID];
                     gu_nb_vx = dUdx_p_bnd[adjID];
                     gu_nb_vy = dUdy_p_bnd[adjID];
                     gu_nb_vz = dUdz_p_bnd[adjID];
                     
//                   gu_nb_vx = dUdxi_p_bnd[adjID][0];//dUdx_p_bnd[adjID];
//                   gu_nb_vy = dUdxi_p_bnd[adjID][1];//dUdy_p_bnd[adjID];
//                   gu_nb_vz = dUdxi_p_bnd[adjID][2];//dUdz_p_bnd[adjID];
                     
                     u_nb = U[adjID];
                 }
                 else
                 {
                     //u_nb     = ghost->getVal(adjID-Nel,0);
                     u_nb     = U[gEl];
                     gu_nb_vx = gu_c_x->getVal(lid,0);
                     gu_nb_vy = gu_c_y->getVal(lid,0);
                     gu_nb_vz = gu_c_z->getVal(lid,0);
                     
                 }
                 
                 Vec3D* nj          = normals[gEl][j];
                 Vec3D* rj          = rvector[gEl][j];
                 
                 double alpha     = DotVec3D(nj,rj);
                 
                 Vec3D* nf_m_arf    = new Vec3D;

                 nf_m_arf->c0=nj->c0-alpha*rj->c0;
                 nf_m_arf->c1=nj->c1-alpha*rj->c1;
                 nf_m_arf->c2=nj->c2-alpha*rj->c2;
                 //std::cout << alpha << std::endl;
                 dphi_dn = alpha * (u_nb - u_c)/dr[gEl][j] +  0.5 * ((gu_nb_vx + gu_c_vx) * nf_m_arf->c0
                                                                  +  (gu_nb_vy + gu_c_vy) * nf_m_arf->c1
                                                                  +  (gu_nb_vz + gu_c_vz) * nf_m_arf->c2);
                 
                 sum_phix = sum_phix+dphi_dn*dxfxc[gEl][j]->c0*dS[gEl][j];
                 sum_phiy = sum_phiy+dphi_dn*dxfxc[gEl][j]->c1*dS[gEl][j];
                 sum_phiz = sum_phiz+dphi_dn*dxfxc[gEl][j]->c2*dS[gEl][j];
                 
                 delete nf_m_arf;
             }
             
             Vol = vol[gEl];
             
             dudx->setVal(i,0,gu_c_x->getVal(i,0));
             dudx->setVal(i,1,gu_c_y->getVal(i,0));
             dudx->setVal(i,2,gu_c_z->getVal(i,0));
             
             gu_c_x->setVal(i,0,1.0/Vol*sum_phix);
             gu_c_y->setVal(i,0,1.0/Vol*sum_phiy);
             gu_c_z->setVal(i,0,1.0/Vol*sum_phiz);
             
             L2normx = L2normx+ sqrt((gu_c_x->getVal(i,0)-dudx->getVal(i,0))*(gu_c_x->getVal(i,0)-dudx->getVal(i,0)));

             L2normy = L2normy+ sqrt((gu_c_y->getVal(i,0)-dudx->getVal(i,1))*(gu_c_y->getVal(i,0)-dudx->getVal(i,1)));

             L2normz = L2normz+ sqrt((gu_c_z->getVal(i,0)-dudx->getVal(i,2))*(gu_c_z->getVal(i,0)-dudx->getVal(i,2)));
         }
        
        dUdx_p_bnd.clear();
        dUdy_p_bnd.clear();
        dUdz_p_bnd.clear();
        
        MPI_Allreduce(&L2normx, &L2normx_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(&L2normy, &L2normy_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        MPI_Allreduce(&L2normz, &L2normz_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        t = clock()-t;
        double tn = ((double)t);
        double tmax = 0.0;
        MPI_Allreduce(&tn, &tmax, 1, MPI_DOUBLE, MPI_MAX, comm);
        tmax=tmax/CLOCKS_PER_SEC;
        if(rank == 0)
        {
            std::cout << it << " " <<  L2normx_max <<","<<L2normy_max<<","<<L2normz_max << " time = " << tmax << std::endl;
        }
        
        if(L2normx_max<1.0e-07 && L2normy_max<1.0e-07 && L2normz_max<1.0e-07)
        {
            break;
        }
        
    }
    
    return dudx;
}

inline Array<double>* Gradients::ComputedUdx_LSQ_US3D_v2(Partition* Pa, std::map<int,double> U,Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   std::vector<Vert> LocalVs             = Pa->getLocalVerts();
   std::map<int,std::vector<int> > gE2lV = Pa->getGlobElem2LocVerts();
   std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
   std::map<int,int> gV2lV               = Pa->getGlobalVert2LocalVert();
   std::map<int,int> gE2lE               = Pa->getGlobalElement2LocalElement();
   std::vector<int> Loc_Elem             = Pa->getLocElem();
   int nLoc_Elem                         = Loc_Elem.size();
   int Nel = Pa->getGlobalPartition()->getNrow();
   Array<int>* ifn = meshTopo->getIFN();
   std::map<int,std::vector<int> > iee_vec = Pa->getIEEpartmap();
   std::vector<std::vector<double> > iee_dist;
   std::vector<double> dist;
   double* Pijk = new double[8*3];
   double* Padj = new double[8*3];

   double d;
   int loc_vid,adjID,elID;
   int cou = 0;
   int offset = Pa->getParallelState()->getOffset(world_rank);
   Vert* Vc = new Vert;
   int lid = 0;
   double u_ijk, u_po;
   dudx = new Array<double>(nLoc_Elem,3);
//
   for(int i=0;i<nLoc_Elem;i++)
   {
       int nadj = 6;
       Array<double>* Vrt_T = new Array<double>(3,nadj);
       Array<double>* Vrt   = new Array<double>(nadj,3);
       Array<double>* b     = new Array<double>(nadj,1);
//
//       std::vector<std::vector<int> > MatVrt;
//
       for(int q=0;q<nadj;q++)
       {
           for(int j=0;j<3;j++)
           {
               Vrt_T->setVal(j,q,0.0);
               Vrt->setVal(q,j,0.0);
           }
       }
//
       int elID = Loc_Elem[i];
       for(int k=0;k<gE2lV[elID].size();k++)
       {
           loc_vid     = gE2lV[elID][k];
           Pijk[k*3+0] = LocalVs[loc_vid].x;
           Pijk[k*3+1] = LocalVs[loc_vid].y;
           Pijk[k*3+2] = LocalVs[loc_vid].z;
       }
       Vert* Vijk = ComputeCenterCoord(Pijk,8);
       u_ijk = U[elID];
       int t = 0;
       for(int j=0;j<6;j++)
       {
             int adjID = iee_vec[elID][j];

           if(adjID<Nel)
           {
               lid = gE2lE[adjID];
               u_po = U[adjID];

               for(int k=0;k<gE2lV[adjID].size();k++)
               {
                   loc_vid     = gE2lV[adjID][k];
                   Padj[k*3+0] = LocalVs[loc_vid].x;
                   Padj[k*3+1] = LocalVs[loc_vid].y;
                   Padj[k*3+2] = LocalVs[loc_vid].z;
               }

               Vert* Vadj = ComputeCenterCoord(Padj,8);

               d = sqrt((Vadj->x-Vijk->x)*(Vadj->x-Vijk->x)+
                        (Vadj->y-Vijk->y)*(Vadj->y-Vijk->y)+
                        (Vadj->z-Vijk->z)*(Vadj->z-Vijk->z));

               Vrt->setVal(t,0,(1.0/d)*(Vadj->x-Vijk->x));
               Vrt->setVal(t,1,(1.0/d)*(Vadj->y-Vijk->y));
               Vrt->setVal(t,2,(1.0/d)*(Vadj->z-Vijk->z));

               b->setVal(t,0,(1.0/d)*(u_po-u_ijk));
               delete Vadj;
               dist.push_back(d);
               t++;
           }
           else
           {
               int fid = gE2gF[elID][j];

               Vc->x = 0.0;Vc->y = 0.0;Vc->z = 0.0;

               for(int s=0;s<4;s++)
               {
                   int gvid = ifn->getVal(fid,s);
                   int lvid = gV2lV[gvid];

                   Vc->x = Vc->x+LocalVs[lvid].x;
                   Vc->y = Vc->y+LocalVs[lvid].y;
                   Vc->z = Vc->z+LocalVs[lvid].z;
               }
               //std::cout << std::endl;

               Vc->x = Vc->x/4.0;
               Vc->y = Vc->y/4.0;
               Vc->z = Vc->z/4.0;

               d = sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                            (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                            (Vc->z-Vijk->z)*(Vc->z-Vijk->z));

               u_po = ghost->getVal(adjID-Nel,0);

               //double u_fpo = bound->getVal(adjID-Nel,0);

               Vrt->setVal(t,0,(1.0/d)*(Vc->x-Vijk->x));
               Vrt->setVal(t,1,(1.0/d)*(Vc->y-Vijk->y));
               Vrt->setVal(t,2,(1.0/d)*(Vc->z-Vijk->z));
               b->setVal(t,0,(1.0/d)*(u_po-u_ijk));
               t++;
               dist.push_back(d);

           }
//
           double* A_cm = new double[nadj*3];
           double sum =0.0;
           for(int s=0;s<nadj;s++)
           {
               for(int j=0;j<3;j++)
               {
                   A_cm[j*nadj+s] = Vrt->getVal(s,j);
                   //std::cout << Vrt->getVal(s,j) << " ";
               }
               //std::cout << std::endl;
           }

           Array<double>* x = SolveQR(A_cm,nadj,3,b);
//
           dudx->setVal(i,0,x->getVal(0,0));
           dudx->setVal(i,1,x->getVal(1,0));
           dudx->setVal(i,2,x->getVal(2,0));
//
           delete[] A_cm;
           delete x;
       }

       delete Vrt_T;
       delete Vrt;
       delete b;

       iee_dist.push_back(dist);
       dist.clear();
   }

   delete Vc;
   delete[] Pijk;
   delete[] Padj;
//
   return dudx;
}


inline Array<double>* Gradients::getdUdXi()
{
    return dudx;
}

#endif
