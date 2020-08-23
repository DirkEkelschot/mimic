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
    if (strcmp(recon_type, "MMG") == 0)
    {
        dudx = ComputedUdx_MGG(Pa,U,meshTopo,ghost,comm);
    }
    if (strcmp(recon_type, "LSQ") == 0)
    {
        std::cout << "We still need to copy the correct routine for the LSQ grad recon." << std::endl;
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
        Array<double>* dudx    = new Array<double>(nLoc_Elem,3);
        
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
        
        for(int it=0;it<100;it++)
        {
            for(int i=0;i<nLoc_Elem;i++)
            {
                dudx->setVal(i,0,gu_c_x->getVal(i,0));
                dudx->setVal(i,1,gu_c_y->getVal(i,0));
                dudx->setVal(i,2,gu_c_z->getVal(i,0));
            }
            
             //communicate grad phi!!!
            
            std::map<int,double> dUdx_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_x, comm);
            std::map<int,double> dUdy_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_y, comm);
            std::map<int,double> dUdz_p_bnd = Pa->CommunicateAdjacentDataUS3D(gu_c_z, comm);
            
            //std::map<int,std::vector<double> > dUdxi_p_bnd = Pa->CommunicateAdjacentDataUS3DNew(gu_c_old, comm);

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
                     //int fid = ief_vec[gEl][j];
                     
                     if(adjID<Nel)
                     {
                         //l_adjid = gE2lE[adjID];
                         gu_nb_vx = dUdx_p_bnd[adjID];
                         gu_nb_vy = dUdy_p_bnd[adjID];
                         gu_nb_vz = dUdz_p_bnd[adjID];
                         
    //                     gu_nb_vx = dUdxi_p_bnd[adjID][0];//dUdx_p_bnd[adjID];
    //                     gu_nb_vy = dUdxi_p_bnd[adjID][1];//dUdy_p_bnd[adjID];
    //                     gu_nb_vz = dUdxi_p_bnd[adjID][2];//dUdz_p_bnd[adjID];
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
                 gu_c_x->setVal(i,0,1.0/Vol*sum_phix);
                 gu_c_y->setVal(i,0,1.0/Vol*sum_phiy);
                 gu_c_z->setVal(i,0,1.0/Vol*sum_phiz);
             }
            
            dUdx_p_bnd.clear();
            dUdy_p_bnd.clear();
            dUdz_p_bnd.clear();
                    
            L2normx = 0.0;
            L2normy = 0.0;
            L2normz = 0.0;
            
            for(int n=0;n<nLoc_Elem;n++)
            {
                L2normx = L2normx+ sqrt((gu_c_x->getVal(n,0)-dudx->getVal(n,0))*(gu_c_x->getVal(n,0)-dudx->getVal(n,0)));

                L2normy = L2normy+ sqrt((gu_c_y->getVal(n,0)-dudx->getVal(n,1))*(gu_c_y->getVal(n,0)-dudx->getVal(n,1)));

                L2normz = L2normz+ sqrt((gu_c_z->getVal(n,0)-dudx->getVal(n,2))*(gu_c_z->getVal(n,0)-dudx->getVal(n,2)));
            }
            
            MPI_Allreduce(&L2normx, &L2normx_max, 1, MPI_DOUBLE, MPI_MAX, comm);
            MPI_Allreduce(&L2normy, &L2normy_max, 1, MPI_DOUBLE, MPI_MAX, comm);
            MPI_Allreduce(&L2normz, &L2normz_max, 1, MPI_DOUBLE, MPI_MAX, comm);
            
            if(rank == 0)
            {
                std::cout << it << " " <<  L2normx_max <<","<<L2normy_max<<","<<L2normz_max << std::endl;
            }
            if(L2normx_max<1.0e-07 && L2normy_max<1.0e-07 && L2normz_max<1.0e-07)
            {
                
                break;
            }
            
        }
        
        return dudx;
}

inline Array<double>* Gradients::getdUdXi()
{
    return dudx;
}

#endif
