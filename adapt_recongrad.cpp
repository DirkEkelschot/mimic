#include "adapt_recongrad.h"

Array<double>* ComputedUdx_LSQ(Partition* P, std::vector<double> U, int Nel, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::vector<Vert> locVerts            = P->getLocalVerts();
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    int loc_vid = 0;
    int np = 8;
    std::map<int, int>::iterator itmap;
    int e = 0;
    int offset  = P->getLocalPartition()->getOffset(rank);
    int* xadj   = P->getXadj();
    int* adjcny = P->getAdjcny();
    int nloc    = P->getLocalPartition()->getNrow();
    int lid;
    Array<double>* grad = new Array<double>(nloc,3);
    std::vector<std::vector<double> > store_coords;
    std::vector<double> upo_stored;
    std::vector<double> upijk_stored;
    std::vector<double> fac_stored;
    for(int i = 0;i<nloc;i++)
    {
        int start = xadj[i];
        int end   = xadj[i+1];
        int nadj  = xadj[i+1]-xadj[i];
        
        Array<double>* Vrt_T = new Array<double>(3,nadj);
        Array<double>* Vrt   = new Array<double>(nadj,3);
        Array<double>* b     = new Array<double>(nadj,1);
        
        for(int q=0;q<nadj;q++)
        {
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,q,0.0);
                Vrt->setVal(q,j,0.0);
            }
        }
        lid = gE2lE[i+offset];
        double u_ijk = U[lid];
        std::vector<int> vijkIDs = gE2lV[i+offset];
        double* Pijk = new double[np*3];
        for(int k=0;k<vijkIDs.size();k++)
        {
            loc_vid     = vijkIDs[k];
            Pijk[k*3+0] = locVerts[loc_vid].x;
            Pijk[k*3+1] = locVerts[loc_vid].y;
            Pijk[k*3+2] = locVerts[loc_vid].z;
        }
        
        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        int tel   = 0;
        std::vector<double> tmp;
       
        for(int j=start;j<end;j++)
        {
            int adjEl_id = adjcny[j];
            
            lid = gE2lE[adjEl_id];
            double u_po = U[lid];
            double* Po = new double[np*3];

            for(int k=0;k<gE2lV[adjEl_id].size();k++)
            {
                loc_vid   = gE2lV[adjEl_id][k];
                Po[k*3+0] = locVerts[loc_vid].x;
                Po[k*3+1] = locVerts[loc_vid].y;
                Po[k*3+2] = locVerts[loc_vid].z;
            }


            Vert* Vpo = ComputeCenterCoord(Po,8);
            //double Volpo = ComputeVolumeHexCell(Po);

            double wi = sqrt((Vpo->x-Vijk->x)*(Vpo->x-Vijk->x)+
                             (Vpo->y-Vijk->y)*(Vpo->y-Vijk->y)+
                             (Vpo->z-Vijk->z)*(Vpo->z-Vijk->z));

            double fac = (wi);

            Vrt->setVal(tel,0,(Vpo->x-Vijk->x));
            Vrt->setVal(tel,1,(Vpo->y-Vijk->y));
            Vrt->setVal(tel,2,(Vpo->z-Vijk->z));
//
            upo_stored.push_back(u_po);
            upijk_stored.push_back(u_ijk);
            fac_stored.push_back(fac);
            delete[] Po;
            delete Vpo;
            tel++;
        }
       
//        double min = *min_element(fac_stored.begin(), fac_stored.end());
//        double max = *max_element(fac_stored.begin(), fac_stored.end());

        std::vector<double> fac_stored_update(nadj);

        for(int s=0;s<nadj;s++)
        {
            fac_stored_update[s] = (1.0/fac_stored[s]);
        }

        double* A_cm = new double[nadj*3];
        double sum =0.0;
        for(int s=0;s<nadj;s++)
        {
            for(int j=0;j<3;j++)
            {
                A_cm[j*nadj+s] = fac_stored_update[s]*Vrt->getVal(s,j);
            }
        }


        for(int s=0;s<nadj;s++)
        {
            b->setVal(s,0,fac_stored_update[s]*(upo_stored[s]-upijk_stored[s]));
        }
        Array<double>* x = SolveQR(A_cm,nadj,3,b);


        grad->setVal(i,0,x->getVal(0,0));
        grad->setVal(i,1,x->getVal(1,0));
        grad->setVal(i,2,x->getVal(2,0));


        upo_stored.clear();
        upijk_stored.clear();
        fac_stored.clear();
        fac_stored_update.clear();
        store_coords.clear();
        delete[] Pijk;
        vijkIDs.clear();
        delete Vrt_T;
        delete Vrt;
        delete b;
        delete[] A_cm;
        delete x;
        delete Vijk;
        
    }
    
    return grad;
}




Array<double>* ComputedUdx_LSQ_US3D_v1(Partition* P, ParallelState* pstate, ParArray<int>* iee, i_part_map* iee_vec,std::map<int,std::vector<int> > ief_vec, Array<int>* ifn, Array<int>* ief, int Nel, std::vector<double> U, Array<double>* ghost, Array<double>* bound, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<Vert> LocalVs = P->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gF = P->getglobElem2globFaces();
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    std::vector<int> ElemPart             = P->getLocElem();
    std::vector<int> Loc_Elem          = P->getLocElem();
//    int* xadj                             = P->getXadj();
//    int* adjcny                           = P->getAdjcny();
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;
    double* Pijk = new double[8*3];
    double* Padj = new double[8*3];
    
    double d;
    int loc_vid,adjID,elID;
    int cou = 0;
    Vert* Vc = new Vert;
    int lid = 0;
    double u_ijk, u_po;

    int nLocElem = Loc_Elem.size();
    Array<double>* grad = new Array<double>(nLocElem,3);

    for(int i=0;i<nLocElem;i++)
    {
        int nadj = 6;
    
       // std::vector<int> adj_el_real;
//        for(int q=0;q<6;q++)
//        {
//            if(iee->getVal(i,q)<Nel)
//            {
//                adj_el_real.push_back(iee->getVal(i,q));
//                nadj++;
//            }
//        }
                
        Array<double>* Vrt_T = new Array<double>(3,nadj);
        Array<double>* Vrt   = new Array<double>(nadj,3);
        Array<double>* b     = new Array<double>(nadj,1);
        
        std::vector<std::vector<int> > MatVrt;
        
        for(int q=0;q<nadj;q++)
        {
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,q,0.0);
                Vrt->setVal(q,j,0.0);
            }
        }
        
        int elID = Loc_Elem[i];
        for(int k=0;k<gE2lV[elID].size();k++)
        {
            loc_vid     = gE2lV[elID][k];
            Pijk[k*3+0] = LocalVs[loc_vid].x;
            Pijk[k*3+1] = LocalVs[loc_vid].y;
            Pijk[k*3+2] = LocalVs[loc_vid].z;
        }
        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        
        lid = gE2lE[elID];
        u_ijk = U[lid];
        int t = 0;
        for(int j=0;j<6;j++)
        {
            
            adjID   = iee_vec->i_map[elID][j];
            //int fid2 = ief_vec[elID][j];
//
            if(adjID<Nel)
            {
                lid = gE2lE[adjID];
                u_po = U[lid];

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
                
                Vrt->setVal(t,0,(Vadj->x-Vijk->x));
                Vrt->setVal(t,1,(Vadj->y-Vijk->y));
                Vrt->setVal(t,2,(Vadj->z-Vijk->z));
                                
                b->setVal(t,0,u_po-u_ijk);
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
                
                Vc->x = Vc->x/4.0;
                Vc->y = Vc->y/4.0;
                Vc->z = Vc->z/4.0;

                d = sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                             (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                             (Vc->z-Vijk->z)*(Vc->z-Vijk->z));

                u_po = ghost->getVal(adjID-Nel,0);
                
                //double u_fpo = bound->getVal(adjID-Nel,0);
 
                Vrt->setVal(t,0,(Vc->x-Vijk->x));
                Vrt->setVal(t,1,(Vc->y-Vijk->y));
                Vrt->setVal(t,2,(Vc->z-Vijk->z));
//
                b->setVal(t,0,u_po-u_ijk);
                
                t++;
                dist.push_back(d);
            }
            
            double* A_cm = new double[nadj*3];
            double sum =0.0;
            for(int s=0;s<nadj;s++)
            {
                for(int j=0;j<3;j++)
                {
                    A_cm[j*nadj+s] = Vrt->getVal(s,j);
                }
            }
            
            Array<double>* x = SolveQR(A_cm,nadj,3,b);

            grad->setVal(i,0,x->getVal(0,0));
            grad->setVal(i,1,x->getVal(1,0));
            grad->setVal(i,2,x->getVal(2,0));
            
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
    delete[] Padj;
    delete[] Pijk;
    return grad;
}

//(Partition* Pa, std::map<int,double> U,
//Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm)

Array<double>* ComputedUdx_LSQ_US3D_v2(Partition* P, ParallelState* pstate, ParArray<int>* iee, Array<int>* ifn, Array<int>* ief, int Nel, std::vector<double> U, Array<double>* ghost, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<Vert> LocalVs = P->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV = P->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gF = P->getglobElem2globFaces();
    std::map<int,int> gV2lV               = P->getGlobalVert2LocalVert();
    std::map<int,int> gE2lE               = P->getGlobalElement2LocalElement();
    std::vector<int> ElemPart             = P->getLocElem();
    
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;
    double* Pijk = new double[8*3];
    double* Padj = new double[8*3];

    double d;
    int loc_vid,adjID,elID;
    int cou = 0;
    int offset = pstate->getOffset(world_rank);
    Vert* Vc = new Vert;
    int lid = 0;
    double u_ijk, u_po;

    
    Array<double>* grad = new Array<double>(iee->getNrow(),3);
    
    for(int i=0;i<iee->getNrow();i++)
    {
        int nadj = 6;
                
        Array<double>* Vrt_T = new Array<double>(3,nadj);
        Array<double>* Vrt   = new Array<double>(nadj,3);
        Array<double>* b     = new Array<double>(nadj,1);
        
        std::vector<std::vector<int> > MatVrt;
        
        for(int q=0;q<nadj;q++)
        {
            for(int j=0;j<3;j++)
            {
                Vrt_T->setVal(j,q,0.0);
                Vrt->setVal(q,j,0.0);
            }
        }
        
        int elID = i+offset;//ElemPart[i];
        for(int k=0;k<gE2lV[elID].size();k++)
        {
            loc_vid     = gE2lV[elID][k];
            Pijk[k*3+0] = LocalVs[loc_vid].x;
            Pijk[k*3+1] = LocalVs[loc_vid].y;
            Pijk[k*3+2] = LocalVs[loc_vid].z;
        }
        Vert* Vijk = ComputeCenterCoord(Pijk,8);
        lid = gE2lE[i+offset];
        u_ijk = U[lid];
        int t = 0;
        for(int j=0;j<6;j++)
        {
            adjID = iee->getVal(i,j);

            if(adjID<Nel)
            {
                lid = gE2lE[adjID];
                u_po = U[lid];

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
                
                Vrt->setVal(t,0,(Vadj->x-Vijk->x));
                Vrt->setVal(t,1,(Vadj->y-Vijk->y));
                Vrt->setVal(t,2,(Vadj->z-Vijk->z));
                
                b->setVal(t,0,u_po-u_ijk);
                delete Vadj;
                dist.push_back(d);
                t++;
            }
            else
            {
                //int fid = gE2gF[elID][j];
                int fid = ief->getVal(elID,j);
                //double* face_adj = new double[4*3];
                Vc->x = 0.0;Vc->y = 0.0;Vc->z = 0.0;
                
                for(int s=0;s<4;s++)
                {
                    int gvid = ifn->getVal(fid,s);
                    int lvid = gV2lV[gvid];
                    
                    Vc->x = Vc->x+LocalVs[lvid].x;
                    Vc->y = Vc->y+LocalVs[lvid].y;
                    Vc->z = Vc->z+LocalVs[lvid].z;
                }
                
                Vc->x = Vc->x/4.0;
                Vc->y = Vc->y/4.0;
                Vc->z = Vc->z/4.0;
                
                d = sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                             (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                             (Vc->z-Vijk->z)*(Vc->z-Vijk->z));

                u_po = ghost->getVal(adjID-Nel,0);
                
 
                Vrt->setVal(t,0,(Vc->x-Vijk->x));
                Vrt->setVal(t,1,(Vc->y-Vijk->y));
                Vrt->setVal(t,2,(Vc->z-Vijk->z));
                b->setVal(t,0,u_po-u_ijk);
                t++;
                dist.push_back(d);
            }
        }
        
        double* A_cm = new double[nadj*3];
        double sum =0.0;
        for(int s=0;s<nadj;s++)
        {
            for(int j=0;j<3;j++)
            {
                    A_cm[j*nadj+s] = Vrt->getVal(s,j);
            }
        }
            
        Array<double>* x = SolveQR(A_cm,nadj,3,b);
        
        grad->setVal(i,0,x->getVal(0,0));
        grad->setVal(i,1,x->getVal(1,0));
        grad->setVal(i,2,x->getVal(2,0));
        
        delete Vrt;
        delete Vrt_T;
        delete[] A_cm;
        delete x;
        
                
        delete Vrt_T;
        delete Vrt;
        delete b;
        
        iee_dist.push_back(dist);
        dist.clear();
    }
    
    delete Vc;
    delete[] Pijk;
    delete[] Padj;
    
    return grad;
}



Array<double>* ComputedUdx_LSQ_US3D_v3(Partition* Pa, std::map<int,double> U,Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm)
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
   i_part_map* iee_vec = Pa->getIEEpartmap();
   std::vector<std::vector<double> > iee_dist;
   std::vector<double> dist;
   double* Pijk = new double[8*3];
   double* Padj = new double[8*3];

   double d;
   int loc_vid,adjID,elID;
   int cou = 0;
   Vert* Vc = new Vert;
   int lid = 0;
   double u_ijk, u_po;
   Array<double>* dudx = new Array<double>(nLoc_Elem,3);
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
           int adjID = iee_vec->i_map[elID][j];

           if(adjID<Nel)
           {
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
      }
       
      double* A_cm = new double[nadj*3];
      double sum =0.0;
      for(int s=0;s<nadj;s++)
      {
          for(int j=0;j<3;j++)
          {
              A_cm[j*nadj+s] = Vrt->getVal(s,j);
          }
      }

       Array<double>* x = SolveQR(A_cm,nadj,3,b);
//
       dudx->setVal(i,0,x->getVal(0,0));
       dudx->setVal(i,1,x->getVal(1,0));
       dudx->setVal(i,2,x->getVal(2,0));
    
       delete[] A_cm;
       delete x;

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



Array<double>* ComputedUdx_MGG(Partition* Pa, std::map<int,double> U, Mesh_Topology* meshTopo, Array<double>* ghost, MPI_Comm comm)
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
    i_part_map* iee_vec = Pa->getIEEpartmap();
    Array<double>* gu_c_old    = new Array<double>(nLoc_Elem,3);
    
    std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
    
    for(int i=0;i<nLoc_Elem;i++)
    {
        gu_c_x->setVal(i,0,0.0);
        gu_c_y->setVal(i,0,0.0);
        gu_c_z->setVal(i,0,0.0);
        
        gu_c_old->setVal(i,0,0.0);
        gu_c_old->setVal(i,1,0.0);
        gu_c_old->setVal(i,2,0.0);
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
    Vec3D* nj;
    Vec3D* rj;
    for(int it=0;it<1000;it++)
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
                 adjID   = iee_vec->i_map[gEl][j];
                 
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
                     
                     u_nb     = 0.0;
                     gu_nb_vx = 0.0;
                     gu_nb_vy = 0.0;
                     gu_nb_vz = 0.0;
                     
                 }
                 
                 nj          = normals[gEl][j];
                 rj          = rvector[gEl][j];
                 
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
             
             gu_c_old->setVal(i,0,gu_c_x->getVal(i,0));
             gu_c_old->setVal(i,1,gu_c_y->getVal(i,0));
             gu_c_old->setVal(i,2,gu_c_z->getVal(i,0));
             
             gu_c_x->setVal(i,0,1.0/Vol*sum_phix);
             gu_c_y->setVal(i,0,1.0/Vol*sum_phiy);
             gu_c_z->setVal(i,0,1.0/Vol*sum_phiz);
             
             L2normx = L2normx+ sqrt((gu_c_x->getVal(i,0)-gu_c_old->getVal(i,0))*(gu_c_x->getVal(i,0)-gu_c_old->getVal(i,0)));

             L2normy = L2normy+ sqrt((gu_c_y->getVal(i,0)-gu_c_old->getVal(i,1))*(gu_c_y->getVal(i,0)-gu_c_old->getVal(i,1)));

             L2normz = L2normz+ sqrt((gu_c_z->getVal(i,0)-gu_c_old->getVal(i,2))*(gu_c_z->getVal(i,0)-gu_c_old->getVal(i,2)));
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
        
        if(L2normx_max<1.0e-06 && L2normy_max<1.0e-06 && L2normz_max<1.0e-06)
        {
            break;
        }
        
    }
    
    return gu_c_old;
}

