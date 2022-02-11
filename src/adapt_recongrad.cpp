#include "adapt_recongrad.h"



std::map<int,Array<double>* > ComputedUdx_LSQ_Vrt_US3D(Partition* Pa, std::map<int,Array<double>* > Ue, std::map<int,Array<double>* > Uv, Mesh_Topology* meshTopo, std::map<int,double> gbMap, MPI_Comm comm)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   std::vector<Vert*> LocalVs                 = Pa->getLocalVerts();
   std::map<int,std::vector<int> > gE2lV      = Pa->getGlobElem2LocVerts();
   std::map<int,std::vector<int> > gE2gF      = Pa->getglobElem2globFaces();
   std::map<int,int> gV2lV                    = Pa->getGlobalVert2LocalVert();
   std::map<int,int> gE2lE                    = Pa->getGlobalElement2LocalElement();
   std::vector<int> Loc_Elem                  = Pa->getLocElem();
   std::map<int,std::vector<int> > scheme_E2V = meshTopo->getScheme_E2V();
   int nLoc_Elem                              = Loc_Elem.size();
   
   
    
   int Nel                      = Pa->getGlobalPartition()->getNrow();
   i_part_map*  ifn_vec         = Pa->getIFNpartmap();
   i_part_map* ief_part_map     = Pa->getIEFpartmap();
   i_part_map*  iee_vec         = Pa->getIEEpartmap();
   i_part_map* if_Nv_part_map   = Pa->getIF_Nvpartmap();

//   std::vector<std::vector<double> > iee_dist;
//   std::vector<double> dist;

   std::map<int,Array<double>* > dudx_map;
   double d;
   int loc_vid,adjID,elID;
   int cou = 0;
   Vert* Vc = new Vert;
   Vert* Vadj = new Vert;
   int lid = 0;
   double u_ijk, u_po;
    
   int el_contr = 1;
   int nadj_el  = 0;

   std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
   std::map<int,int> LocElem2Nv = Pa->getLocElem2Nv();

   for(int i=0;i<nLoc_Elem;i++)
   {
       int bflip    = 0;
       int elID     = Loc_Elem[i];
       int NvPEl    = LocElem2Nv[elID];
       
       if(el_contr == 1)
       {
           nadj_el  = LocElem2Nf[elID];
       }
       
       std::vector<int> vrts = scheme_E2V[elID];
       
       int nadj_vrts  = vrts.size();
       int nadj_tot   = nadj_vrts+nadj_el;

       Array<double>* Vrt_T = new Array<double>(3,nadj_tot);
       Array<double>* Vrt   = new Array<double>(nadj_tot,3);
       Array<double>* b     = new Array<double>(nadj_tot,1);

       for(int q=0;q<nadj_tot;q++)
       {
           for(int j=0;j<3;j++)
           {
               Vrt_T->setVal(j,q,0.0);
               Vrt->setVal(q,j,0.0);
           }
       }
       
       double* Pijk = new double[NvPEl*3];
       
       for(int k=0;k<gE2lV[elID].size();k++)
       {
           loc_vid     = gE2lV[elID][k];
           Pijk[k*3+0] = LocalVs[loc_vid]->x;
           Pijk[k*3+1] = LocalVs[loc_vid]->y;
           Pijk[k*3+2] = LocalVs[loc_vid]->z;
       }
       
       Vert* Vijk   = ComputeCentroidCoord(Pijk,NvPEl);
       
       u_ijk        = Ue[elID]->getVal(0,0);
       int t        = 0;

       if(el_contr == 1)
       {
           for(int j=0;j<nadj_el;j++)
           {
               int adjID = iee_vec->i_map[elID][j];
               double* Padj = new double[gE2lV[adjID].size()*3];

               if(adjID<Nel)
               {
                   u_po = Ue[adjID]->getVal(0,0);
        
                   for(int k=0;k<gE2lV[adjID].size();k++)
                   {
                       loc_vid     = gE2lV[adjID][k];
                       Padj[k*3+0] = LocalVs[loc_vid]->x;
                       Padj[k*3+1] = LocalVs[loc_vid]->y;
                       Padj[k*3+2] = LocalVs[loc_vid]->z;
                   }
                   
                   Vert* Vadj_el = ComputeCentroidCoord(Padj,gE2lV[adjID].size());
                   
                   d = sqrt((Vadj_el->x-Vijk->x)*(Vadj_el->x-Vijk->x)+
                            (Vadj_el->y-Vijk->y)*(Vadj_el->y-Vijk->y)+
                            (Vadj_el->z-Vijk->z)*(Vadj_el->z-Vijk->z));

                   Vrt->setVal(t,0,(1.0/d)*(Vadj_el->x-Vijk->x));
                   Vrt->setVal(t,1,(1.0/d)*(Vadj_el->y-Vijk->y));
                   Vrt->setVal(t,2,(1.0/d)*(Vadj_el->z-Vijk->z));
                   
                   b->setVal(t,0,(1.0/d)*(u_po-u_ijk));
                   delete Vadj_el;
                   //dist.push_back(d);
                   t++;
               }
               else
               {
                   bflip = 1;

                   int fid    = ief_part_map->i_map[elID][j];
                   int NvPerF = if_Nv_part_map->i_map[fid][0];
                   
                   Vc->x = 0.0;
                   Vc->y = 0.0;
                   Vc->z = 0.0;
                   
                   for(int s=0;s<NvPerF;s++)
                   {
                       int gvid = ifn_vec->i_map[fid][s];
                       int lvid = gV2lV[gvid];

                       Vc->x = Vc->x+LocalVs[lvid]->x;
                       Vc->y = Vc->y+LocalVs[lvid]->y;
                       Vc->z = Vc->z+LocalVs[lvid]->z;
                   }

                   Vc->x = Vc->x/NvPerF;
                   Vc->y = Vc->y/NvPerF;
                   Vc->z = Vc->z/NvPerF;

                   d = sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                            (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                            (Vc->z-Vijk->z)*(Vc->z-Vijk->z));

//                   double rtje = sqrt(Vc->x*Vc->x+(Vc->y-0.5)*(Vc->y-0.5)+(Vc->z-0.5)*(Vc->z-0.5));
//                   double Utje = 0.1*tanh(50*(rtje-0.5))+1.0;
                   
                   double Utje = gbMap[elID];
                   //double Utje    = 0.1*sin(50*Vc->x*Vc->z)+atan(0.1/((sin(5.0*Vc->y)-2.0*Vc->x*Vc->z)));
                   
                   u_po = Utje;
                   //u_po = u_ijk;
                   //u_po = U[elID];
                   
                   //double u_fpo = ghost->getVal(adjID-Nel,0);

                   Vrt->setVal(t,0,(1.0/d)*(Vc->x-Vijk->x));
                   Vrt->setVal(t,1,(1.0/d)*(Vc->y-Vijk->y));
                   Vrt->setVal(t,2,(1.0/d)*(Vc->z-Vijk->z));

                   b->setVal(t,0,(1.0/d)*(u_po-u_ijk));
                   t++;
                   
                   //dist.push_back(d);
               }
               
               delete[] Padj;
          }
       }
       
       
       if(bflip==1)
       {
           Array<double>* b_el     = new Array<double>(nadj_el,1);
           double* A_cm_el         = new double[nadj_tot*3];

           for(int q=0;q<nadj_el;q++)
           {
               for(int j=0;j<3;j++)
               {
                   A_cm_el[j*nadj_el+q] = Vrt->getVal(q,j);
               }
               
               b_el->setVal(q,0,b->getVal(q,0));
           }
           
            Array<double>* x = SolveQR(A_cm_el,nadj_el,3,b_el);
            dudx_map[elID] = x;
            delete[] A_cm_el;
            delete b_el;
            //delete x;
       }
       else
       {
           for(int j=0;j<nadj_vrts;j++)
           {
               int gvid    = vrts[j];
               double Uvrt = Uv[gvid]->getVal(0,0);

               int lvid = gV2lV[gvid];

               Vadj->x = LocalVs[lvid]->x;
               Vadj->y = LocalVs[lvid]->y;
               Vadj->z = LocalVs[lvid]->z;
    //
               d = sqrt((Vadj->x-Vijk->x)*(Vadj->x-Vijk->x)+
                        (Vadj->y-Vijk->y)*(Vadj->y-Vijk->y)+
                        (Vadj->z-Vijk->z)*(Vadj->z-Vijk->z));

               Vrt->setVal(t,0,(1.0/d)*(Vadj->x-Vijk->x));
               Vrt->setVal(t,1,(1.0/d)*(Vadj->y-Vijk->y));
               Vrt->setVal(t,2,(1.0/d)*(Vadj->z-Vijk->z));
               b->setVal(t,0,(1.0/d)*(Uvrt-u_ijk));
               //dist.push_back(d);
               t++;
          }
           
          double* A_cm = new double[nadj_tot*3];
          for(int s=0;s<nadj_tot;s++)
          {
              for(int j=0;j<3;j++)
              {
                  A_cm[j*nadj_tot+s] = Vrt->getVal(s,j);
              }
          }

           Array<double>* x = SolveQR(A_cm,nadj_tot,3,b);
           dudx_map[elID] = x;
           delete[] A_cm;

       }
       
       //delete x;
       delete[] Pijk;
       delete Vrt_T;
       delete Vrt;
       delete b;
       delete Vijk;
       //iee_dist.push_back(dist);
       //dist.clear();
   }
    
   delete Vadj;
   delete Vc;
   LocalVs.clear();
   gE2lV.clear();
   gE2lE.clear();
   gV2lV.clear();
   gE2gF.clear();
   scheme_E2V.clear();
   Loc_Elem.clear();
    
   return dudx_map;
}







/*
std::map<int,Array<double>* > ComputedUdx_LSQ_HO_US3D(Partition* Pa, std::map<int,Array<double>* > Ue, Mesh_Topology* meshTopo, std::map<int,double> gbMap, MPI_Comm comm)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   std::vector<Vert*> LocalVs                 = Pa->getLocalVerts();
   std::map<int,std::vector<int> > gE2lV      = Pa->getGlobElem2LocVerts();
   std::map<int,std::vector<int> > gE2gF      = Pa->getglobElem2globFaces();
   std::map<int,int> gV2lV                    = Pa->getGlobalVert2LocalVert();
   std::map<int,int> gE2lE                    = Pa->getGlobalElement2LocalElement();
   std::vector<int> Loc_Elem                  = Pa->getLocElem();
   //std::map<int,std::vector<int> > scheme_E2V = meshTopo->getScheme_E2V();
   int nLoc_Elem                              = Loc_Elem.size();
    Array<int>* pg = Pa->getGlobalPartition();
    
   int Nel                      = Pa->getGlobalPartition()->getNrow();
   i_part_map*  ifn_vec         = Pa->getIFNpartmap();
    
   i_part_map* ief_part_map     = Pa->getIEFpartmap();
//   i_part_map* ief_adj_part_map = Pa->getIEFADJpartmap();
//   i_part_map* ief_adj2_part_map = Pa->getIEFADJ2partmap();
    
   i_part_map*  iee_vec         = Pa->getIEEpartmap();
    i_part_map*  ien_vec         = Pa->getIENpartmap();

    i_part_map*  ien_part_map         = Pa->getIENpartmap();

//   i_part_map*  iee_adj_vec     = Pa->getIEEADJpartmap();
//   i_part_map*  iee_adj2_vec     = Pa->getIEEADJ2partmap();
    
   i_part_map* if_Nv_part_map   = Pa->getIF_Nvpartmap();
    std::map<int,Vert*> ghostvrts = meshTopo->getGhostVerts();
    std::map<int,std::vector<Vert*> > vfvec = meshTopo->getVfacevector();

//   std::vector<std::vector<double> > iee_dist;
//   std::vector<double> dist;

    
   std::map<int,Array<double>* > dudx_map;
   double d;
   int loc_vid,adjID,elID;
   int cou = 0;
   Vert* Vc = new Vert;
   Vert* Vadj = new Vert;
   int lid = 0;
   double u_ijk, u_po;
    
   int el_contr = 1;
   int nadj_el  = 0;
   int cntbnd = 0;
   int cntbnd2 = 0;
   int cntbnd3 = 0;
   std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
   std::set<int> add2set_lay0;
   std::set<int> add2set_layt;
   std::map<int,Vert*> vrt_collect;
   std::map<int,double> sol_collect;
   std::set<int> add2set_lay1;
    int thismany  = 0;
    int thismany2 = 0;
    int thismany3 = 0;
    Vec3D* v0 = new Vec3D;
    Vec3D* v1 = new Vec3D;
    std::vector<Vert*> face;
    double rdotn;
    double orient0;
    Vec3D* n0 = new Vec3D;
   for(int i=0;i<nLoc_Elem;i++)
    {
        int bflip    = 0;
        int elID     = Loc_Elem[i];
        int NvPEl    = ien_part_map->i_map[elID].size();
        
        double* Pijk = new double[NvPEl*3];
        
        for(int k=0;k<NvPEl;k++)
        {
            int glob_vid     = ien_part_map->i_map[elID][k];
            loc_vid          = gV2lV[glob_vid];
            Pijk[k*3+0]      = LocalVs[loc_vid]->x;
            Pijk[k*3+1]      = LocalVs[loc_vid]->y;
            Pijk[k*3+2]      = LocalVs[loc_vid]->z;
        }
        
        Vert* Vijk      = ComputeCentroidCoord(Pijk,NvPEl);
        u_ijk           = Ue[elID]->getVal(0,0);

        int nadj_tot    = ief_part_map->i_map[elID].size();
       
        for(int j=0;j<nadj_tot;j++)
        {
            int adjid   = iee_vec->i_map[elID][j];
               
            if(vrt_collect.find(adjid)==vrt_collect.end() && adjid<Nel)
            {
                int NvPEladj    = ien_part_map->i_map[adjid].size();

                double* Padj = new double[NvPEladj*3];

                for(int k=0;k<NvPEladj;k++)
                {
                    int glob_vid     = ien_part_map->i_map[adjid][k];
                    loc_vid          = gV2lV[glob_vid];
                    Padj[k*3+0] = LocalVs[loc_vid]->x;
                    Padj[k*3+1] = LocalVs[loc_vid]->y;
                    Padj[k*3+2] = LocalVs[loc_vid]->z;
                }

                Vert* Vadj = ComputeCentroidCoord(Padj,NvPEladj);
                
                delete[] Padj;
                
                vrt_collect[adjid]  = Vadj;
                sol_collect[adjid]  = Ue[adjid]->getVal(0,0);
            
                int nadjadj = iee_vec->i_map[adjID].size();
                
                for(int k=0;k<nadjadj;k++)
                {
                    int adjadjID   = iee_vec->i_map[adjID][k];
                    
                    if(vrt_collect.find(adjadjID)==vrt_collect.end() && adjadjID<Nel)
                    {
                        int NvPEladjadj    = ien_part_map->i_map[adjadjID].size();

                        double* Padjadj = new double[NvPEladjadj*3];

                        for(int l=0;l<NvPEladjadj;l++)
                        {
                            int glob_vid     = ien_part_map->i_map[adjadjID][l];
                            loc_vid          = gV2lV[glob_vid];
                            Padjadj[k*3+0] = LocalVs[loc_vid]->x;
                            Padjadj[k*3+1] = LocalVs[loc_vid]->y;
                            Padjadj[k*3+2] = LocalVs[loc_vid]->z;
                        }

                        Vert* Vadjadj = ComputeCentroidCoord(Padj,NvPEladjadj);
                        
                        delete[] Padjadj;
                        
                        vrt_collect[adjadjID]  = Vadjadj;
                        sol_collect[adjadjID]  = Ue[adjadjID]->getVal(0,0);
                        
                        
                        int nadjadjadj = iee_vec->i_map[adjadjID].size();
                        
                        for(int c=0;c<nadjadjadj;c++)
                        {
                            int adjadjadjID   = iee_vec->i_map[adjadjID][c];
                            
                            if(vrt_collect.find(adjadjadjID)==vrt_collect.end() && adjadjadjID>=Nel)
                            {
                                int NvPEladjadjadj = ien_part_map->i_map[adjadjadjID].size();
                                
                                double* Padjadjadj = new double[NvPEladjadjadj*3];

                                for(int l=0;l<NvPEladjadj;l++)
                                {
                                    int glob_vid     = ien_part_map->i_map[adjadjadjID][l];
                                    loc_vid          = gV2lV[glob_vid];
                                    Padjadjadj[k*3+0] = LocalVs[loc_vid]->x;
                                    Padjadjadj[k*3+1] = LocalVs[loc_vid]->y;
                                    Padjadjadj[k*3+2] = LocalVs[loc_vid]->z;
                                }

                                Vert* Vadjadjadj = ComputeCentroidCoord(Padjadj,NvPEladjadjadj);
                                
                                delete[] Padjadjadj;
                                
                                vrt_collect[adjadjID]  = Vadjadjadj;
                                sol_collect[adjadjID]  = Ue[adjadjadjID]->getVal(0,0);
                                
                            }
                            if(vrt_collect.find(adjadjadjID)==vrt_collect.end() && adjadjadjID>=Nel)
                            {
                                int fid    = ief_part_map->i_map[adjadjID][c];
                                int NvPerF = ifn_vec->i_map[fid].size();
                                
                                Vert* Vc = new Vert;
                                Vc->x    = 0.0;
                                Vc->y    = 0.0;
                                Vc->z    = 0.0;
                                
                                for(int s=0;s<NvPerF;s++)
                                {
                                    int gvid = ifn_vec->i_map[fid][s];
                                    int lvid = gV2lV[gvid];
                                    
                                    Vc->x = Vc->x+LocalVs[lvid]->x;
                                    Vc->y = Vc->y+LocalVs[lvid]->y;
                                    Vc->z = Vc->z+LocalVs[lvid]->z;
                                    
                                    Vert* V = new Vert;
                                    V->x    = LocalVs[lvid]->x;
                                    V->y    = LocalVs[lvid]->y;
                                    V->z    = LocalVs[lvid]->z;
                                    face.push_back(V);
                                }

                                Vc->x = Vc->x/NvPerF;
                                Vc->y = Vc->y/NvPerF;
                                Vc->z = Vc->z/NvPerF;
                                
                                Vec3D* r0 = new Vec3D;
                                r0->c0 = (Vc->x-Vijk->x);
                                r0->c1 = (Vc->y-Vijk->y);
                                r0->c2 = (Vc->z-Vijk->z);
                                            
                                
                                if(NvPerF==3)
                                {
                                    v0->c0 = face[1]->x-face[0]->x;
                                    v0->c1 = face[1]->y-face[0]->y;
                                    v0->c2 = face[1]->z-face[0]->z;

                                    v1->c0 = face[2]->x-face[0]->x;
                                    v1->c1 = face[2]->y-face[0]->y;
                                    v1->c2 = face[2]->z-face[0]->z;
                                    
                                    n0 = ComputeSurfaceNormal(v0,v1);
                                    orient0   = DotVec3D(r0,n0);

                                    if(orient0<0.0)
                                    {
                                        NegateVec3D(n0);
                                    }

                                    rdotn = DotVec3D(r0,n0);
                                    //delete rf;
                                }

                                if(NvPerF==4)
                                {
                                    v0->c0 = face[1]->x-face[0]->x;
                                    v0->c1 = face[1]->y-face[0]->y;
                                    v0->c2 = face[1]->z-face[0]->z;

                                    v1->c0 = face[3]->x-face[0]->x;
                                    v1->c1 = face[3]->y-face[0]->y;
                                    v1->c2 = face[3]->z-face[0]->z;
                                    
                                    n0        = ComputeSurfaceNormal(v0,v1);
                                    orient0   = DotVec3D(r0,n0);
                                    
                                    if(orient0<0.0)
                                    {
                                        NegateVec3D(n0);
                                    }
                                    
                                    double rdotn = DotVec3D(r0,n0);
                                }
                                
                                Vec3D* reflect = new Vec3D;
                                reflect->c0 = r0->c0-2.0*(rdotn)*n0->c0;
                                reflect->c1 = r0->c1-2.0*(rdotn)*n0->c1;
                                reflect->c2 = r0->c2-2.0*(rdotn)*n0->c2;
                                
                                Vc->x = Vc->x - reflect->c0;
                                Vc->y = Vc->y - reflect->c1;
                                Vc->z = Vc->z - reflect->c2;
                                
                                double vgx = Vc->x;
                                double vgy = Vc->y;
                                double vgz = Vc->z;
                                double r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                                double Utje = 0.1*tanh(50*(r-0.5))+1.0;
                                
                                double Utje2 = gbMap[adjadjadjID];

                                vrt_collect[adjadjadjID]  = Vc;
                                sol_collect[adjadjadjID]  = Utje;
                                face.clear();
                                cntbnd++;
                            }
                        }
                    }
                    
                    if(vrt_collect.find(adjadjID)==vrt_collect.end() && adjadjID>=Nel)
                    {
                        int fid    = ief_part_map->i_map[adjID][k];
                        int NvPerF = ifn_vec->i_map[fid].size();
                        
                        Vert* Vc = new Vert;
                        Vc->x    = 0.0;
                        Vc->y    = 0.0;
                        Vc->z    = 0.0;
                        
                        for(int s=0;s<NvPerF;s++)
                        {
                            int gvid = ifn_vec->i_map[fid][s];
                            int lvid = gV2lV[gvid];
                            
                            Vc->x = Vc->x+LocalVs[lvid]->x;
                            Vc->y = Vc->y+LocalVs[lvid]->y;
                            Vc->z = Vc->z+LocalVs[lvid]->z;
                            
                            Vert* V = new Vert;
                            V->x    = LocalVs[lvid]->x;
                            V->y    = LocalVs[lvid]->y;
                            V->z    = LocalVs[lvid]->z;
                            face.push_back(V);
                        }

                        Vc->x = Vc->x/NvPerF;
                        Vc->y = Vc->y/NvPerF;
                        Vc->z = Vc->z/NvPerF;
                        
                        Vec3D* r0 = new Vec3D;
                        r0->c0 = (Vc->x-Vijk->x);
                        r0->c1 = (Vc->y-Vijk->y);
                        r0->c2 = (Vc->z-Vijk->z);
                                    
                        
                        if(NvPerF==3)
                        {
                            v0->c0 = face[1]->x-face[0]->x;
                            v0->c1 = face[1]->y-face[0]->y;
                            v0->c2 = face[1]->z-face[0]->z;

                            v1->c0 = face[2]->x-face[0]->x;
                            v1->c1 = face[2]->y-face[0]->y;
                            v1->c2 = face[2]->z-face[0]->z;
                            
                            n0 = ComputeSurfaceNormal(v0,v1);
                            orient0   = DotVec3D(r0,n0);

                            if(orient0<0.0)
                            {
                                NegateVec3D(n0);
                            }

                            rdotn = DotVec3D(r0,n0);
                            //delete rf;
                        }

                        if(NvPerF==4)
                        {
                            v0->c0 = face[1]->x-face[0]->x;
                            v0->c1 = face[1]->y-face[0]->y;
                            v0->c2 = face[1]->z-face[0]->z;

                            v1->c0 = face[3]->x-face[0]->x;
                            v1->c1 = face[3]->y-face[0]->y;
                            v1->c2 = face[3]->z-face[0]->z;
                            
                            n0        = ComputeSurfaceNormal(v0,v1);
                            orient0   = DotVec3D(r0,n0);
                            
                            if(orient0<0.0)
                            {
                                NegateVec3D(n0);
                            }
                            
                            double rdotn = DotVec3D(r0,n0);
                        }
                        
                        Vec3D* reflect = new Vec3D;
                        reflect->c0 = r0->c0-2.0*(rdotn)*n0->c0;
                        reflect->c1 = r0->c1-2.0*(rdotn)*n0->c1;
                        reflect->c2 = r0->c2-2.0*(rdotn)*n0->c2;
                        
                        Vc->x = Vc->x - reflect->c0;
                        Vc->y = Vc->y - reflect->c1;
                        Vc->z = Vc->z - reflect->c2;
                        
                        double vgx = Vc->x;
                        double vgy = Vc->y;
                        double vgz = Vc->z;
                        double r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                        double Utje = 0.1*tanh(50*(r-0.5))+1.0;
                        
                        double Utje2 = gbMap[adjadjID];

                        vrt_collect[adjadjID]  = Vc;
                        sol_collect[adjadjID]  = Utje;
                        face.clear();
                        cntbnd++;
                    }
                    
                }
                
                
            }
            if(vrt_collect.find(adjid)==vrt_collect.end() && adjid>=Nel)
            {
                int fid    = ief_part_map->i_map[elID][j];
                int NvPerF = ifn_vec->i_map[fid].size();
                
                Vert* Vc = new Vert;
                Vc->x    = 0.0;
                Vc->y    = 0.0;
                Vc->z    = 0.0;
                
                for(int s=0;s<NvPerF;s++)
                {
                    int gvid = ifn_vec->i_map[fid][s];
                    int lvid = gV2lV[gvid];
                    
                    Vc->x = Vc->x+LocalVs[lvid]->x;
                    Vc->y = Vc->y+LocalVs[lvid]->y;
                    Vc->z = Vc->z+LocalVs[lvid]->z;
                    
                    Vert* V = new Vert;
                    V->x    = LocalVs[lvid]->x;
                    V->y    = LocalVs[lvid]->y;
                    V->z    = LocalVs[lvid]->z;
                    face.push_back(V);
                }

                Vc->x = Vc->x/NvPerF;
                Vc->y = Vc->y/NvPerF;
                Vc->z = Vc->z/NvPerF;
                
                Vec3D* r0 = new Vec3D;
                r0->c0 = (Vc->x-Vijk->x);
                r0->c1 = (Vc->y-Vijk->y);
                r0->c2 = (Vc->z-Vijk->z);
                            
                
                if(NvPerF==3)
                {
                    v0->c0 = face[1]->x-face[0]->x;
                    v0->c1 = face[1]->y-face[0]->y;
                    v0->c2 = face[1]->z-face[0]->z;

                    v1->c0 = face[2]->x-face[0]->x;
                    v1->c1 = face[2]->y-face[0]->y;
                    v1->c2 = face[2]->z-face[0]->z;
                    
                    n0 = ComputeSurfaceNormal(v0,v1);
                    orient0   = DotVec3D(r0,n0);

                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }

                    rdotn = DotVec3D(r0,n0);
                    //delete rf;
                }

                if(NvPerF==4)
                {
                    v0->c0 = face[1]->x-face[0]->x;
                    v0->c1 = face[1]->y-face[0]->y;
                    v0->c2 = face[1]->z-face[0]->z;

                    v1->c0 = face[3]->x-face[0]->x;
                    v1->c1 = face[3]->y-face[0]->y;
                    v1->c2 = face[3]->z-face[0]->z;
                    
                    n0        = ComputeSurfaceNormal(v0,v1);
                    orient0   = DotVec3D(r0,n0);
                    
                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }
                    
                    double rdotn = DotVec3D(r0,n0);
                }
                
                Vec3D* reflect = new Vec3D;
                reflect->c0 = r0->c0-2.0*(rdotn)*n0->c0;
                reflect->c1 = r0->c1-2.0*(rdotn)*n0->c1;
                reflect->c2 = r0->c2-2.0*(rdotn)*n0->c2;
                
                Vc->x = Vc->x - reflect->c0;
                Vc->y = Vc->y - reflect->c1;
                Vc->z = Vc->z - reflect->c2;
                
//                double rtje = sqrt(Vc->x*Vc->x+(Vc->y-0.5)*(Vc->y-0.5)+(Vc->z-0.5)*(Vc->z-0.5));
//                double Utje = 0.1*tanh(50*(rtje-0.5))+1.0;
                
                //double Utje = gbMap[adjid];
                
                double vgx = Vc->x;
                double vgy = Vc->y;
                double vgz = Vc->z;
                double r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                double Utje = 0.1*tanh(50*(r-0.5))+1.0;
                
                double Utje2 = gbMap[adjid];

                vrt_collect[adjid]  = Vc;
                sol_collect[adjid]  = Utje;
                face.clear();
                cntbnd++;
   
            }
        }
        
        if(vrt_collect.size() > 9)
        {
            int Ndata = vrt_collect.size();
            Array<double>* Vrt_T = new Array<double>(9,Ndata);
            Array<double>* Vrt   = new Array<double>(Ndata,9);
            Array<double>* bvec  = new Array<double>(Ndata,1);

            for(int q=0;q<Ndata;q++)
            {
                for(int g=0;g<9;g++)
                {
                    Vrt_T->setVal(g,q,0.0);
                    Vrt->setVal(q,g,0.0);
                }
            }

            std::map<int,Vert*>::iterator vit;
            int te = 0;

            double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;

            for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
            {
                double di = sqrt((vit->second->x-Vijk->x)*(vit->second->x-Vijk->x)+
                                 (vit->second->y-Vijk->y)*(vit->second->y-Vijk->y)+
                                 (vit->second->z-Vijk->z)*(vit->second->z-Vijk->z));

                a = (vit->second->x - Vijk->x);
                b = (vit->second->y - Vijk->y);
                c = (vit->second->z - Vijk->z);
                
                h00 = 0.5*a*a; h01 = 1.0*a*b; h02 = 1.0*a*c;
                h11 = 0.5*b*b; h12 = 1.0*b*c;
                h22 = 0.5*c*c;

                Vrt->setVal(te,0,1.0/di*a);
                Vrt->setVal(te,1,1.0/di*b);
                Vrt->setVal(te,2,1.0/di*c);

                Vrt->setVal(te,3, 1.0/di*h00);
                Vrt->setVal(te,4, 1.0/di*h01);
                Vrt->setVal(te,5, 1.0/di*h02);
                Vrt->setVal(te,6, 1.0/di*h11);
                Vrt->setVal(te,7, 1.0/di*h12);
                Vrt->setVal(te,8, 1.0/di*h22);

                double Udata = sol_collect[vit->first];

                bvec->setVal(te,0,(1.0/di)*(Udata-u_ijk));

                te++;
            }

            double* A_cm = new double[Ndata*9];
            for(int s=0;s<Ndata;s++)
            {
                for(int g=0;g<9;g++)
                {
                    A_cm[g*Ndata+s] = Vrt->getVal(s,g);
                }
            }

            Array<double>* x = SolveQR(A_cm,Ndata,9,bvec);


            dudx_map[elID] = x;

            delete[] A_cm;
            delete[] Pijk;
            delete Vrt_T;
            delete Vrt;
            delete bvec;


        }
        else
        {
            int Ndata = vrt_collect.size();
            
            
            
            Array<double>* Vrt_T = new Array<double>(3,Ndata);
            Array<double>* Vrt   = new Array<double>(Ndata,3);
            Array<double>* bvec  = new Array<double>(Ndata,1);

            for(int q=0;q<Ndata;q++)
            {
                for(int g=0;g<3;g++)
                {
                    Vrt_T->setVal(g,q,0.0);
                    Vrt->setVal(q,g,0.0);
                }
            }
            
            std::map<int,Vert*>::iterator vit;
            int te = 0;
            
            double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;
            
            for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
            {
                double di = sqrt((vit->second->x-Vijk->x)*(vit->second->x-Vijk->x)+
                                 (vit->second->y-Vijk->y)*(vit->second->y-Vijk->y)+
                                 (vit->second->z-Vijk->z)*(vit->second->z-Vijk->z));

                a = (vit->second->x - Vijk->x);
                b = (vit->second->y - Vijk->y);
                c = (vit->second->z - Vijk->z);
                
                
                Vrt->setVal(te,0,(1.0/di)*a);
                Vrt->setVal(te,1,(1.0/di)*b);
                Vrt->setVal(te,2,(1.0/di)*c);

                double Udata = sol_collect[vit->first];
                
                bvec->setVal(te,0,(1.0/di)*(Udata-u_ijk));
   
                te++;
            }

            double* A_cm = new double[Ndata*3];
            for(int s=0;s<Ndata;s++)
            {
                for(int g=0;g<3;g++)
                {
                    A_cm[g*Ndata+s] = Vrt->getVal(s,g);
                }
            }

            Array<double>* x = SolveQR(A_cm,Ndata,3,bvec);
            std::cout << "Warning:: not enough data points to reconstruct the higher gradients! Number of neigboring points is " << Ndata << " dudx = [" << x->getVal(0,0) << " " << x->getVal(1,0) << " " << x->getVal(2,0) << std::endl;
            dudx_map[elID] = x;
            
            delete[] A_cm;
            delete[] Pijk;
            delete Vrt_T;
            delete Vrt;
            delete bvec;
        }
//
        vrt_collect.clear();
        sol_collect.clear();
   }
    
   return dudx_map;
}
*/

std::map<int,Array<double>* > ComputedUdx_LSQ_HO_US3D(Partition* Pa, std::map<int,Array<double>* > Ue, Mesh_Topology* meshTopo, std::map<int,double> gbMap, MPI_Comm comm)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   std::vector<Vert*> LocalVs                 = Pa->getLocalVerts();
   std::map<int,std::vector<int> > gE2lV      = Pa->getGlobElem2LocVerts();
   std::map<int,std::vector<int> > gE2gF      = Pa->getglobElem2globFaces();
   std::map<int,int> gV2lV                    = Pa->getGlobalVert2LocalVert();
   std::map<int,int> gE2lE                    = Pa->getGlobalElement2LocalElement();
   std::vector<int> Loc_Elem                  = Pa->getLocElem();
   //std::map<int,std::vector<int> > scheme_E2V = meshTopo->getScheme_E2V();
   int nLoc_Elem                              = Loc_Elem.size();
    Array<int>* pg = Pa->getGlobalPartition();
    
   int Nel                      = Pa->getGlobalPartition()->getNrow();
   i_part_map*  ifn_vec         = Pa->getIFNpartmap();
    
   i_part_map* ief_part_map     = Pa->getIEFpartmap();
//   i_part_map* ief_adj_part_map = Pa->getIEFADJpartmap();
//   i_part_map* ief_adj2_part_map = Pa->getIEFADJ2partmap();
    
   i_part_map*  iee_vec         = Pa->getIEEpartmap();
    
//   i_part_map*  iee_adj_vec     = Pa->getIEEADJpartmap();
//   i_part_map*  iee_adj2_vec     = Pa->getIEEADJ2partmap();
    
   i_part_map* if_Nv_part_map   = Pa->getIF_Nvpartmap();
    

//   std::vector<std::vector<double> > iee_dist;
//   std::vector<double> dist;

    
   std::map<int,Array<double>* > dudx_map;
   double d;
   int loc_vid,adjID,elID;
   int cou = 0;
   Vert* Vc = new Vert;
   Vert* Vadj = new Vert;
   int lid = 0;
   double u_ijk, u_po;
    
   int el_contr = 1;
   int nadj_el  = 0;
   int cntbnd = 0;
   int cntbnd2 = 0;
   int cntbnd3 = 0;
   std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
   std::set<int> add2set_lay0;
   std::set<int> add2set_layt;
   std::map<int,Vert*> vrt_collect;
   std::map<int,double> sol_collect;
   std::set<int> add2set_lay1;
    std::vector<Vert*> face;
    double rdotn;
    Vec3D* n0 = new Vec3D;
    Vec3D* v0 = new Vec3D;
    Vec3D* v1 = new Vec3D;
   for(int i=0;i<nLoc_Elem;i++)
    {
        int bflip    = 0;
        int elID     = Loc_Elem[i];
        int NvPEl    = gE2lV[elID].size();
        
        double* Pijk = new double[NvPEl*3];
        
        for(int k=0;k<gE2lV[elID].size();k++)
        {
            loc_vid     = gE2lV[elID][k];
            Pijk[k*3+0] = LocalVs[loc_vid]->x;
            Pijk[k*3+1] = LocalVs[loc_vid]->y;
            Pijk[k*3+2] = LocalVs[loc_vid]->z;
        }
        
        Vert* Vijk   = ComputeCentroidCoord(Pijk,NvPEl);
        u_ijk        = Ue[elID]->getVal(0,0);

        delete[] Pijk;
        
        int nadj_tot   = LocElem2Nf[elID];
       
        for(int j=0;j<nadj_tot;j++)
        {
            int adjid   = iee_vec->i_map[elID][j];
               
            if(vrt_collect.find(adjid)==vrt_collect.end() && adjid<Nel)
            {
                int NvPEladj    = gE2lV[adjid].size();

                double* Padj = new double[NvPEladj*3];

                for(int k=0;k<gE2lV[adjid].size();k++)
                {
                    loc_vid     = gE2lV[adjid][k];
                    Padj[k*3+0] = LocalVs[loc_vid]->x;
                    Padj[k*3+1] = LocalVs[loc_vid]->y;
                    Padj[k*3+2] = LocalVs[loc_vid]->z;
                }

                Vert* Vadj = ComputeCentroidCoord(Padj,NvPEladj);
                
                delete[] Padj;
                
                vrt_collect[adjid]  = Vadj;
                sol_collect[adjid]  = Ue[adjid]->getVal(0,0);
                
                delete Vadj;
                
            }
            if(vrt_collect.find(adjid)==vrt_collect.end() && adjid>=Nel)
            {
                int fid    = ief_part_map->i_map[elID][j];
                int NvPerF = ifn_vec->i_map[fid].size();
                
                Vert* Vc = new Vert;
                Vc->x = 0.0;
                Vc->y = 0.0;
                Vc->z = 0.0;
                
                for(int s=0;s<NvPerF;s++)
                {
                    int gvid = ifn_vec->i_map[fid][s];
                    int lvid = gV2lV[gvid];

                    Vc->x = Vc->x+LocalVs[lvid]->x;
                    Vc->y = Vc->y+LocalVs[lvid]->y;
                    Vc->z = Vc->z+LocalVs[lvid]->z;
                    
                    Vert* V = new Vert;
                    V->x    = LocalVs[lvid]->x;
                    V->y    = LocalVs[lvid]->y;
                    V->z    = LocalVs[lvid]->z;
                    face.push_back(V);
                }

                Vc->x = Vc->x/NvPerF;
                Vc->y = Vc->y/NvPerF;
                Vc->z = Vc->z/NvPerF;
                
                Vec3D* r0 = new Vec3D;
                r0->c0 = (Vc->x-Vijk->x);
                r0->c1 = (Vc->y-Vijk->y);
                r0->c2 = (Vc->z-Vijk->z);
                
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
                
                Vc->x = Vc->x - reflect->c0;
                Vc->y = Vc->y - reflect->c1;
                Vc->z = Vc->z - reflect->c2;
                
//                double rtje = sqrt(Vc->x*Vc->x+(Vc->y-0.5)*(Vc->y-0.5)+(Vc->z-0.5)*(Vc->z-0.5));
//                double Utje = 0.1*tanh(50*(rtje-0.5))+1.0;
                
                double Utje_test = gbMap[adjid];
                double vgx = Vc->x;
                double vgy = Vc->y;
                double vgz = Vc->z;
                double r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                double Utje = 0.1*tanh(50*(r-0.5))+1.0;

                //double Utje    = 0.1*sin(50*Vc->x*Vc->z)+atan(0.1/((sin(5.0*Vc->y)-2.0*Vc->x*Vc->z)));
                vrt_collect[adjid]  = Vc;
                sol_collect[adjid]  = Utje;
                
                if(fabs(Utje_test-Utje)>1.0e-06)
                {
                    std::cout << "failed l1 :: " << Utje_test << " " << Utje << " " << fabs(Utje_test-Utje) << std::endl;
                }
                
                face.clear();
                cntbnd++;
   
            }
            
               
            if(iee_vec->i_map.find(adjid)!=iee_vec->i_map.end())
            {
                int n_adjid     = iee_vec->i_map[adjid].size();
                int NvPElnew    = gE2lV[adjid].size();
                double* Pijknew = new double[NvPElnew*3];
                
                for(int k=0;k<NvPElnew;k++)
                {
                    loc_vid     = gE2lV[adjid][k];
                    Pijknew[k*3+0] = LocalVs[loc_vid]->x;
                    Pijknew[k*3+1] = LocalVs[loc_vid]->y;
                    Pijknew[k*3+2] = LocalVs[loc_vid]->z;
                }
                
                Vert* Vijknew   = ComputeCentroidCoord(Pijknew,NvPElnew);
                
                delete[] Pijknew;
                
                for(int k=0;k<n_adjid;k++)
                {
                    int adjadj = iee_vec->i_map[adjid][k];
                    
                    if(vrt_collect.find(adjadj)==vrt_collect.end() && adjadj<Nel && adjadj!=elID)
                    {
                        int NvPEladjadj    = gE2lV[adjadj].size();
                        double* Padjadj    = new double[NvPEladjadj*3];

                        for(int k=0;k<gE2lV[adjadj].size();k++)
                        {
                            loc_vid            = gE2lV[adjadj][k];
                            Padjadj[k*3+0]     = LocalVs[loc_vid]->x;
                            Padjadj[k*3+1]     = LocalVs[loc_vid]->y;
                            Padjadj[k*3+2]     = LocalVs[loc_vid]->z;
                        }
                        
                        Vert* Vadjadj          = ComputeCentroidCoord(Padjadj,NvPEladjadj);

                        delete[] Padjadj;
                        
                        vrt_collect[adjadj] = Vadjadj;
                        sol_collect[adjadj] = Ue[adjadj]->getVal(0,0);
                        
                    }
                    if(vrt_collect.find(adjadj)==vrt_collect.end() && adjadj>=Nel && adjadj!=elID)
                    {
                        
                        int fid    = ief_part_map->i_map[adjid][k];
                        int NvPerF = ifn_vec->i_map[fid].size();
                        
                        Vert* Vc    = new Vert;
                        Vc->x       = 0.0;
                        Vc->y       = 0.0;
                        Vc->z       = 0.0;
                        
                        for(int s=0;s<NvPerF;s++)
                        {
                            int gvid = ifn_vec->i_map[fid][s];
                            int lvid = gV2lV[gvid];

                            Vc->x = Vc->x+LocalVs[lvid]->x;
                            Vc->y = Vc->y+LocalVs[lvid]->y;
                            Vc->z = Vc->z+LocalVs[lvid]->z;
                            
                            Vert* V = new Vert;
                            V->x    = LocalVs[lvid]->x;
                            V->y    = LocalVs[lvid]->y;
                            V->z    = LocalVs[lvid]->z;
                            face.push_back(V);
                        }

                        Vc->x = Vc->x/NvPerF;
                        Vc->y = Vc->y/NvPerF;
                        Vc->z = Vc->z/NvPerF;
                        
                        Vec3D* r0 = new Vec3D;
                        r0->c0 = (Vc->x-Vijknew->x);
                        r0->c1 = (Vc->y-Vijknew->y);
                        r0->c2 = (Vc->z-Vijknew->z);
                        
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
                        
                        Vc->x = Vc->x - reflect->c0;
                        Vc->y = Vc->y - reflect->c1;
                        Vc->z = Vc->z - reflect->c2;
                        
                        double Utje_test = gbMap[adjadj];
                        double vgx = Vc->x;
                        double vgy = Vc->y;
                        double vgz = Vc->z;
                        double r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                        double Utje = 0.1*tanh(50*(r-0.5))+1.0;

                        vrt_collect[adjadj]  = Vc;
                        sol_collect[adjadj]  = Utje;
                        
                        if(fabs(Utje_test-Utje)>1.0e-06)
                        {
                            std::cout << "failed l2 :: " << Utje_test << " " << Utje << " " << fabs(Utje_test-Utje) << " " << adjadj << std::endl;
                        }

                        
                        cntbnd++;
                        
                        face.clear();
                    }
                    
                    if(iee_vec->i_map.find(adjadj)!=iee_vec->i_map.end())
                    {
                        int n_adjadj = iee_vec->i_map[adjadj].size();
                        int NvPElnewnew    = gE2lV[adjadj].size();
                        double* Pijknewnew = new double[NvPElnewnew*3];
                        
                        for(int k=0;k<NvPElnewnew;k++)
                        {
                            loc_vid     = gE2lV[adjadj][k];
                            Pijknewnew[k*3+0] = LocalVs[loc_vid]->x;
                            Pijknewnew[k*3+1] = LocalVs[loc_vid]->y;
                            Pijknewnew[k*3+2] = LocalVs[loc_vid]->z;
                            
                            //std::cout << "NvPElnewnew " << NvPElnewnew << " " << loc_vid << " " << Pijknewnew[k*3+0] << " " << Pijknewnew[k*3+1] << " " << Pijknewnew[k*3+2] << std::endl;
                        }
                        
                        Vert* Vijknewnew   = ComputeCentroidCoord(Pijknewnew,NvPElnewnew);
                        
                        delete[] Pijknewnew;
                        
                        for(int k=0;k<n_adjadj;k++)
                        {
                            int adjadjadj = iee_vec->i_map[adjadj][k];
                            
                            if(vrt_collect.find(adjadjadj)==vrt_collect.end() && adjadjadj<Nel && adjadjadj!=elID)
                            {
                                int NvPEladjadjadj    = gE2lV[adjadjadj].size();
                                double* Padjadjadj = new double[NvPEladjadjadj*3];

                                for(int k=0;k<gE2lV[adjadjadj].size();k++)
                                {
                                    loc_vid               = gE2lV[adjadjadj][k];
                                    Padjadjadj[k*3+0]     = LocalVs[loc_vid]->x;
                                    Padjadjadj[k*3+1]     = LocalVs[loc_vid]->y;
                                    Padjadjadj[k*3+2]     = LocalVs[loc_vid]->z;
                                }
                                
                                Vert* Vadjadjadj          = ComputeCentroidCoord(Padjadjadj,NvPEladjadjadj);

                                delete[] Padjadjadj;
                                
                                vrt_collect[adjadjadj] = Vadjadjadj;
                                sol_collect[adjadjadj] = Ue[adjadjadj]->getVal(0,0);
                                
                            }
                            if(vrt_collect.find(adjadjadj)==vrt_collect.end() && adjadjadj>=Nel && adjadjadj!=elID)
                            {
                                
                                int fid    = ief_part_map->i_map[adjadj][k];
                                int NvPerF = ifn_vec->i_map[fid].size();
                                
                                Vert* Vc = new Vert;
                                Vc->x = 0.0;
                                Vc->y = 0.0;
                                Vc->z = 0.0;
                                
                                for(int s=0;s<NvPerF;s++)
                                {
                                    int gvid = ifn_vec->i_map[fid][s];
                                    int lvid = gV2lV[gvid];

                                    Vc->x = Vc->x+LocalVs[lvid]->x;
                                    Vc->y = Vc->y+LocalVs[lvid]->y;
                                    Vc->z = Vc->z+LocalVs[lvid]->z;
                                    
                                    Vert* V = new Vert;
                                    V->x    = LocalVs[lvid]->x;
                                    V->y    = LocalVs[lvid]->y;
                                    V->z    = LocalVs[lvid]->z;
                                    face.push_back(V);
                                }

                                Vc->x = Vc->x/NvPerF;
                                Vc->y = Vc->y/NvPerF;
                                Vc->z = Vc->z/NvPerF;
                                
                                Vec3D* r0 = new Vec3D;
                                r0->c0 = (Vc->x-Vijknewnew->x);
                                r0->c1 = (Vc->y-Vijknewnew->y);
                                r0->c2 = (Vc->z-Vijknewnew->z);
                                
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
                                
                                Vc->x = Vc->x - reflect->c0;
                                Vc->y = Vc->y - reflect->c1;
                                Vc->z = Vc->z - reflect->c2;
                                
                                //double Utje = gbMap[adjadjadj];
                                double Utje_test = gbMap[adjadjadj];
                                
                                
                                
                                double vgx = Vc->x;
                                double vgy = Vc->y;
                                double vgz = Vc->z;
                                double r = sqrt(vgx*vgx+(vgy-0.5)*(vgy-0.5)+(vgz-0.5)*(vgz-0.5));
                                double Utje = 0.1*tanh(50*(r-0.5))+1.0;
                                
                                vrt_collect[adjadjadj]  = Vc;
                                sol_collect[adjadjadj]  = Utje;
                                
                                if(fabs(Utje_test-Utje)>1.0e-06)
                                {
                                    std::cout << "failed l3 :: " << Utje_test << " " << Utje << " " << fabs(Utje_test-Utje) << std::endl;
                                }
                                face.clear();
                                
                                cntbnd++;
                            }
                        }
                    }
                }
            }
        }
        
        
        
        
        if(vrt_collect.size() > 9)
        {
            int Ndata = vrt_collect.size();
            Array<double>* Vrt_T = new Array<double>(9,Ndata);
            Array<double>* Vrt   = new Array<double>(Ndata,9);
            Array<double>* bvec  = new Array<double>(Ndata,1);

            for(int q=0;q<Ndata;q++)
            {
                for(int g=0;g<9;g++)
                {
                    Vrt_T->setVal(g,q,0.0);
                    Vrt->setVal(q,g,0.0);
                }
            }

            std::map<int,Vert*>::iterator vit;
            int te = 0;

            double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;
//            Vert->x,Vert->y,Vert-z;
//            std::map<Vert*,double>
            for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
            {
                double di = sqrt((vit->second->x-Vijk->x)*(vit->second->x-Vijk->x)+
                                 (vit->second->y-Vijk->y)*(vit->second->y-Vijk->y)+
                                 (vit->second->z-Vijk->z)*(vit->second->z-Vijk->z));

                a = (vit->second->x - Vijk->x);
                b = (vit->second->y - Vijk->y);
                c = (vit->second->z - Vijk->z);

//                h00 = 0.5*a*a; h01 = 0.5*a*b; h02 = 0.5*a*c;
//                h10 = 0.5*b*a; h11 = 0.5*b*b; h12 = 0.5*b*c;
//                h20 = 0.5*c*a; h21 = 0.5*c*b; h22 = 0.5*c*c;
//
//                Vrt->setVal(te,0,1.0/di*a);
//                Vrt->setVal(te,1,1.0/di*b);
//                Vrt->setVal(te,2,1.0/di*c);
//
//                Vrt->setVal(te,3, 1.0/di*h00);
//                Vrt->setVal(te,4, 1.0/di*h01);
//                Vrt->setVal(te,5, 1.0/di*h02);
//                Vrt->setVal(te,6, 1.0/di*h10);
//                Vrt->setVal(te,7, 1.0/di*h11);
//                Vrt->setVal(te,8, 1.0/di*h12);
//                Vrt->setVal(te,9, 1.0/di*h20);
//                Vrt->setVal(te,10, 1.0/di*h21);
//                Vrt->setVal(te,11, 1.0/di*h22);
                
                h00 = 0.5*a*a; h01 = 1.0*a*b; h02 = 1.0*a*c;
                h11 = 0.5*b*b; h12 = 1.0*b*c;
                h22 = 0.5*c*c;

                Vrt->setVal(te,0,1.0/di*a);
                Vrt->setVal(te,1,1.0/di*b);
                Vrt->setVal(te,2,1.0/di*c);

                Vrt->setVal(te,3, 1.0/di*h00);
                Vrt->setVal(te,4, 1.0/di*h01);
                Vrt->setVal(te,5, 1.0/di*h02);
                Vrt->setVal(te,6, 1.0/di*h11);
                Vrt->setVal(te,7, 1.0/di*h12);
                Vrt->setVal(te,8, 1.0/di*h22);

                double Udata = sol_collect[vit->first];
                //std::cout << " check " << Ndata << " " << Udata << " " << di << std::endl;
                bvec->setVal(te,0,(1.0/di)*(Udata-u_ijk));

                te++;
            }

            double* A_cm = new double[Ndata*9];
            for(int s=0;s<Ndata;s++)
            {
                for(int g=0;g<9;g++)
                {
                    A_cm[g*Ndata+s] = Vrt->getVal(s,g);
                }
            }

            Array<double>* x = SolveQR(A_cm,Ndata,9,bvec);


            dudx_map[elID] = x;

            delete[] A_cm;
            //delete[] Pijk;
            delete Vrt_T;
            delete Vrt;
            delete bvec;


        }
        else
        {
            int Ndata = vrt_collect.size();
            
            std::cout << "Warning:: not enough data points to reconstruct the gradient! Number of neigboring points is " << Ndata << std::endl;
            
            Array<double>* Vrt_T = new Array<double>(3,Ndata);
            Array<double>* Vrt   = new Array<double>(Ndata,3);
            Array<double>* bvec  = new Array<double>(Ndata,1);

            for(int q=0;q<Ndata;q++)
            {
                for(int g=0;g<3;g++)
                {
                    Vrt_T->setVal(g,q,0.0);
                    Vrt->setVal(q,g,0.0);
                }
            }
            
            std::map<int,Vert*>::iterator vit;
            int te = 0;
            
            double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;
            
            for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
            {
                double di = sqrt((vit->second->x-Vijk->x)*(vit->second->x-Vijk->x)+
                                 (vit->second->y-Vijk->y)*(vit->second->y-Vijk->y)+
                                 (vit->second->z-Vijk->z)*(vit->second->z-Vijk->z));

                a = (vit->second->x - Vijk->x);
                b = (vit->second->y - Vijk->y);
                c = (vit->second->z - Vijk->z);
                
                
                Vrt->setVal(te,0,(1.0/di)*a);
                Vrt->setVal(te,1,(1.0/di)*b);
                Vrt->setVal(te,2,(1.0/di)*c);

                double Udata = sol_collect[vit->first];
                
                bvec->setVal(te,0,(1.0/di)*(Udata-u_ijk));
   
                te++;
            }

            double* A_cm = new double[Ndata*3];
            for(int s=0;s<Ndata;s++)
            {
                for(int g=0;g<3;g++)
                {
                    A_cm[g*Ndata+s] = Vrt->getVal(s,g);
                }
            }

            Array<double>* x = SolveQR(A_cm,Ndata,3,bvec);
            
            dudx_map[elID] = x;
            
            delete[] A_cm;
            //delete[] Pijk;
            delete Vrt_T;
            delete Vrt;
            delete bvec;
        }
//
        vrt_collect.clear();
        sol_collect.clear();
   }
    
    
   return dudx_map;
}




std::map<int,Array<double>* > ComputedUdx_LSQ_US3D(Partition* Pa,
                                                   std::map<int,Array<double>* > U,
                                                   Mesh_Topology* meshTopo,
                                                   std::map<int,double> gbMap,
                                                   MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<Vert*> LocalVs            = Pa->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV = Pa->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
    std::map<int,int> gV2lV               = Pa->getGlobalVert2LocalVert();
    std::map<int,int> gE2lE               = Pa->getGlobalElement2LocalElement();
    std::vector<int> Loc_Elem             = Pa->getLocElem();
    int nLoc_Elem                         = Loc_Elem.size();
    
//    std::map<int,vector<Vec3D*> > rvector   = meshTopo->getRvectors();
//    std::map<int,vector<Vec3D*> > dxfxc     = meshTopo->getdXfXc();
//    std::map<int,vector<double> > dS        = meshTopo->getdS();
//    std::map<int,vector<double> > dr        = meshTopo->getdr();
//    std::map<int,double > vol               = meshTopo->getVol();
//    std::map<int,std::vector<Vert*> > vfvec = meshTopo->getVfacevector();
    
    std::map<int,Vert*> ghostvrts = meshTopo->getGhostVerts();
    std::map<int,std::vector<Vert*> > vfvec = meshTopo->getVfacevector();

    int Nel = Pa->getGlobalPartition()->getNrow();
    i_part_map*  ifn_vec         = Pa->getIFNpartmap();
    i_part_map* ief_part_map     = Pa->getIEFpartmap();
    i_part_map*  iee_vec         = Pa->getIEEpartmap();
    i_part_map*  ien_vec         = Pa->getIENpartmap();
    i_part_map* if_Nv_part_map   = Pa->getIF_Nvpartmap();
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;

    std::map<int,Array<double>* > dudx_map;
    double d;
    int loc_vid,adjID,elID;
    int cou = 0;
    Vert* Vc = new Vert;
    Vert* Vadj;
    int lid = 0;
    double u_ijk, u_po;
    //Array<double>* dudx = new Array<double>(nLoc_Elem,3);
    std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
    std::map<int,int> LocElem2Nv = Pa->getLocElem2Nv();
    std::vector<Vert*> face;
    double rdotn;
    Vec3D* v0 = new Vec3D;
    Vec3D* v1 = new Vec3D;
    Vec3D* n0 = new Vec3D;
    double orient0;
    int isZero = 0;
    int isNotZero = 0;
    for(int i=0;i<nLoc_Elem;i++)
    {
       int elID  = Loc_Elem[i];
       int NvPEl = LocElem2Nv[elID];
       int nadj  = LocElem2Nf[elID];
       
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
       double* Pijk = new double[NvPEl*3];
       for(int k=0;k<gE2lV[elID].size();k++)
       {
           loc_vid     = gE2lV[elID][k];
           Pijk[k*3+0] = LocalVs[loc_vid]->x;
           Pijk[k*3+1] = LocalVs[loc_vid]->y;
           Pijk[k*3+2] = LocalVs[loc_vid]->z;
       }
       
       Vert* Vijk   = ComputeCentroidCoord(Pijk,NvPEl);
       u_ijk        = U[elID]->getVal(0,0);
       int t        = 0;
       
       for(int j=0;j<nadj;j++)
       {
           int adjID = iee_vec->i_map[elID][j];

           if(adjID<Nel)
           {
               int nVadj = ien_vec->i_map[adjID].size();

               double* Padj = new double[nVadj*3];

               for(int k=0;k<nVadj;k++)
               {
                   int global_vid  = ien_vec->i_map[adjID][k];
                   loc_vid         = gV2lV[global_vid];
                   Padj[k*3+0] = LocalVs[loc_vid]->x;
                   Padj[k*3+1] = LocalVs[loc_vid]->y;
                   Padj[k*3+2] = LocalVs[loc_vid]->z;
               }
               
               Vadj = ComputeCentroidCoord(Padj,gE2lV[adjID].size());
               
               d = sqrt((Vadj->x-Vijk->x)*(Vadj->x-Vijk->x)+
                        (Vadj->y-Vijk->y)*(Vadj->y-Vijk->y)+
                        (Vadj->z-Vijk->z)*(Vadj->z-Vijk->z));
               
               Vrt->setVal(t,0,(1.0/d)*(Vadj->x-Vijk->x));
               Vrt->setVal(t,1,(1.0/d)*(Vadj->y-Vijk->y));
               Vrt->setVal(t,2,(1.0/d)*(Vadj->z-Vijk->z));
               
               
               u_po = U[adjID]->getVal(0,0);

               b->setVal(t,0,(1.0/d)*(u_po-u_ijk));
               delete Vadj;
               dist.push_back(d);
               delete[] Padj;
               t++;
               
           }
           else
           {
               
               int fid    = ief_part_map->i_map[elID][j];
               int NvPerF = if_Nv_part_map->i_map[fid][0];

               Vc->x = 0.0;
               Vc->y = 0.0;
               Vc->z = 0.0;

               for(int s=0;s<NvPerF;s++)
               {
                   int gvid = ifn_vec->i_map[fid][s];
                   int lvid = gV2lV[gvid];

                   Vc->x = Vc->x+LocalVs[lvid]->x;
                   Vc->y = Vc->y+LocalVs[lvid]->y;
                   Vc->z = Vc->z+LocalVs[lvid]->z;

                   Vert* V = new Vert;
                   V->x    = LocalVs[lvid]->x;
                   V->y    = LocalVs[lvid]->y;
                   V->z    = LocalVs[lvid]->z;

                   face.push_back(V);
               }

               Vc->x = Vc->x/NvPerF;
               Vc->y = Vc->y/NvPerF;
               Vc->z = Vc->z/NvPerF;

               Vec3D* r0 = new Vec3D;
               r0->c0 = (Vc->x-Vijk->x);
               r0->c1 = (Vc->y-Vijk->y);
               r0->c2 = (Vc->z-Vijk->z);

               if(NvPerF==3) // triangle
               {
                   v0->c0 = face[1]->x-face[0]->x;
                   v0->c1 = face[1]->y-face[0]->y;
                   v0->c2 = face[1]->z-face[0]->z;

                   v1->c0 = face[2]->x-face[0]->x;
                   v1->c1 = face[2]->y-face[0]->y;
                   v1->c2 = face[2]->z-face[0]->z;

                   n0 = ComputeSurfaceNormal(v0,v1);
                   orient0   = DotVec3D(r0,n0);

                   if(orient0<0.0)
                   {
                       NegateVec3D(n0);
                   }

                   rdotn = DotVec3D(r0,n0);
               }

               if(NvPerF==4) // triangle
               {
                   v0->c0 = face[1]->x-face[0]->x;
                   v0->c1 = face[1]->y-face[0]->y;
                   v0->c2 = face[1]->z-face[0]->z;

                   v1->c0 = face[3]->x-face[0]->x;
                   v1->c1 = face[3]->y-face[0]->y;
                   v1->c2 = face[3]->z-face[0]->z;

                   n0 = ComputeSurfaceNormal(v0,v1);
                   orient0   = DotVec3D(r0,n0);

                   if(orient0<0.0)
                   {
                       NegateVec3D(n0);
                   }

                   rdotn = DotVec3D(r0,n0);
               }

               Vec3D* reflect = new Vec3D;
               reflect->c0 = r0->c0-2.0*(rdotn)*n0->c0;
               reflect->c1 = r0->c1-2.0*(rdotn)*n0->c1;
               reflect->c2 = r0->c2-2.0*(rdotn)*n0->c2;

               Vert* Vf = new Vert;
               Vf->x = Vc->x;
               Vf->y = Vc->y;
               Vf->z = Vc->z;
               
               Vc->x = Vc->x - reflect->c0;
               Vc->y = Vc->y - reflect->c1;
               Vc->z = Vc->z - reflect->c2;

               Vert* Vc2 = ghostvrts[adjID];
               Vert* Vface_compare = vfvec[elID][j];
               
               
               double sum_error = sqrt((Vc->x-Vc2->x)*(Vc->x-Vc2->x)+(Vc->y-Vc2->y)*(Vc->y-Vc2->y)+(Vc->z-Vc2->z)*(Vc->z-Vc2->z));

               double sum_error_face = sqrt((Vf->x-Vface_compare->x)*(Vf->x-Vface_compare->x)+
                                            (Vf->y-Vface_compare->y)*(Vf->y-Vface_compare->y)+
                                            (Vf->z-Vface_compare->z)*(Vf->z-Vface_compare->z));
//
               if(sum_error_face > 0.0)
               {
                   //std::cout << sum_error_face << std::endl;
                   std::cout << Nel << " " << adjID << " (" << Vf->x << " " << Vf->y << " " << Vf->z << ") " << " (" << Vc->x << " " << Vc->y << " " << Vc->z << ") "<< " (" << Vface_compare->x << " " << Vface_compare->y << " " << Vface_compare->z << ") "<< std::endl;
                   isNotZero++;
               }
               if(sum_error_face == 0.0)
               {
                   isZero++;
               }
               
               double Utje = gbMap[adjID];
               
               d = sqrt((Vc->x-Vijk->x)*(Vc->x-Vijk->x)+
                        (Vc->y-Vijk->y)*(Vc->y-Vijk->y)+
                        (Vc->z-Vijk->z)*(Vc->z-Vijk->z));
               
               u_po = Utje;
               
               Vrt->setVal(t,0,(1.0/d)*(Vc->x-Vijk->x));
               Vrt->setVal(t,1,(1.0/d)*(Vc->y-Vijk->y));
               Vrt->setVal(t,2,(1.0/d)*(Vc->z-Vijk->z));
               
               b->setVal(t,0,(1.0/d)*(u_po-u_ijk));
               
               //delete r0;
               face.clear();
               t++;
           }
           
           
           
      }
       
      double* A_cm = new double[nadj*3];
      for(int s=0;s<nadj;s++)
      {
          for(int j=0;j<3;j++)
          {
              A_cm[j*nadj+s] = Vrt->getVal(s,j);
          }
      }
       
       Array<double>* x = SolveQR(A_cm,nadj,3,b);

       
       dudx_map[elID] = x;
       delete[] A_cm;
       //delete x;
       delete[] Pijk;
       delete Vrt_T;
       delete Vrt;
       delete b;

       //iee_dist.push_back(dist);
       //dist.clear();
    }

    //std::cout << "isZero " << isZero <<   " isNotZero  " << isNotZero << std::endl;
    delete Vc;

    return dudx_map;
}




std::map<int,Array<double>* > ComputedUdx_LSQ_US3D_LargeStencil(Partition* Pa, std::map<int,Array<double>* > Ue, Mesh_Topology* meshTopo, std::map<int,double> gbMap, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<Vert*> LocalVs                 = Pa->getLocalVerts();
    std::map<int,std::vector<int> > gE2lV      = Pa->getGlobElem2LocVerts();
    std::map<int,std::vector<int> > gE2gF      = Pa->getglobElem2globFaces();
    std::map<int,int> gV2lV                    = Pa->getGlobalVert2LocalVert();
    std::map<int,int> gE2lE                    = Pa->getGlobalElement2LocalElement();
    std::vector<int> Loc_Elem                  = Pa->getLocElem();
    //std::map<int,std::vector<int> > scheme_E2V = meshTopo->getScheme_E2V();
    int nLoc_Elem                              = Loc_Elem.size();
     Array<int>* pg = Pa->getGlobalPartition();
     
    int Nel                      = Pa->getGlobalPartition()->getNrow();
    i_part_map*  ifn_vec         = Pa->getIFNpartmap();
     
    i_part_map* ief_part_map     = Pa->getIEFpartmap();
 //   i_part_map* ief_adj_part_map = Pa->getIEFADJpartmap();
 //   i_part_map* ief_adj2_part_map = Pa->getIEFADJ2partmap();
     
    i_part_map*  iee_vec         = Pa->getIEEpartmap();
 //   i_part_map*  iee_adj_vec     = Pa->getIEEADJpartmap();
 //   i_part_map*  iee_adj2_vec     = Pa->getIEEADJ2partmap();
     
    i_part_map* if_Nv_part_map   = Pa->getIF_Nvpartmap();
     

 //   std::vector<std::vector<double> > iee_dist;
 //   std::vector<double> dist;

     
    std::map<int,Array<double>* > dudx_map;
    double d;
    int loc_vid,adjID,elID;
    int cou = 0;
    Vert* Vc = new Vert;
    Vert* Vadj = new Vert;
    int lid = 0;
    double u_ijk, u_po;
     
    int el_contr = 1;
    int nadj_el  = 0;
    int cntbnd = 0;
    int cntbnd2 = 0;
    int cntbnd3 = 0;
    std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
    std::set<int> add2set_lay0;
    std::set<int> add2set_layt;
    std::map<int,Vert*> vrt_collect;
    std::map<int,double> sol_collect;
    std::set<int> add2set_lay1;
    for(int i=0;i<nLoc_Elem;i++)
     {
         int bflip    = 0;
         int elID     = Loc_Elem[i];
         int NvPEl    = gE2lV[elID].size();
         
         double* Pijk = new double[NvPEl*3];
         
         for(int k=0;k<gE2lV[elID].size();k++)
         {
             loc_vid     = gE2lV[elID][k];
             Pijk[k*3+0] = LocalVs[loc_vid]->x;
             Pijk[k*3+1] = LocalVs[loc_vid]->y;
             Pijk[k*3+2] = LocalVs[loc_vid]->z;
         }
         
         Vert* Vijk   = ComputeCentroidCoord(Pijk,NvPEl);
         u_ijk        = Ue[elID]->getVal(0,0);

         int nadj_tot   = LocElem2Nf[elID];
        
         for(int j=0;j<nadj_tot;j++)
         {
             int adjid   = iee_vec->i_map[elID][j];
                
             if(vrt_collect.find(adjid)==vrt_collect.end() && adjid<Nel)
             {
                 int NvPEladj    = gE2lV[adjid].size();

                 double* Padj = new double[NvPEladj*3];

                 for(int k=0;k<gE2lV[adjid].size();k++)
                 {
                     loc_vid     = gE2lV[adjid][k];
                     Padj[k*3+0] = LocalVs[loc_vid]->x;
                     Padj[k*3+1] = LocalVs[loc_vid]->y;
                     Padj[k*3+2] = LocalVs[loc_vid]->z;
                 }

                 Vert* Vadj = ComputeCentroidCoord(Padj,NvPEladj);
                 
                 delete[] Padj;
                 
                 vrt_collect[adjid]  = Vadj;
                 sol_collect[adjid]  = Ue[adjid]->getVal(0,0);
                 
                 delete Vadj;
                 
             }
             if(vrt_collect.find(adjid)==vrt_collect.end() && adjid>=Nel)
             {
                 int fid    = ief_part_map->i_map[elID][j];
                 int NvPerF = if_Nv_part_map->i_map[fid][0];
                 
                 Vert* Vc = new Vert;
                 
                 Vc->x = 0.0;
                 Vc->y = 0.0;
                 Vc->z = 0.0;
                 
                 for(int s=0;s<NvPerF;s++)
                 {
                     int gvid = ifn_vec->i_map[fid][s];
                     int lvid = gV2lV[gvid];

                     Vc->x = Vc->x+LocalVs[lvid]->x;
                     Vc->y = Vc->y+LocalVs[lvid]->y;
                     Vc->z = Vc->z+LocalVs[lvid]->z;
                 }

                 Vc->x = Vc->x/NvPerF;
                 Vc->y = Vc->y/NvPerF;
                 Vc->z = Vc->z/NvPerF;
                 
                 Vec3D* r0 = new Vec3D;
                 r0->c0 = (Vc->x-Vijk->x);
                 r0->c1 = (Vc->y-Vijk->y);
                 r0->c2 = (Vc->z-Vijk->z);
                 
                 NegateVec3D(r0);
                 
                 Vc->x = Vc->x + r0->c0;
                 Vc->y = Vc->y + r0->c1;
                 Vc->z = Vc->z + r0->c2;
                                  
                 double Utje = gbMap[adjid];
                 
                 //double Utje    = 0.1*sin(50*Vc->x*Vc->z)+atan(0.1/((sin(5.0*Vc->y)-2.0*Vc->x*Vc->z)));
                 
                 vrt_collect[adjid]  = Vc;
                 sol_collect[adjid]  = Utje;
                 
                 delete r0;
                 
                 cntbnd++;
    
             }
             
                
             if(iee_vec->i_map.find(adjid)!=iee_vec->i_map.end())
             {
                 int n_adjid = iee_vec->i_map[adjid].size();
                    
                 for(int k=0;k<n_adjid;k++)
                 {
                     int adjadj = iee_vec->i_map[adjid][k];
                     
                     if(vrt_collect.find(adjadj)==vrt_collect.end() && adjadj<Nel && adjadj!=elID)
                     {
                         int NvPEladjadj    = gE2lV[adjadj].size();
                         double* Padjadj = new double[NvPEladjadj*3];

                         for(int k=0;k<gE2lV[adjadj].size();k++)
                         {
                             loc_vid            = gE2lV[adjadj][k];
                             Padjadj[k*3+0]     = LocalVs[loc_vid]->x;
                             Padjadj[k*3+1]     = LocalVs[loc_vid]->y;
                             Padjadj[k*3+2]     = LocalVs[loc_vid]->z;
                         }
                         
                         Vert* Vadjadj          = ComputeCentroidCoord(Padjadj,NvPEladjadj);

                         delete[] Padjadj;
                         
                         vrt_collect[adjadj] = Vadjadj;
                         sol_collect[adjadj] = Ue[adjadj]->getVal(0,0);
                         
                     }
                     if(vrt_collect.find(adjadj)==vrt_collect.end() && adjadj>=Nel && adjadj!=elID)
                     {
                         
                         int fid    = ief_part_map->i_map[adjid][k];
                         int NvPerF = 3;
                         
                         Vert* Vc    = new Vert;
                         Vc->x       = 0.0;
                         Vc->y       = 0.0;
                         Vc->z       = 0.0;
                         
                         for(int s=0;s<NvPerF;s++)
                         {
                             int gvid = ifn_vec->i_map[fid][s];
                             int lvid = gV2lV[gvid];

                             Vc->x = Vc->x+LocalVs[lvid]->x;
                             Vc->y = Vc->y+LocalVs[lvid]->y;
                             Vc->z = Vc->z+LocalVs[lvid]->z;
                         }

                         Vc->x = Vc->x/NvPerF;
                         Vc->y = Vc->y/NvPerF;
                         Vc->z = Vc->z/NvPerF;
                         
                         Vec3D* r0 = new Vec3D;
                         r0->c0 = (Vc->x-Vijk->x);
                         r0->c1 = (Vc->y-Vijk->y);
                         r0->c2 = (Vc->z-Vijk->z);
                         
                         NegateVec3D(r0);
                         
                         Vc->x = Vc->x + r0->c0;
                         Vc->y = Vc->y + r0->c1;
                         Vc->z = Vc->z + r0->c2;
                         
                         double Utje = gbMap[adjadj];

                         vrt_collect[adjadj]  = Vc;
                         sol_collect[adjadj]  = Utje;
                         
                         delete r0;
                         
                         cntbnd++;
                     }
                     
                     if(iee_vec->i_map.find(adjadj)!=iee_vec->i_map.end())
                     {
                         int n_adjadj = iee_vec->i_map[adjadj].size();
                            
                         for(int k=0;k<n_adjadj;k++)
                         {
                             int adjadjadj = iee_vec->i_map[adjadj][k];
                             
                             if(vrt_collect.find(adjadjadj)==vrt_collect.end() && adjadjadj<Nel && adjadjadj!=elID)
                             {
                                 int NvPEladjadjadj    = gE2lV[adjadjadj].size();
                                 double* Padjadjadj = new double[NvPEladjadjadj*3];

                                 for(int k=0;k<gE2lV[adjadjadj].size();k++)
                                 {
                                     loc_vid            = gE2lV[adjadjadj][k];
                                     Padjadjadj[k*3+0]     = LocalVs[loc_vid]->x;
                                     Padjadjadj[k*3+1]     = LocalVs[loc_vid]->y;
                                     Padjadjadj[k*3+2]     = LocalVs[loc_vid]->z;
                                 }
                                 
                                 Vert* Vadjadjadj          = ComputeCentroidCoord(Padjadjadj,NvPEladjadjadj);

                                 delete[] Padjadjadj;
                                 
                                 vrt_collect[adjadjadj] = Vadjadjadj;
                                 sol_collect[adjadjadj] = Ue[adjadjadj]->getVal(0,0);
                                 
                             }
                             if(vrt_collect.find(adjadjadj)==vrt_collect.end() && adjadjadj>=Nel && adjadjadj!=elID)
                             {
                                 
                                 int fid    = ief_part_map->i_map[adjadj][k];
                                 int NvPerF = 3;
                                 
                                 Vert* Vc = new Vert;
                                 Vc->x = 0.0;
                                 Vc->y = 0.0;
                                 Vc->z = 0.0;
                                 
                                 for(int s=0;s<NvPerF;s++)
                                 {
                                     int gvid = ifn_vec->i_map[fid][s];
                                     int lvid = gV2lV[gvid];

                                     Vc->x = Vc->x+LocalVs[lvid]->x;
                                     Vc->y = Vc->y+LocalVs[lvid]->y;
                                     Vc->z = Vc->z+LocalVs[lvid]->z;
                                 }

                                 Vc->x = Vc->x/NvPerF;
                                 Vc->y = Vc->y/NvPerF;
                                 Vc->z = Vc->z/NvPerF;
                                 
                                 Vec3D* r0 = new Vec3D;
                                 r0->c0 = (Vc->x-Vijk->x);
                                 r0->c1 = (Vc->y-Vijk->y);
                                 r0->c2 = (Vc->z-Vijk->z);
                                 
                                 NegateVec3D(r0);
                                 
                                 Vc->x = Vc->x + r0->c0;
                                 Vc->y = Vc->y + r0->c1;
                                 Vc->z = Vc->z + r0->c2;
                                 
                                 double Utje = gbMap[adjadjadj];

                                 vrt_collect[adjadjadj]  = Vc;
                                 sol_collect[adjadjadj]  = Utje;
                                 
                                 delete r0;
                                 
                                 cntbnd++;
                             }
                         }
                     }
                 }
             }
         }
         
         int Ndata = vrt_collect.size();
         Array<double>* Vrt_T = new Array<double>(3,Ndata);
         Array<double>* Vrt   = new Array<double>(Ndata,3);
         Array<double>* bvec  = new Array<double>(Ndata,1);

         for(int q=0;q<Ndata;q++)
         {
             for(int g=0;g<3;g++)
             {
                 Vrt_T->setVal(g,q,0.0);
                 Vrt->setVal(q,g,0.0);
             }
         }

         std::map<int,Vert*>::iterator vit;
         int te = 0;

         double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;
//            Vert->x,Vert->y,Vert-z;
//            std::map<Vert*,double>
         for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
         {
             double di = sqrt((vit->second->x-Vijk->x)*(vit->second->x-Vijk->x)+
                              (vit->second->y-Vijk->y)*(vit->second->y-Vijk->y)+
                              (vit->second->z-Vijk->z)*(vit->second->z-Vijk->z));

             a = (vit->second->x - Vijk->x);
             b = (vit->second->y - Vijk->y);
             c = (vit->second->z - Vijk->z);
             
             Vrt->setVal(te,0,1.0/di*a);
             Vrt->setVal(te,1,1.0/di*b);
             Vrt->setVal(te,2,1.0/di*c);

             double Udata = sol_collect[vit->first];

             bvec->setVal(te,0,(1.0/di)*(Udata-u_ijk));

             te++;
         }

         double* A_cm = new double[Ndata*3];
         for(int s=0;s<Ndata;s++)
         {
             for(int g=0;g<3;g++)
             {
                 A_cm[g*Ndata+s] = Vrt->getVal(s,g);
             }
         }

         Array<double>* x = SolveQR(A_cm,Ndata,3,bvec);

         dudx_map[elID] = x;

         delete[] A_cm;
         delete[] Pijk;
         delete Vrt_T;
         delete Vrt;
         delete bvec;

         vrt_collect.clear();
         sol_collect.clear();
    }
     
     
    return dudx_map;
}



std::map<int,Array<double>* >  ComputedUdx_MGG(Partition* Pa, std::map<int,Array<double>*> U, Mesh_Topology* meshTopo, std::map<int,double> gbMap, MPI_Comm comm)
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
    
    std::map<int,double> gu_c_x_m;
    std::map<int,double> gu_c_y_m;
    std::map<int,double> gu_c_z_m;
    
    std::map<int,Array<double>*> gudxi;
    
    Array<double>* gu_c_x      = new Array<double>(nLoc_Elem,1);
    Array<double>* gu_c_y      = new Array<double>(nLoc_Elem,1);
    Array<double>* gu_c_z      = new Array<double>(nLoc_Elem,1);
    
    i_part_map* iee_vec = Pa->getIEEpartmap();
    Array<double>* gu_c_old    = new Array<double>(nLoc_Elem,3);

    std::map<int,std::vector<int> > gE2gF = Pa->getglobElem2globFaces();
    
    for(int i=0;i<nLoc_Elem;i++)
    {
        int gid  = Loc_Elem[i];
        
        gu_c_x_m[gid] = 1.0;
        gu_c_y_m[gid] = 1.0;
        gu_c_z_m[gid] = 1.0;
        
        gu_c_x->setVal(i,0,1.0);
        gu_c_y->setVal(i,0,1.0);
        gu_c_z->setVal(i,0,1.0);
        
        gu_c_old->setVal(i,0,1.0);
        gu_c_old->setVal(i,1,1.0);
        gu_c_old->setVal(i,2,1.0);
    }
    
    //std::cout << "Computing the MGG " << std::endl;
    std::map<int,vector<Vec3D*> > normals   = meshTopo->getNormals();
    std::map<int,vector<Vec3D*> > rvector   = meshTopo->getRvectors();
    std::map<int,vector<Vec3D*> > dxfxc     = meshTopo->getdXfXc();
    std::map<int,vector<double> > dS        = meshTopo->getdS();
    std::map<int,vector<double> > dr        = meshTopo->getdr();
    std::map<int,double > vol               = meshTopo->getVol();
    std::map<int,std::vector<Vert*> > vfvec = meshTopo->getVfacevector();
    std::map<int,Vert*> ghostvrts = meshTopo->getGhostVerts();

    std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();

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
    for(int it=0;it<1;it++)
    {
        t = clock();
        //communicate grad phi!!!
        
        Pa->AddStateForAdjacentElements(gu_c_x_m, comm);
        Pa->AddStateForAdjacentElements(gu_c_y_m, comm);
        Pa->AddStateForAdjacentElements(gu_c_z_m, comm);
        
//        std::map<int,double> dUdx_p_bnd = Pa->CommunicateStateAdjacentElements(gu_c_x_m, comm);
//        std::map<int,double> dUdy_p_bnd = Pa->CommunicateStateAdjacentElements(gu_c_y_m, comm);
//        std::map<int,double> dUdz_p_bnd = Pa->CommunicateStateAdjacentElements(gu_c_z_m, comm);

        //std::map<int,std::vector<double> > dUdxi_p_bnd = Pa->CommunicateStateAdjacentElementsNew(gu_c_old, comm);
                
        L2normx = 0.0;
        L2normy = 0.0;
        L2normz = 0.0;
        
        for(int i=0;i<nLoc_Elem;i++)
        {
             gEl        = Loc_Elem[i];
             lid        = i;
             int nadj   = LocElem2Nf[gEl];

             u_c = U[gEl]->getVal(0,0);
             
             gu_c_vx = gu_c_x->getVal(lid,0);
             gu_c_vy = gu_c_y->getVal(lid,0);
             gu_c_vz = gu_c_z->getVal(lid,0);

             sum_phix = 0.0;
             sum_phiy = 0.0;
             sum_phiz = 0.0;
            
             if(rvector[gEl].size()!=nadj || dxfxc[gEl].size()!=nadj)
             {
                 std::cout << "Huge error " <<  rvector[gEl].size() << " "  << nadj << " " << dxfxc[gEl].size() << std::endl;
             }
             for(int j=0;j<nadj;j++)
             {
                 adjID   = iee_vec->i_map[gEl][j];
                 
                 if(adjID<Nel)
                 {
                     //l_adjid = gE2lE[adjID];
                     gu_nb_vx = gu_c_x_m[adjID];
                     gu_nb_vy = gu_c_y_m[adjID];
                     gu_nb_vz = gu_c_z_m[adjID];
                     
//                   gu_nb_vx = dUdxi_p_bnd[adjID][0];//dUdx_p_bnd[adjID];
//                   gu_nb_vy = dUdxi_p_bnd[adjID][1];//dUdy_p_bnd[adjID];
//                   gu_nb_vz = dUdxi_p_bnd[adjID][2];//dUdz_p_bnd[adjID];
                     
                     u_nb = U[adjID]->getVal(0,0);
                 }
                 else
                 {
                     Vert* Vc = ghostvrts[adjID];
                     
                     double rtje = sqrt(Vc->x*Vc->x+(Vc->y-0.5)*(Vc->y-0.5)+(Vc->z-0.5)*(Vc->z-0.5));
                     double Utje = 0.1*tanh(50*(rtje-0.5))+1.0;
                                     //double Utje    = 0.1*sin(50*Vc->x*Vc->z)+atan(0.1/((sin(5.0*Vc->y)-2.0*Vc->x*Vc->z)));
                     
                     u_nb     = Utje;
                     gu_nb_vx = gu_c_x->getVal(lid,0);
                     gu_nb_vy = gu_c_y->getVal(lid,0);
                     gu_nb_vz = gu_c_z->getVal(lid,0);

                     
                 }
                 
                 nj          = normals[gEl][j];
                 rj          = rvector[gEl][j];
                 
                 double alpha = DotVec3D(nj,rj);

                 Vec3D* nf_m_arf = new Vec3D;

                 nf_m_arf->c0=nj->c0-alpha*rj->c0;
                 nf_m_arf->c1=nj->c1-alpha*rj->c1;
                 nf_m_arf->c2=nj->c2-alpha*rj->c2;

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
             
             gu_c_x_m[gEl] = 1.0/Vol*sum_phix;
             gu_c_y_m[gEl] = 1.0/Vol*sum_phiy;
             gu_c_z_m[gEl] = 1.0/Vol*sum_phiz;

             L2normx = L2normx+ sqrt((gu_c_x->getVal(i,0)-gu_c_old->getVal(i,0))*(gu_c_x->getVal(i,0)-gu_c_old->getVal(i,0)));

             L2normy = L2normy+ sqrt((gu_c_y->getVal(i,0)-gu_c_old->getVal(i,1))*(gu_c_y->getVal(i,0)-gu_c_old->getVal(i,1)));

             L2normz = L2normz+ sqrt((gu_c_z->getVal(i,0)-gu_c_old->getVal(i,2))*(gu_c_z->getVal(i,0)-gu_c_old->getVal(i,2)));
            
            gudxi[gEl] = new Array<double>(3,1);
            gudxi[gEl]->setVal(0,0,1.0/Vol*sum_phix);
            gudxi[gEl]->setVal(1,0,1.0/Vol*sum_phiy);
            gudxi[gEl]->setVal(2,0,1.0/Vol*sum_phiz);
         }
        
//        dUdx_p_bnd.clear();
//        dUdy_p_bnd.clear();
//        dUdz_p_bnd.clear();
        
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
        
        if(L2normx_max<1.0e-09 && L2normy_max<1.0e-09 && L2normz_max<1.0e-09)
        {
            break;
        }
        
    }
    
    delete gu_c_old;
    delete gu_c_x;
    delete gu_c_y;
    delete gu_c_z;
    
    return gudxi;
}

