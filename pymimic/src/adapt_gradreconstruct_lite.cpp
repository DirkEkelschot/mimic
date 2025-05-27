#include "adapt_gradreconstruct_lite.h"

std::map<int,std::vector<double> > ComputedUdx_LSQ_LS_US3D_Lite_V2(RepartitionObject* RePa,
                                                           std::map<int,std::vector<double> > ghosts,
                                                           int Nel,
                                                           int variable,
                                                           int nvariables,
                                                           MPI_Comm comm,
                                                           int extrap)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();

   std::map<int,std::vector<double> > Ue             = RePa->getElement2DataMap();
   RePa->AddStateVecForAdjacentElements(Ue,nvariables,comm);
    // if(world_rank == 0)
    // {
    //     std::cout << "Running ComputedUdx_LSQ_LS_US3D_Lite_V2 " << std::endl;
    // }

   std::map<int,std::vector<double> > LocalVs           = RePa->getLocalVertsMap();
   std::vector<int> Loc_Elem                            = RePa->getLocElem();
   int nLoc_Elem                                        = Loc_Elem.size();
   std::map<int, std::vector<int> > Face2VertexMap      = RePa->getFace2VertexMap();
   std::map<int, std::vector<int> > Element2FaceMap     = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2FacesMap    = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2VertexMap   = RePa->getElement2VertexMap();
   std::map<int, std::vector<int> > Element2ElementMap  = RePa->getElement2ElementMap();

//    std::cout << "before Usize " << Ue.size() << " " << Element2VertexMap.size() << " " << world_rank << std::endl;
//    RePa->AddStateVecForAdjacentElements(Ue,1,comm);
//    std::cout << "after Usize " << Ue.size() << " " << Element2VertexMap.size() << " " << world_rank << std::endl;

   std::map<int, std::vector<int> >::iterator itmii;


     

   std::map<int,std::vector<double> > dudx_map;
   std::map<int,std::vector<double> > dudx_map2update;
   double d;
   int loc_vid,adjID,elID;
   int cou = 0;

   std::vector<double> Vc(3);
   std::vector<double> Vadj(3);
   int lid = 0;
   double u_ijk, u_po;
    
   
   std::map<int,std::vector<double> > vrt_collect;
   std::map<int,double> sol_collect;
   std::vector<std::vector<double> > face;
   double rdotn;
    
   std::vector<double> n0(3);
   std::vector<double> v0(3);
   std::vector<double> v1(3);

   int rr =0;
   std::set<int> vert_scheme;
   int cc = 0;
   
   for(int i=0;i<nLoc_Elem;i++)
   {

        int ghostelem = 0;
        std::vector<double> x;
        int elID     = Loc_Elem[i];
        int NvPEl    = Element2VertexMap[elID].size();
        int nadj     = Element2ElementMap[elID].size();
        
        std::vector<double> Pijk(NvPEl*3);
        for(int k=0;k<Element2VertexMap[elID].size();k++)
        {
            int global_vid = Element2VertexMap[elID][k];
            Pijk[k*3+0] = LocalVs[global_vid][0];
            Pijk[k*3+1] = LocalVs[global_vid][1];
            Pijk[k*3+2] = LocalVs[global_vid][2];
        }
        
        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,NvPEl);
        u_ijk = Ue[elID][variable];

        for(int j=0;j<nadj;j++)
        {
            int adjid  = Element2ElementMap[elID][j];
            int Faceid = Element2FaceMap[elID][j];
                       
            if(Element2VertexMap.find(adjid)!=Element2VertexMap.end() 
                        && ghosts.find(adjid)==ghosts.end())
            {
                int NvPEladj = Element2VertexMap[adjid].size();

                std::vector<double> Padj(NvPEladj*3);
                for(int k=0;k<Element2VertexMap[adjid].size();k++)
                {
                    int vtag     = Element2VertexMap[adjid][k];
                    Padj[k*3+0] = LocalVs[vtag][0];
                    Padj[k*3+1] = LocalVs[vtag][1];
                    Padj[k*3+2] = LocalVs[vtag][2];
                }

                std::vector<double> Vadj = ComputeCentroidCoord(Padj,NvPEladj);
                            
                vrt_collect[adjid]  = Vadj;
                sol_collect[adjid]  = Ue[adjid][variable];
            }
            else if(ghosts.find(adjid)!=ghosts.end())
            {
                ghostelem++;
                if(vrt_collect.find(adjid)==vrt_collect.end())
                {
                    
                    int NvPerF = Face2VertexMap[Faceid].size();
                    
                    std::vector<double> Vc(3);
                    Vc[0] = 0.0;
                    Vc[1] = 0.0;
                    Vc[2] = 0.0;
                    
                    for(int s=0;s<NvPerF;s++)
                    {
                        int gvid    = Face2VertexMap[Faceid][s];
                        
                        Vc[0]       = Vc[0]+LocalVs[gvid][0];
                        Vc[1]       = Vc[1]+LocalVs[gvid][1];
                        Vc[2]       = Vc[2]+LocalVs[gvid][2];
                        
                        std::vector<double> V(3,0.0);
                        V[0]        = LocalVs[gvid][0];
                        V[1]        = LocalVs[gvid][1];
                        V[2]        = LocalVs[gvid][2];

                        face.push_back(V);    
                    }
                    
                    Vc[0] = Vc[0]/NvPerF;
                    Vc[1] = Vc[1]/NvPerF;
                    Vc[2] = Vc[2]/NvPerF;
                
                    std::vector<double> r0(3);
                    r0[0] = (Vc[0]-Vijk[0]);
                    r0[1] = (Vc[1]-Vijk[1]);
                    r0[2] = (Vc[2]-Vijk[2]);
                    

                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];
                    
                    v1[0] = face[2][0]-face[0][0];
                    v1[1] = face[2][1]-face[0][1];
                    v1[2] = face[2][2]-face[0][2];
                   
                    std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);
                    double orient0   = DotVec3D(r0,n0);
                    
                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }
                    
                    rdotn = DotVec3D(r0,n0);
                    
                    std::vector<double> reflect(3);
                    
                    reflect[0] = r0[0]-2.0*(rdotn)*n0[0];
                    reflect[1] = r0[1]-2.0*(rdotn)*n0[1];
                    reflect[2] = r0[2]-2.0*(rdotn)*n0[2];
                    
                    Vc[0] = Vc[0] - reflect[0];
                    Vc[1] = Vc[1] - reflect[1];
                    Vc[2] = Vc[2] - reflect[2];
                    
                    vrt_collect[adjid]  = Vc;
                    sol_collect[adjid]  = ghosts[adjid][variable];
                
                    face.clear();
                }
            }
        }

        int Ndata = vrt_collect.size();
        //std::cout << "vrt_collect " << Ndata << std::endl;
        std::vector<std::vector<double> > Vrt(Ndata);
        std::vector<double> bvec(Ndata,0);
        
        for(int q=0;q<Ndata;q++)
        {
            std::vector<double> row(3,0);
            Vrt[q] = row;
        }
        std::map<int,std::vector<double> >::iterator vit;
        int te = 0;

        double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;
        
        for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
        {
            double di = sqrt((vit->second[0]-Vijk[0])*(vit->second[0]-Vijk[0])+
                                (vit->second[1]-Vijk[1])*(vit->second[1]-Vijk[1])+
                                (vit->second[2]-Vijk[2])*(vit->second[2]-Vijk[2]));

            a = (vit->second[0] - Vijk[0]);
            b = (vit->second[1] - Vijk[1]);
            c = (vit->second[2] - Vijk[2]);

            h00 = 0.5*a*a; h01 = 1.0*a*b; h02 = 1.0*a*c;
            h11 = 0.5*b*b; h12 = 1.0*b*c;
            h22 = 0.5*c*c;

            Vrt[te][0] = (1.0/di)*a;
            Vrt[te][1] = (1.0/di)*b;
            Vrt[te][2] = (1.0/di)*c;

            // Vrt[te][3] = (1.0/di)*h00;
            // Vrt[te][4] = (1.0/di)*h01;
            // Vrt[te][5] = (1.0/di)*h02;
            // Vrt[te][6] = (1.0/di)*h11;
            // Vrt[te][7] = (1.0/di)*h12;
            // Vrt[te][8] = (1.0/di)*h22;

            double Udata = sol_collect[vit->first];



            bvec[te]   = (1.0/di)*(Udata-u_ijk);

            // std::cout << "te " << te << " " << bvec[te] << std::endl;
            te++;
        }

        std::vector<double> A_cm(Ndata*3,0.0);
        for(int s=0;s<Ndata;s++)
        {
            for(int g=0;g<3;g++)
            {
                A_cm[g*Ndata+s] = Vrt[s][g];
            }
        }

        x = SolveQR_Lite(A_cm,Ndata,3,bvec);

        dudx_map[elID] = x;

        vrt_collect.clear();
        sol_collect.clear();
   }
    
   return dudx_map;
}










std::map<int,std::vector<double> > ComputedUdx_LSQ_LS_US3D_Lite(RepartitionObject* RePa,
                                                           std::map<int,std::vector<double> > ghosts,
                                                           int Nel,
                                                           int variable,
                                                           int nvariables,
                                                           MPI_Comm comm,
                                                           int extrap)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();

   std::map<int,std::vector<double> > Ue             = RePa->getElement2DataMap();
    //    RePa->AddStateVecForAdjacentElements(Ue,nvariables,comm);


   std::map<int,std::vector<double> > LocalVs           = RePa->getLocalVertsMap();
   std::vector<int> Loc_Elem                            = RePa->getLocElem();
   int nLoc_Elem                                        = Loc_Elem.size();
   std::map<int, std::vector<int> > Face2VertexMap      = RePa->getFace2VertexMap();
   std::map<int, std::vector<int> > Element2FaceMap     = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2FacesMap    = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2VertexMap   = RePa->getElement2VertexMap();
   std::map<int, std::vector<int> > Element2ElementMap  = RePa->getElement2ElementMap();

//    std::cout << "before Usize " << Ue.size() << " " << Element2VertexMap.size() << " " << world_rank << std::endl;
//    RePa->AddStateVecForAdjacentElements(Ue,1,comm);
//    std::cout << "after Usize " << Ue.size() << " " << Element2VertexMap.size() << " " << world_rank << std::endl;

   std::map<int, std::vector<int> >::iterator itmii;


     

   std::map<int,std::vector<double> > dudx_map;
   std::map<int,std::vector<double> > dudx_map2update;
   double d;
   int loc_vid,adjID,elID;
   int cou = 0;

   std::vector<double> Vc(3);
   std::vector<double> Vadj(3);
   int lid = 0;
   double u_ijk, u_po;
    
   std::map<int,std::vector<double> > vrt_collect;
   std::map<int,double> sol_collect;
   std::vector<std::vector<double> > face;
   double rdotn;
    
   std::vector<double> n0(3);
   std::vector<double> v0(3);
   std::vector<double> v1(3);

   int rr =0;
   std::set<int> vert_scheme;
   int cc = 0;
   
   for(int i=0;i<nLoc_Elem;i++)
   {
        int ghostelem = 0;
        std::vector<double> x;
        int elID     = Loc_Elem[i];
        int NvPEl    = Element2VertexMap[elID].size();
        int nadj     = Element2ElementMap[elID].size();
        
        std::vector<double> Pijk(NvPEl*3);
        for(int k=0;k<Element2VertexMap[elID].size();k++)
        {
            int global_vid = Element2VertexMap[elID][k];
            Pijk[k*3+0] = LocalVs[global_vid][0];
            Pijk[k*3+1] = LocalVs[global_vid][1];
            Pijk[k*3+2] = LocalVs[global_vid][2];
        }
        
        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,NvPEl);
        u_ijk = Ue[elID][variable];
        //std::cout << "u_ijk " << u_ijk << " " << variable << std::endl;
        for(int j=0;j<nadj;j++)
        {
            int adjid  = Element2ElementMap[elID][j];
            int Faceid = Element2FaceMap[elID][j];
                       
            if(Element2VertexMap.find(adjid)!=Element2VertexMap.end() 
                        && ghosts.find(adjid)==ghosts.end())
            {
                int NvPEladj = Element2VertexMap[adjid].size();

                std::vector<double> Padj(NvPEladj*3);
                for(int k=0;k<Element2VertexMap[adjid].size();k++)
                {
                    int vtag    = Element2VertexMap[adjid][k];
                    Padj[k*3+0] = LocalVs[vtag][0];
                    Padj[k*3+1] = LocalVs[vtag][1];
                    Padj[k*3+2] = LocalVs[vtag][2];
                }

                std::vector<double> Vadj = ComputeCentroidCoord(Padj,NvPEladj);
                            
                vrt_collect[adjid]  = Vadj;
                sol_collect[adjid]  = Ue[adjid][variable];
            }
            else if(ghosts.find(adjid)!=ghosts.end())
            {
                ghostelem++;
                if(vrt_collect.find(adjid)==vrt_collect.end())
                {
                    
                    int NvPerF = Face2VertexMap[Faceid].size();
                    
                    std::vector<double> Vc(3);
                    Vc[0] = 0.0;
                    Vc[1] = 0.0;
                    Vc[2] = 0.0;
                    
                    for(int s=0;s<NvPerF;s++)
                    {
                        int gvid    = Face2VertexMap[Faceid][s];
                        
                        Vc[0]       = Vc[0]+LocalVs[gvid][0];
                        Vc[1]       = Vc[1]+LocalVs[gvid][1];
                        Vc[2]       = Vc[2]+LocalVs[gvid][2];
                        
                        std::vector<double> V(3,0.0);
                        V[0]        = LocalVs[gvid][0];
                        V[1]        = LocalVs[gvid][1];
                        V[2]        = LocalVs[gvid][2];

                        face.push_back(V);    
                    }
                    
                    Vc[0] = Vc[0]/NvPerF;
                    Vc[1] = Vc[1]/NvPerF;
                    Vc[2] = Vc[2]/NvPerF;
                
                    std::vector<double> r0(3);
                    r0[0] = (Vc[0]-Vijk[0]);
                    r0[1] = (Vc[1]-Vijk[1]);
                    r0[2] = (Vc[2]-Vijk[2]);
                    

                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];
                    
                    v1[0] = face[2][0]-face[0][0];
                    v1[1] = face[2][1]-face[0][1];
                    v1[2] = face[2][2]-face[0][2];
                   
                    std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);
                    double orient0   = DotVec3D(r0,n0);
                    
                    if(orient0<0.0)
                    {
                        NegateVec3D(n0);
                    }
                    
                    rdotn = DotVec3D(r0,n0);
                    
                    std::vector<double> reflect(3);
                    
                    reflect[0] = r0[0]-2.0*(rdotn)*n0[0];
                    reflect[1] = r0[1]-2.0*(rdotn)*n0[1];
                    reflect[2] = r0[2]-2.0*(rdotn)*n0[2];
                    
                    Vc[0] = Vc[0] - reflect[0];
                    Vc[1] = Vc[1] - reflect[1];
                    Vc[2] = Vc[2] - reflect[2];
                    
                    vrt_collect[adjid]  = Vc;

                    if(extrap==1)
                    {
                        sol_collect[adjid]  = u_ijk;
                    }
                    else
                    {
                        sol_collect[adjid]  = ghosts[adjid][variable];
                    }
                    
                
                    face.clear();
                }
            }


            
            if(Element2ElementMap.find(adjid)!=Element2ElementMap.end())
            {
                int n_adjid         = Element2ElementMap[adjid].size();
                int NvPElnew        = Element2VertexMap[adjid].size();

                std::vector<double> Pijknew(NvPElnew*3);
                for(int k=0;k<NvPElnew;k++)
                {
                    int glob_vid    = Element2VertexMap[adjid][k];
                    Pijknew[k*3+0]  = LocalVs[glob_vid][0];
                    Pijknew[k*3+1]  = LocalVs[glob_vid][1];
                    Pijknew[k*3+2]  = LocalVs[glob_vid][2];
                }

                std::vector<double> Vijknew = ComputeCentroidCoord(Pijknew,NvPElnew);
                
                for(int k=0;k<n_adjid;k++)
                {
                    int adjadj      = Element2ElementMap[adjid][k];
                    int fadjadj     = Element2FacesMap[adjid][k];

                    if(Element2VertexMap.find(adjadj)!=Element2VertexMap.end() 
                        && ghosts.find(adjadj)==ghosts.end() && adjadj!=elID)
                    {
                        int NvPEladjadj    = Element2VertexMap[adjadj].size();
                        std::vector<double> Padjadj(NvPEladjadj*3);
                        for(int k=0;k<Element2VertexMap[adjadj].size();k++)
                        {
                            int vtag           = Element2VertexMap[adjadj][k];
                            Padjadj[k*3+0]     = LocalVs[vtag][0];
                            Padjadj[k*3+1]     = LocalVs[vtag][1];
                            Padjadj[k*3+2]     = LocalVs[vtag][2];
                        }
                        
                        std::vector<double> Vadjadj = ComputeCentroidCoord(Padjadj,NvPEladjadj);

                        vrt_collect[adjadj] = Vadjadj;
                        sol_collect[adjadj] = Ue[adjadj][variable];
                        
                    }
                    else if(ghosts.find(adjadj)!=ghosts.end())
                    {

                        ghostelem++;
                        if(vrt_collect.find(adjadj)==vrt_collect.end())
                        {
                            int NvPerF = Face2VertexMap[fadjadj].size();
                            
                            std::vector<double> Vc(3);
                            Vc[0]       = 0.0;
                            Vc[1]       = 0.0;
                            Vc[2]       = 0.0;
                            
                            for(int s=0;s<NvPerF;s++)
                            {
                                int gvid = Face2VertexMap[fadjadj][s];

                                Vc[0] = Vc[0]+LocalVs[gvid][0];
                                Vc[1] = Vc[1]+LocalVs[gvid][1];
                                Vc[2] = Vc[2]+LocalVs[gvid][2];
                                
                                std::vector<double> V(3);
                                V[0]    = LocalVs[gvid][0];
                                V[1]    = LocalVs[gvid][1];
                                V[2]    = LocalVs[gvid][2];
                                face.push_back(V);
                            }

                            Vc[0] = Vc[0]/NvPerF;
                            Vc[1] = Vc[1]/NvPerF;
                            Vc[2] = Vc[2]/NvPerF;
                            
                            
                            std::vector<double> r0(3);
                            r0[0] = (Vc[0]-Vijknew[0]);
                            r0[1] = (Vc[1]-Vijknew[1]);
                            r0[2] = (Vc[2]-Vijknew[2]);
                            
                            v0[0] = face[1][0]-face[0][0];
                            v0[1] = face[1][1]-face[0][1];
                            v0[2] = face[1][2]-face[0][2];

                            v1[0] = face[2][0]-face[0][0];
                            v1[1] = face[2][1]-face[0][1];
                            v1[2] = face[2][2]-face[0][2];
                            
                            std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);
                            double orient0   = DotVec3D(r0,n0);
                            
                            if(orient0<0.0)
                            {
                                NegateVec3D(n0);
                            }
                            
                            rdotn = DotVec3D(r0,n0);
                            
                            std::vector<double> reflect(3);
                            
                            reflect[0] = r0[0]-2.0*(rdotn)*n0[0];
                            reflect[1] = r0[1]-2.0*(rdotn)*n0[1];
                            reflect[2] = r0[2]-2.0*(rdotn)*n0[2];
                            
                            Vc[0] = Vc[0] - reflect[0];
                            Vc[1] = Vc[1] - reflect[1];
                            Vc[2] = Vc[2] - reflect[2];
                                            

                            vrt_collect[adjadj]  = Vc;

                            if(extrap==1)
                            {
                                sol_collect[adjid]  = u_ijk;
                            }
                            else
                            {
                                sol_collect[adjadj]  = ghosts[adjadj][variable];
                            }
                            
                            
                            
                            face.clear();
                        }                       
                    }
                    if(Element2ElementMap.find(adjadj)!=Element2ElementMap.end())
                    {
                        int n_adjadj        = Element2ElementMap[adjadj].size();
                        int NvPElnewnew     = Element2VertexMap[adjadj].size();
                        std::vector<double> Pijknewnew(NvPElnewnew*3);
                        for(int k=0;k<NvPElnewnew;k++)
                        {
                            int vtag     = Element2VertexMap[adjadj][k];
                            Pijknewnew[k*3+0] = LocalVs[vtag][0];
                            Pijknewnew[k*3+1] = LocalVs[vtag][1];
                            Pijknewnew[k*3+2] = LocalVs[vtag][2];
                            
                        }
                        
                        std::vector<double> Vijknewnew = ComputeCentroidCoord(Pijknewnew,NvPElnewnew);
                        
                        for(int k=0;k<n_adjadj;k++)
                        {
                            int adjadjadj  = Element2ElementMap[adjadj][k];
                            int fadjadjadj = Element2FaceMap[adjadj][k];

                            if(Element2VertexMap.find(adjadjadj)!=Element2VertexMap.end() 
                                && ghosts.find(adjadjadj)==ghosts.end() && adjadjadj!=elID)
                            {

                                int NvPEladjadjadj    = Element2VertexMap[adjadjadj].size();
                                std::vector<double> Padjadjadj(NvPEladjadjadj*3);
                                for(int k=0;k<Element2VertexMap[adjadjadj].size();k++)
                                {
                                    int vtag              = Element2VertexMap[adjadjadj][k];
                                    Padjadjadj[k*3+0]     = LocalVs[vtag][0];
                                    Padjadjadj[k*3+1]     = LocalVs[vtag][1];
                                    Padjadjadj[k*3+2]     = LocalVs[vtag][2];
                                }
                                
                                std::vector<double> Vadjadjadj = ComputeCentroidCoord(Padjadjadj,NvPEladjadjadj);

                                vrt_collect[adjadjadj] = Vadjadjadj;
                                sol_collect[adjadjadj] = Ue[adjadjadj][variable];

                            }
                            else if(ghosts.find(adjadjadj)!=ghosts.end())
                            {
                                ghostelem++;
                                int NvPerF = Face2VertexMap[fadjadjadj].size();
                                
                                std::vector<double> Vc(3);
                                Vc[0] = 0.0;
                                Vc[1] = 0.0;
                                Vc[2] = 0.0;
                                
                                for(int s=0;s<NvPerF;s++)
                                {
                                    int gvid = Face2VertexMap[fadjadjadj][s];

                                    Vc[0] = Vc[0]+LocalVs[gvid][0];
                                    Vc[1] = Vc[1]+LocalVs[gvid][1];
                                    Vc[2] = Vc[2]+LocalVs[gvid][2];
                                    
                                    std::vector<double> V(3);
                                    V[0]    = LocalVs[gvid][0];
                                    V[1]    = LocalVs[gvid][1];
                                    V[2]    = LocalVs[gvid][2];
                                    face.push_back(V);
                                }

                                Vc[0] = Vc[0]/NvPerF;
                                Vc[1] = Vc[1]/NvPerF;
                                Vc[2] = Vc[2]/NvPerF;
                                
                                std::vector<double> r0(3);
                                r0[0] = (Vc[0]-Vijknewnew[0]);
                                r0[1] = (Vc[1]-Vijknewnew[1]);
                                r0[2] = (Vc[2]-Vijknewnew[2]);
                                
                                v0[0] = face[1][0]-face[0][0];
                                v0[1] = face[1][1]-face[0][1];
                                v0[2] = face[1][2]-face[0][2];

                                v1[0] = face[2][0]-face[0][0];
                                v1[1] = face[2][1]-face[0][1];
                                v1[2] = face[2][2]-face[0][2];
                                
                                std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);
                                double orient0   = DotVec3D(r0,n0);
                                
                                if(orient0<0.0)
                                {
                                    NegateVec3D(n0);
                                }
                                
                                rdotn = DotVec3D(r0,n0);
                                
                                std::vector<double> reflect(3);
                                
                                reflect[0] = r0[0]-2.0*(rdotn)*n0[0];
                                reflect[1] = r0[1]-2.0*(rdotn)*n0[1];
                                reflect[2] = r0[2]-2.0*(rdotn)*n0[2];
                                
                                Vc[0] = Vc[0] - reflect[0];
                                Vc[1] = Vc[1] - reflect[1];
                                Vc[2] = Vc[2] - reflect[2];
    
                                vrt_collect[adjadjadj]  = Vc;
                                if(extrap==1)
                                {
                                    sol_collect[adjid]  = u_ijk;
                                }
                                else
                                {
                                    sol_collect[adjadjadj]  = ghosts[adjadjadj][variable];
                                }
                                
                                face.clear();
                                
                                // cntbnd++;
                                
                            }
                        }
                    }
                }
            }
        }

        //std::cout << "vrt_collect " << vrt_collect.size() << std::endl;
        if(vrt_collect.size() >= 9)
        {
            int Ndata = vrt_collect.size();

            std::vector<std::vector<double> > Vrt(Ndata);
            std::vector<double> bvec(Ndata,0);
           
            for(int q=0;q<Ndata;q++)
            {
                std::vector<double> row(9,0);
                Vrt[q] = row;
            }
            std::map<int,std::vector<double> >::iterator vit;
            int te = 0;

            double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;
            
            for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
            {
                double di = sqrt((vit->second[0]-Vijk[0])*(vit->second[0]-Vijk[0])+
                                 (vit->second[1]-Vijk[1])*(vit->second[1]-Vijk[1])+
                                 (vit->second[2]-Vijk[2])*(vit->second[2]-Vijk[2]));

                a = (vit->second[0] - Vijk[0]);
                b = (vit->second[1] - Vijk[1]);
                c = (vit->second[2] - Vijk[2]);

                h00 = 0.5*a*a; h01 = 1.0*a*b; h02 = 1.0*a*c;
                h11 = 0.5*b*b; h12 = 1.0*b*c;
                h22 = 0.5*c*c;

                Vrt[te][0] = (1.0/di)*a;
                Vrt[te][1] = (1.0/di)*b;
                Vrt[te][2] = (1.0/di)*c;

                Vrt[te][3] = (1.0/di)*h00;
                Vrt[te][4] = (1.0/di)*h01;
                Vrt[te][5] = (1.0/di)*h02;
                Vrt[te][6] = (1.0/di)*h11;
                Vrt[te][7] = (1.0/di)*h12;
                Vrt[te][8] = (1.0/di)*h22;

                double Udata = sol_collect[vit->first];



                bvec[te]   = (1.0/di)*(Udata-u_ijk);

                // std::cout << "te " << te << " " << bvec[te] << std::endl;
                te++;
            }

            std::vector<double> A_cm(Ndata*9,0.0);
            for(int s=0;s<Ndata;s++)
            {
                for(int g=0;g<9;g++)
                {
                    A_cm[g*Ndata+s] = Vrt[s][g];
                }
            }

            x = SolveQR_Lite(A_cm,Ndata,9,bvec);

            dudx_map[elID] = x;
            
            
        }
        else
        {
            int Ndata = vrt_collect.size();
            
            std::cout << "Warning:: not enough data points to reconstruct the gradient! Number of neigboring points is " << Ndata << " number of ghost elements -> " << ghostelem << " for element ID " << elID << std::endl;
            
            std::vector<std::vector<double> > Vrt(Ndata);
            std::vector<double> bvec(Ndata,0);
           
            for(int q=0;q<Ndata;q++)
            {
                std::vector<double> row(3,0);
                Vrt[q] = row;
            }
            std::map<int,std::vector<double> >::iterator vit;
            int te = 0;
            
            double a,b,c,h00,h01,h02,h10,h11,h12,h20,h21,h22;
            
            for(vit=vrt_collect.begin();vit!=vrt_collect.end();vit++)
            {
                double di = sqrt((vit->second[0]-Vijk[0])*(vit->second[0]-Vijk[0])+
                                 (vit->second[1]-Vijk[1])*(vit->second[1]-Vijk[1])+
                                 (vit->second[2]-Vijk[2])*(vit->second[2]-Vijk[2]));

                a = (vit->second[0] - Vijk[0]);
                b = (vit->second[1] - Vijk[1]);
                c = (vit->second[2] - Vijk[2]);
                
                Vrt[te][0]      = (1.0/di)*a;
                Vrt[te][1]      = (1.0/di)*b;
                Vrt[te][2]      = (1.0/di)*c;
                double Udata    = sol_collect[vit->first];
                bvec[te]        = (1.0/di)*(Udata-u_ijk);
                te++;
            }

            std::vector<double> A_cm(Ndata*3,0);
            for(int s=0;s<Ndata;s++)
            {
                for(int g=0;g<3;g++)
                {
                    A_cm[g*Ndata+s] = Vrt[s][g];
                }
            }

            x = SolveQR_Lite(A_cm,Ndata,3,bvec);

            std::vector<double> xnew(9,0.0);
            xnew[0] = x[0];
            xnew[1] = x[1];
            xnew[2] = x[2];
            
            // x.resize(9);
            // x[0] = 0.0;
            // x[1] = 0.0;
            // x[2] = 0.0;
            // x[3] = 0.0;
            // x[4] = 0.0;
            // x[5] = 0.0;
            // x[6] = 0.0;
            // x[7] = 0.0;
            // x[8] = 0.0;
            dudx_map[elID] = xnew;


            // dudx_map2update[elID] = xnew;
            // cc++;
            // delete[] A_cm;
            // ////delete[] Pijk;
            // delete Vrt_T;
            // delete Vrt;
            // delete bvec;
        }

        
        //
        
        vrt_collect.clear();
        sol_collect.clear();
        
   }
   std::map<int,std::vector<double> >::iterator itr;

    if(dudx_map2update.size()!=0)
    {
        std::cout << "We are not updating " << std::endl;
        int inhere = 0;
        int inhere2 = 0;
        std::map<int, std::vector<int> > new_map;
        for(itr=dudx_map2update.begin();itr!=dudx_map2update.end();itr++)
        {
            int elid = itr->first;
            
            if(Element2ElementMap.find(elid)!=Element2ElementMap.end())
            {
                inhere++;
                int nadj = Element2ElementMap[elid].size();
                for(int i=0;i<nadj;i++)
                {
                    int adjid = Element2ElementMap[elid][i];

                    if(dudx_map.find(adjid)!=dudx_map.end())
                    {
                        new_map[elid].push_back(adjid);
                    }

                    if(Element2ElementMap.find(adjid)!=Element2ElementMap.end())
                    {
                        inhere2++;
                        int nadj2 = Element2ElementMap[adjid].size();
                        for(int j=0;j<nadj2;j++)
                        {
                            int adjid2 = Element2ElementMap[adjid][j];
                            if(dudx_map.find(adjid2)!=dudx_map.end())
                            {
                                new_map[elid].push_back(adjid2);
                            }

                            if(Element2ElementMap.find(adjid2)!=Element2ElementMap.end())
                            {
                                int nadj3 = Element2ElementMap[adjid2].size();
                                for(int k=0;k<nadj3;k++)
                                {
                                    int adjid3 = Element2ElementMap[adjid2][k];
                                    if(dudx_map.find(adjid3)!=dudx_map.end())
                                    {
                                        new_map[elid].push_back(adjid3);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        std::cout << "new_map.size() " << new_map.size() << " inhere " << inhere << " " << inhere2 << std::endl;
        std::map<int,std::vector<int> >::iterator itr0;
        for(itr0=new_map.begin();itr0!=new_map.end();itr0++)
        {
                std::cout << " itr0 " << itr0->first << " " << itr0->second.size() << std::endl;
        }
    }
    




   



   for(itr=dudx_map2update.begin();itr!=dudx_map2update.end();itr++)
   {
        int elid = itr->first;
        
        int nadj = Element2ElementMap[elid].size();
        
        std::vector<double> xnew(9,0);
        int tel = 0;
        for(int i=0;i<nadj;i++)
        {
            int adjid = Element2ElementMap[elid][i];
            
            if(dudx_map.find(adjid)!=dudx_map.end())
            {
                xnew[0] = xnew[0] + dudx_map[adjid][0];
                xnew[1] = xnew[1] + dudx_map[adjid][1];
                xnew[2] = xnew[2] + dudx_map[adjid][2];
                xnew[3] = xnew[3] + dudx_map[adjid][3];
                xnew[4] = xnew[4] + dudx_map[adjid][4];
                xnew[5] = xnew[5] + dudx_map[adjid][5];
                xnew[6] = xnew[6] + dudx_map[adjid][6];
                xnew[7] = xnew[7] + dudx_map[adjid][7];
                xnew[8] = xnew[8] + dudx_map[adjid][8];
                tel = tel + 1;
            }

            if(Element2ElementMap.find(adjid)!=Element2ElementMap.end()
                && ghosts.find(adjid)==ghosts.end())
            {
                int nadjadj = Element2ElementMap[adjid].size();

                for(int j=0;j<nadjadj;j++)
                {
                    int adjadjid = Element2ElementMap[adjid][j];

                    if(dudx_map.find(adjadjid)!=dudx_map.end())
                    {
                        xnew[0] = xnew[0] + dudx_map[adjadjid][0];
                        xnew[1] = xnew[1] + dudx_map[adjadjid][1];
                        xnew[2] = xnew[2] + dudx_map[adjadjid][2];
                        xnew[3] = xnew[3] + dudx_map[adjadjid][3];
                        xnew[4] = xnew[4] + dudx_map[adjadjid][4];
                        xnew[5] = xnew[5] + dudx_map[adjadjid][5];
                        xnew[6] = xnew[6] + dudx_map[adjadjid][6];
                        xnew[7] = xnew[7] + dudx_map[adjadjid][7];
                        xnew[8] = xnew[8] + dudx_map[adjadjid][8];
                        tel = tel + 1;
                    }

                     if(Element2ElementMap.find(adjadjid)!=Element2ElementMap.end()
                     && ghosts.find(adjadjid)==ghosts.end())
                     {
                        int nadjadjadj = Element2ElementMap[adjadjid].size();
                        for(int k=0;k<nadjadjadj;k++)
                        {
                            int adjadjadjid = Element2ElementMap[adjadjid][k];
                            if(dudx_map.find(adjadjadjid)!=dudx_map.end())
                            {
                                xnew[0] = xnew[0] + dudx_map[adjadjadjid][0];
                                xnew[1] = xnew[1] + dudx_map[adjadjadjid][1];
                                xnew[2] = xnew[2] + dudx_map[adjadjadjid][2];
                                xnew[3] = xnew[3] + dudx_map[adjadjadjid][3];
                                xnew[4] = xnew[4] + dudx_map[adjadjadjid][4];
                                xnew[5] = xnew[5] + dudx_map[adjadjadjid][5];
                                xnew[6] = xnew[6] + dudx_map[adjadjadjid][6];
                                xnew[7] = xnew[7] + dudx_map[adjadjadjid][7];
                                xnew[8] = xnew[8] + dudx_map[adjadjadjid][8];
                                tel = tel + 1;
                            }

                            if(Element2ElementMap.find(adjadjadjid)!=Element2ElementMap.end()
                            && ghosts.find(adjadjadjid)==ghosts.end())
                            {
                                int nadjadjadjadj = Element2ElementMap[adjadjadjid].size();
                                for(int m=0;m<nadjadjadjadj;m++)
                                {
                                    int adjadjadjadjid = Element2ElementMap[adjadjadjid][m];
                                    if(dudx_map.find(adjadjadjadjid)!=dudx_map.end())
                                    {
                                        xnew[0] = xnew[0] + dudx_map[adjadjadjadjid][0];
                                        xnew[1] = xnew[1] + dudx_map[adjadjadjadjid][1];
                                        xnew[2] = xnew[2] + dudx_map[adjadjadjadjid][2];
                                        xnew[3] = xnew[3] + dudx_map[adjadjadjadjid][3];
                                        xnew[4] = xnew[4] + dudx_map[adjadjadjadjid][4];
                                        xnew[5] = xnew[5] + dudx_map[adjadjadjadjid][5];
                                        xnew[6] = xnew[6] + dudx_map[adjadjadjadjid][6];
                                        xnew[7] = xnew[7] + dudx_map[adjadjadjadjid][7];
                                        xnew[8] = xnew[8] + dudx_map[adjadjadjadjid][8];
                                        tel = tel + 1;
                                    }

                                    if(Element2ElementMap.find(adjadjadjadjid)!=Element2ElementMap.end()
                                        && ghosts.find(adjadjadjadjid)==ghosts.end())
                                    {
                                        int nadjadjadjadjid = Element2ElementMap[adjadjadjadjid].size();

                                        for(int n=0;n<nadjadjadjadjid;n++)
                                        {
                                            int adjadjadjadjadjid = Element2ElementMap[adjadjadjadjid][n];

                                            if(dudx_map.find(adjadjadjadjadjid)!=dudx_map.end())
                                            {
                                                xnew[0] = xnew[0] + dudx_map[adjadjadjadjadjid][0];
                                                xnew[1] = xnew[1] + dudx_map[adjadjadjadjadjid][1];
                                                xnew[2] = xnew[2] + dudx_map[adjadjadjadjadjid][2];
                                                xnew[3] = xnew[3] + dudx_map[adjadjadjadjadjid][3];
                                                xnew[4] = xnew[4] + dudx_map[adjadjadjadjadjid][4];
                                                xnew[5] = xnew[5] + dudx_map[adjadjadjadjadjid][5];
                                                xnew[6] = xnew[6] + dudx_map[adjadjadjadjadjid][6];
                                                xnew[7] = xnew[7] + dudx_map[adjadjadjadjadjid][7];
                                                xnew[8] = xnew[8] + dudx_map[adjadjadjadjadjid][8];
                                                tel = tel + 1;
                                            }
                                        }        
                                    }
                                }
                            }
                        }
                     }
                }
            }
        }

        if(tel != 0)
        {
            xnew[0] = xnew[0]/tel;
            xnew[1] = xnew[1]/tel;
            xnew[2] = xnew[2]/tel;
            xnew[3] = xnew[3]/tel;
            xnew[4] = xnew[4]/tel;
            xnew[5] = xnew[5]/tel;
            xnew[6] = xnew[6]/tel;
            xnew[7] = xnew[7]/tel;
            xnew[8] = xnew[8]/tel;
        }
        else
        {
            //std::cout << "Element ID " << elid <<  " has no proper gradient information." << std::endl;
        }
        
        // std::cout << "number of points found: " << tel << std::endl;
        // std::cout << "===================dUdx       dUdy        dUdz===================" << std::endl;
        // std::cout << "                    " << xnew[0] << "       " << xnew[1] << "       " << xnew[2] << std::endl;
        // std::cout << "===================         Hessian           ===================" << std::endl;
        // std::cout << "                    " << xnew[3] << "       " << xnew[4] << "       " << xnew[5] << std::endl;
        // std::cout << "                    " << xnew[4] << "       " << xnew[6] << "       " << xnew[7] << std::endl;
        // std::cout << "                    " << xnew[5] << "       " << xnew[7] << "       " << xnew[8] << std::endl;
        dudx_map[elid] = xnew;
   }
   
 
    
   return dudx_map;
}











std::map<int,std::vector<double> > ComputedUdx_LSQ_LS_US3D_Lite_Update_Test(RepartitionObject* RePa,
                                                           std::map<int,std::vector<double> > ghosts,
                                                           int Nel,
                                                           int variable,
                                                           int nvariables,
                                                           MPI_Comm comm,
                                                           int extrap)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();

   std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();
    //    RePa->AddStateVecForAdjacentElements(Ue,nvariables,comm);


   std::map<int,std::vector<double> > LocalVs           = RePa->getLocalVertsMap();
   std::vector<int> Loc_Elem                            = RePa->getLocElem();
   int nLoc_Elem                                        = Loc_Elem.size();
   std::map<int, std::vector<int> > Face2VertexMap      = RePa->getFace2VertexMap();
   std::map<int, std::vector<int> > Element2FaceMap     = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2FacesMap    = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2VertexMap   = RePa->getElement2VertexMap();
   std::map<int, std::vector<int> > Element2ElementMap  = RePa->getElement2ElementMap();

   std::map<int,std::vector<double> > ghostface_vrt     = RePa->getGhostFaceData();
   std::map<int, std::set<int> > element2adj            = RePa->getExtendedAdjacencyData();
//    std::map<int, std::set<int> > element2adj         = RePa->getExtendedAdjacencyData(ghostface_vrt);
   
   std::map<int, std::vector<double> > element2centroid = RePa->GetElement2CentroidData();
   std::map<int, std::vector<int> >::iterator itmii;
   std::map<int,std::vector<double> > dudx_map;
   double d;
   int loc_vid,adjID,elID;
   int cou = 0;

   std::vector<double> Vc(3);
   std::vector<double> Vadj(3);
   int lid = 0;
   double u_ijk, u_po;
    
   std::map<int,std::vector<double> > vrt_collect;
   std::map<int,double> sol_collect;
   std::vector<std::vector<double> > face;
   double rdotn;
    
   std::vector<double> n0(3);
   std::vector<double> v0(3);
   std::vector<double> v1(3);

   int rr = 0;
   std::set<int> vert_scheme;
   int cc = 0;

   for(int i=0;i<nLoc_Elem;i++)
   {
        std::vector<double> x(9,0.0);
        int ghostelem               = 0;
        int elID                    = Loc_Elem[i];
        int nadj                    = element2adj[elID].size();
        //std::cout << " nadj " <<  nadj << std::endl;
        std::set<int> adjel         = element2adj[elID];
        int g                       = 0;
        std::vector<double> Vijk    = element2centroid[elID];
        double uijk                 = Ue[elID][variable];
        if(nadj>=9)
        {
            std::vector<double> bvec(nadj,0.0);
            std::vector<double> A_cm(nadj*9,0.0);

            std::set<int>::iterator its;

            for(its=adjel.begin();its!=adjel.end();its++)
            {

                double h00=0.0,h01=0.0,h02=0.0;
                double h11=0.0,h12=0.0,h22=0.0;
                double a=0.0,b=0.0,c=0.0,di=0.0;

                double uadj = 0.0;
                int adjid   = *its;

                if(ghosts.find(adjid)!=ghosts.end())
                {  
                    std::vector<double> Vadj    = ghostface_vrt[adjid];
                    uadj                        = ghosts[adjid][variable];
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                    di                          = sqrt(a*a+b*b+c*c);    
                }
                else
                {
                    std::vector<double> Vadj    = element2centroid[adjid]; 
                    uadj                        = Ue[adjid][variable];                   
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                    di                          = sqrt(a*a+b*b+c*c);
                    
                }

                h00 = 0.5*a*a; h01 = 1.0*a*b; h02 = 1.0*a*c;
                h11 = 0.5*b*b; h12 = 1.0*b*c;
                h22 = 0.5*c*c;

                A_cm[0*nadj+g] = (1.0/di)*a;
                A_cm[1*nadj+g] = (1.0/di)*b;
                A_cm[2*nadj+g] = (1.0/di)*c;

                A_cm[3*nadj+g] = (1.0/di)*h00;
                A_cm[4*nadj+g] = (1.0/di)*h01;
                A_cm[5*nadj+g] = (1.0/di)*h02;
                A_cm[6*nadj+g] = (1.0/di)*h11;
                A_cm[7*nadj+g] = (1.0/di)*h12;
                A_cm[8*nadj+g] = (1.0/di)*h22;

                bvec[g]        = 1.0/di*(uadj-uijk);
                
                g++;
            }

            // if(world_rank == 0)
            // {
            //     for(int h=0;h<nadj;h++)
            //     {
            //         for(int u=0;u<9;u++)
            //         {
            //             std::cout << A_cm[u*nadj+h] << " ";
            //         }
            //         std::cout << std::endl;
            //     }
            //     std::cout << std::endl;
            // }
            

            x = SolveQR_Lite(A_cm,nadj,9,bvec);

            // for(int u=0;u<9;u++)
            // {
            //     std::cout << x[u] << " ";
            // }
            // std::cout << std::endl;

            dudx_map[elID] = x;
            
        }
        else
        {
            std::vector<double> bvec(nadj,0.0);
            std::vector<double> A_cm(nadj*3,0.0);
            std::set<int>::iterator its;
            for(its=adjel.begin();its!=adjel.end();its++)
            {
                int adjid = *its;
                double h00=0.0,h01=0.0,h02=0.0;
                double h11=0.0,h12=0.0,h22=0.0;
                double a=0.0,b=0.0,c=0.0,di=0.0;
                double uadj=0.0;

                if(adjid>Nel)
                {  
                    std::vector<double> Vadj    = ghostface_vrt[adjid];
                    double uadj                 = ghosts[adjid][variable];
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                }
                else
                {
                    std::vector<double> Vadj    = element2centroid[adjid];
                    double uadj                 = Ue[adjid][variable];
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                    
                }

                di                              = sqrt(a*a+b*b+c*c);
                bvec[g]                         = 1.0/di*(uadj-uijk);
                A_cm[0*nadj+g]                  = (1.0/di)*a;
                A_cm[1*nadj+g]                  = (1.0/di)*b;
                A_cm[2*nadj+g]                  = (1.0/di)*c;
                
                g++;
            }

            x = SolveQR_Lite(A_cm,nadj,3,bvec);

            std::vector<double> xnew(9,0.0);
            xnew[0] = x[0];
            xnew[1] = x[1];
            xnew[2] = x[2];

            dudx_map[elID] = xnew;

        }        
   }
   
   return dudx_map;
}






std::map<int,std::vector<double> > ComputedUdx_LSQ_LS_US3D_Lite_Update(RepartitionObject* RePa,
                                                           std::map<int,std::vector<double> > ghosts,
                                                           int Nel,
                                                           int variable,
                                                           int nvariables,
                                                           MPI_Comm comm,
                                                           int extrap)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   // Get the rank of the process
   int world_rank;
   MPI_Comm_rank(comm, &world_rank);
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();
   //std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();

   std::map<int,std::vector<double> > Ue                = RePa->getElement2DataMap();
    //    RePa->AddStateVecForAdjacentElements(Ue,nvariables,comm);


   std::map<int,std::vector<double> > LocalVs           = RePa->getLocalVertsMap();
   std::vector<int> Loc_Elem                            = RePa->getLocElem();
   int nLoc_Elem                                        = Loc_Elem.size();
   std::map<int, std::vector<int> > Face2VertexMap      = RePa->getFace2VertexMap();
   std::map<int, std::vector<int> > Element2FaceMap     = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2FacesMap    = RePa->getElement2FacesMap();
   std::map<int, std::vector<int> > Element2VertexMap   = RePa->getElement2VertexMap();
   std::map<int, std::vector<int> > Element2ElementMap  = RePa->getElement2ElementMap();

   std::map<int,std::vector<double> > ghostface_vrt     = RePa->getGhostFaceData();
   std::map<int, std::set<int> > element2adj            = RePa->getExtendedAdjacencyData();
//    std::map<int, std::set<int> > element2adj         = RePa->getExtendedAdjacencyData(ghostface_vrt);
   
   std::map<int, std::vector<double> > element2centroid = RePa->GetElement2CentroidData();
   std::map<int, std::vector<int> >::iterator itmii;
   std::map<int,std::vector<double> > dudx_map;
   double d;
   int loc_vid,adjID,elID;
   int cou = 0;

   std::vector<double> Vc(3);
   std::vector<double> Vadj(3);
   int lid = 0;
   double u_ijk, u_po;
    
   std::map<int,std::vector<double> > vrt_collect;
   std::map<int,double> sol_collect;
   std::vector<std::vector<double> > face;
   double rdotn;
    
   std::vector<double> n0(3);
   std::vector<double> v0(3);
   std::vector<double> v1(3);

   int rr = 0;
   std::set<int> vert_scheme;
   int cc = 0;

   for(int i=0;i<nLoc_Elem;i++)
   {
        std::vector<double> x(9,0.0);
        int ghostelem               = 0;
        int elID                    = Loc_Elem[i];
        int nadj                    = element2adj[elID].size();
        //std::cout << " nadj " <<  nadj << std::endl;
        std::set<int> adjel         = element2adj[elID];
        int g                       = 0;
        std::vector<double> Vijk    = element2centroid[elID];
        double uijk                 = Ue[elID][variable];
        if(nadj>=9)
        {
            std::vector<double> bvec(nadj,0.0);
            std::vector<double> A_cm(nadj*9,0.0);

            std::set<int>::iterator its;

            for(its=adjel.begin();its!=adjel.end();its++)
            {

                double h00=0.0,h01=0.0,h02=0.0;
                double h11=0.0,h12=0.0,h22=0.0;
                double a=0.0,b=0.0,c=0.0,di=0.0;

                double uadj = 0.0;
                int adjid   = *its;

                if(ghosts.find(adjid)!=ghosts.end())
                {  
                    // std::vector<double> Vadj(3,0.0);
                    // for(int n=0;n<3;n++)
                    // {
                    //     Vadj[n] = Vijk[n];
                    // }
                    std::vector<double> Vadj    = ghostface_vrt[adjid];
                    uadj                        = uijk;//ghosts[adjid][variable];
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                    di                          = sqrt(a*a+b*b+c*c);    
                    // std::cout << "uadj " << uadj << " " << uijk << std::endl;
                }
                else
                {
                    std::vector<double> Vadj    = element2centroid[adjid]; 
                    uadj                        = Ue[adjid][variable];                   
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                    di                          = sqrt(a*a+b*b+c*c);
                    
                }

                h00 = 0.5*a*a; h01 = 1.0*a*b; h02 = 1.0*a*c;
                h11 = 0.5*b*b; h12 = 1.0*b*c;
                h22 = 0.5*c*c;

                A_cm[0*nadj+g] = (1.0/di)*a;
                A_cm[1*nadj+g] = (1.0/di)*b;
                A_cm[2*nadj+g] = (1.0/di)*c;

                A_cm[3*nadj+g] = (1.0/di)*h00;
                A_cm[4*nadj+g] = (1.0/di)*h01;
                A_cm[5*nadj+g] = (1.0/di)*h02;
                A_cm[6*nadj+g] = (1.0/di)*h11;
                A_cm[7*nadj+g] = (1.0/di)*h12;
                A_cm[8*nadj+g] = (1.0/di)*h22;

                bvec[g]        = 1.0/di*(uadj-uijk);
                
                g++;
            }

            // if(world_rank == 0)
            // {
            //     for(int h=0;h<nadj;h++)
            //     {
            //         for(int u=0;u<9;u++)
            //         {
            //             std::cout << A_cm[u*nadj+h] << " ";
            //         }
            //         std::cout << std::endl;
            //     }
            //     std::cout << std::endl;
            // }
            

            x = SolveQR_Lite(A_cm,nadj,9,bvec);

            // for(int u=0;u<9;u++)
            // {
            //     std::cout << x[u] << " ";
            // }
            // std::cout << std::endl;

            dudx_map[elID] = x;
            
        }
        else
        {
            std::vector<double> bvec(nadj,0.0);
            std::vector<double> A_cm(nadj*3,0.0);
            std::set<int>::iterator its;
            for(its=adjel.begin();its!=adjel.end();its++)
            {
                int adjid = *its;
                double h00=0.0,h01=0.0,h02=0.0;
                double h11=0.0,h12=0.0,h22=0.0;
                double a=0.0,b=0.0,c=0.0,di=0.0;
                double uadj=0.0;

                if(adjid>Nel)
                {  
                    // std::vector<double> Vadj(3,0.0);
                    // for(int n=0;n<3;n++)
                    // {
                    //     Vadj[n] = Vijk[n];
                    // }
                    std::vector<double> Vadj    = ghostface_vrt[adjid];
                    double uadj                 = uijk;//ghosts[adjid][variable];
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                }
                else
                {
                    std::vector<double> Vadj    = element2centroid[adjid];
                    double uadj                 = Ue[adjid][variable];
                    a                           = (Vadj[0] - Vijk[0]);
                    b                           = (Vadj[1] - Vijk[1]);
                    c                           = (Vadj[2] - Vijk[2]);
                    
                }

                di                              = sqrt(a*a+b*b+c*c);
                bvec[g]                         = 1.0/di*(uadj-uijk);
                A_cm[0*nadj+g]                  = (1.0/di)*a;
                A_cm[1*nadj+g]                  = (1.0/di)*b;
                A_cm[2*nadj+g]                  = (1.0/di)*c;
                
                g++;
            }

            x = SolveQR_Lite(A_cm,nadj,3,bvec);

            std::vector<double> xnew(9,0.0);
            xnew[0] = x[0];
            xnew[1] = x[1];
            xnew[2] = x[2];

            dudx_map[elID] = xnew;

        }        
   }
   
   return dudx_map;
}






std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Lite(RepartitionObject* RePa,
                                                             std::map<int,std::vector<double> > ghosts,
                                                             int Nel,
                                                             int variable,
                                                             int nvariable,
                                                             MPI_Comm comm,
                                                             int extrap)
{

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process;
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::map<int,std::vector<double> > Ue               = RePa->getElement2DataMap();
    std::map<int,std::vector<double> > dudx_map;
    std::vector<int> Loc_Elem                           = RePa->getLocElem();
    std::map<int, std::vector<double> > LocalVs         = RePa->getLocalVertsMap();
    std::map<int, std::vector<int> > Element2VertexMap  = RePa->getElement2VertexMap();
    std::map<int, std::vector<int> > Element2FacesMap   = RePa->getElement2FacesMap();
    std::map<int, std::vector<int> > Element2ElementMap = RePa->getElement2ElementMap();
    std::map<int, std::vector<int> > Face2VertexMap     = RePa->getFace2VertexMap();
    std::vector<int> Owned_Elem                         = RePa->getLocElem();
    int nLoc_Elem                                       = Owned_Elem.size();

    // std::map<int,std::vector<double> > ghostvrts = meshTopo->getGhostVerts();
    // std::map<int,std::vector<std::vector<double> > > vfvec = meshTopo->getVfacevector();
    
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;

    double d;
    int gvid,adjID,elID;
    int cou = 0;
    std::vector<double> Vc(3);
    std::vector<double> Vadj;
    int lid = 0;
    double u_ijk, u_po;
    //Array<double>* dudx = new Array<double>(nLoc_Elem,3);
    // std::map<int,int> LocElem2Nf = Pa->getLocElem2Nf();
    // std::map<int,int> LocElem2Nv = Pa->getLocElem2Nv();
    std::vector<std::vector<double> > face;
    double rdotn;
    
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    std::vector<double> n0(3);
    int notFound  = 0;
    int notFound2 = 0;
    int notFound3 = 0;
    int Found = 0;
    double orient0;
    int isZero = 0;
    int isNotZero = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    std::set<int> treated;
    int traceelem=0;
    int letmeknow=0;

    
    for(int q = 0;q < Owned_Elem.size();q++)
    {
        std::vector<double> x;
        int elID  = Owned_Elem[q];

        int NvPEl = Element2VertexMap[elID].size();
        int nadj  = Element2ElementMap[elID].size();
   
        std::vector<double> A_cm(nadj*3,0.0);
        std::vector<double> b(nadj,0.0);
        Array<double>* bvec     = new Array<double>(nadj,1);
        std::vector<double> Pijk(NvPEl*3);
        for(int k=0;k<Element2VertexMap[elID].size();k++)
        {
            gvid         = Element2VertexMap[elID][k];
            Pijk[k*3+0]  = LocalVs[gvid][0];
            Pijk[k*3+1]  = LocalVs[gvid][1];
            Pijk[k*3+2]  = LocalVs[gvid][2];
        }
   
        std::vector<double> Vijk     = ComputeCentroidCoord(Pijk,NvPEl);
   
        int t                        = 0;

        if(Ue.find(elID)!=Ue.end())
        {
            u_ijk                    = Ue[elID][variable];
        }
        else
        {
            notFound3++;
        }
        for(int j=0;j<nadj;j++)
        {
            int adjID       = Element2ElementMap[elID][j];
            int faceID      = Element2FacesMap[elID][j];
            // if(Element2FacesMap.find(elID)!=Element2FacesMap.end())
            // {
            //     faceID  = Element2FacesMap[elID][j];
            // }
            // else
            // {
            //     std::cout << "not here "<< std::endl;
            // }

            // if(treated.find(faceID)==treated.end())
            // {
            //     treated.insert(faceID);
            //     if(trace2elements.find(faceID)!=trace2elements.end())
            //     {
            //         Found++;
            //     }
            // }
        
            if(Element2VertexMap.find(adjID)!=Element2VertexMap.end()
                    && ghosts.find(adjID)==ghosts.end())
            {
                int nvpf        = Face2VertexMap[faceID].size();
                int nvpf_real   = Face2VertexMap[faceID].size();
                int nVadj       = Element2VertexMap[adjID].size();
                std::vector<double> Padj(nVadj*3,0.0);

                for(int k=0;k<nVadj;k++)
                {
                    int global_vid   = Element2VertexMap[adjID][k];

                    Padj[k*3+0] = LocalVs[global_vid][0];
                    Padj[k*3+1] = LocalVs[global_vid][1];
                    Padj[k*3+2] = LocalVs[global_vid][2];
                }
                
                Vadj = ComputeCentroidCoord(Padj,Element2VertexMap[adjID].size());
                
                d    = sqrt((Vadj[0]-Vijk[0])*(Vadj[0]-Vijk[0])+
                            (Vadj[1]-Vijk[1])*(Vadj[1]-Vijk[1])+
                            (Vadj[2]-Vijk[2])*(Vadj[2]-Vijk[2]));


                for(int s=0;s<3;s++)
                {
                    A_cm[s*nadj+t] = (1.0/d)*(Vadj[s]-Vijk[s]);
                    
                }

                u_po = Ue[adjID][variable];

                b[t] = (1.0/d)*(u_po-u_ijk);
                dist.push_back(d);
                t++;
                            
            }
            else
            {    
                int ghost_id    = adjID;
                int nvpf_real   = Face2VertexMap[faceID].size();

                std::vector<double> Vface(3,0.0);

                for(int u=0;u<nvpf_real;u++)
                {
                    int global_vid = Face2VertexMap[faceID][u];
                    Vface[0] = Vface[0]+LocalVs[global_vid][0];
                    Vface[1] = Vface[1]+LocalVs[global_vid][1];
                    Vface[2] = Vface[2]+LocalVs[global_vid][2];

                    std::vector<double> V(3,0.0);
                    V[0] = +LocalVs[global_vid][0];
                    V[1] = +LocalVs[global_vid][1];
                    V[2] = +LocalVs[global_vid][2];
                    face.push_back(V);
                }

                Vface[0] = Vface[0]/nvpf_real;
                Vface[1] = Vface[1]/nvpf_real;
                Vface[2] = Vface[2]/nvpf_real;

                std::vector<double> r0(3);
                r0[0] = (Vface[0]-Vijk[0]);
                r0[1] = (Vface[1]-Vijk[1]);
                r0[2] = (Vface[2]-Vijk[2]);

                if(nvpf_real==3) // triangle
                {
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[2][0]-face[0][0];
                    v1[1] = face[2][1]-face[0][1];
                    v1[2] = face[2][2]-face[0][2];
                }

                if(nvpf_real==4) // triangle
                {
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[3][0]-face[0][0];
                    v1[1] = face[3][1]-face[0][1];
                    v1[2] = face[3][2]-face[0][2];
                }


                std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);

                orient0   = DotVec3D(r0,n0);

                if(orient0<0.0)
                {
                    NegateVec3D(n0);
                }

                rdotn = DotVec3D(r0,n0);

                std::vector<double> reflect(3);
               
                reflect[0] = r0[0]-2.0*(rdotn)*n0[0];
                reflect[1] = r0[1]-2.0*(rdotn)*n0[1];
                reflect[2] = r0[2]-2.0*(rdotn)*n0[2];
                
                Vface[0] = Vface[0] - reflect[0];
                Vface[1] = Vface[1] - reflect[1];
                Vface[2] = Vface[2] - reflect[2];

                
                if(extrap==1)
                {
                    u_po = u_ijk;
                }
                else
                {
                    u_po = ghosts[ghost_id][variable];
                }


                d = sqrt((Vface[0]-Vijk[0])*(Vface[0]-Vijk[0])+
                         (Vface[1]-Vijk[1])*(Vface[1]-Vijk[1])+
                         (Vface[2]-Vijk[2])*(Vface[2]-Vijk[2]));
               
                for(int s=0;s<3;s++)
                {
                    A_cm[s*nadj+t] = (1.0/d)*(Vface[s]-Vijk[s]);
                }

                b[t] = (1.0/d)*(u_po-u_ijk);

                t++;
                
                face.clear();
            }    
        }

        // if(t!=4)
        // {
        //     std::cout << "t = " << t << " :: " << elID << ", " << nadj  << " " << Nel << std::endl;
        // }
        

    
        x = SolveQR_Lite(A_cm,nadj,3,b);
        // Array<double>* xold = SolveQR(A_cm.data(),nadj,3,bvec);
        dudx_map[elID] = x;

        

    }   

    
    //std::cout << " world_rank " << Found << " " << world_rank  << " traceelem :: " << traceelem  << "  letmeknow  " << letmeknow << std::endl;

    return dudx_map;
}


















std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Vrt_Lite(RepartitionObject* RePa,
                                                             std::map<int,std::vector<double> > Uval,
                                                             std::map<int,std::vector<double> > ghosts,
                                                             int Nel,
                                                             int variable,
                                                             int nvariable,
                                                             MPI_Comm comm,
                                                             int extrap)

{
    std::map<int,std::vector<double> > dudx_map;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process;
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    std::vector<int> Loc_Elem                           = RePa->getLocElem();
    std::map<int,std::vector<double> > LocalVs          = RePa->getLocalVertsMap();
    std::map<int, std::vector<int> > Element2VertexMap  = RePa->getElement2VertexMap();
    std::map<int, std::vector<int> > Element2FacesMap   = RePa->getElement2FacesMap();
    std::map<int, std::vector<int> > Element2ElementMap = RePa->getElement2ElementMap();
    std::map<int, std::vector<int> > Face2VertexMap     = RePa->getFace2VertexMap();
    std::vector<int> Owned_Elem                         = RePa->getLocElem();
    int nLoc_Elem                                       = Owned_Elem.size();
    //RePa->AddStateVecForAdjacentElements(Uval,nvariable,comm);
    std::map<int,std::map<int, double> > n2el_dist      = RePa->GetNode2ElementMapV2();
    std::map<int,std::vector<double> > node_val         = RePa->ReduceStateVecToVertices(n2el_dist,Uval,nvariable);
    std::map<int,std::map<int, double> > e2n_dist       = RePa->GetElement2NodeMap();
    

    std::map<int, double>::iterator e2dit;
    std::map<int,double> Uvrt_map;
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    std::vector<double> n0(3);
    std::vector<std::vector<double> > face;

    for(int q = 0;q < Owned_Elem.size();q++)
    {

        std::vector<std::vector<double> > vertices;
        std::vector<double> solution;

        int elid    = Owned_Elem[q];
        int NadjE   = Element2ElementMap[elid].size();
        int Nv      = Element2VertexMap[elid].size();
        std::map<int, double> el2node = e2n_dist[elid];
        int NadjV   = el2node.size();
        int nadj    = NadjE + NadjV; 
        std::vector<double> A_cm(nadj*3,0.0);
        std::vector<double> bvec(nadj,0.0);
        std::vector<double> Pijk(Nv*3);
        for(int k=0;k<Nv;k++)
        {
            int gvid     = Element2VertexMap[elid][k];
            Pijk[k*3+0]  = LocalVs[gvid][0];
            Pijk[k*3+1]  = LocalVs[gvid][1];
            Pijk[k*3+2]  = LocalVs[gvid][2];
        }
        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,Nv);
        double u_ijk = Uval[elid][variable];
        int t = 0;
        
        for(int k=0;k<NadjE;k++)
        {
            int adjID       = Element2ElementMap[elid][k];
            int faceID      = Element2FacesMap[elid][k];
            if(Element2VertexMap.find(adjID)!=Element2VertexMap.end()
                    && ghosts.find(adjID)==ghosts.end())
            {
                
                int Nvadj = Element2VertexMap[adjID].size();

                std::vector<double> Padjijk(Nvadj*3);
                for(int k=0;k<Nvadj;k++)
                {
                    int gvid        = Element2VertexMap[adjID][k];
                    Padjijk[k*3+0]  = LocalVs[gvid][0];
                    Padjijk[k*3+1]  = LocalVs[gvid][1];
                    Padjijk[k*3+2]  = LocalVs[gvid][2];
                }
                std::vector<double> Vadjijk = ComputeCentroidCoord(Padjijk,Nvadj);

                vertices.push_back(Vadjijk);
                solution.push_back(Uval[adjID][variable]);
                // double dx = Vadjijk[0]-Vijk[0];
                // double dy = Vadjijk[1]-Vijk[1];
                // double dz = Vadjijk[2]-Vijk[2];

                // double dist = sqrt(dx*dx+dy*dy+dz*dz);

                // A_cm[0*nadj+t]  = (1.0/dist)*(Vadjijk[0]-Vijk[0]);
                // A_cm[1*nadj+t]  = (1.0/dist)*(Vadjijk[1]-Vijk[1]);
                // A_cm[2*nadj+t]  = (1.0/dist)*(Vadjijk[2]-Vijk[2]);
                // double u_po     = Uval[adjID][variable];
                // b[t]            = (1.0/dist)*(u_po-u_ijk);

                t++;
            }
            else
            {
                int ghost_id    = adjID;
                int nvpf_real   = Face2VertexMap[faceID].size();
                std::vector<double> Vface(nvpf_real,0.0);

                for(int u=0;u<nvpf_real;u++)
                {
                    int global_vid = Face2VertexMap[faceID][u];
                    Vface[0] = Vface[0]+LocalVs[global_vid][0];
                    Vface[1] = Vface[1]+LocalVs[global_vid][1];
                    Vface[2] = Vface[2]+LocalVs[global_vid][2];

                    std::vector<double> V(3,0.0);
                    V[0] = +LocalVs[global_vid][0];
                    V[1] = +LocalVs[global_vid][1];
                    V[2] = +LocalVs[global_vid][2];
                    face.push_back(V);
                }

                Vface[0] = Vface[0]/nvpf_real;
                Vface[1] = Vface[1]/nvpf_real;
                Vface[2] = Vface[2]/nvpf_real;

                std::vector<double> r0(3);
                r0[0] = (Vface[0]-Vijk[0]);
                r0[1] = (Vface[1]-Vijk[1]);
                r0[2] = (Vface[2]-Vijk[2]);

                if(nvpf_real==3) // triangle
                {
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[2][0]-face[0][0];
                    v1[1] = face[2][1]-face[0][1];
                    v1[2] = face[2][2]-face[0][2];
                }

                if(nvpf_real==4) // triangle
                {
                    v0[0] = face[1][0]-face[0][0];
                    v0[1] = face[1][1]-face[0][1];
                    v0[2] = face[1][2]-face[0][2];

                    v1[0] = face[3][0]-face[0][0];
                    v1[1] = face[3][1]-face[0][1];
                    v1[2] = face[3][2]-face[0][2];
                }


                std::vector<double> n0 = ComputeSurfaceNormal(v0,v1);

                double orient0   = DotVec3D(r0,n0);

                if(orient0<0.0)
                {
                    NegateVec3D(n0);
                }

                double rdotn = DotVec3D(r0,n0);

                std::vector<double> reflect(3);
               
                reflect[0] = r0[0]-2.0*(rdotn)*n0[0];
                reflect[1] = r0[1]-2.0*(rdotn)*n0[1];
                reflect[2] = r0[2]-2.0*(rdotn)*n0[2];
                
                Vface[0] = Vface[0] - reflect[0];
                Vface[1] = Vface[1] - reflect[1];
                Vface[2] = Vface[2] - reflect[2];

                // double u_po = ghosts[ghost_id][variable];

                // double d = sqrt((Vface[0]-Vijk[0])*(Vface[0]-Vijk[0])+
                //                 (Vface[1]-Vijk[1])*(Vface[1]-Vijk[1])+
                //                 (Vface[2]-Vijk[2])*(Vface[2]-Vijk[2]));
               
                // for(int s=0;s<3;s++)
                // {
                //     A_cm[s*nadj+t] = (1.0/d)*(Vface[s]-Vijk[s]);
                // }

                // b[t] = (1.0/d)*(u_po-u_ijk);

                vertices.push_back(Vface);

                if(extrap==1)
                {
                    solution.push_back(u_ijk);
                }
                else
                {
                    solution.push_back(ghosts[ghost_id][variable]);
                }
                

                t++;
                
                face.clear();
            }
        }

        std::map<int, double>::iterator iter;

        for(iter=el2node.begin();iter!=el2node.end();iter++)
        {
            int gvid        = iter->first;
            double dist     = iter->second;
            double Uvrt     = node_val[gvid][variable];
            std::vector<double> Vvrt(3,0.0);
            Vvrt[0]        = LocalVs[gvid][0];
            Vvrt[1]        = LocalVs[gvid][1];
            Vvrt[2]        = LocalVs[gvid][2];

            // A_cm[0*nadj+t]  = (1.0/dist)*(Vvrt[0]-Vijk[0]);
            // A_cm[1*nadj+t]  = (1.0/dist)*(Vvrt[1]-Vijk[1]);
            // A_cm[2*nadj+t]  = (1.0/dist)*(Vvrt[2]-Vijk[2]);
            // b[t]            = (1.0/dist)*(Uvrt-u_ijk);

            vertices.push_back(Vvrt);
            solution.push_back(node_val[gvid][variable]);

            t++;
        }

        //std::cout << vertices.size() << " " << solution.size() << " " << nadj << std::endl;
        int npoints = vertices.size();

        
        for(int i=0;i<npoints;i++)
        {
            double Usolution = solution[i];
            double a = (vertices[i][0] - Vijk[0]);
            double b = (vertices[i][1] - Vijk[1]);
            double c = (vertices[i][2] - Vijk[2]);

            double dist = sqrt(a*a+b*b+c*c);

            A_cm[0*npoints+i] = (1.0/dist)*a;
            A_cm[1*npoints+i] = (1.0/dist)*b;
            A_cm[2*npoints+i] = (1.0/dist)*c;

            bvec[i]   = (1.0/dist)*(Usolution-u_ijk);
        }
        
        std::vector<double> x = SolveQR_Lite(A_cm,npoints,3,bvec);

        dudx_map[elid] = x;

        // if(world_rank==0)
        // {
        //     for(int q=0;q<nadj;q++)
        //     {
        //         for(int p=0;p<3;p++)
        //         {
        //             std::cout << A_cm[p*nadj+q] << " ";
        //         }
        //         std::cout << std::endl;
        //     }

        //     std::cout << std::endl;

        //     for(int q=0;q<nadj;q++)
        //     {
        //         std::cout << b[q] << std::endl;
        //     }
        //     std::cout << std::endl;
        // }


        
        //std::vector<double> x = SolveQR_Lite(A_cm,nadj,3,b);
        // Array<double>* xold = SolveQR(A_cm.data(),nadj,3,bvec);

        //dudx_map[elid] = x;
    }
    

    return dudx_map;
}