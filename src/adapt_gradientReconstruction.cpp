#include "adapt_gradientReconstruction.h"

std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Lite(PartObject* RePa,
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
    
    std::map<int,std::vector<double> > dudx_map;

    std::map<int,std::vector<int> > m_Elem2Vert     = RePa->getElem2VertMap();
    std::map<int,std::vector<int> > m_Elem2Face     = RePa->getElem2FaceMap();
    std::map<int,std::vector<int> > m_Elem2Elem     = RePa->getElem2ElemMap();
    std::map<int,std::vector<double> > Ue           = RePa->getElem2DataMap();
    std::map<int,int> m_Elem2Rank                   = RePa->getElem2RankMap();
    std::map<int,std::vector<double> > m_VertCoords = RePa->getLocalVertsMap();
    std::map<int,std::vector<int> > m_Face2Vert     = RePa->getFace2VertMap();
    std::map<int,std::vector<int> > m_Face2Elem     = RePa->getFace2ElemMap();
    std::map<int,int> m_partMap                     = RePa->getPartMap();
    std::set<int> m_TraceVertsOnRank                = RePa->getTraceVertsOnRankMap();
    std::set<int> m_TraceFacesOnRank                = RePa->getTraceFacesOnRankMap();
    std::set<int> m_ElemSet                         = RePa->getLocalElemSet();

    
    std::vector<std::vector<double> > iee_dist;
    std::vector<double> dist;

    double d;
    int gvid,adjID,elID;
    int cou = 0;
    std::vector<double> Vc(3);
    std::vector<double> Vadj;
    int lid = 0;
    double u_ijk, u_po;

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

    std::set<int>::iterator its;
    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {
        std::vector<double> x;
        int elID  = *its;

        int NvPEl = m_Elem2Vert[elID].size();
        int nadj  = m_Elem2Elem[elID].size();
   
        std::vector<double> A_cm(nadj*3,0.0);
        std::vector<double> b(nadj,0.0);
        Array<double>* bvec     = new Array<double>(nadj,1);
        std::vector<double> Pijk(NvPEl*3);
        for(int k=0;k<m_Elem2Vert[elID].size();k++)
        {
            gvid         = m_Elem2Vert[elID][k];
            Pijk[k*3+0]  = m_VertCoords[gvid][0];
            Pijk[k*3+1]  = m_VertCoords[gvid][1];
            Pijk[k*3+2]  = m_VertCoords[gvid][2];
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
            int adjID       = m_Elem2Elem[elID][j];
            int faceID      = m_Elem2Face[elID][j];
        
            if(m_Elem2Vert.find(adjID)!=m_Elem2Vert.end()
                    && ghosts.find(adjID)==ghosts.end())
            {
                int nvpf        = m_Face2Vert[faceID].size();
                int nvpf_real   = m_Face2Vert[faceID].size();
                int nVadj       = m_Elem2Vert[adjID].size();
                std::vector<double> Padj(nVadj*3,0.0);

                for(int k=0;k<nVadj;k++)
                {
                    int global_vid   = m_Elem2Vert[adjID][k];

                    Padj[k*3+0] = m_VertCoords[global_vid][0];
                    Padj[k*3+1] = m_VertCoords[global_vid][1];
                    Padj[k*3+2] = m_VertCoords[global_vid][2];
                }
                
                Vadj = ComputeCentroidCoord(Padj,m_Elem2Vert[adjID].size());
                
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
                int nvpf_real   = m_Face2Vert[faceID].size();

                std::vector<double> Vface(3,0.0);

                for(int u=0;u<nvpf_real;u++)
                {
                    int global_vid = m_Face2Vert[faceID][u];
                    Vface[0] = Vface[0]+m_VertCoords[global_vid][0];
                    Vface[1] = Vface[1]+m_VertCoords[global_vid][1];
                    Vface[2] = Vface[2]+m_VertCoords[global_vid][2];

                    std::vector<double> V(3,0.0);
                    V[0] = +m_VertCoords[global_vid][0];
                    V[1] = +m_VertCoords[global_vid][1];
                    V[2] = +m_VertCoords[global_vid][2];
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

    
        x = SolveQR_Lite(A_cm,nadj,3,b);
        dudx_map[elID] = x;

        

    }   

    
    return dudx_map;
}



std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Vrt_Lite(PartObject* RePa,
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


    std::map<int,std::vector<int> > Elem2Vert                   = RePa->getElem2VertMap();
    std::map<int,std::vector<int> > Elem2Face                   = RePa->getElem2FaceMap();
    std::map<int,std::vector<int> > Elem2Elem                   = RePa->getElem2ElemMap();
    std::map<int,std::vector<int> > Face2Vert                   = RePa->getFace2VertMap();

    std::map<int,std::vector<double> > Ue                       = RePa->getElem2DataMap();
    std::map<int,int> m_Elem2Rank                               = RePa->getElem2RankMap();
    std::map<int,std::vector<double> > VertCoords               = RePa->getLocalVertsMap();
    std::set<int> m_ElemSet                                     = RePa->getLocalElemSet();

    std::map<int,std::map<int, double> > Vert2Elem              = RePa->GetNode2ElementMap();
    std::map<int,std::vector<double> > DataAtVert               = RePa->ReduceStateVecToVertices(Vert2Elem,Uval,nvariable);


    std::map<int,std::map<int, double> > Elem2ConnectedVertDist  = RePa->getElem2ConnectedVertMap();
    

    std::map<int, double>::iterator e2dit;
    std::map<int,double> Uvrt_map;
    std::vector<double> v0(3);
    std::vector<double> v1(3);
    std::vector<double> n0(3);
    std::vector<std::vector<double> > face;

    std::set<int>::iterator its;
    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {

        std::vector<std::vector<double> > vertices;
        std::vector<double> solution;

        int elid    = *its;
        int NadjE   = Elem2Elem[elid].size();
        int Nv      = Elem2Vert[elid].size();
        std::map<int, double> Elem2ConnectedVerts = Elem2ConnectedVertDist[elid];
        int NadjV   = Elem2ConnectedVerts.size();
        int nadj    = NadjE + NadjV; 
        std::vector<double> A_cm(nadj*3,0.0);
        std::vector<double> bvec(nadj,0.0);
        std::vector<double> Pijk(Nv*3);
        for(int k=0;k<Nv;k++)
        {
            int gvid     = Elem2Vert[elid][k];
            Pijk[k*3+0]  = VertCoords[gvid][0];
            Pijk[k*3+1]  = VertCoords[gvid][1];
            Pijk[k*3+2]  = VertCoords[gvid][2];
        }
        std::vector<double> Vijk = ComputeCentroidCoord(Pijk,Nv);
        double u_ijk = Uval[elid][variable];
        int t = 0;
        
        for(int k=0;k<NadjE;k++)
        {
            int adjID       = Elem2Elem[elid][k];
            int faceID      = Elem2Face[elid][k];
            if(Elem2Vert.find(adjID)!=Elem2Vert.end()
                    && ghosts.find(adjID)==ghosts.end())
            {
                
                int Nvadj = Elem2Vert[adjID].size();

                std::vector<double> Padjijk(Nvadj*3);
                for(int k=0;k<Nvadj;k++)
                {
                    int gvid        = Elem2Vert[adjID][k];
                    Padjijk[k*3+0]  = VertCoords[gvid][0];
                    Padjijk[k*3+1]  = VertCoords[gvid][1];
                    Padjijk[k*3+2]  = VertCoords[gvid][2];
                }
                std::vector<double> Vadjijk = ComputeCentroidCoord(Padjijk,Nvadj);

                vertices.push_back(Vadjijk);
                solution.push_back(Uval[adjID][variable]);
                t++;
            }
            else
            {
                int ghost_id    = adjID;
                int nvpf_real   = Face2Vert[faceID].size();
                std::vector<double> Vface(nvpf_real,0.0);

                for(int u=0;u<nvpf_real;u++)
                {
                    int global_vid = Face2Vert[faceID][u];
                    Vface[0] = Vface[0]+VertCoords[global_vid][0];
                    Vface[1] = Vface[1]+VertCoords[global_vid][1];
                    Vface[2] = Vface[2]+VertCoords[global_vid][2];

                    std::vector<double> V(3,0.0);
                    V[0] = +VertCoords[global_vid][0];
                    V[1] = +VertCoords[global_vid][1];
                    V[2] = +VertCoords[global_vid][2];
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

        for(iter=Elem2ConnectedVerts.begin();iter!=Elem2ConnectedVerts.end();iter++)
        {
            int gvid        = iter->first;
            double dist     = iter->second;
            double Uvrt     = DataAtVert[gvid][variable];
            std::vector<double> Vvrt(3,0.0);
            Vvrt[0]        = VertCoords[gvid][0];
            Vvrt[1]        = VertCoords[gvid][1];
            Vvrt[2]        = VertCoords[gvid][2];

            vertices.push_back(Vvrt);
            solution.push_back(DataAtVert[gvid][variable]);

            t++;
        }
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

    }
    

    return dudx_map;
}





std::map<int,std::vector<double> > ComputedUdx_LSQ_LS_US3D(PartObject* RePa,
                                                           std::map<int,std::vector<double> > ghosts,
                                                           int Nel,
                                                           int variable,
                                                           int nvariables,
                                                           MPI_Comm comm,
                                                           int extrap)
{

    std::map<int,std::vector<double> > dudx_map;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    std::map<int,std::vector<int> > Elem2Vert                   = RePa->getElem2VertMap();
    std::map<int,std::vector<int> > Elem2Face                   = RePa->getElem2FaceMap();
    std::map<int,std::vector<int> > Elem2Elem                   = RePa->getElem2ElemMap();
    std::map<int,std::vector<int> > Face2Vert                   = RePa->getFace2VertMap();

    std::map<int,std::vector<double> > Ue                       = RePa->getElem2DataMap();
    std::map<int,int> m_Elem2Rank                               = RePa->getElem2RankMap();
    std::map<int,std::vector<double> > VertCoords               = RePa->getLocalVertsMap();
    std::set<int> m_ElemSet                                     = RePa->getLocalElemSet();

    std::map<int,std::map<int, double> > Vert2Elem              = RePa->GetNode2ElementMap();
    std::map<int,std::vector<double> > GhostFaceVert            = RePa->getGhostFaceVert();

    std::map<int,std::map<int, double> > Elem2ConnectedVertDist = RePa->getElem2ConnectedVertMap();
    std::map<int,std::vector<double> > Elem2Centroid            = RePa->getElem2CentroidData();
    std::map<int,std::set<int> > Elem2AdjElem                   = RePa->getExtendedAdjacencyData();

   
    std::map<int, std::vector<int> >::iterator itmii;
    
    
    
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

    std::set<int>::iterator its;
    for(its=m_ElemSet.begin();its!=m_ElemSet.end();its++)
    {
            std::vector<double> x(9,0.0);
            int ghostelem               = 0;
            int elID                    = *its;
            int nadj                    = Elem2AdjElem[elID].size();
            std::set<int> adjel         = Elem2AdjElem[elID];
            int g                       = 0;
            std::vector<double> Vijk    = Elem2Centroid[elID];
            double uijk                 = Ue[elID][variable];

            /* */  
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
                        std::vector<double> Vadj    = GhostFaceVert[adjid];
                        uadj                        = ghosts[adjid][variable];
                        a                           = (Vadj[0] - Vijk[0]);
                        b                           = (Vadj[1] - Vijk[1]);
                        c                           = (Vadj[2] - Vijk[2]);
                        di                          = sqrt(a*a+b*b+c*c);    
                    }
                    else
                    {
                        std::vector<double> Vadj    = Elem2Centroid[adjid]; 
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


                x = SolveQR_Lite(A_cm,nadj,9,bvec);

                dudx_map[elID] = x;
                
            }
            else
            {
                // std::vector<double> bvec(nadj,0.0);
                // std::vector<double> A_cm(nadj*3,0.0);
                // std::set<int>::iterator its;
                // for(its=adjel.begin();its!=adjel.end();its++)
                // {
                //     int adjid = *its;
                //     double h00=0.0,h01=0.0,h02=0.0;
                //     double h11=0.0,h12=0.0,h22=0.0;
                //     double a=0.0,b=0.0,c=0.0,di=0.0;
                //     double uadj=0.0;

                //     if(adjid>Nel)
                //     {  
                //         std::vector<double> Vadj    = GhostFaceVert[adjid];
                //         double uadj                 = ghosts[adjid][variable];
                //         a                           = (Vadj[0] - Vijk[0]);
                //         b                           = (Vadj[1] - Vijk[1]);
                //         c                           = (Vadj[2] - Vijk[2]);
                //     }
                //     else
                //     {
                //         std::vector<double> Vadj    = Elem2Centroid[adjid];
                //         double uadj                 = Ue[adjid][variable];
                //         a                           = (Vadj[0] - Vijk[0]);
                //         b                           = (Vadj[1] - Vijk[1]);
                //         c                           = (Vadj[2] - Vijk[2]);
                        
                //     }

                //     di                              = sqrt(a*a+b*b+c*c);
                //     bvec[g]                         = 1.0/di*(uadj-uijk);
                //     A_cm[0*nadj+g]                  = (1.0/di)*a;
                //     A_cm[1*nadj+g]                  = (1.0/di)*b;
                //     A_cm[2*nadj+g]                  = (1.0/di)*c;
                    
                //     g++;
                // }

                // x = SolveQR_Lite(A_cm,nadj,3,bvec);

                // std::vector<double> xnew(9,0.0);
                // xnew[0] = x[0];
                // xnew[1] = x[1];
                // xnew[2] = x[2];

                // dudx_map[elID] = xnew;

            } 
                
    }
   /**/
    return dudx_map;
}