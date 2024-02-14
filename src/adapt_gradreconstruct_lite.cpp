#include "adapt_gradreconstruct_lite.h"









std::map<int,std::vector<double> > ComputedUdx_LSQ_US3D_Lite(RepartitionObject* RePa,
                                                             PrismTetraTrace* trace,
                                                             std::map<int,std::vector<double> > ghosts,
                                                             int Nel,
                                                             int variable,
                                                             int approxOrder,
                                                             MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    std::map<int,std::vector<double> > dudx_map;
    std::map<int,std::vector<double> > U                = RePa->getElement2DataMap();
    std::vector<int> Loc_Elem                           = RePa->getLocElem();
    std::map<int, std::vector<double> > LocalVs         = RePa->getLocalVertsMap();
    std::map<int, std::vector<int> > Element2VertexMap  = RePa->getElement2VertexMap();
    std::map<int, std::vector<int> > Element2FacesMap   = RePa->getElement2FacesMap();
    std::map<int, std::vector<int> > Element2ElementMap = RePa->getElement2ElementMap();
    std::map<int, std::vector<int> > Face2VertexMap     = RePa->getFace2VertexMap();
    std::map<int, std::vector<int> > Face2NVertexMap    = RePa->getFace2NVertexMap();
    std::vector<int> Owned_Elem                         = RePa->getLocElem();
    int nLoc_Elem                                       = Owned_Elem.size();
    std::map<int,std::map<int,int> > trace2elements     = trace->GetTrace();

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

        if(U.find(elID)!=U.end())
        {
            u_ijk                   = U[elID][variable];
            //std::cout << "u_ijk " << u_ijk << std::endl;
        }
        else
        {
            notFound3++;
        //std::cout << adjID << " adjID not found " << Nel << std::endl;

        }
        // std::cout << "nadj " << nadj << std::endl;
        for(int j=0;j<nadj;j++)
        {
            int adjID       = Element2ElementMap[elID][j];
            int faceID      = -1;
            if(Element2FacesMap.find(elID)!=Element2FacesMap.end())
            {
                faceID  = Element2FacesMap[elID][j];
            }
            else
            {
                letmeknow++;
            }

            if(treated.find(faceID)==treated.end())
            {
                treated.insert(faceID);
                if(trace2elements.find(faceID)!=trace2elements.end())
                {
                    Found++;
                }
            }
        
            if(adjID<Nel)
            {
                int nvpf        = Face2VertexMap[faceID].size();
                int nvpf_real   = Face2NVertexMap[faceID][0];

                if(Element2VertexMap.find(adjID)!=Element2VertexMap.end())
                {
                    int nVadj = Element2VertexMap[adjID].size();
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

                    u_po = U[adjID][1];

                    b[t] = (1.0/d)*(u_po-u_ijk);
                    dist.push_back(d);
                    t++;
                }
                else // Its on the trace between prisms and tets, The face is still internal but has a different element type on each side.
                {
                    traceelem++;
                }                
            }
            else if(approxOrder == 0)
            {    

                int ghost_id    = adjID;
                int nvpf        = Face2VertexMap[faceID].size();
                int nvpf_real   = Face2NVertexMap[faceID][0];

                std::vector<double> Vface(3,0.0);

                for(int u=0;u<nvpf_real;u++)
                {
                    int global_vid = Face2VertexMap[faceID][u];
                    Vface[0] = Vface[0]+LocalVs[global_vid][0];
                    Vface[1] = Vface[1]+LocalVs[global_vid][1];
                    Vface[2] = Vface[2]+LocalVs[global_vid][2];
                }

                Vface[0] = Vface[0]/nvpf_real;
                Vface[1] = Vface[1]/nvpf_real;
                Vface[2] = Vface[2]/nvpf_real;

                u_po = u_ijk;


                d = sqrt((Vface[0]-Vijk[0])*(Vface[0]-Vijk[0])+
                        (Vface[1]-Vijk[1])*(Vface[1]-Vijk[1])+
                        (Vface[2]-Vijk[2])*(Vface[2]-Vijk[2]));
               
                for(int s=0;s<3;s++)
                {
                    A_cm[s*nadj+t] = (1.0/d)*(Vface[s]-Vijk[s]);
                }

                b[t] = (1.0/d)*(u_po-u_ijk);
            }    
            else if(approxOrder == 1)
            {    

                int ghost_id    = adjID;
                int nvpf        = Face2VertexMap[faceID].size();
                int nvpf_real   = Face2NVertexMap[faceID][0];

                std::vector<double> Vface(3,0.0);

                for(int u=0;u<nvpf_real;u++)
                {
                    int global_vid = Face2VertexMap[faceID][u];
                    Vface[0] = Vface[0]+LocalVs[global_vid][0];
                    Vface[1] = Vface[1]+LocalVs[global_vid][1];
                    Vface[2] = Vface[2]+LocalVs[global_vid][2];
                }

                Vface[0] = Vface[0]/nvpf_real;
                Vface[1] = Vface[1]/nvpf_real;
                Vface[2] = Vface[2]/nvpf_real;

                u_po = ghosts[ghost_id][variable];


                d = sqrt((Vface[0]-Vijk[0])*(Vface[0]-Vijk[0])+
                        (Vface[1]-Vijk[1])*(Vface[1]-Vijk[1])+
                        (Vface[2]-Vijk[2])*(Vface[2]-Vijk[2]));
               
                for(int s=0;s<3;s++)
                {
                    A_cm[s*nadj+t] = (1.0/d)*(Vface[s]-Vijk[s]);
                }

                b[t] = (1.0/d)*(u_po-u_ijk);

                t++;
            }    
        }
    
        x = SolveQR_Lite(A_cm,nadj,3,b);
        
        dudx_map[elID] = x;
    }
    
    //std::cout << " world_rank " << Found << " " << world_rank  << " traceelem :: " << traceelem  << "  letmeknow  " << letmeknow << std::endl;

    return dudx_map;
}