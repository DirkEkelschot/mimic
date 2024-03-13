
#include <chrono>
#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
#include "../../src/adapt_prismtetratrace.h"
#include "../../src/adapt_repartition.h"
#include "../../src/adapt_output_vtk.h"
#include "../../src/adapt_meshtopology_lite.h"
#include "../../src/adapt_gradreconstruct_lite.h"
#include "../../src/adapt_runparmmg.h"
#include "../../src/adapt_inputs.h"
#include "../../src/adapt_writeus3ddata.h"
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// typedef CGAL::Simple_cartesian<double> Kernel;
// typedef Kernel::Point_2 Point_2;
// typedef Kernel::Segment_2 Segment_2;



void SetAnalyticalSolution(RepartitionObject* tetra_repart,
                            int Nel_glob,
                            std::map<int,std::vector<double> > &U_map,
                            std::map<int,std::vector<double> > &Usol,
                            std::map<int,std::vector<double> > &gbMap,
                            std::map<int,std::vector<double> > &gbMap_dUdx,
                            std::map<int,std::vector<double> > &gbMap_dUdy,
                            std::map<int,std::vector<double> > &gbMap_dUdz,
                            std::map<int,std::vector<double> > &gbMap_dU2dx2,
                            std::map<int,std::vector<double> > &gbMap_dU2dxy,
                            std::map<int,std::vector<double> > &gbMap_dU2dxz)
{


    std::map<int,std::vector<int> >::iterator itmii;
    std::vector<double> v0(3,0.0);
    std::vector<double> v1(3,0.0);
    std::vector<double> n0(3,0.0);
    std::vector<double> n1(3,0.0);
    std::vector<double> r0(3,0.0);
    double rdotn, U, dUdx, dUdy, dUdz, dU2dx2, dU2dxy, dU2dxz;
    std::vector<std::vector<double> > face;
    int loc_vid;

    std::map<int,std::vector<int> > element2vert_map    = tetra_repart->getElement2VertexMap();
    std::map<int,std::vector<int> > element2element_map = tetra_repart->getElement2ElementMap();
    std::map<int,std::vector<int> > element2face_map    = tetra_repart->getElement2FacesMap();
    std::map<int,std::vector<int> > face2vert_map       = tetra_repart->getFace2VertexMap();
    std::map<int,std::vector<double> > locVs            = tetra_repart->getLocalVertsMap();

    for(itmii=element2vert_map.begin();itmii!=element2vert_map.end();itmii++)
    {
        int gid   = itmii->first;
        int nvrts = element2vert_map[gid].size();

        std::vector<double> Pv(nvrts*3);
        for(int q=0;q<nvrts;q++)
        {
            int gvid = element2vert_map[gid][q];
            Pv[q*3+0] = locVs[gvid][0];
            Pv[q*3+1] = locVs[gvid][1];
            Pv[q*3+2] = locVs[gvid][2];
        }
        
        std::vector<double> Vijk = ComputeCentroidCoord(Pv,nvrts);

        double vmx = Vijk[0];
        double vmy = Vijk[1];
        double vmz = Vijk[2];

        
        double S = 50.0;
        std::vector<double> Uvec(7,0.0);
        std::vector<double> Uvec2(1,0.0);
        Uvec2[0] = sin(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4));
        Uvec[0] = sin(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4));
        Uvec[1] = S*(vmy-0.4)*(vmz-0.4)*cos(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4));
        Uvec[2] = S*(vmx-0.4)*(vmz-0.4)*cos(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4));
        Uvec[3] = S*(vmx-0.4)*(vmy-0.4)*cos(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4));

        Uvec[4] = -(S*S)*(vmy-0.4)*(vmy-0.4)*(vmz-0.4)*(vmz-0.4)*sin(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4));
        Uvec[5] = S*(vmz-0.4)*(cos(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4))-S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4)*sin(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4)));
        Uvec[6] = S*(vmy-0.4)*(cos(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4))-S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4)*sin(S*(vmx-0.4)*(vmy-0.4)*(vmz-0.4)));
        U_map[gid] = Uvec;
        Usol[gid] = Uvec2;
        int nadj = element2element_map[gid].size();

        for(int j=0;j<nadj;j++)
        {
            int adjID = element2element_map[gid][j];
            
            if(adjID>=Nel_glob)
            {
                
                int fid = element2face_map[gid][j];
                int nvf = face2vert_map[fid].size();
                double* Fa = new double[nvf*3];
                std::vector<double> vface(3);

                for(int k=0;k<nvf;k++)
                {
                    int gvid = face2vert_map[fid][k];
                 
                    vface[0] = vface[0] + locVs[gvid][0];
                    vface[1] = vface[1] + locVs[gvid][1];
                    vface[2] = vface[2] + locVs[gvid][2];
                    
                    std::vector<double> V(3);
                    V[0]    = locVs[gvid][0];
                    V[1]    = locVs[gvid][1];
                    V[2]    = locVs[gvid][2];
                    face.push_back(V);
                }
                
                vface[0] = vface[0]/nvf;
                vface[1] = vface[1]/nvf;
                vface[2] = vface[2]/nvf;
                
                
                std::vector<double> r0(3);
                r0[0] = (vface[0]-Vijk[0]);
                r0[1] = (vface[1]-Vijk[1]);
                r0[2] = (vface[2]-Vijk[2]);
                
                v0[0] = face[1][0]-face[0][0];
                v0[1] = face[1][1]-face[0][1];
                v0[2] = face[1][2]-face[0][2];

                v1[0] = face[2][0]-face[0][0];
                v1[1] = face[2][1]-face[0][1];
                v1[2] = face[2][2]-face[0][2];
                
                n0 = ComputeSurfaceNormal(v0,v1);
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
                
                double vgx = vface[0] - reflect[0];
                double vgy = vface[1] - reflect[1];
                double vgz = vface[2] - reflect[2];
                
                
                U      = sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                dUdx   = S*(vgy-0.4)*(vgz-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                dUdy   = S*(vgx-0.4)*(vgz-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                dUdz   = S*(vgx-0.4)*(vgy-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                dU2dx2 = -(S*S)*(vgy-0.4)*(vgy-0.4)*(vgz-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                dU2dxy = S*(vgz-0.4)*(cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4))-S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)));
                dU2dxz = S*(vgy-0.4)*(cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4))-S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)));
                
                std::vector<double> gbMap_vec(1,0.0);
                gbMap_vec[0] = U;
                gbMap[adjID]=gbMap_vec;
                
                std::vector<double> gbMap_dUdx_vec(1,0.0);
                gbMap_dUdx_vec[0] = dUdx;
                gbMap_dUdx[adjID] = gbMap_dUdx_vec;
                std::vector<double> gbMap_dUdy_vec(1,0.0);
                gbMap_dUdy_vec[0] = dUdy;
                gbMap_dUdy[adjID] = gbMap_dUdy_vec;
                std::vector<double> gbMap_dUdz_vec(1,0.0);
                gbMap_dUdz_vec[0] = dUdz;
                gbMap_dUdz[adjID] = gbMap_dUdz_vec;
                
                std::vector<double> gbMap_dU2dx2_vec(1,0.0);
                gbMap_dU2dx2_vec[0] = dU2dx2;
                gbMap_dU2dx2[adjID]=gbMap_dU2dx2_vec;
                std::vector<double> gbMap_dU2dxy_vec(1,0.0);
                gbMap_dU2dxy_vec[0] = dU2dxy;
                gbMap_dU2dxy[adjID] = gbMap_dU2dxy_vec;
                std::vector<double> gbMap_dU2dxz_vec(1,0.0);
                gbMap_dU2dxz_vec[0] = dU2dxz;
                gbMap_dU2dxz[adjID] = gbMap_dU2dxz_vec;
                
                face.clear();
            }
            
            if(element2element_map.find(adjID)!=element2element_map.end())
            {
                int NadjadjID = element2element_map[adjID].size();
                
                int NvPadjadjID = element2vert_map[adjID].size();
                //double* PadjadjID = new double[NvPadjadjID*3];
                std::vector<double> PadjadjID(NvPadjadjID*3);
                for(int k=0;k<NvPadjadjID;k++)
                {
                    loc_vid          = element2vert_map[adjID][k];
                    PadjadjID[k*3+0] = locVs[loc_vid][0];
                    PadjadjID[k*3+1] = locVs[loc_vid][1];
                    PadjadjID[k*3+2] = locVs[loc_vid][2];
                }
                
                std::vector<double> VadjadjID   = ComputeCentroidCoord(PadjadjID,NvPadjadjID);
                                
                for(int k=0;k<NadjadjID;k++)
                {
                    int adjadjID = element2element_map[adjID][k];
                    
                    if(adjadjID>=Nel_glob)
                    {
                        int fid = element2face_map[adjID][k];
                        int nvf = face2vert_map[fid].size();
                        
                        std::vector<double> vface(3,0);
                        for(int k=0;k<nvf;k++)
                        {
                            int gvid = face2vert_map[fid][k];
                            vface[0] = vface[0] + locVs[gvid][0];
                            vface[1] = vface[1] + locVs[gvid][1];
                            vface[2] = vface[2] + locVs[gvid][2];
                            
                            std::vector<double> V(3);
                            V[0]    = locVs[gvid][0];
                            V[1]    = locVs[gvid][1];
                            V[2]    = locVs[gvid][2];
                            face.push_back(V);
                        }
                        
                        vface[0] = vface[0]/nvf;
                        vface[1] = vface[1]/nvf;
                        vface[2] = vface[2]/nvf;
                        
                        std::vector<double> r0(3);
                        r0[0] = (vface[0]-VadjadjID[0]);
                        r0[1] = (vface[1]-VadjadjID[1]);
                        r0[2] = (vface[2]-VadjadjID[2]);
                        
                        v0[0] = face[1][0]-face[0][0];
                        v0[1] = face[1][1]-face[0][1];
                        v0[2] = face[1][2]-face[0][2];

                        v1[0] = face[2][0]-face[0][0];
                        v1[1] = face[2][1]-face[0][1];
                        v1[2] = face[2][2]-face[0][2];
                        
                        n0 = ComputeSurfaceNormal(v0,v1);
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
                        
                        double vgx = vface[0] - reflect[0];
                        double vgy = vface[1] - reflect[1];
                        double vgz = vface[2] - reflect[2];
                        
                        U = sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                        dUdx   = S*(vgy-0.4)*(vgz-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                        dUdy   = S*(vgx-0.4)*(vgz-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                        dUdz   = S*(vgx-0.4)*(vgy-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                        dU2dx2 = -(S*S)*(vgy-0.4)*(vgy-0.4)*(vgz-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                        dU2dxy = S*(vgz-0.4)*(cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4))-S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)));
                        dU2dxz = S*(vgy-0.4)*(cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4))-S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)));
                        
                        
                        std::vector<double> gbMap_vec(1,0.0);
                        gbMap_vec[0] = U;
                        gbMap[adjadjID]=gbMap_vec;
                        
                        std::vector<double> gbMap_dUdx_vec(1,0.0);
                        gbMap_dUdx_vec[0] = dUdx;
                        gbMap_dUdx[adjadjID]=gbMap_dUdx_vec;
                        std::vector<double> gbMap_dUdy_vec(1,0.0);
                        gbMap_dUdy_vec[0] = dUdy;
                        gbMap_dUdy[adjadjID]=gbMap_dUdy_vec;
                        std::vector<double> gbMap_dUdz_vec(1,0.0);
                        gbMap_dUdz_vec[0] = dUdz;
                        gbMap_dUdz[adjadjID]=gbMap_dUdz_vec;
                        
                        std::vector<double> gbMap_dU2dx2_vec(1,0.0);
                        gbMap_dU2dx2_vec[0] = dU2dx2;
                        gbMap_dU2dx2[adjadjID]=gbMap_dU2dx2_vec;
                        std::vector<double> gbMap_dU2dxy_vec(1,0.0);
                        gbMap_dU2dxy_vec[0] = dU2dxy;
                        gbMap_dU2dxy[adjadjID]=gbMap_dU2dxy_vec;
                        std::vector<double> gbMap_dU2dxz_vec(1,0.0);
                        gbMap_dU2dxz_vec[0] = dU2dxz;
                        gbMap_dU2dxz[adjadjID]=gbMap_dU2dxz_vec;
                        
                        face.clear();
                    }
                    
                    if(element2element_map.find(adjadjID)!=element2element_map.end())
                    {
                        int NadjadjadjID = element2element_map[adjadjID].size();
//
                        int NvPadjadjadjID = element2vert_map[adjadjID].size();
                        //double* PadjadjadjID = new double[NvPadjadjadjID*3];
                        std::vector<double> PadjadjadjID(NvPadjadjadjID*3);
                        for(int k=0;k<NvPadjadjadjID;k++)
                        {
                            loc_vid     = element2vert_map[adjadjID][k];
                            PadjadjadjID[k*3+0] = locVs[loc_vid][0];
                            PadjadjadjID[k*3+1] = locVs[loc_vid][1];
                            PadjadjadjID[k*3+2] = locVs[loc_vid][2];
                        }

                        std::vector<double> VadjadjadjID   = ComputeCentroidCoord(PadjadjadjID,NvPadjadjadjID);
//
                        //delete[] PadjadjadjID;

                        for(int f=0;f<NadjadjadjID;f++)
                        {
                            int adjadjadjID = element2element_map[adjadjID][f];

                            if(adjadjadjID>=Nel_glob)
                            {
                                int fid = element2face_map[adjadjID][f];
                                int nvf = face2vert_map[fid].size();
                            
                                std::vector<double> vface(3);
                                for(int ve=0;ve<nvf;ve++)
                                {
                                    int gvid = face2vert_map[fid][ve];
                                    vface[0] = vface[0] + locVs[gvid][0];
                                    vface[1] = vface[1] + locVs[gvid][1];
                                    vface[2] = vface[2] + locVs[gvid][2];

                                    std::vector<double> V(3);
                                    V[0]    = locVs[gvid][0];
                                    V[1]    = locVs[gvid][1];
                                    V[2]    = locVs[gvid][2];
                                    face.push_back(V);
                                }

                                vface[0] = vface[0]/nvf;
                                vface[1] = vface[1]/nvf;
                                vface[2] = vface[2]/nvf;

                                std::vector<double> r0(3);
                                r0[0] = (vface[0]-VadjadjadjID[0]);
                                r0[1] = (vface[1]-VadjadjadjID[1]);
                                r0[2] = (vface[2]-VadjadjadjID[2]);

                                v0[0] = face[1][0]-face[0][0];
                                v0[1] = face[1][1]-face[0][1];
                                v0[2] = face[1][2]-face[0][2];

                                v1[0] = face[2][0]-face[0][0];
                                v1[1] = face[2][1]-face[0][1];
                                v1[2] = face[2][2]-face[0][2];

                                n0 = ComputeSurfaceNormal(v0,v1);
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

                                double vgx = vface[0] - reflect[0];
                                double vgy = vface[1] - reflect[1];
                                double vgz = vface[2] - reflect[2];

                                
                                U = sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                                dUdx   = S*(vgy-0.4)*(vgz-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                                dUdy   = S*(vgx-0.4)*(vgz-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                                dUdz   = S*(vgx-0.4)*(vgy-0.4)*cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                                dU2dx2 = -(S*S)*(vgy-0.4)*(vgy-0.4)*(vgz-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4));
                                dU2dxy = S*(vgz-0.4)*(cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4))-S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)));
                                dU2dxz = S*(vgy-0.4)*(cos(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4))-S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)*sin(S*(vgx-0.4)*(vgy-0.4)*(vgz-0.4)));

                                std::vector<double> gbMap_vec(1,0.0);
                                gbMap_vec[0] = U;
                                gbMap[adjadjadjID]=gbMap_vec;
                                
                                std::vector<double> gbMap_dUdx_vec(1,0.0);
                                gbMap_dUdx_vec[0] = dUdx;
                                gbMap_dUdx[adjadjadjID]=gbMap_dUdx_vec;
                                std::vector<double> gbMap_dUdy_vec(1,0.0);
                                gbMap_dUdy_vec[0] = dUdy;
                                gbMap_dUdy[adjadjadjID]=gbMap_dUdy_vec;
                                std::vector<double> gbMap_dUdz_vec(1,0.0);
                                gbMap_dUdz_vec[0] = dUdz;
                                gbMap_dUdz[adjadjadjID]=gbMap_dUdz_vec;
                                
                                std::vector<double> gbMap_dU2dx2_vec(1,0.0);
                                gbMap_dU2dx2_vec[0] = dU2dx2;
                                gbMap_dU2dx2[adjadjadjID]=gbMap_dU2dx2_vec;
                                std::vector<double> gbMap_dU2dxy_vec(1,0.0);
                                gbMap_dU2dxy_vec[0] = dU2dxy;
                                gbMap_dU2dxy[adjadjadjID]=gbMap_dU2dxy_vec;
                                std::vector<double> gbMap_dU2dxz_vec(1,0.0);
                                gbMap_dU2dxz_vec[0] = dU2dxz;
                                gbMap_dU2dxz[adjadjadjID]=gbMap_dU2dxz_vec;
                                face.clear();
                            }
                        }
                    }
                }
            }
        }
    }
}





int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j,k;
    clock_t t0_met = clock();

    int ier,opt;
    int debug = 1;

    
    std::map<int,std::vector<const char*> > grids;
    grids[10].push_back("inputs/grid101010.h5");
    grids[10].push_back("inputs/conn101010.h5");
    grids[25].push_back("inputs/grid252525.h5");
    grids[25].push_back("inputs/conn252525.h5");
    grids[50].push_back("inputs/grid505050.h5");
    grids[50].push_back("inputs/conn505050.h5");
    // Read in the inputs from metric.xml
    Inputs* inputs = ReadXmlFile(comm, "inputs/metric.xml");


    std::map<int,std::vector<const char*> >::iterator itg;
    for(itg=grids.begin();itg!=grids.end();itg++)
    {
        const char* fn_grid = itg->second[0];
        const char* fn_conn = itg->second[1];

        //===========================================================================
        // Read in the data from grid.h5/conn.h5/data.h5 in parallel using HDF5.
        // the outputted data in meshRead contains uniformly distributed data structures that will have to be partitioned.
        mesh* meshRead = ReadUS3DMesh(fn_conn,fn_grid,
                                        inputs->ReadFromStats,
                                        inputs->StateVar,
                                        comm,info);

        int Nel_loc = meshRead->ien.size();
        //Start building the trace object that contains the information regarding the prism-tetra interfaces.
        // It contains the unique vertex and face information.
        //====================================================================================
        

        clock_t start1, end1, start, end;
        double dur_max1,time_taken1;

        start1 = clock();

        PrismTetraTrace* pttrace = new PrismTetraTrace(comm, 
                                                    meshRead->element2rank, 
                                                    meshRead->ife,
                                                    meshRead->ifn,
                                                    meshRead->iet, 
                                                    meshRead->nElem, 
                                                    meshRead->nFace, 
                                                    meshRead->nVert);

        
        end1 = clock();
        time_taken1 = ( end1 - start1) / (double) CLOCKS_PER_SEC;
        MPI_Allreduce(&time_taken1, &dur_max1, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
            cout << "Time taken to broadcast boundary layer/tetra trace data is : " << fixed
            << dur_max1 << setprecision(16);
            cout << " sec " << endl;
        }


        std::map<int,std::map<int,int> > trace_elem     = pttrace->GetTrace();
        std::map<int,std::vector<int> > trace_verts     = pttrace->GetTraceVerts();
        std::map<int,int> unique_trace_verts2refmap     = pttrace->GetUniqueTraceVerts2RefMap();
        std::map<int,std::vector<int> > leftright_trace = pttrace->GetLeftRightElements();
        std::map<int,int> trace_ref                     = pttrace->GetTraceRef();
        FaceSetPointer FaceTraceRefs                    = pttrace->GetRefTraceFaceSet();

        //Filter out the tetrahedra and prisms into seperate maps from the IO data structures (meshRead).
        //=====================================================================================
        std::map<int,std::vector<int> > tetras_e2v,tetras_e2f,tetras_e2e;
        std::map<int,std::vector<double> > tetras_data;
        std::map<int,std::vector<int> > prisms_e2v,prisms_e2f,prisms_e2e;
        std::map<int,std::vector<double> > prisms_data;

        int ntetra      = meshRead->ntetra;
        int nprism      = meshRead->nprism;
        std::map<int,std::vector<int> >::iterator itmiv;
        int foundte     = 0;
        int foundpr     = 0;

        int ndata    = meshRead->interior.begin()->second.size();
        int nTrcFace = trace_verts.size();

        std::map<int,double> tracePrismData;
        std::map<int,double> traceTetraData;
        std::map<int,int> tetra2type;
        std::map<int,int> prism2type;   

        if(world_rank == 0)
        {
            std::cout << "Start filtering the element types..." << std::endl; 
        }

        for(itmiv=meshRead->ien.begin();itmiv!=meshRead->ien.end();itmiv++)
        {
            int elid   = itmiv->first;
            int eltype = meshRead->iet[elid];

            if(eltype == 2)
            {
                tetras_e2v[elid]  = itmiv->second;
                tetras_e2f[elid]  = meshRead->ief[elid];
                tetras_e2e[elid]  = meshRead->iee[elid];
                std::vector<double> Uvec(1,0.0);
                Uvec[0] = 0.0;
                tetras_data[elid] = Uvec;

                tetra2type[elid]  = eltype;
            }
            else
            {
                prisms_e2v[elid]  = itmiv->second;
                prisms_e2f[elid]  = meshRead->ief[elid];
                prisms_e2e[elid]  = meshRead->iee[elid];
                std::vector<double> Uvec(1,0.0);
                Uvec[0]           = 0.0;
                prisms_data[elid] = Uvec;

                prism2type[elid] = eltype;
                int nf = meshRead->ief[elid].size();

                for(int j=0;j<nf;j++)
                {
                    int faceID = meshRead->ief[elid][j];
                    if(trace_verts.find(faceID)!=trace_verts.end())
                    {
                        if(tracePrismData.find(elid)==tracePrismData.end())
                        {
                            tracePrismData[elid] = meshRead->interior[elid][1];
                        }
                    }
                }
            }   
        }

        bool tetra_ifn = true;
        bool prism_ifn = true;
        if(meshRead->elTypes[3]>0)
        {
            prism_ifn = false;
        }

        if(world_rank == 0)
        {
            std::cout << "Done filtering the element types..." << std::endl; 
        }
        
        //I am adding the prism elements and their data to the ghost map so that that data is in the boundaries data structures.
        std::map<int,double> tracePrismData_glob = AllGatherMap_T(tracePrismData,comm);
        
        std::map<int,double>::iterator itr;
        for(itr=tracePrismData_glob.begin();itr!=tracePrismData_glob.end();itr++)
        {
            int elid    = itr->first;
            double data = itr->second;

            std::vector<double> rowghost(2,0.0);
            rowghost[0] = 0.0;
            rowghost[1] = data;

            if(meshRead->ghost.find(elid)==meshRead->ghost.end())
            {
                meshRead->ghost[elid] = rowghost;
            }
        }
        int nLocTetra  = tetras_e2v.size();
        int nLocPrism  = prisms_e2v.size();
        int nElemsGlob_T = 0;
        int nElemsGlob_P = 0;
        MPI_Allreduce(&nLocTetra, &nElemsGlob_T, 1, MPI_INT, MPI_SUM, comm);
        MPI_Allreduce(&nLocPrism, &nElemsGlob_P, 1, MPI_INT, MPI_SUM, comm);

        

        if(world_rank == 0)
        {
            std::cout << "Done trace operation..." << std::endl; 
        }

        //=========END FILTERING OUT TETRA AND PRISMS FROM IO DATA STRUCTURES===============

        // we need to pass the number of verts per element in case the partition has no elements of this type.

    
        //  You can call it like this : start = time(NULL); 
        // in both the way start contain total time in seconds 
        // since the Epoch. 


        double dur_max,time_taken;

        start = clock();
        RepartitionObject* tetra_repart = new RepartitionObject(meshRead, 
                                                            tetras_e2v, 
                                                            tetras_e2f,
                                                            tetras_e2e,
                                                            tetra2type,
                                                            pttrace, 
                                                            tetras_data,
                                                            2,
                                                            tetra_ifn,
                                                            comm);

        

        tetras_e2v.clear();
        tetras_e2f.clear();
        tetras_e2e.clear();
        tetra2type.clear();
        tetras_data.clear();
                                                            
        end = clock();
        time_taken = ( end - start) / (double) CLOCKS_PER_SEC;

        
        MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
            cout << "Time taken to execute repartioning tetrahedera is : " << fixed 
            << dur_max << setprecision(16); 
            cout << " sec " << endl;
        }

        
        
        tetras_e2v.clear();tetras_e2f.clear();tetras_e2e.clear();

        tetra_repart->buildUpdatedVertexAndFaceNumbering(comm, 
                                                        pttrace, 
                                                        meshRead->ranges_id, 
                                                        meshRead->ranges_ref);

        tetra_repart->buildInteriorSharedAndBoundaryFaceMaps(comm, 
                                                            pttrace,
                                                            meshRead->ranges_id,
                                                            meshRead->ranges_ref);

        int nGlob = meshRead->nElem;
        std::map<int,std::vector<double> > Usol;
        std::map<int,std::vector<double> > U_map;
        std::map<int,std::vector<double> > dUdx_map;
        std::map<int,std::vector<double> > dUdy_map;
        std::map<int,std::vector<double> > dUdz_map;
        std::map<int,std::vector<double> > dU2dx2_map;
        std::map<int,std::vector<double> > dU2dxy_map;
        std::map<int,std::vector<double> > dU2dxz_map;

        std::map<int,std::vector<double> > gbMap;
        std::map<int,std::vector<double> > gbMap_dUdx;
        std::map<int,std::vector<double> > gbMap_dUdy;
        std::map<int,std::vector<double> > gbMap_dUdz;
        std::map<int,std::vector<double> > gbMap_dU2dx2;
        std::map<int,std::vector<double> > gbMap_dU2dxy;
        std::map<int,std::vector<double> > gbMap_dU2dxz;
        
        SetAnalyticalSolution(tetra_repart, nGlob,
                            U_map, Usol,
                            gbMap, gbMap_dUdx, gbMap_dUdy, gbMap_dUdz, gbMap_dU2dx2, gbMap_dU2dxy, gbMap_dU2dxz);
        // 0,1,2,3,4,5
        //std::cout << "compare size " << U_map.size() << " " << tetra_repart->getElement2VertexMap().size() << std::endl;
        tetra_repart->SetStateVec(U_map,7);

        start = clock();
        std::map<int,std::vector<double> > tetra_grad = ComputedUdx_LSQ_LS_US3D_Lite(tetra_repart,
                                                                                pttrace,
                                                                                gbMap,
                                                                                meshRead->nElem,
                                                                                0,
                                                                                1, 
                                                                                1,
                                                                                comm);

        std::map<int,std::vector<double> > tetra_grad_v2 = ComputedUdx_LSQ_US3D_Lite(tetra_repart,
                                                                                pttrace,
                                                                                Usol,
                                                                                gbMap,
                                                                                meshRead->nElem,
                                                                                0,
                                                                                1, 
                                                                                0,
                                                                                comm,
                                                                                0);
        
        std::map<int,std::vector<double> >::iterator iti;
        std::map<int,std::vector<double> > dudx_map;
        std::map<int,std::vector<double> > dudy_map;
        std::map<int,std::vector<double> > dudz_map;
        for(iti=tetra_grad_v2.begin();iti!=tetra_grad_v2.end();iti++)
        {
            int elid = iti->first;
            dudx_map[elid].push_back(iti->second[0]);
            dudy_map[elid].push_back(iti->second[1]);
            dudz_map[elid].push_back(iti->second[2]);
        }

        std::map<int,std::vector<double> > dU2dx2 = ComputedUdx_LSQ_US3D_Lite(tetra_repart,
                                                                                pttrace,
                                                                                dudx_map,
                                                                                gbMap,
                                                                                meshRead->nElem,
                                                                                0,
                                                                                1, 
                                                                                0,
                                                                                comm,
                                                                                0);
        
        // std::map<int,std::vector<double> > dU2dy2 = ComputedUdx_LSQ_US3D_Lite(tetra_repart,
        //                                                                         pttrace,
        //                                                                         dudy_map,
        //                                                                         gbMap,
        //                                                                         meshRead->nElem,
        //                                                                         0,
        //                                                                         1, 
        //                                                                         1,
        //                                                                         comm,
        //                                                                         0);

        // std::map<int,std::vector<double> > dU2dz2 = ComputedUdx_LSQ_US3D_Lite(tetra_repart,
        //                                                                         pttrace,
        //                                                                         dudz_map,
        //                                                                         gbMap,
        //                                                                         meshRead->nElem,
        //                                                                         0,
        //                                                                         1, 
        //                                                                         1,
        //                                                                         comm,
        //                                                                         0);
        
        std::map<int,std::vector<double> > tetra_grad_final;
        double grad_sum0 = 0.0;
        double grad_sum1 = 0.0;
        double grad_sum2 = 0.0;
        for(iti=dudx_map.begin();iti!=dudx_map.end();iti++)
        {
            int elid = iti->first;
            std::vector<double> row(9,0.0);
            row[0] = dudx_map[elid][0];
            //std::cout << "elID " << elid << " " << dudx_map[elid][0] << std::endl;
            row[1] = dudy_map[elid][0];
            row[2] = dudz_map[elid][0];


            grad_sum0 = grad_sum0 + dudx_map[elid][0];
            grad_sum1 = grad_sum1 + dudy_map[elid][0];
            grad_sum2 = grad_sum2 + dudz_map[elid][0];

            row[3] = dU2dx2[elid][0];
            row[4] = dU2dx2[elid][1];
            row[5] = dU2dx2[elid][2];
            tetra_grad_final[elid] = row;

        }

        std::cout << "double grad_sum " << grad_sum0 << " " << grad_sum1 << " " << grad_sum2 << std::endl;

        end = clock();
        time_taken = ( end - start) / (double) CLOCKS_PER_SEC;

        MPI_Allreduce(&time_taken, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
        if(world_rank == 0)
        {
            cout << "Time taken to execute calculating first and second order gradients : " << fixed 
            << dur_max << setprecision(16); 
            cout << " sec " << endl;
        }

        std::map<int,std::vector<double> >::iterator itmivd;
        int var = 0;
        //std::vector<std::vector<double> > L2_errors(tetra_grad.size(),0.0);
        double rec=0.0, exact=0.0;
        std::map<int,std::vector<double> > error_c;
        std::map<int,std::vector<double> > error_c_nobnd;
        int el = 0;
        std::set<int> el2consider;
        std::map<int,std::vector<int> > element2face_map    = tetra_repart->getElement2FacesMap();
        std::map<int,std::vector<int> > element2vert_map    = tetra_repart->getElement2VertexMap();

        std::map<int,std::vector<int> > face2element_map    = tetra_repart->getFace2ElementMap();
        std::map<int,std::vector<double> > locVerts         = tetra_repart->getLocalVertsMap();

        std::map<int,double> error_c2;
        double dUdx_error_cul = 0.0;
        double volume_tot = 0.0;
        for(itmivd=tetra_grad.begin();itmivd!=tetra_grad.end();itmivd++)
        {
            int elid        = itmivd->first;
            int nf          = element2face_map[elid].size();
            int bnd         = 0;

            std::vector<double> Pijk(4*3);
            
            for(int k=0;k<4;k++)
            {
                int global_vid   = element2vert_map[elid][k];
                Pijk[k*3+0]      = locVerts[global_vid][0];
                Pijk[k*3+1]      = locVerts[global_vid][1];
                Pijk[k*3+2]      = locVerts[global_vid][2];
            }
            double volume   = ComputeVolumeTetCell(Pijk);


            for(int j=0;j<nf;j++)
            {
                int fid = element2face_map[elid][j];
                int el0 = face2element_map[fid][0];
                int el1 = face2element_map[fid][1];

                if(el0>nGlob || el1 >nGlob)
                {

                    bnd = 1;
                }
            }

            int nvar = U_map[elid].size();
            std::vector<double> row(nvar-1,0.0);

            for(int i=1;i<nvar;i++)
            {
                
                exact             = U_map[elid][i];
                rec               = itmivd->second[i-1];
                //rec             = tetra_grad_final[elid][i-1];
                row[i-1]          = (rec-exact)*(rec-exact);
            }
            //std::cout << std::endl;
            error_c[elid]         = row;


            if(bnd == 0)
            {
                error_c_nobnd[elid] = row;
                el2consider.insert(elid);    
            }
            volume_tot = volume_tot + volume;
            el++;
        }
        int el_sum = 0;
        MPI_Allreduce(&el, &el_sum, 1, MPI_INT, MPI_SUM, comm);
        double volume_tot_sum = 0.0;
        MPI_Allreduce(&volume_tot, &volume_tot_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
        std::cout << "volume_tot " << volume_tot << " " << volume_tot_sum << std::endl;


        std::vector<double>L2errors(6,0.0);
        for(itmivd=error_c.begin();itmivd!=error_c.end();itmivd++)
        {
            int elid = itmivd->first;
            for(int i=0;i<6;i++)
            {
                L2errors[i] = L2errors[i] + error_c[elid][i];
            }
        }

        std::vector<double> L2errors_reduce(6,0.0);
        std::map<int,std::string > vnames;
        vnames[0]     = "dUdx";
        vnames[1]     = "dUdy";
        vnames[2]     = "dUdz";
        vnames[3]     = "dU2dx2";
        vnames[4]     = "dU2dxy";
        vnames[5]     = "dU2dxz";
        if(world_rank == 0)
        {
            std::cout << "[";
            for(int i=0;i<6;i++)
            {           
                MPI_Allreduce(&L2errors[i], &L2errors_reduce[i], 1, MPI_DOUBLE, MPI_SUM, comm);
                L2errors_reduce[i] = sqrt(L2errors_reduce[i]);

    
                if(world_rank == 0 && i < 5)
                {
                    std::cout << L2errors_reduce[i]/el_sum << ", ";
                }
            
                if(world_rank == 0 && i == 5)
                {
                    std::cout << L2errors_reduce[i]/el_sum;
                }
            }
            std::cout << "]" << std::endl;
        }
        

        std::map<int,std::vector<double> > loc_data_t       = tetra_repart->getElement2DataMap();
        std::map<int,std::vector<int> > gE2lV_t             = tetra_repart->getGlobalElement2LocalVertMap();
        std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
        std::map<int, std::vector<double> > LocalVertsMap_t = tetra_repart->getLocalVertsMap();
        std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();

        std::map<int,std::string > varnamesGrad;

        varnamesGrad[0]     = "dUdx";
        varnamesGrad[1]     = "dUdy";
        varnamesGrad[2]     = "dUdz";
        varnamesGrad[3]     = "dU2dx2";
        varnamesGrad[4]     = "dU2dxy";
        varnamesGrad[5]     = "dU2dxz";
        varnamesGrad[6]     = "dU2dy2";
        varnamesGrad[7]     = "dU2dyz";
        varnamesGrad[8]     = "dU2dz2";
        
        string filename_tg = "tetraRecon"+std::to_string(itg->first)+"_" + std::to_string(world_rank) + ".vtu";

        OutputTetraMeshPartitionVTK(comm,
                                filename_tg, 
                                Owned_Elem_t, 
                                gE2gV_t, 
                                tetra_grad, 
                                varnamesGrad, 
                                LocalVertsMap_t);

        string filename_te = "tetraExact"+std::to_string(itg->first)+"_" + std::to_string(world_rank) + ".vtu";

        std::map<int,std::string > varnamesGraExact;

        varnamesGraExact[0]     = "U";
        varnamesGraExact[1]     = "dUdx";
        varnamesGraExact[2]     = "dUdy";
        varnamesGraExact[3]     = "dUdz";
        varnamesGraExact[4]     = "dU2dx2";
        varnamesGraExact[5]     = "dU2dxy";
        varnamesGraExact[6]     = "dU2dxz";
        
        OutputTetraMeshPartitionVTK(comm,
                                filename_te, 
                                Owned_Elem_t, 
                                gE2gV_t, 
                                U_map, 
                                varnamesGraExact, 
                                LocalVertsMap_t);

        string filename_terror = "tetraError"+std::to_string(itg->first)+"_" + std::to_string(world_rank) + ".vtu";

        std::map<int,std::string > varnamesGrad_error;

        varnamesGrad_error[0]     = "dUdx_error";
        varnamesGrad_error[1]     = "dUdy_error";
        varnamesGrad_error[2]     = "dUdz_error";
        varnamesGrad_error[3]     = "dU2dx_error";
        varnamesGrad_error[4]     = "dU2dxy_error";
        varnamesGrad_error[5]     = "dU2dxz_error";

        OutputTetraMeshPartitionVTK(comm,
                                filename_terror, 
                                Owned_Elem_t, 
                                gE2gV_t, 
                                error_c, 
                                varnamesGrad_error, 
                                LocalVertsMap_t);

                               /* */
    }



    


/**/
    //=======================================================================================
    
    MPI_Finalize();
        
}

