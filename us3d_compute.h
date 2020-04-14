//#include "us3d_datastruct.h"

double ComputeDetJac(double *P0,double *P1,double *P2,double *P3)
{
    
    double *JP1 = new double[9];
    
    JP1[0] = (P1[0]-P0[0]); JP1[1] = (P2[0]-P0[0]); JP1[2] = (P3[0]-P0[0]);
    JP1[3] = (P1[1]-P0[1]); JP1[4] = (P2[1]-P0[1]); JP1[5] = (P3[1]-P0[1]);
    JP1[6] = (P1[2]-P0[2]); JP1[7] = (P2[2]-P0[2]); JP1[8] = (P3[2]-P0[2]);
    
    double DetJ = JP1[0]*(JP1[4]*JP1[8]-JP1[7]*JP1[5])
                 -JP1[1]*(JP1[3]*JP1[8]-JP1[6]*JP1[5])
                 +JP1[2]*(JP1[3]*JP1[7]-JP1[6]*JP1[4]);
    return DetJ;
}




double ComputeJ(double*P, int ElType){
    
    double J = 0.0;
    if (ElType==4)
    {
       double *P0  = new double[3];
       double *P1  = new double[3];
       double *P2  = new double[3];
       double *P3  = new double[3];
       
       P0[0]  = P[0*3+0]; P0[1]  = P[0*3+1]; P0[2] = P[0*3+2];
       P1[0]  = P[1*3+0]; P1[1]  = P[1*3+1]; P1[2] = P[1*3+2];
       P2[0]  = P[3*3+0]; P2[1]  = P[3*3+1]; P2[2] = P[3*3+2];
       P3[0]  = P[4*3+0]; P3[1]  = P[4*3+1]; P3[2] = P[4*3+2];
       double DJP0 = ComputeDetJac(P0,P1,P2,P3);

       P0[0]=P[1*3+0];P0[1]=P[1*3+1];P0[2]=P[1*3+2];
       P1[0]=P[2*3+0];P1[1]=P[2*3+1];P1[2]=P[2*3+2];
       P2[0]=P[0*3+0];P2[1]=P[0*3+1];P2[2]=P[0*3+2];
       P3[0]=P[5*3+0];P3[1]=P[5*3+1];P3[2]=P[5*3+2];
       double DJP1 = ComputeDetJac(P0,P1,P2,P3);
       
       // P2
       P0[0]=P[2*3+0];P0[1]=P[2*3+1];P0[2]=P[2*3+2];
       P1[0]=P[3*3+0];P1[1]=P[3*3+1];P1[2]=P[3*3+2];
       P2[0]=P[1*3+0];P2[1]=P[1*3+1];P2[2]=P[1*3+2];
       P3[0]=P[6*3+0];P3[1]=P[6*3+1];P3[2]=P[6*3+2];
       double DJP2 = ComputeDetJac(P0,P1,P2,P3);
       
       // P3
       P0[0]=P[3*3+0];P0[1]=P[3*3+1];P0[2]=P[3*3+2];
       P1[0]=P[0*3+0];P1[1]=P[0*3+1];P1[2]=P[0*3+2];
       P2[0]=P[2*3+0];P2[1]=P[2*3+1];P2[2]=P[2*3+2];
       P3[0]=P[7*3+0];P3[1]=P[7*3+1];P3[2]=P[7*3+2];
       double DJP3 = ComputeDetJac(P0,P1,P2,P3);
    
       // P4
       P0[0]=P[4*3+0];P0[1]=P[4*3+1];P0[2]=P[4*3+2];
       P1[0]=P[7*3+0];P1[1]=P[7*3+1];P1[2]=P[7*3+2];
       P2[0]=P[5*3+0];P2[1]=P[5*3+1];P2[2]=P[5*3+2];
       P3[0]=P[0*3+0];P3[1]=P[0*3+1];P3[2]=P[0*3+2];
       double DJP4 = ComputeDetJac(P0,P1,P2,P3);
       
       // P5
       P0[0]=P[5*3+0];P0[1]=P[5*3+1];P0[2]=P[5*3+2];
       P1[0]=P[4*3+0];P1[1]=P[4*3+1];P1[2]=P[4*3+2];
       P2[0]=P[6*3+0];P2[1]=P[6*3+1];P2[2]=P[6*3+2];
       P3[0]=P[1*3+0];P3[1]=P[1*3+1];P3[2]=P[1*3+2];
       double DJP5 = ComputeDetJac(P0,P1,P2,P3);
       
       // P6
       P0[0]=P[6*3+0];P0[1]=P[6*3+1];P0[2]=P[6*3+2];
       P1[0]=P[5*3+0];P1[1]=P[5*3+1];P1[2]=P[5*3+2];
       P2[0]=P[7*3+0];P2[1]=P[7*3+1];P2[2]=P[7*3+2];
       P3[0]=P[2*3+0];P3[1]=P[2*3+1];P3[2]=P[2*3+2];
       double DJP6 = ComputeDetJac(P0,P1,P2,P3);

       // P7
       P0[0]=P[7*3+0];P0[1]=P[7*3+1];P0[2]=P[7*3+2];
       P1[0]=P[6*3+0];P1[1]=P[6*3+1];P1[2]=P[6*3+2];
       P2[0]=P[4*3+0];P2[1]=P[4*3+1];P2[2]=P[4*3+2];
       P3[0]=P[3*3+0];P3[1]=P[3*3+1];P3[2]=P[3*3+2];
       double DJP7 = ComputeDetJac(P0,P1,P2,P3);
        
       J = (DJP0+DJP1+DJP2+DJP3+DJP4+DJP5+DJP6+DJP7)/8;

    }
    else if(ElType==6)
    {
        double *P0  = new double[3];
        double *P1  = new double[3];
        double *P2  = new double[3];
        double *P3  = new double[3];
        
        P0[0]  = P[0*3+0]; P0[1]  = P[0*3+1]; P0[2] = P[0*3+2];
        P1[0]  = P[1*3+0]; P1[1]  = P[1*3+1]; P1[2] = P[1*3+2];
        P2[0]  = P[2*3+0]; P2[1]  = P[2*3+1]; P2[2] = P[2*3+2];
        P3[0]  = P[3*3+0]; P3[1]  = P[3*3+1]; P3[2] = P[3*3+2];
        double DJP0 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[1*3+0]; P0[1]  = P[1*3+1]; P0[2] = P[1*3+2];
        P1[0]  = P[2*3+0]; P1[1]  = P[2*3+1]; P1[2] = P[2*3+2];
        P2[0]  = P[0*3+0]; P2[1]  = P[0*3+1]; P2[2] = P[0*3+2];
        P3[0]  = P[4*3+0]; P3[1]  = P[4*3+1]; P3[2] = P[4*3+2];
        double DJP1 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[2*3+0]; P0[1]  = P[2*3+1]; P0[2] = P[2*3+2];
        P1[0]  = P[0*3+0]; P1[1]  = P[0*3+1]; P1[2] = P[0*3+2];
        P2[0]  = P[1*3+0]; P2[1]  = P[1*3+1]; P2[2] = P[1*3+2];
        P3[0]  = P[5*3+0]; P3[1]  = P[5*3+1]; P3[2] = P[5*3+2];
        double DJP2 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[3*3+0]; P0[1]  = P[2*3+1]; P0[2] = P[2*3+2];
        P1[0]  = P[5*3+0]; P1[1]  = P[5*3+1]; P1[2] = P[5*3+2];
        P2[0]  = P[4*3+0]; P2[1]  = P[4*3+1]; P2[2] = P[4*3+2];
        P3[0]  = P[0*3+0]; P3[1]  = P[0*3+1]; P3[2] = P[0*3+2];
        double DJP3 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[4*3+0]; P0[1]  = P[2*3+1]; P0[2] = P[2*3+2];
        P1[0]  = P[3*3+0]; P1[1]  = P[3*3+1]; P1[2] = P[3*3+2];
        P2[0]  = P[5*3+0]; P2[1]  = P[5*3+1]; P2[2] = P[5*3+2];
        P3[0]  = P[1*3+0]; P3[1]  = P[1*3+1]; P3[2] = P[1*3+2];
        double DJP4 = ComputeDetJac(P0,P1,P2,P3);
        
        P0[0]  = P[5*3+0]; P0[1]  = P[5*3+1]; P0[2] = P[5*3+2];
        P1[0]  = P[4*3+0]; P1[1]  = P[4*3+1]; P1[2] = P[4*3+2];
        P2[0]  = P[3*3+0]; P2[1]  = P[3*3+1]; P2[2] = P[3*3+2];
        P3[0]  = P[2*3+0]; P3[1]  = P[2*3+1]; P3[2] = P[2*3+2];
        double DJP5 = ComputeDetJac(P0,P1,P2,P3);
        
        J = (DJP0+DJP1+DJP2+DJP3+DJP4+DJP5)/6;
    }
    return J;
}


inline double ComputeEdgeLength(Vert* v0, Vert* v1)
{
    return sqrt((v0->x - v1->x) * (v0->x - v1->x)+
                (v0->y - v1->y) * (v0->y - v1->y)+
                (v0->z - v1->z) * (v0->z - v1->z));
}



double ComputeVolumeHexCell(double *P)
{
    
    double L01=0.0;
    double L15=0.0;
    double L04=0.0;
    double L45=0.0;
    double L37=0.0;
    double L23=0.0;
    double L26=0.0;
    double L67=0.0;
    
    double b0,b1,b2,b3;
    double H12=0.0,H47=0.0,H30=0.0,H56=0.0;
    
    Vert v0;
    Vert v1;
    
    v0.x = P[0*3+0]; v1.x = P[1*3+0];
    v0.y = P[0*3+1]; v1.y = P[1*3+1];
    v0.z = P[0*3+2]; v1.z = P[1*3+2];
    
    L01 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[1*3+0]; v1.x = P[5*3+0];
    v0.y = P[1*3+1]; v1.y = P[5*3+1];
    v0.z = P[1*3+2]; v1.z = P[5*3+2];
    
    L15 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[1*3+0]; v1.x = P[2*3+0];
    v0.y = P[1*3+1]; v1.y = P[2*3+1];
    v0.z = P[1*3+2]; v1.z = P[2*3+2];
    
    H12 = ComputeEdgeLength(&v0,&v1);
    b0 = 0.5*L01*L15;
    double vol0 = 1.0/3.0*b0*H12;
    
    //==================================================
    
    
    v0.x = P[0*3+0]; v1.x = P[4*3+0];
    v0.y = P[0*3+1]; v1.y = P[4*3+1];
    v0.z = P[0*3+2]; v1.z = P[4*3+2];
    
    L04 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[4*3+0]; v1.x = P[5*3+0];
    v0.y = P[4*3+1]; v1.y = P[5*3+1];
    v0.z = P[4*3+2]; v1.z = P[5*3+2];
    
    L45 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[4*3+0]; v1.x = P[7*3+0];
    v0.y = P[4*3+1]; v1.y = P[7*3+1];
    v0.z = P[4*3+2]; v1.z = P[7*3+2];
    
    H47 = ComputeEdgeLength(&v0,&v1);
    b1 = 0.5*L04*L45;
    double vol1 = 1.0/3.0*b1*H47;
    
    //==================================================
    
    v0.x = P[3*3+0]; v1.x = P[7*3+0];
    v0.y = P[3*3+1]; v1.y = P[7*3+1];
    v0.z = P[3*3+2]; v1.z = P[7*3+2];
    
    L37 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[2*3+0]; v1.x = P[3*3+0];
    v0.y = P[2*3+1]; v1.y = P[3*3+1];
    v0.z = P[2*3+2]; v1.z = P[3*3+2];
    
    L23 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[3*3+0]; v1.x = P[0*3+0];
    v0.y = P[3*3+1]; v1.y = P[0*3+1];
    v0.z = P[3*3+2]; v1.z = P[0*3+2];
   
    H30 = ComputeEdgeLength(&v0,&v1);
    b2 = 0.5*L37*L23;
    double vol2 = 1.0/3.0*b2*H30;
    
    //==================================================
    
    v0.x = P[2*3+0]; v1.x = P[6*3+0];
    v0.y = P[2*3+1]; v1.y = P[6*3+1];
    v0.z = P[2*3+2]; v1.z = P[6*3+2];
    
    L26 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[6*3+0]; v1.x = P[7*3+0];
    v0.y = P[6*3+1]; v1.y = P[7*3+1];
    v0.z = P[6*3+2]; v1.z = P[7*3+2];
    
    L67 = ComputeEdgeLength(&v0,&v1);
    
    v0.x = P[5*3+0]; v1.x = P[6*3+0];
    v0.y = P[5*3+1]; v1.y = P[6*3+1];
    v0.z = P[5*3+2]; v1.z = P[6*3+2];
    
    H56 = ComputeEdgeLength(&v0,&v1);
    b3 = 0.5*L26*L67;
    double vol3 = 1.0/3.0*b3*H56;
    
    return vol0+vol1+vol2+vol3;
}

// This function outputs J as an array of 9 values where the matrix is defined as:

/*
 Jac = [J[0], J[1], J[2]
        J[3], J[4], J[5]
        J[6], J[7], J[8]]
*/
// J is computed using the 8-point isoparametric mapping for a hex. The 8-point rule should be sufficient since everything is linear anyways.

double* ComputeJAtCenter(double*P, int np)
{
    //Vert V;
    double * J = new double[9];
    for(int i=0;i<9;i++)
    {
        J[i] = 0.0;
    }
    int * ref = new int[np*3];
        
    // Allocate the arrays for the mapping function and its derivatives.
    double * N      = new double[np];
    double * dNdeta = new double[np];
    double * dNdmu  = new double[np];
    double * dNdksi = new double[np];
    
    for(int i=0;i<np;i++)
    {
        N[i] = 0.0;
        dNdeta[i] = 0.0;
        dNdmu[i]  = 0.0;
        dNdksi[i] = 0.0;
    }
    
    // Define the reference element as [-1,1]^3
    ref[0*3+0] = -1;    ref[0*3+1] = -1;    ref[0*3+2] = -1;
    ref[1*3+0] =  1;    ref[1*3+1] = -1;    ref[1*3+2] = -1;
    ref[2*3+0] =  1;    ref[2*3+1] =  1;    ref[2*3+2] = -1;
    ref[3*3+0] = -1;    ref[3*3+1] =  1;    ref[3*3+2] = -1;
            
    ref[4*3+0] = -1;    ref[4*3+1] = -1;    ref[4*3+2] = 1;
    ref[5*3+0] =  1;    ref[5*3+1] = -1;    ref[5*3+2] = 1;
    ref[6*3+0] =  1;    ref[6*3+1] =  1;    ref[6*3+2] = 1;
    ref[7*3+0] = -1;    ref[7*3+1] =  1;    ref[7*3+2] = 1;
     
    // We want to compute the Jacobian at the center of the cell
    // So we set eta, mu and ksi equal to 0 which is the center
    // of the reference cell.
    double eta = 0;
    double mu  = 0;
    double ksi = 0;
    double xphys = 0.0;double yphys = 0.0;double zphys = 0.0;

    // The basis functions for the 8 point isoparametric mapping for a hex
    // looks like the following.
    
    /*
    N[0] = 1.0/8.0*(1-eta)*(1-mu)*(1-ksi);
    N[1] = 1.0/8.0*(1+eta)*(1-mu)*(1-ksi);
    N[2] = 1.0/8.0*(1+eta)*(1+mu)*(1-ksi);
    N[3] = 1.0/8.0*(1-eta)*(1+mu)*(1-ksi);
    N[4] = 1.0/8.0*(1-eta)*(1-mu)*(1+ksi);
    N[5] = 1.0/8.0*(1+eta)*(1-mu)*(1+ksi);
    N[6] = 1.0/8.0*(1+eta)*(1+mu)*(1+ksi);
    N[7] = 1.0/8.0*(1-eta)*(1+mu)*(1+ksi);
    */
        
    for(int i = 0; i < np; i++)
    {
        N[i] = 1.0/8.0*(1+ref[i*3+0]*eta)*(1+ref[i*3+1]*mu)*(1+ref[i*3+2]*ksi);
        
        dNdeta[i] = (1.0/8.0 * (1+ref[i*3+1]*mu)  * (1+ref[i*3+2]*ksi))*ref[i*3+0];
        dNdmu[i]  = (1.0/8.0 * (1+ref[i*3+0]*eta) * (1+ref[i*3+2]*ksi))*ref[i*3+1];
        dNdksi[i] = (1.0/8.0 * (1+ref[i*3+0]*eta) * (1+ref[i*3+1]*mu))*ref[i*3+2];
                        
        xphys = xphys+N[i]*P[i*3+0];
        yphys = yphys+N[i]*P[i*3+1];
        zphys = zphys+N[i]*P[i*3+2];
            
        J[0] = J[0]+dNdeta[i]*P[i*3+0];
        J[1] = J[1]+dNdeta[i]*P[i*3+1];
        J[2] = J[2]+dNdeta[i]*P[i*3+2];
        
        J[3] = J[3]+dNdmu[i]*P[i*3+0];
        J[4] = J[4]+dNdmu[i]*P[i*3+1];
        J[5] = J[5]+dNdmu[i]*P[i*3+2];
            
        J[6] = J[6]+dNdksi[i]*P[i*3+0];
        J[7] = J[7]+dNdksi[i]*P[i*3+1];
        J[8] = J[8]+dNdksi[i]*P[i*3+2];
            
    }
    
    return J;
}

double ComputeDeterminantJ(double*P, int np)
{
    double* JP1 = ComputeJAtCenter(P, np);
    
    double DetJ = JP1[0]*(JP1[4]*JP1[8]-JP1[7]*JP1[5])
    -JP1[1]*(JP1[3]*JP1[8]-JP1[6]*JP1[5])
    +JP1[2]*(JP1[3]*JP1[7]-JP1[6]*JP1[4]);
    
    return DetJ;
    
}


double* ComputeDeterminantofJacobian(Array<double>* xcn, Array<int>* ien, int nloc, int offset, Array<double>* detJ)
{
    
    std::cout << "offset = " << offset << std::endl;
    int np = 8;
    double* Jac = new double[nloc];
    double* P = new double[np*3];
    int Vid;
    for(int i=0;i<nloc;i++)
    {
        np = 8;
        
        for(int j=0;j<np;j++)
        {
            Vid = ien->getVal(i,j+1)-1;
            
            P[j*3+0] = xcn->getVal(Vid+offset,0);
            P[j*3+1] = xcn->getVal(Vid+offset,1);
            P[j*3+2] = xcn->getVal(Vid+offset,2);
        }
        
        double dJ = ComputeDeterminantJ(P,8);
        
        detJ->setVal(i,0,dJ);
        Jac[i] = dJ;
        
        //std::cout << "================================" << std::endl;
        //std::cout << J[0] << " " << J[1] << " " << J[2] << std::endl;
        //std::cout << J[3] << " " << J[4] << " " << J[5] << std::endl;
        //std::cout << J[6] << " " << J[7] << " " << J[8] << std::endl;
        //std::cout << "================================" << std::endl;
    }
    return Jac;
}


double* ComputeVolumeCells(Array<double>* xcn, Array<int>* ien, MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    
    int nloc     = int(ien->nglob/world_size) + ( world_rank < ien->nglob%world_size );
    //  compute offset of rows for each proc;
    int offset   = world_rank*int(ien->nglob/world_size) + MIN(world_rank, ien->nglob%world_size);
    
    int Nelements = nloc;
    double * vol_cells = new double[Nelements];
    int np = 8;
    double* P = new double[np*3];
    
    
    double vhex = 0.0;
    
    int Vid;
    for(int i=0;i<Nelements;i++)
    {
        for(int j=0;j<np;j++)
        {
            Vid = ien->getVal(i+offset,j+1)-1;
            
            P[j*3+0] = xcn->getVal(Vid,0);
            P[j*3+1] = xcn->getVal(Vid,1);
            P[j*3+2] = xcn->getVal(Vid,2);
        }
        
        vhex = ComputeVolumeHexCell(P);
    
        vol_cells[i] = vhex;
        
    }
    
    return vol_cells;
    
}







double* ComputeVolumeCellsReducedToVerts(Array<double>* xcn, Array<int>* ien)
{
    
    
    int Nelements = ien->nglob;
    int Nnodes = xcn->nglob;
    
    double * vol_cells = new double[Nelements];
    double * vert_cnt   = new double[Nnodes];
    double * sum_vol    = new double[Nnodes];
    double * vol_verts  = new double[Nnodes];
    
    int np = 8;
    double* P = new double[np*3];
    
    double vhex = 0.0;
    
    int Vid;
    for(int i=0;i<Nelements;i++)
    {
        for(int j=0;j<np;j++)
        {
            Vid      = ien->getVal(i,j+1)-1;
            
            P[j*3+0] = xcn->getVal(Vid,0);
            P[j*3+1] = xcn->getVal(Vid,1);
            P[j*3+2] = xcn->getVal(Vid,2);
        }
        
        vhex = ComputeVolumeHexCell(P);
        
        for(int j=0;j<np;j++)
        {
            Vid = ien->getVal(i,j+1)-1;
            
            vert_cnt[Vid] = vert_cnt[Vid] + 1;
            sum_vol[Vid]  = sum_vol[Vid] + vhex;
        }
    }
    
    for(int i=0;i<Nnodes;i++)
    {
        vol_verts[i] = sum_vol[i]/vert_cnt[i];
    }
    
    return vol_verts;
}





/*
double* ComputeVolumeCells(Array<double>* xcn, Array<int>* ien)
{
    int Nelements = ien->nrow;
    double * vol_cells = new double[Nelements];
    int np = 8;
    double* P = new double[np*3];
    
    double L01=0.0;
    double L15=0.0;
    double L04=0.0;
    double L45=0.0;
    double L37=0.0;
    double L23=0.0;
    double L26=0.0;
    double L67=0.0;
    
    double b0,b1,b2,b3;
    double v0,v1,v2,v3,vhex;
    double H12,H47,H30,H56;
    
    int Vid;
    for(int i=0;i<Nelements;i++)
    {
        for(int j=0;j<np;j++)
        {
            Vid = ien->getVal(i,j+1)-1;
            P[j*3+0] = xcn->getVal(Vid,0);
            P[j*3+1] = xcn->getVal(Vid,1);
            P[j*3+2] = xcn->getVal(Vid,2);
        }
        
        L01 = sqrt((P[0*3+0]-P[1*3+0])*(P[0*3+0]-P[1*3+0])+
                   (P[0*3+1]-P[1*3+1])*(P[0*3+1]-P[1*3+1])+
                   (P[0*3+2]-P[1*3+2])*(P[0*3+2]-P[1*3+2]));
        
        L15 = sqrt((P[1*3+0]-P[5*3+0])*(P[1*3+0]-P[5*3+0])+
                   (P[1*3+1]-P[5*3+1])*(P[1*3+1]-P[5*3+1])+
                   (P[1*3+2]-P[5*3+2])*(P[1*3+2]-P[5*3+2]));
        
        H12 = sqrt((P[1*3+0]-P[2*3+0])*(P[1*3+0]-P[2*3+0])+
                   (P[1*3+1]-P[2*3+1])*(P[1*3+1]-P[2*3+1])+
                   (P[1*3+2]-P[2*3+2])*(P[1*3+2]-P[2*3+2]));
        
        
        b0 = 0.5*L01*L15;
        v0 = 1.0/3.0*b0*H12;
        
        L04 = sqrt((P[0*3+0]-P[4*3+0])*(P[0*3+0]-P[4*3+0])+
                   (P[0*3+1]-P[4*3+1])*(P[0*3+1]-P[4*3+1])+
                   (P[0*3+2]-P[4*3+2])*(P[0*3+2]-P[4*3+2]));
        
        L45 = sqrt((P[4*3+0]-P[5*3+0])*(P[4*3+0]-P[5*3+0])+
                   (P[4*3+1]-P[5*3+1])*(P[4*3+1]-P[5*3+1])+
                   (P[4*3+2]-P[5*3+2])*(P[4*3+2]-P[5*3+2]));
        
        H47 = sqrt((P[4*3+0]-P[7*3+0])*(P[4*3+0]-P[7*3+0])+
                   (P[4*3+1]-P[7*3+1])*(P[4*3+1]-P[7*3+1])+
                   (P[4*3+2]-P[7*3+2])*(P[4*3+2]-P[7*3+2]));
        
        b1 = 0.5*L04*L45;
        v1 = 1.0/3.0*b1*H47;
        
        L37 = sqrt((P[3*3+0]-P[7*3+0])*(P[3*3+0]-P[7*3+0])+
                   (P[3*3+1]-P[7*3+1])*(P[3*3+1]-P[7*3+1])+
                   (P[3*3+2]-P[7*3+2])*(P[3*3+2]-P[7*3+2]));
        
        L23 = sqrt((P[2*3+0]-P[3*3+0])*(P[2*3+0]-P[3*3+0])+
                   (P[2*3+1]-P[3*3+1])*(P[2*3+1]-P[3*3+1])+
                   (P[2*3+2]-P[3*3+2])*(P[2*3+2]-P[3*3+2]));
        
        H30 = sqrt((P[3*3+0]-P[0*3+0])*(P[3*3+0]-P[0*3+0])+
                   (P[3*3+1]-P[0*3+1])*(P[3*3+1]-P[0*3+1])+
                   (P[3*3+2]-P[0*3+2])*(P[3*3+2]-P[0*3+2]));
        
        b2 = 0.5*L37*L23;
        v2 = 1.0/3.0*b2*H30;
        
        L26 = sqrt((P[2*3+0]-P[6*3+0])*(P[2*3+0]-P[6*3+0])+
                   (P[2*3+1]-P[6*3+1])*(P[2*3+1]-P[6*3+1])+
                   (P[2*3+2]-P[6*3+2])*(P[2*3+2]-P[6*3+2]));
        
        L67 = sqrt((P[6*3+0]-P[7*3+0])*(P[6*3+0]-P[7*3+0])+
                   (P[6*3+1]-P[7*3+1])*(P[6*3+1]-P[7*3+1])+
                   (P[6*3+2]-P[7*3+2])*(P[6*3+2]-P[7*3+2]));
        
        H56 = sqrt((P[5*3+0]-P[6*3+0])*(P[5*3+0]-P[6*3+0])+
                   (P[5*3+1]-P[6*3+1])*(P[5*3+1]-P[6*3+1])+
                   (P[5*3+2]-P[6*3+2])*(P[5*3+2]-P[6*3+2]));
        
        b3 = 0.5*L26*L67;
        v3 = 1.0/3.0*b3*H56;
        vhex = v0+v1+v2+v3;
    
        vol_cells[i] = vhex;
        
    }
    return vol_cells;
}
*/
