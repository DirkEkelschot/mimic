#include "adapt_metricdistance.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>


MetricDistance::MetricDistance()
{

}

MetricDistance::MetricDistance(std::vector<std::vector<double> > edge, std::vector<std::vector<double> > metrics)
{
    double Ax = edge[0][0];
    double Ay = edge[0][1];
    double Az = edge[0][2];

    double Bx = edge[1][0];
    double By = edge[1][1];
    double Bz = edge[1][2];

    std::vector<double> S(3,0.0);

    double dSxdxi = Bx-Ax;
    double dSydxi = By-Ay;
    double dSzdxi = Bz-Az;

    gsl_interp_accel *accel;
    gsl_interp *interpX;
    gsl_interp *interpY;

    std::vector<double> V(3,0.0);
    std::vector<double> tmp(3,0.0);
    V[0] = Bx - Ax;
    V[1] = By - Ay;
    V[2] = Bz - Az;

    // std::cout << V[0] << " " << V[1] << " " << V[2] << std::endl;

    int nq = 40;

    std::vector<std::vector<double> > Mi(3);
    for(int k=0;k<3;k++)
    {
        std::vector<double> rowMi(3,0.0);
        Mi[k] = rowMi;
    }

    double porder   = 2;
    double lp       = 0.0;

    for(int i=0;i<nq;i++)
    {
        double xq  = xquad[i];
        double wq  = wquad[i];

        double jac = 1.0;

        S[0] = Ax*(1-xq)+Bx*xq;
        S[1] = Ay*(1-xq)+By*xq;
        S[2] = Az*(1-xq)+Bz*xq;

        // Linear interpolation of the metric field
        Mi[0][0] = metrics[0][0]*(1-xq)+metrics[1][0]*xq;   Mi[0][1] = metrics[0][1]*(1-xq)+metrics[1][1]*xq;   Mi[0][2] = metrics[0][2]*(1-xq)+metrics[1][2]*xq;
        Mi[1][0] = metrics[0][1]*(1-xq)+metrics[1][1]*xq;   Mi[1][1] = metrics[0][3]*(1-xq)+metrics[1][3]*xq;   Mi[1][2] = metrics[0][4]*(1-xq)+metrics[1][4]*xq;
        Mi[2][0] = metrics[0][2]*(1-xq)+metrics[1][2]*xq;   Mi[2][1] = metrics[0][4]*(1-xq)+metrics[1][4]*xq;   Mi[2][2] = metrics[0][5]*(1-xq)+metrics[1][5]*xq;

        // tmp = M*S
        tmp[0] = Mi[0][0]*V[0] + Mi[0][1]*V[1] + Mi[0][2]*V[2];
        tmp[1] = Mi[1][0]*V[0] + Mi[1][1]*V[1] + Mi[1][2]*V[2];
        tmp[2] = Mi[2][0]*V[0] + Mi[2][1]*V[1] + Mi[2][2]*V[2];

        // tmp[0] = Mi[0][0]*S[0] + Mi[0][1]*S[1] + Mi[0][2]*S[2];
        // tmp[1] = Mi[1][0]*S[0] + Mi[1][1]*S[1] + Mi[1][2]*S[2];
        // tmp[2] = Mi[2][0]*S[0] + Mi[2][1]*S[1] + Mi[2][2]*S[2];

        double lpNorm = pow(tmp[0],porder)+pow(tmp[1],porder)+pow(tmp[2],porder);

        lp = lp + wq*pow(lpNorm,1.0/porder)*jac;
    
    }

    std::cout << "lp " << lp << std::endl;

    

}

MetricDistance::~MetricDistance()
{
}

double MetricDistance::GetLength()
{

}
