#include "adapt_io.h"

US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, MPI_Comm comm, MPI_Info info)
{
    
    US3D* us3d = new US3D;
    ParArray<double>* xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);

    ParArray<int>* ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
    ParArray<int>* iee = ReadDataSetFromFileInParallel<int>(fn_conn,"iee",comm,info);

    Array<int>* ifn = ReadDataSetFromFile<int>(fn_grid,"ifn");
    Array<int>* ife = ReadDataSetFromFile<int>(fn_conn,"ife");
    
    int Nel = ien->getNglob();
    
    ParArray<double>* interior = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","interior",0,Nel,comm,info);
    Array<double>* ghost = ReadUS3DGhostCellsFromRun<double>(fn_data,"run_1","interior",Nel);
    
    int i,j;
    int nglob = ien->getNglob();
    int nrow  = ien->getNrow();
    int ncol  = ien->getNcol()-1;
    //
    ParArray<int>* ien_copy = new ParArray<int>(nglob,ncol,comm);
    //
    for(i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol;j++)
        {
            ien_copy->setVal(i,j,ien->getVal(i,j+1)-1);
        }
    }
    delete ien;
    //
    int ncol_ief = ief->getNcol()-1;
    ParArray<int>* ief_copy = new ParArray<int>(nglob,ncol_ief,comm);

    for(i=0;i<nrow;i++)
    {
        for(j=0;j<ncol_ief;j++)
        {
            ief_copy->setVal(i,j,fabs(ief->getVal(i,j+1))-1);
        }
    }
    delete ief;
    
    int ncol_iee = iee->getNcol()-1;
    ParArray<int>* iee_copy = new ParArray<int>(nglob,6,comm);
    
    for(i=0;i<nrow;i++)
    {
        for(j=0;j<ncol_iee;j++)
        {
            iee_copy->setVal(i,j,iee->getVal(i,j+1)-1);
        }
    }
    delete iee;
    
    int nrow_ifn = ifn->getNrow();

    int ncol_ifn = 4;
    Array<int>* ifn_copy = new Array<int>(nrow_ifn,ncol_ifn);
    for(i=0;i<nrow_ifn;i++)
    {
        for(j=0;j<ncol_ifn;j++)
        {
            ifn_copy->setVal(i,j,ifn->getVal(i,j+1)-1);
        }
    }
    delete ifn;
    
    int nrow_ife = ife->getNrow();
    int ncol_ife = 2;
    Array<int>* ife_copy = new Array<int>(nrow_ife,ncol_ife);
    for(i=0;i<nrow_ife;i++)
    {
        for(j=0;j<ncol_ife;j++)
        {
            ife_copy->setVal(i,j,ife->getVal(i,j)-1);
        }
    }
    delete ife;
    
    us3d->xcn = xcn;
    
    us3d->ien = ien_copy;
    us3d->ief = ief_copy;
    us3d->iee = iee_copy;
    
    us3d->ifn = ifn_copy;
    us3d->ife = ife_copy;
    
    us3d->interior = interior;
    us3d->ghost    = ghost;
    
    return us3d;
}
