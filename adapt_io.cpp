#include "adapt_io.h"
#include "adapt_output.h"
US3D* ReadUS3DData(const char* fn_conn, const char* fn_grid, const char* fn_data, MPI_Comm comm, MPI_Info info)
{
    int size,rank;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    MPI_Comm_rank(comm, &rank);
    
    US3D* us3d = new US3D;
    ParArray<double>* xcn = ReadDataSetFromFileInParallel<double>(fn_grid,"xcn",comm,info);

    ParArray<int>* ien = ReadDataSetFromFileInParallel<int>(fn_conn,"ien",comm,info);
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);
    ParArray<int>* iee = ReadDataSetFromFileInParallel<int>(fn_conn,"iee",comm,info);

    Array<int>* ifn = ReadDataSetFromFile<int>(fn_grid,"ifn");
    Array<int>* ife = ReadDataSetFromFile<int>(fn_conn,"ife");
    
    int Nel = ien->getNglob();
    
    ParArray<double>* interior  = ReadDataSetFromRunInFileInParallel<double>(fn_data,"run_1","interior",0,Nel,comm,info);
    Array<double>* ghost        = ReadUS3DGhostCellsFromRun<double>(fn_data,"run_1","interior",Nel);
    
    Array<int>*    zdefs        = ReadDataSetFromGroupFromFile<int>(fn_grid,"zones","zdefs");
    Array<char>*  znames        = ReadDataSetFromGroupFromFile<char>(fn_grid,"zones","znames");
    
    // Collect boundary data;
    std::vector<int> bnd_m;
    int t=0;
    for(int i=4;i<zdefs->getNrow();i++)
    {
        bnd_m.push_back(zdefs->getVal(i,3));
    }
    bnd_m.push_back(zdefs->getVal(zdefs->getNrow()-1,4));
    
    PlotBoundaryData(znames,zdefs,comm);
    
    int nBnd = zdefs->getNrow()-3;
    int* bnd_map = new int[zdefs->getNrow()-3];
    for(int i=3;i<zdefs->getNrow();i++)
    {
        bnd_map[i-3] = zdefs->getVal(i,3)-1;
        //bnd_m.push_back(zdefs->getVal(i,3));
    }
    int fint = bnd_map[0];
    if(rank == 0)
    {
        for(int i=0;i<zdefs->getNrow()-3;i++)
        {
            std::cout << "bnd_map " << bnd_map[i] << " " << nBnd << std::endl;
        }
    }
    
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
    Array<int>* ifn_ref  = new Array<int>(nrow_ifn,1);
    int ref;
    std::ofstream myfile20;
    myfile20.open("ifn_ref.dat");
    std::map<std::set<int>,int> tria_ref_map;
    std::set<int> tria0;
    std::set<int> tria1;
    std::vector<int*> tria_ref;
    
    for(i=0;i<nrow_ifn;i++)
    {
        for(j=0;j<ncol_ifn;j++)
        {
            ifn_copy->setVal(i,j,ifn->getVal(i,j+1)-1);
            
            if(i<fint) // identify the internal face;
            {
                ref = 0;
                ifn_ref->setVal(i,0,ref);
                myfile20 << ref << std::endl;
            }
            else // identify the boundary interface and based on bnd_map, determine the reference value.
            {
                ref = FindBoundaryID(bnd_map,nBnd,i)+1;
                ifn_ref->setVal(i,0,ref);
                myfile20 << ref << std::endl;
            }
        }
        
        tria0.insert(ifn->getVal(i,0+1)-1);
        tria0.insert(ifn->getVal(i,1+1)-1);
        tria0.insert(ifn->getVal(i,2+1)-1);
        
        tria1.insert(ifn->getVal(i,0+1)-1);
        tria1.insert(ifn->getVal(i,1+1)-1);
        tria1.insert(ifn->getVal(i,3+1)-1);
        
        if(tria_ref_map.find(tria0)==tria_ref_map.end() && ref!=0)
        {
            
            tria_ref_map[tria0] = ref;
            int* tria_tmp   = new int[4];
            tria_tmp[0]     = ifn->getVal(i,0+1)-1;
            tria_tmp[1]     = ifn->getVal(i,1+1)-1;
            tria_tmp[2]     = ifn->getVal(i,2+1)-1;
            tria_tmp[3]     = ref;
            tria_ref.push_back(tria_tmp);
        }
        if(tria_ref_map.find(tria1)==tria_ref_map.end() && ref!=0)
        {
            tria_ref_map[tria1] = ref;
            int* tria_tmp   = new int[4];
            tria_tmp[0]     = ifn->getVal(i,0+1)-1;
            tria_tmp[1]     = ifn->getVal(i,1+1)-1;
            tria_tmp[2]     = ifn->getVal(i,3+1)-1;
            tria_tmp[3]     = ref;
            tria_ref.push_back(tria_tmp);
        }
        
        tria0.clear();
        tria1.clear();
        
    }
    myfile20.close();
    
    std::ofstream myfile21;
    myfile21.open("tria_ref.dat");
    for(int i=0;i<tria_ref.size();i++)
    {
        myfile21 << tria_ref[i][0] << " " << tria_ref[i][1] << " " << tria_ref[i][2] << " " << tria_ref[i][3]  << std::endl;
    }
    myfile21.close();
    
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
    
    us3d->ien           = ien_copy;
    us3d->ief           = ief_copy;
    us3d->iee           = iee_copy;
    
    us3d->ifn           = ifn_copy;
    us3d->ifn_ref       = ifn_copy;
    us3d->ife           = ife_copy;
    
    us3d->interior      = interior;
    us3d->ghost         = ghost;
    
    us3d->zdefs         = zdefs;
    us3d->bnd_m         = bnd_m;
    us3d->bnd_map       = bnd_map;
    us3d->nBnd          = nBnd;
    us3d->tria_ref_map  = tria_ref_map;
    us3d->tria_ref      = tria_ref;
    
    return us3d;
}
