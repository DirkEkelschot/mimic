#include "adapt_output.h"

using namespace std;

void OutputZone(Partition* part, MPI_Comm comm)
{
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    int nloc = part->loc_elem2verts_loc->getNrow();
    int ncol = part->loc_elem2verts_loc->getNcol();
    string filename = "quantity_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"rho\"" << std::endl;
    int nvert =  part->Verts.size();

    myfile <<"ZONE N = " << nvert << ", E = " << nloc << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    std::cout << rank << " number of nodes -> " << nvert << std::endl;
    for(int i=0;i<nvert;i++)
    {
       myfile << part->Verts[i].x << "   " << part->Verts[i].y << "   " << part->Verts[i].z << "   " << part->variable[i] << std::endl;
    }
    

    
    for(int i=0;i<nloc;i++)
    {
       myfile << part->loc_elem2verts_loc->getVal(i,0)+1 << "  " <<
                 part->loc_elem2verts_loc->getVal(i,1)+1 << "  " <<
                 part->loc_elem2verts_loc->getVal(i,2)+1 << "  " <<
                 part->loc_elem2verts_loc->getVal(i,3)+1 << "  " <<
                 part->loc_elem2verts_loc->getVal(i,4)+1 << "  " <<
                 part->loc_elem2verts_loc->getVal(i,5)+1 << "  " <<
                 part->loc_elem2verts_loc->getVal(i,6)+1 << "  " <<
                 part->loc_elem2verts_loc->getVal(i,7)+1 << std::endl;
    }
    
    
    myfile.close();
}

void OutputQuantityPartition(Partition_old* pa, Array<double>* Quan, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nrow = pa->ien->getNrow();
    int ncol = pa->ien->getNcol();
    int nloc = nrow;
    
    int gid;
    int lid;
    Vert V;
    
    //set<int> gid_set;
    std::map<int,Vert> vert_out;
    std::map<int,double> quan_out;
    int v=0;int u=0;int el=0;
    int* l_vert_id = new int[nrow*(ncol-1)];
    map< int, int > gid_set;

    double Q;
    for(int i=0;i<nloc;i++)
    {
        Q = Quan->getVal(i,0);
        for(int j=0;j<ncol-1;j++)
        {
            gid = pa->ien->getVal(i,j+1)-1;
            lid = pa->glob2loc_Vmap[gid];
            
            if ( gid_set.find( gid ) != gid_set.end() )
            {
                l_vert_id[el*(ncol-1)+j]=gid_set[gid];
            }
            else
            {
                l_vert_id[el*(ncol-1)+j]=v;
                
                V.x = pa->Verts->getVal(lid,0);
                V.y = pa->Verts->getVal(lid,1);
                V.z = pa->Verts->getVal(lid,2);
                
                vert_out[u] = V;
                quan_out[u] = Q;
                u++;
            }
            v++;
        }
        el++;
    }
    
    string filename = "quantity_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"dJ\"" << std::endl;
    myfile <<"ZONE N = " << vert_out.size() << ", E = " << nrow << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    
    std::cout << "number of nodes -> " << vert_out.size() << std::endl;
    
    for(int i=0;i<vert_out.size();i++)
    {
       myfile << vert_out[(i)].x << "   " << vert_out[(i)].y << "   " << vert_out[(i)].z << "   " << quan_out[(i)] << std::endl;
    }

    for(int i=0;i<nrow;i++)
    {
       myfile << l_vert_id[i*8+0]+1 << "  " <<
                 l_vert_id[i*8+1]+1 << "  " <<
                 l_vert_id[i*8+2]+1 << "  " <<
                 l_vert_id[i*8+3]+1 << "  " <<
                 l_vert_id[i*8+4]+1 << "  " <<
                 l_vert_id[i*8+5]+1 << "  " <<
                 l_vert_id[i*8+6]+1 << "  " <<
                 l_vert_id[i*8+7]+1 << std::endl;
    }
    
    
    myfile.close();
    delete[] l_vert_id;

}


void OutputPartionVolumes(ParArray<int>* ien, Array<double>* xcn_on_root, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    Partition_old* pv = CollectVerticesPerRank(ien,xcn_on_root,comm);
    
    int nrow = ien->getNrow();
    int ncol = ien->getNcol();
    
    int gid;
    int lid;
    Vert V;
    
    //set<int> gid_set;
    std::map<int,Vert> vert_out;
    int v=0;int u=0;int el=0;
    int* l_vert_id = new int[nrow*(ncol-1)];
    map< int, int > gid_set;
    for(int i=0;i<nrow;i++)
    {
        for(int j=0;j<ncol-1;j++)
        {
            gid = ien->getVal(i,j+1)-1;
            lid = pv->glob2loc_Vmap[gid];
            
            if ( gid_set.find( gid ) != gid_set.end() )
            {
                l_vert_id[el*(ncol-1)+j]=gid_set[gid];
            }
            else
            {
                l_vert_id[el*(ncol-1)+j]=v;
                
                V.x = pv->Verts->getVal(lid,0);
                V.y = pv->Verts->getVal(lid,1);
                V.z = pv->Verts->getVal(lid,2);
                
                vert_out[u]=V;
                
                u++;
            }
            v++;
        }
        el++;
    }
    
    string filename = "volume_per_rank_" + std::to_string(rank) + ".dat";
    ofstream myfile;
    myfile.open(filename);
    myfile << "TITLE=\"volume_part_"  + std::to_string(rank) +  ".tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << vert_out.size() << ", E = " << nrow << ", DATAPACKING = POINT, ZONETYPE = FEBRICK" << std::endl;
    
    std::cout << "number of nodes -> " << vert_out.size() << std::endl;
    
    for(int i=0;i<vert_out.size();i++)
    {
       myfile << vert_out[(i)].x << "   " << vert_out[(i)].y << "   " << vert_out[(i)].z << std::endl;
    }

    for(int i=0;i<nrow;i++)
    {
       myfile << l_vert_id[i*8+0]+1 << "  " <<
                 l_vert_id[i*8+1]+1 << "  " <<
                 l_vert_id[i*8+2]+1 << "  " <<
                 l_vert_id[i*8+3]+1 << "  " <<
                 l_vert_id[i*8+4]+1 << "  " <<
                 l_vert_id[i*8+5]+1 << "  " <<
                 l_vert_id[i*8+6]+1 << "  " <<
                 l_vert_id[i*8+7]+1 << std::endl;
    }
    myfile.close();
    delete[] l_vert_id;
}



void OutputPartitionFaces()
{
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    const char* fn_conn="grids/adept/conn.h5";
    const char* fn_grid="grids/adept/grid.h5";
    ParArray<int>* ief = ReadDataSetFromFileInParallel<int>(fn_conn,"ief",comm,info);

    std::vector<int> partfaces = GetAdjacencyForUS3D_V4(ief, comm);
    
    if(rank == 0)
    {
        Array<double>* xcn = ReadDataSetFromFile<double>(fn_grid,"xcn");
        Array<int>* ifn = ReadDataSetFromFile<int>(fn_grid,"ifn");
        int n_bc_faces = partfaces.size();
        Vert V;
        int* Loc = new int[n_bc_faces*4];
        map< int, int > Loc2GlobBound;
        map< int, Vert> P_verts;
        int cnt = 0;
        int tel = 0;
        int teller = 0;

        for(int j=0;j<partfaces.size();j++)
        {
            int face = partfaces[j];
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(face,k+1);

                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    
                    P_verts[cnt] = V;
                    
                    cnt++;
                }
                
                teller=teller+1;
            }
            tel++;
        }
        
        ofstream myfile;
        
        string filename = "partfaces_" + std::to_string(rank) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"partitionfaces.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << P_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        
        std::cout << "number of boundary nodes -> " << P_verts.size() << std::endl;
        
        for(int i=0;i<P_verts.size();i++)
        {
           myfile << P_verts[(i)].x << "   " << P_verts[(i)].y << "   " << P_verts[(i)].z << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "  "
                  << Loc[i*4+1]+1 << "  "
                  << Loc[i*4+2]+1 << "  "
                  << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
        
        delete[] Loc;
        delete ifn;
        delete xcn;
     }
}



void WriteBoundaryDataInSerial3(Array<double>* xcn)
{
    string filename = "boundary_nodes.dat";
    ofstream myfile;
    myfile.open(filename);
    
    
    for(int i=0;i<3097156;i++)
    {
        myfile << xcn->getVal(i,0) << " " << xcn->getVal(i,1) << " " << xcn->getVal(i,2) << std::endl;
    }
    myfile.close();
}


void WriteBoundaryDataInSerial(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* detJ_verts, double* vol_verts, double* Jnorm_verts)
{
    for(int bc=3;bc<zdefs->getNrow();bc++)
    {
        map< int, int > Loc2GlobBound;
        map< int, Vert> BC_verts;
        map< int, double> J_verts;
        map< int, double> V_verts;
        map< int, double> Jno_verts;
        int n_bc_faces = zdefs->getVal(bc,4)-zdefs->getVal(bc,3);
        int* Loc = new int[n_bc_faces*4];
        
        int b_start = zdefs->getVal(bc,3)-1;
        int b_end   = zdefs->getVal(bc,4)-1;
        int cnt = 0;
        Vert V;
        int tel = 0;
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(j,k+1);
                
                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    BC_verts[cnt] = V;
                    //std::cout << detJ_verts->getVal(val-1,0)*1.0e08 << std::endl;
                    
                    J_verts[cnt]  = detJ_verts[val-1];
                    V_verts[cnt]  = vol_verts[val-1];
                    Jno_verts[cnt]  = Jnorm_verts[val-1];
                    cnt++;
                }
            }
            tel++;
        }
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"dJ\", \"Vol\", \"dJ/Vol\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        for(int i=0;i<BC_verts.size();i++)
        {
           myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << "   " << J_verts[(i)] << "   " << V_verts[(i)] << "   " << Jno_verts[(i)] << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
        
        delete[] Loc;
    }
}



void WriteBoundaryDataInSerial2(Array<double>* xcn, Array<int>* zdefs, Array<int>* ifn, double* vol_verts)
{
    for(int bc=3;bc<zdefs->getNrow();bc++)
    {
        map< int, int > Loc2GlobBound;
        map< int, Vert> BC_verts;
        map< int, double> J_verts;
        map< int, double> V_verts;
        map< int, double> Jno_verts;
        int n_bc_faces = zdefs->getVal(bc,4)-zdefs->getVal(bc,3);
        int* Loc = new int[n_bc_faces*4];
        
        int b_start = zdefs->getVal(bc,3)-1;
        int b_end   = zdefs->getVal(bc,4)-1;
        int cnt = 0;
        Vert V;
        int tel = 0;
        //int* pltJ = new int[n_bc_faces*4];
        int teller = 0;
        for(int j=b_start;j<b_end;j++)
        {
            for(int k=0;k<4;k++)
            {
                int val = ifn->getVal(j,k+1);
                //pltJ[teller] = val;
                //std::cout << val << " ";
                if ( Loc2GlobBound.find( val ) != Loc2GlobBound.end() )
                {
                    Loc[tel*4+k]=Loc2GlobBound[val];
                }
                else
                {
                    Loc2GlobBound[val] = cnt;
                    Loc[tel*4+k]=cnt;
                    V.x = xcn->getVal(val-1,0);
                    V.y = xcn->getVal(val-1,1);
                    V.z = xcn->getVal(val-1,2);
                    BC_verts[cnt] = V;
                    //std::cout << detJ_verts->getVal(val-1,0)*1.0e08 << std::endl;
                    
                    V_verts[cnt]   = vol_verts[val-1];
                    cnt++;
                }
                
                teller=teller+1;
            }
            //std::cout << std::endl;
            tel++;
        }
        //cout << "\nlargest id = " << largest(pltJ,n_bc_faces*4) << std::endl;
        
        ofstream myfile;
        
        string filename = "boundary_" + std::to_string(bc) + ".dat";
        
        myfile.open(filename);
        myfile << "TITLE=\"boundary.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\", \"Vol\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile <<"ZONE N = " << BC_verts.size() << ", E = " << n_bc_faces << ", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL" << std::endl;
        
        std::cout << "number of boundary nodes -> " << BC_verts.size() << std::endl;
        
        for(int i=0;i<BC_verts.size();i++)
        {
           myfile << BC_verts[(i)].x << "   " << BC_verts[(i)].y << "   " << BC_verts[(i)].z << "   " << V_verts[i] << std::endl;
        }
        
        for(int i=0;i<n_bc_faces;i++)
        {
           myfile << Loc[i*4+0]+1 << "    " << Loc[i*4+1]+1 << "   " << Loc[i*4+2]+1 << "  " << Loc[i*4+3]+1 << std::endl;
        }
        myfile.close();
        
        delete[] Loc;
    }
}
