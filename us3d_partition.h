#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


ParVar* CreateParallelData(int N, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int nloc             = int(N/size) + ( rank < N%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(N/size) + MIN(rank, N%size);
    
    
    int* proc_nlocs                 = new int[size];
    int* proc_offset                = new int[size];
    int* red_proc_nlocs             = new int[size];
    int* red_proc_offset            = new int[size];
    
    for(int i=0;i<size;i++)
    {
        red_proc_nlocs[i] = 0;
        red_proc_offset[i] = 0;
        
        if(i==rank)
        {
            proc_nlocs[i]  = nloc;
            proc_offset[i] = offset;
        }
        else
        {
            proc_nlocs[i]  = 0;
            proc_offset[i] = 0;
        }
    }
    
    MPI_Allreduce(proc_nlocs,  red_proc_nlocs,  size, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(proc_offset, red_proc_offset, size, MPI_INT, MPI_SUM, comm);
    
    ParVar* pv = new ParVar;
    
    pv->size    = size;
    pv->nlocs   = red_proc_nlocs;
    pv->offsets = red_proc_offset;
    
    
    return pv;
}



// This function computes the eptr and eind array which are required for most ParMetis APIs.
// This is done based on the global element2node (e2n) array, the number of elements, the communicator
// and the type of elements we are currently using. For now the type of elements is fixed to an integer
// assuming that we use one type of element throughout the whole mesh. However this needs to become
// an int* in order to allow for hybrid meshes.
// e2n has the Nvert per element stored consecutively for each element. Hence this array is Nel*NvertPerElement long.

ParVar_ParMetis* CreateParallelDataParmetis(Array<int>* e2n, MPI_Comm comm, int type)
{
    int size;
    MPI_Comm_size(comm, &size);

    // Get the rank of the process;
    int rank;
    MPI_Comm_rank(comm, &rank);
    int Nel = e2n->nloc;
    //std::cout << "number of elements = " << Nel;
    int nloc             = int(Nel/size) + ( rank < Nel%size );
    //  compute offset of rows for each proc;
    int offset           = rank*int(Nel/size) + MIN(rank, Nel%size);
    
    int npo_loc = 0;
    for(int i=0;i<nloc;i++)
    {
        npo_loc += type;
    }
    
    int* nlocs                 = new int[size];
    int* red_nlocs             = new int[size];
    int* npo_locs              = new int[size];
    int* red_npo_locs          = new int[size];

    for(int i=0;i<size;i++)
    {
        red_nlocs[i]        = 0;
        red_npo_locs[i]     = 0;
        
        if(i==rank)
        {
            nlocs[i]        = nloc;
            npo_locs[i]     = npo_loc;
        }
        else
        {
            nlocs[i]        = 0;
            npo_locs[i]     = 0;
        }
    }
    
    int* elm_dist              = new int[size+1];
    int* npo_offset            = new int[size+1];
    int* red_elm_dist          = new int[size+1];
    int* red_npo_offset        = new int[size+1];
    
    for(int i=0;i<size+1;i++)
    {
        red_elm_dist[i]   = 0;
        red_npo_offset[i] = 0;
        if(i==rank)
        {
            elm_dist[i]   = offset;
            npo_offset[i] = offset*8;
        }
        else
        {
            elm_dist[i]  = 0;
            npo_offset[i] = 0;
        }
    }
    
    
    MPI_Allreduce(nlocs,        red_nlocs,      size,   MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_locs,     red_npo_locs,   size,   MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(elm_dist,     red_elm_dist,   size+1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(npo_offset,   red_npo_offset, size+1,   MPI_INT, MPI_SUM, comm);

    red_elm_dist[size] = Nel;
    red_npo_offset[size] = Nel*8;

    
    int* eptr = new int[nloc+1];
    int* eind = new int[npo_loc];
    eptr[0]  = 0;
    for(int i=0;i<nloc;i++)
    {
        eptr[i+1]  = eptr[i]+type;
        
        for(int j=eptr[i];j<eptr[i+1];j++)
        {
            eind[j] = e2n->data[npo_offset[rank]+j];
        }
    }
    
    if (rank == 0)
    {
        for(int i=0;i<size+1;i++)
        {
            std::cout << red_elm_dist[i] << " " << red_npo_offset[i] << std::endl;
        }
    }
    //std::cout << rank << " npo " << npo_loc << std::endl;
    
    ParVar_ParMetis* pv_parmetis = new ParVar_ParMetis;
    
    pv_parmetis->size        =  size;
    pv_parmetis->nlocs       =  red_nlocs;
    pv_parmetis->elmdist     =  red_elm_dist;
    pv_parmetis->npo_locs    =  red_npo_locs;
    pv_parmetis->npo_offset  =  red_npo_offset;
    pv_parmetis->eptr        =  eptr;
    pv_parmetis->eind        =  eind;
    
    return pv_parmetis;
}
