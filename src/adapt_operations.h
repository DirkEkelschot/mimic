#include "adapt.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_datastruct.h"
#include "adapt_distri_parstate.h"

#ifndef ADAPT_OPERATIONS_H
#define ADAPT_OPERATIONS_H


template <typename T>
std::map<int,std::vector<T> > GatherJaggedGlobalMapOnRoot_T(std::map<int,std::vector<T> > mappie, MPI_Comm mpi_comm)
{
    std::map<int,std::vector<T> > mappie_glob;
    int world_size;
    MPI_Comm_size(mpi_comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(mpi_comm, &world_rank);

    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    if constexpr (std::is_same_v<T, int>) 
    {
        mpi_type = MPI_INT;
    }
    if constexpr (std::is_same_v<T, double>) 
    {
        mpi_type = MPI_DOUBLE;
    }

    int mapSizeLoc                      = mappie.size();

    typename std::map<int,std::vector<T> >::iterator itt;
    int mapvec_size = 0;
    for(itt=mappie.begin();itt!=mappie.end();itt++)
    {
        mapvec_size = mapvec_size + itt->second.size();
    }

    // int rowsize                          = mappie.begin()->second.size();
    DistributedParallelState* distrimap     = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot                          = distrimap->getNel();
    DistributedParallelState* distrEntrymap = new DistributedParallelState(mapvec_size,mpi_comm);
    int mapEntrySizeTot                     = distrEntrymap->getNel();

    std::vector<int> key_loc(mapSizeLoc,0);
    std::vector<int> size_loc(mapSizeLoc,0);
    std::vector<T> val_loc(mapvec_size,0);
    std::vector<int> key_tot;
    std::vector<int> size_tot;
    std::vector<T> val_tot;

    if(world_rank == 0)
    {
        key_tot.resize(mapSizeTot);
        size_tot.resize(mapSizeTot);
        val_tot.resize(mapEntrySizeTot);
    }

    int i = 0;
    
    typename std::map<int,std::vector<T> >::iterator itred;
    int offset = 0;
    for(itred=mappie.begin();itred!=mappie.end();itred++)
    {
        key_loc[i]  = itred->first;
        size_loc[i] = itred->second.size();
        int rowsize = itred->second.size();

        for(int q=0;q<rowsize;q++)
        {
            val_loc[offset+q] = itred->second[q];
        }
        
        offset = offset + rowsize;
        i++;
    }

    DistributedParallelState* distMapKey = new DistributedParallelState(mapSizeLoc,mpi_comm);
    DistributedParallelState* distMapVal = new DistributedParallelState(mapvec_size,mpi_comm);
    
    MPI_Gatherv(&key_loc.data()[0],
            mapSizeLoc,
            MPI_INT,
            &key_tot.data()[0],
            distMapKey->getNlocs(),
            distMapKey->getOffsets(),
            MPI_INT, 0, mpi_comm);
    
    MPI_Gatherv(&size_loc.data()[0],
            mapSizeLoc,
            MPI_INT,
            &size_tot.data()[0],
            distMapKey->getNlocs(),
            distMapKey->getOffsets(),
            MPI_INT, 0, mpi_comm);
    
    MPI_Gatherv(&val_loc.data()[0],
            val_loc.size(),
            mpi_type,
            &val_tot.data()[0],
            distMapVal->getNlocs(),
            distMapVal->getOffsets(),
            mpi_type, 0, mpi_comm);


    if(world_rank == 0)
    {
        int duple = 0;
        int key;
        T val;
        int offset = 0;
        for(int i=0;i<mapSizeTot;i++)
        {
            key         = key_tot[i];
            int rowsize = size_tot[i];
            std::vector<T> valrow(rowsize);
            for(int q=0;q<rowsize;q++)
            {
                valrow[q] = val_tot[offset+q];
            }
            
            if(mappie_glob.find(key)==mappie_glob.end())
            {
                mappie_glob[key] = valrow;
            }
            else
            {
                duple++;
            }
            offset = offset + rowsize;
        }

        std::cout << "duplicates " << duple << std::endl;
    }

    /**/
    return mappie_glob;
}




template <typename T>
std::map<int,std::vector<T> > GatherGlobalMapOnRoot_T(std::map<int,std::vector<T> > mappie, MPI_Comm mpi_comm)
{
    int world_size;
    MPI_Comm_size(mpi_comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(mpi_comm, &world_rank);

    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    if constexpr (std::is_same_v<T, int>) 
    {
        mpi_type = MPI_INT;
    }
    if constexpr (std::is_same_v<T, double>) 
    {
        mpi_type = MPI_DOUBLE;
    }

    int mapSizeLoc = mappie.size();
    int rowsize = mappie.begin()->second.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot = distrimap->getNel();

    std::vector<int> key_loc(mapSizeLoc,0);
    std::vector<T> val_loc(mapSizeLoc*rowsize,0);
    std::vector<int> key_tot;
    std::vector<T> val_tot;
    if(world_rank == 0)
    {
        key_tot.resize(mapSizeTot);
        val_tot.resize(mapSizeTot*rowsize);
    }

    int i = 0;
    
    typename std::map<int,std::vector<T> >::iterator itred;
    for(itred=mappie.begin();itred!=mappie.end();itred++)
    {
        key_loc[i] = itred->first;

        for(int q=0;q<rowsize;q++)
        {
            val_loc[i*rowsize+q] = itred->second[q];
        }
        i++;
    }

    DistributedParallelState* distMapKey = new DistributedParallelState(mapSizeLoc,mpi_comm);
    DistributedParallelState* distMapVal = new DistributedParallelState(mapSizeLoc*rowsize,mpi_comm);
    MPI_Gatherv(&key_loc.data()[0],
            mapSizeLoc,
            MPI_INT,
            &key_tot.data()[0],
            distMapKey->getNlocs(),
            distMapKey->getOffsets(),
            MPI_INT, 0, mpi_comm);

    MPI_Gatherv(&val_loc.data()[0],
            val_loc.size(),
            mpi_type,
            &val_tot.data()[0],
            distMapVal->getNlocs(),
            distMapVal->getOffsets(),
            mpi_type, 0, mpi_comm);

    std::map<int,std::vector<T> > mappie_glob;

    if(world_rank == 0)
    {
        int key;
        T val;
        for(int i=0;i<mapSizeTot;i++)
        {
            key = key_tot[i];
            std::vector<T> valrow(rowsize);
            for(int q=0;q<rowsize;q++)
            {
                valrow[q] = val_tot[i*rowsize+q];
            }
            
            if(mappie_glob.find(key)==mappie_glob.end())
            {
                mappie_glob[key] = valrow;
            }
        }
    }


    return mappie_glob;
}




template <typename T>
std::vector<std::vector<T> > GatherJaggedGlobalVectorOnRoot_T(std::vector<std::vector<T> > mappie, MPI_Comm mpi_comm)
{
    std::vector<std::vector<T> > mappie_glob;
    int world_size;
    MPI_Comm_size(mpi_comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(mpi_comm, &world_rank);

    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    if constexpr (std::is_same_v<T, int>) 
    {
        mpi_type = MPI_INT;
    }
    if constexpr (std::is_same_v<T, double>) 
    {
        mpi_type = MPI_DOUBLE;
    }

    int mapSizeLoc                      = mappie.size();

    int mapvec_size = 0;
    for(int q=0;q<mappie.size();q++)
    {
        mapvec_size = mapvec_size + mappie[q].size();
    }

    // int rowsize                         = mappie.begin()->second.size();
    DistributedParallelState* distrimap     = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot                          = distrimap->getNel();
    DistributedParallelState* distrEntrymap = new DistributedParallelState(mapvec_size,mpi_comm);
    int mapEntrySizeTot                     = distrEntrymap->getNel();

    std::vector<int> size_loc(mapSizeLoc,0);
    std::vector<T> val_loc(mapvec_size,0);
    std::vector<int> size_tot;
    std::vector<T> val_tot;

    if(world_rank == 0)
    {
        size_tot.resize(mapSizeTot);
        val_tot.resize(mapEntrySizeTot);
    }

    int i = 0;
    
    typename std::vector<std::vector<T> >::iterator itred;
    int offset = 0;
    for(int q=0;q<mappie.size();q++)
    {
        
        size_loc[i] = mappie[q].size();
        int rowsize = mappie[q].size();

        for(int p=0;p<rowsize;p++)
        {
            val_loc[offset+p] = mappie[q][p];
        }
        
        offset = offset + rowsize;
        i++;
    }

    DistributedParallelState* distMapKey = new DistributedParallelState(mapSizeLoc,mpi_comm);
    DistributedParallelState* distMapVal = new DistributedParallelState(mapvec_size,mpi_comm);
    
    // MPI_Gatherv(&key_loc.data()[0],
    //         mapSizeLoc,
    //         MPI_INT,
    //         &key_tot.data()[0],
    //         distMapKey->getNlocs(),
    //         distMapKey->getOffsets(),
    //         MPI_INT, 0, mpi_comm);
    
    MPI_Gatherv(&size_loc.data()[0],
            mapSizeLoc,
            MPI_INT,
            &size_tot.data()[0],
            distMapKey->getNlocs(),
            distMapKey->getOffsets(),
            MPI_INT, 0, mpi_comm);
    
    MPI_Gatherv(&val_loc.data()[0],
            val_loc.size(),
            mpi_type,
            &val_tot.data()[0],
            distMapVal->getNlocs(),
            distMapVal->getOffsets(),
            mpi_type, 0, mpi_comm);


    if(world_rank == 0)
    {
        int duple = 0;
        int key;
        T val;
        int offset = 0;
        for(int i=0;i<mapSizeTot;i++)
        {
            int rowsize = size_tot[i];
            std::vector<T> valrow(rowsize);
            for(int q=0;q<rowsize;q++)
            {
                valrow[q] = val_tot[offset+q];
            }
            
            mappie_glob.push_back(valrow);
            // if(mappie_glob.find(key)==mappie_glob.end())
            // {
            //     mappie_glob[key] = valrow;
            // }
            // else
            // {
            //     duple++;
            // }
            offset = offset + rowsize;
        }

        std::cout << "duplicates " << duple << std::endl;
    }

    /**/
    return mappie_glob;
}




template <typename T>
std::map<int,T> AllGatherMap_T(std::map<int,T> mappie, MPI_Comm mpi_comm)
{
    int world_size;
    MPI_Comm_size(mpi_comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(mpi_comm, &world_rank);

    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    if constexpr (std::is_same_v<T, int>) 
    {
        mpi_type = MPI_INT;
    }
    if constexpr (std::is_same_v<T, double>) 
    {
        mpi_type = MPI_DOUBLE;
    }

    int mapSizeLoc = mappie.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot = distrimap->getNel();
    
    int* key_loc    = new int[mapSizeLoc];
    T* val_loc      = new T[mapSizeLoc];
    int* key_tot    = new int[mapSizeTot];
    T* val_tot      = new T[mapSizeTot];
    int i = 0;
    
    typename std::map<int,T>::iterator itred;
    for(itred=mappie.begin();itred!=mappie.end();itred++)
    {
        key_loc[i] = itred->first;
        val_loc[i] = itred->second;
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    MPI_Allgatherv(key_loc,
                   mapSizeLoc,
                   MPI_INT,
                   key_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    
    MPI_Allgatherv(val_loc,
                   mapSizeLoc,
                   mpi_type,
                   val_tot,
                   nlocs,
                   offsets,
                   mpi_type, mpi_comm);
    
    int key;
    T val;
    std::map<int,T> mappie_glob;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = key_tot[i];
        val = val_tot[i];
        if(mappie_glob.find(key)==mappie_glob.end())
        {
        	mappie_glob[key] = val;
        }
    }
    
    delete[] key_loc;
    delete[] val_loc;
    delete[] key_tot;
    delete[] val_tot;

    return mappie_glob;
}


template <typename T>
std::set<T> AllGatherSet(std::set<T> set_tmp, MPI_Comm mpi_comm)
{
    int world_size;
    MPI_Comm_size(mpi_comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(mpi_comm, &world_rank);

    MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
    if constexpr (std::is_same_v<T, int>) 
    {
        mpi_type = MPI_INT;
    }
    if constexpr (std::is_same_v<T, double>) 
    {
        mpi_type = MPI_DOUBLE;
    }

    int mapSizeLoc = set_tmp.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot = distrimap->getNel();
    
    // T* key_loc    = new int[mapSizeLoc];
    // T* key_tot    = new int[mapSizeTot];
    int i = 0;
    
    std::vector<T> key_loc(mapSizeLoc);
    std::vector<T> key_tot(mapSizeTot);

    typename std::set<T>::iterator itred;
    for(itred=set_tmp.begin();itred!=set_tmp.end();itred++)
    {
        key_loc[i] = *itred;
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    MPI_Allgatherv(&key_loc.data()[0],
                   mapSizeLoc,
                   mpi_type,
                   &key_tot.data()[0],
                   nlocs,
                   offsets,
                   mpi_type, mpi_comm);
    
    T key;
    std::set<T> set_glob;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = key_tot[i];
        if(set_glob.find(key)==set_glob.end())
        {
        	set_glob.insert(key);
        }
    }
    
    key_loc.clear();
    key_tot.clear();

    return set_glob;
}

std::map<int,int> AllGatherMap_I(std::map<int,int> mappie, MPI_Comm mpi_comm);
std::map<int,double> AllGatherMap_D(std::map<int,double> mappie, MPI_Comm mpi_comm);

int FindRank(int* arr, int size, int val);

int FindBoundaryID(int* arr, int size, int val);

std::vector<int> FindDuplicates(std::vector<int> arr);

std::vector<int> FindDuplicatesInParallel(int* arr, int loc_size, int glob_size, MPI_Comm comm);

//InteriorPartitionEntity* FindDuplicatesInParallel_Vec(std::vector<int> arr, int arr_size, int glob_size, MPI_Comm comm);

std::vector<int> FindDuplicatesInParallel_VecV2(std::vector<int> arr, int arr_size, int glob_size, MPI_Comm comm);

int compare (const void * a, const void * b);

int* merge(int* a, int* b, int* merged, int size);

std::vector<int> merge_vec(std::vector<int> a, std::vector<int> b);

int* mergeSort(int height, int id, int* localArray, int size, MPI_Comm comm, int* globalArray);

std::vector<int> mergeSort_vec(int height, int rank, std::vector<int> localArray, int size, MPI_Comm comm, std::vector<int> globalArray);

int binarySearch(int* arr, int low, int high, int key);

int largest(int arr[], int n);

void TestFindRank(MPI_Comm comm);

void mergeNew(int *, int *, int, int, int);

void mergeSortNew(int *, int *, int, int);



#endif
