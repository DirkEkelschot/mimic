#include "us3d_ops.h"

using namespace std;
std::vector<int> FindDuplicates(std::vector<int> arr)
{
    int N = arr.size();
    sort(arr.begin(),arr.end());
    std::vector<int> res;
    std::set<int> check;
    for(int i=0;i<N;i++)
    {
        if(arr[i+1]==arr[i])
        {
            if(check.find(arr[i])==check.end())
            {
                check.insert(arr[i]);
                res.push_back(arr[i]);
            }
        }
    }
    
    return res;
}


std::vector<int> FindDuplicatesInParallel(int* arr, int arr_size, int glob_size, MPI_Comm comm)
{
    int size;
    MPI_Comm_size(comm, &size);
    // Get the rank of the process
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    int levels = log2(size);
    
    int N = arr_size;

    int* glob_arr;
    if (rank == 0)
    {
        glob_arr = new int[glob_size];
    }
    
    int* sorted = mergeSort(levels, rank, arr, N, comm, glob_arr);
    
    std::vector<int> res;
    if(rank == 0)
    {
        
        std::set<int> check;
        for(int i=0;i<glob_size;i++)
        {
            if(sorted[i+1]==sorted[i])
            {
                if(check.find(sorted[i])==check.end())
                {
                    check.insert(sorted[i]);
                    res.push_back(sorted[i]);
                }
            }
        }
    }

    return res;
}


int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int* merge(int* a, int* b, int size_a, int size_b, int* merged, int size)
{

    int size_res = size_a+size_b;
    
    int i=0;
    int j=0;
    int k=0;
    
    while(i<=size_a-1 && j<=size_b-1)
    {
        if(a[i] <= b[j])
        {
            merged[k++] = a[i++];
        }
        else
        {
            merged[k++] = b[j++];
        }
    }
    while(i <= size_a-1)
    {
        merged[k++] = a[i++];
    }
    while(j <= size_b-1)
    {
        merged[k++] = b[j++];
    }
    
    return merged;
}


std::vector<int> merge_vec(std::vector<int> a, std::vector<int> b)
{
    int n = a.size();
    int m = b.size();
    int size = n+m;
    
    std::vector<int> merged(size);
    int i=0;
    int j=0;
    int k=0;
    
    while(i<=n-1 && j<=m-1)
    {
        if(a[i] <= b[j])
        {
            merged[k++] = a[i++];
        }
        else
        {
            merged[k++] = b[j++];
        }
    }
    while(i <= n-1)
    {
        merged[k++] = a[i++];
    }
    while(j <= m-1)
    {
        merged[k++] = b[j++];
    }
    
    return merged;
}


int* mergeSort(int height, int id, int* localArray, int size, MPI_Comm comm, int* globalArray)
{
    int parent, rightChild, myHeight;
    int *half1, *half2, *mergeResult;

    myHeight = 0;
    qsort(localArray, size, sizeof(int), compare); // sort local array
    half1 = localArray;  // assign half1 to localArray
    int size_half1, size_half2;
    
    while (myHeight < height) { // not yet at top
        parent = (id & (~(1 << myHeight)));

        if (parent == id) { // left child
            rightChild = (id | (1 << myHeight));

              // allocate memory and receive array of right child
              
              MPI_Recv(&size_half2, 1, MPI_INT, rightChild, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
              half2 = new int[size_half2];
            
              MPI_Recv(half2, size_half2, MPI_INT, rightChild, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
              // allocate memory for result of merge
              
              // merge half1 and half2 into mergeResult
              int size_half1 = size;
//              if (size_half1!=size_half2)
//              {
//                  std::cout << size_half2<< " " << size_half1 << std::endl;
//              }
              
              //mergeResult = (int*) malloc (size_half2+size_half1);
              mergeResult = new int[size_half2+size_half1];
              mergeResult = merge(half1, half2, size_half1, size_half2, mergeResult, size);
              // reassign half1 to merge result
            half1 = mergeResult;
            size = size_half1+size_half2;  // double size
            
            delete[] half2;
            mergeResult = NULL;

            myHeight++;

        }
        else
        {
            // right child
            // send local array to parent
            MPI_Send(&size,    1, MPI_INT, parent, 0, MPI_COMM_WORLD);
            MPI_Send(half1, size, MPI_INT, parent, 0, MPI_COMM_WORLD);
            if(myHeight != 0)
            {
                delete[] half1;
            }
            myHeight = height;
        }
    }

    if(id == 0){
        globalArray = half1;   // reassign globalArray to half1
    }
    return globalArray;
}



std::vector<int> mergeSort_vec(int height, int id, std::vector<int> localArray, int size, MPI_Comm comm, std::vector<int> globalArray){
    
    int parent, rightChild, myHeight;

    myHeight = 0;
    sort(localArray.begin(),localArray.end());
    //qsort(localArray, size, sizeof(int), compare); // sort local array
    std::vector<int> half1 = localArray;  // assign half1 to localArray
    int size_half1, size_half2;
    while (myHeight < height) { // not yet at top
        parent = (id & (~(1 << myHeight)));

        if (parent == id) { // left child
            rightChild = (id | (1 << myHeight));

              // allocate memory and receive array of right child
              //half2 = (int*) malloc (size * sizeof(int));
            MPI_Recv(&size_half2, 1, MPI_INT, rightChild, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::vector<int> half2(size_half2);
            MPI_Recv(&half2[0], size_half2, MPI_INT, rightChild, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // allocate memory for result of merge
            //mergeResult = (int*) malloc (size * 2 * sizeof(int));
              
            // merge half1 and half2 into mergeResult
            std::vector<int> mergeResult = merge_vec(half1, half2);
            // reassign half1 to merge result
            half1 = mergeResult;
            size = size * 2;  // double size
            half2.erase(half2.begin(),half2.end());
            //free(half2);
            mergeResult.erase(mergeResult.begin(),mergeResult.end());

            myHeight++;

        } else { // right child
              // send local array to parent
            MPI_Send(&size,    1, MPI_INT, parent, 0, MPI_COMM_WORLD);

            MPI_Send(&half1[0], size, MPI_INT, parent, 0, MPI_COMM_WORLD);
            if(myHeight != 0)
            {
                half1.erase(half1.begin(),half1.end());
            }
            myHeight = height;
        }
    }

    if(id == 0){
        globalArray = half1;   // reassign globalArray to half1
    }
    return globalArray;
}
