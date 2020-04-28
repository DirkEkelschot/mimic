#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


ParVar* CreateParallelData(int N, MPI_Comm comm);

// This function computes the eptr and eind array which are required for most ParMetis APIs.
// This is done based on the global element2node (e2n) array, the number of elements, the communicator
// and the type of elements we are currently using. For now the type of elements is fixed to an integer
// assuming that we use one type of element throughout the whole mesh. However this needs to become
// an int* in order to allow for hybrid meshes.
// e2n has the Nvert per element stored consecutively for each element. Hence this array is Nel*NvertPerElement long.

ParVar_ParMetis* CreateParallelDataParmetis(Array<int>* e2n, MPI_Comm comm, int type);
