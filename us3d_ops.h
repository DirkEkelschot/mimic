#include "adapt.h"

std::vector<int> FindDuplicates(std::vector<int> arr);

int compare (const void * a, const void * b);

int* merge(int* a, int* b, int* merged, int size);

std::vector<int> merge_vec(std::vector<int> a, std::vector<int> b);

int* mergeSort(int height, int id, int* localArray, int size, MPI_Comm comm, int* globalArray);

std::vector<int> mergeSort_vec(int height, int id, std::vector<int> localArray, int size, MPI_Comm comm, std::vector<int> globalArray);
