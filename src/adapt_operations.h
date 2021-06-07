#include "adapt.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_datastruct.h"

#ifndef ADAPT_OPERATIONS_H
#define ADAPT_OPERATIONS_H

int FindRank(int* arr, int size, int val);

int FindBoundaryID(int* arr, int size, int val);

std::vector<int> FindDuplicates(std::vector<int> arr);

std::vector<int> FindDuplicatesInParallel(int* arr, int loc_size, int glob_size, MPI_Comm comm);

InteriorPartitionEntity* FindDuplicatesInParallel_Vec(std::vector<int> arr, int arr_size, int glob_size, MPI_Comm comm);

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
