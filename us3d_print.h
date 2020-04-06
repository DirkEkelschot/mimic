#include <iostream>
using namespace std;

template <class T>
void printArray(Array<T> A)
{
    int m = A.nrow;
    int n = A.ncol;
    std::cout << " " << std::endl;
    std::cout << "[";
    for(int i=0;i<m;i++)
    {
        std::cout << "[";
        for(int j=0;j<n;j++)
        {
            std::cout << A.data[i*n+j] << ", ";
        }
        
        if (i == m-1)
        {
            std::cout << "]";
        }
        else
        {
            std::cout << "]," << std::endl;
        }
    }
    std::cout << "]" << std::endl;
    std::cout << " " << std::endl;
}




//template <class T>
void PrintUS3dData(US3dData Vec, char tag)
{
    if (tag == 'g')
    {
        std::cout << " ===============Grid Coordinate Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_grd;i++)
        {
            for(int j=0;j<Vec.cols_grd;j++)
            {
                std::cout << Vec.Coordinates[i*Vec.cols_grd+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
    if (tag == 'c')
    {
        std::cout << " ===============Connection Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_conn;i++)
        {
            for(int j=0;j<Vec.cols_conn;j++)
            {
                std::cout << Vec.Connection[i*Vec.cols_conn+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
    if (tag == 'b')
    {
        std::cout << " ===============Boundary Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_bound;i++)
        {
            for(int j=0;j<Vec.cols_bound;j++)
            {
                std::cout << Vec.Boundaries[i*Vec.cols_bound+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
    if (tag == 'i')
    {
        std::cout << " ===============Interior Values ======= " << std::endl;
        std::cout << " ====================================== " << std::endl;
        for(int i=0;i<Vec.rows_interior;i++)
        {
            for(int j=0;j<Vec.cols_interior;j++)
            {
                std::cout << Vec.Interior[i*Vec.cols_interior+j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << " ====================================== " << std::endl;
    }
    

}

