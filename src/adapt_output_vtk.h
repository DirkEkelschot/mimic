#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkHexahedron.h>
#include <vtkHexagonalPrism.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include "adapt.h"
#include "NekFace.h"
#include <string>
using namespace std;


void OutputHexahedralMeshOnRootVTK(MPI_Comm comm,
                                string filename, 
                                std::set<int> OwnedElem,
                                std::map<int,std::vector<int> > gE2lV,
                                std::map<int,std::vector<double> > loc_data,
                                std::map<int,std::string > varnames,
                                std::map<int, std::vector<double> > LocalVerts);


void OutputTetraMeshOnRootVTK(MPI_Comm comm,
                                string filename, 
                                std::set<int> OwnedElem,
                                std::map<int,std::vector<int> > gE2lV,
                                std::map<int,std::vector<double> > loc_data,
                                std::map<int,std::string > varnames,
                                std::map<int, std::vector<double> > LocalVerts);


void OutputHexahedralMeshPartitionVTK(MPI_Comm comm,
                                string filename, 
                                std::set<int> OwnedElem,
                                std::map<int,std::vector<int> > gE2lV,
                                std::map<int,std::vector<double> > loc_data,
                                std::map<int,std::string > varnames,
                                std::map<int, std::vector<double> > LocalVerts);

void OutputTetraMeshPartitionVTK(MPI_Comm comm,
                            string filename, 
                            std::set<int> OwnedElem,
                            std::map<int,std::vector<int> > gE2lV,
                            std::map<int,std::vector<double> > loc_data,
                            std::map<int,std::string > varnames,
                            std::map<int, std::vector<double> > LocalVerts);

void OutputTetraMeshNoSolutionPartitionVTK(MPI_Comm comm,
                            string filename, 
                            std::set<int> OwnedElem,
                            std::map<int,std::vector<int> > gE2lV,
                            std::map<int, std::vector<double> > LocalVerts);

void OutputTetraMeshNoSolutionOnRootVTK(MPI_Comm comm,
                                string filename, 
                                std::set<int> OwnedElem,
                                std::map<int,std::vector<int> > gE2lV,
                                std::map<int, std::vector<double> > LocalVerts);

void OutputTriMeshPartitionVTK(MPI_Comm comm,
                            string filename, 
                            FaceSetPointer FaceMap,
                            std::map<int, std::vector<double> > LocalVerts);


void OutputPrismMeshPartitionVTK(string filename,
							std::vector<int> Owned_Elem, 
							std::map<int,std::vector<int> > gE2lV,
							std::map<int,std::vector<double> > loc_data,
							std::map<int,std::string > varnames,
							std::map<int, std::vector<double> > LocalVerts);