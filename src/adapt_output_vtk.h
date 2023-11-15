#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkTetra.h>
#include <vtkHexagonalPrism.h>
#include <vtkCellArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include "adapt.h"
#include <string>
using namespace std;

void OutputMeshPartitionVTK(string filename, std::map<int,std::vector<int> > gE2lV,
							std::map<int,std::vector<double> > loc_data,
							std::map<int,std::string > varnames,
							std::vector<std::vector<double> > LocalVerts);


void OutputPrismMeshPartitionVTK(string filename, std::map<int,std::vector<int> > gE2lV,
							std::map<int,std::vector<double> > loc_data,
							std::map<int,std::string > varnames,
							std::map<int, std::vector<double> > LocalVerts);