#include "adapt_output_vtk.h"



void OutputTetraMeshPartitionVTK(string filename, 
                            std::vector<int> OwnedElem,
                            std::map<int,std::vector<int> > gE2lV,
                            std::map<int,std::vector<double> > loc_data,
                            std::map<int,std::string > varnames,
                            std::map<int, std::vector<double> > LocalVerts)
{
	//==================================================================================
    //==================================================================================
    //==================================================================================
    
    const int length = filename.length();
    char* filename_char = new char[length + 1];
    strcpy(filename_char, filename.c_str());


    std::map<int,char *> varnames_new;
    std::map<int,std::string>::iterator itis;

    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        const int length = itis->second.length();
        char* varname_char = new char[length + 1];
        strcpy(varname_char, itis->second.c_str());
        varnames_new[itis->first]=varname_char;
    }


    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkCellArray> cellArray =
            vtkSmartPointer<vtkCellArray>::New();
    vtkNew<vtkPoints> points;

    std::map<int,int> g2lv;
    int lvid = 0;
    int numberOfTetrahedra = gE2lV.size();

    if(loc_data.begin()->second.size() != varnames.size())
    {
        std::cout << "Warning :: the length of the variable name map does not correspond with the data size." << std::endl;
    }

    std::map<int,vtkDoubleArray*> mapVars;
    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        vtkDoubleArray* VArray = vtkDoubleArray::New();    
        VArray->SetNumberOfComponents(itis->first);
        VArray->SetName(varnames_new[itis->first]);
        varnames_new[itis->first];

        mapVars[itis->first] = VArray;
    }

    vtkDoubleArray* TArray = vtkDoubleArray::New();    
    TArray->SetNumberOfComponents(1);
    TArray->SetName("Temperature");
    vtkDoubleArray* TKEArray = vtkDoubleArray::New();    
    TKEArray->SetNumberOfComponents(0);
    TKEArray->SetName("TKE");

    vtkIdType ielement = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    // for(itmiv=gE2lV.begin();itmiv!=gE2lV.end();itmiv++)
    for(int k=0;k<OwnedElem.size();k++)
    {
        int elid   = OwnedElem[k];
        vtkNew<vtkTetra> tetra;
        
        for(int q=0;q<gE2lV[elid].size();q++)
        {
            int gvid = gE2lV[elid][q];
            tetra->GetPointIds()->SetId(q, gvid);
            
            if(g2lv.find(gvid)==g2lv.end())
            {
                g2lv[gvid]=lvid;
                points->InsertNextPoint(LocalVerts[gvid][0],
                                        LocalVerts[gvid][1],
                                        LocalVerts[gvid][2]);

                tetra->GetPointIds()->SetId(q, lvid);

                lvid++;
            }
            else
            {
                int lvidn = g2lv[gvid];
                tetra->GetPointIds()->SetId(q, lvidn);

            }

        }
        cellArray->InsertNextCell(tetra);
        
        TKEArray->InsertNextTuple(&loc_data[elid][0]);
        TArray->InsertNextTuple(&loc_data[elid][1]);

        for(itis=varnames.begin();itis!=varnames.end();itis++)
        {
            mapVars[itis->first]->InsertNextTuple(&loc_data[elid][itis->first]);
        }
        ielement++;
    }

    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        vtkmesh->GetCellData()->AddArray(mapVars[itis->first]);
    }

    vtkmesh->SetPoints(points);
    vtkmesh->SetCells(VTK_TETRA, cellArray);
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename_char);
    writer->SetHeaderTypeToUInt64();  // Use UInt64 header format
    writer->SetDataModeToBinary();      // Write in ASCII mode
    writer->SetInputData(vtkmesh);
    writer->Write();
}





void OutputPrismMeshPartitionVTK(string filename, 
                                std::vector<int> OwnedElem,
                                std::map<int,std::vector<int> > gE2lV,
                                std::map<int,std::vector<double> > loc_data,
                                std::map<int,std::string > varnames,
                                std::map<int, std::vector<double> > LocalVerts)
{
    //==================================================================================
    //==================================================================================
    //==================================================================================
    
    const int length = filename.length();
    char* filename_char = new char[length + 1];
    strcpy(filename_char, filename.c_str());


    std::map<int,char *> varnames_new;
    std::map<int,std::string>::iterator itis;

    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        const int length = itis->second.length();
        char* varname_char = new char[length + 1];
        strcpy(varname_char, itis->second.c_str());
        varnames_new[itis->first]=varname_char;
    }


    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkCellArray> cellArray =
            vtkSmartPointer<vtkCellArray>::New();
    vtkNew<vtkPoints> points;

    std::map<int,int> g2lv;
    int llvid = 0;
    int numberOfTetrahedra = OwnedElem.size();

    if(loc_data.begin()->second.size() != varnames.size())
    {
        std::cout << "loc_data.begin()->second.size()  " << loc_data.begin()->second.size()  << std::endl;
        std::cout << "Warning :: the length of the variable name map does not correspond with the data size." << std::endl;
    }

    std::map<int,vtkDoubleArray*> mapVars;
    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        vtkDoubleArray* VArray = vtkDoubleArray::New();    
        VArray->SetNumberOfComponents(itis->first);
        VArray->SetName(varnames_new[itis->first]);
        varnames_new[itis->first];

        mapVars[itis->first] = VArray;
    }




    vtkIdType ielement = 0;
    std::map<int,std::vector<int> >::iterator itmiv;

    
    //for(itmiv=gE2lV.begin();itmiv!=gE2lV.end();itmiv++)
    for(int k=0;k<OwnedElem.size();k++)
    {
        int elid   = OwnedElem[k];
        vtkSmartPointer<vtkHexagonalPrism> prism = vtkSmartPointer<vtkHexagonalPrism>::New();
        
        for(int q=0;q<gE2lV[elid].size();q++)
        {
            int lvid = gE2lV[elid][q];
            
            if(g2lv.find(lvid)==g2lv.end())
            {
                g2lv[lvid]=llvid;

                if(LocalVerts.find(lvid)!=LocalVerts.end())
                {
                    points->InsertNextPoint(LocalVerts[lvid][0],
                                        LocalVerts[lvid][1],
                                        LocalVerts[lvid][2]);
                }
                else
                {
                    std::cout << "mot found " << lvid << " :: " << std::endl;
                }
                

                prism->GetPointIds()->SetId(q, llvid);

                llvid++;
            }
            else
            {
                int lvidn = g2lv[lvid];
                prism->GetPointIds()->SetId(q, lvidn);
            }
        }


        vtkmesh->InsertNextCell(VTK_WEDGE, prism->GetPointIds());


        for(itis=varnames.begin();itis!=varnames.end();itis++)
        {
            //std::cout << "wqwqf " << itis->first << " " << loc_data[elid][itis->first] << std::endl;
            mapVars[itis->first]->InsertNextTuple(&loc_data[elid][itis->first]);
        }

        ielement++;
    }
    
    
    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        vtkmesh->GetCellData()->AddArray(mapVars[itis->first]);
    }

    //vtkmesh->GetCellData()->AddArray(TArray);
    //vtkmesh->GetCellData()->AddArray(TKEArray);
    vtkmesh->SetPoints(points);
    //vtkmesh->SetCells(VTK_TETRA, cellArray);
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename_char);
    writer->SetHeaderTypeToUInt64();  // Use UInt64 header format
    writer->SetDataModeToBinary();      // Write in ASCII mode
    writer->SetInputData(vtkmesh);
    writer->Write();/**/
}