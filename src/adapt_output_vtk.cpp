#include "adapt_output_vtk.h"
#include "adapt_distri_parstate.h"
#include "adapt_operations.h"


void OutputHexahedralMeshOnRootVTK(MPI_Comm comm,
    string filename, 
    std::set<int> OwnedElem,
    std::map<int,std::vector<int> > gE2lV,
    std::map<int,std::vector<double> > loc_data,
    std::map<int,std::string > varnames,
    std::map<int, std::vector<double> > LocalVerts)
{
    //==================================================================================
    //==================================================================================
    //==================================================================================

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

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

    std::map<int,std::vector<int> > gE2lV_tmp;
    std::set<int>::iterator its;
    for(its=OwnedElem.begin();its!=OwnedElem.end();its++)
    {
        int elid = *its;//OwnedElem[i];
        if(gE2lV.find(elid)!=gE2lV_tmp.end() && gE2lV_tmp.find(elid)==gE2lV_tmp.end())
        {
            gE2lV_tmp[elid] = gE2lV[elid];
        }
    }

    //===============================================================================
    // Pack data;

    int nlocelem = OwnedElem.size();
    std::vector<int> OwnedElem_vec(OwnedElem.begin(), OwnedElem.end());
    
    DistributedParallelState* distElem = new DistributedParallelState(nlocelem,comm);
    int nElem    = distElem->getNel();
    std::vector<int> ElemOnRoot;
    if(world_rank == 0)
    {
        ElemOnRoot.resize(nElem,0);
    }
    MPI_Gatherv(&OwnedElem_vec.data()[0],
    nlocelem,
    MPI_INT,
    &ElemOnRoot.data()[0],
    distElem->getNlocs(),
    distElem->getOffsets(),
    MPI_INT, 0, comm);

    std::map<int,std::vector<int> > gE2lVOnRoot     = GatherGlobalMapOnRoot_T(gE2lV_tmp,comm);
    std::map<int,std::vector<double> > glob_data    = GatherGlobalMapOnRoot_T(loc_data,comm);
    std::map<int,std::vector<double> > VertsOnRoot  = GatherGlobalMapOnRoot_T(LocalVerts,comm);

    //===============================================================================

    if(world_rank == 0)
    {
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
            std::cout << "Warning :: the length of the variable name map does not correspond with the data size. ---> " << filename << " " << loc_data.begin()->second.size() << " " << varnames.size() << std::endl;
        }

        std::map<int,vtkDoubleArray*> mapVars;
        for(itis=varnames.begin();itis!=varnames.end();itis++)
        {
            vtkDoubleArray* VArray = vtkDoubleArray::New();    
            VArray->SetNumberOfComponents(1);
            VArray->SetName(varnames_new[itis->first]);
            varnames_new[itis->first];
            mapVars[itis->first] = VArray;
        }

        // vtkDoubleArray* TArray = vtkDoubleArray::New();    
        // TArray->SetNumberOfComponents(1);
        // TArray->SetName("Temperature");
        // vtkDoubleArray* TKEArray = vtkDoubleArray::New();    
        // TKEArray->SetNumberOfComponents(0);
        // TKEArray->SetName("TKE");

        vtkIdType ielement = 0;
        std::map<int,std::vector<int> >::iterator itmiv;

        //for(int k=0;k<ElemOnRoot.size();k++)
        std::map<int,std::vector<int> >::iterator iter;

        for(iter=gE2lVOnRoot.begin();iter!=gE2lVOnRoot.end();iter++)
        {
            //int elid   = ElemOnRoot[k];
            int elid = iter->first;

            vtkNew<vtkHexahedron> hex;


            for(int q=0;q<gE2lVOnRoot[elid].size();q++)
            {
                int gvid = gE2lVOnRoot[elid][q];
                //hex->GetPointIds()->SetId(q, gvid);

                if(g2lv.find(gvid)==g2lv.end())
                {
                    g2lv[gvid]=lvid;
                    points->InsertNextPoint(VertsOnRoot[gvid][0],
                                    VertsOnRoot[gvid][1],
                                    VertsOnRoot[gvid][2]);


                    hex->GetPointIds()->SetId(q, lvid);

                    lvid++;
                }
                else
                {
                    int lvidn = g2lv[gvid];
                    hex->GetPointIds()->SetId(q, lvidn);
                }
            }

            cellArray->InsertNextCell(hex);

            // TKEArray->InsertNextTuple(&loc_data[elid][0]);
            // TArray->InsertNextTuple(&loc_data[elid][1]);

            for(itis=varnames.begin();itis!=varnames.end();itis++)
            {
                mapVars[itis->first]->InsertNextTuple(&glob_data[elid][itis->first]);
            }

            ielement++;
        }

        for(itis=varnames.begin();itis!=varnames.end();itis++)
        {
            vtkmesh->GetCellData()->AddArray(mapVars[itis->first]);
        }

        vtkmesh->SetPoints(points);
        vtkmesh->SetCells(VTK_HEXAHEDRON, cellArray);
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(filename_char);
        writer->SetHeaderTypeToUInt64();  // Use UInt64 header format
        writer->SetDataModeToBinary();      // Write in ASCII mode
        writer->SetInputData(vtkmesh);
        writer->Write();

    }

/**/


}





void OutputTetraMeshOnRootVTK(MPI_Comm comm,
                                string filename, 
                                std::set<int> OwnedElem,
                                std::map<int,std::vector<int> > gE2lV,
                                std::map<int,std::vector<double> > loc_data,
                                std::map<int,std::string > varnames,
                                std::map<int, std::vector<double> > LocalVerts)
{
    //==================================================================================
    //==================================================================================
    //==================================================================================
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

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

    std::map<int,std::vector<int> > gE2lV_tmp;
    std::set<int>::iterator its;
    for(its=OwnedElem.begin();its!=OwnedElem.end();its++)
    {
        int elid = *its;//OwnedElem[i];
        if(gE2lV.find(elid)!=gE2lV_tmp.end() && gE2lV_tmp.find(elid)==gE2lV_tmp.end())
        {
            gE2lV_tmp[elid] = gE2lV[elid];
        }
    }
    
    //===============================================================================
    // Pack data;

    int nlocelem = OwnedElem.size();
    std::vector<int> OwnedElem_vec(OwnedElem.begin(), OwnedElem.end());
   \
    DistributedParallelState* distElem = new DistributedParallelState(nlocelem,comm);
    int nElem    = distElem->getNel();
    std::vector<int> ElemOnRoot;
    if(world_rank == 0)
    {
        ElemOnRoot.resize(nElem,0);
    }
    MPI_Gatherv(&OwnedElem_vec.data()[0],
            nlocelem,
            MPI_INT,
            &ElemOnRoot.data()[0],
            distElem->getNlocs(),
            distElem->getOffsets(),
            MPI_INT, 0, comm);
    
    std::map<int,std::vector<int> > gE2lVOnRoot     = GatherGlobalMapOnRoot_T(gE2lV_tmp,comm);
    std::map<int,std::vector<double> > glob_data    = GatherGlobalMapOnRoot_T(loc_data,comm);
    std::map<int,std::vector<double> > VertsOnRoot  = GatherGlobalMapOnRoot_T(LocalVerts,comm);
    
    //===============================================================================
    
    if(world_rank == 0)
    {

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
            std::cout << "Warning :: the length of the variable name map does not correspond with the data size. ---> " << filename << " " << loc_data.begin()->second.size() << " " << varnames.size() << std::endl;
        }
        
        std::map<int,vtkDoubleArray*> mapVars;
        for(itis=varnames.begin();itis!=varnames.end();itis++)
        {
            vtkDoubleArray* VArray = vtkDoubleArray::New();    
            VArray->SetNumberOfComponents(1);
            VArray->SetName(varnames_new[itis->first]);
            varnames_new[itis->first];
            mapVars[itis->first] = VArray;
        }

        // vtkDoubleArray* TArray = vtkDoubleArray::New();    
        // TArray->SetNumberOfComponents(1);
        // TArray->SetName("Temperature");
        // vtkDoubleArray* TKEArray = vtkDoubleArray::New();    
        // TKEArray->SetNumberOfComponents(0);
        // TKEArray->SetName("TKE");

        vtkIdType ielement = 0;
        std::map<int,std::vector<int> >::iterator itmiv;
        
        //for(int k=0;k<ElemOnRoot.size();k++)
        std::map<int,std::vector<int> >::iterator iter;
        
        for(iter=gE2lVOnRoot.begin();iter!=gE2lVOnRoot.end();iter++)
        {
            //int elid   = ElemOnRoot[k];
            int elid = iter->first;

            vtkNew<vtkTetra> tetra;
            
            
            for(int q=0;q<gE2lVOnRoot[elid].size();q++)
            {
                
                int gvid = gE2lVOnRoot[elid][q];
                tetra->GetPointIds()->SetId(q, gvid);
                
                if(g2lv.find(gvid)==g2lv.end())
                {
                    g2lv[gvid]=lvid;
                    points->InsertNextPoint(VertsOnRoot[gvid][0],
                                            VertsOnRoot[gvid][1],
                                            VertsOnRoot[gvid][2]);
                    

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
            
            // TKEArray->InsertNextTuple(&loc_data[elid][0]);
            // TArray->InsertNextTuple(&loc_data[elid][1]);

            for(itis=varnames.begin();itis!=varnames.end();itis++)
            {
                mapVars[itis->first]->InsertNextTuple(&glob_data[elid][itis->first]);
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
}






void OutputTetraMeshNoSolutionOnRootVTK(MPI_Comm comm,
                                string filename, 
                                std::set<int> OwnedElem,
                                std::map<int,std::vector<int> > gE2lV,
                                std::map<int, std::vector<double> > LocalVerts)
{
    //==================================================================================
    //==================================================================================
    //==================================================================================
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    const int length = filename.length();
    char* filename_char = new char[length + 1];
    strcpy(filename_char, filename.c_str());

    std::map<int,char *> varnames_new;
    std::map<int,std::string>::iterator itis;
    
    // for(itis=varnames.begin();itis!=varnames.end();itis++)
    // {
    //     const int length = itis->second.length();
    //     char* varname_char = new char[length + 1];
    //     strcpy(varname_char, itis->second.c_str());
    //     varnames_new[itis->first]=varname_char;
    // }

    std::map<int,std::vector<int> > gE2lV_tmp;
    std::set<int>::iterator its;
    for(its=OwnedElem.begin();its!=OwnedElem.end();its++)
    {
        int elid = *its;//OwnedElem[i];
        if(gE2lV.find(elid)!=gE2lV_tmp.end() && gE2lV_tmp.find(elid)==gE2lV_tmp.end())
        {
            gE2lV_tmp[elid] = gE2lV[elid];
        }
    }
    
    //===============================================================================
    // Pack data;

    int nlocelem = OwnedElem.size();
    std::vector<int> OwnedElem_vec(OwnedElem.begin(), OwnedElem.end());
   \
    DistributedParallelState* distElem = new DistributedParallelState(nlocelem,comm);
    int nElem    = distElem->getNel();
    std::vector<int> ElemOnRoot;
    if(world_rank == 0)
    {
        ElemOnRoot.resize(nElem,0);
    }
    MPI_Gatherv(&OwnedElem_vec.data()[0],
            nlocelem,
            MPI_INT,
            &ElemOnRoot.data()[0],
            distElem->getNlocs(),
            distElem->getOffsets(),
            MPI_INT, 0, comm);
    
    std::map<int,std::vector<int> > gE2lVOnRoot     = GatherGlobalMapOnRoot_T(gE2lV_tmp,comm);
    //std::map<int,std::vector<double> > glob_data    = GatherGlobalMapOnRoot_T(loc_data,comm);
    std::map<int,std::vector<double> > VertsOnRoot  = GatherGlobalMapOnRoot_T(LocalVerts,comm);
    
    //===============================================================================
    //std::cout << "gE2lV_tmp " << gE2lV_tmp.size() << std::endl;
    if(world_rank == 0)
    {

        vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkCellArray> cellArray =
                vtkSmartPointer<vtkCellArray>::New();
        vtkNew<vtkPoints> points;

        std::map<int,int> g2lv;
        int lvid = 0;
        int numberOfTetrahedra = gE2lV.size();

        // if(loc_data.begin()->second.size() != varnames.size())
        // {
        //     std::cout << "Warning :: the length of the variable name map does not correspond with the data size. ---> " << filename << " " << loc_data.begin()->second.size() << " " << varnames.size() << std::endl;
        // }
        
        // std::map<int,vtkDoubleArray*> mapVars;
        // for(itis=varnames.begin();itis!=varnames.end();itis++)
        // {
        //     vtkDoubleArray* VArray = vtkDoubleArray::New();    
        //     VArray->SetNumberOfComponents(1);
        //     VArray->SetName(varnames_new[itis->first]);
        //     varnames_new[itis->first];k
        //     mapVars[itis->first] = VArray;
        // }

        // vtkDoubleArray* TArray = vtkDoubleArray::New();    
        // TArray->SetNumberOfComponents(1);
        // TArray->SetName("Temperature");
        // vtkDoubleArray* TKEArray = vtkDoubleArray::New();    
        // TKEArray->SetNumberOfComponents(0);
        // TKEArray->SetName("TKE");

        vtkIdType ielement = 0;
        std::map<int,std::vector<int> >::iterator itmiv;
        
        //for(int k=0;k<ElemOnRoot.size();k++)
        std::map<int,std::vector<int> >::iterator iter;
        std::cout << "gE2lVOnRoot " << gE2lVOnRoot.size() << std::endl;
        for(iter=gE2lVOnRoot.begin();iter!=gE2lVOnRoot.end();iter++)
        {
            //int elid   = ElemOnRoot[k];
            int elid = iter->first;

            vtkNew<vtkTetra> tetra;
            
            
            for(int q=0;q<gE2lVOnRoot[elid].size();q++)
            {
                
                int gvid = gE2lVOnRoot[elid][q];
                tetra->GetPointIds()->SetId(q, gvid);
                
                if(g2lv.find(gvid)==g2lv.end())
                {
                    g2lv[gvid]=lvid;
                    points->InsertNextPoint(VertsOnRoot[gvid][0],
                                            VertsOnRoot[gvid][1],
                                            VertsOnRoot[gvid][2]);
                    

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
            
            // TKEArray->InsertNextTuple(&loc_data[elid][0]);
            // TArray->InsertNextTuple(&loc_data[elid][1]);

            // for(itis=varnames.begin();itis!=varnames.end();itis++)
            // {
            //     mapVars[itis->first]->InsertNextTuple(&glob_data[elid][itis->first]);
            // }
            
            ielement++;
        }
        
        // for(itis=varnames.begin();itis!=varnames.end();itis++)
        // {
        //     vtkmesh->GetCellData()->AddArray(mapVars[itis->first]);
        // }
        
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
}







void OutputTetraMeshPartitionVTK(MPI_Comm comm,
                            string filename, 
                            std::set<int> OwnedElem,
                            std::map<int,std::vector<int> > gE2lV,
                            std::map<int,std::vector<double> > loc_data,
                            std::map<int,std::string > varnames,
                            std::map<int, std::vector<double> > LocalVerts)
{
	//==================================================================================
    //==================================================================================
    //==================================================================================
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

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

    //std::cout << "numberOfTetrahedra " << numberOfTetrahedra << std::endl;
    if(loc_data.begin()->second.size() != varnames.size())
    {
        std::cout << "Warning :: the length of the variable name map does not correspond with the data size ---> " << filename << " " << loc_data.begin()->second.size() << " " << varnames.size() << std::endl;
    }

    std::map<int,vtkDoubleArray*> mapVars;
    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        vtkDoubleArray* VArray = vtkDoubleArray::New();    
        VArray->SetNumberOfComponents(1);
        VArray->SetName(varnames_new[itis->first]);
        varnames_new[itis->first];
        mapVars[itis->first] = VArray;
    }

    // vtkDoubleArray* TArray = vtkDoubleArray::New();    
    // TArray->SetNumberOfComponents(1);
    // TArray->SetName("Temperature");
    // vtkDoubleArray* TKEArray = vtkDoubleArray::New();    
    // TKEArray->SetNumberOfComponents(0);
    // TKEArray->SetName("TKE");

    vtkIdType ielement = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    std::set<int>::iterator its;
    for(its=OwnedElem.begin();its!=OwnedElem.end();its++)
    {
        int elid   = *its;//OwnedElem[k];
        
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
        
        // TKEArray->InsertNextTuple(&loc_data[elid][0]);
        // TArray->InsertNextTuple(&loc_data[elid][1]);

        for(itis=varnames.begin();itis!=varnames.end();itis++)
        {
            mapVars[itis->first]->InsertNextTuple(&loc_data[elid][itis->first]);
        }

        ielement++;
        /**/
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








void OutputTetraMeshNoSolutionPartitionVTK(MPI_Comm comm,
                            string filename, 
                            std::set<int> OwnedElem,
                            std::map<int,std::vector<int> > gE2lV,
                            std::map<int, std::vector<double> > LocalVerts)
{
	//==================================================================================
    //==================================================================================
    //==================================================================================
    
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

    const int length = filename.length();
    char* filename_char = new char[length + 1];
    strcpy(filename_char, filename.c_str());


    // std::map<int,char *> varnames_new;
    // std::map<int,std::string>::iterator itis;

    // for(itis=varnames.begin();itis!=varnames.end();itis++)
    // {
    //     const int length = itis->second.length();
    //     char* varname_char = new char[length + 1];
    //     strcpy(varname_char, itis->second.c_str());
    //     varnames_new[itis->first]=varname_char;
    // }


    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkCellArray> cellArray =
            vtkSmartPointer<vtkCellArray>::New();
    vtkNew<vtkPoints> points;

    std::map<int,int> g2lv;
    int lvid = 0;
    int numberOfTetrahedra = gE2lV.size();

    // if(loc_data.begin()->second.size() != varnames.size())
    // {
    //     std::cout << "Warning :: the length of the variable name map does not correspond with the data size ---> " << filename << " " << loc_data.begin()->second.size() << " " << varnames.size() << std::endl;
    // }

    // std::map<int,vtkDoubleArray*> mapVars;
    // for(itis=varnames.begin();itis!=varnames.end();itis++)
    // {
    //     vtkDoubleArray* VArray = vtkDoubleArray::New();    
    //     VArray->SetNumberOfComponents(1);
    //     VArray->SetName(varnames_new[itis->first]);
    //     varnames_new[itis->first];
    //     mapVars[itis->first] = VArray;
    // }

    // vtkDoubleArray* TArray = vtkDoubleArray::New();    
    // TArray->SetNumberOfComponents(1);
    // TArray->SetName("Temperature");
    // vtkDoubleArray* TKEArray = vtkDoubleArray::New();    
    // TKEArray->SetNumberOfComponents(0);
    // TKEArray->SetName("TKE");

    vtkIdType ielement = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    std::set<int>::iterator its;
    for(its=OwnedElem.begin();its!=OwnedElem.end();its++)
    {
        int elid   = *its;//OwnedElem[k];

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
        
        // TKEArray->InsertNextTuple(&loc_data[elid][0]);
        // TArray->InsertNextTuple(&loc_data[elid][1]);

        // for(itis=varnames.begin();itis!=varnames.end();itis++)
        // {
        //     mapVars[itis->first]->InsertNextTuple(&loc_data[elid][itis->first]);
        // }

        ielement++;
    }

    // for(itis=varnames.begin();itis!=varnames.end();itis++)
    // {
    //     vtkmesh->GetCellData()->AddArray(mapVars[itis->first]);
    // }
    
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





void OutputTriMeshPartitionVTK(MPI_Comm comm,
                            string filename, 
                            FaceSetPointer FaceMap,
                            std::map<int, std::vector<double> > LocalVerts)
{

    FaceSetPointer::iterator itit;
    int world_size, world_rank;
    MPI_Comm_size(comm, &world_size);
    MPI_Comm_rank(comm, &world_rank);

    // Convert filename to char array
    const int length = filename.length();
    char* filename_char = new char[length + 1];
    strcpy(filename_char, filename.c_str());

    // Create VTK data structures
    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkNew<vtkPoints> points;

    // Create data arrays for variables
    // std::map<int,vtkDoubleArray*> mapVars;
    // for(auto& [key, name] : varnames) {
    //     vtkDoubleArray* array = vtkDoubleArray::New();
    //     array->SetNumberOfComponents(1);
    //     array->SetName(name.c_str());
    //     mapVars[key] = array;
    // }

    // Vertex mapping and population
    std::map<int,int> g2lv;
    int lvid = 0;
    
    for(itit=FaceMap.begin();itit!=FaceMap.end();itit++)
    {
        vtkNew<vtkTriangle> triangle;
        std::vector<int> vrts = (*itit)->GetEdgeIDs();

        //std::cout << (*itit)->GetFaceRef() << " reffie " << (*itit)->GetFaceID() << std::endl;
        // Set triangle vertices (3 points instead of 4)
        for(int q = 0; q < 3; q++) {  // Triangles have 3 vertices
            int gvid = vrts[q];
            
            if(g2lv.find(gvid) == g2lv.end()) {
                g2lv[gvid] = lvid;
                points->InsertNextPoint(LocalVerts[gvid][0],
                                        LocalVerts[gvid][1],
                                        LocalVerts[gvid][2]);
                lvid++;
            }
            triangle->GetPointIds()->SetId(q, g2lv[gvid]);
        }
        
        cellArray->InsertNextCell(triangle);
        
        // // Add cell data
        // for(auto& [var_id, array] : mapVars) {
        //     array->InsertNextTuple(&loc_data[elid][var_id]);
        // }
    }

    // Assemble the mesh
    vtkmesh->SetPoints(points);
    vtkmesh->SetCells(VTK_TRIANGLE, cellArray);  // Changed to triangle cell type
    
    // Add data arrays
    // for(auto& [var_id, array] : mapVars) {
    //     vtkmesh->GetCellData()->AddArray(array);
    //     array->Delete();
    // }

    // Write output
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = 
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename_char);
    writer->SetHeaderTypeToUInt64();
    writer->SetDataModeToBinary();
    writer->SetInputData(vtkmesh);
    writer->Write();

    // Cleanup
    delete[] filename_char;
}



void OutputHexahedralMeshPartitionVTK(MPI_Comm comm,
    string filename, 
    std::set<int> OwnedElem,
    std::map<int,std::vector<int> > gE2lV,
    std::map<int,std::vector<double> > loc_data,
    std::map<int,std::string > varnames,
    std::map<int, std::vector<double> > LocalVerts)
{
    //==================================================================================
    //==================================================================================
    //==================================================================================

    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);

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
        std::cout << "Warning :: the length of the variable name map does not correspond with the data size ---> " << filename << " " << loc_data.begin()->second.size() << " " << varnames.size() << std::endl;
    }

    std::map<int,vtkDoubleArray*> mapVars;
    for(itis=varnames.begin();itis!=varnames.end();itis++)
    {
        vtkDoubleArray* VArray = vtkDoubleArray::New();    
        VArray->SetNumberOfComponents(1);
        VArray->SetName(varnames_new[itis->first]);
        varnames_new[itis->first];
        mapVars[itis->first] = VArray;
    }

    // vtkDoubleArray* TArray = vtkDoubleArray::New();    
    // TArray->SetNumberOfComponents(1);
    // TArray->SetName("Temperature");
    // vtkDoubleArray* TKEArray = vtkDoubleArray::New();    
    // TKEArray->SetNumberOfComponents(0);
    // TKEArray->SetName("TKE");

    vtkIdType ielement = 0;
    std::map<int,std::vector<int> >::iterator itmiv;
    std::set<int>::iterator its;
    for(its=OwnedElem.begin();its!=OwnedElem.end();its++)
    {
        int elid   = *its;//OwnedElem[k];

        vtkNew<vtkHexahedron> hex;

        for(int q=0;q<gE2lV[elid].size();q++)
        {
            int gvid = gE2lV[elid][q];
            hex->GetPointIds()->SetId(q, gvid);

            if(g2lv.find(gvid)==g2lv.end())
            {
                g2lv[gvid]=lvid;
                points->InsertNextPoint(LocalVerts[gvid][0],
                                LocalVerts[gvid][1],
                                LocalVerts[gvid][2]);

                hex->GetPointIds()->SetId(q, lvid);

                lvid++;
            }
            else
            {
                int lvidn = g2lv[gvid];
                hex->GetPointIds()->SetId(q, lvidn);
            }

        }
        cellArray->InsertNextCell(hex);

        // TKEArray->InsertNextTuple(&loc_data[elid][0]);
        // TArray->InsertNextTuple(&loc_data[elid][1]);

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
    vtkmesh->SetCells(VTK_HEXAHEDRON, cellArray);
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
        std::cout << "Warning :: the length of the variable name map does not correspond with the data size ---> " << filename << " " << loc_data.begin()->second.size() << " " << varnames.size() << std::endl;
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