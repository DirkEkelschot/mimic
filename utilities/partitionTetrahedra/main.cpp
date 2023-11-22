#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
#include "../../src/NekFace.h"
#include "../../src/adapt_prismtetratrace.h"
#include "../../src/adapt_repartition.h"
#include "../../src/adapt_output_vtk.h"
#include "../../src/adapt_meshtopology_lite.h"
#include "../../src/adapt_gradreconstruct_lite.h"




#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))


void ParseEquals(const std::string &line, std::string &lhs,
                                std::string &rhs)
{
    /// Pull out lhs and rhs and eliminate any spaces.
    size_t beg = line.find_first_not_of(" ");
    size_t end = line.find_first_of("=");
    // Check for no parameter name
    if (beg == end)
        throw 1;
    // Check for no parameter value
    if (end != line.find_last_of("="))
        throw 1;
    // Check for no equals sign
    if (end == std::string::npos)
        throw 1;

    lhs = line.substr(line.find_first_not_of(" "), end - beg);
    lhs = lhs.substr(0, lhs.find_last_not_of(" ") + 1);
    rhs = line.substr(line.find_last_of("=") + 1);
    rhs = rhs.substr(rhs.find_first_not_of(" "));
    rhs = rhs.substr(0, rhs.find_last_not_of(" ") + 1);
}


struct Inputs{
    double hgrad;
    double hmin;
    double hmax;
    double MetScale;
    double hausd;
    int ReadFromStats;
    int RunWakRefinement;
    double hwake;
    int niter;
    int recursive;
    int extended;
    int StateVar;
};


Inputs* ReadXmlFile(const char* filename)
{
    TiXmlDocument *m_xmlDoc = new TiXmlDocument;
    TiXmlDocument doc( filename );
    Inputs* inp = new Inputs;
    doc.LoadFile();
    
    TiXmlHandle hDoc(&doc);
    
//    TiXmlHandle docHandle(m_xmlDoc);
    
//    TiXmlElement *e;
//
//    e = doc->FirstChildElement("METRIC").Element();
//
//    TiXmlElement *parametersElement =
//        conditions->FirstChildElement("PARAMETERS");
    
    TiXmlElement *xmlMetric = doc.FirstChildElement("MIMIC");
    
    
    TiXmlElement *xmlParam = xmlMetric->FirstChildElement("PARAMETERS");
    
    std::map<std::string,double> param_map;
    if (xmlParam)
    {
        TiXmlElement *parameter = xmlParam->FirstChildElement("P");
        
        while (parameter)
        {
            TiXmlNode *node = parameter->FirstChild();
            
            std::string line = node->ToText()->Value(), lhs, rhs;
            
            try
            {
                ParseEquals(line, lhs, rhs);
            }
            catch (...)
            {
                std::cout << "Error reading metric.xml " << std::endl;
            }
            
            if (!lhs.empty() && !rhs.empty())
            {
                double value = std::stod(rhs);
                param_map[lhs] = value;
                
            }
            parameter = parameter->NextSiblingElement();
        }
    }
    
    if(param_map.find("hMinimum")!=param_map.end())
    {
        inp->hmin = param_map["hMinimum"];
    }
    else
    {
        std::cout << "Error: hMinimum is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("hMaximum")!=param_map.end())
    {
        inp->hmax = param_map["hMaximum"];
    }
    else
    {
        std::cout << "Error: hMaximum is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("hGradation")!=param_map.end())
    {
        inp->hgrad = param_map["hGradation"];
    }
    else
    {
        std::cout << "Error: hGradation is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("Scaling")!=param_map.end())
    {
        inp->MetScale = param_map["Scaling"];
    }
    else
    {
        std::cout << "Error: Scaling is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("HausDorff")!=param_map.end())
    {
        inp->hausd = param_map["HausDorff"];
    }
    else
    {
        std::cout << "Error: HausDorff is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("nIterations")!=param_map.end())
    {
        inp->niter = param_map["nIterations"];
    }
    else
    {
        std::cout << "Error: nIterations is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("RecursiveReconstruction")!=param_map.end())
    {
        inp->recursive = param_map["RecursiveReconstruction"];
    }
    else
    {
        std::cout << "Error: RecursiveReconstruction is not defined in metric.xml." << std::endl;
    }
    
    if(param_map.find("ExtendedScheme")!=param_map.end())
    {
        inp->extended = param_map["ExtendedScheme"];
    }
    else
    {
        std::cout << "Error: RecursiveReconstruction is not defined in metric.xml." << std::endl;
    }
    
    if(param_map.find("UseStatistics")!=param_map.end())
    {
        inp->ReadFromStats = param_map["UseStatistics"];
    }
    else
    {
        std::cout << "Error: UseStatistics is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("WakeRefinement")!=param_map.end())
    {
        inp->RunWakRefinement = param_map["WakeRefinement"];
    }
    else
    {
        std::cout << "Error: WakeRefinement is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("hWake")!=param_map.end())
    {
        inp->hwake = param_map["hWake"];
    }
    else
    {
        std::cout << "Error: hWake is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("StateVariable")!=param_map.end())
    {
        inp->StateVar = param_map["StateVariable"];
    }
    else
    {
        std::cout << "Error: StateVariable is not defined in metric.xml." << std::endl;
    }
    
    
    return inp;
}



int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);
    FILE            *inm;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int i,j,k;
    clock_t t0_met = clock();

    int ier,opt;
    int debug = 1;
    
    
    const char* fn_grid="inputs/grid.h5";
    const char* fn_conn="inputs/conn.h5";
    const char* fn_data="inputs/data.h5";

    Inputs* inputs = ReadXmlFile("inputs/metric.xml");

    
    
    
    
//    TiXmlElement *parametersElement =
//        conditions->FirstChildElement("PARAMETERS");
    

    
//    TiXmlHandle docHandle(m_xmlDoc);
//    TiXmlElement *e;
//    e = docHandle.FirstChildElement("NEKTAR")
//            .FirstChildElement("CONDITIONS")
//            .Element();
    
    
    
    //std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);

    
    //===========================================================================
//    int StateVar = 0;
//    double hgrad         = metric_inputs[0];
//    double hmin          = metric_inputs[1];
//    double hmax          = metric_inputs[2];
//    double MetScale      = metric_inputs[3];
//    double hausd         = metric_inputs[4];
//    int ReadFromStats    = metric_inputs[5];
//    int RunWakRefinement = metric_inputs[6];
//    double hwake         = metric_inputs[7];
//    int niter            = metric_inputs[8];
//    int recursive	     = metric_inputs[9];
//    int extended         = metric_inputs[10];
//    StateVar         = metric_inputs[11];
    if(world_rank == 0)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "============== Metric parameters ==================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "Nproc	    = " << world_size << std::endl;
        std::cout << "hgrad     = " << inputs->hgrad << std::endl;
        std::cout << "hmin      = " << inputs->hmin << std::endl;
        std::cout << "hmax      = " << inputs->hmax << std::endl;
        std::cout << "MetScale  = " << inputs->MetScale << std::endl;
        std::cout << "Hausdorff = " << inputs->hausd << std::endl;
        std::cout << "NiterPart = " << inputs->niter << std::endl;
        if(inputs->ReadFromStats == 0)
        {
            std::cout << "Reading statistics? -> NO (5th entry in the metric.inp file is set to 0.)" << std::endl;
            std::cout << "The metric is reconstructed based on instantaneous Mach number"<<std::endl;
        }
        if(inputs->ReadFromStats == 1)
        {
            std::cout << "Reading statistics? -> YES (5th entry in the metric.inp file is set to 1.)" << std::endl;
            std::cout << "The metric is reconstructed based on the mean of Mach number."<<std::endl;

        }
        if(inputs->RunWakRefinement==0)
        {
            std::cout << "Wake refinement is switch OFF. (6th entry in the metric.inp file is set to 0. hwake, the 7th entry defined in the metric.inp file, is being ignored)" << std::endl;
            
        }
        if(inputs->RunWakRefinement==1)
        {
            std::cout << "Wake refinement is switch ON with hwake = " << inputs->hwake << "(6th entry in the metric.inp file is set to 1 and hwake is set equal to the 7th entry defined in the metric.inp file.) " << std::endl;
        }
        if(inputs->StateVar == 0)
	{
	    std::cout << "We are adapting based on the Mach number."<<std::endl;
	}
	if(inputs->StateVar == 1)
	{
	    std::cout << "We are adapting based on the static Temperature." << std::endl;
        }
        std::cout << "===================================================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "  " << std::endl;
    }
    //===========================================================================
    
    mesh* meshRead = ReadUS3DMeshData(fn_conn,fn_grid,fn_data,
                                        inputs->ReadFromStats,
                                        inputs->StateVar,
                                       comm,info);

    int Nel_loc = meshRead->ien.size();

    PrismTetraTrace* pttrace = new PrismTetraTrace(comm, 
                                                   meshRead->element2rank, 
                                                   meshRead->ife, 
                                                   meshRead->iet, 
                                                   meshRead->nElem, 
                                                   meshRead->nFace, 
                                                   meshRead->nVert);

    std::map<int,std::vector<int> >::iterator itmiv;

    std::map<int,std::vector<int> > tetras_e2v;
    std::map<int,std::vector<int> > tetras_e2f;
    std::map<int,std::vector<int> > tetras_e2e;

    std::map<int,std::vector<int> > prisms_e2v;
    std::map<int,std::vector<int> > prisms_e2f;
    std::map<int,std::vector<int> > prisms_e2e;

    std::map<int,std::vector<double> > tetras_data;
    std::map<int,std::vector<double> > prisms_data;

    std::map<int,int> loc2glob_prismv;
    std::map<int,int> glob2loc_prismv;

    std::map<int,int> loc2glob_tetrav;
    std::map<int,int> glob2loc_tetrav;

    int tetrav_loc = 0;
    int prismv_loc = 0;
    int tetra_id   = 0;
    int prism_id   = 0;
    int ntetra     = meshRead->ntetra;

    for(itmiv=meshRead->ien.begin();itmiv!=meshRead->ien.end();itmiv++)
    {
        int elid   = itmiv->first;
        int eltype = meshRead->iet[elid];

        if(eltype == 2)
        {
            tetras_e2v[elid]  = itmiv->second;
            tetras_e2f[elid]  = meshRead->ief[elid];
            tetras_e2e[elid]  = meshRead->iee[elid];
            tetras_data[elid] = meshRead->interior[elid];
            tetra_id++;
            
        }
        if(eltype == 6)
        {
            prisms_e2v[elid]  = itmiv->second;
            prisms_e2f[elid]  = meshRead->ief[elid];
            prisms_e2e[elid]  = meshRead->iee[elid];
            if(meshRead->iee[elid].size()!=5)
            {
                std::cout << "meshRead->iee[elid]; " << meshRead->iee[elid].size() << std::endl;
            }
            prisms_data[elid] = meshRead->interior[elid];
            prism_id++;
        }
    }

    // we need to pass the number of verts per element in case the partition has no elements of this type.

    //RedistributeMeshtThroughRoot(tetras,4,comm);

    RepartitionObject* tetra_repart = new RepartitionObject(meshRead, 
                                                        tetras_e2v, 
                                                        tetras_e2f,
                                                        tetras_e2e,
                                                        pttrace, 
                                                        tetras_data,
                                                        comm);



    tetras_e2v.clear();
    tetras_e2f.clear();
    tetras_e2e.clear();

    std::map<int,std::vector<double> > loc_data_t       = tetra_repart->getElement2DataMap();
    std::map<int,std::vector<int> > gE2lV_t             = tetra_repart->getGlobalElement2LocalVertMap();
    std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
    std::map<int, std::vector<double> > LocalVertsMap_t = tetra_repart->getLocalVertsMap();
    std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();


    std::map<int,std::string > varnamesGrad;

    varnamesGrad[0]     = "dUdx";
    varnamesGrad[1]     = "dUdy";
    varnamesGrad[2]     = "dUdz";
    string filename_t   = "tetra_" + std::to_string(world_rank) + ".vtu";

    std::map<int,std::string > varnames;
    varnames[0]         = "TKE";
    varnames[1]         = "Temperature";



    OutputTetraMeshPartitionVTK(filename_t, Owned_Elem_t, gE2gV_t, loc_data_t, varnames, LocalVertsMap_t);


    std::map<int,std::vector<double> > tetra_grad = ComputedUdx_LSQ_US3D_Lite(tetra_repart, 
                                                                              pttrace,
                                                                              meshRead->ghost,
                                                                              meshRead->nElem,
                                                                              1, 
                                                                              1,
                                                                              comm);

    string filename_tg = "tetraGrad_" + std::to_string(world_rank) + ".vtu";

    OutputTetraMeshPartitionVTK(filename_tg, 
                                Owned_Elem_t, 
                                gE2gV_t, 
                                tetra_grad, 
                                varnamesGrad, 
                                LocalVertsMap_t);

    tetras_e2v.clear();
    tetras_e2f.clear();
    tetras_e2e.clear();

    //=======================================================================================
    //=======================================================================================
    //=======================================================================================
    //=======================================================================================
    //=======================================================================================

    
    RepartitionObject* prism_repart = new RepartitionObject(meshRead, 
                                                            prisms_e2v,
                                                            prisms_e2f,
                                                            prisms_e2e, 
                                                            pttrace,
                                                            prisms_data,
                                                            comm);
    prisms_e2v.clear();
    prisms_e2f.clear();
    prisms_e2e.clear();

    std::map<int,std::vector<double> > loc_data_p       = prism_repart->getElement2DataMap();
    std::map<int,std::vector<int> > gE2lV_p             = prism_repart->getGlobalElement2LocalVertMap();
    std::map<int,std::vector<int> > gE2gV_p             = prism_repart->getElement2VertexMap();
    std::map<int, std::vector<double> > LocalVertsMap_p = prism_repart->getLocalVertsMap();
    std::vector<int> Owned_Elem_p                       = prism_repart->getLocElem();

    string filename_p = "prism_" + std::to_string(world_rank) + ".vtu";

    OutputPrismMeshPartitionVTK(filename_p, Owned_Elem_p, gE2gV_p, loc_data_p, varnames, LocalVertsMap_p);

    std::map<int,std::vector<double> > prism_grad = ComputedUdx_LSQ_US3D_Lite(prism_repart, 
                                                                              pttrace,
                                                                              meshRead->ghost,
                                                                              meshRead->nElem,
                                                                              1,
                                                                              1, 
                                                                              comm);

    string filename_pg = "prismGrad_" + std::to_string(world_rank) + ".vtu";

    OutputPrismMeshPartitionVTK(filename_pg, 
                                Owned_Elem_p, 
                                gE2gV_p, 
                                prism_grad, 
                                varnamesGrad, 
                                LocalVertsMap_p);


    //=======================================================================================
    //=======================================================================================
    //=======================================================================================
    //=======================================================================================
    //=======================================================================================






    //Mesh_Topology_Lite* meshTopo = new Mesh_Topology_Lite(tetra_repart,comm);
    /*
    char* filename = "mesh.vtk";

    vtkNew<vtkPoints> points;
    points->InsertNextPoint(0, 0, 0);
    points->InsertNextPoint(1, 0, 0);
    points->InsertNextPoint(1, 1, 0);
    points->InsertNextPoint(0, 1, 1);

    vtkNew<vtkTetra> tetra;
    tetra->GetPointIds()->SetId(0, 0);
    tetra->GetPointIds()->SetId(1, 1);
    tetra->GetPointIds()->SetId(2, 2);
    tetra->GetPointIds()->SetId(3, 3);

    vtkSmartPointer<vtkUnstructuredGrid> vtkmesh =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

        vtkSmartPointer<vtkCellArray> cellArray =
    vtkSmartPointer<vtkCellArray>::New();
    cellArray->InsertNextCell(tetra);
    vtkmesh->SetPoints(points);
    vtkmesh->SetCells(VTK_TETRA, cellArray);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(vtkmesh);
    writer->Write();



    
    int Nve       = us3d->xcn->getNglob();
    
    int Nel_part  = us3d->ien->getNrow();
    
    Array<double>* Ui = new Array<double>(Nel_part,1);
    Array<double>* TKEi;
    int varia = 4;
    double TKE, MState;
    
    if(inputs->ReadFromStats==0)
    {
        for(int i=0;i<Nel_part;i++)
        {
            MState   = us3d->interior->getVal(i,0);
            Ui->setVal(i,0,MState);
        }
    }
    
    if(inputs->ReadFromStats==1)
    {
        TKEi = new Array<double>(Nel_part,1);

        for(int i=0;i<Nel_part;i++)
        {
            TKE      = us3d->interior->getVal(i,0);
            MState   = us3d->interior->getVal(i,1);
            Ui->setVal(i,0,MState);
            TKEi->setVal(i,0,TKE);
        }
    }
    
    
    delete us3d->interior;
 
    Array<double>* gB = new Array<double>(us3d->ghost->getNrow(),1);
    for(int i=0;i<us3d->ghost->getNrow();i++)
    {
        gB->setVal(i,0,us3d->ghost->getVal(i,varia));
    	//std::cout << "ghost " << i << " " << us3d->ghost->getVal(i,varia) << std::endl;
    }
    int ngho = us3d->ghost->getNrow();
    int ngval;
    MPI_Allreduce(&ngho, &ngval, 1, MPI_INT, MPI_MAX, comm);
       
    ParallelState* ien_pstate               = new ParallelState(us3d->ien->getNglob(),comm);
    ParallelState* ife_pstate               = new ParallelState(us3d->ifn->getNglob(),comm);
    
    ParallelState_Parmetis* parmetis_pstate = new ParallelState_Parmetis(us3d->ien,us3d->elTypes,us3d->ie_Nv,comm);
    ParallelState* xcn_pstate               = new ParallelState(us3d->xcn->getNglob(),comm);
    
    clock_t t;
    double tn = 0.0;
    t = clock();
    
    
    

    RedistributePartitionObject* tetra_distri = new RedistributePartitionObject(us3d,
    																			tetrahedra,
																				iferank_map,
                                                                                ief_map,
                                                                                ifn_map,
                                                                                ife_map,
                                                                                ifref_map,
                                                                                ushell,
                                                                                hess_vmap_new, comm);
    
    
    Array<int>* element2node                        = tetra_distri->GetElement2NodeMap();
    std::map<int,Array<double>* > metric            = tetra_distri->GetVert2MetricMap();
    int** ifc_tria_glob                             = tetra_distri->GetFace2GlobalNode();
    int** ifc_tria_loc                              = tetra_distri->GetFace2LocalNode();
    int nFaces                                      = tetra_distri->GetNBoundaryFaces();
    std::vector<std::vector<double> > locVs                        = tetra_distri->GetLocalVertices();
    std::vector<int> faces4parmmg                   = tetra_distri->GetFaces4ParMMG();
    std::map<int,std::vector<int> > face2node                    = tetra_distri->GetFace2NodeMap();
    std::map<int,std::vector<int> > face2element    = tetra_distri->GetFace2ElementMap();
    std::map<int,int> globV2locV                    = tetra_distri->GetGlobalVert2LocalVertMap();
    std::map<int,int> locV2globV                    = tetra_distri->GetLocalVert2GlobalVertMap();
    int ncomm                                       = tetra_distri->GetNcomm();
    int* color_face                                 = tetra_distri->GetColorFace();
    //int** face2globnode                           = tetra_distri->GetFace2GlobalNode();
    int *ntifc                                      = tetra_distri->GetNFacesPerColor();
    std::map<int,int> locShF2globShF                = tetra_distri->GetLocalSharedFace2GlobalSharedFace();
    std::map<int,int> face2ref                      = tetra_distri->GetFace2RefMap();
    std::map<int,int> shell_tet2hybF                = tetra_distri->GetShellTet2HybFaceMap();
    std::map<int,int> shellvert2ref                 = tetra_distri->GetShellVert2RefMap_Global();
    std::map<int,std::set<int> > shellface2vertref  = tetra_distri->GetShellFace2VertRefMap();
    std::map<int,int> shellvert2ref_local           = tetra_distri->GetShellVert2RefMap_Local();
    std::map<int,int> tetF2hybF                     = tetra_distri->GetTetF2HybFMap();
    std::map<int,int> tetV2tagV						= tetra_distri->GetTet2TagVertMap();
    std::map<int,int> shellvertOriginalTag2ref_Glob = tetra_distri->GetShellVertTag2RefMap_Global();
    std::map<int,std::vector<int> > bndref2face     = tetra_distri->GetBndRef2FaceMap();
    std::map<int,std::vector<double> > shellVertCoord2Ref          = tetra_distri->GetShellVertCoords2RefMap_Global();
    std::map<int,int> shellvertTag2ref;
    std::map<int,int> shellvertTag2ref2;

    */
    MPI_Finalize();
    
}

