
#include <chrono>
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
    std::cout << "tracing " << std::endl;
    PrismTetraTrace* pttrace = new PrismTetraTrace(comm, 
                                                   meshRead->element2rank, 
                                                   meshRead->ife,
                                                   meshRead->ifn, 
                                                   meshRead->iet, 
                                                   meshRead->nElem, 
                                                   meshRead->nFace, 
                                                   meshRead->nVert);

    std::cout << " end tracing " << std::endl;
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

    time_t start, end; 
 
    /* You can call it like this : start = time(NULL); 
    in both the way start contain total time in seconds 
    since the Epoch. */
    time(&start); 
    RepartitionObject* tetra_repart = new RepartitionObject(meshRead, 
                                                        tetras_e2v, 
                                                        tetras_e2f,
                                                        tetras_e2e,
                                                        pttrace, 
                                                        tetras_data,
                                                        comm);
    time(&end); 
    double time_taken = double(end - start); 
    cout << "Time taken to repartion tetrahedera is : " << fixed 
        << time_taken << setprecision(16); 
    cout << " sec " << endl;

    tetras_e2v.clear();
    tetras_e2f.clear();
    tetras_e2e.clear();

    std::map<int,std::vector<double> > loc_data_t       = tetra_repart->getElement2DataMap();
    std::map<int,std::vector<int> > gE2lV_t             = tetra_repart->getGlobalElement2LocalVertMap();
    std::map<int,std::vector<int> > gE2gV_t             = tetra_repart->getElement2VertexMap();
    std::map<int, std::vector<double> > LocalVertsMap_t = tetra_repart->getLocalVertsMap();
    std::vector<int> Owned_Elem_t                       = tetra_repart->getLocElem();

    std::set<int> test_set;
    for(int i=0;i<Owned_Elem_t.size();i++)
    {
        int elid = Owned_Elem_t[i];
        for(int j=0;j<4;j++)
        {
            int vid = gE2gV_t[elid][j];

            if(test_set.find(vid)==test_set.end())
            {
                test_set.insert(vid);
            }

        }
    }

    std::cout << "nvert tet = " << test_set.size() << std::endl;

    std::map<int,std::string > varnamesGrad;

    varnamesGrad[0]     = "dUdx";
    varnamesGrad[1]     = "dUdy";
    varnamesGrad[2]     = "dUdz";
    string filename_t   = "tetra_" + std::to_string(world_rank) + ".vtu";

    std::map<int,std::string > varnames;
    varnames[0]         = "TKE";
    varnames[1]         = "Temperature";


    if(world_rank == 0)
    {
        std::cout << "++++++++++++++++++++ BEFORE ++++++++++++++++++++" << world_rank << std::endl;
        std::cout << "              Owned_Elem_t " << Owned_Elem_t.size() << std::endl;
        std::cout << "              gE2gV_t " << gE2gV_t.size() << std::endl;
        std::cout << "              loc_data_t " << loc_data_t.size() << std::endl;
        std::cout << "              LocalVertsMap_t " << loc_data_t.size() << std::endl;  
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }


    

    OutputTetraMeshPartitionVTK(filename_t, Owned_Elem_t, gE2gV_t, loc_data_t, varnames, LocalVertsMap_t);


    time(&start); 
    // std::map<int,std::vector<double> > tetra_grad = ComputedUdx_LSQ_US3D_Lite(tetra_repart, 
    //                                                                           pttrace,
    //                                                                           meshRead->ghost,
    //                                                                           meshRead->nElem,
    //                                                                           1, 
    //                                                                           1,
    //                                                                           comm);
    // time(&end); 
    // time_taken = double(end - start); 
    // cout << "Time taken to calculate gradients for tetrahedra is : " << fixed 
    //     << time_taken << setprecision(16); 
    // cout << " sec " << endl;
    // string filename_tg = "tetraGrad_" + std::to_string(world_rank) + ".vtu";

    // OutputTetraMeshPartitionVTK(filename_tg, 
    //                             Owned_Elem_t, 
    //                             gE2gV_t, 
    //                             tetra_grad, 
    //                             varnamesGrad, 
    //                             LocalVertsMap_t);

    tetras_e2v.clear();
    tetras_e2f.clear();
    tetras_e2e.clear();

    //=======================================================================================
    //=======================================================================================
    //=======================================================================================
    //=======================================================================================
    //=======================================================================================

    
    // RepartitionObject* prism_repart = new RepartitionObject(meshRead, 
    //                                                         prisms_e2v,
    //                                                         prisms_e2f,
    //                                                         prisms_e2e, 
    //                                                         pttrace,
    //                                                         prisms_data,
    //                                                         comm);
    // prisms_e2v.clear();
    // prisms_e2f.clear();
    // prisms_e2e.clear();

    // std::map<int,std::vector<double> > loc_data_p       = prism_repart->getElement2DataMap();
    // std::map<int,std::vector<int> > gE2lV_p             = prism_repart->getGlobalElement2LocalVertMap();
    // std::map<int,std::vector<int> > gE2gV_p             = prism_repart->getElement2VertexMap();
    // std::map<int, std::vector<double> > LocalVertsMap_p = prism_repart->getLocalVertsMap();
    // std::vector<int> Owned_Elem_p                       = prism_repart->getLocElem();

    // prism_repart->buildInteriorSharedAndBoundaryFacesMaps(comm,pttrace, meshRead->ranges_id);
    // std::map<int,std::vector<int> > sharedFaces_prisms   = prism_repart->getSharedFaceMap();
    // std::map<int,std::vector<int> > interiorFaces_prisms = prism_repart->getInteriorFaceMap();
    // std::map<int,std::vector<int> > boundaryFaces_prisms = prism_repart->getBoundaryFaceMap();
    // // prism_repart->buildTag2GlobalElementMapsAndFace2LeftRightElementMap(comm,pttrace, meshRead->ranges_id);
    // tetra_repart->buildInteriorSharedAndBoundaryFacesMaps(comm,pttrace, meshRead->ranges_id);
    // std::map<int,std::vector<int> > sharedFaces_tetras   = tetra_repart->getSharedFaceMap();
    // std::map<int,std::vector<int> > interiorFaces_tetras = tetra_repart->getInteriorFaceMap();
    // std::map<int,std::vector<int> > boundaryFaces_tetras = tetra_repart->getBoundaryFaceMap();

    // string filename_p = "prism_" + std::to_string(world_rank) + ".vtu";

    // OutputPrismMeshPartitionVTK(filename_p, Owned_Elem_p, gE2gV_p, loc_data_p, varnames, LocalVertsMap_p);

    // std::map<int,std::vector<double> > prism_grad = ComputedUdx_LSQ_US3D_Lite(prism_repart, 
    //                                                                           pttrace,
    //                                                                           meshRead->ghost,
    //                                                                           meshRead->nElem,
    //                                                                           1,
    //                                                                           1, 
    //                                                                           comm);

    // string filename_pg = "prismGrad_" + std::to_string(world_rank) + ".vtu";

    // OutputPrismMeshPartitionVTK(filename_pg, 
    //                             Owned_Elem_p, 
    //                             gE2gV_p, 
    //                             prism_grad, 
    //                             varnamesGrad, 
    //                             LocalVertsMap_p);

    prisms_e2v.clear();
    prisms_e2f.clear();
    prisms_e2e.clear();
    tetra_repart->buildInteriorSharedAndBoundaryFacesMaps(comm,pttrace, meshRead->ranges_id);
    tetra_repart->buildTag2GlobalElementMapsAndFace2LeftRightElementMap(comm,pttrace, meshRead->ranges_id);

    std::map<int,std::map<int,int> > trace_elem         = pttrace->GetTrace();
    std::map<int,std::vector<int> > trace_verts         = pttrace->GetTraceVerts();
    std::map<int,int> unique_trace_verts2refmap         = pttrace->GetUniqueTraceVerts2RefMap();
    std::map<int,std::vector<int> > leftright_trace     = pttrace->GetLeftRightElements();
    std::map<int,int> trace_ref                         = pttrace->GetTraceRef();
    // std::map<int,int> tag2element_trace_Prisms          = prism_repart->getTag2ElementTrace();
    // Note that the boundaryFaces_tetra and the boundaryFaces2Ref_tetra data DOES include the trace data.
    std::map<int,int> boundaryFaces_tetra               = tetra_repart->getBoundaryFaces();
    std::map<int,int> boundaryFaces2Ref_tetra           = tetra_repart->getBoundaryFaces2Ref();
    // Note that this boundaryFacesMap_tetra data DOES NOT include the trace data.
    std::map<int,std::vector<int> > boundaryFacesMap_tetra = tetra_repart->getBoundaryFaceMap();

    FaceSetPointer m_PMMG_RefsOnShell_2_Prism;
    int nshell = 0;
    int foundU = 0;


    // for(itmiv=trace_verts.begin();itmiv!=trace_verts.end();itmiv++)
    // {
    //     int fhyb       = itmiv->first;
    //     int Etettag    = leftright_trace[fhyb][0];
    //     int Eprismtag  = leftright_trace[fhyb][1];

    //     int EprismNew  = tag2element_trace_Prisms[Eprismtag];
    
    //     if(trace_verts.find(fhyb)!=trace_verts.end())
    //     {
    //         std::vector<int>::iterator its;
    //         std::vector<int> refs(trace_verts[fhyb].size());
    //         int c = 0;
    //         for(its=trace_verts[fhyb].begin();its!=trace_verts[fhyb].end();its++)
    //         {
    //             int vid = *its;
    //             int ref = unique_trace_verts2refmap[vid];
    //             refs[c] = ref;
    //             c++;
    //         }
    //         FaceSharedPtr RefFacePointer = std::shared_ptr<NekFace>(new NekFace(refs));
    //         pair<FaceSetPointer::iterator, bool> testInsPointer;
    //         testInsPointer = m_PMMG_RefsOnShell_2_Prism.insert(RefFacePointer);
            
    //         if(testInsPointer.second)
    //         {
    //             (*testInsPointer.first)->SetFaceLeftElement(EprismNew);
    //         }
    //     }
    // }

    std::map<int,int> lv2gvID  = tetra_repart->getLocalVertex2GlobalVertexID();
    std::map<int,int> gv2lvID  = tetra_repart->getGlobalVertex2LocalVertexID();
    std::map<int,int> le2geID  = tetra_repart->getTag2GlobElementID();
    std::map<int,int> tag2leID = tetra_repart->getTag2LocalElementID();
    std::map<int,int> le2tagID = tetra_repart->getLocal2TagElementID();

    tetra_repart->UpdateGlobalIDs(comm,pttrace);

    tetra_repart->GetFace2RankMesh(comm);

    std::vector<int> face4parmmg = tetra_repart->getFace4ParMMG();
    std::map<int,int> global2tagF = tetra_repart->getGlobal2TagFMap();
    std::cout << "huh " << boundaryFaces_tetra.size() << " " << face4parmmg.size() << " rank " << world_rank << std::endl; 
    string filename_pg = "tetraData_" + std::to_string(world_rank) + ".vtu";


    if(world_rank == 0)
    {
        std::cout << "++++++++++++++++++++ AFTER ++++++++++++++++++++" << world_rank << std::endl;
        std::cout << "              Owned_Elem_t " << Owned_Elem_t.size() << std::endl;
        std::cout << "              gE2gV_t " << gE2gV_t.size() << std::endl;
        std::cout << "              loc_data_t " << loc_data_t.size() << std::endl;
        std::cout << "              LocalVertsMap_t " << loc_data_t.size() << std::endl;  
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }

    OutputTetraMeshPartitionVTK(filename_pg, 
                                Owned_Elem_t, 
                                gE2gV_t, 
                                loc_data_t, 
                                varnames, 
                                LocalVertsMap_t);



    int nVertices   = lv2gvID.size();
    int nTetrahedra = le2geID.size();
    int nEdges      = 0;
    int nTriangles  = face4parmmg.size();

    PMMG_pParMesh   parmesh;
    PMMG_Init_parMesh(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_pMesh,PMMG_ARG_pMet,
                      PMMG_ARG_dim,3,PMMG_ARG_MPIComm,MPI_COMM_WORLD,
                      PMMG_ARG_end);

    if ( PMMG_Set_meshSize(parmesh,nVertices,nTetrahedra,0,nTriangles,0,nEdges) != 1 ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }

    //PMMG_Set_metSize(PMMG_pParMesh parmesh,int typEntity,int np,int typSol)
    if ( PMMG_Set_metSize(parmesh,MMG5_Vertex,nVertices,MMG5_Tensor) != 1 ) exit(EXIT_FAILURE);

    int vrefmax = -1;
    // std::map<int,int>::iterator itm;
    // for(itm=lv2gvID.begin();itm!=lv2gvID.end();itm++)
    // {
    //     std::cout << "gvid =  " << itm->first << " " << itm->second << " " << nVertices << std::endl;
    // }

    for ( k=0; k<nVertices; ++k )
    {
        int gvid  = lv2gvID[k];

        //std::cout << "gvid " << gvid << std::endl;
        double vx = LocalVertsMap_t[gvid][0];
        double vy = LocalVertsMap_t[gvid][1];
        double vz = LocalVertsMap_t[gvid][2];

        int vref = 86;

        if ( PMMG_Set_vertex(parmesh,vx,vy,vz, vref, k+1) != 1 )
        {
        MPI_Finalize();
        exit(EXIT_FAILURE);
        }

        std::vector<double> tensor(6,0);
        tensor[0] = 1.0;
        tensor[1] = 0.0;
        tensor[2] = 0.0;
        tensor[3] = 1.0;
        tensor[4] = 0.0;
        tensor[5] = 1.0;

        if(PMMG_Set_tensorMet(parmesh,tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],k+1)!=1)
        {
         MPI_Finalize();
         exit(EXIT_FAILURE);
        }
    }


        int v0,v1,v2,v3;
    int v0l,v1l,v2l,v3l;
    int teller = 0;
    int refer  = 0;
    int cref36 = 0;
    //double* c0 = new double[3];
    int iref;
    int c13 = 0;
    int suc = 0;
    int buggi = 0;
    int shelfound = 0;
    int shelfound2 = 0;
    int vertref = 86;
    int flippie = -1;
    
    int locs = 0;
    
    std::map<int,int> shell_g2l;
    std::vector<std::vector<double> > unshellVin;
    std::vector<std::vector<int> > unshellTin;
    int outflowbc  = 0;
    int inflowbc   = 0;
    int tracebc    = 0;
    int symmetrybc = 0;
    int wallbc = 0;
    int nothere = 0;

    std::cout << "nTriangles set into ParMMG structure " << nTriangles << " " << face4parmmg.size()<< std::endl;

    std::map<int,std::vector<int> > f2vmap = tetra_repart->getFace2VertexMap();
    std::map<int,std::vector<int> > f2refmap = tetra_repart->getFace2RefMap();
    int ref2 = 0;
    for ( k=0; k<nTriangles; ++k )
    {
        
        int faceID      = face4parmmg[k];
        int tagFaceID   = global2tagF[faceID];
        int ref = f2refmap[tagFaceID][0];
        
        if(ref==13)
        {
            for(int s=0;s<3;s++)
            {
                int ref = 13;
                int gvid  = trace_verts[tagFaceID][s];

                // int globalVid = global2tagV[gvid];

                double vx = LocalVertsMap_t[gvid][0];
                double vy = LocalVertsMap_t[gvid][1];
                double vz = LocalVertsMap_t[gvid][2];

                int lvid = gv2lvID[gvid];
                
                if(unique_trace_verts2refmap.find(gvid)!=unique_trace_verts2refmap.end())
                {
                    vertref = unique_trace_verts2refmap[gvid];
                }
                else
                {
                    std::cout << "Warning:: reference value on trace vertID " << gvid << " is wrong!" << std::endl;
                }
            }
            
            int gvid0   = trace_verts[tagFaceID][0];
            int gvid1   = trace_verts[tagFaceID][1];
            int gvid2   = trace_verts[tagFaceID][2];

            v0l         = gv2lvID[gvid0];
            v1l         = gv2lvID[gvid1];
            v2l         = gv2lvID[gvid2];
        
            if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,ref,k+1) != 1 )
            {
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }

            PMMG_Set_requiredTriangle( parmesh, k+1 );

            tracebc++;
        }
        else
        {
            std::cout << "ref  = " << ref << std::endl;
            if(f2vmap.find(tagFaceID)!=f2vmap.end())
            {
                
                std::cout << tagFaceID << std::endl;
            }
            else
            {
                std::cout << "tagFaceID not found " << tagFaceID << std::endl;
            }


            if(ref==2)
            {
                ref2++;
            }
            int gvid0    = f2vmap[tagFaceID][0];
            int gvid1    = f2vmap[tagFaceID][1];
            int gvid2    = f2vmap[tagFaceID][2];
       
            //std::cout << "local vs " <<  v0l << " " << v1l << " " << v2l << std::endl;
            v0l         = gv2lvID[gvid0];
            v1l         = gv2lvID[gvid1];
            v2l         = gv2lvID[gvid2];
        
            if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,ref,k+1) != 1 )
            {
                MPI_Finalize();
                exit(EXIT_FAILURE);
            }
        }
        // if(trace_verts.find(tagFaceID)!=trace_verts.end())// Handle prism-tetra interfaces/traces differently
        // {
        //     for(int s=0;s<3;s++)
        //     {
        //         int ref = 13;
        //         int gvid  = trace_verts[tagFaceID][s];

        //         // int globalVid = global2tagV[gvid];

        //         double vx = LocalVertsMap_t[gvid][0];
        //         double vy = LocalVertsMap_t[gvid][1];
        //         double vz = LocalVertsMap_t[gvid][2];

        //         int lvid = gv2lvID[gvid];
                
        //         if(unique_trace_verts2refmap.find(gvid)!=unique_trace_verts2refmap.end())
        //         {
        //             vertref = unique_trace_verts2refmap[gvid];
        //         }
        //         else
        //         {
        //             std::cout << "Warning:: reference value on trace vertID " << gvid << " is wrong!" << std::endl;
        //         }
        //     }
            
        //     int gvid0   = trace_verts[tagFaceID][0];
        //     int gvid1   = trace_verts[tagFaceID][1];
        //     int gvid2   = trace_verts[tagFaceID][2];

        //     v0l         = gv2lvID[gvid0];
        //     v1l         = gv2lvID[gvid1];
        //     v2l         = gv2lvID[gvid2];
        
        //     if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,ref,k+1) != 1 )
        //     {
        //         MPI_Finalize();
        //         exit(EXIT_FAILURE);
        //     }

        //     PMMG_Set_requiredTriangle( parmesh, k+1 );

        //     tracebc++;
        // }
        
        // else
        // {
            
            //std::cout << "ref  = " << ref << std::endl;
            // if(f2vmap.find(tagFaceID)!=f2vmap.end())
            // {
                
            //     std::cout << tagFaceID << std::endl;
            // }
            // else
            // {
            //     std::cout << "tagFaceID not found " << tagFaceID << std::endl;
            // }


            // if(ref==2)
            // {
            //     ref2++;
            // }
            // int gvid0    = f2vmap[tagFaceID][0];
            // int gvid1    = f2vmap[tagFaceID][1];
            // int gvid2    = f2vmap[tagFaceID][2];
       
            // //std::cout << "local vs " <<  v0l << " " << v1l << " " << v2l << std::endl;
            // v0l         = gv2lvID[gvid0];
            // v1l         = gv2lvID[gvid1];
            // v2l         = gv2lvID[gvid2];
        
            // if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,ref,k+1) != 1 )
            // {
            //     MPI_Finalize();
            //     exit(EXIT_FAILURE);
            // }

        //}
        
        // if(ref==3)
        // {
        //     wallbc++;
        // }
        // if(ref==36)
        // {
        //     outflowbc++;
        // }
        // if(ref==7)
        // {
        //     symmetrybc++;
        // }
        // if(ref==10)
        // {
        //     inflowbc++;
        // }

        // if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,refer,k+1) != 1 )
        // {
        //     MPI_Finalize();
        //     exit(EXIT_FAILURE);
        // }
        
    }

    std::cout << "nothere " << nothere << "  " << ref2 << " "  <<tracebc<< std::endl; 
    //std::cout << "Set triangles is  " << k << std::endl;
    
    // int outflowbc_sum   = 0;
    // int inflowbc_sum    = 0;
    // int tracebc_sum     = 0;
    // int wallbc_sum      = 0;
    // int symmetrybc_sum  = 0;

    // MPI_Allreduce(&outflowbc, &outflowbc_sum, 1, MPI_INT, MPI_SUM, comm);
    // MPI_Allreduce(&inflowbc, &inflowbc_sum, 1, MPI_INT, MPI_SUM, comm);
    // MPI_Allreduce(&tracebc, &tracebc_sum, 1, MPI_INT, MPI_SUM, comm);
    // MPI_Allreduce(&wallbc, &wallbc_sum, 1, MPI_INT, MPI_SUM, comm);
    // MPI_Allreduce(&symmetrybc, &symmetrybc_sum, 1, MPI_INT, MPI_SUM, comm);

    // if(world_rank == 0)
    // {
    //     std::cout << "=============Boundary Face Stats for Tetrahedra ============" << std::endl;
    //     std::cout << "                  Outflow: " << outflowbc_sum << std::endl;
    //     std::cout << "                  Inflow: " << inflowbc_sum << std::endl;
    //     std::cout << "                  Trace: " << tracebc_sum << std::endl;
    //     std::cout << "                  Wall: " << wallbc_sum << std::endl;
    //     std::cout << "                  Symmetry: " << symmetrybc_sum << std::endl;
    //     std::cout << "============================================================" << std::endl;
    // }

    int API_mode = 0;
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };


    // tetra_repart->buildParMMGCommunicationMaps(comm);

    int** ifc_tria_glob = tetra_repart->getParMMGCommFace2GlobalVertMap();
    int** ifc_tria_loc  = tetra_repart->getParMMGCommFace2LocalVertMap();
    int* color_face     = tetra_repart->getParMMGCommColorFace();
    int *ntifc          = tetra_repart->getParMMGCommNFacesPerColor();
    int ncomm           = tetra_repart->getParMMGNComm();


    std::cout << " ncomm " << ncomm << std::endl;
    
    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);

   
    for(int icomm=0; icomm<ncomm; icomm++ )
    {
      // Set nb. of entities on interface and rank of the outward proc
     std::cout << "icomm " << icomm << " " << ntifc[icomm] << std::endl;
      ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                             color_face[icomm],
                                             ntifc[icomm]);

        std::cout << "ntifc[icomm] " << ntifc[icomm] << " rank " << world_rank << " ";
        //Set local and global index for each entity on the interface
      ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                               ifc_tria_loc[icomm],
                                               ifc_tria_glob[icomm], 1);
    }

    
    std::cout << std::endl;
    
    // //==========================Clear from here==============================
    
    std::map<int,std::vector<int> >::iterator ittet;
    k = 0;
    std::cout << "nTetrahedra  " << nTetrahedra << " " << le2geID.size() << std::endl;

    for ( int t = 0;t < nTetrahedra; t++  )
    {
        //From local element ID get the tag ID (initial global ID).
        int gEid = le2tagID[t];
        //Using the tag Element ID, get the tag Vertex ID.
        v0 = gE2gV_t[gEid][0];
        v1 = gE2gV_t[gEid][1];
        v2 = gE2gV_t[gEid][2];
        v3 = gE2gV_t[gEid][3];
        //From the tag vertex ID, get the local vertex ID.
        v0l = gv2lvID[v0];
        v1l = gv2lvID[v1];
        v2l = gv2lvID[v2];
        v3l = gv2lvID[v3];
        
        std::vector<double> P(4*3);
        P[0*3+0]=LocalVertsMap_t[v0][0];   P[0*3+1]=LocalVertsMap_t[v0][1];    P[0*3+2]=LocalVertsMap_t[v0][2];
        P[1*3+0]=LocalVertsMap_t[v1][0];   P[1*3+1]=LocalVertsMap_t[v1][1];    P[1*3+2]=LocalVertsMap_t[v1][2];
        P[2*3+0]=LocalVertsMap_t[v2][0];   P[2*3+1]=LocalVertsMap_t[v2][1];    P[2*3+2]=LocalVertsMap_t[v2][2];
        P[3*3+0]=LocalVertsMap_t[v3][0];   P[3*3+1]=LocalVertsMap_t[v3][1];    P[3*3+2]=LocalVertsMap_t[v3][2];

        double Vtet = GetQualityTetrahedra(P);
        //std::cout << "Vtet " << Vtet << std::endl; 
        if(Vtet<0.0)
        {
            std::cout << " negative volume in Element " << t << " on rank " << world_rank  <<std::endl;
        }
        if ( PMMG_Set_tetrahedron(parmesh,v0l+1,v1l+1,v2l+1,v3l+1,1.0,t+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    }
    
    
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_niter, inputs->niter ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hausd, inputs->hausd) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    
    if ( PMMG_Set_dparameter(parmesh,PMMG_DPARAM_hgrad, inputs->hgrad) != 1 )
    {
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    
    if( !PMMG_Set_dparameter( parmesh,  PMMG_DPARAM_hgradreq , -1.0 ) ){
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
 
    int nVerticesIN   = 0;
    int nTetrahedraIN = 0;
    int nTrianglesIN  = 0;
    int nEdgesIN      = 0;
    
    if ( PMMG_Get_meshSize(parmesh,&nVerticesIN,&nTetrahedraIN,NULL,&nTrianglesIN,NULL,
                           &nEdgesIN) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }
    
    // std::cout << "nVertices input " << nVerticesIN << std::endl;
    // std::cout << "nTetrahedra input " << nTetrahedraIN << std::endl;
    // std::cout << "nTriangles input " << nTrianglesIN << std::endl;
    // std::cout << "nEdges input " << nEdgesIN << std::endl;
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_globalNum, 1 ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };

    DistributedParallelState* distTetraB = new DistributedParallelState(nTetrahedra,comm);
    int ToTElements_bef                  = distTetraB->getNel();
    
    if(world_rank == 0)
    {
        std::cout << "Total elements = "  << ToTElements_bef << std::endl;
    }

    int ierlib = PMMG_parmmglib_distributed( parmesh );
    
    



    MPI_Finalize();
    
}

