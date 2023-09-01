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
#include <iomanip>

#define MAX2(a,b)      (((a) > (b)) ? (a) : (b))
#define MAX4(a,b,c,d)  (((MAX2(a,b)) > (MAX2(c,d))) ? (MAX2(a,b)) : (MAX2(c,d)))

// This is basically textbook recursive merge sort using std::merge_inplace
// but it considers the offsets of segments that are already sorted




std::vector<std::vector<double> > ReadReferenceData(int world_rank)
{
    std::vector<std::vector<double> > output;
    std::ifstream fin_v;
    fin_v.open("testdata/compareValues_"+ std::to_string(world_rank) + ".txt");
    
    std::vector<double> row_v(6);
    while(fin_v >> row_v[0] >> row_v[1] >> row_v[2] >> row_v[3] >> row_v[4] >> row_v[5])
    {
        output.push_back(row_v);
    }
    fin_v.close();
    
    return output;
    
}

std::map<int,Array<double>*> ReadReferenceData2(int world_rank)
{
    std::map<int,Array<double>*> output;
    std::ifstream fin_v;
    fin_v.open("testdata/UVariaValues_"+ std::to_string(world_rank) + ".txt");
    
    std::vector<double> row_v(2);
    while(fin_v >> row_v[0] >> row_v[1])
    {
        Array<double>* entry = new Array<double>(1,1);
        entry->setVal(0,0,row_v[1]);
        int key = (int) row_v[0];
        output[key] = entry;
        
    }
    fin_v.close();
    
    return output;
    
}






void OutputMesh_PMMG(int nV, double* VertOUT, int nE, int* tetraOUT, string fname)
{
    int pos;
    
    std::ofstream myfile;
    myfile.open(fname);
    myfile << "TITLE=\"new_volume.tec\"" << std::endl;
    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile <<"ZONE N = " << nV << ", E = " << nE << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

    for(int i=0;i<nV;i++)
    {
        pos = 3*i;
        myfile << VertOUT[pos] << " " <<VertOUT[pos+1] << " " << VertOUT[pos+2] <<  std::endl;
    }
    for(int i=0;i<nE;i++)
    {
        pos=4*i;
        
        myfile << tetraOUT[pos] << " " << tetraOUT[pos+1]  << " " << tetraOUT[pos+2]  << " " << tetraOUT[pos+3]  << std::endl;
        
    }
    myfile.close();

}






void OutputMesh_PMMG_V2(int nV, double* VertOUT, int nE, int* tetraOUT, string fname)
{
    int pos;
    
    std::vector<std::vector<int> > plotnodes;
    for(int i=0;i<nE;i++)
    {
        pos=4*i;
        if(VertOUT[(tetraOUT[pos]-1)*3+2]<0.0 &&
           VertOUT[(tetraOUT[pos+1]-1)*3+2]<0.0 &&
           VertOUT[(tetraOUT[pos+2]-1)*3+2]<0.0 &&
           VertOUT[(tetraOUT[pos+3]-1)*3+2]<0.0)
        {
            std::vector<int> roww(4);
            roww[0] = tetraOUT[pos];
            roww[1] = tetraOUT[pos+1];
            roww[2] = tetraOUT[pos+2];
            roww[3] = tetraOUT[pos+3];
            
            plotnodes.push_back(roww);
        }
        
    }
        
    if(plotnodes.size() > 0 )
    {
        std::ofstream myfile;
        myfile.open(fname);
        myfile << "TITLE=\"new_volume.tec\"" << std::endl;
        myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
        myfile <<"ZONE N = " << nV << ", E = " << plotnodes.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

        for(int i=0;i<nV;i++)
        {
            pos = 3*i;
            myfile << VertOUT[pos] << " " <<VertOUT[pos+1] << " " << VertOUT[pos+2] <<  std::endl;
        }
        for(int i=0;i<plotnodes.size();i++)
        {
            pos=4*i;
    //        if(VertOUT[tetraOUT[pos]*3+2]>0.25 &&
    //           VertOUT[tetraOUT[pos+1]*3+2]>0.25 &&
    //           VertOUT[tetraOUT[pos+2]*3+2]>0.25 &&
    //           VertOUT[tetraOUT[pos+3]*3+2]>0.25)
    //        {
            
                //myfile << tetraOUT[pos] << " " << tetraOUT[pos+1]  << " " << tetraOUT[pos+2]  << " " << tetraOUT[pos+3]  << std::endl;
            myfile << plotnodes[i][0] << " " << plotnodes[i][1]  << " " << plotnodes[i][2]  << " " << plotnodes[i][3]  << std::endl;
    //        }
            
        }
        myfile.close();
    }
    
}


std::map<int,int> AllGatherMap(std::map<int,int> mappie, MPI_Comm mpi_comm)
{
    int mapSizeLoc = mappie.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    int mapSizeTot = distrimap->getNel();
    
    int* key_loc = new int[mapSizeLoc];
    int* val_loc = new int[mapSizeLoc];
    int* key_tot = new int[mapSizeTot];
    int* val_tot = new int[mapSizeTot];
    int i = 0;
    
    std::map<int,int>::iterator itred;
    for(itred=mappie.begin();itred!=mappie.end();itred++)
    {
        key_loc[i] = itred->first;
        val_loc[i] = itred->second;
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    
    MPI_Allgatherv(key_loc,
                   mapSizeLoc,
                   MPI_INT,
                   key_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    
    MPI_Allgatherv(val_loc,
                   mapSizeLoc,
                   MPI_INT,
                   val_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    int key,val;
    std::map<int,int> mappie_glob;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = key_tot[i];
        val = val_tot[i];
        
        if(mappie_glob.find(key)==mappie_glob.end())
        {
        	mappie_glob[key] = val;
        }
    }
    
    return mappie_glob;
}



std::map<int,std::vector<double> > AllGatherMapDoubleVec(std::map<int,std::vector<double> > mappie, MPI_Comm mpi_comm)
{
    int mapSizeLoc = mappie.size();
    DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,mpi_comm);
    
    DistributedParallelState* distrimapVal = new DistributedParallelState(mapSizeLoc*3,mpi_comm);
    
    int mapSizeTot = distrimap->getNel();
    int* key_loc = new int[mapSizeLoc];
    double* val_loc = new double[mapSizeLoc*3];
    int* key_tot = new int[mapSizeTot];
    double* val_tot = new double[mapSizeTot*3];
    int i = 0;
    
    std::map<int,std::vector<double> >::iterator itred;
    for(itred=mappie.begin();itred!=mappie.end();itred++)
    {
        key_loc[i] = itred->first;
        int nrow   = itred->second.size();
        for(int q=0;q<nrow;q++)
        {
            val_loc[i*3+q] = itred->second[q];
            //std::cout << "itred->second[q] " << itred->second[q] << std::endl;
        }
        
        i++;
    }
    
    int* offsets = distrimap->getOffsets();
    int* nlocs   = distrimap->getNlocs();
    
    
    MPI_Allgatherv(key_loc,
                   mapSizeLoc,
                   MPI_INT,
                   key_tot,
                   nlocs,
                   offsets,
                   MPI_INT, mpi_comm);
    
    int* offsetsVal = distrimapVal->getOffsets();
    int* nlocsVal   = distrimapVal->getNlocs();
    
    MPI_Allgatherv(val_loc,
                   mapSizeLoc*3,
                   MPI_DOUBLE,
                   val_tot,
                   nlocsVal,
                   offsetsVal,
                   MPI_DOUBLE, mpi_comm);
    
    int key,val;
    std::map<int,std::vector<double> > mappie_glob;
    for(int i=0;i<mapSizeTot;i++)
    {
        key = key_tot[i];
        
        std::vector<double> values(3);
        for(int q=0;q<3;q++)
        {
            values[q] = val_tot[i*3+q];
        }
        
        if(mappie_glob.find(key)==mappie_glob.end())
        {
            
            mappie_glob[key] = values;
            //std::cout << "itred->second[q] " << val[0] << " " << val[1] << " " << val[2] << std::endl;
        }
    }
    
    return mappie_glob;
}







//void OutputTetrahedralMeshOnPartition(TetrahedraMesh* tmesh, MPI_Comm comm)
//{
//
//    int world_size;
//    MPI_Comm_size(comm, &world_size);
//    // Get the rank of the process
//    int world_rank;
//    MPI_Comm_rank(comm, &world_rank);
//
//    std::vector<int> lverts;
//    std::map<int,int> lpartv2gv_v2;
//    std::map<int,int> gv2lpv2;
//
//    std::set<int> gv_set;
//    int lcv2 = 0;
//    Array<int>* ien_part_tetra     = tmesh->ien_part_tetra;
//    Array<int>* ien_part_hybrid    = tmesh->ien_part_hybrid;
//    std::vector<Vert*> locVs       = tmesh->LocalVerts;
//    int nElonRank = ien_part_tetra->getNrow();
//
//    Array<int>* locelem2locnode= new Array<int>(nElonRank,4);
//
//    std::vector<Vert*> printVs;
//
//    for(int i=0;i<ien_part_tetra->getNrow();i++)
//    {
//        for(int q=0;q<ien_part_tetra->getNcol();q++)
//        {
//            int gv = ien_part_tetra->getVal(i,q);
//            int lvv = tmesh->globV2locV[gv];
//
//            if(gv_set.find(gv)==gv_set.end())
//            {
//                gv_set.insert(gv);
//                lverts.push_back(lvv);
//                lpartv2gv_v2[lvv]=gv;
//                gv2lpv2[gv]=lcv2;
//                locelem2locnode->setVal(i,q,lcv2);
//
//                printVs.push_back(locVs[lvv]);
//
//                lcv2=lcv2+1;
//            }
//            else
//            {
//                int lcv_u = gv2lpv2[gv];
//                locelem2locnode->setVal(i,q,lcv_u);
//            }
//        }
//    }
//
//    std::vector<Vert*> lv = tmesh->LocalVerts;
//    std::string filename = "checkPart_" + std::to_string(world_rank) + ".dat";
//    std::ofstream myfile;
//    myfile.open(filename);
//    myfile << "TITLE=\"volume_part_"  + std::to_string(world_rank) +  ".tec\"" << std::endl;
//    myfile <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
//    myfile <<"ZONE N = " << printVs.size() << ", E = " << nElonRank << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;
//
//    for(int i=0;i<printVs.size();i++)
//    {
//        myfile << printVs[i][0] << " " << printVs[i][1] << " " << printVs[i][2] << std::endl;
//    }
//    int gv0,gv1,gv2,gv3,gv4,gv5,gv6,gv7;
//    int lv0,lv1,lv2,lv3,lv4,lv5,lv6,lv7;
//    for(int i=0;i<ien_part_hybrid->getNrow();i++)
//    {
//        myfile <<   locelem2locnode->getVal(i,0)+1 << "  " <<
//        locelem2locnode->getVal(i,1)+1 << "  " <<
//        locelem2locnode->getVal(i,2)+1 << "  " <<
//        locelem2locnode->getVal(i,3)+1 << "  " << std::endl;
//    }
//
//
//    myfile.close();
//}




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
    
    // US3D* us3d              = ReadUS3DData(fn_conn,fn_grid,fn_data,inputs->ReadFromStats,inputs->StateVar,comm,info);
    mesh* meshRead = ReadUS3DMeshData(fn_conn,fn_grid,fn_data,
                                        inputs->ReadFromStats,
                                        inputs->StateVar,
                                       comm,info);
    int Nel_loc = meshRead->ien.size();

    std::cout << "Nel_loc " << Nel_loc << std::endl;
    // for(int i=0;i<Nel_loc;i++)
    // {
    //     for(int j=0;j<meshRead->ien[0].size();j++)
    //     {
    //         std::cout << meshRead->ien[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    PrismTetraTrace* pttrace = new PrismTetraTrace(comm, meshRead->ife, meshRead->iet, 
                                            meshRead->nElem, meshRead->nFace, meshRead->nVert);

    std::map<int,std::vector<int> > trace = pttrace->GetTrace();

    std::cout << "trace " << trace.size() << std::endl;

    std::map<int,std::vector<int> >::iterator itmiv;

    std::map<int,std::vector<int> > tetras;
    std::map<int,std::vector<int> > prisms;
    std::map<int,std::vector<double> > tetras_data;

    std::map<int,int> loc2glob_prismv;
    std::map<int,int> glob2loc_prismv;

    std::map<int,int> loc2glob_tetrav;
    std::map<int,int> glob2loc_tetrav;

    int tetrav_loc = 0;
    int prismv_loc = 0;
    int tetra_id   = 0;
    int prism_id   = 0;
    int ntetra     = meshRead->ntetra;
    int* ien_tetra = new int[ntetra*4];

    for(itmiv=meshRead->ien.begin();itmiv!=meshRead->ien.end();itmiv++)
    {
        int elid   = itmiv->first;
        int eltype = meshRead->iet[elid];

        if(eltype == 2)
        {
            int nv = itmiv->second.size();
            std::vector<int> tetra(nv,0);
            
            for(int q=0;q<itmiv->second.size();q++)
            {
                ien_tetra[tetra_id*4+q] = itmiv->second[q];
            }
            
            tetras[elid] = itmiv->second;
            tetra_id++;
            
        }
        if(eltype == 6)
        {
            int nv = itmiv->second.size();
            std::vector<int> prism(nv,0);
            for(int i=0;i<nv;i++)
            {
                
                if(glob2loc_prismv.find(itmiv->second[i])==glob2loc_prismv.end())
                {
                    glob2loc_prismv[itmiv->second[i]]   = prismv_loc;
                    loc2glob_prismv[prismv_loc]         = itmiv->second[i];
                    prism[i]                            = prismv_loc;
                    prismv_loc++;
                }
                else
                {
                    int prismv_new                      = glob2loc_prismv[itmiv->second[i]];
                    prism[i]                            = prismv_new;
                }
            }

            prisms[prism_id] = prism;
            
            prism_id++;
        }
    }

    // we need to pass the number of verts per element in case the partition has no elements of this type.

    //RedistributeMeshtThroughRoot(tetras,4,comm);

    RepartitionObject* tetra_repart = new RepartitionObject(meshRead, 
                                                        tetras, 
                                                        trace, 
                                                        tetras_data, 
                                                        comm);

    RepartitionObject* prism_repart = new RepartitionObject(meshRead, 
                                                        prisms, 
                                                        trace, 
                                                        tetras_data, 
                                                        comm);

    

    /*
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

