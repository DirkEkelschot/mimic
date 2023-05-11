#include "src/adapt_recongrad.h"
#include "src/adapt_io.h"
#include "src/adapt_parops.h"
#include "src/adapt_output.h"
#include "src/adapt_boundary.h"
#include "src/adapt_distri_parstate.h"
#include "src/adapt_redistribute.h"
#include "src/adapt_DefinePrismMesh.h"
#include "src/adapt_prismaticlayer.h"
#include "src/NekFace.h"
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
//        myfile << printVs[i]->x << " " << printVs[i]->y << " " << printVs[i]->z << std::endl;
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
    
    US3D* us3d    = ReadUS3DData(fn_conn,fn_grid,fn_data,inputs->ReadFromStats,inputs->StateVar,comm,info);
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
    std::cout << "Start of partitioning..."<<std::endl;
    Partition* P = new Partition(us3d->ien, us3d->iee, us3d->ief, us3d->ie_Nv , us3d->ie_Nf,
                                 us3d->ifn, us3d->ife, us3d->if_ref, us3d->if_Nv,
                                 parmetis_pstate, ien_pstate, ife_pstate,
                                 us3d->xcn, xcn_pstate, Ui, us3d->ie_tetCnt, comm);
    
    std::cout << "End of partitioning..."<<std::endl;
    std::vector<int> LocElem                    = P->getLocElem();
    std::vector<double> Uvaria                  = P->getLocElemVaria();
    std::map<int,Array<double>*> Uvaria_map     = P->getLocAndAdjElemVaria();
    std::map<int,std::map<int,double> >  n2n    = P->getNode2NodeMap();

    std::map<int,Array<double>*> Uvaria_map2;
    double UvariaV              = 0.0;
    for(int i=0;i<LocElem.size();i++)
    {
        int gid   = LocElem[i];
        UvariaV   = Uvaria[i];

        Array<double>* Uarr = new Array<double>(1,1);
        Uarr->setVal(0,0,UvariaV);
        Uvaria_map2[gid] = Uarr;
    }
    
    P->AddStateVecForAdjacentElements(Uvaria_map2,1,comm);
    
//    std::map<int,Array<double>*> Uvaria_ref = ReadReferenceData2(world_rank);
//
//    std::map<int,Array<double>* >::iterator pltV;
//    int correctVal = 0;
//    int incorrectVal = 0;
//    int nothere = 0;
//    for(pltV=Uvaria_map.begin();pltV!=Uvaria_map.end();pltV++)
//    {
//        if(Uvaria_map2.find(pltV->first)!=Uvaria_map2.end())
//        {
//            double err = fabs(pltV->second->getVal(0,0)-Uvaria_map2[pltV->first]->getVal(0,0));
//            if(err > 1.0e-05)
//            {
//                std::cout << std::setprecision(16) << " Uvaria_map "<< pltV->first << " " << pltV->second->getVal(0,0) << " " << Uvaria_map2[pltV->first]->getVal(0,0)  << " " << err << std::endl;
//                incorrectVal++;
//            }
//            else
//            {
//                correctVal++;
//            }
//        }
//        else
//        {
//            nothere++;
//        }
//    }
//
//
//    std::cout << "correct vs incorrect = " << world_rank << " " << correctVal << " " << nothere << " " << incorrectVal << "( " << Uvaria_map2.size() << " " << Uvaria_map.size() << std::endl;
//
    //std::cout << "world rank " << world_rank <<  " " << Uvaria_map.size() << " " << Uvaria.size() << std::endl;
    
    std::map<int,Array<double>* > Mvar_vmap = P->ReduceStateVecToAllVertices_V2(Uvaria_map,1);

    
    std::map<int,Array<double>*> Utke_map;
    std::map<int,Array<double>* >  TKE_vmap;
    

    if(inputs->ReadFromStats==1)
    {
        Utke_map = P->PartitionAuxilaryData(TKEi, comm);
        P->AddStateVecForAdjacentElements(Utke_map,1,comm);
        TKE_vmap = P->ReduceStateVecToAllVertices_V2(Utke_map,1);
    }
    
    std::map<int,double> gbMap;
    std::vector<int> LocAndAdjElem = P->getLocAndAdjElem();
    i_part_map* iee_pmap            = P->getIEEpartmap();
    int Nel_glob  = us3d->ien->getNglob();
    
    for(int i=0;i<gB->getNrow();i++)
    {
        int eid = Nel_glob+i;
        gbMap[eid]=gB->getVal(i,0);
    }
    
    delete gB;
        
    delete us3d->ghost;
    delete us3d->ien;
    delete us3d->iee;
    delete us3d->ief;
    delete us3d->ie_Nv;
    delete us3d->ie_Nf;
//  delete us3d->ifn;
    delete us3d->ife;
    //delete us3d->if_ref;
    delete us3d->if_Nv;
//  delete us3d->xcn;
    
    std::map<int,Array<double>* > dudx_vmap;
    std::map<int,Array<double>* > dudy_vmap;
    std::map<int,Array<double>* > dudz_vmap;
    
    std::map<int,Array<double>* >::iterator itgg;
    std::map<int,Array<double>* > hess_vmap;
    
    Mesh_Topology* meshTopo      = new Mesh_Topology(P,comm);

    if(inputs->recursive == 0)
    {
        if(world_rank == 0)
        {
            std::cout << "We are running extended WLSqGR..." << std::endl;
        }
        
        std::map<int,Array<double>* > Hess_map = ComputedUdx_LSQ_HO_US3D(P,Uvaria_map,meshTopo,gbMap,comm);
        
        P->AddStateVecForAdjacentElements(Hess_map,9,comm);

        hess_vmap = P->ReduceStateVecToAllVertices_V2(Hess_map,9);
        
        for(itgg = hess_vmap.begin();itgg!=hess_vmap.end();itgg++)
        {
            Array<double>* dudx_v = new Array<double>(1,1);
            dudx_v->setVal(0,0,itgg->second->getVal(0,0));
            Array<double>* dudy_v = new Array<double>(1,1);
            dudy_v->setVal(0,0,itgg->second->getVal(1,0));
            Array<double>* dudz_v = new Array<double>(1,1);
            dudz_v->setVal(0,0,itgg->second->getVal(2,0));
            
            dudx_vmap[itgg->first] = dudx_v;
            dudy_vmap[itgg->first] = dudy_v;
            dudz_vmap[itgg->first] = dudz_v;
        }
        
        double po = 6.0;
        if(inputs->RunWakRefinement == 0)
        {
            ComputeMetric(P,comm,hess_vmap,1.0,po,inputs->recursive,inputs->extended,inputs->hmin,inputs->hmax,inputs->MetScale);

        }
        if(inputs->RunWakRefinement == 1)
        {
            ComputeMetricWithWake(P, comm, TKE_vmap, hess_vmap, 1.0, po, inputs->hwake, inputs->recursive,inputs->hmin,inputs->hmax,inputs->MetScale);
        }

        std::map<int,Array<double>* >::iterator itgg;
        
        
            
            
        for(itgg=Hess_map.begin();itgg!=Hess_map.end();itgg++)
        {
            delete itgg->second;
        }
        
        
    }
    

    if(inputs->recursive == 1)
    {
        
        
        std::map<int,double> Volumes_tmp = meshTopo->getVol();
        std::map<int,double> Volumes;
        std::map<int,double>::iterator itc;
        for(itc=Volumes_tmp.begin();itc!=Volumes_tmp.end();itc++)
        {
            int elid = itc->first;
            int vol  = itc->second;
            Volumes[elid] = vol;
        }
        
        std::map<int,Array<double>* > dUdXi;

        if(inputs->extended == 1)
        {
            if(world_rank == 0)
            {
                std::cout << "We are running WLSqGR recursively with an extended reconstruction scheme..." << std::endl;
            }
            
            dUdXi = ComputedUdx_LSQ_LS_US3D(P,Uvaria_map,meshTopo,gbMap,comm);
        }
        if(inputs->extended == 0)
        {
            if(world_rank == 0)
            {
                std::cout << "We are running WLSqGR recursively with an conventional reconstruction scheme..." << std::endl;
            }
            
            dUdXi = ComputedUdx_LSQ_US3D(P,Uvaria_map2,meshTopo,gbMap,comm);
//            dUdXi = ComputedUdx_LSQ_Vrt_US3D(P,Uvaria_map,Mvar_vmap,meshTopo,gbMap,comm);

        }
        
        P->AddStateVecForAdjacentElements(dUdXi,3,comm);
        
        std::map<int,Array<double>* >::iterator grit;
        std::map<int,Array<double>* > dUidxi_map;
        std::map<int,Array<double>* > dUidyi_map;
        std::map<int,Array<double>* > dUidzi_map;
        std::map<int,Array<double>* > dUidXi_map;
        for(grit=dUdXi.begin();grit!=dUdXi.end();grit++)
        {
            Array<double>* dUdx_E = new Array<double>(1,1);
            dUdx_E->setVal(0,0,grit->second->getVal(0,0));
            Array<double>* dUdy_E = new Array<double>(1,1);
            dUdy_E->setVal(0,0,grit->second->getVal(1,0));
            Array<double>* dUdz_E = new Array<double>(1,1);
            dUdz_E->setVal(0,0,grit->second->getVal(2,0));
            
            dUidxi_map[grit->first]=dUdx_E;
            dUidyi_map[grit->first]=dUdy_E;
            dUidzi_map[grit->first]=dUdz_E;
            delete grit->second;
        }
        
        dudx_vmap = P->ReduceStateVecToAllVertices_V2(dUidxi_map,1);
        dudy_vmap = P->ReduceStateVecToAllVertices_V2(dUidyi_map,1);
        dudz_vmap = P->ReduceStateVecToAllVertices_V2(dUidzi_map,1);
        
        std::map<int,Array<double>* > dU2dXi2;
        std::map<int,Array<double>* > dU2dYi2;
        std::map<int,Array<double>* > dU2dZi2;
        
        
        if(inputs->extended == 1)
        {
            dU2dXi2 = ComputedUdx_LSQ_LS_US3D(P,dUidxi_map,meshTopo,gbMap,comm);
            dU2dYi2 = ComputedUdx_LSQ_LS_US3D(P,dUidyi_map,meshTopo,gbMap,comm);
            dU2dZi2 = ComputedUdx_LSQ_LS_US3D(P,dUidzi_map,meshTopo,gbMap,comm);
        }
        if(inputs->extended == 0)
        {
            
            dU2dXi2 = ComputedUdx_LSQ_US3D(P,dUidxi_map,meshTopo,gbMap,comm);
            dU2dYi2 = ComputedUdx_LSQ_US3D(P,dUidyi_map,meshTopo,gbMap,comm);
            dU2dZi2 = ComputedUdx_LSQ_US3D(P,dUidzi_map,meshTopo,gbMap,comm);
            
//            dU2dXi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidxi_map,dudx_vmap,meshTopo,gbMap,comm);
//            dU2dYi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidyi_map,dudy_vmap,meshTopo,gbMap,comm);
//            dU2dZi2 = ComputedUdx_LSQ_Vrt_US3D(P,dUidzi_map,dudz_vmap,meshTopo,gbMap,comm);
        }
        
        
        std::map<int,Array<double>* > Hess_map;
        
        for(itgg=dU2dXi2.begin();itgg!=dU2dXi2.end();itgg++)
        {
            int gid = itgg->first;
            
            Array<double>* Hess = new Array<double>(6,1);
            
            Hess->setVal(0,0,dU2dXi2[gid]->getVal(0,0));
            Hess->setVal(1,0,dU2dXi2[gid]->getVal(1,0));
            Hess->setVal(2,0,dU2dXi2[gid]->getVal(2,0));

            Hess->setVal(3,0,dU2dYi2[gid]->getVal(1,0));
            Hess->setVal(4,0,dU2dYi2[gid]->getVal(2,0));
            Hess->setVal(5,0,dU2dZi2[gid]->getVal(2,0));
            
            Hess_map[gid] = Hess;
            
            delete dU2dXi2[gid];
            delete dU2dYi2[gid];
            delete dU2dZi2[gid];
            
            t++;
        }
        
        for(grit=dUidxi_map.begin();grit!=dUidxi_map.end();grit++)
        {
            delete grit->second;
            delete dUidyi_map[grit->first];
            delete dUidzi_map[grit->first];
        }
        
        P->AddStateVecForAdjacentElements(Hess_map,6,comm);
        hess_vmap = P->ReduceStateVecToAllVertices_V2(Hess_map,6);
        
        double po = 6.0;
        
        if(inputs->RunWakRefinement == 0)
        {
            ComputeMetric(P,comm,hess_vmap,1.0,po,inputs->recursive,inputs->extended,inputs->hmin,inputs->hmax,inputs->MetScale);
        }
        if(inputs->RunWakRefinement == 1)
        {
            ComputeMetricWithWake(P, comm, TKE_vmap, hess_vmap, 1.0, po, inputs->hwake, inputs->recursive,inputs->hmin,inputs->hmax,inputs->MetScale);
        }
                   
        for(itgg=Hess_map.begin();itgg!=Hess_map.end();itgg++)
        {
            delete itgg->second;
        }
    }
    
    std::map<int,Array<double>* > hess_vmap_new;// = hess_vmap;
    double m00,m01,m02,m10,m11,m12,m20,m21,m22;
    double sum_dist;
    int vid;
    double di;
    std::map<int,Array<double>* >::iterator hessit;
    std::map<int,double>::iterator itn2n;
    int howmany = 0;
    //std::cout << "hess_vmap size =  " << hess_vmap.size() << std::endl; 
    for(hessit=hess_vmap.begin();hessit!=hess_vmap.end();hessit++)
    {
        int gvid       = hessit->first;

        //sum        = 0.0;
        //du         = 0.0;

        std::map<int,double>::iterator n2dit;

        if(n2n.find(gvid)!=n2n.end())
        {
            sum_dist   = 0.0;
            m00 = 0.0;m01 = 0.0;m02 = 0.0;
            m12 = 0.0;m11 = 0.0;m12 = 0.0;
            m20 = 0.0;m12 = 0.0;m22 = 0.0;

            std::map<int,double> n2nmap = n2n[gvid];
            //std::cout << "n2nmap.size " << n2nmap.size() << std::endl;
            for(itn2n=n2nmap.begin();itn2n!=n2nmap.end();itn2n++)
            {
                vid       = itn2n->first;
                di        = itn2n->second;

                if(hess_vmap.find(vid)!=hess_vmap.end())
                {
                    sum_dist  = sum_dist+di;

                    m00       = m00+hess_vmap[vid]->getVal(0,0)*di;
                    m01       = m01+hess_vmap[vid]->getVal(0,1)*di;
                    m02       = m02+hess_vmap[vid]->getVal(0,2)*di;
                    m10       = m10+hess_vmap[vid]->getVal(1,0)*di;
                    m11       = m11+hess_vmap[vid]->getVal(1,1)*di;
                    m12       = m12+hess_vmap[vid]->getVal(1,2)*di;
                    m20       = m20+hess_vmap[vid]->getVal(2,0)*di;
                    m21       = m21+hess_vmap[vid]->getVal(2,1)*di;
                    m22       = m22+hess_vmap[vid]->getVal(2,2)*di;
                }
            }
	    //std::cout << "di " << gvid << " " << sum_dist << " " << n2nmap.size() << " --->" << m00 << " " << m01 << " " << m02 << " " << m11 << " " << m12 << " " << m22 << std::endl;
            Array<double>* HabsNew = new Array<double>(3,3);

            HabsNew->setVal(0,0,m00/sum_dist);
            HabsNew->setVal(0,1,m01/sum_dist);
            HabsNew->setVal(0,2,m02/sum_dist);
            HabsNew->setVal(1,0,m10/sum_dist);
            HabsNew->setVal(1,1,m11/sum_dist);
            HabsNew->setVal(1,2,m12/sum_dist);
            HabsNew->setVal(2,0,m20/sum_dist);
            HabsNew->setVal(2,1,m21/sum_dist);
            HabsNew->setVal(2,2,m22/sum_dist);

            hess_vmap_new[gvid] = HabsNew;
            //hess_vmap_new[gvid] = hess_vmap[gvid];
        }
        else
        {
            hess_vmap_new[gvid] = hess_vmap[gvid];
            howmany++;
        }
    }
    
    

    //std::cout << "WorldRank " << world_rank << " " << howmany << " " << hess_vmap_new.size() << std::endl;


//    for(itgg=hess_vmap.begin();itgg!=hess_vmap.end();itgg++)
//    {
//        delete itgg->second;
//    }
    
    if(world_rank == 0)
    {
        std::cout << " Done computing the metric " << std::endl;
    }
    //std::cout << "Mvar_vmap " << Mvar_vmap.size() << std::endl;
    
    //std::vector<int> LocElem     = P->getLocElem();
    //================================================================================
    //================================================================================
    //================================================================================
    
    
    Domain* pDom = P->getPartitionDomain();
    
    std::map<int,std::vector<int> > tetraLoc   = pDom->Tetras;
    std::map<int,std::vector<int> > tetraEl    = pDom->GTetras;
    std::map<int,std::vector<int> > prismEl    = pDom->GPrisms;
    std::map<int,std::vector<int> > pyramidEl  = pDom->GPyramids;
    std::map<int,std::vector<int> > ushell_o   = pDom->ushell;
    std::map<int,Vert* > ushell_cen            = pDom->ushell_centroid;
    i_part_map* ief_part_map                   = P->getIEFpartmap();
    i_part_map* if_Nv_part_map                 = P->getIF_Nvpartmap();
    i_part_map* ifn_part_map                   = P->getIFNpartmap();
    i_part_map* ife_part_map                   = P->getIFEpartmap();
    i_part_map* if_ref_part_map                = P->getIFREFpartmap();
    i_part_map* if_Erank_part_map			   = P->getIFERankpartmap();
    std::map<int,int> tag2locV_map      	   = P->getGlobalVert2LocalVert();
    std::vector<Vert*> LocVerts_part     	   = P->getLocalVerts();
    std::map<int,int> lpartv2gv                = pDom->lpartv2gv;
    std::vector<int> loc_part_verts            = pDom->loc_part_verts;
    //================================================================================
    //================================================================================
    //================================================================================
    // Copy maps relevant maps so that we can destroy Partition*
    std::map<int,std::vector<int> >::iterator itmm;
    std::map<int,int >::iterator itr;
    std::map<int,Vert* >::iterator itrV;

    // ==== Copy required facemaps ====; 
    std::map<int,std::vector<int> > ifn_map;
    std::map<int,std::vector<int> > ife_map;
    std::map<int,std::vector<int> > iferank_map;
    std::map<int,int> ifref_map;
    std::map<int,int> if_Nv_map;
    
    for(itmm=ifn_part_map->i_map.begin();itmm!=ifn_part_map->i_map.end();itmm++)
    {
        int faceid        = itmm->first;
        int Nfn           = itmm->second.size();
        int fref          = if_ref_part_map->i_map[faceid][0];
        int fNv			  = if_Nv_part_map->i_map[faceid][0];
        
        std::vector<int>fn(Nfn);
        std::vector<int>fe(2);
        std::vector<int>ferank(2);

        for(int j=0;j<Nfn;j++)
        {
            fn[j]=ifn_part_map->i_map[faceid][j];
            if(j<2)
            {
                fe[j]     = ife_part_map->i_map[faceid][j];
                ferank[j] = if_Erank_part_map->i_map[faceid][j];
            }
        }
        ifn_map[faceid]     = fn;
        ife_map[faceid]     = fe;
        iferank_map[faceid] = ferank;
        ifref_map[faceid]   = fref;
        if_Nv_map[faceid]   = fNv;
    }
    // ==== Copy required elementmaps ====
    std::map<int,std::vector<int> > ief_map;
    
    for(itmm=ief_part_map->i_map.begin();itmm!=ief_part_map->i_map.end();itmm++)
	{
		int elemid        = itmm->first;
		int Nf            = itmm->second.size();

		std::vector<int>ef(Nf);

		for(int j=0;j<Nf;j++)
		{
			ef[j]=ief_part_map->i_map[elemid][j];
		}
		
		ief_map[elemid] 	= ef;
		
	}
    
    // ==== Copy required element2node maps for tetrahedra ====
    std::map<int,std::vector<int> > tetrahedra;
	for(itmm=tetraEl.begin();itmm!=tetraEl.end();itmm++)
	{
		int elemid        = itmm->first;
		int nn            = itmm->second.size();
		std::vector<int> tetraNodes(nn);
		for(int q=0;q<nn;q++)
		{
			tetraNodes[q] = itmm->second[q];
		}
		tetrahedra[elemid] = tetraNodes;
	}
	
    // ==== Copy required element2node maps for prisms ====
	std::map<int,std::vector<int> > prisms;
	for(itmm=prismEl.begin();itmm!=prismEl.end();itmm++)
	{
		int elemid        = itmm->first;
		int nn            = itmm->second.size();
		std::vector<int> prismNodes(nn);
		for(int q=0;q<nn;q++)
		{
			prismNodes[q] = itmm->second[q];
		}
		prisms[elemid] = prismNodes;
	}
	
	std::map<int,std::vector<int> > pyramids;
	for(itmm=pyramidEl.begin();itmm!=pyramidEl.end();itmm++)
	{
		int elemid        = itmm->first;
		int nn            = itmm->second.size();
		std::vector<int> pyramidNodes(nn);
		for(int q=0;q<nn;q++)
		{
			pyramidNodes[q] = itmm->second[q];
		}
		pyramids[elemid] = pyramidNodes;
	}
	
    // ==== Copy required element2node maps for prisms ====
	std::map<int,std::vector<int> > ushell;
	for(itmm=ushell_o.begin();itmm!=ushell_o.end();itmm++)
	{
		int elemid        = itmm->first;
		int nn            = itmm->second.size();
		std::vector<int> shellnodes(nn);
		for(int q=0;q<nn;q++)
		{
			shellnodes[q] = itmm->second[q];
		}
		ushell[elemid] = shellnodes;
	}
		
	
    std::map<int,int> tag2locVrtMap;
    for(itr=tag2locV_map.begin();itr!=tag2locV_map.end();itr++)
    {
    	int key = itr->first;
    	int val = itr->second;
    	
    	tag2locVrtMap[key] = val;
    }
    
    
    std::vector<Vert*> LocVerts(LocVerts_part.size());
    for(int q=0;q<LocVerts_part.size();q++)
    {
        Vert* val = new Vert;
        val->x = LocVerts_part[q]->x;
        val->y = LocVerts_part[q]->y;
        val->z = LocVerts_part[q]->z;
        
        LocVerts[q] = val;
    }

    //std::cout << "Copied all relevant data over " << std::endl;
	
    //====================================================================

    if(tetraLoc.size()!=0 && debug == 1)
    {
        std::ofstream myfilet;
        myfilet.open("metricFieldGrad_" + std::to_string(world_rank) + ".dat");
        myfilet << "TITLE=\"new_volume.tec\"" << std::endl;
        myfilet <<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"G0\", \"G1\", \"G2\", \"M00\", \"M01\", \"M02\", \"M11\", \"M12\", \"M22\"" << std::endl;
        //myfilet <<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\"" << std::endl;
        myfilet <<"ZONE N = " << loc_part_verts.size() << ", E = " << tetraLoc.size() << ", DATAPACKING = POINT, ZONETYPE = FETETRAHEDRON" << std::endl;

        for(int i=0;i<loc_part_verts.size();i++)
        {
            int loc_vid  = loc_part_verts[i];
            int glob_vid = lpartv2gv[loc_vid];
            myfilet << LocVerts[loc_vid]->x << " " <<
                       LocVerts[loc_vid]->y << " " <<
                       LocVerts[loc_vid]->z << " " <<
                       Mvar_vmap[glob_vid]->getVal(0,0)  << " " <<
                       dudx_vmap[glob_vid]->getVal(0,0) << " " <<
                       dudy_vmap[glob_vid]->getVal(0,0) << " " <<
                       dudz_vmap[glob_vid]->getVal(0,0) << " " <<
                       hess_vmap_new[glob_vid]->getVal(0,0) << " " <<
                       hess_vmap_new[glob_vid]->getVal(0,1) << " " <<
                       hess_vmap_new[glob_vid]->getVal(0,2) << " " <<
                       hess_vmap_new[glob_vid]->getVal(1,1) << " " <<
                       hess_vmap_new[glob_vid]->getVal(1,2) << " " <<
                       hess_vmap_new[glob_vid]->getVal(2,2) << std::endl;
        }

        std::map<int,std::vector<int> >::iterator itmm2;
        for(itmm2=tetraLoc.begin();itmm2!=tetraLoc.end();itmm2++)
        {
            myfilet << itmm2->second[0]+1 << " " << itmm2->second[1]+1 << " "
                    << itmm2->second[2]+1 << " " << itmm2->second[3]+1 <<  std::endl;
        }

        myfilet.close();
        
        delete P;
    }
    else
    {
        delete P;
    }
   //====================================================================
    
    
    clock_t t0_redis = clock();

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
    std::vector<Vert*> locVs                        = tetra_distri->GetLocalVertices();
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
    std::map<int,Vert*> shellVertCoord2Ref          = tetra_distri->GetShellVertCoords2RefMap_Global();
    std::map<int,int> shellvertTag2ref;
    std::map<int,int> shellvertTag2ref2;

    

    std::map<int,std::vector<int> >::iterator bndtest;
    
    std::map<int,int>::iterator itt;
    for(itt=shellvert2ref.begin();itt!=shellvert2ref.end();itt++)
    {
    	int vtet = itt->first;
    	int vhyb = tetV2tagV[vtet];
    	int ref  = itt->second;
    	shellvertTag2ref[vhyb] = ref;
    }
    
    //int icomm;
    //Based on the new local tetrahedra mesh, we output a tecplot file per processor that has the geometry of the computational domain that is owned by world_rank.
    
    //OutputTetrahedralMeshOnPartition(tmesh,comm);
    
    Array<int>* ien_part_tetra   = element2node;
    int nTetrahedra              = element2node->getNrow();
    int nVertices                = locVs.size();
    int nTriangles               = nFaces;
    int nEdges                   = 0;
    int nPrisms                  = 0;
    int nQuadrilaterals          = 0;
    
    PrismaticLayer* prsmLyr = new PrismaticLayer(prisms, ief_map, ifn_map, ife_map, ifref_map,
                                                 if_Nv_map, iferank_map, ushell, tag2locVrtMap,
												 LocVerts, shellvertOriginalTag2ref_Glob, comm);
    
    nPrisms                                                 = prisms.size();
    DistributedParallelState* distPrismIN                   = new DistributedParallelState(nPrisms,comm);
    int nPrismsTot                                          = distPrismIN->getNel();
    
    std::map<int,int> tagE2gE                               = prsmLyr->getTag2GlobalElementMap();
    std::map<int,int> gE2tagE                               = prsmLyr->getGlobal2TagElementMap();
    std::map<int,int> rhp                                   = prsmLyr->getRightElementGlobalIDForFaces();
    std::map<int,int> lhp                                   = prsmLyr->getLeftElementGlobalIDForFaces();
    std::map<int,std::vector<int> > pbcmap                  = prsmLyr->getBoundaryCondition2FaceID();
    std::map<int,int> tag2element_shell                     = prsmLyr->getTag2Element4TetPrismInterface();
    std::map<int,std::vector<int> > shared_face2node_prism  = prsmLyr->getOwnedSharedFace2NodeMap();
    std::map<int,std::vector<int> > int_face2node_prism     = prsmLyr->getInternalFace2NodeMap();
    std::map<int,std::vector<int> > bc_face2node_prism      = prsmLyr->getBoundaryFace2NodeMap();
    std::map<int,int> tag2glob_prism                        = prsmLyr->getVertexTag2GlobalMap();
    Array<double>* xcn_prisms_int                           = prsmLyr->getInternalCoordinates();
    Array<double>* xcn_prisms_shared                        = prsmLyr->getSharedCoordinates();
    std::map<int,int> SharedVertsNotOwned                   = prsmLyr->getNotOwnedSharedVerticesMap();
    Array<int>* parmmg_iet_prisms                           = prsmLyr->getElementType();
    //std::map<int,int> sharedVmap                            = prsmLyr->getSharedVertexMap();

    
    DistributedParallelState* distPrismIntVerts = new DistributedParallelState(xcn_prisms_int->getNrow(),comm);
    DistributedParallelState* distPrismShaVerts = new DistributedParallelState(xcn_prisms_shared->getNrow(),comm);

    int foundU = 0;
    
    //std::map<std::set<int>, int > vertref2shell_prism;
    FaceSetPointer m_PMMG_RefsOnShell_2_Prism;
    int nshell = 0;
    for(itmm=ushell.begin();itmm!=ushell.end();itmm++)
    {
        int fhyb      = itmm->first;
        int Etettag   = itmm->second[0];
        int Eprismtag = itmm->second[1];
        int EprismNew = tag2element_shell[Eprismtag];
    
        if(shellface2vertref.find(fhyb)!=shellface2vertref.end())
        {
            std::set<int>::iterator its;
            std::vector<int> refs(shellface2vertref[fhyb].size());
            int c = 0;
            for(its=shellface2vertref[fhyb].begin();
                its!=shellface2vertref[fhyb].end();its++)
            {
                refs[c] = *its;
                c++;
            }
            FaceSharedPtr RefFacePointer = std::shared_ptr<NekFace>(new NekFace(refs));
            pair<FaceSetPointer::iterator, bool> testInsPointer;
            testInsPointer = m_PMMG_RefsOnShell_2_Prism.insert(RefFacePointer);
            
            if(testInsPointer.second)
            {
                (*testInsPointer.first)->SetFaceLeftElement(EprismNew);
            }
            
            //vertref2shell_prism[shellface2vertref[fhyb]] = EprismNew;
            foundU++;
        }
        else
        {
           std::cout << "Fhyb NotFound-> " << fhyb << " " << world_rank << std::endl;
        }
        
        nshell++;
    }
    
    
    clock_t t1_redis = clock();
    double duration_redis = ( t1_redis - t0_redis) / (double) CLOCKS_PER_SEC;
    double dur_max_redis;
    MPI_Allreduce(&duration_redis, &dur_max_redis, 1, MPI_DOUBLE, MPI_MAX, comm);
    
    
    //std::cout << "WOR " << world_rank << " " << vertref2shell_prism.size() << " " << ushell.size() << " shellface2vertref " << shellface2vertref.size() << " " << foundU <<std::endl;
    
    //===============================================================================================
    //===============================================================================================
    //===============================================================================================
    
    
    
    clock_t t0_input = clock();

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
    for ( k=0; k<nVertices; ++k )
    {
        double vx = locVs[k]->x;
        double vy = locVs[k]->y;
        double vz = locVs[k]->z;
        
        int vert = locV2globV[k];

        int vref = 86;

        if ( PMMG_Set_vertex(parmesh,vx,vy,vz, vref, k+1) != 1 )
        {
        MPI_Finalize();
        exit(EXIT_FAILURE);
        }

        double* tensor = new double[6];


        tensor[0] = metric[vert]->getVal(0,0);
        tensor[1] = metric[vert]->getVal(1,0);
        tensor[2] = metric[vert]->getVal(2,0);
        tensor[3] = metric[vert]->getVal(3,0);
        tensor[4] = metric[vert]->getVal(4,0);
        tensor[5] = metric[vert]->getVal(5,0);

        if(PMMG_Set_tensorMet(parmesh,tensor[0],tensor[1],tensor[2],tensor[3],tensor[4],tensor[5],k+1)!=1)
        {
         MPI_Finalize();
         exit(EXIT_FAILURE);
        }
        delete[] tensor;
        //delete metric[vert];
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
    std::vector<Vert*> unshellVin;
    std::vector<std::vector<int> > unshellTin;
    int outflowFound = 0;
    int inflowFound = 0;
    int shellFound = 0;
    for ( k=0; k<nTriangles; ++k )
    {
        int faceID  = faces4parmmg[k];
        int facetag = tetF2hybF[faceID];
    
        v0      = face2node[faceID][0];
        v1      = face2node[faceID][1];
        v2      = face2node[faceID][2];
        
        v0l     = globV2locV[v0];
        v1l     = globV2locV[v1];
        v2l     = globV2locV[v2];
        
        if(face2ref.find(faceID)!=face2ref.end())
        {
            refer   = face2ref[faceID];
        }
        else
        {
            refer   = 0;
        }
        
        if(refer == 36)
        {
            outflowFound++;
        }
        if(refer == 10)
        {
            inflowFound++;
        }
        if(refer == 13)
        {
            shellFound++;
        }
        
        
        if ( PMMG_Set_triangle(parmesh,v0l+1,v1l+1,v2l+1,refer,k+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        
        
        
        
        if(refer == 13)
        {
            //flippie = -1;
            double testerr;
            double testerrStored;
            if(shell_tet2hybF.find(faceID)!=shell_tet2hybF.end())
            {
                if(shell_tet2hybF[faceID]!=facetag)
                {
                    std::cout << "Error:: Tet2Hyb maps are incorrect." << std::endl;
                }
                else
                {
                    std::set<int>::iterator its;
                    std::set<int> sf = shellface2vertref[facetag];
                    int inde = -1;
                    std::set<int> sft_t;
                    //flippie = -1;
                    std::vector<int> row(3);
                    for(int u=0;u<3;u++)
                    {
                        
                        int vid = face2node[faceID][u];
                        row[u] = vid;
                        int vidg = globV2locV[vid];
                        //double* c0 = new double[3];
                        std::vector<double> c0(3);
                        c0[0] = locVs[vidg]->x;
                        c0[1] = locVs[vidg]->y;
                        c0[2] = locVs[vidg]->z;
                        vertref = 86;
                        if(shellvert2ref.find(vid)!=shellvert2ref.end())
                        {
                            vertref = shellvert2ref[vid];
                        }
                        //vertref = shellvert2ref[vid];
                        
                        
                        if(shell_g2l.find(vid)==shell_g2l.end())
                        {
                            Vert* vsh = new Vert;
                            vsh->x = c0[0];
                            vsh->y = c0[1];
                            vsh->z = c0[2];
                            unshellVin.push_back(vsh);
                            shell_g2l[vid] = locs;
                            locs++;
                        }
                        
                        if ( PMMG_Set_vertex(parmesh, c0[0], c0[1], c0[2], vertref, vidg+1) != 1 )
                        {
                        MPI_Finalize();
                        exit(EXIT_FAILURE);
                        }
                        
                        
                        //PMMG_Set_vertex( parmesh, c0[0], c0[1], c0[2], vertref, vidg+1 );

                        sft_t.insert(vertref);
                        
                    }
                    
                    unshellTin.push_back(row);
                    
                    sft_t.clear();
                }
                suc++;
            }
            PMMG_Set_requiredTriangle( parmesh, k+1 );
        }
    }
    
    int outflowFound_sum = 0;
    int inflowFound_sum = 0;
    int shellFound_sum = 0;
    MPI_Allreduce(&outflowFound, &outflowFound_sum, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&inflowFound, &inflowFound_sum, 1, MPI_INT, MPI_SUM, comm);
    MPI_Allreduce(&shellFound, &shellFound_sum, 1, MPI_INT, MPI_SUM, comm);
    if(world_rank == 0)
    {
        std::cout << " Outflow: " << outflowFound_sum << " Inflow: " << inflowFound_sum << " Shell: " << shellFound_sum << std::endl;

    }
    
    if(unshellTin.size()!=0 && debug == 1)
    {
        ofstream myfile_INP;

        string filename_shellINPUT = "inputShell_" + std::to_string(world_rank) + ".dat";

        myfile_INP.open(filename_shellINPUT);
        myfile_INP << "TITLE=\"inputShell.tec\"" << std::endl;
        myfile_INP <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
        //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
        myfile_INP <<"ZONE N = " << unshellVin.size() << ", E = " << unshellTin.size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
        for(int i=0;i<unshellVin.size();i++)
        {
          myfile_INP << unshellVin[i]->x << "   " << unshellVin[i]->y << "   " << unshellVin[i]->z << std::endl;
        }

        for(int i=0;i<unshellTin.size();i++)
        {
            int g0=unshellTin[i][0];
            int g1=unshellTin[i][1];
            int g2=unshellTin[i][2];
            
            int l0=shell_g2l[g0];
            int l1=shell_g2l[g1];
            int l2=shell_g2l[g2];
            
          myfile_INP << l0+1 << "    " << l1+1 << "   " << l2+1 << std::endl;
        }
        myfile_INP.close();
    }
    
//    DistributedParallelState* sucDist = new DistributedParallelState(suc,comm);
    
//    double* vertIN = new double[nVertices*3];
//    int* refIN = new int[nVertices];
//    int *requiredIN = (int*)calloc(MAX4(nVertices,nTetrahedra,nTriangles,nEdges),sizeof(int));
//
//    int *cornerIN = (int*)calloc(nVertices,sizeof(int));
//    int foundRrefIN    = 0;
//    std::vector<int> rfIN;
//    std::map<int,int> frfIN;
//    std::map<int,std::vector<double> > ref2coordinates;
//    for ( k=0; k<nVertices; k++ )
//    {
//          int pos = 3*k;
//          if ( PMMG_Get_vertex(parmesh,&(vertIN[pos]),&(vertIN[pos+1]),&(vertIN[pos+2]),
//                               &(refIN[k]),&(cornerIN[k]),&(requiredIN[k])) != 1 ) {
//            fprintf(inm,"Unable to get mesh vertex %d \n",k);
//            ier = PMMG_STRONGFAILURE;
//          }
//        int refi = refIN[k];
//        if(refi != -20 && refi != 0)
//        {
//            std::vector<double> coords(3);
//            coords[0] = vertIN[pos];
//            coords[1] = vertIN[pos+1];
//            coords[2] = vertIN[pos+2];
//            ref2coordinates[refi] = coords;
//        }
//    }
        
    
    //std::map<int,std::vector<double> > ref2coordsAll = AllGatherMapDoubleVec(ref2coordinates,comm);

    
    
    int API_mode = 0;
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_APImode, API_mode ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    ier = PMMG_Set_numberOfFaceCommunicators(parmesh, ncomm);
    
    for(int icomm=0; icomm<ncomm; icomm++ )
    {
      // Set nb. of entities on interface and rank of the outward proc
     
      ier = PMMG_Set_ithFaceCommunicatorSize(parmesh, icomm,
                                             color_face[icomm],
                                             ntifc[icomm]);

    //Set local and global index for each entity on the interface
      ier = PMMG_Set_ithFaceCommunicator_faces(parmesh, icomm,
                                               ifc_tria_loc[icomm],
                                               ifc_tria_glob[icomm], 1);
    }
    
    
    std::map<int,std::vector<int> >::iterator ittet;
    k = 0;
    for ( int t = 0;t < nTetrahedra; t++  )
    {
        v0 = ien_part_tetra->getVal(t,0);
        v1 = ien_part_tetra->getVal(t,1);
        v2 = ien_part_tetra->getVal(t,2);
        v3 = ien_part_tetra->getVal(t,3);
        
        v0l = globV2locV[v0];
        v1l = globV2locV[v1];
        v2l = globV2locV[v2];
        v3l = globV2locV[v3];
        
        double* P = new double[4*3];
        
        P[0*3+0]=locVs[v0l]->x;   P[0*3+1]=locVs[v0l]->y;    P[0*3+2]=locVs[v0l]->z;
        P[1*3+0]=locVs[v1l]->x;   P[1*3+1]=locVs[v1l]->y;    P[1*3+2]=locVs[v1l]->z;
        P[2*3+0]=locVs[v2l]->x;   P[2*3+1]=locVs[v2l]->y;    P[2*3+2]=locVs[v2l]->z;
        P[3*3+0]=locVs[v3l]->x;   P[3*3+1]=locVs[v3l]->y;    P[3*3+2]=locVs[v3l]->z;

        double Vtet = GetQualityTetrahedra(P);
        
        if(Vtet<0.0)
        {
            std::cout << " negative volume in Element " << t << " on rank " << world_rank  <<std::endl;
        }
        if ( PMMG_Set_tetrahedron(parmesh,v0l+1,v1l+1,v2l+1,v3l+1,1.0,t+1) != 1 )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        
        delete[] P;
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
    
    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_globalNum, 1 ) ) {
      MPI_Finalize();
      exit(EXIT_FAILURE);
    };
    
    
//    if( !PMMG_Set_iparameter( parmesh, PMMG_IPARAM_mem, 2000 ) ) {
//      MPI_Finalize();
//      exit(EXIT_FAILURE);
//    };
    
    
    
    
//    //======================================================================
//    //======================================================================
//    //======================================================================
//
//    // Destructor of redistribute functions;
//
//    delete element2node;
//    std::map<int,Array<double>* >::iterator itdes;
//    for(itdes = metric.begin();itdes != metric.end(); itdes++ )
//    {
//        delete itdes->second;
//    }
//
//    for(int q=0;q<locVs.size();q++)
//    {
//        delete locVs[q];
//    }
//
//    std::map<int,int*>::iterator itpdes;
//    for(itpdes = face2node.begin();itpdes != face2node.end(); itpdes++ )
//    {
//        delete[] itpdes->second;
//    }
////    delete[] color_face;
////    delete[] ntifc;
////
////
////    shellface2vertref.clear();
//    bndref2face.clear(); //not necessary
//    face2element.clear(); //not necessary
//    faces4parmmg.clear();
//    globV2locV.clear();
//    locV2globV.clear();
//    locShF2globShF.clear();
//    face2ref.clear();
//    shell_tet2hybF.clear();
//    shellvert2ref.clear();
//    shellface2vertref.clear();
//    shellvert2ref_local.clear();
//    tetF2hybF.clear();
//    tetV2tagV.clear();
////    shellvertOriginalTag2ref_Glob.clear();
//    bndref2face.clear();
//
//    //delete prsmLyr;
//    //======================================================================
//    //======================================================================
//    //======================================================================

    DistributedParallelState* distTetraB = new DistributedParallelState(nTetrahedra,comm);
    DistributedParallelState* distPrismB = new DistributedParallelState(nPrisms,comm);
    int ToTElements_prism_bef            = distPrismB->getNel();
    int ToTElements_bef                  = distTetraB->getNel();
    
    if(world_rank == 0)
    {
        std::cout << "Total elements = " << ToTElements_prism_bef << " " << ToTElements_bef << std::endl;
    }
    //std::cout << "number of tets per rank  " << world_rank << " ("<<nTetrahedra<<", " << ToTElements_bef << ") (" << nPrisms <<", " << ToTElements_prism_bef <<")" << std::endl;
    
    clock_t t1_input = clock();
    double duration_input = ( t1_input - t0_input) / (double) CLOCKS_PER_SEC;
    double dur_max_input;
    MPI_Allreduce(&duration_input, &dur_max_input, 1, MPI_DOUBLE, MPI_MAX, comm);
        
    
    
    delete tetra_distri;
    
    clock_t t0_adapt = clock();

    int ierlib = PMMG_parmmglib_distributed( parmesh );
   
    clock_t t1_adapt = clock();
    double duration_adapt = ( t1_adapt - t0_adapt) / (double) CLOCKS_PER_SEC;
    double dur_max_adapt;
    MPI_Allreduce(&duration_adapt, &dur_max_adapt, 1, MPI_DOUBLE, MPI_MAX, comm);
    
    
    clock_t t0_ijk = clock();
    
    if(ierlib==0 && world_rank == 0)
    {
        std::cout << "SUCCESFULLY adapted the mesh in parallel!" << std::endl;
    }
    else if(ierlib==1 && world_rank == 0)
    {
        std::cout << "FAILED to adapt the mesh in parallel! "<< std::endl;
    }
    
    
    
    int nVerticesOUT   = 0;
    int nTetrahedraOUT = 0;
    int nTrianglesOUT  = 0;
    int nEdgesOUT      = 0;
    
    if ( PMMG_Get_meshSize(parmesh,&nVerticesOUT,&nTetrahedraOUT,NULL,&nTrianglesOUT,NULL,
                           &nEdgesOUT) !=1 )
    {
        ier = PMMG_STRONGFAILURE;
    }
    
    int             nodeGloNumber,nodeOwner;
    std::map<int,int> loc2globVid;
    std::map<int,int> glob2locVid;
    std::vector<int> globIDs;
    for( k = 1; k <= nVerticesOUT; k++ )
    {
        if( !PMMG_Get_vertexGloNum( parmesh, &nodeGloNumber, &nodeOwner ) )
        {
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
    
        loc2globVid[k]=nodeGloNumber;
        glob2locVid[nodeGloNumber]=k;
        globIDs.push_back(nodeGloNumber);
    }
    
    int *required = (int*)calloc(MAX4(nVerticesOUT,
                                      nTetrahedraOUT,
                                      nTrianglesOUT,
                                      nEdgesOUT),sizeof(int));
    
    int *refOUT = (int*)calloc(MAX4(nVerticesOUT,
                                    nTetrahedraOUT,
                                    nTrianglesOUT,
                                    nEdgesOUT),sizeof(int));
    
    int *corner = (int*)calloc(nVerticesOUT,sizeof(int));
    int pos;
    //int *refOUT = new int[nVerticesOUT];

    std::vector<std::vector<int> > outT;
    
    double *vertOUT = (double*)calloc((nVerticesOUT)*3,sizeof(double));
    //std::map<int,std::vector<double> > ref2coordinatesOUT;
    for ( k=0; k<nVerticesOUT; k++ )
    {
          pos = 3*k;
          if ( PMMG_Get_vertex(parmesh,&(vertOUT[pos]),&(vertOUT[pos+1]),&(vertOUT[pos+2]),
                               &(refOUT[k]),&(corner[k]),&(required[k])) != 1 ) {
            fprintf(inm,"Unable to get mesh vertex %d \n",k);
            ier = PMMG_STRONGFAILURE;
          }
    }
    

    
    if(inputs->niter==0)
    {
        std::cout << "Check the input and outputted shard faces." << std::endl;
        
        int **out_tria_loc;
        int *nitem_face_comm;
        int next_face_comm;
        ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
        int *color_node_out,*color_face_out;
        color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
        nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                 &color_face_out[icomm],
                                                 &nitem_face_comm[icomm]);
        
        out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
        
        // Check matching of input interface nodes with the set ones
        
        // Get input triangle nodes
        int** faceNodes2 = (int **) malloc(ncomm*sizeof(int *));
        for( int icomm = 0; icomm < ncomm; icomm++ ) {
          faceNodes2[icomm] = (int *) malloc(3*ntifc[icomm]*sizeof(int));
          for( i = 0; i < ntifc[icomm]; i++ ) {
              
            int faceID = ifc_tria_loc[icomm][i]-1;
            int faceID2 = locShF2globShF[faceID];

            v0 = face2node[faceID2][0];
            v1 = face2node[faceID2][1];
            v2 = face2node[faceID2][2];
    
            v0l = globV2locV[v0];
            v1l = globV2locV[v1];
            v2l = globV2locV[v2];

            //pos = ifc_tria_loc[icomm][i];
            faceNodes2[icomm][3*i]     = v0l+1; // tria_vert[3*(pos-1)];
            faceNodes2[icomm][3*i+1]   = v1l+1; // tria_vert[3*(pos-1)+1];
            faceNodes2[icomm][3*i+2]   = v2l+1; // tria_vert[3*(pos-1)+2];
          }
        }

        // Check matching of input interface triangles with the set ones
        if( !PMMG_Check_Set_FaceCommunicators(parmesh,ncomm,ntifc,
                                              color_face,faceNodes2) ) {
          printf("### FAILED:: Wrong set face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: Set the correct face communicators\n");
        }
        
        int *ref2       = (int*)calloc(nTriangles,sizeof(int));
        int *required2  = (int*)calloc(nTriangles,sizeof(int));
        int *triaNodes2 = (int*)calloc(3*nTriangles,sizeof(int));

        if ( PMMG_Get_triangles(parmesh,triaNodes2,ref2,required2) != 1 ) {
          fprintf(stderr,"FAILED:: Unable to get mesh triangles\n");
          ier = PMMG_STRONGFAILURE;
        }
        else
        {
            printf("### SUCCES:: retrieved all mesh triangles\n");
        }

        int** faceNodes_out = (int **) malloc(next_face_comm*sizeof(int *));
        for( int icomm = 0; icomm < next_face_comm; icomm++ )
        {
              faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
              for( i = 0; i < nitem_face_comm[icomm]; i++ )
              {
                  int pos = out_tria_loc[icomm][i];
                  faceNodes_out[icomm][3*i]   = triaNodes2[3*(pos-1)];
                  faceNodes_out[icomm][3*i+1] = triaNodes2[3*(pos-1)+1];
                  faceNodes_out[icomm][3*i+2] = triaNodes2[3*(pos-1)+2];
              }
        }

        
        // Check matching of input interface triangles with the output ones
        if( !PMMG_Check_Get_FaceCommunicators(parmesh,ncomm,ntifc,
                                              color_face,faceNodes2,
                                              next_face_comm,nitem_face_comm,
                                              color_face_out,faceNodes_out) )
        {
          printf("### FAILED:: Wrong retrieved face communicators!\n");
          MPI_Finalize();
          exit(EXIT_FAILURE);
        }
        else
        {
            printf("### SUCCES:: retrieved the correct face communicators\n");
        }
        
        int lvt = 0;
        int lvt_o = 0;
        std::map<int,int> sharedVrts_Owned;
        std::map<int,int> sharedVert;
        for( int icomm = 0; icomm < next_face_comm; icomm++ )
        {
            faceNodes_out[icomm] = (int *) malloc(3*nitem_face_comm[icomm]*sizeof(int));
            int nPartFace                 = nPartFace + nitem_face_comm[icomm];
            //rank2icomm[color_face_out[icomm]] = icomm;

            if(world_rank < color_face_out[icomm])
            {
                int nTshared_owned = nTshared_owned + nitem_face_comm[icomm];

                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft   = out_tria_loc[icomm][i];
                    int reff = ref2[ft-1];

                    for(int k=0;k<3;k++)
                    {
                        int vt = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];

                        if(sharedVrts_Owned.find(vt)==sharedVrts_Owned.end())
                        {
                            sharedVrts_Owned[vt] = lvt_o;
                            sharedVert[vt]       = lvt;

                            lvt++;
                            lvt_o++;
                        }
                    }
                }
            }
            else
            {
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft = out_tria_loc[icomm][i];

                    for(int k=0;k<3;k++)
                    {
                        int vt3 = triaNodes2[3*(ft-1)+k];
                        faceNodes_out[icomm][3*i+k] = triaNodes2[3*(ft-1)+k];
                        if(sharedVert.find(vt3)==sharedVert.end())
                        {
                            sharedVert[vt3] = lvt;
                            lvt++;
                        }
                    }
                }
            }
        }
    }
    else
    {
        //std::cout << "Outputting the new triangles..." << std::endl;
        //===================================================================
        int **out_tria_loc, **out_node_loc;
        int *nitem_face_comm,*nitem_node_comm;
        int next_face_comm, next_node_comm;
        int *color_node_out,*color_face_out;
        
        ier = PMMG_Get_numberOfNodeCommunicators(parmesh,&next_node_comm);
        
        color_node_out  = (int *) malloc(next_node_comm*sizeof(int));
        nitem_node_comm = (int *) malloc(next_node_comm*sizeof(int));
        for( int icomm=0; icomm<next_node_comm; icomm++ )
          ier = PMMG_Get_ithNodeCommunicatorSize(parmesh, icomm,
                                                 &color_node_out[icomm],
                                                 &nitem_node_comm[icomm]);
        
        // Get IDs of nodes on each interface //
        out_node_loc = (int **) malloc(next_node_comm*sizeof(int *));
        for( int icomm=0; icomm<next_node_comm; icomm++ )
          out_node_loc[icomm] = (int *) malloc(nitem_node_comm[icomm]*sizeof(int));
        ier = PMMG_Get_NodeCommunicator_nodes(parmesh, out_node_loc);
        
        //===================================================================
        ier = PMMG_Get_numberOfFaceCommunicators(parmesh,&next_face_comm);
        color_face_out  = (int *) malloc(next_face_comm*sizeof(int));
        nitem_face_comm = (int *) malloc(next_face_comm*sizeof(int));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          ier = PMMG_Get_ithFaceCommunicatorSize(parmesh, icomm,
                                                 &color_face_out[icomm],
                                                 &nitem_face_comm[icomm]);
        
        
        
        out_tria_loc = (int **) malloc(next_face_comm*sizeof(int *));
        for( int icomm=0; icomm<next_face_comm; icomm++ )
          out_tria_loc[icomm] = (int *) malloc(nitem_face_comm[icomm]*sizeof(int));
        ier = PMMG_Get_FaceCommunicator_faces(parmesh, out_tria_loc);
        
        // Check matching of input interface nodes with the set ones
        int *ref2       = (int*)calloc(nTrianglesOUT,sizeof(int));
        int *required2  = (int*)calloc(nTrianglesOUT,sizeof(int));
        int *triaNodes2 = (int*)calloc(3*nTrianglesOUT,sizeof(int));
        
//        if ( PMMG_Get_triangles(parmesh,triaNodes2,ref2,required2) != 1 ) {
//          fprintf(stderr,"FAILED:: Unable to get mesh triangles\n");
//          ier = PMMG_STRONGFAILURE;
//        }
        
        int nreq2   = 0;
        int pos     = 0;
        int cnt_bnd = 0;
        int cnt_th  = 0;
        int cnt_int = 0;
        
        FaceSetPointer m_PMMG_Face2RefPointer;
        
        //std::map<std::set<int>, int> PMMG_Shell2Prism;
        FaceSetPointer m_PMMG_ShellFace2PrismPointer;

        int c36 = 0;
        int gv0,gv1,gv2;
        int prismFound = 0;
        FaceSetPointer m_PMMG_ShellFacePointer;
        
        
        int faceCovered = 0;
        double tolerance = 1.0e-16;
        std::vector<int> vref1;
        
        std::map<int,int> outshell_g2l;
        std::vector<Vert*> outshellVerts;
        std::vector<std::vector<int> > outshellT;
        int locv = 0;
        for ( k=0; k<nTrianglesOUT; k++ )
        {
            pos = 3*k;
            if ( PMMG_Get_triangle(parmesh,&(triaNodes2[pos]),
                                         &(triaNodes2[pos+1]),
                                         &(triaNodes2[pos+2]),
                                         &(ref2[k]),
                                         &(required2[k])) != 1 )
            {
            fprintf(inm,"Unable to get mesh triangle %d \n",k);
            ier = PMMG_STRONGFAILURE;
            }

            
            gv0 = loc2globVid[triaNodes2[pos]];
            gv1 = loc2globVid[triaNodes2[pos+1]];
            gv2 = loc2globVid[triaNodes2[pos+2]];
        
            
            std::set<int> faceSh;
            faceSh.insert(gv0);
            faceSh.insert(gv1);
            faceSh.insert(gv2);
            
            std::vector<int> faceShVec(3);
            faceShVec[0] = gv0;
            faceShVec[1] = gv1;
            faceShVec[2] = gv2;
            

            if(ref2[k]!=0 && ref2[k]!=13)
            {
                
                FaceSharedPtr Face2RefPointer = std::shared_ptr<NekFace>(new NekFace(faceShVec));
                pair<FaceSetPointer::iterator, bool> testFace2RefPointer;
                testFace2RefPointer = m_PMMG_Face2RefPointer.insert(Face2RefPointer);
                
                if(testFace2RefPointer.second)
                {
                    (*testFace2RefPointer.first)->SetFaceRef(ref2[k]);
                }
                
                
            }
            if(ref2[k]==13)
            {
                std::vector<int> row(3);
                
                for(int vr=0;vr<3;vr++)
                {
                    row[vr] = triaNodes2[pos+vr];
                    int vertid = triaNodes2[pos+vr];
                    
                    if(outshell_g2l.find(vertid)==outshell_g2l.end())
                    {
                        outshell_g2l[vertid] = locv;
                        
                        Vert* vout = new Vert;
                        vout->x = vertOUT[(vertid-1)*3];
                        vout->y = vertOUT[(vertid-1)*3+1];
                        vout->z = vertOUT[(vertid-1)*3+2];
                        outshellVerts.push_back(vout);
                        
                        locv++;
                    }
                }
                
                outshellT.push_back(row);
                
                FaceSharedPtr ShellFacePointer = std::shared_ptr<NekFace>(new NekFace(faceShVec));
                pair<FaceSetPointer::iterator, bool> testInsPointer;
                testInsPointer = m_PMMG_ShellFacePointer.insert(ShellFacePointer);
                
                if(testInsPointer.second)
                {
                    (*testInsPointer.first)->SetFaceID(ref2[k]);
                }
                
                
                int flipper = 0;

                std::vector<int> shFace(3);
                shFace[0] = triaNodes2[pos];
                shFace[1] = triaNodes2[pos+1];
                shFace[2] = triaNodes2[pos+2];
                
                
                
                
                if(refOUT[triaNodes2[pos]-1]==86)
                {
                    clock_t t0,t1;
                    t0 = clock();
                    
                    flipper = 1;
                    int vertid = triaNodes2[pos]-1;
                    Vert* vcon = new Vert;
                    vcon->x = vertOUT[vertid*3];
                    vcon->y = vertOUT[vertid*3+1];
                    vcon->z = vertOUT[vertid*3+2];
                                        
                    std::map<int,Vert*>::iterator brutus;
                    for(brutus=shellVertCoord2Ref.begin();brutus!=shellVertCoord2Ref.end();brutus++)
                    {
                        int vrtref = brutus->first;
                        Vert* vrtCoord = brutus->second;
                        
                        double errx = fabs(vrtCoord->x - vcon->x);
                        double erry = fabs(vrtCoord->y - vcon->y);
                        double errz = fabs(vrtCoord->z - vcon->z);

                        if(errx < tolerance && erry < tolerance && errz < tolerance)
                        {
                            refOUT[triaNodes2[pos]-1] = vrtref;
                        }
                    }
                    
                    
                    
                    t1 = clock();
                    double duration = ( t1 - t0) / (double) CLOCKS_PER_SEC;
                    std::cout << "Warning:: There is an issue with the ref value of vertex " << vertid << " on rank " << world_rank << std::endl;
                    std::cout << "It is swapped back to the correct value by brute force which took " << duration << " seconds to find out of a set of " << shellVertCoord2Ref.size() << " vertices." << std::endl;
                    
                    
                }
                if(refOUT[triaNodes2[pos+1]-1]==86)
                {
                    clock_t t0,t1;
                    t0 = clock();
                    
                    flipper = 1;

                    int vertid = triaNodes2[pos+1]-1;
                    Vert* vcon = new Vert;
                    vcon->x = vertOUT[vertid*3];
                    vcon->y = vertOUT[vertid*3+1];
                    vcon->z = vertOUT[vertid*3+2];

                    std::map<int,Vert*>::iterator brutus;
                    for(brutus=shellVertCoord2Ref.begin();brutus!=shellVertCoord2Ref.end();brutus++)
                    {
                        int vrtref = brutus->first;
                        Vert* vrtCoord = brutus->second;
                        
                        double errx = fabs(vrtCoord->x - vcon->x);
                        double erry = fabs(vrtCoord->y - vcon->y);
                        double errz = fabs(vrtCoord->z - vcon->z);
                        if(errx < tolerance && erry < tolerance && errz < tolerance)
                        {
                            refOUT[triaNodes2[pos+1]-1] = vrtref;
                        }
                    }
                    
                    t1 = clock();
                    double duration = ( t1 - t0) / (double) CLOCKS_PER_SEC;
                    std::cout << "Warning:: There is an issue with the ref value of vertex " << vertid << " on rank " << world_rank << std::endl;
                    std::cout << "It is swapped back to the correct value by brute force which took " << duration << " seconds to find out of a set of " << shellVertCoord2Ref.size() << " vertices." << std::endl;
                    
                }
                if(refOUT[triaNodes2[pos+2]-1]==86)
                {
                    
                    clock_t t0,t1;
                    t0 = clock();
                    
                    
                    flipper = 1;

                    
                    int vertid = triaNodes2[pos+2]-1;
                    Vert* vcon = new Vert;
                    vcon->x = vertOUT[vertid*3];
                    vcon->y = vertOUT[vertid*3+1];
                    vcon->z = vertOUT[vertid*3+2];

                    std::map<int,Vert*>::iterator brutus;
                    for(brutus=shellVertCoord2Ref.begin();brutus!=shellVertCoord2Ref.end();brutus++)
                    {
                        int vrtref = brutus->first;
                        Vert* vrtCoord = brutus->second;
                        
                        double errx = fabs(vrtCoord->x - vcon->x);
                        double erry = fabs(vrtCoord->y - vcon->y);
                        double errz = fabs(vrtCoord->z - vcon->z);
                        
                        
                        if(errx < tolerance && erry < tolerance && errz < tolerance)
                        {
                            refOUT[triaNodes2[pos+2]-1] = vrtref;
                        }
                    }
                    
                    t1 = clock();
                    double duration = ( t1 - t0) / (double) CLOCKS_PER_SEC;
                    std::cout << "Warning:: There is an issue with the ref value of vertex " << vertid << " on rank " << world_rank << std::endl;
                    std::cout << "It is swapped back to the correct value by brute force which took " << duration << " seconds to find out of a set of " << shellVertCoord2Ref.size() << " vertices." << std::endl;
                }
                
                
                std::vector<int> refOnFace(3);
                refOnFace[0] = refOUT[triaNodes2[pos]-1];
                refOnFace[1] = refOUT[triaNodes2[pos+1]-1];
                refOnFace[2] = refOUT[triaNodes2[pos+2]-1];
                
                FaceSharedPtr RefOnFacePointer              = std::shared_ptr<NekFace>(new NekFace(refOnFace));
                FaceSetPointer::iterator RefsOnShellPointer = m_PMMG_RefsOnShell_2_Prism.find(RefOnFacePointer);

                // Mechanism to match up the shell faces from the adapted tetrahedra to the fixed prisms.
                if(RefsOnShellPointer!=m_PMMG_RefsOnShell_2_Prism.end())
                {
                    int PrismID = (*RefsOnShellPointer)->GetFaceLeftElement();
                    
                    
                    FaceSharedPtr FacePointer = std::shared_ptr<NekFace>(new NekFace(faceShVec));
                    
                    pair<FaceSetPointer::iterator, bool> Shell2PrismPointer;
                    Shell2PrismPointer = m_PMMG_ShellFace2PrismPointer.insert(FacePointer);
                    
                    if(Shell2PrismPointer.second)
                    {
                        (*Shell2PrismPointer.first)->SetFaceLeftElement(PrismID);
                    }
                }

                
                faceCovered++;
                // Need to catch the references here //
            }
         
            faceSh.clear();
    
            if ( required2 && required2[k] )  nreq2++;
        }
        
        
        if(outshellT.size()!=0 && debug == 1)
        {
            ofstream myfile_OUT;
                       
            string filename_shellOUT = "outputShell_" + std::to_string(world_rank) + ".dat";

            myfile_OUT.open(filename_shellOUT);
            myfile_OUT << "TITLE=\"outputShell.tec\"" << std::endl;
            myfile_OUT <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
            //ZONE N = 64, E = 48, DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL
            myfile_OUT <<"ZONE N = " << outshellVerts.size() << ", E = " << outshellT.size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
            for(int i=0;i<outshellVerts.size();i++)
            {
              myfile_OUT << outshellVerts[i]->x << "   " << outshellVerts[i]->y << "   " << outshellVerts[i]->z << std::endl;
            }

            for(int i=0;i<outshellT.size();i++)
            {
                int g0=outshellT[i][0];
                int g1=outshellT[i][1];
                int g2=outshellT[i][2];
                
                int l0=outshell_g2l[g0];
                int l1=outshell_g2l[g1];
                int l2=outshell_g2l[g2];
                
              myfile_OUT << l0+1 << "    " << l1+1 << "   " << l2+1 << std::endl;
            }
            myfile_OUT.close();
        }
        
        
        
        
        DistributedParallelState* pmmg_shellfacedist = new DistributedParallelState(m_PMMG_ShellFacePointer.size(),comm);
        DistributedParallelState* PMMG_Shell2Prismdist = new DistributedParallelState(m_PMMG_ShellFace2PrismPointer.size(),comm);
        
        if(m_PMMG_ShellFacePointer.size()!=m_PMMG_ShellFace2PrismPointer.size())
        {
            std::cout << world_rank << " m_PMMG_ShellFacePointer size " << m_PMMG_ShellFacePointer.size() << " pmmg_shellfacedist->getNel() " << pmmg_shellfacedist->getNel() << " PMMG_Shell2Prismdist->getNel() = "<< PMMG_Shell2Prismdist->getNel() << " prismFound-> " << prismFound << " " << faceCovered << " " << m_PMMG_ShellFace2PrismPointer.size() <<  std::endl;
        }
        
        int itt2            = 0;
        int nTshared_owned  = 0;

        int vt,ft,gvt,gft;
        
        FaceSetPointer m_PMMG_SharedFacePointer;
        FaceSetPointer m_PMMG_OwnedSharedFacePointer;

        std::set<int> PMMG_SharedVertsOwned;
        int nPartFace = 0;
        std::map<int,int> rank2icomm;
        std::map<int,std::vector<int> > LocateOppositeFace;
        std::map<int,int> vertonrank;
        int locID_NotShVrt = 0;
        std::vector<int> PMMG_SharedVertsOwned_vec;
        std::vector<int> PMMG_SharedVertsOwnedRank_vec;
        std::set<int> PMMG_SharedVertices;
        std::map<int,std::vector<int> > Color_SharedOwned;

        for(int icomm = 0; icomm < next_face_comm; icomm++ )
        {
            nPartFace                           = nPartFace + nitem_face_comm[icomm];
            rank2icomm[color_face_out[icomm]] = icomm;
                    
            if(world_rank < color_face_out[icomm])
            {
                nTshared_owned = nTshared_owned + nitem_face_comm[icomm];
                
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft       = out_tria_loc[icomm][i];
                    int face_ref = ref2[ft-1];
                    //Color_SharedOwned[icomm].push_back(i);
                    std::set<int> faceSh;
                    std::vector<int> faceShVec(3);
                    for(int k=0;k<3;k++)
                    {
                        int vt  = triaNodes2[3*(ft-1)+k];
                        int gvt = loc2globVid[vt];
                        faceSh.insert(gvt);
                        faceShVec[k] = gvt;
                        
                        if(PMMG_SharedVertices.find(gvt)==PMMG_SharedVertices.end())
                        {
                            PMMG_SharedVertices.insert(gvt);
                        }
                        
                        if(PMMG_SharedVertsOwned.find(gvt)==PMMG_SharedVertsOwned.end())
                        {
                            PMMG_SharedVertsOwned.insert(gvt);
                            PMMG_SharedVertsOwned_vec.push_back(gvt);
                            PMMG_SharedVertsOwnedRank_vec.push_back(world_rank);
                        }
                    }
                    
                    Color_SharedOwned[color_face_out[icomm]].push_back(i);
                    
                    FaceSharedPtr sharedFacePointer = std::shared_ptr<NekFace>(new NekFace(faceShVec));
                    pair<FaceSetPointer::iterator, bool> SharedFPointer;
                    SharedFPointer      = m_PMMG_SharedFacePointer.insert(sharedFacePointer);
                    pair<FaceSetPointer::iterator, bool> OwnedSharedFPointer;
                    OwnedSharedFPointer = m_PMMG_OwnedSharedFacePointer.insert(sharedFacePointer);
                    
                    if(SharedFPointer.second)
                    {
                        (*SharedFPointer.first)->SetFaceID(ft);
                    }
                    
                    if(OwnedSharedFPointer.second)
                    {
                        (*OwnedSharedFPointer.first)->SetFaceID(ft);
                    }
                    
                    faceSh.clear();
                }
            }
            else
            {
                for( i = 0; i < nitem_face_comm[icomm]; i++ )
                {
                    int ft       = out_tria_loc[icomm][i];
                    
                    std::set<int> faceSh;
                    std::vector<int> faceShVec(3);
                    for(int k=0;k<3;k++)
                    {
                        int vt  = triaNodes2[3*(ft-1)+k];
                        int gvt = loc2globVid[vt];
                        faceSh.insert(gvt);
                        faceShVec[k]=gvt;
                        if(PMMG_SharedVertices.find(gvt)==PMMG_SharedVertices.end())
                        {
                            PMMG_SharedVertices.insert(gvt);
                        }
                    }
//                    if(PMMG_SharedFaces.find(faceSh)==PMMG_SharedFaces.end())
//                    {
//                        PMMG_SharedFaces[faceSh]=ft;
//                    }
                    
                    
                    FaceSharedPtr sharedFacePointer = std::shared_ptr<NekFace>(new NekFace(faceShVec));
                    pair<FaceSetPointer::iterator, bool> testInsPointer;
                    testInsPointer = m_PMMG_SharedFacePointer.insert(sharedFacePointer);
                    if(testInsPointer.second)
                    {
                        (*testInsPointer.first)->SetFaceID(ft);
                    }
                    
                    faceSh.clear();
                }
            }
        }
       
        
        //std::cout << world_rank << " m_PMMG_SharedFacePointer " << m_PMMG_SharedFacePointer.size() << " " << PMMG_SharedFaces.size() << std::endl;
        
        
        int nLocallyOwnedVerts                      = PMMG_SharedVertsOwned.size();
        DistributedParallelState* locallyOwnedVerts = new DistributedParallelState(nLocallyOwnedVerts,comm);
        int* OwnedVertsDistri                       = new int[locallyOwnedVerts->getNel()];
        int* OwnedVertsDistriRank                   = new int[locallyOwnedVerts->getNel()];
        int* PMMG_SharedVertsOwned_arr              = new int[nLocallyOwnedVerts];
        int* PMMG_SharedVertsOwnedRank_arr          = new int[nLocallyOwnedVerts];
        
        for(int i=0;i<nLocallyOwnedVerts;i++)
        {
            PMMG_SharedVertsOwned_arr[i]            = PMMG_SharedVertsOwned_vec[i];
            PMMG_SharedVertsOwnedRank_arr[i]        = PMMG_SharedVertsOwnedRank_vec[i];
        }
        
        MPI_Allgatherv(&PMMG_SharedVertsOwned_arr[0],
                       nLocallyOwnedVerts,
                       MPI_INT,
                       OwnedVertsDistri,
                       locallyOwnedVerts->getNlocs(),
                       locallyOwnedVerts->getOffsets(),
                       MPI_INT, comm);
        
        MPI_Allgatherv(&PMMG_SharedVertsOwnedRank_arr[0],
                       nLocallyOwnedVerts,
                       MPI_INT,
                       OwnedVertsDistriRank,
                       locallyOwnedVerts->getNlocs(),
                       locallyOwnedVerts->getOffsets(),
                       MPI_INT, comm);
        
        std::map<int,int> ActualOwnedVertDistr;
        std::map<int,std::vector<int> > ActualOwnedVertDistr_map;
        int ownedID;
        
        for(int i=0;i<locallyOwnedVerts->getNel();i++)
        {
            int gvd = OwnedVertsDistri[i];
            
            if(ActualOwnedVertDistr.find(gvd)==ActualOwnedVertDistr.end())
            {
                ActualOwnedVertDistr_map[gvd].push_back(OwnedVertsDistriRank[i]);

                ActualOwnedVertDistr[gvd] = OwnedVertsDistriRank[i];
                ownedID = OwnedVertsDistriRank[i];
            }
            else
            {
                if(OwnedVertsDistriRank[i]<ActualOwnedVertDistr[gvd])
                {
                    ActualOwnedVertDistr[gvd] = OwnedVertsDistriRank[i];
                    ownedID = OwnedVertsDistriRank[i];
                }
                else
                {
                    ActualOwnedVertDistr[gvd] = OwnedVertsDistriRank[i];
                    ownedID = ActualOwnedVertDistr[gvd];
                }
            }
        }
        
        delete[] OwnedVertsDistri;
        delete[] OwnedVertsDistriRank;
        delete[] PMMG_SharedVertsOwned_arr;
        delete[] PMMG_SharedVertsOwnedRank_arr;
        
        std::map<int,std::vector<int> > ActualOwnedVertDistr_map_update;

        std::map<int,std::vector<int> >::iterator itm;
        std::map<int,int>::iterator iitm;
        int tel = 0;
        int* tells = new int[world_size];
        for(int u=0;u<world_size;u++)
        {
            tells[u] = 0;
        }
        
        std::map<int,int> LocationSharedVert;
        for(itm=ActualOwnedVertDistr_map.begin();itm!=ActualOwnedVertDistr_map.end();itm++)
        {
            int gv = itm->first;
            int ra = *min_element(itm->second.begin(), itm->second.end());
            tells[ra] = tells[ra]+1;
            ActualOwnedVertDistr_map_update[ra].push_back(gv);
            LocationSharedVert[gv]=ra;
            tel++;
        }

        int* tetraOUT = new int[nTetrahedraOUT*4];
        std::vector<std::vector<int> > tetrasOUT;
        std::vector<int> NonSharedVrts_vec;

        std::map<int,std::vector<double> > coords_int;
        std::map<int,int> lhshown;
        std::map<int,int> lhbnd;
        std::map<int,int> lhshell;
        std::map<int,int> rhshell;
        std::map<int,int> lhsh;
        std::map<int,int> lh;
        std::map<int,int> rh;
        
        DistributedParallelState* distTetraOut = new DistributedParallelState(nTetrahedraOUT,comm);
        int* TetraOUT_offsets = distTetraOut->getOffsets();
        
        std::map<int,std::vector<int> > ienOUT;
        
        
        std::map<int,FaceSharedPtr> BoundaryFaces;
        std::map<int,FaceSharedPtr> InternalFaces;
        std::map<int,FaceSharedPtr> SharedFaces;
        std::map<int,std::vector<int> > fmShell;
        std::map<int,std::vector<int> > fm;
        
        FaceSetPointer m_FaceSetPointer;
        std::map<int,std::vector<int> > bcmap;
        
        std::map<int,std::vector<int> > pfmInt;
        std::map<int,std::vector<int> > pfmSha;
        std::map<int,std::vector<int> > pfm;
        
        Array<int>* parmmg_iet = new Array<int>(nTetrahedraOUT,1);
        int q;
        
        std::map<int,int> gvid2shid;
        int fid  = 0;
        int lshf = 0;
        int tetra_faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
        
        int NoNShFaces = 0;
        
        std::map<int,int> locShF2globShF;
        std::map<int,int> locBF2globBF;
        std::map<int,int> NonSharedVrts;
        int ngv = 0;
        
        std::vector<int> Elvrts(4);
        int fv0,fv1,fv2;
        Vec3D* v0 = new Vec3D;
        Vec3D* v1 = new Vec3D;
        Vec3D* r0 = new Vec3D;

        Vert* Vface = new Vert;
        int negit = 0;
        
        int curElID = TetraOUT_offsets[world_rank]+1+nPrismsTot;

        int* refTET = new int[nTetrahedraOUT];
        
        std::map<int,int> tag2shelltag;
        std::map<int,int> shelltag2tag;
        int hellofound = 0;
        
        //std::cout << "Start ID for tets on rank = " << curElID << std::endl;
        
        for ( k=0; k<nTetrahedraOUT; k++ )
        {
            parmmg_iet->setVal(k,0,2);
            pos = 4*k;
            if ( PMMG_Get_tetrahedron(parmesh,
                                        &(tetraOUT[pos  ]),&(tetraOUT[pos+1]),
                                        &(tetraOUT[pos+2]),&(tetraOUT[pos+3]),
                                        &(refTET[k]),&(required[k])) != 1 )
            {
                fprintf(inm,"Unable to get mesh tetra %d \n",k);
                ier = PMMG_STRONGFAILURE;
            }
            
            Elvrts[0] = tetraOUT[pos];
            Elvrts[1] = tetraOUT[pos+1];
            Elvrts[2] = tetraOUT[pos+2];
            Elvrts[3] = tetraOUT[pos+3];
            
            double* P = new double[4*3];
            
            P[0*3+0]=vertOUT[(Elvrts[0]-1)*3];
            P[0*3+1]=vertOUT[(Elvrts[0]-1)*3+1];
            P[0*3+2]=vertOUT[(Elvrts[0]-1)*3+2];
            
            P[1*3+0]=vertOUT[(Elvrts[1]-1)*3];
            P[1*3+1]=vertOUT[(Elvrts[1]-1)*3+1];
            P[1*3+2]=vertOUT[(Elvrts[1]-1)*3+2];
            
            P[2*3+0]=vertOUT[(Elvrts[2]-1)*3];
            P[2*3+1]=vertOUT[(Elvrts[2]-1)*3+1];
            P[2*3+2]=vertOUT[(Elvrts[2]-1)*3+2];
            
            P[3*3+0]=vertOUT[(Elvrts[3]-1)*3];
            P[3*3+1]=vertOUT[(Elvrts[3]-1)*3+1];
            P[3*3+2]=vertOUT[(Elvrts[3]-1)*3+2];
            
            ienOUT[curElID] = Elvrts;
            tetrasOUT.push_back(Elvrts);
            double Vtet = GetQualityTetrahedra(P);
            Vert* vCenter = ComputeCentroidCoord(P, 4);
            
            if(Vtet<0.0)
            {
                std::cout << "Error " << Vtet << std::endl;
            }
            
            for(int u=0;u<4;u++)
            {
                fv0 = loc2globVid[Elvrts[tetra_faces[u][0]]];
                fv1 = loc2globVid[Elvrts[tetra_faces[u][1]]];
                fv2 = loc2globVid[Elvrts[tetra_faces[u][2]]];
                
                std::set<int> Face;
                Face.insert(fv0);
                Face.insert(fv1);
                Face.insert(fv2);
                
                std::vector<int> FaceVec(3);
                FaceVec[0] = fv0;
                FaceVec[1] = fv1;
                FaceVec[2] = fv2;
                FaceSharedPtr f2ePointer = std::shared_ptr<NekFace>(new NekFace(FaceVec));
                pair<FaceSetPointer::iterator, bool> testInsPointer;
                testInsPointer = m_FaceSetPointer.insert(f2ePointer);
                
                if(testInsPointer.second)
                {
                    (*testInsPointer.first)->SetFaceID(fid);
                    std::vector<int> fce(3);
                    
                    fce[0]  = fv0;
                    fce[1]  = fv1;
                    fce[2]  = fv2;
                    
                    double v0x = vertOUT[(Elvrts[tetra_faces[u][0]]-1)*3];
                    double v0y = vertOUT[(Elvrts[tetra_faces[u][0]]-1)*3+1];
                    double v0z = vertOUT[(Elvrts[tetra_faces[u][0]]-1)*3+2];
                    
                    double v1x = vertOUT[(Elvrts[tetra_faces[u][1]]-1)*3];
                    double v1y = vertOUT[(Elvrts[tetra_faces[u][1]]-1)*3+1];
                    double v1z = vertOUT[(Elvrts[tetra_faces[u][1]]-1)*3+2];
                    
                    double v2x = vertOUT[(Elvrts[tetra_faces[u][2]]-1)*3];
                    double v2y = vertOUT[(Elvrts[tetra_faces[u][2]]-1)*3+1];
                    double v2z = vertOUT[(Elvrts[tetra_faces[u][2]]-1)*3+2];
                    
                    Vface->x = (v0x+v1x+v2x)/3.0;
                    Vface->y = (v0y+v1y+v2y)/3.0;
                    Vface->z = (v0z+v1z+v2z)/3.0;
                    
                    r0->c0 = (Vface->x-vCenter->x);///Lr;
                    r0->c1 = (Vface->y-vCenter->y);///Lr;
                    r0->c2 = (Vface->z-vCenter->z);///Lr;
                    
                    v0->c0 = v1x-v0x;
                    v0->c1 = v1y-v0y;
                    v0->c2 = v1z-v0z;

                    v1->c0 = v2x-v0x;
                    v1->c1 = v2y-v0y;
                    v1->c2 = v2z-v0z;
                    
                    Vec3D* n0        = ComputeSurfaceNormal(v0,v1);
                    double orient0   = DotVec3D(r0,n0);
                    if(orient0<0.0)
                    {
                        negit++;
                    }
                    
                    fm[fid] = fce;
                    
                    FaceSetPointer::iterator SharedFPointer         = m_PMMG_SharedFacePointer.find(f2ePointer);
                    FaceSetPointer::iterator OwnedSharedFPointer    = m_PMMG_OwnedSharedFacePointer.find(f2ePointer);
                    FaceSetPointer::iterator testShellPointer       = m_PMMG_ShellFacePointer.find(f2ePointer);
                    FaceSetPointer::iterator testFace2RefPointer    = m_PMMG_Face2RefPointer.find(f2ePointer);
                    FaceSetPointer::iterator testShell2PrismPointer = m_PMMG_ShellFace2PrismPointer.find(f2ePointer);

                    if(SharedFPointer != m_PMMG_SharedFacePointer.end())
                    {
                        lshf                     = (*SharedFPointer)->GetFaceID();
                        locShF2globShF[lshf]     = fid;
                        lhsh[fid]                = curElID;
                    }
                    
                    
                    if(SharedFPointer == m_PMMG_SharedFacePointer.end()
                       && testFace2RefPointer == m_PMMG_Face2RefPointer.end()
                       && testShellPointer == m_PMMG_ShellFacePointer.end())
                    {
                        InternalFaces[fid] = (*testInsPointer.first);
                        //fmInt[fid]  = fce;
                        lh[fid]     = curElID;
                        
                        for(int s=0;s<3;s++)
                        {
                            int gvm2 = fce[s];
                            if(NonSharedVrts.find(gvm2)==NonSharedVrts.end() &&
                                    PMMG_SharedVertices.find(gvm2)==PMMG_SharedVertices.end())
                            {
                                NonSharedVrts[gvm2] = ngv;
                                ngv++;
                            }
                        }
                    }
                    
                    
                    
                    if(OwnedSharedFPointer != m_PMMG_OwnedSharedFacePointer.end())
                    {
                        
                        SharedFaces[fid]         = (*testInsPointer.first);
                        //fmSha[fid]               = fce;
                        lhshown[fid]             = curElID;
                    }
                    
                    
                    if(testShellPointer != m_PMMG_ShellFacePointer.end())
                    {
                        fmShell[fid]        = fce;
                        lhshell[fid]        = curElID;
                       
                        
                        if(testShell2PrismPointer != m_PMMG_ShellFace2PrismPointer.end())
                        {
                            rhshell[fid] = (*testShell2PrismPointer)->GetFaceLeftElement();
                        
                            for(int y=0;y<3;y++)
                            {
                                int shelltag = refOUT[Elvrts[tetra_faces[u][y]]-1];
                                int vtag     = fce[y];
                                
                                if(tag2shelltag.find(vtag)==tag2shelltag.end())
                                {
                                    tag2shelltag[vtag]=shelltag;
                                    shelltag2tag[shelltag]=vtag;
                                }
                            }
                            
                            hellofound++;
                        }
                    }
                    
                    
                    
                    if(testFace2RefPointer != m_PMMG_Face2RefPointer.end())
                    {
                        int FaceRef         = (*testFace2RefPointer)->GetFaceRef();
                        
                        BoundaryFaces[fid] = f2ePointer;
                        //fmBnd[fid]          = fce;
                        lhbnd[fid]          = curElID;
                        bcmap[FaceRef].push_back(fid);
                    }
                }
                else
                {
                    
                    int fid_n     = (*testInsPointer.first)->GetFaceID();
                    rh[fid_n]     = curElID;
                }
                
                    
                fid++;
                
                Face.clear();
            }
            
            curElID++;

            delete[] P;
        }
       
        delete[] refTET;
        
        PMMG_Free_all(PMMG_ARG_start,
                      PMMG_ARG_ppParMesh,&parmesh,
                      PMMG_ARG_end);
        
        
        if(nTetrahedraOUT!=0 && debug==1)
        {
            string filename = "pmmg_" + std::to_string(world_rank) + ".dat";
            OutputMesh_PMMG_V2(nVerticesOUT,vertOUT,nTetrahedraOUT,tetraOUT,filename);
        }

        
        //std::cout << "BoundaryFaces " << BoundaryFaces.size() << std::endl;
        
        
        
        if(m_PMMG_ShellFace2PrismPointer.size()!=m_PMMG_ShellFacePointer.size())
        {
            std::cout << "NOTIFIED !!! " << world_rank << " --> " << rhshell.size() << " " << fmShell.size() << " " << hellofound << " " << m_PMMG_ShellFace2PrismPointer.size() << " " << m_PMMG_ShellFacePointer.size() << std::endl;
        }
       
        
        std::map<int,int> tag2shelltag_glob = AllGatherMap(tag2shelltag,comm);
        std::map<int,int> shelltag2tag_glob = AllGatherMap(shelltag2tag,comm);

        int nLocIntVrts         = NonSharedVrts.size();
        int nLocShVrts          = ActualOwnedVertDistr_map_update[world_rank].size();
        int nLocTotVrts         = nLocIntVrts+nLocShVrts;
        DistributedParallelState* distnLocTotVrts = new DistributedParallelState(nLocTotVrts,comm);

        int vert = 0;//distnLocTotVrts->getOffsets()[world_rank];
        int ToTVrts =  distnLocTotVrts->getNel();
        int* ToTVrts_offsets = distnLocTotVrts->getOffsets();
        int ToTVrts_offset = ToTVrts_offsets[world_rank];
        Array<double>* xcn_parmmg = new Array<double>(nLocTotVrts,3);
        std::map<int,int> tag2glob;
        std::map<int,int> glob2tag;
        
        
        
        int nPrismVerts_tmp = distPrismIntVerts->getNel()+distPrismShaVerts->getNel();
        

        for(iitm=NonSharedVrts.begin();iitm!=NonSharedVrts.end();iitm++)
        {
            int tag = iitm->first;
            int lvert = glob2locVid[tag];

            double xc = vertOUT[(lvert-1)*3+0];
            double yc = vertOUT[(lvert-1)*3+1];
            double zc = vertOUT[(lvert-1)*3+2];

            
            xcn_parmmg->setVal(vert,0,xc);
            xcn_parmmg->setVal(vert,1,yc);
            xcn_parmmg->setVal(vert,2,zc);
            
            tag2glob[tag] = vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp;
            glob2tag[vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp] = tag;
            
            vert++;
        }
        
        
        for(int i=0;i<nLocShVrts;i++)
        {
            int tag = ActualOwnedVertDistr_map_update[world_rank][i];
            int lvert = glob2locVid[tag];

            double xc = vertOUT[(lvert-1)*3+0];
            double yc = vertOUT[(lvert-1)*3+1];
            double zc = vertOUT[(lvert-1)*3+2];

            xcn_parmmg->setVal(vert,0,xc);
            xcn_parmmg->setVal(vert,1,yc);
            xcn_parmmg->setVal(vert,2,zc);

            tag2glob[tag] = vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp;
            glob2tag[vert+distnLocTotVrts->getOffsets()[world_rank]+1+nPrismVerts_tmp] = tag;

            vert++;
        }
        
        int* updateGlobSharedVrtID       = new int[LocationSharedVert.size()];
        int* updateGlobSharedVrtID_red   = new int[LocationSharedVert.size()];
        
        int* originalGlobSharedVrtID     = new int[LocationSharedVert.size()];
        int* originalGlobSharedVrtID_red = new int[LocationSharedVert.size()];
        int ig = 0;

        for(iitm=LocationSharedVert.begin();iitm!=LocationSharedVert.end();iitm++)
        {
            int tag = iitm->first;
            int ra  = iitm->second;
            
            originalGlobSharedVrtID_red[ig] = 0;
            updateGlobSharedVrtID_red[ig]   = 0;
            
            if(ra == world_rank)
            {
                originalGlobSharedVrtID[ig] = tag;
                updateGlobSharedVrtID[ig]   = tag2glob[tag];
            }
            else
            {
                originalGlobSharedVrtID[ig] = 0;
                updateGlobSharedVrtID[ig]   = 0;
                
            }
            ig++;
        }
        
        //std::cout << "Performing an AllReduce of the global vertex IDs of the shared Vertices between partitions." << std::endl;
        
        MPI_Allreduce(originalGlobSharedVrtID,
                      originalGlobSharedVrtID_red,
                      LocationSharedVert.size(),
                      MPI_INT, MPI_SUM, comm);
        
        MPI_Allreduce(updateGlobSharedVrtID,
                      updateGlobSharedVrtID_red,
                      LocationSharedVert.size(),
                      MPI_INT, MPI_SUM, comm);
        
        std::map<int,int> LocationSharedVert_update;
        for(int i=0;i<LocationSharedVert.size();i++)
        {
            LocationSharedVert_update[originalGlobSharedVrtID_red[i]] = updateGlobSharedVrtID_red[i];
        }
    
        int nBndFaces       = m_PMMG_Face2RefPointer.size();
        int nShaFaces       = m_PMMG_OwnedSharedFacePointer.size();
        int nIntFaces       = InternalFaces.size();
        int nLocFaceNBnd    = InternalFaces.size()+SharedFaces.size()+fmShell.size();
        Array<int>* ifnOUT  = new Array<int>(nLocFaceNBnd,8);
        std::map<int,int> updateID;
        int ftot            =  0;
        std::vector<int> fq(3);
        int flag            = -1;
        int flagged         =  0;
        
        //std::cout <<  "right-handed " << world_rank << " " << rh.size() << " " << nShaFaces << " VS " << InternalFaces.size() << " " << SharedFaces.size() << " " << fmShell.size() << std::endl;
        
        int nothere2 = 0;
        
        int testElem = 9;
        int testComp = 10;
        
        
        Vert* Vijk = new Vert;
        Vert* VcF = new Vert;

        std::map<int,FaceSharedPtr>::iterator itf;
        //        for(itm=fmInt.begin();itm!=fmInt.end();itm++)

        for(itf=InternalFaces.begin();itf!=InternalFaces.end();itf++)
        {
            fid             = itf->first;
            updateID[fid]   = ftot;
            ifnOUT->setVal(ftot,0,3);
            flag            = -1;
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            
            std::vector<int> edges = InternalFaces[fid]->GetEdgeIDs();
            for(int q=0;q<3;q++)
            {
                int vertexID = InternalFaces[fid]->GetEdgeIDs()[q];

                int lvert = glob2locVid[vertexID];
                Vert* Vf = new Vert;
                Vf->x  = vertOUT[(lvert-1)*3+0];
                Vf->y  = vertOUT[(lvert-1)*3+1];
                Vf->z  = vertOUT[(lvert-1)*3+2];
                
                VcF->x = VcF->x + Vf->x;
                VcF->y = VcF->y + Vf->y;
                VcF->z = VcF->z + Vf->z;
                
                Vfaces.push_back(Vf);
                
                
                if(LocationSharedVert_update.find(vertexID)!=LocationSharedVert_update.end())
                {
                    fq[q] = LocationSharedVert_update[vertexID];
                    flag  = q;
                }
                else
                {
                    fq[q] = tag2glob[vertexID];
                }
            }
            
            VcF->x = VcF->x/3.0;
            VcF->y = VcF->y/3.0;
            VcF->z = VcF->z/3.0;
            
            int gv0 = fq[0];
            int gv1 = fq[1];
            int gv2 = fq[2];
            
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            //ienOUT[curElID] = Elvrts
            for(int u=0;u<ienOUT[lh[fid]].size();u++)
            {
                
                Vijk->x = Vijk->x + vertOUT[(ienOUT[lh[fid]][u]-1)*3];
                Vijk->y = Vijk->y + vertOUT[(ienOUT[lh[fid]][u]-1)*3+1];
                Vijk->z = Vijk->z + vertOUT[(ienOUT[lh[fid]][u]-1)*3+2];
                
            }
            
            Vijk->x = Vijk->x/ienOUT[lh[fid]].size();
            Vijk->y = Vijk->y/ienOUT[lh[fid]].size();
            Vijk->z = Vijk->z/ienOUT[lh[fid]].size();
            
            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            if(orient0 < 0.0)
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv2);
                ifnOUT->setVal(ftot,3,gv1);
                ifnOUT->setVal(ftot,4,0);
            }
            else
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv1);
                ifnOUT->setVal(ftot,3,gv2);
                ifnOUT->setVal(ftot,4,0);
            }
            
            ifnOUT->setVal(ftot,5,rh[fid]);
            ifnOUT->setVal(ftot,6,lh[fid]);
            ifnOUT->setVal(ftot,7,2);
            
            if(rh[fid] == 0)
            {
            	std::cout <<"Found the face -> " << fid << " " << world_rank << " " << lh[fid] << std::endl;
            }
                ftot++;
            }
        
        //Color_SharedOwned
        
        ScheduleObj* ish_schedule = DoScheduling(Color_SharedOwned,comm);
        std::map<int,std::vector<int> > recv_ids;
        std::map<int,std::vector<int> >::iterator it;
          
        for(int q=0;q<world_size;q++)
        {
            if(world_rank==q)
            {
                int i=0;
                for (it = Color_SharedOwned.begin(); it != Color_SharedOwned.end(); it++)
                {
                    int n_req           = it->second.size();
                    int dest            = it->first;

                    MPI_Send(&n_req, 1, MPI_INT, dest, 6798+78*dest, comm);
                    MPI_Send(&it->second[0], n_req, MPI_INT, dest, 14876+dest, comm);
                    i++;
                }
            }
            else if (ish_schedule->SendFromRank2Rank[q].find( world_rank ) != ish_schedule->SendFromRank2Rank[q].end())
            {
                int n_reqstd_ids;
                MPI_Recv(&n_reqstd_ids, 1, MPI_INT, q, 6798+78*world_rank, comm, MPI_STATUS_IGNORE);

                std::vector<int> recv_reqstd_ids(n_reqstd_ids);
                
                MPI_Recv(&recv_reqstd_ids[0], n_reqstd_ids, MPI_INT, q, 14876+world_rank, comm, MPI_STATUS_IGNORE);
                recv_ids[q] = recv_reqstd_ids;
            }
        }
        
        std::map<int,std::vector<int> > sendEl;
        std::map<int,std::vector<int> >::iterator rcvit;
        for(rcvit=recv_ids.begin();rcvit!=recv_ids.end();rcvit++)
        {
            int frank = rcvit->first;
            
            int icomm = rank2icomm[frank];
            int nF    = rcvit->second.size();
            
            for(int j=0;j<nF;j++)
            {
                int fidInt = rcvit->second[j];
                int ofid   = out_tria_loc[icomm][fidInt];
                int nfid   = locShF2globShF[ofid];
                sendEl[frank].push_back(lhsh[nfid]);
                
                if(lhsh[nfid] == 0)
                {
                    std::cout << "yep commi not correct " << std::endl;
                }
            }
        }
        
        ScheduleObj* ishBack_schedule = DoScheduling(sendEl,comm);

        std::map<int,std::vector<int> > adj_ids;
        for(int q=0;q<world_size;q++)
        {
            if(world_rank==q)
            {
                int i=0;
                for (it = sendEl.begin(); it != sendEl.end(); it++)
                {
                    int n_req           = it->second.size();
                    int dest            = it->first;

                    MPI_Send(&n_req, 1,
                            MPI_INT, dest,
                            6798+78000*dest, comm);
                    
                    MPI_Send(&it->second[0],
                            n_req, MPI_INT,
                            dest, 14876000+dest, comm);

                    i++;
                }
            }
            else if (ishBack_schedule->SendFromRank2Rank[q].find( world_rank ) != ishBack_schedule->SendFromRank2Rank[q].end())
            {
                int n_reqstd_ids;
                
                MPI_Recv(&n_reqstd_ids,
                        1, MPI_INT, q,
                        6798+78000*world_rank,
                        comm, MPI_STATUS_IGNORE);

                std::vector<int> recv_reqstd_ids(n_reqstd_ids);
                
                MPI_Recv(&recv_reqstd_ids[0],
                        n_reqstd_ids,
                        MPI_INT, q,
                        14876000+world_rank,
                        comm, MPI_STATUS_IGNORE);
                
                adj_ids[q] = recv_reqstd_ids;

            }
        }
        
        DistributedParallelState* rhbefore = new DistributedParallelState(rh.size(),comm);
        DistributedParallelState* lhbefore = new DistributedParallelState(lh.size(),comm);

        std::map<int,int> adjElements;
        int fid_loc,fid_glo;
        int telli = 0;
        std::set<int> frh;
        int outside = 0;
        for(itm=Color_SharedOwned.begin();itm!=Color_SharedOwned.end();itm++)
        {
            int rrank = itm->first;
            
            for(int j=0;j<itm->second.size();j++)
            {
                fid_loc              = itm->second[j];
                int icomm            = rank2icomm[rrank];
                int ft               = out_tria_loc[icomm][fid_loc];
                fid_glo              = locShF2globShF[ft];
                rh[fid_glo]          = adj_ids[itm->first][j];
                
                if(frh.find(fid_glo)==frh.end())
                {
                    adjElements[fid_glo] = adj_ids[itm->first][j];
                    
                    frh.insert(fid_glo);
                    outside++;
                }
                
                telli++;
            }
        }
        
        //DistributedParallelState* rhafter  = new DistributedParallelState(rh.size(),comm);
        //DistributedParallelState* lhafter  = new DistributedParallelState(lh.size(),comm);
        //DistributedParallelState* adjafter = new DistributedParallelState(adjElements.size(),comm);
        
        int nothere =  0;
        int elLh    = -1;
        int elRh    = -1;
        int buthere =  0;
        
        int ftot_inter = ftot;
        //for(itm=fmSha.begin();itm!=fmSha.end();itm++)

        for(itf=SharedFaces.begin();itf!=SharedFaces.end();itf++)
        {
            fid             = itf->first;
            updateID[fid]   = ftot;
            ifnOUT->setVal(ftot,0,3);
            flag = -1;
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            for(int q=0;q<3;q++)
            {
                int lvert = glob2locVid[itf->second->GetEdgeIDs()[q]];
                Vert* Vf = new Vert;
                Vf->x  = vertOUT[(lvert-1)*3+0];
                Vf->y  = vertOUT[(lvert-1)*3+1];
                Vf->z  = vertOUT[(lvert-1)*3+2];
                
                VcF->x = VcF->x + Vf->x;
                VcF->y = VcF->y + Vf->y;
                VcF->z = VcF->z + Vf->z;
                
                Vfaces.push_back(Vf);
                
                
                if(LocationSharedVert_update.find(itf->second->GetEdgeIDs()[q])!=LocationSharedVert_update.end())
                {
                    fq[q] = LocationSharedVert_update[itf->second->GetEdgeIDs()[q]];
                    flag = q;
                }
                else
                {
                    fq[q] = tag2glob[itf->second->GetEdgeIDs()[q]];
                }
            }
            
            VcF->x = VcF->x/3.0;
            VcF->y = VcF->y/3.0;
            VcF->z = VcF->z/3.0;
            
            int gv0 = fq[0];
            int gv1 = fq[1];
            int gv2 = fq[2];
            
            if(lhshown.find(fid)!=lhshown.end())
            {
                elLh = lhshown[fid];
                elRh = adjElements[fid];
            }
            
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            //ienOUT[curElID] = Elvrts
            for(int u=0;u<ienOUT[elLh].size();u++)
            {
                
                Vijk->x = Vijk->x + vertOUT[(ienOUT[elLh][u]-1)*3];
                Vijk->y = Vijk->y + vertOUT[(ienOUT[elLh][u]-1)*3+1];
                Vijk->z = Vijk->z + vertOUT[(ienOUT[elLh][u]-1)*3+2];
                
            }
            
            Vijk->x = Vijk->x/ienOUT[elLh].size();
            Vijk->y = Vijk->y/ienOUT[elLh].size();
            Vijk->z = Vijk->z/ienOUT[elLh].size();
            
            
            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            
            if(orient0 < 0.0)
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv2);
                ifnOUT->setVal(ftot,3,gv1);
                ifnOUT->setVal(ftot,4,0);
            }
            else
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv1);
                ifnOUT->setVal(ftot,3,gv2);
                ifnOUT->setVal(ftot,4,0);
            }
            
            ifnOUT->setVal(ftot,5,elRh);
            ifnOUT->setVal(ftot,6,elLh);
            ifnOUT->setVal(ftot,7,2);
	    
	    if(elRh == 0)
            {    
                std::cout <<"Found the face in shared -> " << fid << " " << world_rank << " " << elLh << std::endl;
            } 

            ftot++;
        }
        
        
        std::map<int,int> shelltag2glob;
        int alhere      = 0;
        int nothere3    = 0;
        int tetshe      = 0;
        int sna         = 0;
        int notsh       = 0;
        int tellie      = 0;
        int ormin       = 0;
        int tellie2     = 0;
        for(itm=fmShell.begin();itm!=fmShell.end();itm++)
        {
            fid             = itm->first;
            updateID[fid]   = ftot;
            ifnOUT->setVal(ftot,0,3);
            flag = -1;
            
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            
            std::vector<Vert*> Vfaces;
            std::vector<int> reference(3);
            for(int q=0;q<3;q++)
            {
                int lvert = glob2locVid[itm->second[q]];
                Vert* Vf = new Vert;
                Vf->x  = vertOUT[(lvert-1)*3+0];
                Vf->y  = vertOUT[(lvert-1)*3+1];
                Vf->z  = vertOUT[(lvert-1)*3+2];
                
                reference[q] = refOUT[(lvert-1)];
                
                VcF->x = VcF->x + Vf->x;
                VcF->y = VcF->y + Vf->y;
                VcF->z = VcF->z + Vf->z;
                
                Vfaces.push_back(Vf);
                
//                if(reference[q] == 86)
//                {
//                    int pri = glob2locVid[itm->second[q]];
//                    std::cout << "HuH ! -> " << fq[q] << " " << pri << " " << std::setprecision(16) << " ("<< Vf->x << ", " << Vf->y << ", " << Vf->z << ") " << std::endl;
//                }
                
                if(LocationSharedVert_update.find(itm->second[q])!=LocationSharedVert_update.end())
                {
                    fq[q] = LocationSharedVert_update[itm->second[q]];
                    
                    int shelltag = tag2shelltag_glob[itm->second[q]];
                    
                    
                    
                    if(shelltag2glob.find(shelltag)==shelltag2glob.end())
                    {
                        shelltag2glob[shelltag] = fq[q];
                    }
                                        
                    flag = q;
                    
                }
                else
                {
                    fq[q] = tag2glob[itm->second[q]];
                    
                    int shelltag = tag2shelltag_glob[itm->second[q]];
                    
                    if(shelltag2glob.find(shelltag)==shelltag2glob.end())
                    {
                        shelltag2glob[shelltag] = fq[q];
                    }
                    else
                    {
                        alhere++;
                    }
                }
            }
            
            VcF->x = VcF->x/3.0;
            VcF->y = VcF->y/3.0;
            VcF->z = VcF->z/3.0;
            
            
            int gv0 = fq[0];
            int gv1 = fq[1];
            int gv2 = fq[2];
            
            
            if(rhshell.find(fid)!=rhshell.end())
            {
                elRh = rhshell[fid];
                elLh = lhshell[fid];
		
            }
            else
            {
                elRh = rh[fid];
                elLh = lhshell[fid];
                std::cout << "NotFound! -> " << rh[fid] << " " << lhshell[fid] << " " << rhshell.size() << " " << fmShell.size() << " " << std::endl;
                for(int b=0;b<3;b++)
                {
                    int lvert = glob2locVid[itm->second[b]];
                	std::cout << " ("  << reference[b] << ", " <<  fq[b] << ") ";
                	if(reference[b] == 86)
                	{

                		std::cout << std::setprecision(16) << "not correct " << " ("<< vertOUT[(lvert-1)*3+0] << ", " << vertOUT[(lvert-1)*3+1] << ", " << vertOUT[(lvert-1)*3+2] << ") " << std::endl;
                	}
                    else
                    {
                        std::cout << std::setprecision(16) << "sup correct " << " ("<< vertOUT[(lvert-1)*3+0] << ", " << vertOUT[(lvert-1)*3+1] << ", " << vertOUT[(lvert-1)*3+2] << ") " << std::endl;
                    }

                }
                std::cout << std::endl;
                std::cout << "NotFound!" << rh[fid] << " " << lhshell[fid] << " " << rhshell.size() << " " << fmShell.size() << std::endl;
                notsh++;
            }
            
            
            Vert* Vijk   = new Vert;
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            
            // compute element center;
            //ienOUT[curElID] = Elvrts
            for(int u=0;u<ienOUT[elLh].size();u++)
            {
                Vijk->x = Vijk->x + vertOUT[(ienOUT[elLh][u]-1)*3];
                Vijk->y = Vijk->y + vertOUT[(ienOUT[elLh][u]-1)*3+1];
                Vijk->z = Vijk->z + vertOUT[(ienOUT[elLh][u]-1)*3+2];
            }
            
            Vijk->x = Vijk->x/ienOUT[elLh].size();
            Vijk->y = Vijk->y/ienOUT[elLh].size();
            Vijk->z = Vijk->z/ienOUT[elLh].size();
            
            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            
            
            if(orient0 < 0.0)
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv2);
                ifnOUT->setVal(ftot,3,gv1);
                ifnOUT->setVal(ftot,4,0);
            }
            else
            {
                ifnOUT->setVal(ftot,1,gv0);
                ifnOUT->setVal(ftot,2,gv1);
                ifnOUT->setVal(ftot,3,gv2);
                ifnOUT->setVal(ftot,4,0);
            }
            
            ifnOUT->setVal(ftot,5,elRh);
            ifnOUT->setVal(ftot,6,elLh);
            ifnOUT->setVal(ftot,7,2);
            if(elRh == 0)
            {     
                std::cout <<"Found the face in shell -> " << fid << " " << world_rank << " " << elLh << std::endl;
            } 

            Vfaces.clear();

            ftot++;
        }
      
        if(notsh != 0)
	{
		std::cout << world_rank << " found notsh " << std::endl;
	} 
        
        std::map<int,int> shelltag2glob_glob = AllGatherMap(shelltag2glob,comm);
        
        int mapSizeLoc = shelltag2glob.size();
        DistributedParallelState* distrimap = new DistributedParallelState(mapSizeLoc,comm);
        int mapSizeTot = distrimap->getNel();
        int* shelltag_loc = new int[mapSizeLoc];
        int* tag_loc = new int[mapSizeLoc];
        int* shelltag_tot = new int[mapSizeTot];
        int* tag_tot = new int[mapSizeTot];

        int i = 0;
        
        std::map<int,int>::iterator itred;
        
        for(itred=shelltag2glob.begin();itred!=shelltag2glob.end();itred++)
        {
            shelltag_loc[i] = itred->first;
            tag_loc[i] = itred->second;
            i++;
        }
        
        int* offsets = distrimap->getOffsets();
        int* nlocs   = distrimap->getNlocs();
        
        
        MPI_Allgatherv(shelltag_loc,
                       mapSizeLoc,
                       MPI_INT,
                       shelltag_tot,
                       nlocs,
                       offsets,
                       MPI_INT, comm);
        
        
        MPI_Allgatherv(tag_loc,
                       mapSizeLoc,
                       MPI_INT,
                       tag_tot,
                       nlocs,
                       offsets,
                       MPI_INT, comm);
        
        int key,val;
        std::map<int,int> shelltag2glob_global;
        for(int i=0;i<mapSizeTot;i++)
        {
            key = tag_tot[i];
            val = shelltag_tot[i];
            
            if(shelltag2glob_global.find(val)==shelltag2glob_global.end())
            {
                shelltag2glob_global[val] = key;
            }
        }
        
        
        std::map<int,std::vector<int> >::iterator prit;
        int pid   = 0;
        int pfid  = 0;
        int f13   = 0;
        int fptot = 0; //   IntFprims_offsets[world_rank];
        int llvid = 0;
        std::map<int,int> fftell;
        int ftell = 0;
        int gftel = 0;
        int lbb = 1;
        std::map<int,int> Tag2globPrimsVid;
        std::map<int,int> allPnodes;
        int ggvid;
        
        int cc       = 0;
        int notFound = 0;
        int tagnew   = 0;
        std::map<int,int> t2g_tmp;
        std::map<int,int> testmap;
        int notfo  = 0;
        int fo     = 0;
        int notfo2 = 0;
        int redflag = 0;
    
        std::map<int,int> collect;
        std::map<int,int> collect2;
        int cn = 0;
        int fshell = 0;
        
        Array<int>* ifnOUT_prism = new Array<int>(int_face2node_prism.size()+shared_face2node_prism.size(),8);
        
        for( prit=int_face2node_prism.begin();prit!=int_face2node_prism.end();prit++)
        {
            int gfid    = prit->first;
            int npf     = prit->second.size();
            std::vector<int> fce(npf);
            ifnOUT_prism->setVal(fptot,0,npf);
            
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            fshell = 0;
            for(int g=0;g<npf;g++)
            {
                int oldtag = prit->second[g];
                
                if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end() &&
                   shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = LocVerts[lvp]->x;
                    Vf->y = LocVerts[lvp]->y;
                    Vf->z = LocVerts[lvp]->z;
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = tag2glob_prism[oldtag];
                    fce[g]     = globid;
                
                }
                else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end() &&
                        shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = LocVerts[lvp]->x;
                    Vf->y = LocVerts[lvp]->y;
                    Vf->z = LocVerts[lvp]->z;
                    
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = SharedVertsNotOwned[oldtag];
                    fce[g]     = globid;
                    
                }
                else if(shellvertOriginalTag2ref_Glob.find(oldtag)!=shellvertOriginalTag2ref_Glob.end())
                {
                    int ref     = shellvertOriginalTag2ref_Glob[oldtag];
                    int globid  = shelltag2glob_global[ref];
                    fce[g]      = globid;
                    fshell = 1;
                    if(globid>nPrismVerts_tmp)
                    {
                        if(collect2.find(globid)==collect2.end())
                        {
                            collect2[globid] = cn;
                            cn++;
                        }
                    }
                    if(collect.find(globid)==collect.end())
                    {
                        collect[globid] = g;
                    }
                }
            }
            
            VcF->x = VcF->x/npf;
            VcF->y = VcF->y/npf;
            VcF->z = VcF->z/npf;
            
//            if(npf == 3)
//            {
//                ifnOUT_prism->setVal(fptot,4,0);
//            }
            
            if(rhp[gfid] == 0 || lhp[gfid] == 0)
            {
                std::cout << world_rank << " rhp[gfid] and lhp[gfid] are zero " << std::endl;
            }
            
            int leftEl  = lhp[gfid];
            int leftTag = gE2tagE[leftEl];
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            int nvp = prisms[leftTag].size();

            for(int q=0;q<nvp;q++)
            {
                int tag  = prisms[leftTag][q];
                int lvp  = tag2locV_map[tag];

                Vijk->x = Vijk->x + LocVerts[lvp]->x;
                Vijk->y = Vijk->y + LocVerts[lvp]->y;
                Vijk->z = Vijk->z + LocVerts[lvp]->z;
            }

            Vijk->x = Vijk->x/nvp;
            Vijk->y = Vijk->y/nvp;
            Vijk->z = Vijk->z/nvp;

            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
            
            if(orient0 < 0.0)
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[2]);
                    ifnOUT_prism->setVal(fptot,3,fce[1]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[3]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[1]);
                }
                
            }
            else
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[3]);
                }
            }

            
            ifnOUT_prism->setVal(fptot,5,rhp[gfid]);
            ifnOUT_prism->setVal(fptot,6,lhp[gfid]);
            ifnOUT_prism->setVal(fptot,7,2);
            
            if(rhp[gfid] == 0)
            {     
                std::cout <<"Found the face -> prism " << gfid << " " << world_rank << " " << lhp[gfid] << std::endl;
            } 
            fptot++;
            pid++;
        }
    

        int cnttt     = 0;
        int notany    = 0;
        for( prit=shared_face2node_prism.begin();prit!=shared_face2node_prism.end();prit++)
        {
            int gfid  = prit->first;
            int npf   = prit->second.size();
            
            ifnOUT_prism->setVal(fptot,0,npf);
            std::vector<int> fce(npf);
            std::vector<Vert*> Vfaces;
            VcF->x = 0.0;
            VcF->y = 0.0;
            VcF->z = 0.0;
            for(int g=0;g<npf;g++)
            {
                int oldtag = prit->second[g];
                if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end() &&
                   shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = LocVerts[lvp]->x;
                    Vf->y = LocVerts[lvp]->y;
                    Vf->z = LocVerts[lvp]->z;
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = tag2glob_prism[oldtag];
                    fce[g]     = globid;
                    
                }
                else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end() &&
                        shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
                {
                    int lvp  = tag2locV_map[oldtag];
                    Vert* Vf = new Vert;
                    Vf->x = LocVerts[lvp]->x;
                    Vf->y = LocVerts[lvp]->y;
                    Vf->z = LocVerts[lvp]->z;
                    VcF->x = VcF->x + Vf->x;
                    VcF->y = VcF->y + Vf->y;
                    VcF->z = VcF->z + Vf->z;
                    
                    Vfaces.push_back(Vf);
                    int globid = SharedVertsNotOwned[oldtag];
                    fce[g]     = globid;
                    
                }
                else if(shellvertOriginalTag2ref_Glob.find(oldtag)!=shellvertOriginalTag2ref_Glob.end())
                {
                    int ref     = shellvertOriginalTag2ref_Glob[oldtag];
                    int globid = shelltag2glob_global[ref];
                    fce[g]     = globid;
                    if(globid>nPrismVerts_tmp)
                    {
                        if(collect2.find(globid)==collect2.end())
                        {
                            collect2[globid] = cn;
                            cn++;
                        }
                    }
                    if(collect.find(globid)==collect.end())
                    {
                        collect[globid] = g;
                    }
                }
            }
            
            VcF->x = VcF->x/npf;
            VcF->y = VcF->y/npf;
            VcF->z = VcF->z/npf;
            if(npf == 3)
            {
                ifnOUT_prism->setVal(fptot,4,0);
            }
            
            int leftEl  = lhp[gfid];
            int leftTag = gE2tagE[leftEl];
            
            Vijk->x = 0.0;
            Vijk->y = 0.0;
            Vijk->z = 0.0;
            // compute element center;
            int nvp = prisms[leftTag].size();

            for(int q=0;q<nvp;q++)
            {
                int tag  = prisms[leftTag][q];
                int lvp  = tag2locV_map[tag];

                Vijk->x = Vijk->x + LocVerts[lvp]->x;
                Vijk->y = Vijk->y + LocVerts[lvp]->y;
                Vijk->z = Vijk->z + LocVerts[lvp]->z;
            }

            Vijk->x = Vijk->x/nvp;
            Vijk->y = Vijk->y/nvp;
            Vijk->z = Vijk->z/nvp;

            double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
//
            if(orient0 < 0.0)
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[2]);
                    ifnOUT_prism->setVal(fptot,3,fce[1]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[3]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[1]);
                }
            }
            else
            {
                if(npf == 3)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,0);
                }
                if(npf == 4)
                {
                    ifnOUT_prism->setVal(fptot,1,fce[0]);
                    ifnOUT_prism->setVal(fptot,2,fce[1]);
                    ifnOUT_prism->setVal(fptot,3,fce[2]);
                    ifnOUT_prism->setVal(fptot,4,fce[3]);
                }
            }
            
            ifnOUT_prism->setVal(fptot,5,rhp[gfid]);
            ifnOUT_prism->setVal(fptot,6,lhp[gfid]);
            ifnOUT_prism->setVal(fptot,7,2);
            
	    if(rhp[gfid] == 0)
            {
                std::cout <<"Found the face -> prism shared " << gfid << " " << world_rank << " " << lhp[gfid] << std::endl;
            }

            fptot++;
            pid++;
        }
        

        
        // End reducing the map.
        
        DistributedParallelState* distftot   = new DistributedParallelState(ftot,comm);
        DistributedParallelState* distfptot  = new DistributedParallelState(fptot,comm);

        int nTotInteriorFaces          = distftot->getNel();
        int* TotIntFaces_offsets       = distftot->getOffsets();
        
        int nTotInteriorFaces_prism    = distfptot->getNel();
        int* TotIntFaces_offsets_prism = distfptot->getOffsets();
        
        std::map<int,std::vector<int> >::iterator bit;
        std::set<int> sorted_BCid;
        std::map<int,int> sorted_BCid_map;
        std::map<int,int> sorted_NBCid_map;
        int ii = 0;
        int Nbf;
        int nloc_bcs_p = pbcmap.size();
        std::set<int> bcids_tot;
        
        std::vector<int> Lbcs;
        
        for(bit=bcmap.begin();bit!=bcmap.end();bit++)
        {
            if(bcids_tot.find(bit->first)==bcids_tot.end())
            {
                bcids_tot.insert(bit->first);
                Lbcs.push_back(bit->first);
            }
        }
        
        for(bit=pbcmap.begin();bit!=pbcmap.end();bit++)
		{
			if(bcids_tot.find(bit->first)==bcids_tot.end())
			{
				bcids_tot.insert(bit->first);
				Lbcs.push_back(bit->first);
			}
		}
        
        int nloc_bcs   = Lbcs.size();

        DistributedParallelState* distLocBCs = new DistributedParallelState(nloc_bcs,comm);
        
        int Nt_BCs                 = distLocBCs->getNel();
        int* BCs_offsets           = distLocBCs->getOffsets();
        int* BCs_nlocs             = distLocBCs->getNlocs();
        int* BCs_arr               = new int[Nt_BCs];
        
        MPI_Allgatherv(&Lbcs[0],
                       nloc_bcs,
                       MPI_INT,
                       BCs_arr,
                       BCs_nlocs,
                       BCs_offsets,
                       MPI_INT, comm);

        std::set<int> bcsToT;
        std::map<int,int> bcentry;
        for(int i=0;i<Nt_BCs;i++)
        {
            if(bcsToT.find(BCs_arr[i])==bcsToT.end())
            {
                bcsToT.insert(BCs_arr[i]);
            }
        }
        
        int* bcid = new int[bcsToT.size()];
        int* nlbc = new int[bcsToT.size()];
        std::set<int>::iterator entry;
        int cnt = 0;

        int lac;
        q = 0;
        for(entry=bcsToT.begin();entry!=bcsToT.end();entry++)
        {
            int ee = *entry;
            int Nbft = 0;
            int Nbfp = 0;
            
            if(bcmap.find(ee)!=bcmap.end())
            {
                Nbft = bcmap[ee].size();
            }
            if(pbcmap.find(ee)!=pbcmap.end())
            {
            	Nbfp = pbcmap[ee].size();
            }

            bcid[q]  = ee;
            nlbc[q]  = Nbft+Nbfp;
            q++;
        }
        
        
        
        std::vector<Array<int>* > bcArrays;
        std::map<int,int> bcsizing;
        std::vector<int> bci_offsets;
        std::vector<int> bciTot_offsets;
        int nTotBCFaces = 0;
        int nTotBCFaces_offset = 0;
    
        int failbc = 0;
        int globalVid;
        

        for(int i=0;i<bcsToT.size();i++)
        {
            int bc_id = bcid[i];
            DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);
            int NelLoc_bci = nlbc[i];
            int NelTot_bci = distBCi->getNel();
            
            Array<int>* ifn_bc_i = new Array<int>(NelLoc_bci,8);
            int offsetbci        = distBCi->getOffsets()[world_rank];
            int fbc  = 0;
            int Nbft = 0;
			int Nbfp = 0;

            if(bcmap.find(bc_id)!=bcmap.end())
            {
                Nbft = bcmap[bc_id].size();
            }
            if(pbcmap.find(bc_id)!=pbcmap.end())
            {
                Nbfp = pbcmap[bc_id].size();
            }
            
            if(Nbfp!=0)
            {
                int sk = 0;
                
                for(int q=0;q<Nbfp;q++)
				{
					int bcface = pbcmap[bc_id][q];
					int flag = -1;
					
					int nppf = bc_face2node_prism[bcface].size();
					
					ifn_bc_i->setVal(fbc,0,nppf);
                    std::vector<int> face_tmp(nppf);
                    std::vector<int> fce(nppf);
                    std::vector<Vert*> Vfaces;
                    VcF->x = 0.0;
                    VcF->y = 0.0;
                    VcF->z = 0.0;
                    
                    for(int g=0;g<nppf;g++)
                    {
                        int oldtag = bc_face2node_prism[bcface][g];
                        
                        if(tag2glob_prism.find(oldtag)!=tag2glob_prism.end() &&
                           shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
                        {
                            int lvp  = tag2locV_map[oldtag];
                            Vert* Vf = new Vert;
                            Vf->x = LocVerts[lvp]->x;
                            Vf->y = LocVerts[lvp]->y;
                            Vf->z = LocVerts[lvp]->z;
                            VcF->x = VcF->x + Vf->x;
                            VcF->y = VcF->y + Vf->y;
                            VcF->z = VcF->z + Vf->z;
                            
                            Vfaces.push_back(Vf);
                            int globid = tag2glob_prism[oldtag];
                            fce[g]     = globid;
                        }
                        else if(SharedVertsNotOwned.find(oldtag)!=SharedVertsNotOwned.end() &&
                                shellvertOriginalTag2ref_Glob.find(oldtag)==shellvertOriginalTag2ref_Glob.end())
                        {
                            int lvp  = tag2locV_map[oldtag];
                            Vert* Vf = new Vert;
                            Vf->x = LocVerts[lvp]->x;
                            Vf->y = LocVerts[lvp]->y;
                            Vf->z = LocVerts[lvp]->z;
                            VcF->x = VcF->x + Vf->x;
                            VcF->y = VcF->y + Vf->y;
                            VcF->z = VcF->z + Vf->z;
                            
                            Vfaces.push_back(Vf);
                            int globid = SharedVertsNotOwned[oldtag];
                            fce[g]     = globid;
                        }
                        else if(shellvertOriginalTag2ref_Glob.find(oldtag)!=shellvertOriginalTag2ref_Glob.end())
                        {
                            int ref     = shellvertOriginalTag2ref_Glob[oldtag];
                            int globid  = shelltag2glob_global[ref];
                            fce[g]      = globid;
                            if(globid>nPrismVerts_tmp)
                            {
                                if(collect2.find(globid)==collect2.end())
                                {
                                    collect2[globid] = cn;
                                    cn++;
                                }
                            }
                            if(collect.find(globid)==collect.end())
                            {
                                collect[globid] = g;
                            }
                        }
                    }
                    
                    VcF->x = VcF->x/nppf;
                    VcF->y = VcF->y/nppf;
                    VcF->z = VcF->z/nppf;
                    
                    int leftEl  = lhp[bcface];
                    int leftTag = gE2tagE[leftEl];
                    
                    Vijk->x = 0.0;
                    Vijk->y = 0.0;
                    Vijk->z = 0.0;
                    // compute element center;
                    int nvp = prisms[leftTag].size();

                    for(int q=0;q<nvp;q++)
                    {
                        int tag  = prisms[leftTag][q];
                        int lvp  = tag2locV_map[tag];

                        Vijk->x = Vijk->x + LocVerts[lvp]->x;
                        Vijk->y = Vijk->y + LocVerts[lvp]->y;
                        Vijk->z = Vijk->z + LocVerts[lvp]->z;
                    }

                    Vijk->x = Vijk->x/nvp;
                    Vijk->y = Vijk->y/nvp;
                    Vijk->z = Vijk->z/nvp;

                    double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);
                    
                    if(orient0 < 0.0)
                    {
                        if(nppf == 3)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[2]);
                            ifn_bc_i->setVal(fbc,3,fce[1]);
                            ifn_bc_i->setVal(fbc,4,0);
                        }
                        if(nppf == 4)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[3]);
                            ifn_bc_i->setVal(fbc,3,fce[2]);
                            ifn_bc_i->setVal(fbc,4,fce[1]);
                        }
                        
                    }
                    else
                    {
                        if(nppf == 3)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[1]);
                            ifn_bc_i->setVal(fbc,3,fce[2]);
                            ifn_bc_i->setVal(fbc,4,0);
                        }
                        if(nppf == 4)
                        {
                            ifn_bc_i->setVal(fbc,1,fce[0]);
                            ifn_bc_i->setVal(fbc,2,fce[1]);
                            ifn_bc_i->setVal(fbc,3,fce[2]);
                            ifn_bc_i->setVal(fbc,4,fce[3]);
                        }
                    }
                    
					ifn_bc_i->setVal(fbc,5,0);
					ifn_bc_i->setVal(fbc,6,lhp[bcface]);
					ifn_bc_i->setVal(fbc,7,bc_id);
					
					fbc++;
				}
                
                
            }
            
            if(Nbft!=0 )
            {
                for(int q=0;q<Nbft;q++)
                {
                    int bcface = bcmap[bc_id][q];
                    ifn_bc_i->setVal(fbc,0,3);

                    std::vector<Vert*> Vfaces;
                    
                    VcF->x = 0.0;
                    VcF->y = 0.0;
                    VcF->z = 0.0;
                    flag   = -1;

                    for(int w=0;w<3;w++)
                    {
                        int vertexID = BoundaryFaces[bcface]->GetEdgeIDs()[w];
                        int lvert = glob2locVid[vertexID];
                        
                        Vert* Vf = new Vert;
                        Vf->x  = vertOUT[(lvert-1)*3+0];
                        Vf->y  = vertOUT[(lvert-1)*3+1];
                        Vf->z  = vertOUT[(lvert-1)*3+2];
                        
                        VcF->x = VcF->x + Vf->x;
                        VcF->y = VcF->y + Vf->y;
                        VcF->z = VcF->z + Vf->z;
                        
                        Vfaces.push_back(Vf);
                        
                        if(LocationSharedVert_update.find(vertexID)!=LocationSharedVert_update.end())
                        {
                            fq[w] = LocationSharedVert_update[vertexID];
                            ifn_bc_i->setVal(fbc,w+1,LocationSharedVert_update[vertexID]);
                            if(ifn_bc_i->getVal(fbc,w+1)==0)
                            {
                                std::cout << " TET Boundary V zero in the tet local LocationSharedVert_update " << world_rank << " " << fbc << std::endl;
                            }
                            flag = w;
                        }
                        else
                        {
                            fq[w] = tag2glob[vertexID];
                            ifn_bc_i->setVal(fbc,w+1,tag2glob[vertexID]);
                            if(ifn_bc_i->getVal(fbc,w+1)==0)
                            {
                                std::cout << " TET Boundary V zero in the tet local tag2glob " << world_rank << " " << fbc << std::endl;
                            }
                        }
                    }
                    
                    //std::cout << std::endl;
                    int elLh = lhbnd[bcface];
                    
                    Vijk->x = 0.0;
                    Vijk->y = 0.0;
                    Vijk->z = 0.0;
                    // compute element center;
                    //ienOUT[curElID] = Elvrts
                    for(int u=0;u<ienOUT[elLh].size();u++)
                    {
                        
                        Vijk->x = Vijk->x + vertOUT[(ienOUT[elLh][u]-1)*3];
                        Vijk->y = Vijk->y + vertOUT[(ienOUT[elLh][u]-1)*3+1];
                        Vijk->z = Vijk->z + vertOUT[(ienOUT[elLh][u]-1)*3+2];
                        
                    }
                    
                    Vijk->x = Vijk->x/ienOUT[elLh].size();
                    Vijk->y = Vijk->y/ienOUT[elLh].size();
                    Vijk->z = Vijk->z/ienOUT[elLh].size();
                    
                    double orient0 = CheckFaceOrientation(VcF,Vfaces,Vijk);

                    if(orient0 < 0.0)
                    {
                        std::cout << "Weve got negative faces " << orient0 << " " << bc_id << std::endl;
                    }
                    
                    Vfaces.clear();
                    
                    ifn_bc_i->setVal(fbc,4,0);
                    ifn_bc_i->setVal(fbc,5,0);
                    ifn_bc_i->setVal(fbc,6,lhbnd[bcface]);
                    ifn_bc_i->setVal(fbc,7,bc_id);
                    
                    fbc++;
                }
            }
            
            
            int nbt = Nbft+Nbfp;
            
            bcsizing[bc_id] = NelTot_bci;
            bci_offsets.push_back(offsetbci);
            bciTot_offsets.push_back(nTotBCFaces_offset);
            bcArrays.push_back(ifn_bc_i);
            
            nTotBCFaces_offset = nTotBCFaces_offset + NelTot_bci;
            nTotBCFaces        = nTotBCFaces + NelTot_bci;
            
//            if(world_rank == 0 && bc_id == 3)
//            {
//                std::cout << "bcID and nBCFaces -> " << bc_id << " " << NelTot_bci << " " << pmmg_shellfacedist->getNel() << " " << PMMG_Shell2Prismdist->getNel() << " " << sucDist->getNel() << std::endl;
//            }
            
        }

        int nPrismOUT = parmmg_iet_prisms->getNrow();

        DistributedParallelState* distTetraVerts = new DistributedParallelState(xcn_parmmg->getNrow(),comm);
        DistributedParallelState* distPrismVerts = new DistributedParallelState(xcn_prisms_int->getNrow()+xcn_prisms_shared->getNrow(),comm);
        
        
        DistributedParallelState* distTetra      = new DistributedParallelState(nTetrahedraOUT,comm);
        DistributedParallelState* distPrism      = new DistributedParallelState(nPrismOUT,comm);
        int ToTElements_prism           = distPrism->getNel();
        int ToTElements                 = distTetra->getNel();

        int ToTElements_offset_prism    = distPrism->getOffsets()[world_rank];
        int ToTElements_offset          = distTetra->getOffsets()[world_rank];
        int nTotTetraVerts_v2           = distTetraVerts->getNel();

        int nTotElements                = ToTElements_prism+ToTElements;
        int nTotFaces                   = nTotInteriorFaces_prism + nTotInteriorFaces + nTotBCFaces;
        int nTotIntFaces                = nTotInteriorFaces_prism + nTotInteriorFaces;
        
        int nTotPrismVerts_v2           = distPrismVerts->getNel();
        int nTotPrismIntVerts_v2        = distPrismIntVerts->getNel();
        int nTotPrismShaVerts_v2        = distPrismShaVerts->getNel();
    
        int TotPrismVerts_offset_int    = distPrismIntVerts->getOffsets()[world_rank];
        int TotPrismVerts_offset_sha    = distPrismShaVerts->getOffsets()[world_rank];

        int TotPrismVerts_offset        = distPrismVerts->getOffsets()[world_rank];
        
        int nTotVertsPrismTetra = nTotPrismVerts_v2+nTotTetraVerts_v2;
        int nbo = bcArrays.size();
        //std::cout << "-- Constructing the zdefs array..."<<std::endl;
        Array<int>* adapt_zdefs = new Array<int>(3+nbo,7);
        // Collect node data (10) . Starting index-ending index Nodes
        adapt_zdefs->setVal(0,0,10);
        adapt_zdefs->setVal(0,1,-1);
        adapt_zdefs->setVal(0,2,1);
        adapt_zdefs->setVal(0,3,1);
        adapt_zdefs->setVal(0,4,nTotVertsPrismTetra);
        adapt_zdefs->setVal(0,5,us3d->zdefs->getVal(0,5));
        adapt_zdefs->setVal(0,6,us3d->zdefs->getVal(0,6));
        // Collect element data (12) . Starting index-ending index Element
        adapt_zdefs->setVal(1,0,12);
        adapt_zdefs->setVal(1,1,-1);
        adapt_zdefs->setVal(1,2,2);
        adapt_zdefs->setVal(1,3,1);
        adapt_zdefs->setVal(1,4,nTotElements);
        adapt_zdefs->setVal(1,5,us3d->zdefs->getVal(1,5));
        adapt_zdefs->setVal(1,6,2);
        // Collect internal face data (13) . Starting index-ending index internal face.
        adapt_zdefs->setVal(2,0,13);
        adapt_zdefs->setVal(2,1,-1);
        adapt_zdefs->setVal(2,2, 3);
        adapt_zdefs->setVal(2,3, 1);
        adapt_zdefs->setVal(2,4,nTotIntFaces);
        adapt_zdefs->setVal(2,5,us3d->zdefs->getVal(2,5));
        adapt_zdefs->setVal(2,6,2);
        
        int qq  = 1;
        int nb  = 0;
        int face_start = nTotIntFaces+1;
        int face_end;
        std::map<int,int>::iterator itr;
        for(itr=bcsizing.begin();itr!=bcsizing.end();itr++)
        {
        	
            int bnd_ref  = itr->first;
            int bnd_size = itr->second;
            face_end = face_start+bnd_size-1;
            adapt_zdefs->setVal(3+nb,0,13);
            adapt_zdefs->setVal(3+nb,1,-1);
            adapt_zdefs->setVal(3+nb,2,3+qq);
            adapt_zdefs->setVal(3+nb,3,face_start);
            adapt_zdefs->setVal(3+nb,4,face_end);
            adapt_zdefs->setVal(3+nb,5,bnd_ref);
            adapt_zdefs->setVal(3+nb,6,2);
            face_start = face_end+1;

            nb++;
            qq++;
        }
        
        
        
        if(world_rank == 0)
        {
            PlotBoundaryData(us3d->znames,adapt_zdefs);
        }
        

        
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        //===================================================================================
        
        hid_t ret;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, comm, info);
        hid_t file_id = H5Fcreate("grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

        hid_t    dset_id;
        hid_t    filespace;
        hid_t    memspace;
        hid_t    status;
        hsize_t     dimsf[2];
        hsize_t     countH5[2];
        hsize_t     offsetH5[2];
        
        hsize_t dimsf_att = 1;
        hid_t att_space = H5Screate_simple(1, &dimsf_att, NULL);
        hid_t type =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type, 14);
        ret = H5Tset_strpad(type,H5T_STR_SPACEPAD);
        hid_t attr_id   = H5Acreate (file_id, "filetype", type, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri[] = "US3D Grid File";
        status = H5Awrite(attr_id, type, &stri);
        H5Aclose(attr_id);
        
        hsize_t dimsf_att2 = 1;
         att_space = H5Screate_simple(1, &dimsf_att2, NULL);
        hid_t type2 =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type2, 5);
        ret = H5Tset_strpad(type2,H5T_STR_SPACEPAD);
        attr_id   = H5Acreate (file_id, "filevers", type2, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri2[] = "1.1.8";
        status = H5Awrite(attr_id, type2, &stri2);
        H5Aclose(attr_id);
        
        hid_t group_info_id  = H5Gcreate(file_id, "info", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

        hsize_t dimsf_att3 = 1;
        att_space = H5Screate_simple(1, &dimsf_att3, NULL);
        hid_t type3 =  H5Tcopy (H5T_C_S1);
        ret = H5Tset_size (type3, 10);
        ret = H5Tset_strpad(type3,H5T_STR_SPACEPAD);
        attr_id   = H5Acreate (group_info_id, "date", type3, att_space, H5P_DEFAULT, H5P_DEFAULT);
        char stri3[] = "27-05-1987";
        status = H5Awrite(attr_id, type3, &stri3);
        H5Aclose(attr_id);
        
        hid_t group_grid_id  = H5Gcreate(group_info_id, "grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        dimsf_att = 1;
        att_space = H5Screate_simple(1, &dimsf_att, NULL);
        attr_id   = H5Acreate (group_grid_id, "nc", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        int value = nTotElements;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "nf", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotFaces;// nTotInteriorFaces_prism+nTotInteriorFaces+nTotBCFaces;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "ng", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotBCFaces;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        attr_id   = H5Acreate (group_grid_id, "nn", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = nTotPrismVerts_v2+nTotTetraVerts_v2;//ToTVrts;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        //std::cout << " NN " << nTotVertsPrismTetra << std::endl;
        //====================================================================================
        // Add iet map to the grid.h5 file
        //====================================================================================
        dimsf[0] = nTotElements;
        dimsf[1] = parmmg_iet_prisms->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);

        dset_id = H5Dcreate(file_id, "iet",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = parmmg_iet_prisms->getNrow();
        countH5[1]  = parmmg_iet_prisms->getNcol();
        
        offsetH5[0] = ToTElements_offset_prism;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, parmmg_iet_prisms->data);
        delete parmmg_iet_prisms;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        countH5[0]  = parmmg_iet->getNrow();
        countH5[1]  = parmmg_iet->getNcol();

        offsetH5[0] = ToTElements_prism+ToTElements_offset;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        //filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, parmmg_iet->data);
        delete parmmg_iet;
        //====================================================================================
        // Add xcn map to the grid.h5 file
        //====================================================================================
        

        
        
        dimsf[0] = nTotPrismIntVerts_v2+nTotPrismShaVerts_v2+nTotTetraVerts_v2;
        dimsf[1] = xcn_prisms_int->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);
        
        dset_id = H5Dcreate(file_id, "xcn",
                            H5T_NATIVE_DOUBLE, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = xcn_prisms_int->getNrow();
        countH5[1]  = xcn_prisms_int->getNcol();
        
        offsetH5[0] = TotPrismVerts_offset_int;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_prisms_int->data);
        delete xcn_prisms_int;
    
//        dimsf[0] = nTotPrismShaVerts_v2;//+nTotTetraVerts_v2;
//        dimsf[1] = xcn_prisms_shared->getNcol();
//        filespace = H5Screate_simple(2, dimsf, NULL);
        
//        dset_id = H5Dcreate(file_id, "xcn",
//                            H5T_NATIVE_DOUBLE, filespace,
//                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        countH5[0]  = xcn_prisms_shared->getNrow();
        countH5[1]  = xcn_prisms_shared->getNcol();
        
        
        offsetH5[0] = nTotPrismIntVerts_v2+TotPrismVerts_offset_sha;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);

        filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        
        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_prisms_shared->data);
        //std::cout << "world " << nTotPrismIntVerts_v2 << " " << TotPrismVerts_offset_sha << " " << xcn_prisms_shared->getNrow() << std::endl;
        delete xcn_prisms_shared;
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
        countH5[0]  = xcn_parmmg->getNrow();
        countH5[1]  = xcn_parmmg->getNcol();

        offsetH5[0] = nTotPrismIntVerts_v2+nTotPrismShaVerts_v2+ToTVrts_offset;
        offsetH5[1] = 0;
        memspace = H5Screate_simple(2, countH5, NULL);
        //filespace = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

        status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xcn_parmmg->data);
        delete xcn_parmmg;
        
       
        //===================================================================================
//        int nTotInteriorFaces          = distftot->getNel();
//        int* TotIntFaces_offsets       = distftot->getOffsets();
//
//        int nTotInteriorFaces_prism    = distfptot->getNel();
//        int* TotIntFaces_offsets_prism = distfptot->getOffsets();
        
        dimsf[0]  = nTotInteriorFaces_prism+nTotInteriorFaces+nTotBCFaces;
        dimsf[1]  = ifnOUT_prism->getNcol();
        
        filespace = H5Screate_simple(2, dimsf, NULL);
        dset_id   = H5Dcreate(file_id, "ifn",
                            H5T_NATIVE_INT, filespace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//
        countH5[0]  = ifnOUT_prism->getNrow();
        countH5[1]  = dimsf[1];
        
        offsetH5[0] = TotIntFaces_offsets_prism[world_rank];
        offsetH5[1] = 0;
        
        memspace      = H5Screate_simple(2, countH5, NULL);
        filespace     = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id      = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
//
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, ifnOUT_prism->data);
        
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        countH5[0]  = ifnOUT->getNrow();
        countH5[1]  = dimsf[1];

        offsetH5[0] = nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank];
        offsetH5[1] = 0;
        //std::cout << "nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank ] " << TotIntFaces_offsets_prism[world_rank] << " " << ifnOUT_prism->getNrow() << " " << nTotInteriorFaces_prism << " " << TotIntFaces_offsets[world_rank ] << " " << ifnOUT->getNrow() << world_rank << " -> " <<nTotInteriorFaces_prism << " + " << nTotInteriorFaces << " " << ftot << " " << fptot  << " " << world_rank << " " << nTotInteriorFaces_prism+nTotInteriorFaces << std::endl;

        memspace     = H5Screate_simple(2, countH5, NULL);
        filespace     = H5Dget_space(dset_id);
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);
        plist_id     = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
//
        status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                      plist_id, ifnOUT->data);

        //std::cout << "world " << nTotInteriorFaces_prism+TotIntFaces_offsets[world_rank ]+ifnOUT->getNrow() << " --> " << nTotInteriorFaces_prism+nTotInteriorFaces << " +++ "<< world_rank << std::endl;
        //===================================================================================

        for(int i=0;i<bcsToT.size();i++)
        {
            
            
            int bc_id = bcid[i];
            DistributedParallelState* distBCi = new DistributedParallelState(nlbc[i],comm);

            int NelLoc_bci = nlbc[i];
            int NelTot_bci = distBCi->getNel();

            Array<int>* ifn_bc_i = bcArrays[i];

            countH5[0]    = ifn_bc_i->getNrow();
            countH5[1]    = dimsf[1];
            
            offsetH5[0]   = nTotInteriorFaces_prism+nTotInteriorFaces+bciTot_offsets[i]+bci_offsets[i];
            offsetH5[1]   = 0;
            memspace      = H5Screate_simple(2, countH5, NULL);
            filespace     = H5Dget_space(dset_id);

            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetH5, NULL, countH5, NULL);

            plist_id     = H5Pcreate(H5P_DATASET_XFER);
            H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    //
            status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace,
                          plist_id, ifn_bc_i->data);

        }
        
        // Create group;
        //====================================================================================
        hid_t group_zones_id  = H5Gcreate(file_id, "zones", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        // Add attribute to group:
        //====================================================================================
        dimsf_att = 1;
        att_space = H5Screate_simple(1, &dimsf_att, NULL);
        attr_id   = H5Acreate (group_zones_id, "nz", H5T_STD_I32BE, att_space, H5P_DEFAULT, H5P_DEFAULT);
        value = 3+nbo;
        status = H5Awrite(attr_id, H5T_NATIVE_INT, &value);
        H5Aclose(attr_id);
        //====================================================================================
        dimsf[0] = adapt_zdefs->getNrow();
        dimsf[1] = adapt_zdefs->getNcol();
        filespace = H5Screate_simple(2, dimsf, NULL);
        hid_t dset_zdefs_id = H5Dcreate(group_zones_id, "zdefs", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Sclose(filespace);
        
        countH5[0]  = dimsf[0];
        countH5[1]  = dimsf[1];
        offsetH5[0] = 0;
        offsetH5[1] = 0;
        memspace  = H5Screate_simple(2, countH5, NULL);
        filespace = H5Dget_space(dset_zdefs_id);
        
        //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

        status = H5Dwrite(dset_zdefs_id, H5T_NATIVE_INT, memspace, filespace, plist_id, adapt_zdefs->data);
        //====================================================================================
        
        dimsf_att = us3d->znames->getNrow();
        filespace = H5Screate_simple(1, &dimsf_att, NULL);
        type =  H5Tcopy (H5T_C_S1);
        ret  = H5Tset_size (type, 20);
        ret  = H5Tset_strpad(type, H5T_STR_SPACEPAD);
        hid_t dset_znames_id = H5Dcreate(group_zones_id, "znames", type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Sclose(filespace);

        hsize_t cnt2 = us3d->znames->getNrow();
        memspace  = H5Screate_simple(1, &cnt2, NULL);
        filespace = H5Dget_space(dset_znames_id);

        status = H5Dwrite(dset_znames_id, type, memspace, filespace, plist_id, us3d->znames->data);

        //===================================================================================
        //===================================================================================
        //===================================================================================
    }
    
    clock_t t1_ijk = clock();
    double duration = ( t1_ijk - t0_ijk) / (double) CLOCKS_PER_SEC;
    double dur_max;
    MPI_Allreduce(&duration, &dur_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    
    if(world_rank==0)
    {
        //std::cout << std::setprecision(16) << "Computing the metric takes " << dur_max_met << " seconds using " << world_size << " procs. " << std::endl;
        std::cout << std::setprecision(16) << "Redistributing the tetrahedra takes " << dur_max_redis << " seconds using " << world_size << " procs. " << std::endl;
        std::cout << std::setprecision(16) << "Adapting the tetrahedra takes " << dur_max_adapt << " seconds using " << world_size << " procs. " << std::endl;
        std::cout << std::setprecision(16) << "Writing out the grid takes " << dur_max << " seconds using " << world_size << "procs. " << std::endl;
    	std::cout << "Finalizing process" << std::endl;     
    }
    
    /**/
    
    MPI_Finalize();
    
}

