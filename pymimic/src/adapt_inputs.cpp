#include "adapt_inputs.h"


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





Inputs* ReadXmlFile(MPI_Comm comm, const char* filename)
{


    int world_size;
    MPI_Comm_size(comm, &world_size);
    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(comm, &world_rank);


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
    if(param_map.find("MetricProvided")!=param_map.end())
    {
        inp->MetricProvided = param_map["MetricProvided"];
    }
    else
    {
        std::cout << "Error: MetricProvided is not defined in metric.xml." << std::endl;
    }
    if(param_map.find("RunNumber")!=param_map.end())
    {
        inp->RunNumber = param_map["RunNumber"];
    }
    else
    {
        std::cout << "Error: RunNumber is not defined in metric.xml." << std::endl;
    }




    if(world_rank == 0)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << "============== Metric parameters ==================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "Nproc             = " << world_size << std::endl;
        std::cout << "hgrad             = " << inp->hgrad << std::endl;
        std::cout << "hmin              = " << inp->hmin << std::endl;
        std::cout << "hmax              = " << inp->hmax << std::endl;
        std::cout << "MetScale          = " << inp->MetScale << std::endl;
        std::cout << "Hausdorff         = " << inp->hausd << std::endl;
        std::cout << "NiterPart         = " << inp->niter << std::endl;
        std::cout << "StateVariable     = " << inp->StateVar << std::endl;
        std::cout << "MetricProvided    = " << inp->MetricProvided << std::endl;
        std::cout << "RunNumber         = " << inp->RunNumber << std::endl;


        if(inp->ReadFromStats == 0)
        {
            std::cout << "Reading statistics? -> NO (5th entry in the metric.inp file is set to 0.)" << std::endl;
            std::cout << "The metric is reconstructed based on instantaneous Mach number"<<std::endl;
        }
        if(inp->ReadFromStats == 1)
        {
            std::cout << "Reading statistics? -> YES (5th entry in the metric.inp file is set to 1.)" << std::endl;
            std::cout << "The metric is reconstructed based on the mean of Mach number."<<std::endl;

        }
        if(inp->RunWakRefinement==0)
        {
            std::cout << "Wake refinement is switch OFF. (6th entry in the metric.inp file is set to 0. hwake, the 7th entry defined in the metric.inp file, is being ignored)" << std::endl;
            
        }
        if(inp->RunWakRefinement==1)
        {
            std::cout << "Wake refinement is switch ON with hwake = " << inp->hwake << "(6th entry in the metric.inp file is set to 1 and hwake is set equal to the 7th entry defined in the metric.inp file.) " << std::endl;
        }
        if(inp->StateVar == 0)
	{
	    std::cout << "We are adapting based on the Mach number."<<std::endl;
	}
	if(inp->StateVar == 1)
	{
	    std::cout << "We are adapting based on the static Temperature." << std::endl;
        }
        std::cout << "===================================================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "  " << std::endl;
    }




    return inp;
}