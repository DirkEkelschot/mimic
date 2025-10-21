
#include "adapt.h"
#include "adapt_datastruct.h"
#include <iomanip>

void ParseEquals(const std::string &line, std::string &lhs,
                                std::string &rhs);


Inputs* ReadXmlFile(MPI_Comm comm, const char* filename);