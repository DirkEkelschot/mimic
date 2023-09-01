#include "adapt.h"
#include "adapt_datastruct.h"

#ifndef ADAPT_PRISMTETRATRACE_H
#define ADAPT_PRISMTETRATRACE_H

class PrismTetraTrace{
    public:
        PrismTetraTrace(){};
        PrismTetraTrace(MPI_Comm comm,
                        std::map<int,std::vector<int> > ife, 
                        std::map<int,int> iet,
                        int Nelem,
                        int Nface,
                        int Nvert);
        ~PrismTetraTrace();
        std::map<int,std::vector<int> > GetTrace();

    private:
        std::map<int,std::vector<int> > trace;
};
#endif
