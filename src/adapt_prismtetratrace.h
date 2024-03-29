#include "adapt.h"
#include "adapt_datastruct.h"
#include "NekFace.h"

#ifndef ADAPT_PRISMTETRATRACE_H
#define ADAPT_PRISMTETRATRACE_H

class PrismTetraTrace{
    public:
        PrismTetraTrace(){};
        PrismTetraTrace(MPI_Comm comm,
                        std::vector<int> element2rank,
                        std::map<int,std::vector<int> > ife,
                        std::map<int,std::vector<int> > ifn,  
                        std::map<int,int> iet,
                        int Nelem,
                        int Nface,
                        int Nvert);

        ~PrismTetraTrace();
        std::map<int,int> GetUniqueTraceVerts2RefMap();

    private:
        std::map<int,int> unique_trace_verts;
};
#endif
