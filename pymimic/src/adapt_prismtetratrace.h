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

        FaceSetPointer GetRefTraceFaceSet();
        std::map<int,std::map<int,int> > GetTrace();
        std::map<int,std::vector<int> > GetTraceVerts();
        std::map<int,int> GetUniqueTraceVerts2RefMap();
        std::map<int,std::vector<int> > GetLeftRightElements();
        std::map<int,int> GetTraceRef();

    private:

        FaceSetPointer m_RefTraceFaces;
        std::map<int,int> unique_trace_verts;
        std::map<int,std::map<int,int> > trace_elems;
        std::map<int,std::vector<int> > trace_verts;
        std::map<int,std::vector<int> > trace_LR_elem;
        std::map<int,int> trace_ref;
};
#endif
