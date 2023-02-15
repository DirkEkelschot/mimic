#include "adapt.h"
#include <memory>
#include "adapt_array.h"
#include "adapt_hashutils.h"

class NekFace {
    public:
        NekFace();
        NekFace(std::vector<int> eid);
        ~NekFace();
        int GetFaceID() const;
        std::vector<int> GetEdgeIDs() const;
        void SetFaceID(int idje);
        int GetFaceRef() const;
        void SetFaceRef(int ref);
        std::size_t ComputeFaceHash(std::vector<int> Face) const;
        std::size_t GetFaceHash();
        void SetFaceLeftElement(int le);
        void SetFaceRightElement(int re);
        int GetFaceLeftElement() const;
        int GetFaceRightElement() const;
        bool GetHandled() const;
        void SetHandled(bool H);
        bool operator==(const NekFace& otherFace) const;
//        bool operator<(const NekFace& otherFace);
        
      
     private:
        std::vector<int> m_eId;
        int m_id;
        int m_seed;
        int m_ref;
        int m_le;
        int m_re;
        bool handled;
};

typedef std::shared_ptr<NekFace> FaceSharedPtr;



bool operator==(FaceSharedPtr const &p1,
                FaceSharedPtr const &p2);

//bool operator<(FaceSharedPtr const &p1,
//               FaceSharedPtr const &p2);



struct FaceHash : std::unary_function<NekFace, std::size_t>
{
    std::size_t operator()(const NekFace &vf) const
    {
        std::vector<int> eid = vf.GetEdgeIDs();
        unsigned int nVert = eid.size();
        std::size_t seed = 0;
        std::vector<unsigned int> ids(nVert);

        for (int i = 0; i < nVert; ++i)
        {
            ids[i] = eid[i];
        }
        
        std::sort(ids.begin(), ids.end());
        hash_range(seed, ids.begin(), ids.end());

        return seed;
    }
};


struct FaceHashPointer : std::unary_function<FaceSharedPtr, std::size_t>
{
    std::size_t operator()(const FaceSharedPtr &vf) const
    {
        std::vector<int> eid = vf->GetEdgeIDs();
        unsigned int nVert = eid.size();
        std::size_t seed = 0;
        std::vector<unsigned int> ids(nVert);

        for (int i = 0; i < nVert; ++i)
        {
            ids[i] = eid[i];
        }
        
        std::sort(ids.begin(), ids.end());
        hash_range(seed, ids.begin(), ids.end());

        return seed;
    }
};

typedef std::unordered_set<NekFace, FaceHash> FaceSet;
//typedef std::map<NekFace, FaceHash, int> FaceMap;
typedef std::unordered_set<FaceSharedPtr, FaceHashPointer> FaceSetPointer;
