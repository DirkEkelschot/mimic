#include "NekFace.h"

NekFace::NekFace(){}
NekFace::NekFace(std::vector<int> eid)
{
    int nEdges = eid.size();
    
    for(int i=0;i<eid.size();i++)
    {
        m_eId.push_back(eid[i]);
    }
    
    handled = false;
}
NekFace::~NekFace()
{
}
std::vector<int> NekFace::GetEdgeIDs() const
{
    return m_eId;
}
int NekFace::GetFaceID() const
{
    return m_id;
}
void NekFace::SetFaceID(int idje)
{
    m_id = idje;
}
void NekFace::SetFaceRef(int ref)
{
    m_ref = ref;
}
int NekFace::GetFaceRef() const
{
    return m_ref;
}
void NekFace::SetFaceLeftElement(int le)
{
    m_le = le;
}
void NekFace::SetFaceRightElement(int re)
{
    m_re = re;
}
int NekFace::GetFaceLeftElement() const
{
    return m_le;
}
int NekFace::GetFaceRightElement() const
{
    return m_re;
}
void NekFace::SetHandled(bool H)
{
    handled = H;
}
bool NekFace::GetHandled() const
{
    return handled;
}

std::size_t NekFace::ComputeFaceHash(std::vector<int> Face) const
{
    sort(Face.begin(),Face.end());
    std::size_t seed = 0;
    std::hash<int> hasher;
    for (int i : Face) {
        seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
    return seed;
}

bool NekFace::operator==(const NekFace& otherFace) const
{

    int l0 = m_eId.size();
    int l1 = otherFace.m_eId.size();
    std::vector<int> Face0(l0);
    std::vector<int> Face1(l1);
    
    for(int i=0;i<l0;i++)
    {
        Face0[i] = m_eId[i];
        Face1[i] = otherFace.m_eId[i];
    }
    if(l0!=l1)
    {
        return false;
    }
    else
    {
        for (std::vector<int>::iterator it1=Face0.begin();it1!=Face0.end();it1++)
        {
            if (find(Face1.begin(), Face1.end(),
                     *it1) == Face1.end())
            {
                return false;
            }
        }
        
        return true;
        
    }
}


bool operator==(FaceSharedPtr const &p1, FaceSharedPtr const &p2)
{
    int l0 = p1->GetEdgeIDs().size();
    int l1 = p2->GetEdgeIDs().size();
    std::vector<int> Face0(l0);
    std::vector<int> Face1(l1);
    
    for(int i=0;i<l0;i++)
    {
        Face0[i] = p1->GetEdgeIDs()[i];
        Face1[i] = p2->GetEdgeIDs()[i];
    }
    if(l0!=l1)
    {
        return false;
    }
    else
    {
        for (std::vector<int>::iterator it1=Face0.begin();it1!=Face0.end();it1++)
        {
            if (find(Face1.begin(), Face1.end(),
                     *it1) == Face1.end())
            {
                return false;
            }
        }
        return true;
    }
}

std::size_t NekFace::GetFaceHash()
{
    return m_seed;
}

