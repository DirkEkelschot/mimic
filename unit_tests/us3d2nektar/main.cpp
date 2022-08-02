#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
//#include "/Users/dekelsch/Software/boost_1_71_0/boost/core/ignore_unused.hpp"
#include <iomanip>
#include <sstream>


template<typename K>
void sort3(K& x, K& y, K& z)
{
#define SWAP(a,b) if (a > b) std::swap(a,b);
    SWAP(y, z);
    SWAP(x, z);
    SWAP(x, y);
#undef SWAP
}

//void ReadXmlFile(const char*  filename)
//{
//    TiXmlDocument doc( filename );
//    bool loadOkay = doc.LoadFile();
//
//    if ( !loadOkay )
//    {
//        std::cout<< "Could not load test file 'demotest.xml'. Error=. Exiting" << doc.ErrorDesc() << std::endl;
//        exit( 1 );
//    }
//    else
//    {
//        std::cout << "Read succesfully!" << std::endl;
//    }
//}



void WriteXmlFile(const char*  filename,Array<double>* xcn,Array<int>* iet,Array<int>* ien,Array<int>* if_Nv,Array<int>* ifn,std::map<int,std::vector<int> > edgeMap,std::map<int,std::vector<int> > element_map,Array<int>* zdefs,std::map<int,std::vector<int> > faceMap,Array<int>* ief)
{
    TiXmlDocument doc(filename);
    TiXmlDeclaration *decl = new TiXmlDeclaration("1.0", "utf-8", "");
    doc.LinkEndChild(decl);

    TiXmlElement *root = new TiXmlElement("NEKTAR");
    doc.LinkEndChild(root);

    TiXmlElement *geom = new TiXmlElement("GEOMETRY");
    root->LinkEndChild(geom);

    geom->SetAttribute("DIM", 3);
    geom->SetAttribute("SPACE", 3);

    // Add Vertices

    TiXmlElement *vertTag = new TiXmlElement("VERTEX");

    for(int i=0;i<xcn->getNrow();i++)
    {
        std::stringstream s;
        s << scientific << setprecision(8) << xcn->getVal(i,0) << " " << xcn->getVal(i,1) << " " << xcn->getVal(i,2);
        TiXmlElement *t = new TiXmlElement("V");
        t->SetAttribute("ID",i);
        TiXmlText *vList = new TiXmlText(s.str().c_str());
        t->LinkEndChild(vList);
        vertTag->LinkEndChild(t);
    }

    geom->LinkEndChild(vertTag);


    // Add Vertices

    TiXmlElement *edgeTag = new TiXmlElement("EDGE");
    std::map<int,std::vector<int> >::iterator itedge;
    for(itedge=edgeMap.begin();itedge!=edgeMap.end();itedge++)
    {
        std::stringstream s;

        s << itedge->second[0]<< " " << itedge->second[1];

        TiXmlElement *e = new TiXmlElement("E");
        e->SetAttribute("ID",itedge->first);
        TiXmlText *vList = new TiXmlText(s.str().c_str());
        e->LinkEndChild(vList);
        edgeTag->LinkEndChild(e);
    }

    geom->LinkEndChild(edgeTag);

    // Add Faces

    TiXmlElement *faceTag = new TiXmlElement("FACE");
    std::map<int,std::vector<int> >::iterator itface;

    for(itface=faceMap.begin();itface!=faceMap.end();itface++)
    {
        if(itface->second.size()==3)
        {
            std::stringstream s;
            s << itface->second[0] << " " << itface->second[1] << " " << itface->second[2];
            TiXmlElement *f = new TiXmlElement("T");
            f->SetAttribute("ID",itface->first);
            TiXmlText *vList = new TiXmlText(s.str().c_str());
            f->LinkEndChild(vList);
            faceTag->LinkEndChild(f);
        }

        if(itface->second.size()==4)
        {
            std::stringstream s;
            s << itface->second[0] << " " << itface->second[1] << " " << itface->second[2] << " " << itface->second[3];
            TiXmlElement *f = new TiXmlElement("Q");
            f->SetAttribute("ID",itface->first);
            TiXmlText *vList = new TiXmlText(s.str().c_str());
            f->LinkEndChild(vList);
            faceTag->LinkEndChild(f);
        }

    }

    geom->LinkEndChild(faceTag);

    // Add Elements

    TiXmlElement *elemTag = new TiXmlElement("ELEMENT");

    std::vector<int> prisms;
    std::vector<int> tets;
    std::map<int,std::vector<int> >::iterator itelem;

    for(itelem=element_map.begin();itelem!=element_map.end();itelem++)
    {
        if(itelem->second.size()==4)
        {
            std::stringstream s;
            s << itelem->second[0] << " " << itelem->second[1] << " " << itelem->second[2] << " " << itelem->second[3];
            TiXmlElement *e = new TiXmlElement("A");
            e->SetAttribute("ID",itelem->first);
            TiXmlText *vList = new TiXmlText(s.str().c_str());
            e->LinkEndChild(vList);
            elemTag->LinkEndChild(e);
            tets.push_back(itelem->first);
        }
        if(itelem->second.size()==5)
        {
            std::stringstream s;
            s << itelem->second[0] << " " << itelem->second[1] << " " << itelem->second[2] << " " << itelem->second[3] << " " << itelem->second[4];
            TiXmlElement *e = new TiXmlElement("R");
            e->SetAttribute("ID",itelem->first);
            TiXmlText *vList = new TiXmlText(s.str().c_str());
            e->LinkEndChild(vList);
            elemTag->LinkEndChild(e);
            prisms.push_back(itelem->first);
        }

    }

    std::cout << "A " << tets[0] << " " << tets[tets.size()-1] << std::endl;
    std::cout << "R " << prisms[0] << " " << prisms[prisms.size()-1] << std::endl;
    geom->LinkEndChild(elemTag);

    TiXmlElement *compTag = new TiXmlElement("COMPOSITE");
    TiXmlElement *cT = new TiXmlElement("C");

    std::stringstream ssT;
//    for(int i=0;i<tets.size();i++)
//    {
//        if(i<(tets.size()-1))
//        {
//            if(i==0)
//            {
//                ssT << "[" << tets[i]<< " ";
//            }
//            else
//            {
//                ssT << tets[i]<< " ";
//            }
//        }
//        if(i==(tets.size()-1))
//        {
//            ssT << tets[i]<< "]";
//        }
//    }

    ssT << " A[" << tets[0] << "-" << tets[tets.size()-1] << "]";

    cT->SetAttribute("ID",0);
    TiXmlText *vListT = new TiXmlText(ssT.str().c_str());
    cT->LinkEndChild(vListT);
    compTag->LinkEndChild(cT);

    TiXmlElement *cP = new TiXmlElement("C");

    std::stringstream ssP;
//    for(int i=0;i<prisms.size();i++)
//    {
//        if(i<(prisms.size()-1))
//        {
//            if(i==0)
//            {
//                ssP << "[" << prisms[i]<< " ";
//            }
//            else
//            {
//                ssP << prisms[i]<< " ";
//            }
//        }
//        if(i==(prisms.size()-1))
//        {
//            ssP << prisms[i]<< "]";
//        }
//    }

    ssP << " R[" << prisms[0] << "-" << prisms[prisms.size()-1] << "]";

    cP->SetAttribute("ID",1);
    TiXmlText *vListP = new TiXmlText(ssP.str().c_str());
    cP->LinkEndChild(vListP);
    compTag->LinkEndChild(cP);

    for(int i=3;i<zdefs->getNrow();i++)
    {
        std::stringstream ssB;
        ssB << " F[" <<  zdefs->getVal(i,3)-1 << "-" << zdefs->getVal(i,4)-1 << "]";
        TiXmlElement *cB = new TiXmlElement("C");
        cB->SetAttribute("ID",i);
        TiXmlText *vListB = new TiXmlText(ssB.str().c_str());
        cB->LinkEndChild(vListB);
        compTag->LinkEndChild(cB);
    }

    geom->LinkEndChild(compTag);






    if ( doc.Error() )
    {
        std::cout << "Error in "<< doc.Value()<< ": " <<  doc.ErrorDesc() << std::endl;
        exit( 1 );
    }
    doc.SaveFile();

}



//vector<int> PrismReordering(std::vector<int> pe)
//{
//    const int order = 1;
//    const int n     = order - 1;
//
//    int i;
//    vector<int> mapping(6);
//
//    // To get from Gmsh -> Nektar++ prism, coordinates axes are
//    // different; need to mirror in the triangular faces, and then
//    // reorder vertices to make ordering anticlockwise on base quad.
//    static int nekToGmshVerts[6] = {3, 4, 1, 0, 5, 2};
//    // Inverse (gmsh vert -> Nektar++ vert): 3->0, 4->1, 1->2, 0->3, 5->4, 2->5
//
//    for (i = 0; i < 6; ++i)
//    {
//        mapping[i] = nekToGmshVerts[i];
//    }
//
//    if (order == 1)
//    {
//        return mapping;
//    }
//
//    // Curvilinear edges.
//    mapping.resize(6 + 9 * n);
//
//    // Nektar++ edges:
//    //   0 =    1 =    2 =    3 =    4 =    5 =    6 =    7 =    8 =
//    //   {0,1}, {1,2}, {3,2}, {0,3}, {0,4}, {1,4}, {2,5}, {3,5}, {4,5}
//    // Gmsh edges (from Geo/MPrism.h of Gmsh source):
//    //   {0,1}, {0,2}, {0,3}, {1,2}, {1,4}, {2,5}, {3,4}, {3,5}, {4,5}
//    // Apply inverse of gmshToNekVerts map:
//    //   {3,2}, {3,5}, {3,0}, {2,5}, {2,1}, {5,4}, {0,1}, {0,4}, {1,4}
//    // Nektar++ mapping (negative indicates reverse orientation):
//
//    //   = 2    = 7    = -3   = 6    = -1   = -8   = 0    = 4    = 5
//    static int gmshToNekEdge[9] = {2, 7, 3, 6, 1, 8, 0, 4, 5};
//    static int gmshToNekRev[9]  = {0, 0, 1, 0, 1, 1, 0, 0, 0};
//
//    // Reorder edges.
//    int offset, cnt = 6;
//    for (i = 0; i < 9; ++i)
//    {
//        offset = 6 + n * gmshToNekEdge[i];
//
//        if (gmshToNekRev[i])
//        {
//            for (int j = 0; j < n; ++j)
//            {
//                mapping[offset + n - j - 1] = cnt++;
//            }
//        }
//        else
//        {
//            for (int j = 0; j < n; ++j)
//            {
//                mapping[offset + j] = cnt++;
//            }
//        }
//    }
//
////    if (conf.m_faceNodes == false)
////    {
////        return mapping;
////    }
//
//    int nTriInt  = n * (n - 1) / 2;
//    int nQuadInt = n * n;
//
//    // Nektar++ faces:
//    //   0 = {0,1,2,3}, 1 = {0,1,4}, 2 = {1,2,5,4}, 3 = {3,2,5}, 4 = {0,3,5,4}
//    // Gmsh faces (from Geo/MPrism.h of Gmsh source):
//    //   {0,2,1}, {3,4,5}, {0,1,4,3}, {0,3,5,2}, {1,2,5,4}
//    // Apply inverse of gmshToNekVerts map:
//    //   {3,5,2}, {0,1,4}, {3,2,1,0}, {3,0,4,5}, {2,5,4,1}
//    //   = 3      = 1      = 0        = 4        = 2
//    // This gives gmsh -> Nektar++ faces:
//    static int gmshToNekFace[5] = {3, 1, 0, 4, 2};
//
//    // Face offsets
//    vector<int> offsets(5), offsets2(5);
//    offset = 6 + 9 * n;
//
//    // Offsets in the gmsh order: face ordering is TTQQQ
//    offsets[0] = offset;
//    offsets[1] = offset + nTriInt;
//    offsets[2] = offset + 2 * nTriInt;
//    offsets[3] = offset + 2 * nTriInt + nQuadInt;
//    offsets[4] = offset + 2 * nTriInt + 2 * nQuadInt;
//
//    // Offsets in the Nektar++ order: face ordering is QTQTQ
//    offsets2[0] = offset;
//    offsets2[1] = offset + nQuadInt;
//    offsets2[2] = offset + nQuadInt + nTriInt;
//    offsets2[3] = offset + 2 * nQuadInt + nTriInt;
//    offsets2[4] = offset + 2 * nQuadInt + 2 * nTriInt;
//
//    mapping.resize(6 + 9 * n + 3 * nQuadInt + 2 * nTriInt);
//
//    offset = 6 + 9 * n;
//
//    vector<int> triVertId(3);
//    triVertId[0] = 0;
//    triVertId[1] = 1;
//    triVertId[2] = 2;
//
//    vector<int> quadVertId(4);
//    quadVertId[0] = 0;
//    quadVertId[1] = 1;
//    quadVertId[2] = 2;
//    quadVertId[3] = 3;
//
//    for (i = 0; i < 5; ++i)
//    {
//        int face    = gmshToNekFace[i];
//        int offset2 = offsets[i];
//        offset      = offsets2[face];
//
//        bool tri = i < 2;
//        int nFacePts = tri ? nTriInt : nQuadInt;
//
//        if (nFacePts == 0)
//        {
//            continue;
//        }
//
//        // Create a list of interior face nodes for this face only.
//        vector<int> faceNodes(nFacePts);
//        vector<int> toAlign(tri ? 3 : 4);
//        for (int j = 0; j < nFacePts; ++j)
//        {
//            faceNodes[j] = offset2 + j;
//        }
//
//        if (tri)
//        {
//            // Now get the reordering of this face, which puts Gmsh
//            // recursive ordering into Nektar++ row-by-row order.
//            vector<int> tmp = triTensorNodeOrdering(faceNodes, n - 1);
//            HOTriangle<int> hoTri(triVertId, tmp);
//
//            // Apply reorientation
//            if (i == 0)
//            {
//                // Triangle verts {0,2,1} --> {0,1,2}
//                toAlign[0] = 0;
//                toAlign[1] = 2;
//                toAlign[2] = 1;
//                hoTri.Align(toAlign);
//            }
//
//            // Fill in mapping.
//            for (int j = 0; j < nTriInt; ++j)
//            {
//                mapping[offset + j] = hoTri.surfVerts[j];
//            }
//        }
//        else
//        {
//            vector<int> tmp = quadTensorNodeOrdering(faceNodes, n);
//            HOQuadrilateral<int> hoQuad(quadVertId, tmp);
//
//            // Apply reorientation
//            if (i == 2)
//            {
//                toAlign[0] = 3;
//                toAlign[1] = 2;
//                toAlign[2] = 1;
//                toAlign[3] = 0;
//            }
//            else if (i == 3)
//            {
//                toAlign[0] = 1;
//                toAlign[1] = 0;
//                toAlign[2] = 3;
//                toAlign[3] = 2;
//            }
//            else if (i == 4)
//            {
//                toAlign[0] = 3;
//                toAlign[1] = 0;
//                toAlign[2] = 1;
//                toAlign[3] = 2;
//            }
//
//            hoQuad.Align(toAlign);
//
//            // Fill in mapping.
//            for (int j = 0; j < nQuadInt; ++j)
//            {
//                mapping[offset + j] = hoQuad.surfVerts[j];
//            }
//        }
//    }
//
////    if (conf.m_volumeNodes == false)
////    {
////        return mapping;
////    }
//
//    // Interior points
//    offset = offsets[4] + nQuadInt;
//    vector<int> intPoints, tmp;
//
//    for (int i = offset; i < (order+1) * (order+1) * (order+2) / 2; ++i)
//    {
//        intPoints.push_back(i);
//    }
//
//    // Reorder interior points
//    tmp = prismTensorNodeOrdering(intPoints, order - 1, log);
//    mapping.insert(mapping.end(), tmp.begin(), tmp.end());
//
//    return mapping;
//}




//vector<int> TetReordering(std::vector<int> )
//{
//    const int order = 1;
//    const int n     = order - 1;
//    const int n2    = n * (n - 1) / 2;
//
//    int i, j;
//    vector<int> mapping(4);
//
//    // Copy vertices.
//    for (i = 0; i < 4; ++i)
//    {
//        mapping[i] = i;
//    }
//
//    if (order == 1)
//    {
//        return mapping;
//    }
//
//    // Curvilinear edges.
//    mapping.resize(4 + 6 * n);
//
//    // Curvilinear edges.
//    static int gmshToNekEdge[6] = {0, 1, 2, 3, 5, 4};
//    static int gmshToNekRev[6]  = {0, 0, 1, 1, 1, 1};
//
//    // Reorder edges.
//    int offset, cnt = 4;
//    for (i = 0; i < 6; ++i)
//    {
//        offset = 4 + n * gmshToNekEdge[i];
//
//        if (gmshToNekRev[i])
//        {
//            for (int j = 0; j < n; ++j)
//            {
//                mapping[offset + n - j - 1] = cnt++;
//            }
//        }
//        else
//        {
//            for (int j = 0; j < n; ++j)
//            {
//                mapping[offset + j] = cnt++;
//            }
//        }
//    }
//
////    if (conf.m_faceNodes == false || n2 == 0)
////    {
////        return mapping;
////    }
//
//    // Curvilinear faces.
//    mapping.resize(4 + 6 * n + 4 * n2);
//
//    static int gmshToNekFace[4] = {0, 1, 3, 2};
//
//    vector<int> triVertId(3);
//    triVertId[0] = 0;
//    triVertId[1] = 1;
//    triVertId[2] = 2;
//
//    // Loop over Gmsh faces
//    for (i = 0; i < 4; ++i)
//    {
//        int face    = gmshToNekFace[i];
//        int offset2 = 4 + 6 * n + i * n2;
//        offset      = 4 + 6 * n + face * n2;
//
//        // Create a list of interior face nodes for this face only.
//        vector<int> faceNodes(n2);
//        vector<int> toAlign(3);
//        for (j = 0; j < n2; ++j)
//        {
//            faceNodes[j] = offset2 + j;
//        }
//
//        // Now get the reordering of this face, which puts Gmsh
//        // recursive ordering into Nektar++ row-by-row order.
//        vector<int> tmp = triTensorNodeOrdering(faceNodes, n - 1);
//        HOTriangle<int> hoTri(triVertId, tmp);
//
//        // Apply reorientation
//        if (i == 0 || i == 2)
//        {
//            // Triangle verts {0,2,1} --> {0,1,2}
//            toAlign[0] = 0;
//            toAlign[1] = 2;
//            toAlign[2] = 1;
//            hoTri.Align(toAlign);
//        }
//        else if (i == 3)
//        {
//            // Triangle verts {1,2,0} --> {0,1,2}
//            toAlign[0] = 1;
//            toAlign[1] = 2;
//            toAlign[2] = 0;
//            hoTri.Align(toAlign);
//        }
//
//        // Fill in mapping.
//        for (j = 0; j < n2; ++j)
//        {
//            mapping[offset + j] = hoTri.surfVerts[j];
//        }
//    }
//
////    if (conf.m_volumeNodes == false)
////    {
////        return mapping;
////    }
//
//    const int nInt = (order - 3) * (order - 2) * (order - 1) / 6;
//    if (nInt <= 0)
//    {
//        return mapping;
//    }
//
//    if (nInt == 1)
//    {
//        mapping.push_back(mapping.size());
//        return mapping;
//    }
//
//    int ntot = (order+1)*(order+2)*(order+3)/6;
//    vector<int> interior;
//
//    for (int i = 4 + 6 * n + 4 * n2; i < ntot; ++i)
//    {
//        interior.push_back(i);
//    }
//
//    if (interior.size() > 0)
//    {
//        interior = tetTensorNodeOrdering(interior, order-3);
//    }
//
//    mapping.insert(mapping.end(), interior.begin(), interior.end());
//
//    return mapping;
//}


//inline void hash_combine(std::size_t& seed)
//{
//    boost::ignore_unused(seed);
//}
//
//
//template <typename T, typename... Args>
//inline void hash_combine(std::size_t& seed, const T& v, Args... args)
//{
//    std::hash<T> hasher;
//    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//    hash_combine(seed, args...);
//}
//
//template <typename T, typename... Args>
//inline std::size_t hash_combine(const T& v, Args... args)
//{
//    boost::ignore_unused(v);
//    std::size_t seed = 0;
//    hash_combine(seed, args...);
//    return seed;
//}
//
//template<typename Iter>
//std::size_t hash_range(Iter first, Iter last)
//{
//    std::size_t seed = 0;
//    for(; first != last; ++first)
//    {
//        hash_combine(seed, *first);
//    }
//    return seed;
//}
//
//template<typename Iter>
//void hash_range(std::size_t &seed, Iter first, Iter last)
//{
//    hash_combine(seed, hash_range(first, last));
//}
//
//struct EnumHash
//{
//    template <typename T>
//    std::size_t operator()(T t) const
//    {
//        return static_cast<std::size_t>(t);
//    }
//};
//
//struct PairHash {
//    template <class T1, class T2>
//    std::size_t operator()(const std::pair<T1, T2> &p) const
//    {
//        std::size_t seed = 0;
//        auto h1 = std::hash<T1>{}(p.first);
//        auto h2 = std::hash<T2>{}(p.second);
//        hash_combine(seed, h1, h2);
//        return seed;
//    }
//};
//
//
//struct FaceHash : std::unary_function<std::vector<int>, std::size_t>
//{
//    std::size_t operator()(std::vector<int> const &vf) const
//    {
//        unsigned int nVert = vf.size();
//        std::cout << nVert << std::endl;
//        std::size_t seed = 0;
//        std::vector<unsigned int> ids(nVert);
//
//        for (int i = 0; i < nVert; ++i)
//        {
//            ids[i] = vf[i];
//        }
//
//        std::sort(ids.begin(), ids.end());
//        hash_range(seed, ids.begin(), ids.end());
//
//        return seed;
//    }
//};
//typedef std::unordered_set<std::vector<int>, FaceHash> FaceSet;
//
//
//struct FaceHashSet : std::unary_function<std::vector<int>, std::size_t>
//{
//    std::size_t operator()(std::vector<int> const &vf) const
//    {
//        unsigned int nVert = vf.size();
//        std::cout << nVert << std::endl;
//        std::size_t seed = 0;
//        std::vector<unsigned int> ids(nVert);
//
//        for (int i = 0; i < nVert; ++i)
//        {
//            ids[i] = vf[i];
//        }
//
//        std::sort(ids.begin(), ids.end());
//        hash_range(seed, ids.begin(), ids.end());
//
//        return seed;
//    }
//};
//typedef std::unordered_set<std::set<int>, FaceHashSet> FaceSetSet;

struct NekElement{
    std::vector<int> nodes;
    std::vector<std::vector<int> > m_edges;
    //std::vector<std::vector<int> > m_faces;
    std::vector<std::vector<int> > m_faceVertices;
    std::vector<std::vector<std::vector<int> > > m_faceEdges;
};



NekElement* SetTetrahedron(std::vector<int> nodes, std::vector<std::vector<double> > nodesCoords)
{
    NekElement* tet = new NekElement;
    
    int m_orientationMap[4];
    int m_origVertMap[4];
    int m_faceVertMap[4][3] = {
        {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}
    };
    
    int m_edgeVertMap[6][2] = {
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}
    };
    
    int m_faceEdgeMap[4][3] = {
        {0, 1, 2}, {0, 4, 3}, {1, 5, 4}, {2, 5, 3}
    };
    
    // Create a copy of the original vertex ordering. This is used to
    // construct a mapping, #orientationMap, which maps the original
    // face ordering to the new face ordering.
    int orig_faces[4][3];
    for (int i = 0; i < 4; ++i)
    {
        int v0id = nodes[m_faceVertMap[i][0]];
        int v1id = nodes[m_faceVertMap[i][1]];
        int v2id = nodes[m_faceVertMap[i][2]];
        sort3(v0id, v1id, v2id);
        orig_faces[i][0] = v0id;
        orig_faces[i][1] = v1id;
        orig_faces[i][2] = v2id;
    }

    // Store a copy of the original vertex ordering so we can create a
    // permutation map later.
    
    std::vector<int> origVert(nodes.size());
    
    for(int i=0;i<nodes.size();i++)
    {
        origVert[i] = nodes[i];
    }
    
//    vector<NodeSharedPtr> origVert = nodes;

    // Order vertices with highest global vertex at top degenerate
    // point. Place second highest global vertex at base degenerate
    // point.
    sort(nodes.begin(), nodes.end());

    
    
    
    // Calculate a.(b x c) if negative, reverse order of
    // non-degenerate points to correctly orientate the tet.

    double ax = nodesCoords[1][0] - nodesCoords[0][0];
    double ay = nodesCoords[1][1] - nodesCoords[0][1];
    double az = nodesCoords[1][2] - nodesCoords[0][2];
    double bx = nodesCoords[2][0] - nodesCoords[0][0];
    double by = nodesCoords[2][1] - nodesCoords[0][1];
    double bz = nodesCoords[2][2] - nodesCoords[0][2];
    double cx = nodesCoords[3][0] - nodesCoords[0][0];
    double cy = nodesCoords[3][1] - nodesCoords[0][1];
    double cz = nodesCoords[3][2] - nodesCoords[0][2];

    double nx   = (ay * bz - az * by);
    double ny   = (az * bx - ax * bz);
    double nz   = (ax * by - ay * bx);
    double nmag = sqrt(nx * nx + ny * ny + nz * nz);
    
    nx /= nmag;
    ny /= nmag;
    nz /= nmag;

    double area = 0.5 * nmag;

    // distance of top vertex from base
    double dist = cx * nx + cy * ny + cz * nz;

    if (fabs(dist) / area <= 1e-4 )
    {
        cerr << "Warning: degenerate tetrahedron, 3rd vertex is = " << dist
             << " from face" << endl;
    }
    
    if (dist < 0)
    {
        std::swap(nodes[0], nodes[1]);
    }
    
    nx   = (ay * cz - az * cy);
    ny   = (az * cx - ax * cz);
    nz   = (ax * cy - ay * cx);
    nmag = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= nmag;
    ny /= nmag;
    nz /= nmag;

    area = 0.5 * nmag;

    // distance of top vertex from base
    dist = bx * nx + by * ny + bz * nz;

    if (fabs(dist) / area <= 1e-4)
    {
        cerr << "Warning: degenerate tetrahedron, 2nd vertex is = " << dist
             << " from face" << endl;
    }

    nx   = (by * cz - bz * cy);
    ny   = (bz * cx - bx * cz);
    nz   = (bx * cy - by * cx);
    nmag = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= nmag;
    ny /= nmag;
    nz /= nmag;

    area = 0.5 * nmag;

    // distance of top vertex from base
    dist = ax * nx + ay * ny + az * nz;

    if (fabs(dist) / area <= 1e-4)
    {
        cerr << "Warning: degenerate tetrahedron, 1st vertex is = " << dist
             << " from face" << endl;
    }

    // Search for the face in the original set of face nodes. Then use
    // this to construct the #orientationMap.
    for (int i = 0; i < 4; ++i)
    {
        int v0id = nodes[m_faceVertMap[i][0]];
        int v1id = nodes[m_faceVertMap[i][1]];
        int v2id = nodes[m_faceVertMap[i][2]];
        sort3(v0id, v1id, v2id);
        for (int j = 0; j < 4; ++j)
        {
            if (v0id == orig_faces[j][0] && v1id == orig_faces[j][1] &&
                v2id == orig_faces[j][2])
            {
                m_orientationMap[j] = i;
                break;
            }
        }

        for (int j = 0; j < 4; ++j)
        {
            if (nodes[i] == origVert[j])
            {
                m_origVertMap[j] = i;
                break;
            }
        }
    }
    // Create edges (with corresponding set of edge points). Apply orientation
    // logic to get the right interior points for each edge.
    for (int i = 0; i < 6; ++i)
    {
        //std::vector<NodeSharedPtr> edgeNodes(n);

        int origEdge = -1;
        bool rev = false;
        for (int j = 0; j < 6; ++j)
        {
            if (m_edgeVertMap[i][0] == m_origVertMap[m_edgeVertMap[j][0]] &&
                m_edgeVertMap[i][1] == m_origVertMap[m_edgeVertMap[j][1]])
            {
                origEdge = j;
                break;
            }
            else if (m_edgeVertMap[i][0] == m_origVertMap[m_edgeVertMap[j][1]] &&
                     m_edgeVertMap[i][1] == m_origVertMap[m_edgeVertMap[j][0]])
            {
                origEdge = j;
                rev = true;
                break;
            }
        }
        
        if (rev)
        {
            std::vector<int> fedge(2);
            fedge[0] = nodes[m_edgeVertMap[i][1]];
            fedge[1] = nodes[m_edgeVertMap[i][0]];
            tet->m_edges.push_back(fedge);
            
        }
        else
        {
            std::vector<int> fedge(2);
            fedge[0] = nodes[m_edgeVertMap[i][0]];
            fedge[1] = nodes[m_edgeVertMap[i][1]];
            tet->m_edges.push_back(fedge);
        }
    }
    
    // Create faces
    for (int j = 0; j < 4; ++j)
    {
        std::vector<int> fnodes(3);
        std::vector<std::vector<int> >faceEdges;
        for (int k = 0; k < 3; ++k)
        {
            fnodes[k] = nodes[m_faceVertMap[j][k]];

            std::vector<int> fedge(2);
            fedge[0] = tet->m_edges[m_faceEdgeMap[j][k]][0];
            fedge[1] = tet->m_edges[m_faceEdgeMap[j][k]][1];
            faceEdges.push_back(fedge);
        }
        tet->m_faceVertices.push_back(fnodes);
        tet->m_faceEdges.push_back(faceEdges);
    }
    
    for(int i=0;i<nodes.size();i++)
    {
        tet->nodes.push_back(nodes[i]);
    }
    
    return tet;
}





NekElement* SetPrism(std::vector<int> nodes, std::vector<std::vector<double> > nodesCoords)
{
    NekElement* prism = new NekElement;
    int m_orientation = 0;
    int m_faceVertMap[5][4] = {
        {0, 1, 2, 3}, {0, 1, 4, -1}, {1, 2, 5, 4}, {3, 2, 5, -1}, {0, 3, 5, 4}
    };
    
    int m_edgeVertMap[9][2] = {{0, 1},
                               {1, 2},
                               {3, 2},
                               {0, 3},
                               {0, 4},
                               {1, 4},
                               {2, 5},
                               {3, 5},
                               {4, 5}};
    
    //int m_edge[9][2];
    
    std::vector<std::vector<int> > m_edges;
    for(int i=0;i<9;i++)
    {
        std::vector<int> edge;
        
        for(int j=0;j<2;j++)
        {
            //m_edge[i][j] = nodes[m_edgeVertMap[i][j]];
            edge.push_back(nodes[m_edgeVertMap[i][j]]);
        }
        m_edges.push_back(edge);
    }
    
    
    // linear element;
    int n = 0;

    int eid = 0;
    // Create edges (with corresponding set of edge points)
    
    int lid[6], gid[6];
    // Re-order vertices.
    for (int i = 0; i < 6; ++i)
    {
        lid[i] = i;
        gid[i] = nodes[i];
    }

    gid[0] = gid[3] = max(gid[0], gid[3]);
    gid[1] = gid[2] = max(gid[1], gid[2]);
    gid[4] = gid[5] = max(gid[4], gid[5]);
    
    for (int i = 1; i < 6; ++i)
    {
        if (gid[0] < gid[i])
        {
            swap(gid[i], gid[0]);
            swap(lid[i], lid[0]);
        }
    }
    
    if (lid[0] == 4 || lid[0] == 5)
    {
        m_orientation = 0;
    }
    
    else if (lid[0] == 1 || lid[0] == 2)
    {
        // Rotate prism clockwise in p-r plane
        vector<int> vertexmap(6);
        vertexmap[0]  = nodes[4];
        vertexmap[1]  = nodes[0];
        vertexmap[2]  = nodes[3];
        vertexmap[3]  = nodes[5];
        vertexmap[4]  = nodes[1];
        vertexmap[5]  = nodes[2];
        nodes         = vertexmap;
        m_orientation = 1;
    }
    else if (lid[0] == 0 || lid[0] == 3)
    {
        // Rotate prism counter-clockwise in p-r plane
        vector<int> vertexmap(6);
        vertexmap[0]  = nodes[1];
        vertexmap[1]  = nodes[4];
        vertexmap[2]  = nodes[5];
        vertexmap[3]  = nodes[2];
        vertexmap[4]  = nodes[0];
        vertexmap[5]  = nodes[3];
        nodes         = vertexmap;
        m_orientation = 2;
    }
    else
    {
        cerr << "Warning: possible prism orientation problem." << endl;
    }
    
    // Create faces
    int face_edges[5][4];
    int face_offset[5];
    face_offset[0] = 6 + 9 * n;
    for (int j = 0; j < 4; ++j)
    {
        int facenodes      = j % 2 == 0 ? n * n : n * (n - 1) / 2;
        face_offset[j + 1] = face_offset[j] + facenodes;
    }
    
    for (int j = 0; j < 5; ++j)
    {
        vector<int> faceVertices;
        vector<std::vector<int> > faceEdges;
//        vector<NodeSharedPtr> faceNodes;
        int nEdge = 3 - (j % 2 - 1);

        for (int k = 0; k < nEdge; ++k)
        {
            faceVertices.push_back(nodes[m_faceVertMap[j][k]]);
            
//            NodeSharedPtr a = nodes[m_faceIds[j][k]];
//            NodeSharedPtr b = nodes[m_faceIds[j][(k + 1) % nEdge]];
            
            int avid = nodes[m_faceVertMap[j][k]];
            int bvid = nodes[m_faceVertMap[j][(k + 1) % nEdge]];
            
            unsigned int i;
            for (i = 0; i < m_edges.size(); ++i)
            {
                if ((m_edges[i][0] == avid &&
                     m_edges[i][1] == bvid) ||
                    (m_edges[i][0] == bvid &&
                     m_edges[i][1] == avid))
                {
                    faceEdges.push_back(m_edges[i]);
                    face_edges[j][k] = i;
                    break;
                }
            }

            if (i == m_edges.size())
            {
                face_edges[j][k] = -1;
            }
        }
        
        prism->m_faceEdges.push_back(faceEdges);
        prism->m_faceVertices.push_back(faceVertices);
    
    }
    
    
    for(int i=0;i<nodes.size();i++)
    {
        prism->nodes.push_back(nodes[i]);
    }

    return prism;
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

    int ier,opt;
    int debug = 1;
    
    int tet_edgeVertMap[6][2] = {
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}
    };

    /// Local vertices that make up each tetrahedral face.
    int tet_faceVertMap[4][3] = {
        {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}
    };

    /// Local edges that make up each tetrahedral face.
    int tet_faceEdgeMap[4][3] = {
        {0, 1, 2}, {0, 4, 3}, {1, 5, 4}, {2, 5, 3}
    };
    
//    int prism_faceVertMap[5][4] = {
//        {0, 1, 2, 3}, {0, 1, 4, -1}, {1, 2, 5, 4}, {3, 2, 5, -1}, {0, 3, 5, 4}
//    };
    
    
//    int prism_edgeVertMap[9][2] = {{0, 1},
//                                   {1, 2},
//                                   {3, 2},
//                                   {0, 3},
//                                   {0, 4},
//                                   {1, 4},
//                                   {2, 5},
//                                   {3, 5},
//                                   {4, 5}};
//  const char* fn_grid="../test_mesh/cylinder_hybrid/grid.h5";
//  const char* fn_conn="../test_mesh/cylinder_hybrid/conn.h5";
//  const char* fn_data="../test_mesh/cylinder_hybrid/data.h5";
    
    const char* fn_grid     =    "inputs/grid.h5";
    const char* fn_conn     =    "inputs/conn.h5";
    const char* fn_metric   =    "inputs/metric.inp";
    
    std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
    
    int ReadFromStats    = metric_inputs[4];

    US3D* us3d                        = ReadUS3DGrid(fn_conn,fn_grid,ReadFromStats,comm,info);

    const char*  filename ="extrude.xml";
    
    std::set<std::set<int> > edges;
    std::map<int,std::vector<int> > edgeMap;
    std::map<std::set<int>,int > edgeMap_inv;

    int edgecnt = 0;
    std::map<int,std::vector<int> > elemMap;
    std::map<int,std::vector<int> > element_map;
    std::map<std::set<int>,int > edge_set;
    std::map<int,std::vector<int> > edge_map;
    std::map<std::set<int>, int > qface_set;
    std::map<std::set<int>, int > face_set;
    std::map<std::set<int>,int > face_set2;
    std::map<int,std::vector<int> > face_map;
    
    std::map<std::set<int>,int > pface_set;
    std::map<int,std::vector<int> > pface_map;
    std::set<int> m_face;
    std::set<int> m_edges;
    std::set<std::set<int> > edge_setset;
    int eId   = 0;
    int eIdn  = 0;
    int fId   = 0;
    int fIdn  = 0;
    int pfId  = 0;
    int pfIdn = 0;
    
    
    
    //FaceSet m_faceSet;
    std::map<std::size_t,std::vector<int> > hshCheck;
    std::set<int> hashedFaces;
    std::set<int> phashedFaces;
    
    int n = 0;

    std::map<std::set<int>, int> edgeSet_map;
    
    int edgeGid = 0;
    int nEdge   = 9;
    
//    int tet_faceVertMap[4][3] = {
//        {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}
//    };
    
    
    std::vector<std::vector<int> > tet_Face2VertMap;
    std::vector<int> face0(3);
    face0[0] = 0;
    face0[1] = 1;
    face0[2] = 2;
    tet_Face2VertMap.push_back(face0);
    std::vector<int> face1(3);
    face1[0] = 0;
    face1[1] = 1;
    face1[2] = 3;
    tet_Face2VertMap.push_back(face1);
    std::vector<int> face2(3);
    face2[0] = 1;
    face2[1] = 2;
    face2[2] = 3;
    tet_Face2VertMap.push_back(face2);
    std::vector<int> face3(3);
    face3[0] = 0;
    face3[1] = 2;
    face3[2] = 3;
    tet_Face2VertMap.push_back(face3);
    //    int prism_edgeVertMap[9][2] = {{0, 1},
    //                                   {1, 2},
    //                                   {3, 2},
    //                                   {0, 3},
    //                                   {0, 4},
    //                                   {1, 4},
    //                                   {2, 5},
    //                                   {3, 5},
    //                                   {4, 5}};
    std::vector<std::vector<int> > prism_edgeVertMap;
    std::vector<int> edge0(2);
    edge0[0] = 0;
    edge0[1] = 1;
    prism_edgeVertMap.push_back(edge0);
    std::vector<int> edge1(2);
    edge1[0] = 1;
    edge1[1] = 2;
    prism_edgeVertMap.push_back(edge1);
    std::vector<int> edge2(2);
    edge2[0] = 3;
    edge2[1] = 2;
    prism_edgeVertMap.push_back(edge2);
    std::vector<int> edge3(2);
    edge3[0] = 0;
    edge3[1] = 3;
    prism_edgeVertMap.push_back(edge3);
    std::vector<int> edge4(2);
    edge4[0] = 0;
    edge4[1] = 4;
    prism_edgeVertMap.push_back(edge4);
    std::vector<int> edge5(2);
    edge5[0] = 1;
    edge5[1] = 4;
    prism_edgeVertMap.push_back(edge5);
    std::vector<int> edge6(2);
    edge6[0] = 2;
    edge6[1] = 5;
    prism_edgeVertMap.push_back(edge6);
    std::vector<int> edge7(2);
    edge7[0] = 3;
    edge7[1] = 5;
    prism_edgeVertMap.push_back(edge7);
    std::vector<int> edge8(2);
    edge8[0] = 4;
    edge8[1] = 5;
    prism_edgeVertMap.push_back(edge8);
    
    
    
//    int prism_faceVertMap[5][4] = {
//        {0, 1, 2, 3}, {0, 1, 4, -1}, {1, 2, 5, 4}, {3, 2, 5, -1}, {0, 3, 5, 4}
//    };
    std::vector<std::vector<int> > prism_faceVertMap;
    std::vector<int> pface0(4);
    pface0[0] = 0;
    pface0[1] = 1;
    pface0[2] = 2;
    pface0[3] = 3;
    prism_faceVertMap.push_back(pface0);
    std::vector<int> pface1(3);
    pface1[0] = 0;
    pface1[1] = 1;
    pface1[2] = 4;
    prism_faceVertMap.push_back(pface1);
    std::vector<int> pface2(4);
    pface2[0] = 1;
    pface2[1] = 2;
    pface2[2] = 5;
    pface2[3] = 4;
    prism_faceVertMap.push_back(pface2);
    std::vector<int> pface3(3);
    pface3[0] = 3;
    pface3[1] = 2;
    pface3[2] = 5;
    prism_faceVertMap.push_back(pface3);
    std::vector<int> pface4(4);
    pface4[0] = 0;
    pface4[1] = 3;
    pface4[2] = 5;
    pface4[3] = 4;
    prism_faceVertMap.push_back(pface4);
    
//    std::vector<int> te(4);
//    std::vector<int> pr(6);
    int s = 0;
    // Orient the elements correctly:
    int nEl = us3d->iet->getNrow();
    std::vector<NekElement*> mesh(nEl);
    for(int i=0;i<nEl;i++)
    {
        if(us3d->iet->getVal(i,0)==2)
        {
            std::vector<int> nodes(4);
            nodes[0] = us3d->ien->getVal(i,0);
            nodes[1] = us3d->ien->getVal(i,1);
            nodes[2] = us3d->ien->getVal(i,2);
            nodes[3] = us3d->ien->getVal(i,3);
            std::vector<std::vector<double> > nodesCoords;
            for (int j = 0; j < 4; ++j)
            {
                int vid  = nodes[j];
                std::vector<double> Coords(3);
                Coords[0] = us3d->xcn->getVal(vid,0);
                Coords[1] = us3d->xcn->getVal(vid,1);
                Coords[2] = us3d->xcn->getVal(vid,2);
                
                nodesCoords.push_back(Coords);

            }
            
            NekElement* tet = SetTetrahedron(nodes,nodesCoords);
            
            mesh[i] = tet;
        }
        
        
        if(us3d->iet->getVal(i,0)==6)
        {
            std::vector<int> nodes(6);
            nodes[3] = us3d->ien->getVal(i,0);
            nodes[2] = us3d->ien->getVal(i,2);
            nodes[5] = us3d->ien->getVal(i,1);
            nodes[0] = us3d->ien->getVal(i,3);
            nodes[1] = us3d->ien->getVal(i,5);
            nodes[4] = us3d->ien->getVal(i,4);
            
            std::vector<std::vector<double> > nodesCoords;
            for (int j = 0; j < 6; ++j)
            {
                int vid  = nodes[j];
                std::vector<double> Coords(3);
                Coords[0] = us3d->xcn->getVal(vid,0);
                Coords[1] = us3d->xcn->getVal(vid,1);
                Coords[2] = us3d->xcn->getVal(vid,2);
                
                nodesCoords.push_back(Coords);

            }
            
            NekElement* prism = SetPrism(nodes,nodesCoords);
            
            mesh[i] = prism;
        }
        
        
    }
    
    /*
    
    
    
    
    std::map<int,std::map<int,int> > el2_l2g_edge;
    for(int i=0;i<us3d->iet->getNrow();i++)
    {
        if(us3d->iet->getVal(i,0)==2)
        {
            std::vector<int> te(4);
            te[0] = us3d->ien->getVal(i,0);
            te[1] = us3d->ien->getVal(i,1);
            te[2] = us3d->ien->getVal(i,2);
            te[3] = us3d->ien->getVal(i,3);
            // Redistribute is not necessary for tets.
            std::vector<int> gedges(6);
            std::map<int,int> Tet_l2g_edge;
            for(s=0;s<6;s++)
            {
                std::vector<int> edef(2);
                std::set<int> eset;
                for(int k=0;k<2;k++)
                {
                    eset.insert(te[tet_edgeVertMap[s][k]]);
                    edef[k] = te[tet_edgeVertMap[s][k]];
                }
                if(edge_set.find(eset)==edge_set.end())
                {
                    edge_set[eset]  = eId;
                    edge_map[eId]   = edef;
                    gedges[s]       = eId;
                    Tet_l2g_edge[s] = eId;
                    eId++;
                }
                else
                {
                    int eIdn        = edge_set[eset];
                    gedges[s]       = eIdn;
                    Tet_l2g_edge[s] = eIdn;

                }
            }
            
            el2_l2g_edge[i] = Tet_l2g_edge;
        }
        
        if(us3d->iet->getVal(i,0)==6)
        {
            std::vector<int> pr(6);

            pr[3] = us3d->ien->getVal(i,0);
            pr[2] = us3d->ien->getVal(i,2);
            pr[5] = us3d->ien->getVal(i,1);

            pr[0] = us3d->ien->getVal(i,3);
            pr[1] = us3d->ien->getVal(i,5);
            pr[4] = us3d->ien->getVal(i,4);
            
            std::vector<int> gpedges(9);
            std::map<std::set<int> ,int> loc_Edges;
            std::map<int ,int> loc2glob_Edges;
            std::map<int,int> Prism_l2g_edge;
            for(int s=0;s<9;s++)
            {
                std::vector<int> edef(2);
                std::set<int> edgeset;
                for(int k=0;k<2;k++)
                {
                    int vvid = pr[prism_edgeVertMap[s][k]];
                    
                    edgeset.insert(vvid);

                    edef[k] = pr[prism_edgeVertMap[s][k]];
                }
                
                if(edge_set.find(edgeset)==edge_set.end())
                {
                    edge_set[edgeset]  = eId;
                    edge_map[eId]      = edef;
                    gpedges[s]         = eId;
                    Prism_l2g_edge[s]  = eId;
                    eId++;
                }
                else
                {
                    int eIdn           = edge_set[edgeset];
                    gpedges[s]         = eIdn;
                    Prism_l2g_edge[s]  = eIdn;
                }
            }
            el2_l2g_edge[i] = Prism_l2g_edge;
        }
    }
    
    for(int i=0;i<us3d->iet->getNrow();i++)
    {
        if(us3d->iet->getVal(i,0)==2)
        {
            std::vector<int> tetr(4);
            tetr[0] = us3d->ien->getVal(i,0);
            tetr[1] = us3d->ien->getVal(i,2);
            tetr[2] = us3d->ien->getVal(i,1);
            tetr[3] = us3d->ien->getVal(i,3);
            std::vector<int> gfaces(4);

            for(int j=0;j<4;j++)
            {
                std::vector<int> fdef(3);
                std::set<int> fset_test;
                
                for(int k=0;k<3;k++)
                {
                    int loc_edge_id  = tet_faceEdgeMap[j][k];
                    int glob_edge_id = el2_l2g_edge[i][loc_edge_id];
                    fset_test.insert(tetr[tet_Face2VertMap[j][k]]);
                    //fset.insert(glob_edge_id);
                    fdef[k] = glob_edge_id;
                    
                    
                }
                if(face_set.find(fset_test)==face_set.end())
                {
                    face_set[fset_test]  = fId;
                    face_map[fId]        = fdef;
                    gfaces[j]            = fId;
                    if(fId==11 || fId==1)
                    {
                        std::cout << fId << " a tet " << i << " " << fdef[0] << " " << fdef[1] << " " << fdef[2] << std::endl;
                    }
                    fId++;
                }
                else
                {
                    int fIdn        = face_set[fset_test];
                    if(fId==11 || fId==1)
                    {
                        std::cout << fId << " b tet " << i << " " << fdef[0] << " " << fdef[1] << " " << fdef[2] << std::endl;
                    }
                    gfaces[j]       = fIdn;
                }
            }
            
            element_map[i] = gfaces;
        }
        if(us3d->iet->getVal(i,0)==6)
        {
            
            int face_edges[5][4];
            std::vector<int> pr(6);
            pr[3] = us3d->ien->getVal(i,0);
            pr[2] = us3d->ien->getVal(i,2);
            pr[5] = us3d->ien->getVal(i,1);

            pr[0] = us3d->ien->getVal(i,3);
            pr[1] = us3d->ien->getVal(i,5);
            pr[4] = us3d->ien->getVal(i,4);
            
            std::vector<int> gpfaces(5);
            std::map<int,std::vector<int> > face2edge;
            
            for(int j=0;j<5;j++)
            {
                
                //===============================================================
                // since we dont have the face2edge map from nektar++ we need to do this first.
                //===============================================================


                int nEdge = 3 - (j % 2 - 1);
                
                std::vector<int> f2edge(nEdge);

               
                for (int k = 0; k < nEdge; ++k)
                {
                    std::set<int> edgeCheck;
                    
                    edgeCheck.insert(pr[prism_faceVertMap[j][k]]);
                    edgeCheck.insert(pr[prism_faceVertMap[j][(k + 1) % nEdge]]);
                    
                    int glo_eid = edge_set[edgeCheck];
                    
                    f2edge[k]   = glo_eid;
                }
                
                face2edge[j] = f2edge;

                std::set<int> pfset;
                int sizero = prism_faceVertMap[j].size();
                for(int k=0;k<sizero;k++)
                {
                    pfset.insert(pr[prism_faceVertMap[j][k]]);
                    //pedef.push_back(face_edges[j][k]);
                    
                }
                //===============================================================
                //===============================================================
                if(face_set.find(pfset)==face_set.end())
                {
                    face_set[pfset]   = fId;
                    face_map[fId]      = f2edge;
                    gpfaces[j]         = fId;
                    if(fId==11 || fId==1)
                    {
                        std::cout << fId << " a pri " << i << " " << f2edge[0] << " " << f2edge[1] << " " << f2edge[2] << std::endl;
                    }
//
                    fId++;
                }
                if(face_set.find(pfset)!=face_set.end())
                {
                    int pfIdn        = face_set[pfset];
                    if(fId==11 || fId==1)
                    {
                        std::cout << fId << " b pri " << i << " " << f2edge[0] << " " << f2edge[1] << " " << f2edge[2] << std::endl;
                    }
                    gpfaces[j]       = pfIdn;
                }
                //pfdef.clear();
                //pfset.clear();
                
            }
            element_map[i] = gpfaces;
        }
    }
        

    
    std::cout << "number of edges " << edge_map.size() << std::endl;
    std::cout << "number of faces " << face_map.size() << " " << face_set.size() << std::endl;
//
//    ReadXmlFile(filename);
    WriteXmlFile("writefile2.xml",us3d->xcn,us3d->iet,us3d->ien,us3d->if_Nv,us3d->ifn,edge_map,element_map,us3d->zdefs,face_map,us3d->ief);
    */
    MPI_Finalize();
    
}

