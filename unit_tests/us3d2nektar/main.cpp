#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
#include <boost/core/ignore_unused.hpp>
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


std::vector<double> ComputeCentroid(std::vector<std::vector<double> > elemcoords)
{
    std::vector<double> cent(3);
    cent[0] = 0.0;cent[1] = 0.0;cent[2] = 0.0;
    for(int i = 0;i<elemcoords.size();i++)
    {
        cent[0] = cent[0] + elemcoords[i][0];
        cent[1] = cent[1] + elemcoords[i][1];
        cent[2] = cent[2] + elemcoords[i][2];
    }
    
    cent[0] = cent[0] / elemcoords.size();
    cent[1] = cent[1] / elemcoords.size();
    cent[2] = cent[2] / elemcoords.size();
    
    return cent;
}



void WriteXmlFile(const char*  filename,Array<double>* xcn,std::map<int,std::vector<int> > edgeMap,std::map<int,std::vector<int> > element_map,Array<int>* zdefs,std::map<int,std::vector<int> > faceMap,Array<int>* ief,std::map<int,std::vector<int> > BoundaryComposites)
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

        s << itedge->second[0] << " " << itedge->second[1];

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
            if(itelem->first == 602 || itelem->first == 175238)
            {
                std::cout << "hiero = " << itelem->second[0] << " " << itelem->second[1] << " " << itelem->second[2] << " " << itelem->second[3] << " " << itelem->second[4] << std::endl;
            }
            if(itelem->second[0]==1914 && itelem->second[1]==1915 && itelem->second[2]==1916 && itelem->second[3]==1907 && itelem->second[4]==1917)
            {
                std::cout << "YESSIR " << itelem->first << std::endl;
            }
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

    ssT << " A[" << tets[0] << "-" << tets[tets.size()-1] << "]";

    cT->SetAttribute("ID",0);
    TiXmlText *vListT = new TiXmlText(ssT.str().c_str());
    cT->LinkEndChild(vListT);
    compTag->LinkEndChild(cT);

    TiXmlElement *cP = new TiXmlElement("C");

    std::stringstream ssP;


    ssP << " R[" << prisms[0] << "-" << prisms[prisms.size()-1] << "]";

    cP->SetAttribute("ID",1);
    TiXmlText *vListP = new TiXmlText(ssP.str().c_str());
    cP->LinkEndChild(vListP);
    compTag->LinkEndChild(cP);

    
    std::map<int,std::vector<int> >::iterator biter;
    
    int bcnt = 2;
    for(biter=BoundaryComposites.begin();biter!=BoundaryComposites.end();biter++)
    {
        int bid = biter->first;
        std::cout << "boundary info " << bid << " " << biter->second.size() << std::endl;
        std::stringstream ssB;
        
        ssB << " F[";
        for(int q=0;q<biter->second.size();q++)
        {
            if(q<(biter->second.size()-1))
            {
                ssB <<  biter->second[q] << ",";
            }
            if(q==(biter->second.size()-1))
            {
                ssB <<  biter->second[q];
            }
        }
        ssB << "]";
        
        TiXmlElement *cB = new TiXmlElement("C");
        cB->SetAttribute("ID",bcnt);
        TiXmlText *vListB = new TiXmlText(ssB.str().c_str());
        cB->LinkEndChild(vListB);
        compTag->LinkEndChild(cB);
        
        
        bcnt++;
    }
    
    
    
//    for(int i=3;i<zdefs->getNrow();i++)
//    {
//        std::stringstream ssB;
//        ssB << " F[" <<  zdefs->getVal(i,3)-1 << "-" << zdefs->getVal(i,4)-1 << "]";
//        TiXmlElement *cB = new TiXmlElement("C");
//        cB->SetAttribute("ID",i);
//        TiXmlText *vListB = new TiXmlText(ssB.str().c_str());
//        cB->LinkEndChild(vListB);
//        compTag->LinkEndChild(cB);
//    }

    geom->LinkEndChild(compTag);
    
    
    TiXmlElement *domTag = new TiXmlElement("DOMAIN");
    TiXmlElement *cD = new TiXmlElement("D");

    std::stringstream ssD;
//    for(int i=0;i<tets.size();i++)
//    {

    ssD << " C[0-1]";

    cD->SetAttribute("ID",0);
    TiXmlText *vListD = new TiXmlText(ssD.str().c_str());
    cD->LinkEndChild(vListD);
    domTag->LinkEndChild(cD);

    geom->LinkEndChild(domTag);



    TiXmlElement *expTag = new TiXmlElement("EXPANSIONS");
    TiXmlElement *exp0 = new TiXmlElement("E");
    exp0->SetAttribute("COMPOSITE",
                      "C[0]");
    exp0->SetAttribute("NUMMODES", 4);
    exp0->SetAttribute("TYPE", "MODIFIED");
    exp0->SetAttribute("FIELDS", "rho,rhou,rhov,rhow,E");
    
    expTag->LinkEndChild(exp0);
    TiXmlElement *exp1 = new TiXmlElement("E");
    exp1->SetAttribute("COMPOSITE",
                      "C[1]");
    exp1->SetAttribute("NUMMODES", 4);
    exp1->SetAttribute("TYPE", "MODIFIED");
    exp1->SetAttribute("FIELDS", "rho,rhou,rhov,rhow,E");
    expTag->LinkEndChild(exp1);
    
    root->LinkEndChild(expTag);
    
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


inline void hash_combine(std::size_t& seed)
{
    boost::ignore_unused(seed);
}
//
//
template <typename T, typename... Args>
inline void hash_combine(std::size_t& seed, const T& v, Args... args)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    hash_combine(seed, args...);
}

template <typename T, typename... Args>
inline std::size_t hash_combine(const T& v, Args... args)
{
    boost::ignore_unused(v);
    std::size_t seed = 0;
    hash_combine(seed, args...);
    return seed;
}
//
template<typename Iter>
std::size_t hash_range(Iter first, Iter last)
{
    std::size_t seed = 0;
    for(; first != last; ++first)
    {
        hash_combine(seed, *first);
    }
    return seed;
}
//
template<typename Iter>
void hash_range(std::size_t &seed, Iter first, Iter last)
{
    hash_combine(seed, hash_range(first, last));
}
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

std::size_t ComputeFaceHash(std::vector<int> Face)
{
    sort(Face.begin(),Face.end());
    std::size_t seed = 0;
    std::hash<int> hasher;
    for (int i : Face) {
        seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
    return seed;
}


std::size_t ComputeFaceHash_V2(std::vector<int> Face)
{
    sort(Face.begin(),Face.end());
    
    size_t hash = 0;
    hash_range(hash,Face.begin(),Face.end());
//    if(Face.size()==3)
//    {
//        hash_combine(hash, Face[0], Face[1], Face[2]);
//    }
//
//    if(Face.size()==4)
//    {
//        hash_combine(hash, Face[0], Face[1], Face[2], Face[3]);
//    }
    
    return hash;
}


struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        
        //sort(v.begin(),v.end());
        
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};
typedef std::unordered_set<std::vector<int>, VectorHash> MySet;


struct FaceHash : std::unary_function<std::vector<int>, std::size_t>
{
    std::size_t operator()(std::vector<int> const &vf) const
    {
        unsigned int nVert = vf.size();
        std::size_t seed = 0;
        std::vector<unsigned int> ids(nVert);

        for (int i = 0; i < nVert; ++i)
        {
            ids[i] = vf[i];
        }

        //std::sort(ids.begin(), ids.end());
        hash_range(seed, ids.begin(), ids.end());

        return seed;
    }
};
typedef std::unordered_set<std::vector<int>, FaceHash> FaceSet;
//

struct NekElement{
    std::map<int,std::vector<double> > nodesCoordMap;
    std::vector<std::vector<double> > nodesCoords;
    std::vector<int> nodes;
    std::vector<std::vector<int> > m_edges;
    std::vector<std::vector<int> > m_edgeMap;
    std::vector<std::vector<int> > m_faceVertices;
    std::vector<std::vector<std::vector<int> > > m_faceEdges;
    std::vector<int> gfaces;
    std::vector<int> m_eorient;
    std::map<int,int> globE2locE;
    std::map<int,int> globV2locV;
};

struct NekNode{
    int gid;
    double x;
    double y;
    double z;
};


bool CheckTetRotation(NekElement* elem, int id)
{
    bool RotationOK = true;
    double abx, aby, abz;
    
    std::vector<Vert> v(4);

    v[0].x = elem->nodesCoords[0][0];
    v[0].y = elem->nodesCoords[0][1];
    v[0].z = elem->nodesCoords[0][2];

    v[1].x = elem->nodesCoords[1][0];
    v[1].y = elem->nodesCoords[1][1];
    v[1].z = elem->nodesCoords[1][2];

    v[2].x = elem->nodesCoords[2][0];
    v[2].y = elem->nodesCoords[2][1];
    v[2].z = elem->nodesCoords[2][2];

    v[3].x = elem->nodesCoords[3][0];
    v[3].y = elem->nodesCoords[3][1];
    v[3].z = elem->nodesCoords[3][2];

    // cross product of edge 0 and 2
    abx = (v[1].y - v[0].y) * (v[2].z - v[0].z) -
          (v[1].z - v[0].z) * (v[2].y - v[0].y);
    aby = (v[1].z - v[0].z) * (v[2].x - v[0].x) -
          (v[1].x - v[0].x) * (v[2].z - v[0].z);
    abz = (v[1].x - v[0].x) * (v[2].y - v[0].y) -
          (v[1].y - v[0].y) * (v[2].x - v[0].x);

    //std::cout << "Value = " << ((v[3].x - v[0].x) * abx + (v[3].y - v[0].y) * aby +
    //                            (v[3].z - v[0].z) * abz) << std::endl;
    // inner product of cross product with respect to edge 3 should be positive
    if (((v[3].x - v[0].x) * abx + (v[3].y - v[0].y) * aby +
         (v[3].z - v[0].z) * abz) < 0.0)
    {
        cerr << "ERROR: Element " << id + 1 << " /val/ "<<((v[3].x - v[0].x) * abx + (v[3].y - v[0].y) * aby +
                                                           (v[3].z - v[0].z) * abz) << " is NOT counter-clockwise\n"
             << endl;
        RotationOK = false;
    }
    
    
    // Check face rotation
    if (elem->m_faceVertices[0][2] != elem->nodes[2])
    {
        cout << "ERROR: Face " << elem->gfaces[0]
             << " (vert "
             << elem->m_faceVertices[0][2]
             << ") is not aligned with base vertex of Tet "
             << id << " (vert "
             << elem->nodes[2] << ")" << endl;
        RotationOK = false;
    }

    for (int i = 1; i < 4; ++i)
    {
        if(elem->m_faceVertices[i][2] != elem->nodes[3])
        {
            cout << "ERROR: Face " << elem->gfaces[i]
                 << " is not aligned with top Vertex of Tet "
                 << id << endl;
           RotationOK = false;
        }
    }
    
    return RotationOK;
}





bool CheckPrismRotation(NekElement* elem, int id)
{
    bool RotationOK = true;

    if(elem->m_faceVertices[1][2]!=elem->nodes[4])
    {
        cout << "ERROR: Face " << elem->gfaces[1]
             << " (vert " << elem->m_faceVertices[1][2]
             << ") not aligned to face 1 singular vert of Prism "
             << id << " (vert "
             << elem->nodes[4] << ")" << endl;
        RotationOK = false;
    }

    // Check face rotation
    if (elem->m_faceVertices[3][2] !=elem->nodes[5])
    {
        cout << "ERROR: Face " << elem->gfaces[3]
             << " (vert " << elem->m_faceVertices[3][2]
             << ") not aligned to face 3 singular vert of Prism "
             << id << " (vert "
             << elem->nodes[5] << ")" << endl;
        RotationOK = false;
    }
    
    return RotationOK;
}


NekElement* SetTetrahedron(std::vector<int> nodes, std::vector<std::vector<double> > nodesCoords, int id)
{
    NekElement* tet = new NekElement;
    
    int m_orientationMap[4];
    for (int i = 0; i < 4; ++i)
    {
        m_orientationMap[i] = i;
    }
    
    int m_origVertMap[4];
    
    for (int i = 0; i < 4; ++i)
    {
        m_origVertMap[i] = i;
    }
    
    int m_faceVertMap[4][3] = {
        {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}
    };
    
    int m_edgeVertMap[6][2] = {
        {0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}
    };
    
    int m_faceEdgeMap[4][3] = {
        {0, 1, 2}, {0, 4, 3}, {1, 5, 4}, {2, 5, 3}
    };
    
    std::vector<int> origVert(nodes.size());
    std::vector<std::vector<double> > origCoords(nodes.size());
    std::map<int,std::vector<double> > nodeMap;
    for(int i=0;i<nodes.size();i++)
    {
        origVert[i]   = nodes[i];
        
        origCoords[i] = nodesCoords[i];
        nodeMap[nodes[i]] = nodesCoords[i];
        
        
        
    }
    
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
    
//    std::vector<int> origVert(nodes.size());
//    std::vector<std::vector<double> > origCoords(nodes.size());
//    for(int i=0;i<nodes.size();i++)
//    {
//        origVert[i]   = nodes[i];
//        origCoords[i] = nodesCoords[i];
//    }
    
//    vector<NodeSharedPtr> origVert = nodes;

    // Order vertices with highest global vertex at top degenerate
    // point. Place second highest global vertex at base degenerate
    // point.
    sort(nodes.begin(), nodes.end());
    
    // Calculate a.(b x c) if negative, reverse order of
    // non-degenerate points to correctly orientate the tet.

    double ax = nodeMap[nodes[1]][0] - nodeMap[nodes[0]][0];
    double ay = nodeMap[nodes[1]][1] - nodeMap[nodes[0]][1];
    double az = nodeMap[nodes[1]][2] - nodeMap[nodes[0]][2];
    double bx = nodeMap[nodes[2]][0] - nodeMap[nodes[0]][0];
    double by = nodeMap[nodes[2]][1] - nodeMap[nodes[0]][1];
    double bz = nodeMap[nodes[2]][2] - nodeMap[nodes[0]][2];
    double cx = nodeMap[nodes[3]][0] - nodeMap[nodes[0]][0];
    double cy = nodeMap[nodes[3]][1] - nodeMap[nodes[0]][1];
    double cz = nodeMap[nodes[3]][2] - nodeMap[nodes[0]][2];

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
    /**/
    
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
        tet->nodesCoords.push_back(nodeMap[nodes[i]]);
        
        //std::cout << m_origVertMap[i] << " ";
    }
    
    
    std::vector<int> m_eorient(6);
    for(int i=0;i<6;i++)
    {
        if(tet->m_edges[i][0] == nodes[m_edgeVertMap[i][0]])
        {
            m_eorient[i] = 1;
        }
        else if(tet->m_edges[i][0] == nodes[m_edgeVertMap[i][1]])
        {
            m_eorient[i] = -1;
        }
    }
    
    tet->m_eorient = m_eorient;
    
    //std::cout << std::endl;
    
    return tet;
}





NekElement* SetPrism(std::vector<int> nodes, std::vector<std::vector<double> > nodesCoords,int elid)
{
    
    
    
    
    NekElement* prism = new NekElement;
    int m_orientation = 0;
    std::map<int,int> glob2loc;
    
    int m_faceVertMap[5][4] = {{0, 1, 2,  3},
                               {0, 1, 4, -1},
                               {1, 2, 5,  4},
                               {3, 2, 5, -1},
                               {0, 3, 5,  4}};
    
    int m_edgeVertMap[9][2] = {{0, 1},
                               {1, 2},
                               {3, 2},
                               {0, 3},
                               {0, 4},
                               {1, 4},
                               {2, 5},
                               {3, 5},
                               {4, 5}};
    
    std::vector<int> m_eorient(9);
    //int m_edge[9][2];
    
    
    
    
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

    for(int i=0;i<6;i++)
    {
        glob2loc[nodes[i]] = i;
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
    
    
    //std::vector<std::vector<int> > m_edges;
    for(int i=0;i<9;i++)
    {
        std::vector<int> edge;
        
        for(int j=0;j<2;j++)
        {
            //m_edge[i][j] = nodes[m_edgeVertMap[i][j]];
            edge.push_back(nodes[m_edgeVertMap[i][j]]);
            
        }
        prism->m_edges.push_back(edge);
//prism->m_edgeMap2LocID[edge] = i;
//        m_edgeMap
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
            for (i = 0; i < prism->m_edges.size(); ++i)
            {
                if ((prism->m_edges[i][0] == avid &&
                     prism->m_edges[i][1] == bvid) ||
                    (prism->m_edges[i][0] == bvid &&
                     prism->m_edges[i][1] == avid))
                {
                    faceEdges.push_back(prism->m_edges[i]);
                    face_edges[j][k] = i;
                    break;
                }
            }

            if (i == prism->m_edges.size())
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
    
    if(face_edges[0][0] == -1) {std::cout << "face_edges[0][0] == -1" << std::endl;}
    if(face_edges[0][1] == -1) {std::cout << "face_edges[0][1] == -1" << std::endl;}
    if(face_edges[0][2] == -1) {std::cout << "face_edges[0][2] == -1" << std::endl;}
    if(face_edges[0][3] == -1) {std::cout << "face_edges[0][3] == -1" << std::endl;}
    if(face_edges[1][2] == -1) {std::cout << "face_edges[1][2] == -1" << std::endl;}
    if(face_edges[1][1] == -1) {std::cout << "face_edges[1][1] == -1" << std::endl;}
    if(face_edges[2][2] == -1) {std::cout << "face_edges[2][2] == -1" << std::endl;}
    if(face_edges[3][2] == -1) {std::cout << "face_edges[3][2] == -1" << std::endl;}
    if(face_edges[4][2] == -1) {std::cout << "face_edges[4][2] == -1" << std::endl;}
//    std::cout << "m_orientation " << m_orientation << std::endl;
    int tmp[9][2];
    
    tmp[0][0] = prism->m_edges[face_edges[0][0]][0];
    tmp[0][1] = prism->m_edges[face_edges[0][0]][1];

    tmp[1][0] = prism->m_edges[face_edges[0][1]][0];
    tmp[1][1] = prism->m_edges[face_edges[0][1]][1];

    tmp[2][0] = prism->m_edges[face_edges[0][2]][0];
    tmp[2][1] = prism->m_edges[face_edges[0][2]][1];

    tmp[3][0] = prism->m_edges[face_edges[0][3]][0];
    tmp[3][1] = prism->m_edges[face_edges[0][3]][1];

    tmp[4][0] = prism->m_edges[face_edges[1][2]][0];
    tmp[4][1] = prism->m_edges[face_edges[1][2]][1];

    tmp[5][0] = prism->m_edges[face_edges[1][1]][0];
    tmp[5][1] = prism->m_edges[face_edges[1][1]][1];

    tmp[6][0] = prism->m_edges[face_edges[2][1]][0];
    tmp[6][1] = prism->m_edges[face_edges[2][1]][1];

    tmp[7][0] = prism->m_edges[face_edges[3][2]][0];
    tmp[7][1] = prism->m_edges[face_edges[3][2]][1];

    tmp[8][0] = prism->m_edges[face_edges[4][2]][0];
    tmp[8][1] = prism->m_edges[face_edges[4][2]][1];
//
    prism->m_edges[0][0] = tmp[0][0];
    prism->m_edges[0][1] = tmp[0][1];

    prism->m_edges[1][0] = tmp[1][0];
    prism->m_edges[1][1] = tmp[1][1];

    prism->m_edges[2][0] = tmp[2][0];
    prism->m_edges[2][1] = tmp[2][1];

    prism->m_edges[3][0] = tmp[3][0];
    prism->m_edges[3][1] = tmp[3][1];

    prism->m_edges[4][0] = tmp[4][0];
    prism->m_edges[4][1] = tmp[4][1];

    prism->m_edges[5][0] = tmp[5][0];
    prism->m_edges[5][1] = tmp[5][1];

    prism->m_edges[6][0] = tmp[6][0];
    prism->m_edges[6][1] = tmp[6][1];

    prism->m_edges[7][0] = tmp[7][0];
    prism->m_edges[7][1] = tmp[7][1];

    prism->m_edges[8][0] = tmp[8][0];
    prism->m_edges[8][1] = tmp[8][1];
    
    
    
    for(int i=0;i<9;i++)
    {
        if(prism->m_edges[i][0] == nodes[m_edgeVertMap[i][0]])
        {
            m_eorient[i] = 1;
        }
        else if(prism->m_edges[i][0] == nodes[m_edgeVertMap[i][1]])
        {
            m_eorient[i] = -1;
        }
    }
    
    std::map<int,std::vector<double> > nodesCoordMap;
    std::vector<std::vector<double> > nodesCoords_copy;
//    std::cout << glob2loc[nodes[0]] << " " <<  glob2loc[nodes[1]]  << " " <<  glob2loc[nodes[2]]  << " " <<  glob2loc[nodes[3]]  << " " <<  glob2loc[nodes[4]]  << " " <<  glob2loc[nodes[5]] << std::endl;
    std::vector<double> node0(3);
    node0[0] = nodesCoords[glob2loc[nodes[0]]][0];
    node0[1] = nodesCoords[glob2loc[nodes[0]]][1];
    node0[2] = nodesCoords[glob2loc[nodes[0]]][2];
    nodesCoords_copy.push_back(node0);
    std::vector<double> node1(3);
    node1[0] = nodesCoords[glob2loc[nodes[1]]][0];
    node1[1] = nodesCoords[glob2loc[nodes[1]]][1];
    node1[2] = nodesCoords[glob2loc[nodes[1]]][2];
    nodesCoords_copy.push_back(node1);
    std::vector<double> node2(3);
    node2[0] = nodesCoords[glob2loc[nodes[2]]][0];
    node2[1] = nodesCoords[glob2loc[nodes[2]]][1];
    node2[2] = nodesCoords[glob2loc[nodes[2]]][2];
    nodesCoords_copy.push_back(node2);
    std::vector<double> node3(3);
    node3[0] = nodesCoords[glob2loc[nodes[3]]][0];
    node3[1] = nodesCoords[glob2loc[nodes[3]]][1];
    node3[2] = nodesCoords[glob2loc[nodes[3]]][2];
    nodesCoords_copy.push_back(node3);
    std::vector<double> node4(3);
    node4[0] = nodesCoords[glob2loc[nodes[4]]][0];
    node4[1] = nodesCoords[glob2loc[nodes[4]]][1];
    node4[2] = nodesCoords[glob2loc[nodes[4]]][2];
    nodesCoords_copy.push_back(node4);
    std::vector<double> node5(3);
    node5[0] = nodesCoords[glob2loc[nodes[5]]][0];
    node5[1] = nodesCoords[glob2loc[nodes[5]]][1];
    node5[2] = nodesCoords[glob2loc[nodes[5]]][2];
    nodesCoords_copy.push_back(node5);

    nodesCoordMap[nodes[0]] = node0;
    nodesCoordMap[nodes[1]] = node1;
    nodesCoordMap[nodes[2]] = node2;
    nodesCoordMap[nodes[3]] = node3;
    nodesCoordMap[nodes[4]] = node4;
    nodesCoordMap[nodes[5]] = node5;
    
    prism->m_eorient     = m_eorient;
    prism->nodesCoords   = nodesCoords_copy;
    prism->nodesCoordMap = nodesCoordMap;
    std::map<int,int> globV2locV;
    
    for(int i=0;i<prism->nodes.size();i++)
    {
        globV2locV[nodes[i]] = i;
    }
    prism->globV2locV = globV2locV;
    return prism;
}

struct PrismLines{
    std::map<int,std::map<int,std::vector<int> > > ElementLinesNodes;
    std::map<int,std::vector<int> > ElementLines;
    std::map<int,int> Prism2WallFace;
    std::map<int,int> Face2WallFace;
    std::map<int,std::vector<int> > Element2Faces;
    std::map<int,std::vector<int> > FaceLines;
    //std::map<int,std::vector<int> > facesPerElement;
};




PrismLines* GetPrismLines_V2(US3D* us3d, std::vector<int> wallfaces)
{
    PrismLines* plines = new PrismLines;
    // Get prism lines:
    int el_cur = -1;
    int nEl = us3d->iet->getNrow();
    int fid;
//    std::map<int,std::vector<int> > ElementLines;
//    std::map<int,int> Prism2WallFace;
//    std::map<int,std::vector<int> > FaceLines;
    std::map<int,std::map<int,std::vector<int> > > Lnodes;
    
    
    for(int wf=0;wf<wallfaces.size();wf++)
    {
        int wfid     = wallfaces[wf];
        int elid0    = us3d->ife->getVal(wfid,0);
        int elid1    = us3d->ife->getVal(wfid,1);
        
        std::vector<int> PrismLine;
        std::vector<int> FaceLine;
        std::map<int, std::vector<int> > PrismLineNodes;
        if(elid0<nEl)
        {
            el_cur = elid0;
        }
        else
        {
            el_cur = elid1;
        }
        
        
        std::vector<std::vector<double> > Coords;
        
        for(int j=0;j<6;j++)
        {
            int vid = us3d->ien->getVal(el_cur,j);
            
            std::vector<double> coord(3);
            
            coord[0] = us3d->xcn->getVal(vid,0);
            coord[1] = us3d->xcn->getVal(vid,1);
            coord[2] = us3d->xcn->getVal(vid,2);
            
            Coords.push_back(coord);
            
        }
        
        std::vector<double> cent = ComputeCentroid(Coords);
        
        int Nv = us3d->if_Nv->getVal(wfid,0);
        std::vector<std::vector<double> > Face;
        std::vector<double> FaceCentroid(3);
        FaceCentroid[0] = 0.0;
        FaceCentroid[1] = 0.0;
        FaceCentroid[2] = 0.0;
        std::vector<int> faceVids;
        std::vector<int> faceVids_copy;
        
        for(int k = 0;k<Nv;k++)
        {
            std::vector<double> Fcoord(3);
            int fvid = us3d->ien->getVal(el_cur,k);
            Fcoord[0] = us3d->xcn->getVal(fvid,0);
            Fcoord[1] = us3d->xcn->getVal(fvid,1);
            Fcoord[2] = us3d->xcn->getVal(fvid,2);
            faceVids.push_back(fvid);
            faceVids_copy.push_back(fvid);
            Face.push_back(Fcoord);
            FaceCentroid[0] = FaceCentroid[0] + Fcoord[0];
            FaceCentroid[1] = FaceCentroid[1] + Fcoord[1];
            FaceCentroid[2] = FaceCentroid[2] + Fcoord[2];
        }
        
        FaceCentroid[0] = FaceCentroid[0]/Nv;
        FaceCentroid[1] = FaceCentroid[1]/Nv;
        FaceCentroid[2] = FaceCentroid[2]/Nv;
        
        Vec3D* r0 = new Vec3D;
        r0->c0 = (FaceCentroid[0]-cent[0]);
        r0->c1 = (FaceCentroid[1]-cent[1]);
        r0->c2 = (FaceCentroid[2]-cent[2]);
        
        Vec3D* v0 = new Vec3D;
        v0->c0 = Face[1][0]-Face[0][0];
        v0->c1 = Face[1][1]-Face[0][1];
        v0->c2 = Face[1][2]-Face[0][2];
        Vec3D* v1 = new Vec3D;
        v1->c0 = Face[2][0]-Face[0][0];
        v1->c1 = Face[2][1]-Face[0][1];
        v1->c2 = Face[2][2]-Face[0][2];

        Vec3D* nbf = ComputeSurfaceNormal(v0,v1);
        
        double orient0 = DotVec3D(r0,nbf);
        
        std::vector<int> prismsVerts(6);
        
        prismsVerts[0] = us3d->ien->getVal(el_cur,0);
        prismsVerts[1] = us3d->ien->getVal(el_cur,1);
        prismsVerts[2] = us3d->ien->getVal(el_cur,2);
        prismsVerts[3] = us3d->ien->getVal(el_cur,3);
        prismsVerts[4] = us3d->ien->getVal(el_cur,4);
        prismsVerts[5] = us3d->ien->getVal(el_cur,5);
        
        
//        if(orient0<0.0)
//        {
//            NegateVec3D(nbf);
//            prismsVerts[0] = us3d->ien->getVal(el_cur,0);
//            prismsVerts[1] = us3d->ien->getVal(el_cur,2);
//            prismsVerts[2] = us3d->ien->getVal(el_cur,1);
//            prismsVerts[3] = us3d->ien->getVal(el_cur,3);
//            prismsVerts[4] = us3d->ien->getVal(el_cur,5);
//            prismsVerts[5] = us3d->ien->getVal(el_cur,4);
//
//        }
//        else
//        {
//            prismsVerts[0] = us3d->ien->getVal(el_cur,0);
//            prismsVerts[1] = us3d->ien->getVal(el_cur,1);
//            prismsVerts[2] = us3d->ien->getVal(el_cur,2);
//            prismsVerts[3] = us3d->ien->getVal(el_cur,3);
//            prismsVerts[4] = us3d->ien->getVal(el_cur,4);
//            prismsVerts[5] = us3d->ien->getVal(el_cur,5);
//        }
        
        
        int tetFound = 0;
        int tel = 0;
        int startV;
        
        while(tetFound==0)
        {
            PrismLine.push_back(el_cur);
            FaceLine.push_back(wfid);
            plines->Prism2WallFace[el_cur] = wfid;
            
            std::set<int> fnewSet0;
            fnewSet0.insert(prismsVerts[0]);
            fnewSet0.insert(prismsVerts[1]);
            fnewSet0.insert(prismsVerts[2]);
            
            std::set<int> fnewSet1;
            fnewSet1.insert(prismsVerts[3]);
            fnewSet1.insert(prismsVerts[4]);
            fnewSet1.insert(prismsVerts[5]);
            
            
            int testFaceId;
            for(int j=0;j<6;j++)
            {
                testFaceId = us3d->ief->getVal(el_cur,j);

                if(testFaceId == wfid)
                {
                    fid = j;
                    break;
                }
            }
            
            int nextLocalFaceId = 0;

            if(fid==0)
            {
                nextLocalFaceId = 1;
            }
            if(fid==1)
            {
                nextLocalFaceId = 0;
            }
            if(fid==2)
            {
                nextLocalFaceId = 3;
            }
            if(fid==3)
            {
                nextLocalFaceId = 2;
            }
            if(fid==4)
            {
                nextLocalFaceId = 5;
            }
            if(fid==5)
            {
                nextLocalFaceId = 4;
            }

            int nextFaceId = us3d->ief->getVal(el_cur,nextLocalFaceId);
            
            int Nv_next = us3d->if_Nv->getVal(nextFaceId,0);
            std::vector<std::vector<double> > Face_next;
            std::vector<double> FaceCentroid_next(3);
            FaceCentroid_next[0] = 0.0;
            FaceCentroid_next[1] = 0.0;
            FaceCentroid_next[2] = 0.0;
            std::vector<int> faceVids_next;
            std::vector<int> faceVids_next_copy;
            std::set<int> ftest;
            for(int k = 0;k<Nv_next;k++)
            {
                int fvid_next  = us3d->ifn->getVal(nextFaceId,k);
                ftest.insert(fvid_next);
            }
            
            if(ftest==fnewSet0)
            {
                startV = 0;
            }
            if(ftest==fnewSet1)
            {
                startV = 3;
            }
            
            std::set<int>::iterator its;
            
//            for(its=ftest.begin();its!=ftest.end();its++)
//            {
//                std::cout << "ftest " << *its << std::endl;
//            }
//
//            for(its=fnewSet0.begin();its!=fnewSet0.end();its++)
//            {
//                std::cout << "fnewSet0 " << *its << std::endl;
//            }
//
//            for(its=fnewSet1.begin();its!=fnewSet1.end();its++)
//            {
//                std::cout << "fnewSet1 " << *its << std::endl;
//            }
//
//            std::cout << "startV " << startV << std::endl;
            
            for(int k = 0;k<Nv_next;k++)
            {
                std::vector<double> Fcoord_next(3);
                int fvid_next  = prismsVerts[startV+k];
                Fcoord_next[0] = us3d->xcn->getVal(fvid_next,0);
                Fcoord_next[1] = us3d->xcn->getVal(fvid_next,1);
                Fcoord_next[2] = us3d->xcn->getVal(fvid_next,2);
                faceVids_next.push_back(fvid_next);
                faceVids_next_copy.push_back(fvid_next);

                Face_next.push_back(Fcoord_next);
                FaceCentroid_next[0] = FaceCentroid_next[0] + Fcoord_next[0];
                FaceCentroid_next[1] = FaceCentroid_next[1] + Fcoord_next[1];
                FaceCentroid_next[2] = FaceCentroid_next[2] + Fcoord_next[2];
            }
            
            FaceCentroid_next[0] = FaceCentroid_next[0]/Nv;
            FaceCentroid_next[1] = FaceCentroid_next[1]/Nv;
            FaceCentroid_next[2] = FaceCentroid_next[2]/Nv;
            
            Vec3D* r0_next = new Vec3D;
            r0_next->c0 = (FaceCentroid_next[0]-cent[0]);
            r0_next->c1 = (FaceCentroid_next[1]-cent[1]);
            r0_next->c2 = (FaceCentroid_next[2]-cent[2]);
            Vec3D* v0_next = new Vec3D;
            v0_next->c0 = Face_next[1][0]-Face_next[0][0];
            v0_next->c1 = Face_next[1][1]-Face_next[0][1];
            v0_next->c2 = Face_next[1][2]-Face_next[0][2];
            Vec3D* v1_next = new Vec3D;
            v1_next->c0 = Face_next[2][0]-Face_next[0][0];
            v1_next->c1 = Face_next[2][1]-Face_next[0][1];
            v1_next->c2 = Face_next[2][2]-Face_next[0][2];
            
            Vec3D* nbf_next = ComputeSurfaceNormal(v0_next,v1_next);
            
            double orient0_test = DotVec3D(r0,nbf);
            double orient0_next = DotVec3D(nbf,nbf_next);
            
            std::cout << "orientations " << el_cur << " " << orient0 << " " << orient0_test << " " << orient0_next << std::endl;
            
            
            if(orient0_next>0.0)
            {
                std::cout << "Both face normals are pointing in the same direction " << std::endl;
            }
            
            PrismLineNodes[el_cur] = prismsVerts;
            
            tel++;
            int nextElemId;
            int nextElemId_0 = us3d->ife->getVal(nextFaceId,0);
            int nextElemId_1 = us3d->ife->getVal(nextFaceId,1);
            
            if(nextElemId_0 == el_cur)
            {
                nextElemId = nextElemId_1;
            }
            else
            {
                nextElemId = nextElemId_0;
            }

            int ntype = us3d->iet->getVal(nextElemId,0);

            if(ntype==2)
            {
                tetFound = 1;
            }
            else
            {
                prismsVerts[0] = us3d->ien->getVal(nextElemId,0);
                prismsVerts[1] = us3d->ien->getVal(nextElemId,1);
                prismsVerts[2] = us3d->ien->getVal(nextElemId,2);
                
                
                prismsVerts[3] = us3d->ien->getVal(nextElemId,3);
                prismsVerts[4] = us3d->ien->getVal(nextElemId,4);
                prismsVerts[5] = us3d->ien->getVal(nextElemId,5);

            }
            
            plines->Face2WallFace[nextFaceId] = wfid;
            std::vector<int> facesPerElement(2);
            
            facesPerElement[0] = wfid;
            facesPerElement[1] = nextFaceId;
            plines->Element2Faces[el_cur] = facesPerElement;
            el_cur = nextElemId;
            wfid   = nextFaceId;
            nbf->c0 = -nbf_next->c0;
            nbf->c1 = -nbf_next->c1;
            nbf->c2 = -nbf_next->c2;
            
            ftest.clear();
            fnewSet0.clear();
            fnewSet1.clear();
            
        }
        
        Lnodes[wallfaces[wf]]               = PrismLineNodes;
        plines->ElementLines[wallfaces[wf]] = PrismLine;
        plines->FaceLines[wallfaces[wf]]    = FaceLine;
        
        //std::cout << "PrismLine = " << PrismLine.size() << std::endl;
        //std::cout << std::endl;
    }
    
    plines->ElementLinesNodes = Lnodes;
    
    return plines;
}



PrismLines* GetPrismLines(US3D* us3d, std::vector<int> wallfaces)
{
    PrismLines* plines = new PrismLines;
    // Get prism lines:
    int el_cur = -1;
    int nEl = us3d->iet->getNrow();
//    int prismTris[2][3] = {{0, 1, 4}, {3, 2, 5}};
//    std::map<int,int> vertex_map;
//    int nodeId = 0;
    int fid;
//    std::map<int,std::vector<int> > ElementLines;
//    std::map<int,int> Prism2WallFace;
//    std::map<int,std::vector<int> > FaceLines;
    for(int wf=0;wf<wallfaces.size();wf++)
    {
        int wfid     = wallfaces[wf];
        int elid0    = us3d->ife->getVal(wfid,0);
        int elid1    = us3d->ife->getVal(wfid,1);
        
        std::vector<int> PrismLine;
        std::vector<int> FaceLine;

        if(elid0<nEl)
        {
            el_cur = elid0;
        }
        else
        {
            el_cur = elid1;
        }
        
        int tetFound = 0;
        
        while(tetFound==0)
        {
            PrismLine.push_back(el_cur);
            FaceLine.push_back(wfid);
            plines->Prism2WallFace[el_cur] = wfid;
            
            
            int testFaceId;
            for(int j=0;j<6;j++)
            {
                testFaceId = us3d->ief->getVal(el_cur,j);

                if(testFaceId == wfid)
                {
                    fid = j;
                    break;
                }
            }
            
            std::vector<std::vector<double> > Coords;
            for(int j=0;j<6;j++)
            {
                int vid = us3d->ien->getVal(el_cur,j);
                
                std::vector<double> coord(3);
                
                coord[0] = us3d->xcn->getVal(vid,0);
                coord[1] = us3d->xcn->getVal(vid,1);
                coord[2] = us3d->xcn->getVal(vid,2);
                
                Coords.push_back(coord);
                
            }
            
            std::vector<double> cent = ComputeCentroid(Coords);
            
            int nextLocalFaceId = 0;

            if(fid==0)
            {
                nextLocalFaceId = 1;
            }
            if(fid==1)
            {
                nextLocalFaceId = 0;
            }
            if(fid==2)
            {
                nextLocalFaceId = 3;
            }
            if(fid==3)
            {
                nextLocalFaceId = 2;
            }
            if(fid==4)
            {
                nextLocalFaceId = 5;
            }
            if(fid==5)
            {
                nextLocalFaceId = 4;
            }

            int nextFaceId = us3d->ief->getVal(el_cur,nextLocalFaceId);
            
            int nextElemId;
            int nextElemId_0 = us3d->ife->getVal(nextFaceId,0);
            int nextElemId_1 = us3d->ife->getVal(nextFaceId,1);
            
            if(nextElemId_0 == el_cur)
            {
                nextElemId = nextElemId_1;
            }
            else
            {
                nextElemId = nextElemId_0;
            }

            int ntype = us3d->iet->getVal(nextElemId,0);

            if(ntype==2)
            {
                tetFound = 1;
            }
            
            plines->Face2WallFace[nextFaceId] = wfid;
            std::vector<int> facesPerElement(2);
            std::cout << "Face2Element " << el_cur << " " << nextElemId << std::endl;
            facesPerElement[0] = wfid;
            facesPerElement[1] = nextFaceId;
            plines->Element2Faces[el_cur] = facesPerElement;
            el_cur = nextElemId;
            wfid   = nextFaceId;
        }
        
        plines->ElementLines[wallfaces[wf]] = PrismLine;
        plines->FaceLines[wallfaces[wf]] = FaceLine;
        
    }
    
    return plines;
}


struct FaceMap{
    std::map<int,std::vector<int> > Faces;
    std::map<int,std::vector<int> > collisions;
};

void HashFaceToFaceMap(FaceMap* &fmap, std::vector<int> face)
{
    int dupl = 0;
    int addnew = 0;
    int nv = face.size();
    std::vector<int> tmp(nv);
    for(int i=0;i<nv;i++)
    {
        tmp[i] = face[i];

    }
    sort(tmp.begin(),tmp.end());

    size_t HA = ComputeFaceHash_V2(tmp);
    
    if(fmap->Faces.find(HA) == fmap->Faces.end())
    {
        fmap->Faces[HA] = face;
    }
    else
    {
        std::vector<int> triCom = fmap->Faces[HA];
        int err = 0;
        for(int i = 0;i < face.size();i++)
        {
            err = err + (face[i]-triCom[i]);
        }

        if(err==0)
        {
            dupl++;
        }
        else
        {
            if(fmap->collisions.find(HA)==fmap->collisions.end())
            {
                fmap->collisions[HA] = face;
                addnew++;
            }
        }
    }
}

std::vector<int> GetFaceFromMap(FaceMap* fmap, size_t key)
{
    std::vector<int> face;
    
    if(fmap->Faces.find(key)!=fmap->Faces.end())
    {
        
        
    }
    
    
    if(fmap->collisions.find(key)!=fmap->collisions.end())
    {
        face = fmap->collisions[key];
    }
    else
    {
        if(fmap->Faces.find(key)!=fmap->Faces.end())
        {
            face = fmap->Faces[key];
        }
        else
        {
            std::cout << "Warning: " << key << " is not in FaceMap" << std::endl;
        }
    }
    
    return face;
}

std::vector<double> NodeCurl(std::vector<double> nodeOr, std::vector<double> node)
{
    double xOr = nodeOr[0];
    double yOr = nodeOr[1];
    double zOr = nodeOr[2];
    
    double x = node[0];
    double y = node[1];
    double z = node[2];
    
    std::vector<double> curl(3);
    curl[0] = yOr*z - zOr*y;
    curl[1] = zOr*x - xOr*z;
    curl[2] = xOr*y - yOr*x;
    return curl;
    
}

std::vector<int> SortFaceNodes(std::map<int,std::vector<double> > coordmap,
                               std::vector<std::vector<int> > Face2NodeMap)
{
    std::vector<int> returnVal(coordmap.size());
        
    if(Face2NodeMap.size()==4)
    {
        std::vector<int> f0 = Face2NodeMap[0];
        //std::cout << "f0 " << f0[0] << " " << f0[1] << " " << f0[2] << std::endl;
        int indx0 = f0[0];
        int indx1 = f0[1];
        int indx2 = f0[2];
        int indx3 = -1;
        
        std::vector<double> a(3);
        a[0] = coordmap[indx1][0] - coordmap[indx0][0];
        a[1] = coordmap[indx1][1] - coordmap[indx0][1];
        a[2] = coordmap[indx1][2] - coordmap[indx0][2];
        
        
        std::vector<double> b(3);
        b[0] = coordmap[indx2][0] - coordmap[indx0][0];
        b[1] = coordmap[indx2][1] - coordmap[indx0][1];
        b[2] = coordmap[indx2][2] - coordmap[indx0][2];
        
        
        std::vector<int> f1 = Face2NodeMap[1];
        //std::cout << "f1 " << f1[0] << " " << f1[1] << " " << f1[2] << std::endl;
        for(int i=0;i<3;i++)
        {
            if ((f1[i] != indx0) && (f1[i] != indx1) && (f1[i] != indx2))
            {
                indx3 = f1[i];
                break;
            }
        }
        
        std::vector<double> c(3);
        
        c[0] = coordmap[indx3][0] - coordmap[indx0][0];
        c[1] = coordmap[indx3][1] - coordmap[indx0][1];
        c[2] = coordmap[indx3][2] - coordmap[indx0][2];

        std::vector<double> acurlb = NodeCurl(a,b);
        
        double acurlb_dotc = acurlb[0]*c[0] + acurlb[1]*c[1] + acurlb[2]*c[2];
        
        if (acurlb_dotc < 0.0)
        {
            returnVal[0] = indx0;
            returnVal[1] = indx1;
            returnVal[2] = indx2;
            returnVal[3] = indx3;
        }
        else
        {
            returnVal[0] = indx1;
            returnVal[1] = indx0;
            returnVal[2] = indx2;
            returnVal[3] = indx3;
        }
    }
    else if(Face2NodeMap.size()==5)
    {
        std::vector<std::vector<int> > triFaces;
        std::vector<std::vector<int> > quadFaces;
        for(int f=0;f<5;f++)
        {
            if(Face2NodeMap[f].size() == 3)
            {
                triFaces.push_back(Face2NodeMap[f]);
            }
            if(Face2NodeMap[f].size() == 4)
            {
                quadFaces.push_back(Face2NodeMap[f]);
            }
        }
        
        //std::cout << "Failed 1" << std::endl;

        int indx0, indx1, indx2, indx3, indx4;

        std::vector<int> tF = triFaces[0];
        std::vector<int> qF = quadFaces[0];
        
        indx0 = indx1 = indx2 = indx3 = indx4 = -1;
        int j;
        for (int i = 0; i < 4; ++i)
        {
            for (j = 0; j < 3; ++j)
            {
                if (tF[j] == qF[i])
                {
                    break; // same node break
                }
            }

            if (j == 3) // Vertex not in quad face
            {
                if (indx2 == -1)
                {
                    indx2 = qF[i];
                }
                else if (indx3 == -1)
                {
                    indx3 = qF[i];
                }
                else
                {
                    std::cout << "More than two vertices do not match "
                                 << "triangular face" << endl;
                }
            }
            else // if found match then set indx0,indx1;
            {
                if (indx0 == -1)
                {
                    indx0 = qF[i];
                }
                else
                {
                    indx1 = qF[i];
                }
            }
        }
        //std::cout << "Failed 2 " << std::endl;

        // Finally check for top vertex
        for (int i = 0; i < 3; ++i)
        {
            //std::cout << "Failed 2.25111 "<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;

            if (tF[i] != indx0 && tF[i] != indx1 &&
                tF[i] != indx2)
            {
                indx4 = tF[i];
                
                break;
            }
        }
        //std::cout << "Failed 2.25111 22222"<< " " << quadFaces.size() << " " << triFaces.size() << " " << coordmap.size() << " " << indx0 << " " << indx1 << " " << indx2 << " " << indx3 << std::endl;
        std::vector<double> a(3);
        a[0] = coordmap[indx1][0] - coordmap[indx0][0];
        a[1] = coordmap[indx1][1] - coordmap[indx0][1];
        a[2] = coordmap[indx1][2] - coordmap[indx0][2];
        
        //std::cout << "Failed 2.25 "<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;

        std::vector<double> b(3);
        a[0] = coordmap[indx4][0] - coordmap[indx0][0];
        a[1] = coordmap[indx4][1] - coordmap[indx0][1];
        a[2] = coordmap[indx4][2] - coordmap[indx0][2];
        
        std::vector<double> c(3);
        c[0] = coordmap[indx2][0] - coordmap[indx0][0];
        c[1] = coordmap[indx2][1] - coordmap[indx0][1];
        c[2] = coordmap[indx2][2] - coordmap[indx0][2];
        
        std::vector<double> acurlb = NodeCurl(a,b);
        //std::cout << "Failed 2.5 "<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;

        double acurlb_dotc = acurlb[0]*c[0] + acurlb[1]*c[1] + acurlb[2]*c[2];

        if (acurlb_dotc < 0.0)
        {
            returnVal[0] = indx0;
            returnVal[1] = indx1;
            returnVal[4] = indx4;
        }
        else
        {
            returnVal[0] = indx1;
            returnVal[1] = indx0;
            returnVal[4] = indx4;
        }
        
        //std::cout << "Failed 3"<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;
        
        std::vector<int> qF1 = quadFaces[1];
        std::vector<int> qF2 = quadFaces[2];
        int cnt = 0;
        for (int i = 0; i < 4; ++i)
        {
            if (qF1[i] == returnVal[1] ||
                qF1[i] == indx2)
            {
                cnt++;
            }
        }
        if (cnt == 2) // have two matching vertices
        {
            returnVal[2] = indx2;
            returnVal[3] = indx3;
        }
        else
        {
            cnt = 0;
            for (int i = 0; i < 4; ++i)
            {
                if (qF2[i] == returnVal[1] ||
                    qF2[i] == indx2)
                {
                    cnt++;
                }
            }

            if (cnt != 2) // neither of the other faces has two matching
                          // nodes so reverse
            {
                returnVal[2] = indx3;
                returnVal[3] = indx2;
            }
            else // have two matching vertices
            {
                returnVal[2] = indx2;
                returnVal[3] = indx3;
            }
        }
        
        //std::cout << "Failed 4"<<std::endl;

        // finally need to find last vertex from second triangular face.
        std::vector<int> tF1 = triFaces[1];
        for (int i = 0; i < 3; ++i)
        {
            if (tF1[i] != indx2 && tF1[i] != indx3)
            {
                returnVal[5] = tF1[i];
                break;
            }
        }
    }
    
    
    
    return returnVal;
}


std::map<int,std::vector<int> > ResetNodes(US3D* us3d,
                                           std::map<int,std::vector<int> > Element2Nodes_reset,
                                           std::map<int,int> new2old_v,
                                           std::map<int,int> TetTransferMap,
                                           std::vector<int> wallfaces,
                                           std::map<int,int> &vertex_map,
                                           std::map<int,int> &vertex_map_inv,
                                           std::map<int,std::vector<std::vector<int> > > &ElementFace2NodeMap)
{
    
    std::map<int,std::vector<int> > ElementNodes;
    
    //std::map<int,std::vector<std::vector<int> > > ElementFace2NodeMap;
    
    PrismLines* plines = GetPrismLines(us3d,wallfaces);
    int i,j,k;
    std::map<int,std::vector<int> >::iterator itm;
    std::set<int> vert_set;

    int prismTris[2][3] = {{0, 1, 4}, {3, 2, 5}};
    std::set<std::set<int> > facesDone;
    int nodeId = 0;

    //std::map<int,std::vector<int> > newPrismsOrigNodes;
    std::map<int,std::vector<int> > newPrismsNodes;
    std::map<int,std::map<int,std::vector<double> > > newPrismsCoords;
    //std::map<int,std::map<int,int> > newPrismsOldVid2NewVid;

    int prismtel0=0;
    std::map<int,int> oVid2nVid_global;


    std::cout << "Starting off resetting the vertices..." << std::endl;
    
    for(itm=plines->ElementLines.begin();itm!=plines->ElementLines.end();itm++)
    {
        int nprismOnLine = itm->second.size();
        for(int q=0;q<nprismOnLine;q++)
        {
            int elId = itm->second[q];
            std::vector<int> nodes(6);
            std::vector<int> orig_nodes(6);
            //std::vector<int> lnodes = lineNodes[elId];

            
            orig_nodes[0] = us3d->ien->getVal(elId,3);
            orig_nodes[1] = us3d->ien->getVal(elId,4);
            orig_nodes[2] = us3d->ien->getVal(elId,1);
            orig_nodes[3] = us3d->ien->getVal(elId,0);
            orig_nodes[4] = us3d->ien->getVal(elId,5);
            orig_nodes[5] = us3d->ien->getVal(elId,2);
            
            
            orig_nodes[0] = Element2Nodes_reset[elId][3];
            orig_nodes[1] = Element2Nodes_reset[elId][4];
            orig_nodes[2] = Element2Nodes_reset[elId][1];
            orig_nodes[3] = Element2Nodes_reset[elId][0];
            orig_nodes[4] = Element2Nodes_reset[elId][5];
            orig_nodes[5] = Element2Nodes_reset[elId][2];
            
//            orig_nodes[0] = us3d->ien->getVal(elId,3);
//            orig_nodes[1] = us3d->ien->getVal(elId,4);
//            orig_nodes[2] = us3d->ien->getVal(elId,1);
//            orig_nodes[3] = us3d->ien->getVal(elId,0);
//            orig_nodes[4] = us3d->ien->getVal(elId,5);
//            orig_nodes[5] = us3d->ien->getVal(elId,2);
            
//            orig_nodes[0] = us3d->ien->getVal(elId,0);
//            orig_nodes[1] = us3d->ien->getVal(elId,1);
//            orig_nodes[2] = us3d->ien->getVal(elId,4);
//            orig_nodes[3] = us3d->ien->getVal(elId,3);
//            orig_nodes[4] = us3d->ien->getVal(elId,2);
//            orig_nodes[5] = us3d->ien->getVal(elId,5);
            
            std::map<int,std::vector<double> > coordmap;
            std::map<int,std::vector<double> > coordmap_orig;
            std::vector<std::vector<double> > Coords_orig;
            
            for(int v=0;v<6;v++)
            {
                int vertexid = new2old_v[orig_nodes[v]];
                
                std::vector<double> Coords(3);
                Coords[0] = us3d->xcn->getVal(vertexid,0);
                Coords[1] = us3d->xcn->getVal(vertexid,1);
                Coords[2] = us3d->xcn->getVal(vertexid,2);
                Coords_orig.push_back(Coords);
                coordmap_orig[vertexid] = Coords;
            }
            if(elId == 4 || elId == 6 || elId == 28 || elId == 30)
            {
                std::cout << "elID " << elId << std::endl;
                for(int v=0;v<6;v++)
                {
                    std::cout << v << " " << orig_nodes[v] << " " << Coords_orig[v][0] << " " << Coords_orig[v][1] << " " << Coords_orig[v][2] << std::endl;
                    
                }
                std::cout << std::endl;
            }
            
            
//            std::vector<std::vector<int> > ElementFace2Node;
//            for(int f=0;f<5;f++)
//            {
//                int fid = us3d->ief->getVal(elId,f);
//                int Nv = us3d->if_Nv->getVal(fid,0);
//
//                std::vector<int> fv(Nv);
//
//                for(int ff=0;ff<Nv;ff++)
//                {
//                    fv[ff] = us3d->ifn->getVal(elId,ff);
//                }
//
//                ElementFace2Node.push_back(fv);
//            }
            
            NekElement* prismOrig = SetPrism(orig_nodes,Coords_orig,elId);
//
            std::vector<std::vector<int> > ElementFace2Node = prismOrig->m_faceVertices;
//            std::map<int,std::vector<double> > coordMap     = prismOrig->nodesCoordMap;

//            std::vector<int> NodesReOrder = SortFaceNodes(coordMap, ElementFace2Node);
            
            
            if(elId == 4 || elId == 6 || elId == 28 || elId == 30)
            {
                std::cout << "elID after 1" << elId << std::endl;
                for(int v=0;v<6;v++)
                {
                    std::cout << v << " " << prismOrig->nodes[v] << " " << Coords_orig[v][0] << " " << Coords_orig[v][1] << " " << Coords_orig[v][2] << std::endl;
                    
                }
                std::cout << std::endl;
            }
            
            
            std::vector<int> NodesReOrder = prismOrig->nodes;
//
//            std::map<int,std::vector<double> > coordmap_origReorder;
//            std::vector<std::vector<double> > Coords_origReorder;
//
//
//
//            for(int v=0;v<6;v++)
//            {
//                int vertexid = NodesReOrder[v];
//
//                std::vector<double> Coords_Reorder(3);
//                Coords_Reorder[0] = us3d->xcn->getVal(vertexid,0);
//                Coords_Reorder[1] = us3d->xcn->getVal(vertexid,1);
//                Coords_Reorder[2] = us3d->xcn->getVal(vertexid,2);
//                Coords_origReorder.push_back(Coords_Reorder);
//                coordmap_origReorder[vertexid] = Coords_Reorder;
//
//            }
//
//            NekElement* prismOrigReorder = SetPrism(NodesReOrder,Coords_origReorder,elId);
//
//
//          NekElement* prismOrigReorder = SetPrism(orig_nodes,Coords_orig,elId);

            
            std::set<int> face0;
            std::vector<int> face0v;
            face0.insert(NodesReOrder[0]);
            face0.insert(NodesReOrder[4]);
            face0.insert(NodesReOrder[1]);
            face0v.push_back(NodesReOrder[0]);
            face0v.push_back(NodesReOrder[4]);
            face0v.push_back(NodesReOrder[1]);
            
            std::set<int> face1;
            std::vector<int> face1v;
            face1.insert(NodesReOrder[3]);
            face1.insert(NodesReOrder[5]);
            face1.insert(NodesReOrder[2]);
            face1v.push_back(NodesReOrder[3]);
            face1v.push_back(NodesReOrder[5]);
            face1v.push_back(NodesReOrder[2]);
            
            //std::vector<double> cent = ComputeCentroid(prismOrigReorder->nodesCoords);
            
            std::vector<int> facespelem = plines->Element2Faces[elId];
            
            std::map<int,int> oVid2nVid;
            
            if(facesDone.find(face0)!=facesDone.end() &&
               facesDone.find(face1)!=facesDone.end())
            {
                std::cout << "ERROR: Both faces are process already..." << std::endl;
            }
            
            if(facesDone.find(face0)==facesDone.end() &&
               facesDone.find(face1)==facesDone.end())
            {
                for (j = 0; j < 2; ++j)
                {
                    for (k = 0; k < 3; ++k)
                    {
                        int vid = NodesReOrder[prismTris[j][k]];
                        
                        //std::vector<double> coords = coordmap_orig[vid];
                        
                        if(vertex_map.find(vid)==vertex_map.end())
                        {
                            vertex_map[vid]          = nodeId;
                            vertex_map_inv[nodeId]   = vid;
                            nodes[prismTris[j][k]]   = nodeId;
                            nodeId++;
                        }
                        else
                        {
                            int nodeIdn              = vertex_map[vid];
                            nodes[prismTris[j][k]]   = nodeIdn;
                        }
                    }
                }

                facesDone.insert(face0);
                facesDone.insert(face1);
            }
            
            if(facesDone.find(face0)!=facesDone.end() &&
               facesDone.find(face1)==facesDone.end())
            {
                int tmp1[3] = {NodesReOrder[prismTris[0][0]],
                    NodesReOrder[prismTris[0][1]],
                    NodesReOrder[prismTris[0][2]]};
                                
                int tmp2[3] = {0, 1, 2};
                
                if (tmp1[0] > tmp1[1])
                {
                    swap(tmp1[0], tmp1[1]);
                    swap(tmp2[0], tmp2[1]);
                }

                if (tmp1[1] > tmp1[2])
                {
                    swap(tmp1[1], tmp1[2]);
                    swap(tmp2[1], tmp2[2]);
                }

                if (tmp1[0] > tmp1[2])
                {
                    swap(tmp1[0], tmp1[2]);
                    swap(tmp2[0], tmp2[2]);
                }
                
                for(j = 0; j < 3; ++j)
                {
                    int vid = NodesReOrder[prismTris[0][tmp2[j]]];
                    nodes[prismTris[0][tmp2[j]]] = vertex_map[vid];
                }
                
                // Renumber this face so that highest ID matches.
                for (j = 0; j < 3; ++j)
                {
                    int vid = NodesReOrder[prismTris[1][tmp2[j]]];
                    std::vector<double> coords = coordmap_orig[vid];
                    
                    if(vertex_map.find(vid)==vertex_map.end())
                    {
                        vertex_map[vid]                 = nodeId;
                        vertex_map_inv[nodeId]          = vid;
                        nodes[prismTris[1][tmp2[j]]]    = nodeId;
                        nodeId++;
                    }
                    else
                    {
                        int nodeIdn                     = vertex_map[vid];
                        nodes[prismTris[1][tmp2[j]]]    = nodeIdn;
                    }
                }

                facesDone.insert(face1);
                
            }
            if(facesDone.find(face0)==facesDone.end() &&
               facesDone.find(face1)!=facesDone.end())
            {
                int tmp1[3] = {NodesReOrder[prismTris[1][0]],
                    NodesReOrder[prismTris[1][1]],
                    NodesReOrder[prismTris[1][2]]};
            
                
                int tmp2[3] = {0, 1, 2};
                
                if (tmp1[0] > tmp1[1])
                {
                    swap(tmp1[0], tmp1[1]);
                    swap(tmp2[0], tmp2[1]);
                }

                if (tmp1[1] > tmp1[2])
                {
                    swap(tmp1[1], tmp1[2]);
                    swap(tmp2[1], tmp2[2]);
                }

                if (tmp1[0] > tmp1[2])
                {
                    swap(tmp1[0], tmp1[2]);
                    swap(tmp2[0], tmp2[2]);
                }
                                
                for(j = 0; j < 3; ++j)
                {
                    int vid = NodesReOrder[prismTris[1][tmp2[j]]];
                    nodes[prismTris[1][tmp2[j]]] = vertex_map[vid];
                }
                
                // Renumber this face so that highest ID matches.
                // Renumber this face so that highest ID matches.
                
                for (j = 0; j < 3; ++j)
                {
                    int vid = NodesReOrder[prismTris[0][tmp2[j]]];
                    std::vector<double> coords = coordmap_orig[vid];
                    
                    if(vertex_map.find(vid)==vertex_map.end())
                    {
                        vertex_map[vid]              = nodeId;
                        vertex_map_inv[nodeId]       = vid;
                        nodes[prismTris[0][tmp2[j]]] = nodeId;
                        nodeId++;
                    }
                    else
                    {
                        int nodeIdn                  = vertex_map[vid];
                        nodes[prismTris[0][tmp2[j]]] = nodeIdn;
                    }
                }

                facesDone.insert(face0);
            }
            
            
            if(elId == 4 || elId == 6 || elId == 28 || elId == 30)
            {
                std::cout << "elID after 2 " << elId << std::endl;
                for(int v=0;v<6;v++)
                {
                    std::cout << vertex_map_inv[nodes[v]] << std::endl;
                    
                }
                std::cout << std::endl;
            }
            
            
            ElementNodes[elId]               = nodes;
            ElementFace2NodeMap[elId]        = ElementFace2Node;
            prismtel0++;
             
        }
    }
    
    // fill in any unset nodes at from other shapes
    std::map<int,int>::iterator itelm2;
    for(itelm2=TetTransferMap.begin();itelm2!=TetTransferMap.end();itelm2++)
    {
        int elId = itelm2->first;
        std::vector<int> tnodes(4);
        std::vector<int> nodes(4);
        tnodes[0] = us3d->ien->getVal(itelm2->first,0);
        tnodes[1] = us3d->ien->getVal(itelm2->first,1);
        tnodes[2] = us3d->ien->getVal(itelm2->first,2);
        tnodes[3] = us3d->ien->getVal(itelm2->first,3);
        
        tnodes[0] = Element2Nodes_reset[elId][0];
        tnodes[1] = Element2Nodes_reset[elId][1];
        tnodes[2] = Element2Nodes_reset[elId][2];
        tnodes[3] = Element2Nodes_reset[elId][3];
        
        std::cout << Element2Nodes_reset[elId][0] << " " << Element2Nodes_reset[elId][1] << " " << Element2Nodes_reset[elId][2] << " " << Element2Nodes_reset[elId][3] << std::endl;
        
        int oldid = itelm2->first;
        
        for(int t=0;t<4;t++)
        {
            int vid = tnodes[t];
            if(vertex_map.find(vid)==vertex_map.end())
            {
                vertex_map[vid]         = nodeId;
                vertex_map_inv[nodeId]  = vid;
                nodes[t]                = nodeId;
                nodeId++;
            }
            else
            {
                int nodeIdn = vertex_map[vid];
                nodes[t] = nodeIdn;
            }
        }
        
        //ElementFace2NodeMap[oldid] = ElementFace2Node;
        ElementNodes[oldid]        = nodes;
    }
    /* */
    return ElementNodes;
}


std::map<int,std::vector<int> > ActualResetNodes(US3D* us3d, std::map<int,int> PrismTransferMap,std::map<int,int> TetTransferMap,std::map<int,int> &new2old_v)
{
    std::map<int,int>::iterator itr;
    std::map<int,std::vector<int> > ElementNodeReset;
    std::map<int,int> vreset;
    std::map<int,int> vreset_inv;
    int vidn = 0;
    for(itr=PrismTransferMap.begin();itr!=PrismTransferMap.end();itr++)
    {
        int elId = itr->first;
        std::vector<int> nodes(6);
        
        for(int v=0;v<6;v++)
        {
            int vid = us3d->ien->getVal(elId,v);
            
            if(vreset.find(vid)==vreset.end())
            {
                vreset[vid]      = vidn;
                new2old_v[vidn]  = vid;
                nodes[v]         = vidn;
                vidn++;
            }
            else
            {
                int vidn2 = vreset[vid];
                nodes[v]  = vidn2;
            }
        }
        
        ElementNodeReset[elId]=nodes;
        
    }
    
    
    std::map<int,int>::iterator itelm2;
    for(itelm2=TetTransferMap.begin();itelm2!=TetTransferMap.end();itelm2++)
    {
        std::vector<int> tnodes(4);
        std::vector<int> nodes(4);
        tnodes[0] = us3d->ien->getVal(itelm2->first,0);
        tnodes[1] = us3d->ien->getVal(itelm2->first,1);
        tnodes[2] = us3d->ien->getVal(itelm2->first,2);
        tnodes[3] = us3d->ien->getVal(itelm2->first,3);
        
        int oldid = itelm2->first;
        
        for(int t=0;t<4;t++)
        {
            int vid = tnodes[t];
            if(vreset.find(vid)==vreset.end())
            {
                vreset[vid]         = vidn;
                new2old_v[vidn]     = vid;
                nodes[t]            = vidn;
                vidn++;
            }
            else
            {
                int nodeIdn = vreset[vid];
                nodes[t] = nodeIdn;
            }
        }
        
        //ElementFace2NodeMap[oldid] = ElementFace2Node;
        ElementNodeReset[oldid]      = nodes;
    }
    
    
    return ElementNodeReset;
    
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
    
    int nEl = us3d->iet->getNrow();
    
    std::map<int,int> TetTransferMap;
    std::map<int,int> TetTransferMap_Inv;
    int tetcnt = 0;
    for(int el=0;el<us3d->iet->getNrow();el++)
    {
        if(us3d->iet->getVal(el,0)==2)
        {
            TetTransferMap[el] = tetcnt;
            TetTransferMap_Inv[tetcnt] = el;
            tetcnt++;
        }
        
    }
    std::map<int,int> PrismTransferMap;
    std::map<int,int> PrismTransferMap_Inv;
    int prismcnt = tetcnt;
    for(int el=0;el<us3d->iet->getNrow();el++)
    {
        if(us3d->iet->getVal(el,0)==6)
        {
            PrismTransferMap[el] = prismcnt;
            PrismTransferMap_Inv[prismcnt] = el;
            prismcnt++;
        }
        
    }
    
    std::cout << "PrismTransferMap  " << PrismTransferMap.size() << " " << TetTransferMap.size() << std::endl;
    std::cout << "TetTransferMap  " << TetTransferMap.size() << " " << TetTransferMap.size() << std::endl;

    FaceSet m_ShellFaceSet;

    const char*  filename ="extrude.xml";
    
    std::set<std::set<int> > edges;
    std::map<int,std::vector<int> > edgeMap;
    std::map<std::set<int>,int > edgeMap_inv;

    int edgecnt = 0;
    std::map<int,std::vector<int> > element_map;
    std::map<std::set<int>,int > edge_set;
    std::map<int,std::vector<int> > edge_map;
    std::map<std::set<int>, int > face_set;
    std::map<int,std::vector<int> > face_map;
    std::map<std::set<int>, int> refFaceMap;
    std::vector<int> wallfaces;
    int facetel = 0;
    for(int i=0;i<us3d->ifn->getNrow();i++)
    {
        std::set<int> fset;
        
        int fref = us3d->if_ref->getVal(i,0);
        int lId  = us3d->ife->getVal(i,0);
        int rId  = us3d->ife->getVal(i,1);
        if(fref != 2)
        {
            int fNv  = us3d->if_Nv->getVal(i,0);
            
            for(int j=0;j<fNv;j++)
            {
                fset.insert(us3d->ifn->getVal(i,j));
            }
            
            if(refFaceMap.find(fset)==refFaceMap.end())
            {
                refFaceMap[fset] = fref;
            }
        }
        
        if(lId < nEl && rId < nEl)
        {
            int lnType = us3d->iet->getVal(lId,0);
            int rnType = us3d->iet->getVal(rId,0);
            
            if(lnType != rnType)
            {
                int fNv  = us3d->if_Nv->getVal(i,0);
                std::vector<int> faceVec(fNv);
                for(int j=0;j<fNv;j++)
                {
                    faceVec[j] = us3d->ifn->getVal(i,j);
                }
                pair<FaceSet::iterator, bool> testIns;
                testIns = m_ShellFaceSet.insert(faceVec);
                
                if(testIns.second)
                {
                    facetel++;
                }
            }
        }
//
//        pair<FaceSet::iterator, bool> testIns;
//        testIns = m_mesh->m_faceSet.insert(elmt[i]->GetFace(j));
//
        if(fref == 3)
        {
            wallfaces.push_back(i);
        }
    }
    
    
    std::map<int,int> vert_map;
    std::map<int,int> vert_map_inv;
    std::map<int,std::vector<std::vector<int> > > ElementFace2NodeMap;
    std::map<int,int> new2old_v;
    std::map<int,std::vector<int> > ResettedElementNodes = ActualResetNodes(us3d,PrismTransferMap,TetTransferMap,new2old_v);
    
    
    std::map<int,std::vector<int> > ElementNodes = ResetNodes(us3d,ResettedElementNodes,new2old_v,TetTransferMap,wallfaces,vert_map,vert_map_inv,ElementFace2NodeMap);
    
    
    
    
    std::map<int,NekElement*> mesh;
    std::map<int,std::vector<int> >::iterator itelem;
    for(itelem=ElementNodes.begin();itelem!=ElementNodes.end();itelem++)
    {
        if(itelem->second.size()==4)
        {
            int newTetID = TetTransferMap[itelem->first];
            std::vector<int> Nodes = itelem->second;
            std::vector<std::vector<double> > coords;
            std::map<int,std::vector<double> > coordMap;
            
            for(int i=0;i<Nodes.size();i++)
            {
                int ovid = new2old_v[vert_map_inv[Nodes[i]]];
                std::vector<double> crd(3);
                crd[0] = us3d->xcn->getVal(ovid,0);
                crd[1] = us3d->xcn->getVal(ovid,1);
                crd[2] = us3d->xcn->getVal(ovid,2);
                
                coords.push_back(crd);
                
                coordMap[Nodes[i]] = crd;
            }
            
            NekElement* tet = SetTetrahedron(Nodes,coords,newTetID);
            
            mesh[newTetID] = tet;
        }
        if(itelem->second.size()==6)
        {
            int newPrismID = PrismTransferMap[itelem->first];
            
            std::vector<int> Nodes = itelem->second;
            std::vector<std::vector<double> > coords;
            std::map<int,std::vector<double> > coordMap;
            for(int i=0;i<Nodes.size();i++)
            {
                int ovid = new2old_v[vert_map_inv[Nodes[i]]];
                std::vector<double> crd(3);
                crd[0] = us3d->xcn->getVal(ovid,0);
                crd[1] = us3d->xcn->getVal(ovid,1);
                crd[2] = us3d->xcn->getVal(ovid,2);
                coordMap[Nodes[i]] = crd;
                coords.push_back(crd);
            }
            
            if(itelem->first == 4 || itelem->first == 6 || itelem->first == 28 || itelem->first == 30)
            {
                std::cout << "el " << itelem->first << " :: " << std::endl;
                for(int i=0;i<Nodes.size();i++)
                {
                    std::cout << i << " " << Nodes[i] << " " << coords[i][0] << " " << coords[i][1] << " " << coords[i][2] << std::endl;
                }
            }
            
            /*
            std::cout << "ElementFace2NodeMap[itelem->first].size() " << ElementFace2NodeMap[itelem->first].size() << std::endl;
            std::vector<std::vector<int> > UpdateElementFace2Node;
            
            int nfaces = ElementFace2NodeMap[itelem->first].size();
            for(int q=0;q<nfaces;q++)
            {
                std::vector<int> f2n(ElementFace2NodeMap[itelem->first][q].size());
                for(int p=0;p<ElementFace2NodeMap[itelem->first][q].size();p++)
                {
                    int nvid = vert_map[ElementFace2NodeMap[itelem->first][q][p]];
                    f2n[p] = nvid;
                }
                
                UpdateElementFace2Node.push_back(f2n);
            }
            
            std::vector<int> NodesReOrder = SortFaceNodes(coordMap, UpdateElementFace2Node);
//
            for(int i=0;i<NodesReOrder.size();i++)
            {
                std::cout << "compare nodes " << NodesReOrder[i] << " " << Nodes[i] << std::endl;

            }
            */
            
            NekElement* prism = SetPrism(Nodes,coords,newPrismID);
            
            mesh[newPrismID] = prism;
            
        }
    }
    
    int eId   = 0;
    int eIdn  = 0;
    int fId   = 0;
    int fIdn  = 0;
    int pfId  = 0;
    int pfIdn = 0;
    
    std::map<int,NekElement*>::iterator itmesh;
    std::map<int,std::map<int,int> > el2_l2g_edge;
    
    
    // ===================== Setting up the edges =====================
    
    for(itmesh=mesh.begin();itmesh!=mesh.end();itmesh++)
    {
        NekElement* elem = itmesh->second;
        int nid = itmesh->first;
        
        if(elem->nodes.size()==4)
        {
            std::vector<std::vector<int> > edges = elem->m_edges;
            std::map<int,int> Tet_l2g_edge;
            
            for(int s=0;s<6;s++)
            {
                std::set<int> eset;
                for(int k=0;k<2;k++)
                {
                    eset.insert(edges[s][k]);
                }
                if(edge_set.find(eset)==edge_set.end())
                {
                    edge_set[eset]  = eId;
                    edge_map[eId]   = edges[s];
                    Tet_l2g_edge[s] = eId;
                    eId++;
                }
                else
                {
                    int eIdn        = edge_set[eset];
                    Tet_l2g_edge[s] = eIdn;
                }
            }
            el2_l2g_edge[nid] = Tet_l2g_edge;
            
        }
        if(elem->nodes.size()==6)
        {
            std::vector<std::vector<int> > edges = elem->m_edges;
            std::map<int,int> Prism_l2g_edge;
            
            for(int s=0;s<9;s++)
            {
                std::set<int> edgeset;
                for(int k=0;k<2;k++)
                {
                    edgeset.insert(edges[s][k]);
                }
                
                if(edge_set.find(edgeset)==edge_set.end())
                {
                    edge_set[edgeset]  = eId;
                    edge_map[eId]      = edges[s];
                    Prism_l2g_edge[s]  = eId;
                    eId++;
                }
                else
                {
                    int eIdn           = edge_set[edgeset];
                    Prism_l2g_edge[s]  = eIdn;
                }
            }
            
            el2_l2g_edge[nid] = Prism_l2g_edge;
        }
    }
    // =================== Done setting up the edges ==================
    
    
    
    // ===================== Setting up the faces =====================
    
    
    int gfid = 0;
    std::map<int,std::vector<int> > BoundaryComposites;
    int BadPrisms = 0;
    int BadTets = 0;
    
    for(itmesh=mesh.begin();itmesh!=mesh.end();itmesh++)
    {
        int elid = itmesh->first;
        
        std::vector<std::vector<std::vector<int> > > faceEdges = itmesh->second->m_faceEdges;
        std::vector<std::vector<int> > faceVerts = itmesh->second->m_faceVertices;

        int nfaces = faceEdges.size();
        std::vector<int> gfaces_copy(nfaces);
        
        if(elid == 28 || elid == 30)
        {
            for(int q=0;q<nfaces;q++)
            {
                int nEdges = faceEdges[q].size();
                int nVrts  = faceVerts[q].size();
                
                for(int u=0;u<nVrts;u++)
                {
                    std::cout << "(" << faceVerts[q][u] << ", " << itmesh->second->globV2locV[faceVerts[q][u]] << ") ";
                }
                std::cout << std::endl;
            }
        }
        
        for(int q=0;q<nfaces;q++)
        {
            int nEdges = faceEdges[q].size();
            int nVrts  = faceVerts[q].size();
            std::vector<int> face2edge(nEdges);
            std::set<int> face2edge_set;
            std::set<int> faceverts_set;
            std::set<int> faceOverts_set;
            std::vector<std::vector<double> > faceCoords;
            
            for(int p=0;p<nVrts;p++)
            {
                
                int ovid = new2old_v[vert_map_inv[faceVerts[q][p]]];
                faceOverts_set.insert(ovid);
                
                std::vector<double> coords(3);
                coords[0] = us3d->xcn->getVal(ovid,0);
                coords[1] = us3d->xcn->getVal(ovid,1);
                coords[2] = us3d->xcn->getVal(ovid,2);
                
                faceCoords.push_back(coords);
            }
            
            
            std::map<int,int> gEdge2lEdgeMap;
            for(int p=0;p<nEdges;p++)
            {
                std::set<int> edge;
                std::vector<int> edge_vec;
                for(int r=0;r<2;r++)
                {
                    edge.insert(faceEdges[q][p][r]);
                    
                }
                
                int geid = edge_set[edge];
                
//                if(elid == 28 || elid == 30)
                if(itmesh->second->nodes.size()==6)
                {
                    std::cout << "Element " << elid << " PrismTransferMap_Inv " << PrismTransferMap_Inv[elid] << " Face " << q << " :: ";
                    std::cout << "edgie " << geid << " " << p
                    << " --> (" <<  faceEdges[q][p][0] << ", " << itmesh->second->globV2locV[faceEdges[q][p][0]] << ") " << " "
                    << itmesh->second->nodesCoordMap[faceEdges[q][p][0]][0] << " "
                    << itmesh->second->nodesCoordMap[faceEdges[q][p][0]][1] << " "
                    << itmesh->second->nodesCoordMap[faceEdges[q][p][0]][2] << " "
                    << " --> (" <<  faceEdges[q][p][1] << ", " << itmesh->second->globV2locV[faceEdges[q][p][1]] << ") " << " "
                    << itmesh->second->nodesCoordMap[faceEdges[q][p][1]][0] << " "
                    << itmesh->second->nodesCoordMap[faceEdges[q][p][1]][1] << " "
                    << itmesh->second->nodesCoordMap[faceEdges[q][p][1]][2] << " " << std::endl;
                    
                }
                
                
                
                
                face2edge[p] = geid;
                face2edge_set.insert(geid);
                gEdge2lEdgeMap[geid] = p;
            }
//            if(elid == 28 || elid == 30)
            if(itmesh->second->nodes.size()==6)
            {
                std::cout << std::endl;
                
            }
            if(face_set.find(face2edge_set)==face_set.end())
            {
                face_set[face2edge_set] = gfid;
                face_map[gfid] = face2edge;
                gfaces_copy[q] = gfid;
                
                if(refFaceMap.find(faceOverts_set)!=refFaceMap.end())
                {
                    int faceRef = refFaceMap[faceOverts_set];
                    BoundaryComposites[faceRef].push_back(gfid);
                }
                
                std::vector<double> cent = ComputeCentroid(faceCoords);

                
                if(fabs(cent[0]+0.1253763782378859)<1.0e-06 &&
                   fabs(cent[1]+0.4810824626539176)<1.0e-06 &&
                   fabs(cent[2]-0.3333333333333333)<1.0e-06)
                {
                    std::cout << "ELLIE " << elid << " PrismTransferMap_Inv " << PrismTransferMap_Inv[elid] << " " << " " << mesh[elid]->nodes.size() << " " << face2edge[0] << " " << face2edge[1] << " " << face2edge[2] << std::endl;
                    
                    for(int p=0;p<nEdges;p++)
                    {
                        std::cout << p << " :: " << face2edge[p] << " ";
                        for(int r=0;r<2;r++)
                        {
                            std::cout << faceEdges[q][p][r] << " ";
                        }
                        std::cout << std::endl;

                    }
                    std::cout << std::endl;
                }
                
                
                gfid++;
            }
            else
            {
                int gfidn = face_set[face2edge_set];
                std::vector<double> cent = ComputeCentroid(faceCoords);
                
                if(fabs(cent[0]+0.1253763782378859)<1.0e-06 &&
                   fabs(cent[1]+0.4810824626539176)<1.0e-06 &&
                   fabs(cent[2]-0.3333333333333333)<1.0e-06)
                {
                    std::cout << "ELLIE " << elid << " PrismTransferMap_Inv " << PrismTransferMap_Inv[elid] << " " << mesh[elid]->nodes.size() << " " << face2edge[0] << " " << face2edge[1] << " " << face2edge[2] << std::endl;
                    
                    for(int p=0;p<nEdges;p++)
                    {
                        std::cout << p << " :: " << face2edge[p] << " ";
                        for(int r=0;r<2;r++)
                        {
                            std::cout << faceEdges[q][p][r] << " ";
                        }
                        std::cout << std::endl;

                    }
                    std::cout << std::endl;
                }
                
                
                
                gfaces_copy[q] = gfidn;
            }
        }
        
        element_map[elid] = gfaces_copy;
        
        if(itmesh->second->nodes.size()==4)
        {
            bool checkTet = CheckTetRotation(itmesh->second,elid);
            if(!checkTet)
            {
                std::cout << "The orientation of tetrahedron with ID " << elid << " is not correct." << std::endl;
                BadTets++;
            }
        }

        if(itmesh->second->nodes.size()==6)
        {
            bool checkPrism = CheckPrismRotation(itmesh->second,elid);

            if(!checkPrism)
            {
                std::cout << "The orientation of prism with ID " << elid << " is not correct." << std::endl;

                BadPrisms++;
            }
        }
        
    }
    
    // =================== Done setting up the faces ==================

    
    Array<double>* xcn_update = new Array<double>(us3d->xcn->getNrow(),3);
    
    std::map<int,int>::iterator itmp;
    
    for(itmp=vert_map.begin();itmp!=vert_map.end();itmp++)
    {
        int oldvid = new2old_v[itmp->first];
        int newvid = itmp->second;
        
        //std::cout << oldvid << " " << newvid << " " << us3d->xcn->getVal(oldvid,0) << " " << us3d->xcn->getVal(oldvid,1) << " " << us3d->xcn->getVal(oldvid,2) << std::endl;
        
        xcn_update->setVal(newvid,0,us3d->xcn->getVal(oldvid,0));
        xcn_update->setVal(newvid,1,us3d->xcn->getVal(oldvid,1));
        xcn_update->setVal(newvid,2,us3d->xcn->getVal(oldvid,2));

    }
    
//    std::cout << "Done checking Prisms " << badRotPrism << " " << TetWRong << " " << facesDone.size() << std::endl;
//
//    std::cout << "number of edges " << edge_map.size() << std::endl;
//    std::cout << "number of faces " << face_map.size() << " " << face_set.size()  << " " << us3d->ifn->getNrow() << std::endl;
//
//    ReadXmlFile(filename);
    WriteXmlFile("nektarpp.xml",xcn_update,edge_map,element_map,us3d->zdefs,face_map,us3d->ief,BoundaryComposites);
    /**/
    MPI_Finalize();
    
}

