#include "../../src/adapt_recongrad.h"
#include "../../src/adapt_io.h"
#include "../../src/adapt_parops.h"
#include "../../src/adapt_output.h"
#include "../../src/adapt_boundary.h"
#include "../../src/adapt_distri_parstate.h"
#include "../../src/adapt_redistribute.h"
#include "../../src/adapt_DefinePrismMesh.h"
#include "../../src/adapt_prismaticlayer.h"
#include "../../src/NekFace.h"
//#include <boost/core/ignore_unused.hpp>
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
//        // std::cout << "Read succesfully!" << std::endl;
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
                // std::cout << "hiero = " << itelem->second[0] << " " << itelem->second[1] << " " << itelem->second[2] << " " << itelem->second[3] << " " << itelem->second[4] << std::endl;
            }
            if(itelem->second[0]==1914 && itelem->second[1]==1915 && itelem->second[2]==1916 && itelem->second[3]==1907 && itelem->second[4]==1917)
            {
                // std::cout << "YESSIR " << itelem->first << std::endl;
            }
            TiXmlElement *e = new TiXmlElement("R");
            e->SetAttribute("ID",itelem->first);
            TiXmlText *vList = new TiXmlText(s.str().c_str());
            e->LinkEndChild(vList);
            elemTag->LinkEndChild(e);
            prisms.push_back(itelem->first);
        }

    }

    // std::cout << "A " << tets[0] << " " << tets[tets.size()-1] << std::endl;
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

    if(prisms.size()!=0)
    {
        // std::cout << "R " << prisms[0] << " " << prisms[prisms.size()-1] << std::endl;

        TiXmlElement *cP = new TiXmlElement("C");

        std::stringstream ssP;

        ssP << " R[" << prisms[0] << "-" << prisms[prisms.size()-1] << "]";

        cP->SetAttribute("ID",1);
        TiXmlText *vListP = new TiXmlText(ssP.str().c_str());
        cP->LinkEndChild(vListP);
        compTag->LinkEndChild(cP);
    }
    

    
    std::map<int,std::vector<int> >::iterator biter;
    
    int bcnt = 2;
    for(biter=BoundaryComposites.begin();biter!=BoundaryComposites.end();biter++)
    {
        int bid = biter->first;
        // std::cout << "boundary info " << bid << " " << biter->second.size() << std::endl;
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

    if(prisms.size()!=0)
    {
        std::stringstream ssD;


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
            // std::cout << "Error in "<< doc.Value()<< ": " <<  doc.ErrorDesc() << std::endl;
            exit( 1 );
        }
    }
    else
    {
        std::stringstream ssD;


        ssD << " C[0]";

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
        
        
        root->LinkEndChild(expTag);
        
        if ( doc.Error() )
        {
            // std::cout << "Error in "<< doc.Value()<< ": " <<  doc.ErrorDesc() << std::endl;
            exit( 1 );
        }
    }

    doc.SaveFile();

}





//std::size_t ComputeFaceHash(std::vector<int> Face)
//{
//    sort(Face.begin(),Face.end());
//    std::size_t seed = 0;
//    std::hash<int> hasher;
//    for (int i : Face) {
//        seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
//    }
//    return seed;
//}
//
//
//std::size_t ComputeFaceHash_V2(std::vector<int> Face)
//{
//    sort(Face.begin(),Face.end());
//
//    size_t hash = 0;
//    hash_range(hash,Face.begin(),Face.end());
////    if(Face.size()==3)
////    {
////        hash_combine(hash, Face[0], Face[1], Face[2]);
////    }
////
////    if(Face.size()==4)
////    {
////        hash_combine(hash, Face[0], Face[1], Face[2], Face[3]);
////    }
//
//    return hash;
//}




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

    //// std::cout << "Value = " << ((v[3].x - v[0].x) * abx + (v[3].y - v[0].y) * aby +
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
        
        //// std::cout << m_origVertMap[i] << " ";
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
    
    //// std::cout << std::endl;
    std::map<int,int> glob2loc;
    for(int i=0;i<6;i++)
    {
        glob2loc[tet->nodes[i]] = i;
    }
//    std::map<int,std::vector<double> > nodesCoordMap;
//    std::vector<std::vector<double> > nodesCoords_copy;
////    // std::cout << glob2loc[nodes[0]] << " " <<  glob2loc[nodes[1]]  << " " <<  glob2loc[nodes[2]]  << " " <<  glob2loc[nodes[3]]  << " " <<  glob2loc[nodes[4]]  << " " <<  glob2loc[nodes[5]] << std::endl;
//    std::vector<double> node0(3);
//    //// std::cout << "Whats up tet " << id << " " << nodesCoords.size() << " " << glob2loc[nodes[0]] << " " << nodes[0] << std::endl;
//    node0[0] = nodesCoords[glob2loc[nodes[0]]][0];
//    node0[1] = nodesCoords[glob2loc[nodes[0]]][1];
//    node0[2] = nodesCoords[glob2loc[nodes[0]]][2];
//    nodesCoords_copy.push_back(node0);
//    std::vector<double> node1(3);
//    node1[0] = nodesCoords[glob2loc[nodes[1]]][0];
//    node1[1] = nodesCoords[glob2loc[nodes[1]]][1];
//    node1[2] = nodesCoords[glob2loc[nodes[1]]][2];
//    nodesCoords_copy.push_back(node1);
//    std::vector<double> node2(3);
//    node2[0] = nodesCoords[glob2loc[nodes[2]]][0];
//    node2[1] = nodesCoords[glob2loc[nodes[2]]][1];
//    node2[2] = nodesCoords[glob2loc[nodes[2]]][2];
//    nodesCoords_copy.push_back(node2);
//    std::vector<double> node3(3);
//    node3[0] = nodesCoords[glob2loc[nodes[3]]][0];
//    node3[1] = nodesCoords[glob2loc[nodes[3]]][1];
//    node3[2] = nodesCoords[glob2loc[nodes[3]]][2];
//    nodesCoords_copy.push_back(node3);
//    std::vector<double> node4(3);
//    node4[0] = nodesCoords[glob2loc[nodes[4]]][0];
//    node4[1] = nodesCoords[glob2loc[nodes[4]]][1];
//    node4[2] = nodesCoords[glob2loc[nodes[4]]][2];
//    nodesCoords_copy.push_back(node4);
//    std::vector<double> node5(3);
//    node5[0] = nodesCoords[glob2loc[nodes[5]]][0];
//    node5[1] = nodesCoords[glob2loc[nodes[5]]][1];
//    node5[2] = nodesCoords[glob2loc[nodes[5]]][2];
//    nodesCoords_copy.push_back(node5);
//
//    nodesCoordMap[nodes[0]] = node0;
//    nodesCoordMap[nodes[1]] = node1;
//    nodesCoordMap[nodes[2]] = node2;
//    nodesCoordMap[nodes[3]] = node3;
//    nodesCoordMap[nodes[4]] = node4;
//    nodesCoordMap[nodes[5]] = node5;
    std::map<int,int> globV2locV;
    for(int i=0;i<tet->nodes.size();i++)
    {
        globV2locV[nodes[i]] = i;
    }
    tet->globV2locV = globV2locV;
    
    tet->nodesCoordMap = nodeMap;
    return tet;
}



NekElement* SetPrism_V2(std::vector<int> nodes, std::vector<std::vector<double> > nodesCoords, int elid)
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

    for(int i=0;i<9;i++)
    {
        std::vector<int> edge;
        
        for(int j=0;j<2;j++)
        {
            edge.push_back(nodes[m_edgeVertMap[i][j]]);
            
        }
        prism->m_edges.push_back(edge);
    }
    
    for(int i=0;i<6;i++)
    {
        glob2loc[nodes[i]] = i;
    }
    
    // Create faces
    int face_edges[5][4];
   // int face_offset[5];
//    face_offset[0] = 6 + 9 * n;
//    for (int j = 0; j < 4; ++j)
//    {
//        int facenodes      = j % 2 == 0 ? n * n : n * (n - 1) / 2;
//        face_offset[j + 1] = face_offset[j] + facenodes;
//    }
//
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
    
    
    std::map<int,std::vector<double> > nodesCoordMap;
    std::vector<std::vector<double> > nodesCoords_copy;
//    // std::cout << glob2loc[nodes[0]] << " " <<  glob2loc[nodes[1]]  << " " <<  glob2loc[nodes[2]]  << " " <<  glob2loc[nodes[3]]  << " " <<  glob2loc[nodes[4]]  << " " <<  glob2loc[nodes[5]] << std::endl;
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
    
    //prism->m_eorient     = m_eorient;
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


void OrientPrism(std::vector<int> &nodes)
{
    // Create edges (with corresponding set of edge points)
    int m_orientation = 0;
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
}

NekElement* SetPrism(std::vector<int> nodes, std::vector<std::vector<double> > nodesCoords,int elid, int orient)
{
    
    NekElement* prism = new NekElement;
    int m_orientation = 0;
    
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
    
    
    std::map<int,int> glob2loc;

    for(int i=0;i<6;i++)
    {
        glob2loc[nodes[i]] = i;
    }
    
    if(orient==1)
    {
        OrientPrism(nodes);
    }
    
    
    
    
    //std::vector<std::vector<int> > m_edges;
    //// std::cout << "m_edgeVertMap[i][j] " << std::endl;
    for(int i=0;i<9;i++)
    {
        std::vector<int> edge;
        
        for(int j=0;j<2;j++)
        {
            //m_edge[i][j] = nodes[m_edgeVertMap[i][j]];
            edge.push_back(nodes[m_edgeVertMap[i][j]]);
            
            //// std::cout << m_edgeVertMap[i][j] << " ";
            
        }
        
        if(edge[0] == 0 && edge[1] == 0)
        {
            //// std::cout << "ANOTHER ERROR " << " " << elid << std::endl;
            
//            for(int q=0;q<nodes.size();q++)
//            {
//                // std::cout << nodes[q] << std::endl;
//            }
//            // std::cout << std::endl;
        }
        //// std::cout << std::endl;
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
        //// std::cout << "face_edges ";
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
                //// std::cout << std::endl;
            }

            if (i == prism->m_edges.size())
            {
                face_edges[j][k] = -1;
                //// std::cout << face_edges[j][k] << " ";
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
//    // std::cout << "m_orientation " << m_orientation << std::endl;
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
            m_eorient[i] =  1;
        }
        else if(prism->m_edges[i][0] == nodes[m_edgeVertMap[i][1]])
        {
            m_eorient[i] = -1;
        }
    }
    
    std::map<int,std::vector<double> > nodesCoordMap;
    std::vector<std::vector<double> > nodesCoords_copy;

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



NekElement* SetPrismV3(std::vector<int> nodes, std::vector<std::vector<double> > nodesCoords,int elid, int orient)
{
    
    NekElement* prism = new NekElement;
    int m_orientation = 0;
    
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
    
    
    std::map<int,int> glob2loc;

    for(int i=0;i<6;i++)
    {
        glob2loc[nodes[i]] = i;
    }
    
    if(orient==1)
    {
        OrientPrism(nodes);
    }
    
    
    
    
    //std::vector<std::vector<int> > m_edges;
    //// std::cout << "m_edgeVertMap[i][j] " << std::endl;
    for(int i=0;i<9;i++)
    {
        std::vector<int> edge;
        
        for(int j=0;j<2;j++)
        {
            //m_edge[i][j] = nodes[m_edgeVertMap[i][j]];
            edge.push_back(nodes[m_edgeVertMap[i][j]]);
            
            //// std::cout << m_edgeVertMap[i][j] << " ";
            
        }
        
//        if(edge[0] == 0 && edge[1] == 0)
//        {
//            // std::cout << "ANOTHER ERROR " << " " << elid << std::endl;
//
//            for(int q=0;q<nodes.size();q++)
//            {
//                // std::cout << nodes[q] << std::endl;
//            }
//            // std::cout << std::endl;
//        }
        //// std::cout << std::endl;
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
        //// std::cout << "face_edges ";
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
                    //// std::cout << prism->m_edges[i][0] << " " << prism->m_edges[i][1] << " ";
                    
                    if(prism->m_edges[i][0] == 0 && prism->m_edges[i][1] == 0)
                    {
                        // std::cout << "Errror "<< std::endl;
                    }
                    break;
                }
                //// std::cout << std::endl;
            }

            if (i == prism->m_edges.size())
            {
                face_edges[j][k] = -1;
                // std::cout << face_edges[j][k] << " ";
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
//    // std::cout << "m_orientation " << m_orientation << std::endl;
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
//    // std::cout << glob2loc[nodes[0]] << " " <<  glob2loc[nodes[1]]  << " " <<  glob2loc[nodes[2]]  << " " <<  glob2loc[nodes[3]]  << " " <<  glob2loc[nodes[4]]  << " " <<  glob2loc[nodes[5]] << std::endl;
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
    std::map<std::set<int>,std::vector<int> > shellFaces2Element;
    FaceSet faceSet_shell;
    FaceSetPointer FaceSetPointer_shell;
    std::map<int,std::vector<int> > shellF_map;
    std::map<int,int > shellFaces2WallFace;
    std::map<int,int> Element2ShellFace;
    std::map<int,std::vector<int> > WallFace2ShellFace;
    std::map<int,std::vector<int> > WallFace2WallVert;
    std::map<int,int > ShellFace2WallFace;
    std::map<int, std::vector<int> > SinglePrismVertsPerWallFace;
    std::map<int,int> wallFace2PrismID;
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
    int shellFid = 0;
    
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
//                // std::cout << "ftest " << *its << std::endl;
//            }
//
//            for(its=fnewSet0.begin();its!=fnewSet0.end();its++)
//            {
//                // std::cout << "fnewSet0 " << *its << std::endl;
//            }
//
//            for(its=fnewSet1.begin();its!=fnewSet1.end();its++)
//            {
//                // std::cout << "fnewSet1 " << *its << std::endl;
//            }
//
//            // std::cout << "startV " << startV << std::endl;
            
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
            
            //// std::cout << "orientations " << el_cur << " " << orient0 << " " << orient0_test << " " << orient0_next << std::endl;
            
            
            if(orient0_next>0.0)
            {
                // std::cout << "Both face normals are pointing in the same direction " << std::endl;
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
                
                std::vector<int> ShellElem(2);
                
                ShellElem[0] = el_cur;
                ShellElem[1] = nextElemId;
                std::set<int> shellface;
                shellface.insert(us3d->ifn->getVal(nextFaceId,0));
                shellface.insert(us3d->ifn->getVal(nextFaceId,1));
                shellface.insert(us3d->ifn->getVal(nextFaceId,2));
                plines->shellFaces2Element[shellface] = ShellElem;

                
                
                std::vector<int> shellfaceVec(3);
                shellfaceVec[0] = us3d->ifn->getVal(nextFaceId,0);
                shellfaceVec[1] = us3d->ifn->getVal(nextFaceId,1);
                shellfaceVec[2] = us3d->ifn->getVal(nextFaceId,2);
                NekFace f2e(shellfaceVec);
                pair<FaceSet::iterator, bool> testIns;
                testIns = plines->faceSet_shell.insert(f2e);
                
                if(testIns.second)
                {
                    f2e.SetFaceID(shellFid);
                    plines->shellF_map[shellFid] = ShellElem;
                    plines->shellFaces2WallFace[shellFid] = wallfaces[wf];

                    shellFid++;
                }
                
                
//                WallFace2ShellFace[wallfaces[wf]] = nextFaceId;
//                ShellFace2WallFace[shellface]     = wallfaces[wf];
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
        
        plines->ElementLinesNodes[wallfaces[wf]]   = PrismLineNodes;
        plines->ElementLines[wallfaces[wf]]        = PrismLine;
        plines->FaceLines[wallfaces[wf]]           = FaceLine;
        
        //// std::cout << "PrismLine = " << PrismLine.size() << std::endl;
        //// std::cout << std::endl;
    }
    
    plines->ElementLinesNodes = Lnodes;
    
    return plines;
}



PrismLines* GetPrismLines(US3D* us3d,
                          std::vector<int> wallfaces,
                          int nTets,
                          std::map<int,std::set<int> > ref_vertices,
                          std::map<std::set<int>, int> &refFaceMap)
{
    std::ofstream myfile1;
    myfile1.open("SinglePrisms.dat");
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
    int newPrismElmID = nTets;
    std::map<int,std::map<int,std::vector<int> > > Lnodes;
    int shellFid = 0;
    int shellFid2 = 0;
    int LocQuadMap[3][4] = {{0,3,5,2},
                            {2,5,4,1},
                            {1,4,3,0}};
    
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
        
        int tetFound = 0;
        
        std::vector<int> SinglePrismVerts(6);
        
        SinglePrismVerts[0] = us3d->ifn->getVal(wfid,0);
        SinglePrismVerts[1] = us3d->ifn->getVal(wfid,1);
        SinglePrismVerts[2] = us3d->ifn->getVal(wfid,2);
        
//        SinglePrismVerts[0] = us3d->ien->getVal(el_cur,0);
//        SinglePrismVerts[1] = us3d->ien->getVal(el_cur,1);
//        SinglePrismVerts[2] = us3d->ien->getVal(el_cur,2);
        
        std::vector<int> wallfaceVec(3);
        wallfaceVec[0] = us3d->ifn->getVal(wfid,0);
        wallfaceVec[1] = us3d->ifn->getVal(wfid,1);
        wallfaceVec[2] = us3d->ifn->getVal(wfid,2);
        
        while(tetFound==0)
        {
            PrismLine.push_back(el_cur);
            FaceLine.push_back(wfid);
            plines->Prism2WallFace[el_cur] = wfid;
            
            std::vector<int> prismsVerts(6);
            
            prismsVerts[0] = us3d->ien->getVal(el_cur,0);
            prismsVerts[1] = us3d->ien->getVal(el_cur,1);
            prismsVerts[2] = us3d->ien->getVal(el_cur,2);
            prismsVerts[3] = us3d->ien->getVal(el_cur,3);
            prismsVerts[4] = us3d->ien->getVal(el_cur,4);
            prismsVerts[5] = us3d->ien->getVal(el_cur,5);
            
            PrismLineNodes[el_cur] = prismsVerts;
            
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
 
            //// std::cout << "ntype " << ntype << std::endl;
            
            if(ntype==2)
            {
                tetFound = 1;
                
                std::vector<int> ShellElem(2);
                
                ShellElem[0] = el_cur;
                ShellElem[1] = nextElemId;
                std::set<int> shellface;
                shellface.insert(us3d->ifn->getVal(nextFaceId,0));
                shellface.insert(us3d->ifn->getVal(nextFaceId,1));
                shellface.insert(us3d->ifn->getVal(nextFaceId,2));
                plines->shellFaces2Element[shellface] = ShellElem;
                
                std::vector<int> shellfaceVec(3);
//              shellfaceVec[0] = us3d->ien->getVal(el_cur,3);
//              shellfaceVec[1] = us3d->ien->getVal(el_cur,4);
//              shellfaceVec[2] = us3d->ien->getVal(el_cur,5);
                
                shellfaceVec[0] = us3d->ifn->getVal(nextFaceId,0);
                shellfaceVec[1] = us3d->ifn->getVal(nextFaceId,2);
                shellfaceVec[2] = us3d->ifn->getVal(nextFaceId,1);
                
                SinglePrismVerts[3] = us3d->ifn->getVal(nextFaceId,0);
                SinglePrismVerts[4] = us3d->ifn->getVal(nextFaceId,1);
                SinglePrismVerts[5] = us3d->ifn->getVal(nextFaceId,2);
                
                std::vector<double> SinglePrismVertsCopy(3);
                SinglePrismVertsCopy[0] = SinglePrismVerts[3];
                SinglePrismVertsCopy[1] = SinglePrismVerts[4];
                SinglePrismVertsCopy[2] = SinglePrismVerts[5];
                
//                SinglePrismVerts[3] = us3d->ien->getVal(el_cur,3);
//                SinglePrismVerts[4] = us3d->ien->getVal(el_cur,4);
//                SinglePrismVerts[5] = us3d->ien->getVal(el_cur,5);
//
                std::vector<std::vector<double> > Coords;
                for(int j=0;j<6;j++)
                {
                    int vid = SinglePrismVerts[j];
                    
                    std::vector<double> coord(3);
                    
                    coord[0] = us3d->xcn->getVal(vid,0);
                    coord[1] = us3d->xcn->getVal(vid,1);
                    coord[2] = us3d->xcn->getVal(vid,2);
                    
                    Coords.push_back(coord);
                    
                }
                std::vector<double> cent = ComputeCentroid(Coords);
                
                
                std::vector<double> FaceCentroid0(3);
                FaceCentroid0[0] = 0.0;
                FaceCentroid0[1] = 0.0;
                FaceCentroid0[2] = 0.0;
                std::vector<double> FaceCentroid1(3);
                FaceCentroid1[0] = 0.0;
                FaceCentroid1[1] = 0.0;
                FaceCentroid1[2] = 0.0;
                std::vector<std::vector<double> > Face0;
                std::vector<std::vector<double> > Face1;
                for(int p=0;p<3;p++)
                {
                    vector<double> F0Coord(3);
                    F0Coord[0] = us3d->xcn->getVal(SinglePrismVerts[0*3+p],0);
                    F0Coord[1] = us3d->xcn->getVal(SinglePrismVerts[0*3+p],1);
                    F0Coord[2] = us3d->xcn->getVal(SinglePrismVerts[0*3+p],2);
                    FaceCentroid0[0] = FaceCentroid0[0] + us3d->xcn->getVal(SinglePrismVerts[0*3+p],0);
                    FaceCentroid0[1] = FaceCentroid0[1] + us3d->xcn->getVal(SinglePrismVerts[0*3+p],1);
                    FaceCentroid0[2] = FaceCentroid0[2] + us3d->xcn->getVal(SinglePrismVerts[0*3+p],2);
                    Face0.push_back(F0Coord);
                    vector<double> F1Coord(3);
                    F1Coord[0] = us3d->xcn->getVal(SinglePrismVerts[1*3+p],0);
                    F1Coord[1] = us3d->xcn->getVal(SinglePrismVerts[1*3+p],1);
                    F1Coord[2] = us3d->xcn->getVal(SinglePrismVerts[1*3+p],2);
                    FaceCentroid1[0] = FaceCentroid1[0] + us3d->xcn->getVal(SinglePrismVerts[1*3+p],0);
                    FaceCentroid1[1] = FaceCentroid1[1] + us3d->xcn->getVal(SinglePrismVerts[1*3+p],1);
                    FaceCentroid1[2] = FaceCentroid1[2] + us3d->xcn->getVal(SinglePrismVerts[1*3+p],2);
                    Face1.push_back(F1Coord);
                }
                
                //=============TESTING=================
                
                std::vector<std::vector<Vec3D*> > orF0;
                std::vector<std::vector<Vec3D*> > orF1;
                std::vector<double> midist0;
                
                for(int p=0;p<3;p++)
                {
                    double F0vxRef = Face0[p][0];
                    double F0vyRef = Face0[p][1];
                    double F0vzRef = Face0[p][2];
                    
                    double F1vxRef = Face1[p][0];
                    double F1vyRef = Face1[p][1];
                    double F1vzRef = Face1[p][2];
                    std::vector<Vec3D*> oris0;
                    std::vector<Vec3D*> oris1;
                    
                    std::vector<double> dists0(3);
                    std::vector<double> dists1(3);
                    for(int s=0;s<3;s++)
                    {
                        double F1vx = Face1[s][0];
                        double F1vy = Face1[s][1];
                        double F1vz = Face1[s][2];
                        
                        Vec3D* or0 = new Vec3D;
                        or0->c0 = (F0vxRef-F1vx);
                        or0->c1 = (F0vyRef-F1vy);
                        or0->c2 = (F0vzRef-F1vz);
                        dists0[s] = sqrt(or0->c0*or0->c0+or0->c1*or0->c1+or0->c2*or0->c2);
                        
                        double F0vx = Face0[s][0];
                        double F0vy = Face0[s][1];
                        double F0vz = Face0[s][2];
                        
                        Vec3D* or1 = new Vec3D;
                        or1->c0 = (F1vxRef-F0vx);
                        or1->c1 = (F1vyRef-F0vy);
                        or1->c2 = (F1vzRef-F0vz);
                        dists1[s] = sqrt(or1->c0*or1->c0+or1->c1*or1->c1+or1->c2*or1->c2);
                        
                        oris0.push_back(or0);
                        oris1.push_back(or1);
                    }
                    std::vector<double>::iterator result0 = std::min_element(dists0.begin(),dists0.end());
                    int mindist0_index = std::distance(dists0.begin(), result0);
                    std::vector<double>::iterator result1 = std::min_element(dists1.begin(),dists1.end());
                    int mindist1_index = std::distance(dists1.begin(), result1);
                    //// std::cout << "mindist " << p << " " << mindist0_index << " " << mindist1_index << std::endl;
                    orF0.push_back(oris0);
                    orF1.push_back(oris1);
                    midist0.push_back(mindist0_index);
                    
                }
                //// std::cout << std::endl;
                
                for(int p=0;p<3;p++)
                {
                    SinglePrismVerts[3+p] = SinglePrismVertsCopy[midist0[p]];
                }
                
                
                for(int s=0;s<3;s++)
                {
                    int newRefFace = 2;
                    std::map<int,std::vector<int> > refVmap;
                    std::map<int,std::vector<int> >::iterator rit;
                    std::set<int> RebuiltFace;
                    
                    for(int t=0;t<4;t++)
                    {
                        int locid             = LocQuadMap[s][t];
                        int vid               = SinglePrismVerts[locid];
                        RebuiltFace.insert(vid);
                        std::set<int> refV = ref_vertices[vid];
                        std::set<int>::iterator its;
                        
                        for(its=refV.begin();its!=refV.end();its++)
                        {
                            refVmap[*its].push_back(vid);
                        }
                    }

                    int refmax = 2;
                    for(rit=refVmap.begin();rit!=refVmap.end();rit++)
                    {
                        int ref     = rit->first;
                        int refsize = rit->second.size();
                            
//                           if(refsize>2)//The face needs to have more than 2 vertices that hacve a boundary reference attached to it for the face to be considered to be a boundary face.
//                            {
                        if(refsize>refmax && refsize>2)
                        {
                            refmax  = refsize;
                            newRefFace = ref;
                        }
                        
                    }
                    
                    
                    if(refFaceMap.find(RebuiltFace)==refFaceMap.end())
                    {
                        //// std::cout << "newRefFace  " << newRefFace << std::endl;
                        
                        if(newRefFace != 2)
                        {
                            refFaceMap[RebuiltFace] = newRefFace;
                        }
                        
                    }
                    
                    
                    
                }
//
//                for(int p=0;p<3;p++)
//                {
//                    // std::cout << "S " << SinglePrismVerts[3+p] << " " << SinglePrismVertsCopy[p] << std::endl;
//                }
//
                
                //=============TESTING=================
                int fon = 0;
                for(int p=0;p<3;p++)
                {
                    std::vector<Vec3D*> ors0 = orF0[p];
                    
                    for(int r=0;r<3;r++)
                    {
                        Vec3D* orRef = ors0[r];
                        
                        for(int s=0;s<3;s++)
                        {
                            std::vector<Vec3D*> ors1 = orF1[s];
                            
                            for(int t=0;t<3;t++)
                            {
                                Vec3D* orTest  = ors1[t];
                                double orient0 = DotVec3D(orRef,orTest);
                                //// std::cout << "oreint test " << orient0 << std::endl;
//                                if(fabs(orient0) < 1.0e-04)
//                                {
//                                    // std::cout << "Found it" << std::endl;
//                                    fon++;
//                                }
                            }
                        }
                    }
                }
                //// std::cout << "fon  " << fon<< std::endl;
                
                
                
                
                
                
                FaceCentroid0[0] = FaceCentroid0[0]/3;
                FaceCentroid0[1] = FaceCentroid0[1]/3;
                FaceCentroid0[2] = FaceCentroid0[2]/3;
                
                FaceCentroid1[0] = FaceCentroid1[0]/3;
                FaceCentroid1[1] = FaceCentroid1[1]/3;
                FaceCentroid1[2] = FaceCentroid1[2]/3;
                
                Vec3D* r0 = new Vec3D;
                r0->c0 = (FaceCentroid0[0]-cent[0]);
                r0->c1 = (FaceCentroid0[1]-cent[1]);
                r0->c2 = (FaceCentroid0[2]-cent[2]);
                
                Vec3D* r1 = new Vec3D;
                r1->c0 = (FaceCentroid1[0]-cent[0]);
                r1->c1 = (FaceCentroid1[1]-cent[1]);
                r1->c2 = (FaceCentroid1[2]-cent[2]);
                
                Vec3D* v0_f0 = new Vec3D;
                v0_f0->c0 = Face0[1][0]-Face0[0][0];
                v0_f0->c1 = Face0[1][1]-Face0[0][1];
                v0_f0->c2 = Face0[1][2]-Face0[0][2];
                Vec3D* v1_f0 = new Vec3D;
                v1_f0->c0 = Face0[2][0]-Face0[0][0];
                v1_f0->c1 = Face0[2][1]-Face0[0][1];
                v1_f0->c2 = Face0[2][2]-Face0[0][2];
                Vec3D* nbf_f0 = ComputeSurfaceNormal(v0_f0,v1_f0);
                double orient0 = DotVec3D(r0,nbf_f0);
                
                
                Vec3D* v0_f1 = new Vec3D;
                v0_f1->c0 = Face1[1][0]-Face1[0][0];
                v0_f1->c1 = Face1[1][1]-Face1[0][1];
                v0_f1->c2 = Face1[1][2]-Face1[0][2];
                Vec3D* v1_f1 = new Vec3D;
                v1_f1->c0 = Face1[2][0]-Face1[0][0];
                v1_f1->c1 = Face1[2][1]-Face1[0][1];
                v1_f1->c2 = Face1[2][2]-Face1[0][2];
                Vec3D* nbf_f1 = ComputeSurfaceNormal(v0_f1,v1_f1);
                double orient1 = DotVec3D(r1,nbf_f1);
                
                //// std::cout << "orientations " << orient0 << " " << orient1 << std::endl;
//                // std::cout << "(" << us3d->ien->getVal(el_cur,3) << ", " << us3d->ifn->getVal(nextFaceId,0) << ") " << "( " << us3d->ien->getVal(el_cur,4) << ", " << us3d->ifn->getVal(nextFaceId,1) << ") " << "(" << us3d->ien->getVal(el_cur,5) << ", " << us3d->ifn->getVal(nextFaceId,2) << ") " << std::endl;
                
//                SinglePrismVerts[3] = us3d->ifn->getVal(nextFaceId,0);
//                SinglePrismVerts[4] = us3d->ifn->getVal(nextFaceId,1);
//                SinglePrismVerts[5] = us3d->ifn->getVal(nextFaceId,2);
                
               
                for(int m=0;m<6;m++)
                {
                    myfile1 << std::setprecision(16) << us3d->xcn->getVal(SinglePrismVerts[m],0) << " " << us3d->xcn->getVal(SinglePrismVerts[m],1) << " " << us3d->xcn->getVal(SinglePrismVerts[m],2) << std::endl;
                }
                /**/
//                NekFace f2e(shellfaceVec);
//                pair<FaceSet::iterator, bool> testIns;
//                testIns = plines->faceSet_shell.insert(f2e);
//
//                if(testIns.second)
//                {
//                    f2e.SetFaceID(shellFid);
//                    plines->shellF_map[shellFid] = ShellElem;
//
//                    shellFid++;
//                }
                
                
                FaceSharedPtr f2ePointer = std::shared_ptr<NekFace>(new NekFace(shellfaceVec));
                
                pair<FaceSetPointer::iterator, bool> testInsPointer;
                testInsPointer = plines->FaceSetPointer_shell.insert(f2ePointer);
                
                if(testInsPointer.second)
                {
                    (*testInsPointer.first)->SetFaceID(shellFid2);
                    plines->shellFaces2WallFace[shellFid2] = wallfaces[wf];
                    
                    plines->shellF_map[shellFid2] = ShellElem;
                    
                    shellFid2++;
                }
                plines->WallFace2ShellFace[wallfaces[wf]] = shellfaceVec;
                plines->WallFace2WallVert[wallfaces[wf]]  = wallfaceVec;
//              WallFace2ShellFace[wallfaces[wf]] = nextFaceId;
//              ShellFace2WallFace[shellface]     = wallfaces[wf];
                
            }
            
            plines->Face2WallFace[nextFaceId] = wfid;
            std::vector<int> facesPerElement(2);

            facesPerElement[0] = wfid;
            facesPerElement[1] = nextFaceId;

            plines->Element2Faces[el_cur] = facesPerElement;
            plines->Element2Faces[el_cur] = facesPerElement;

            el_cur = nextElemId;
            wfid   = nextFaceId;
            
            
        }
        //// std::cout << "PrismLineNodes " << PrismLineNodes.size() << std::endl;
        plines->ElementLinesNodes[wallfaces[wf]] = PrismLineNodes;
        plines->ElementLines[wallfaces[wf]] = PrismLine;
        plines->FaceLines[wallfaces[wf]] = FaceLine;
        plines->wallFace2PrismID[wallfaces[wf]] = newPrismElmID;
        plines->SinglePrismVertsPerWallFace[wallfaces[wf]] = SinglePrismVerts;

        newPrismElmID++;
        //plines>wallFace2ShellFace[wallfaces[wf]] = shellfaceVec;
        
        
    }
        
    
    // std::cout << " plines->shellFaces2Element " << plines->shellFaces2Element.size() << " " << newPrismElmID-nTets<< std::endl;
    myfile1.close();
    return plines;
}


PrismLines* ConvertPrismLines2PrismLayer(US3D* us3d, std::vector<int> wallfaces)
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
    
    std::map<int,std::map<int,std::vector<int> > > Lnodes;
    int shellFid = 0;
    int shellFid2 = 0;
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
        
        int tetFound = 0;
        
        while(tetFound==0)
        {
            PrismLine.push_back(el_cur);
            FaceLine.push_back(wfid);
            plines->Prism2WallFace[el_cur] = wfid;
            
            std::vector<int> prismsVerts(6);
            
            prismsVerts[0] = us3d->ien->getVal(el_cur,0);
            prismsVerts[1] = us3d->ien->getVal(el_cur,1);
            prismsVerts[2] = us3d->ien->getVal(el_cur,2);
            prismsVerts[3] = us3d->ien->getVal(el_cur,3);
            prismsVerts[4] = us3d->ien->getVal(el_cur,4);
            prismsVerts[5] = us3d->ien->getVal(el_cur,5);
            
            PrismLineNodes[el_cur] = prismsVerts;
            
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
                
                std::vector<int> ShellElem(2);
                
                ShellElem[0] = el_cur;
                ShellElem[1] = nextElemId;
                std::set<int> shellface;
                shellface.insert(us3d->ifn->getVal(nextFaceId,0));
                shellface.insert(us3d->ifn->getVal(nextFaceId,1));
                shellface.insert(us3d->ifn->getVal(nextFaceId,2));
                plines->shellFaces2Element[shellface] = ShellElem;
                
                std::vector<int> shellfaceVec(3);
                shellfaceVec[0] = us3d->ifn->getVal(nextFaceId,0);
                shellfaceVec[1] = us3d->ifn->getVal(nextFaceId,1);
                shellfaceVec[2] = us3d->ifn->getVal(nextFaceId,2);
//                NekFace f2e(shellfaceVec);
//                pair<FaceSet::iterator, bool> testIns;
//                testIns = plines->faceSet_shell.insert(f2e);
//
//                if(testIns.second)
//                {
//                    f2e.SetFaceID(shellFid);
//                    plines->shellF_map[shellFid] = ShellElem;
//
//                    shellFid++;
//                }
                
                
                FaceSharedPtr f2ePointer = std::shared_ptr<NekFace>(new NekFace(shellfaceVec));
                
                pair<FaceSetPointer::iterator, bool> testInsPointer;
                testInsPointer = plines->FaceSetPointer_shell.insert(f2ePointer);
                
                if(testInsPointer.second)
                {
                    (*testInsPointer.first)->SetFaceID(shellFid2);
                    plines->shellFaces2WallFace[shellFid2] = wallfaces[wf];

                    plines->shellF_map[shellFid2] = ShellElem;
                    
                    shellFid2++;
                }
                
                
                
//                WallFace2ShellFace[wallfaces[wf]] = nextFaceId;
//                ShellFace2WallFace[shellface]     = wallfaces[wf];
            }
            
            plines->Face2WallFace[nextFaceId] = wfid;
            std::vector<int> facesPerElement(2);
            //// std::cout << "Face2Element " << el_cur << " " << nextElemId << std::endl;
            facesPerElement[0] = wfid;
            facesPerElement[1] = nextFaceId;
            plines->Element2Faces[el_cur] = facesPerElement;
            el_cur = nextElemId;
            wfid   = nextFaceId;
            
            
        }
        //// std::cout << "PrismLineNodes " << PrismLineNodes.size() << std::endl;
        plines->ElementLinesNodes[wallfaces[wf]] = PrismLineNodes;
        plines->ElementLines[wallfaces[wf]] = PrismLine;
        plines->FaceLines[wallfaces[wf]] = FaceLine;
        
    }
    // std::cout << " plines->shellFaces2Element " << plines->shellFaces2Element.size() << std::endl;
    return plines;
}


//struct FaceMap{
//    std::map<int,std::vector<int> > Faces;
//    std::map<int,std::vector<int> > collisions;
//};
//
//void HashFaceToFaceMap(FaceMap* &fmap, std::vector<int> face)
//{
//    int dupl = 0;
//    int addnew = 0;
//    int nv = face.size();
//    std::vector<int> tmp(nv);
//    for(int i=0;i<nv;i++)
//    {
//        tmp[i] = face[i];
//
//    }
//    sort(tmp.begin(),tmp.end());
//
//    size_t HA = ComputeFaceHash_V2(tmp);
//
//    if(fmap->Faces.find(HA) == fmap->Faces.end())
//    {
//        fmap->Faces[HA] = face;
//    }
//    else
//    {
//        std::vector<int> triCom = fmap->Faces[HA];
//        int err = 0;
//        for(int i = 0;i < face.size();i++)
//        {
//            err = err + (face[i]-triCom[i]);
//        }
//
//        if(err==0)
//        {
//            dupl++;
//        }
//        else
//        {
//            if(fmap->collisions.find(HA)==fmap->collisions.end())
//            {
//                fmap->collisions[HA] = face;
//                addnew++;
//            }
//        }
//    }
//}

//std::vector<int> GetFaceFromMap(FaceMap* fmap, size_t key)
//{
//    std::vector<int> face;
//
//    if(fmap->Faces.find(key)!=fmap->Faces.end())
//    {
//
//
//    }
//
//
//    if(fmap->collisions.find(key)!=fmap->collisions.end())
//    {
//        face = fmap->collisions[key];
//    }
//    else
//    {
//        if(fmap->Faces.find(key)!=fmap->Faces.end())
//        {
//            face = fmap->Faces[key];
//        }
//        else
//        {
//            // std::cout << "Warning: " << key << " is not in FaceMap" << std::endl;
//        }
//    }
//
//    return face;
//}

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
        //// std::cout << "f0 " << f0[0] << " " << f0[1] << " " << f0[2] << std::endl;
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
        //// std::cout << "f1 " << f1[0] << " " << f1[1] << " " << f1[2] << std::endl;
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
        
        //// std::cout << "Failed 1" << std::endl;

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
                    std::cout << "More than two vertices do not match " << "triangular face" << endl;
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
        //// std::cout << "Failed 2 " << std::endl;

        // Finally check for top vertex
        for (int i = 0; i < 3; ++i)
        {
            //// std::cout << "Failed 2.25111 "<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;

            if (tF[i] != indx0 && tF[i] != indx1 &&
                tF[i] != indx2)
            {
                indx4 = tF[i];
                
                break;
            }
        }
        //// std::cout << "Failed 2.25111 22222"<< " " << quadFaces.size() << " " << triFaces.size() << " " << coordmap.size() << " " << indx0 << " " << indx1 << " " << indx2 << " " << indx3 << std::endl;
        std::vector<double> a(3);
        a[0] = coordmap[indx1][0] - coordmap[indx0][0];
        a[1] = coordmap[indx1][1] - coordmap[indx0][1];
        a[2] = coordmap[indx1][2] - coordmap[indx0][2];
        
        //// std::cout << "Failed 2.25 "<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;

        std::vector<double> b(3);
        a[0] = coordmap[indx4][0] - coordmap[indx0][0];
        a[1] = coordmap[indx4][1] - coordmap[indx0][1];
        a[2] = coordmap[indx4][2] - coordmap[indx0][2];
        
        std::vector<double> c(3);
        c[0] = coordmap[indx2][0] - coordmap[indx0][0];
        c[1] = coordmap[indx2][1] - coordmap[indx0][1];
        c[2] = coordmap[indx2][2] - coordmap[indx0][2];
        
        std::vector<double> acurlb = NodeCurl(a,b);
        //// std::cout << "Failed 2.5 "<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;

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
        
        //// std::cout << "Failed 3"<< " " << quadFaces.size() << " " << triFaces.size() << std::endl;
        
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
        
        //// std::cout << "Failed 4"<<std::endl;

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
                                           PrismLines* plines,
                                           std::map<int,int> &vertex_map,
                                           std::map<int,int> &vertex_map_inv,
                                           std::map<int,std::vector<std::vector<int> > > &ElementFace2NodeMap,
                                           std::map<int,std::set<int> > ref_vertices_update,
                                           std::map<int,std::set<int> > &ref_vertices_update2)
{
    
    std::map<int,std::vector<int> > ElementNodes;
    
    //std::map<int,std::vector<std::vector<int> > > ElementFace2NodeMap;
    int nGlobElm = us3d->ien->getNrow();
    
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


    // std::cout << "Starting off resetting the vertices... " << Element2Nodes_reset.size() << std::endl;
    

    for(itm=plines->SinglePrismVertsPerWallFace.begin();itm!=plines->SinglePrismVertsPerWallFace.end();itm++)
    {
        int nprismOnLine = itm->second.size();
        std::vector<int> PrismRead = itm->second;
        int elId = plines->wallFace2PrismID[itm->first];
        std::vector<int> nodes(6);
        std::vector<int> orig_nodes(6);
            
//      orig_nodes[0] = PrismRead[3];
//      orig_nodes[1] = PrismRead[4];
//      orig_nodes[2] = PrismRead[1];
//      orig_nodes[3] = PrismRead[0];
//      orig_nodes[4] = PrismRead[5];
//      orig_nodes[5] = PrismRead[2];
        
        orig_nodes[0] = Element2Nodes_reset[elId][3];
        orig_nodes[1] = Element2Nodes_reset[elId][4];
        orig_nodes[2] = Element2Nodes_reset[elId][1];
        orig_nodes[3] = Element2Nodes_reset[elId][0];
        orig_nodes[4] = Element2Nodes_reset[elId][5];
        orig_nodes[5] = Element2Nodes_reset[elId][2];
        
//        if(elId==55736)
//        {
//            // std::cout << "wrong ";
//            for(int u=0;u<6;u++)
//            {
//                // std::cout << orig_nodes[u] << " ";
//            }
//            // std::cout << std::endl;
//        }
        
//        // std::cout << "++++++++++++++++++" << std::endl;
//        // std::cout << "prism 1 " << PrismRead[3] << " " << PrismRead[4] << " " << PrismRead[1] << " " << PrismRead[0] << " " << PrismRead[5] << " " << PrismRead[2] << std::endl;
//
//        // std::cout << "prism 2 " << Element2Nodes_reset[elId][3] << " " << Element2Nodes_reset[elId][4] << " " << Element2Nodes_reset[elId][1] << " " << Element2Nodes_reset[elId][0] << " " << Element2Nodes_reset[elId][5] << " " << Element2Nodes_reset[elId][2] << std::endl;
//        // std::cout << "===================" << std::endl;
 
        std::map<int,std::vector<double> > coordmap;
        std::map<int,std::vector<double> > coordmap_orig;
        std::vector<std::vector<double> > Coords_orig;
        int cntz = 0;
        for(int v=0;v<6;v++)
        {
            
            int vertexid = new2old_v[orig_nodes[v]];
            //// std::cout << vertexid << ", " << orig_nodes[v] << " ";
//            if(elId==55736)
//            {
//                // std::cout << vertexid << " ";
//            }
            
            
            if(orig_nodes[v]==0)
            {
                cntz++;
            }
            std::vector<double> Coords(3);
            Coords[0] = us3d->xcn->getVal(vertexid,0);
            Coords[1] = us3d->xcn->getVal(vertexid,1);
            Coords[2] = us3d->xcn->getVal(vertexid,2);
            Coords_orig.push_back(Coords);
            coordmap_orig[orig_nodes[v]] = Coords;
        }
        //// std::cout << std::endl;
        
        if(cntz>=2)
        {
            // std::cout << "WRONG " <<elId <<std::endl;
            for(int v=0;v<6;v++)
            {
                // std::cout << orig_nodes[v] << " " << Element2Nodes_reset[elId][v] << std::endl;

            }
        }
        
//        if(elId==55736)
//        {
//            // std::cout << std::endl;
//        }
            
        NekElement* prismOrig = SetPrismV3(orig_nodes,Coords_orig,elId,0);
        std::vector<std::vector<int> > ElementFace2Node = prismOrig->m_faceVertices;
        std::map<int,std::vector<double> > coordMap     = prismOrig->nodesCoordMap;


        std::vector<int> NodesReOrder = prismOrig->nodes;
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
            // std::cout << "ERROR: Both faces are process already..." << std::endl;
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
                        ref_vertices_update2[nodeId] = ref_vertices_update[vid];
                        //// std::cout << "ef ";
                        std::set<int>::iterator its;
//                        for(its=ref_vertices_update[vid].begin();its!=ref_vertices_update[vid].end();its++)
//                        {
//                            // std::cout << *its << " ";
//                        }
//                        // std::cout << std::endl;
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
                    ref_vertices_update2[nodeId] = ref_vertices_update[vid];
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
                    ref_vertices_update2[nodeId] = ref_vertices_update[vid];
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
        
        ElementNodes[nGlobElm+prismtel0] = nodes;

        prismtel0++;
    }

    
    // fill in any unset nodes at from other shapes
    std::map<int,int>::iterator itelm2;
    for(itelm2=TetTransferMap.begin();itelm2!=TetTransferMap.end();itelm2++)
    {
        int elId = itelm2->second;
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
        
        int oldid = itelm2->second;
        
        for(int t=0;t<4;t++)
        {
            int vid = tnodes[t];
            if(vertex_map.find(vid)==vertex_map.end())
            {
                vertex_map[vid]         = nodeId;
                vertex_map_inv[nodeId]  = vid;
                nodes[t]                = nodeId;
                ref_vertices_update2[nodeId] = ref_vertices_update[vid];
                nodeId++;
            }
            else
            {
                int nodeIdn = vertex_map[vid];
                nodes[t] = nodeIdn;
            }
        }
        
        ElementNodes[oldid]        = nodes;
        
    }
    /* */
    return ElementNodes;
}


std::map<int,std::vector<int> > ActualResetNodes(US3D* us3d, PrismLines* plines,std::map<int,int> TetTransferMap,std::map<int,int> &new2old_v,std::map<int,int> &old2new_v, std::map<int,std::set<int> > ref_vertices, std::map<int,std::set<int> > &ref_vertices_update)
{
    std::map<int,int>::iterator itr;
    std::map<int,std::vector<int> >::iterator itrvec;
    std::map<int,std::vector<int> > ElementNodeReset;
    std::map<int,int> vreset;
    std::map<int,int> vreset_inv;
    int vidn = 0;
    int nTets = TetTransferMap.size();
    
    for(itrvec =plines->SinglePrismVertsPerWallFace.begin();
        itrvec!=plines->SinglePrismVertsPerWallFace.end();
        itrvec++)
    {
        std::vector<int> nodes(6);
        
        for(int v=0;v<6;v++)
        {
            //int vid = us3d->ien->getVal(elId,v);
            int vid = itrvec->second[v];
            
            if(vreset.find(vid)==vreset.end())
            {
                vreset[vid]               = vidn;
                new2old_v[vidn]           = vid;
                old2new_v[vid]            = vidn;
                nodes[v]                  = vidn;
                ref_vertices_update[vidn] = ref_vertices[vid];
                vidn++;
            }
            else
            {
                int vidn2 = vreset[vid];
                nodes[v]  = vidn2;
            }
        }
        int pelid = plines->wallFace2PrismID[itrvec->first];
        ElementNodeReset[pelid]=nodes;
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
        
        int oldid = itelm2->second;
        
        for(int t=0;t<4;t++)
        {
            int vid = tnodes[t];
            if(vreset.find(vid)==vreset.end())
            {
                vreset[vid]               = vidn;
                new2old_v[vidn]           = vid;
                old2new_v[vid]            = vidn;
                nodes[t]                  = vidn;
                ref_vertices_update[vidn] = ref_vertices[vid];
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


//struct NekElem{
//    std::vector<FaceSharedPtr>;
//};


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
    
    
    int m_edgeVertMap[9][2] = {{0, 1},
                               {1, 2},
                               {3, 2},
                               {0, 3},
                               {0, 4},
                               {1, 4},
                               {2, 5},
                               {3, 5},
                               {4, 5}};
    
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
    
    const char* fn_grid     =    "inputs_new/grid.h5";
    const char* fn_conn     =    "inputs_new/conn.h5";
    //const char* fn_metric   =    "inputs_new/metric.inp";
    
    //std::vector<double> metric_inputs = ReadMetricInputs(fn_metric);
    
//    int ReadFromStats    = metric_inputs[4];

    US3D* us3d                        = ReadUS3DGrid(fn_conn,fn_grid,0,comm,info);
    
    int nEl = us3d->iet->getNrow();
    
    

//    FaceSet m_ShellFaceSet;

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
    
    std::map<int,std::set<int> > ref_vertices;
    
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
                
                ref_vertices[us3d->ifn->getVal(i,j)].insert(fref);
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
//                pair<FaceSet::iterator, bool> testIns;
//                testIns = m_ShellFaceSet.insert(faceVec);
//
//                if(testIns.second)
//                {
//                    facetel++;
//                }
            }
        }
//
//        pair<FaceSet::iterator, bool> testIns;
//        testIns = m_mesh->m_faceSet.insert(elmt[i]->GetFace(j));
//
        if(fref == 3 || fref == 14)
        {
            wallfaces.push_back(i);
        }
    }
    
    std::map<int,int> TetTransferMap;
    std::map<int,int> TetTransferMap_Inv;
    int tetcnt = 0;
    int nPrisms = wallfaces.size();
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
    
    // std::cout << "PrismTransferMap  " << PrismTransferMap.size() << " " << TetTransferMap.size() << std::endl;
    // std::cout << "TetTransferMap  " << TetTransferMap.size() << " " << TetTransferMap.size() << std::endl;
    
    std::map<int,int> vert_map;
    std::map<int,int> vert_map_inv;
    std::map<int,std::vector<std::vector<int> > > ElementFace2NodeMap;
    std::map<int,int> new2old_v;
    std::map<int,int> old2new_v;
    std::map<int,std::set<int> > ref_vertices_update;
    
    int nTets = TetTransferMap.size();
    PrismLines* plines = GetPrismLines(us3d,wallfaces,nTets,ref_vertices,refFaceMap);
    
    std::map<int,std::vector<int> > ResettedElementNodes = ActualResetNodes(us3d,plines,TetTransferMap,new2old_v,old2new_v,ref_vertices, ref_vertices_update);
    std::map<int,std::set<int> > ref_vertices_update2;
    std::map<int,std::vector<int> > ElementNodes = ResetNodes(us3d,ResettedElementNodes,new2old_v,TetTransferMap,plines,vert_map,vert_map_inv,ElementFace2NodeMap,ref_vertices_update,ref_vertices_update2);
    
    // std::cout << "old2new_v " << old2new_v.size() << " " << vert_map.size() << std::endl;
    // Printing
    std::vector<std::vector<int> > u_tris_loc;
    std::set<int> unique_shell_vert;
    std::vector<int>  U_shell_vert;
    std::map<int,int> glob2loc_shell_vert;
    std::map<int,int> loc2glob_shell_vert;
    std::vector<int> u_tri_vec(3);
    std::vector<int> u_tri_vec_loc(3);
    int bb = 0;
    int teller = 0;
    int loc_s_v=0;
    std::map<int,std::vector<int> >::iterator itm2;
    // std::cout << "Starting off resetting the vertices..." << std::endl;
    for(itm2=plines->WallFace2ShellFace.begin();itm2!=plines->WallFace2ShellFace.end();itm2++)
    {
        int bb = 0;
        //std::cout << bb " -- ";
        for(int q=0;q<itm2->second.size();q++)
        {
            int ve = itm2->second[q];
            u_tri_vec[bb]=ve;
            
            //std::cout << ve << ", ";

            if(unique_shell_vert.find(ve)==unique_shell_vert.end())
            {
                //u_tri_vec_loc.push_back(loc_s_v);
                unique_shell_vert.insert(ve);
                U_shell_vert.push_back(ve);
                loc2glob_shell_vert[loc_s_v]=ve;
                glob2loc_shell_vert[ve]=loc_s_v;
                u_tri_vec_loc[bb]= loc_s_v;
                loc_s_v++;
            }
            else
            {
                int local_s_vId = glob2loc_shell_vert[ve];
                u_tri_vec_loc[bb]=local_s_vId;
            }
            bb++;
        }
        //std::cout << std::endl;

        u_tris_loc.push_back(u_tri_vec_loc);
        teller++;
    }
    
    std::ofstream myfile_shell;
    myfile_shell.open("shell.dat");
    myfile_shell << "TITLE=\"new_volume.tec\"" << std::endl;
    myfile_shell <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
    myfile_shell <<"ZONE N = " << U_shell_vert.size() << ", E = " << teller << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
    for(int i=0;i<U_shell_vert.size();i++)
    {
        myfile_shell << us3d->xcn->getVal(loc2glob_shell_vert[i],0)
        << " " << us3d->xcn->getVal(loc2glob_shell_vert[i],1)
        << " " << us3d->xcn->getVal(loc2glob_shell_vert[i],2) << std::endl;
    }
    for(int i=0;i<teller;i++)
    {
        myfile_shell << u_tris_loc[i][0]+1
        << " " << u_tris_loc[i][1]+1
        << " " << u_tris_loc[i][2]+1 << std::endl;
    }
    myfile_shell.close();
    
    // std::cout << "Done writing shell" << std::endl;
    
    
    
    
    std::vector<std::vector<int> > u_tris_loc_wall;
    std::set<int> unique_wall_vert;
    std::vector<int>  U_wall_vert;
    std::map<int,int> glob2loc_wall_vert;
    std::map<int,int> loc2glob_wall_vert;

    bb = 0;
    teller = 0;
    loc_s_v=0;
    std::map<int,std::vector<int> >::iterator itm3;
    // std::cout << "Starting off resetting the vertices..." << std::endl;
    for(itm3=plines->WallFace2WallVert.begin();itm3!=plines->WallFace2WallVert.end();itm3++)
    {
        
        int bb = 0;
        for(int q=0;q<itm3->second.size();q++)
        {
            int ve = itm3->second[q];
            u_tri_vec[bb]=ve;
            
            //// std::cout << "ve = " << ve << std::endl;

            if(unique_wall_vert.find(ve)==unique_wall_vert.end())
            {
                //u_tri_vec_loc.push_back(loc_s_v);
                unique_wall_vert.insert(ve);
                U_wall_vert.push_back(ve);
                loc2glob_wall_vert[loc_s_v]=ve;
                glob2loc_wall_vert[ve]=loc_s_v;
                u_tri_vec_loc[bb]= loc_s_v;
                loc_s_v++;
            }
            else
            {
                int local_s_vId = glob2loc_wall_vert[ve];
                u_tri_vec_loc[bb]=local_s_vId;
            }
            bb++;
        }
        
        u_tris_loc_wall.push_back(u_tri_vec_loc);
        teller++;
        
    }
    
//    std::ofstream myfile_wall;
//    myfile_wall.open("wall.dat");
//    myfile_wall << "TITLE=\"new_volume.tec\"" << std::endl;
//    myfile_wall <<"VARIABLES = \"X\", \"Y\", \"Z\"" << std::endl;
//    myfile_wall <<"ZONE N = " << U_wall_vert.size() << ", E = " << teller << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE" << std::endl;
//    for(int i=0;i<U_wall_vert.size();i++)
//    {
//        myfile_wall << us3d->xcn->getVal(loc2glob_wall_vert[i],0)
//        << " " << us3d->xcn->getVal(loc2glob_wall_vert[i],1)
//        << " " << us3d->xcn->getVal(loc2glob_wall_vert[i],2) << std::endl;
//    }
//    for(int i=0;i<teller;i++)
//    {
//        myfile_wall << u_tris_loc_wall[i][0]+1
//        << " " << u_tris_loc_wall[i][1]+1
//        << " " << u_tris_loc_wall[i][2]+1 << std::endl;
//    }
//    myfile_wall.close();
    
    //// std::cout << "Done writing wall" << std::endl;
    
    
    
    std::map<int,NekElement*> mesh;
    std::map<int,NekElement*> prisms;
    std::map<int,NekElement*> tetrahedra;
    
    std::map<int,std::vector<int> >::iterator itelem;
    for(itelem=ElementNodes.begin();itelem!=ElementNodes.end();itelem++)
    {
        if(itelem->second.size()==4)
        {
            int newTetID = itelem->first;
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
            
            tetrahedra[newTetID] = tet;
        }
        if(itelem->second.size()==6)
        {
            int newPrismID = itelem->first;
            
            std::vector<int> Nodes = itelem->second;
            std::vector<std::vector<double> > coords;
            std::map<int,std::vector<double> > coordMap;
            //// std::cout << "Nodes ";
//            if(newPrismID == 55667)
//            {
//                // std::cout << "gaatie "<< std::endl;
//            }
            for(int i=0;i<Nodes.size();i++)
            {
                //// std::cout << Nodes[i] << " ";
                int ovid = new2old_v[vert_map_inv[Nodes[i]]];
                
                //// std::cout << ovid << " ";
                //if(newPrismID == 55667)
                //{
                    //// std::cout << Nodes[i] << " " << ovid << std::endl;
                //}
                std::vector<double> crd(3);
                crd[0] = us3d->xcn->getVal(ovid,0);
                crd[1] = us3d->xcn->getVal(ovid,1);
                crd[2] = us3d->xcn->getVal(ovid,2);
                coordMap[Nodes[i]] = crd;
                coords.push_back(crd);
            }
            //// std::cout << std::endl;
//            if(newPrismID == 55667)
//            {
                //// std::cout << std::endl;
//            }
            //// std::cout << std::endl;
            NekElement* prism = SetPrism(Nodes,coords,newPrismID,1);
            
            mesh[newPrismID] = prism;
            //prisms[newPrismID] = prism;
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
    
//    for(itmesh=mesh.begin();itmesh!=mesh.end();itmesh++)
//    {
//        NekElement* elem = itmesh->second;
//        int nid = itmesh->first;
//
//        if(elem->nodes.size()==4)
//        {
//            std::vector<std::vector<int> > edges = elem->m_edges;
//            std::map<int,int> Tet_l2g_edge;
//
//            for(int s=0;s<6;s++)
//            {
//                std::set<int> eset;
//                for(int k=0;k<2;k++)
//                {
//                    eset.insert(edges[s][k]);
//                }
//                if(edge_set.find(eset)==edge_set.end())
//                {
//                    edge_set[eset]  = eId;
//                    edge_map[eId]   = edges[s];
//                    Tet_l2g_edge[s] = eId;
//                    eId++;
//                }
//                else
//                {
//                    int eIdn        = edge_set[eset];
//                    Tet_l2g_edge[s] = eIdn;
//                }
//            }
//            el2_l2g_edge[nid] = Tet_l2g_edge;
//
//        }
//        if(elem->nodes.size()==6)
//        {
//
//            std::vector<std::vector<int> > edges = elem->m_edges;
//            std::map<int,int> Prism_l2g_edge;
//
//            for(int s=0;s<9;s++)
//            {
//                std::set<int> edgeset;
//                for(int k=0;k<2;k++)
//                {
//                    edgeset.insert(edges[s][k]);
//                }
//
//                if(edge_set.find(edgeset)==edge_set.end())
//                {
//                    edge_set[edgeset]  = eId;
//                    edge_map[eId]      = edges[s];
//                    Prism_l2g_edge[s]  = eId;
//                    eId++;
//                }
//                else
//                {
//                    int eIdn           = edge_set[edgeset];
//                    Prism_l2g_edge[s]  = eIdn;
//                }
//            }
//
//            el2_l2g_edge[nid] = Prism_l2g_edge;
//        }
//    }
    
    
    for(itmesh=tetrahedra.begin();itmesh!=tetrahedra.end();itmesh++)
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
    teller = 0;
    FaceSet m_FaceSet;
    FaceSetPointer m_FaceSetPointer;
    int cntshell = 0;
    std::map<int,std::vector<int> > face_map2;
    std::map<int,int> WallFace2PrismElem;
    std::map<int,std::vector<int> > shellMapOrientVerts;
    std::map<int,std::vector<int> > shellMapOrientEdges;
    int fidid = 0;
    int fcid = 0;
    int FId = 0;
    std::map<int,NekFace*> face_map4;
    std::map<int,NekFace*> face_map5;
    std::set<std::vector<int> > face_map3;
    int dupl = 0;
    int cntdu = 0;
    std::set<int> hashset;
    int foundd =0;
    int facF = 0;
    for(itmesh=tetrahedra.begin();itmesh!=tetrahedra.end();itmesh++)
    {
        int elid = itmesh->first;

        std::vector<std::vector<std::vector<int> > > faceEdges = itmesh->second->m_faceEdges;
        std::vector<std::vector<int> > faceVerts = itmesh->second->m_faceVertices;

        int nfaces = faceEdges.size();
        std::vector<int> gfaces_copy(nfaces);

        for(int q=0;q<nfaces;q++)
        {
            int nEdges = faceEdges[q].size();
            int nVrts  = faceVerts[q].size();
            std::vector<int> face2edge(nEdges);
            std::vector<int> face2edge_copy(nEdges);
            std::set<int> face2edge_set;
            std::set<int> faceverts_set;
            std::set<int> faceOverts_set;
            std::vector<int> faceOverts_vec;
            std::vector<std::vector<double> > faceCoords;

            for(int p=0;p<nVrts;p++)
            {
                int ovid = new2old_v[vert_map_inv[faceVerts[q][p]]];
                faceOverts_set.insert(ovid);
                faceOverts_vec.push_back(ovid);
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


                face2edge[p] = geid;
                face2edge_copy[p] = geid;
                face2edge_set.insert(geid);
                gEdge2lEdgeMap[geid] = p;
            }
            
//            // std::cout << "face2edge " << face2edge.size() << std::endl;
            NekFace* nf = new NekFace(face2edge);
            int fhash   = nf->GetFaceHash();
//            // std::cout << "fhas " << fhash << std::endl;
            
            if(face_map2.find(fhash)==face_map2.end())
            {
                face_map2[fidid] = face2edge;
                nf->SetFaceID(fidid);
                fidid++;
            }
            
            if(face_map4.find(fhash)==face_map4.end())
            {
                face_map4[fhash] = nf;
                nf->SetFaceID(FId);
                face_map5[FId] = nf;
                FId++;
            }
            else
            {
                NekFace* nfc = face_map4[fhash];
                
                std::vector<int> f_c1(nf->GetEdgeIDs().size());
                std::vector<int> f_c2(nfc->GetEdgeIDs().size());
                
                for(int u=0;u<nf->GetEdgeIDs().size();u++)
                {
                    f_c1[u] = nf->GetEdgeIDs()[u];
                    f_c2[u] = nfc->GetEdgeIDs()[u];
                }
                
                sort(f_c1.begin(),f_c1.end());
                sort(f_c2.begin(),f_c2.end());
                
                if((f_c1[0]==f_c2[0]) &&
                   (f_c1[1]==f_c2[1]) &&
                   (f_c1[2]==f_c2[2]))
                {
                    dupl++;
                }
                else
                {
//                    // std::cout << " 1 " << f_c1[0] << " " << f_c1[1] << " " << f_c1[2] << std::endl;
//                    // std::cout << " 2 " << f_c2[0] << " " << f_c2[1] << " " << f_c2[2] << std::endl;
                    
                    nf->SetFaceID(FId);
                    face_map5[FId] = nf;
                    
                    cntdu++;
                    FId++;
                }
            }
            
            
            
            FaceSharedPtr sharedFPointer = std::shared_ptr<NekFace>(new NekFace(face2edge));
            
            
//            NekFace* sharedF = new NekFace(face2edge);
            NekFace sharedF(face2edge);
            pair<FaceSet::iterator, bool> testIns;
            testIns = m_FaceSet.insert(sharedF);
            
            pair<FaceSetPointer::iterator, bool> testInsPointer;
            testInsPointer = m_FaceSetPointer.insert(sharedFPointer);

//            if(m_FaceSet.find(sharedF)!=m_FaceSet.end())
//            {
//                facF++;
//            }
//            if(testIns.second)
//            {
//                foundd++;
//            }
            
//            NekFace* sharedF = new NekFace(face2edge);
//            NekFace sharedF(face2edge);
//            sort(face2edge_copy.begin(),face2edge_copy.end());
//
//            if(face_map3.find(face2edge_copy)==face_map3.end())
//            {
//                face_map3.insert(face2edge_copy);
//                fcid++;
//            }

            
//            FaceSharedPtr sharedF2 = std::shared_ptr<NekFace>(new NekFace(face2edge_copy));
            
//            NekFace* sharedF2 = new NekFace(face2edge_copy);


//            pair<FaceSet::iterator, bool> testIns;
//            testIns = m_FaceSet.insert(sharedF2);
//
//            int fhash2 = sharedF2->GetFaceHash();
//            if(hashset.find(fhash2)==hashset.end())
//            {
//                hashset.insert(fhash2);
//            }
//
//            if(testIns.second)
//            {
////              // std::cout << "FaceHash " << sharedF->GetFaceHash() << std::endl;
//                sharedF->SetFaceID(facetel);
//                facetel++;
//            }
//            else
//            {
//
//            }
            
            
//          else
//          {
//              // std::cout << "Duplicate " << sharedF->GetFaceID() << std::endl;
//          }
            
            
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
                
                
                
                NekFace F(face2edge);
                if(m_FaceSet.find(F)!=m_FaceSet.end())
                {
                    facF++;
                }
                
//                FaceSharedPtr f2ePointer = std::shared_ptr<NekFace>(new NekFace(f2e));
//
//                pair<FaceSetPointer::iterator, bool> testInsPointer;
//                testInsPointer = plines->FaceSetPointer_shell.insert(f2ePointer);
//
//                if(testInsPointer.second)
//                {
//                    (*testInsPointer.first)->SetFaceID(shellFid);
//
//                    plines->shellF_map[shellFid2] = ShellElem;
//
//                    shellFid2++;
//                }
                
                FaceSharedPtr f2ePointer = std::shared_ptr<NekFace>(new NekFace(faceOverts_vec));
                FaceSetPointer::iterator testInsPointer = plines->FaceSetPointer_shell.find(f2ePointer);
                
                
//                if(plines->shellFaces2Element.find(faceOverts_set)!=plines->shellFaces2Element.end())
                if(testInsPointer!=plines->FaceSetPointer_shell.end())
                {
                    int shellFid = (*testInsPointer)->GetFaceID();
                    
                    //// std::cout << "shellFid " << shellFid << " " << plines->FaceSetPointer_shell.size() << std::endl;
                    
                    std::vector<int> ElPrismConnect = plines->shellF_map[shellFid];
                    
                    std::vector<int> ElPrismConnect2 = plines->shellFaces2Element[faceOverts_set];
                    int wallface                     = plines->shellFaces2WallFace[shellFid];
                    
                    int pelem = ElPrismConnect[0];
                    
                    WallFace2PrismElem[wallface] = pelem;
                    std::vector<int> verts(3);
                    verts[0] = faceVerts[q][0];
                    verts[1] = faceVerts[q][1];
                    verts[2] = faceVerts[q][2];
                    
                    std::vector<int> edges(3);
                    edges[0] = face2edge[0];
                    edges[1] = face2edge[1];
                    edges[2] = face2edge[2];
                    
                    shellMapOrientVerts[wallface] = verts;
                    shellMapOrientEdges[wallface] = edges;
                    
                    cntshell++;
                }

                std::vector<double> cent = ComputeCentroid(faceCoords);

                gfid++;
            }
            else
            {
                int gfidn = face_set[face2edge_set];
                std::vector<double> cent = ComputeCentroid(faceCoords);
                
                

                gfaces_copy[q] = gfidn;
            }
        }

        element_map[elid] = gfaces_copy;

        if(itmesh->second->nodes.size()==4)
        {
            bool checkTet = CheckTetRotation(itmesh->second,elid);
            if(!checkTet)
            {
                // std::cout << "The orientation of tetrahedron with ID " << elid << " is not correct." << std::endl;
                BadTets++;
            }
        }

        if(itmesh->second->nodes.size()==6)
        {
            bool checkPrism = CheckPrismRotation(itmesh->second,elid);
            
            // std::cout << "We are never here " << std::endl;
            
            if(!checkPrism)
            {
                // std::cout << "The orientation of prism with ID " << elid << " is not correct." << std::endl;

                BadPrisms++;
            }
        }
        teller++;
    }
    
    //// std::cout << "m_FaceSet " << m_FaceSet.size() << " " << face_map.size() << " " << face_map2.size() << " " << facetel << " " << face_map3.size() << " " << face_map4.size() << " " << face_map5.size() << " " << dupl << " " << cntdu  << " " << hashset.size() << " found " << foundd << " " << facF << " " << m_FaceSetPointer.size() << std::endl;
    
    //// std::cout << "cntshell " << cntshell << std::endl;
    
    std::map<int,std::vector<int> >::iterator itm;
    int prmsID = TetTransferMap.size();
    int cntr=0;
    int f1cnt=0;
    int f2cnt=0;
    for(itm=plines->SinglePrismVertsPerWallFace.begin();
        itm!=plines->SinglePrismVertsPerWallFace.end();
        itm++)
    {
        int wf = itm->first;
        
        int nprismOnLine = itm->second.size();
        
        int elidStartTest = WallFace2PrismElem[wf];
        
        std::vector<int> testF(3);
        testF[0] = shellMapOrientVerts[wf][0];
        testF[1] = shellMapOrientVerts[wf][1];
        testF[2] = shellMapOrientVerts[wf][2];
        
        //for(int q=0;q<nprismOnLine;q++)
        //{
        //    int elidStart = itm->second[nprismOnLine-1-q];
            
        //    std::map<int,std::vector<int> > LineNodes = plines->ElementLinesNodes[wf];
        
//        if(LineNodes.find(elidStart)!=LineNodes.end())
//        {
            std::vector<int> Lnodes = itm->second;//LineNodes[elidStart];
            
            
            std::cout << Lnodes[0] << " " << Lnodes[1] << " " << Lnodes[2] 
            					  << " " << Lnodes[3] << " " << Lnodes[4] << " " << Lnodes[5] << std::endl;
            
            std::map<int,int> oppositeNodes;
            oppositeNodes[vert_map[old2new_v[Lnodes[0]]]] = vert_map[old2new_v[Lnodes[3]]];
            oppositeNodes[vert_map[old2new_v[Lnodes[1]]]] = vert_map[old2new_v[Lnodes[4]]];
            oppositeNodes[vert_map[old2new_v[Lnodes[2]]]] = vert_map[old2new_v[Lnodes[5]]];
            oppositeNodes[vert_map[old2new_v[Lnodes[3]]]] = vert_map[old2new_v[Lnodes[0]]];
            oppositeNodes[vert_map[old2new_v[Lnodes[4]]]] = vert_map[old2new_v[Lnodes[1]]];
            oppositeNodes[vert_map[old2new_v[Lnodes[5]]]] = vert_map[old2new_v[Lnodes[2]]];
            std::vector<int> prismNodes(6);
            
            prismNodes[0] = testF[0];
            prismNodes[1] = testF[1];
            prismNodes[2] = oppositeNodes[testF[1]];

            prismNodes[3] = oppositeNodes[testF[0]];
            prismNodes[4] = testF[2];
            prismNodes[5] = oppositeNodes[testF[2]];
        
            std::vector<std::vector<double> > coords;
            std::map<int,std::vector<double> > coordMap;
            
            for(int i=0;i<prismNodes.size();i++)
            {
                int ovid = new2old_v[vert_map_inv[prismNodes[i]]];
                
                std::vector<double> crd(3);
                crd[0] = us3d->xcn->getVal(ovid,0);
                crd[1] = us3d->xcn->getVal(ovid,1);
                crd[2] = us3d->xcn->getVal(ovid,2);
                coordMap[prismNodes[i]] = crd;
                coords.push_back(crd);
            }
            if(prmsID==55667)
            {
                // std::cout <<std::endl;
            }
            // int newPrismID = PrismTransferMap[elidStart];
                            
            NekElement* prism   = SetPrism(prismNodes,coords,prmsID,0);
            prisms[prmsID]      = prism;
            
            std::vector<int> testFCopy(3);
            testFCopy[0] = oppositeNodes[testF[0]];
            testFCopy[1] = oppositeNodes[testF[1]];
            testFCopy[2] = oppositeNodes[testF[2]];
            
            testF[0] = testFCopy[0];
            testF[1] = testFCopy[1];
            testF[2] = testFCopy[2];
            
            prmsID++;
            
//        }
        cntr++;
        //}
    }
    
   
    //// std::cout << "Fs  " << f1cnt << " " << f2cnt << std::endl;
    int pit = 0;
    for(itmesh=prisms.begin();itmesh!=prisms.end();itmesh++)
    {
        NekElement* elem = itmesh->second;
        int nid = itmesh->first;

//        if(elem->nodes.size()==4)
//        {
//            std::vector<std::vector<int> > edges = elem->m_edges;
//            std::map<int,int> Tet_l2g_edge;
//
//            for(int s=0;s<6;s++)
//            {
//                std::set<int> eset;
//                for(int k=0;k<2;k++)
//                {
//                    eset.insert(edges[s][k]);
//                }
//                if(edge_set.find(eset)==edge_set.end())
//                {
//                    edge_set[eset]  = eId;
//                    edge_map[eId]   = edges[s];
//                    Tet_l2g_edge[s] = eId;
//                    eId++;
//                }
//                else
//                {
//                    int eIdn        = edge_set[eset];
//                    Tet_l2g_edge[s] = eIdn;
//                }
//            }
//
//            el2_l2g_edge[nid] = Tet_l2g_edge;
//
//        }
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
                    
//                    if(pit==0)
//                    {
//                        // std::cout << "eId " << eId << " " <<  eId << std::endl;
//                    }
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
        
        pit++;
    }
    
    // std::cout << "prisms " << prisms.size() << std::endl;
    
    for(itmesh=prisms.begin();itmesh!=prisms.end();itmesh++)
    {
        int elid = itmesh->first;

        std::vector<std::vector<std::vector<int> > > faceEdges = itmesh->second->m_faceEdges;
        std::vector<std::vector<int> > faceVerts = itmesh->second->m_faceVertices;

        int nfaces = faceEdges.size();
        std::vector<int> gfaces_copy(nfaces);

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

//                if(geid == 69686)
//                {
//                    // std::cout << "69686  " << p << " " << nEdges << std::endl;
//
//                    for(int r=0;r<2;r++)
//                    {
//                        // std::cout << faceEdges[q][p][r] << " ";
//                    }
//                    // std::cout << std::endl;
//                }
                face2edge[p] = geid;
                face2edge_set.insert(geid);
                gEdge2lEdgeMap[geid] = p;
            }

            if(face_set.find(face2edge_set)==face_set.end())
            {
                face_set[face2edge_set] = gfid;
                face_map[gfid] = face2edge;
                gfaces_copy[q] = gfid;
                
//                if(gfid == 114639 || gfid == 114640 || gfid == 114641)
//                {
//                    int nv = face2edge.size();
//                    // std::cout << "GFID = " << gfid << " :: and elid " << elid << " >> ";
//                    for(int u=0;u<nv;u++)
//                    {
//                        // std::cout << face2edge[u] << " ";
//                    }
//
//                    // std::cout << std::endl;
//
//                }

                if(refFaceMap.find(faceOverts_set)!=refFaceMap.end())
                {
                    int faceRef = refFaceMap[faceOverts_set];
                    BoundaryComposites[faceRef].push_back(gfid);
                }
                
                
                
                
                
                std::vector<double> cent = ComputeCentroid(faceCoords);

                gfid++;
            }
            else
            {
                int gfidn = face_set[face2edge_set];
                std::vector<double> cent = ComputeCentroid(faceCoords);

                gfaces_copy[q] = gfidn;
            }
        }

        element_map[elid] = gfaces_copy;

        if(itmesh->second->nodes.size()==4)
        {
            bool checkTet = CheckTetRotation(itmesh->second,elid);
            if(!checkTet)
            {
                // std::cout << "The orientation of tetrahedron with ID " << elid << " is not correct." << std::endl;
                BadTets++;
            }
        }

        if(itmesh->second->nodes.size()==6)
        {
            bool checkPrism = CheckPrismRotation(itmesh->second,elid);

            if(!checkPrism)
            {
                // std::cout << "The orientation of prism with ID " << elid << " is not correct." << std::endl;

                BadPrisms++;
            }
        }

        teller++;
    }
    
    
    /**/
    
    
    
    // =================== Done setting up the faces ==================

    
    Array<double>* xcn_update = new Array<double>(us3d->xcn->getNrow(),3);
    
    std::map<int,int>::iterator itmp;
    
    for(itmp=vert_map.begin();itmp!=vert_map.end();itmp++)
    {
        int oldvid = new2old_v[itmp->first];
        int newvid = itmp->second;
        
        //// std::cout << oldvid << " " << newvid << " " << us3d->xcn->getVal(oldvid,0) << " " << us3d->xcn->getVal(oldvid,1) << " " << us3d->xcn->getVal(oldvid,2) << std::endl;
        
        xcn_update->setVal(newvid,0,us3d->xcn->getVal(oldvid,0));
        xcn_update->setVal(newvid,1,us3d->xcn->getVal(oldvid,1));
        xcn_update->setVal(newvid,2,us3d->xcn->getVal(oldvid,2));

    }
    
//    // std::cout << "Done checking Prisms " << badRotPrism << " " << TetWRong << " " << facesDone.size() << std::endl;
//
//    // std::cout << "number of edges " << edge_map.size() << std::endl;
//    // std::cout << "number of faces " << face_map.size() << " " << face_set.size()  << " " << us3d->ifn->getNrow() << std::endl;
//
//    ReadXmlFile(filename);
    WriteXmlFile("nektarpp.xml",xcn_update,edge_map,element_map,us3d->zdefs,face_map,us3d->ief,BoundaryComposites);
    
    /**/
    MPI_Finalize();
    
}

