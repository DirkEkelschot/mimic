#include "adapt_elements.h"


std::vector<std::vector<int> > getTetraFaceMap()
{
    std::vector<std::vector<int> > tetra_faces(4);
    std::vector<int>  tetra_f0(3,0);
    tetra_f0[0] = 0; tetra_f0[1] = 1; tetra_f0[2] = 2;
    tetra_faces[0]  =  tetra_f0;
    std::vector<int>  tetra_f1(3,0);
    tetra_f1[0] = 0; tetra_f1[1] = 1; tetra_f1[2] = 3;
    tetra_faces[1]  =  tetra_f1;
    std::vector<int>  tetra_f2(3,0);
    tetra_f2[0] = 1; tetra_f2[1] = 2; tetra_f2[2] = 3;
    tetra_faces[2]  =  tetra_f2;
    std::vector<int>  tetra_f3(3,0);
    tetra_f3[0] = 0; tetra_f3[1] = 3; tetra_f3[2] = 2;
    tetra_faces[3]  =  tetra_f3;

    return tetra_faces;
}   
    

std::vector<std::vector<int> > getPrismFaceMap()
{  
    std::vector<std::vector<int> > prism_faces(5);
    std::vector<int> prism_f0(3,0);
    prism_f0[0] = 0;prism_f0[1] = 1;prism_f0[2] = 2;
    prism_faces[0]  = prism_f0;
    std::vector<int> prism_f1(3,0);
    prism_f1[0] = 3;prism_f1[1] = 4;prism_f1[2] = 5;
    prism_faces[1]  = prism_f1;
    std::vector<int> prism_f2(4,0);
    prism_f2[0] = 0;prism_f2[1] = 1;prism_f2[2] = 4;prism_f2[3] = 3;
    prism_faces[2]  = prism_f2;
    std::vector<int> prism_f3(4,0);
    prism_f3[0] = 1;prism_f3[1] = 2;prism_f3[2] = 5;prism_f3[3] = 4;
    prism_faces[3]  = prism_f3;
    std::vector<int> prism_f4(4,0);
    prism_f4[0] = 0;prism_f4[1] = 3;prism_f4[2] = 5;prism_f4[3] = 2;
    prism_faces[4]  = prism_f4;

    return prism_faces;

}

std::vector<std::vector<int> > getPyramidFaceMap()
{   
    std::vector<std::vector<int> > pyramid_faces(5);
    std::vector<int> pyramid_f0(4,0);
    pyramid_f0[0] = 0;pyramid_f0[1] = 1;pyramid_f0[2] = 2;pyramid_f0[3] = 3;
    pyramid_faces[0]  = pyramid_f0;
    std::vector<int> pyramid_f1(3,0);
    pyramid_f1[0] = 0;pyramid_f1[1] = 1;pyramid_f1[2] = 4;
    pyramid_faces[1]  = pyramid_f1;
    std::vector<int> pyramid_f2(3,0);
    pyramid_f2[0] = 1;pyramid_f2[1] = 2;pyramid_f2[2] = 4;
    pyramid_faces[2]  = pyramid_f2;
    std::vector<int> pyramid_f3(3,0);
    pyramid_f3[0] = 2;pyramid_f3[1] = 3;pyramid_f3[2] = 4;
    pyramid_faces[3]  = pyramid_f3;
    std::vector<int> pyramid_f4(3,0);
    pyramid_f4[0] = 0;pyramid_f4[1] = 4;pyramid_f4[2] = 3;
    pyramid_faces[4]  = pyramid_f4;

    return pyramid_faces;
}



std::vector<std::vector<int> > getPyramidFaceMap_vtk()
{   
    std::vector<std::vector<int> > pyramid_faces(5);
    std::vector<int> pyramid_f0(4,0);
    pyramid_f0[0] = 0;pyramid_f0[1] = 1;pyramid_f0[2] = 2;pyramid_f0[3] = 3;
    pyramid_faces[0]  = pyramid_f0;
    std::vector<int> pyramid_f1(3,0);
    pyramid_f1[0] = 0;pyramid_f1[1] = 1;pyramid_f1[2] = 4;
    pyramid_faces[1]  = pyramid_f1;
    std::vector<int> pyramid_f2(3,0);
    pyramid_f2[0] = 1;pyramid_f2[1] = 2;pyramid_f2[2] = 4;
    pyramid_faces[2]  = pyramid_f2;
    std::vector<int> pyramid_f3(3,0);
    pyramid_f3[0] = 2;pyramid_f3[1] = 3;pyramid_f3[2] = 4;
    pyramid_faces[3]  = pyramid_f3;
    std::vector<int> pyramid_f4(3,0);
    pyramid_f4[0] = 0;pyramid_f4[1] = 4;pyramid_f4[2] = 3;
    pyramid_faces[4]  = pyramid_f4;

    return pyramid_faces;
}

std::vector<std::vector<int> > getHexFaceMap()
{  
    std::vector<std::vector<int> > hex_faces(6);
    std::vector<int> f0(4,0);
    f0[0] = 0;f0[1] = 1;f0[2] = 2;f0[3] = 3;
    hex_faces[0]  = f0;
    std::vector<int> f1(4,0);
    f1[0] = 4;f1[1] = 5;f1[2] = 6;f1[3] = 7;
    hex_faces[1]  = f1;
    std::vector<int> f2(4,0);
    f2[0] = 2;f2[1] = 3;f2[2] = 7;f2[3] = 6;
    hex_faces[2]  = f2;
    std::vector<int> f3(4,0);
    f3[0] = 0;f3[1] = 1;f3[2] = 5;f3[3] = 4;
    hex_faces[3]  = f3;
    std::vector<int> f4(4,0);
    f4[0] = 0;f4[1] = 4;f4[2] = 7;f4[3] = 3;
    hex_faces[4]  = f4;
    std::vector<int> f5(4,0);
    f5[0] = 1;f5[1] = 2;f5[2] = 6;f5[3] = 5;
    hex_faces[5]  = f5;

    return hex_faces;
}



std::vector<std::vector<int> > getHexFaceMap_vtk()
{  
    std::vector<std::vector<int> > hex_faces(6);
    std::vector<int> f0(4,0);
    f0[0] = 0;f0[1] = 4;f0[2] = 7;f0[3] = 3;
    hex_faces[0]  = f0;
    std::vector<int> f1(4,0);
    f1[0] = 1;f1[1] = 2;f1[2] = 6;f1[3] = 5;
    hex_faces[1]  = f1;
    std::vector<int> f2(4,0);
    f2[0] = 0;f2[1] = 1;f2[2] = 5;f2[3] = 4;
    hex_faces[2]  = f2;
    std::vector<int> f3(4,0);
    f3[0] = 3;f3[1] = 7;f3[2] = 6;f3[3] = 2;
    hex_faces[3]  = f3;
    std::vector<int> f4(4,0);
    f4[0] = 0;f4[1] = 3;f4[2] = 2;f4[3] = 1;
    hex_faces[4]  = f4;
    std::vector<int> f5(4,0);
    f5[0] = 4;f5[1] = 5;f5[2] = 6;f5[3] = 7;
    hex_faces[5]  = f5;

    return hex_faces;
}