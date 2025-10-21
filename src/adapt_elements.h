
#include "adapt.h"

#ifndef ADAPT_ELEMENTS_H
#define ADAPT_ELEMENTS_H


std::vector<std::vector<int> > getTetraFaceMap();

std::vector<std::vector<int> > getPrismFaceMap();

std::vector<std::vector<int> > getPyramidFaceMap();

std::vector<std::vector<int> > getHexFaceMap();

#endif