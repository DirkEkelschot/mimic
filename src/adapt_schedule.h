#include "adapt.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate.h"
#include "adapt_operations.h"

#ifndef ADAPT_SCHEDULE_H
#define ADAPT_SCHEDULE_H

struct ScheduleObj
{
    std::map<int, std::set<int> >    SendFromRank2Rank;
    std::map<int, std::set<int> >    RecvRankFromRank;
};


ScheduleObj* DoScheduling(std::map<int,std::vector<int> > Rank2RequestEntity, MPI_Comm comm, int rank, int size);

#endif
