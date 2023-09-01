#include "adapt.h"
#include "adapt_schedule.h"
#include "adapt_parstate.h"
#include "adapt_array.h"
#include "adapt_parmetisstate_lite.h"
#include "adapt_datastruct.h"
#include "adapt_array.h"
#include "adapt_distri_parstate.h"


#ifndef ADAPT_REPARTITION_H
#define ADAPT_REPARTITION_H

class RepartitionObject{
        public:
                RepartitionObject(){};
                RepartitionObject(mesh* meshInput,
                                std::map<int,std::vector<int> > elements,
                                std::map<int,std::vector<int> > trace,
                                std::map<int,std::vector<double> > data,
                                MPI_Comm comm);
                ~RepartitionObject();

                void GetSharedTraces(std::map<int,std::vector<int> > elements,
                                        std::map<int,std::vector<int> > ife,
                                        std::map<int,int > if_ref,
                                        std::map<int,int > iet,
                                        std::vector<int> element2rank,
                                        MPI_Comm comm);

                std::map<int,std::vector<int> > GetOptimalDistributionSchedule(std::map<int,std::vector<int> > elements,MPI_Comm comm);

                void DeterminePartitionLayout(std::map<int,std::vector<int> > elements, MPI_Comm comm);

                std::vector<int> getGlobalElement2Rank();
        private:
                std::vector<int> part;
                std::vector<int> part_global;

};

#endif


