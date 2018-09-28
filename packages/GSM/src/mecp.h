#ifndef MECP_H
#define MECP_H


#include "base.h"


class Mecp : public Base
{
 private:
    double V0;

 public:
    void driver();
    void calc_V0(); 
    Mecp(int run, int nprocs,int NNODES);/*, float DQMAG_SSM_MAX, float DQMAG_SSM_MIN,
         int MAX_OPT_ITERS, int STEP_OPT_ITERS, int INITIAL_OPT,
         float SSM_DQMAX, float CONV_TOL, float ADD_NODE_TOL,
         int BOND_FRAGMENTS, int NNODES);*/

};

#endif
