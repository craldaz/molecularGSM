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
    Mecp(int name,int nprocs);

};

#endif
