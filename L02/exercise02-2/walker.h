#include<memory>
#include <iostream>
#include <cmath>
#include "random.h"

class Walker
{
public:
    Walker();
    Walker(Random*, double, double, double, double);

    double x;
    double y;
    double z;
    double step_len;
    Random* rnd;
    
    double R2(); //squared norm of position vector
    void set_step(double);
    void set_pos(double, double, double);
    void lattice_evol();
    void cont_evol();

};
