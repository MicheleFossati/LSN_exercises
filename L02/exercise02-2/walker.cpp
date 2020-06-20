#include "walker.h"
using namespace std;

Walker::Walker(Random* r, double m_step, double m_x = 0, double m_y = 0, double m_z = 0): x(m_x), y(m_y), z(m_z), step_len(m_step), rnd(r)
{   }

double
Walker::R2()
{
    return x*x + y*y + z*z;
}

void
Walker::set_step(double input_spacing)
{
    if(input_spacing <= 0)
        cout << "insert valid step length" << endl;
    else step_len = input_spacing;
}

void
Walker::set_pos(double m_x, double m_y, double m_z)
{
    x = m_x;
    y = m_y;
    z = m_z;
}
void
Walker::lattice_evol()
{
    double r = rnd->Rannyu();
    if(r<1./6.)
        x += step_len;
    else if (1./6 <= r && r < 2./6)
        x += -step_len;
    else if (2./6 <=r && r < 3./6)
        y += step_len;
    else if (3./6 <=r && r < 4./6)
        y += -step_len;
    else if (4./6 <=r && r < 5./6)
        z += step_len;
    else
        z += -step_len;
}

void
Walker::cont_evol()
{
    double phi = 2*M_PI*rnd->Rannyu();
    double theta = rnd->Sin();
    x += step_len*sin(theta)*cos(phi);
    y += step_len*sin(theta)*sin(phi);
    z += step_len*cos(theta);
}
