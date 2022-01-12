#ifndef GREENCHEETAH
#define GREENCHEETAH

#include <armadillo>

using namespace std;
using namespace arma;

const float hbar = 1.054571817E-34;  // J.s
const float Planck = 6.62607015E-34; // J.s
const float mass = 9.1093837015E-31; // kg
const float exrg = 1.602176634E-19;  // C
const float eta = 1E-15;

class GreensFunction {
    private:
        float   hex, hey;
        float   dx, dy;
        int     nr, nc;
        mat     Potential;
        
        cx_fmat leadGreensFunction(float);
        cx_fmat Hamiltonian(int);
                
    public:
        float   Transmission(float, float);
        vec     EnergySweep(float, float, int,float);
        mat     DOS(float, float);
        
        GreensFunction(mat, float, float, float);
        GreensFunction(void);
};
        
#endif        
