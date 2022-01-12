#include "green.hpp"

int main(void) {
    mat V;
    vec t;
    GreensFunction G;
    
    
    V = zeros(50,50);
    G = GreensFunction(V,1E-6,1E-6,0.041);
    
    t = G.EnergySweep(0,1E-3,10,0);
    
    cout << t << endl;
}
    
