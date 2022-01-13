#include "green.hpp"

// Constructor
GreensFunction::GreensFunction(mat V, float height, float width, float efm) {
    nc = V.n_cols;
    nr = V.n_rows;
    
    Potential = V;
    
    dx = width/nc;
    dy = height/nr;
    
    hex = ((hbar/pow(dx,2))*(hbar/(efm*mass)))/exrg;
    hey = ((hbar/pow(dy,2))*(hbar/(efm*mass)))/exrg;
}

GreensFunction::GreensFunction(void) {
    nc = 0;
    nr = 0;
}


// 2D Lead Green's function
cx_fmat GreensFunction::leadGreensFunction(float E) {
    float                   Vx, Vy, dL;
    complex<float>          p,q;
    const complex<float>    img = complex<float>(0,1);
    const float             pi = 3.14159265358979323;
    
    cx_fmat G, U;
    
    G = zeros<cx_fmat>(nr,nc);
    U = zeros<cx_fmat>(nr,nc);
    
    dL = 1.0/(Potential.n_rows+1);
    Vx = -hex;
    Vy = -hey;
        
    q = 2*abs(Vx);
    for(int i = 1; i <= nr; i ++) {
        p = E + complex<float>(2,0)*(Vx+Vy) + complex<float>(2,0)*abs(Vy)*cos(pi*i*dL) + eta*img;
        G(i-1,i-1) = (complex<float>(2,0)*p/pow(q,2))*(complex<float>(1,0)-sqrt(complex<float>(1,0)-q*q/(p*p)));
    }
    
    for(int j = 0; j < nc; j ++) 
        for(int i = 0; i < nr; i ++) 
            U(i,j) = sqrt(2.0*dL)*sin((j*pi*dL)*i);
    
    return U*G*U.t();
}

// Hamiltonian
cx_fmat GreensFunction::Hamiltonian(int n) {
    cx_fmat H;
    
    H = zeros<cx_fmat>(Potential.n_rows, Potential.n_rows);
    H.diag(1) = -hey*ones<cx_fmat>(nr-1);
    H.diag(-1) = -hey*ones<cx_fmat>(nr-1);
    H.diag(0) = 2*(hex+hey) + conv_to<cx_fmat>::from(Potential.col(n));
    
    return H;
}

// General transmission function
float GreensFunction::Transmission(float Energy, float B) {
    cx_fmat V, G_lead, G_nn, G_Ln, Sigma, P, Gama, G;
    
    cx_fmat I = eye<cx_fmat>(nr,nr);
    complex<float> img = complex<float>(0,1);
    complex<float> VB;
    
    // Coupling matrix
    VB = img*exrg*B*dx*dy/hbar;
    V = eye<cx_fmat>(nr,nr);
    for(int n = 0; n < nr; n ++)
        V(n,n) = -hex*exp(complex<float>(n+1,0)*VB);
    
    // Begin with left lead
    G_lead = leadGreensFunction(Energy);
    G_nn = G_lead;
    G_Ln = G_lead;
    
    /* Loop over columns
     * Sigma = V+.G_nn.V
     * G_nn = [(E+i.eta).I - H - S]^(-1)
     * G_Ln = G_Ln.V.G_nn */
    for(int n = 0; n < nc; n ++) {
        Sigma = conj(V)*G_nn*V;
        P = (Energy+img*eta)*I-Hamiltonian(n) - Sigma;
        G_nn = inv(P);
        G_Ln = G_Ln*V*G_nn;
    }
    
    /* Add right lead
     * [G_lead.V.G_nn.V).G_nn = G_lead
     * G = G_Ln.V.G_nn */
    Sigma = conj(V)*G_nn*V;
    P = I - G_lead*Sigma;
    G_nn = solve(P,G_lead);
    G = G_Ln*V*G_nn;
    
    /* Compute Gamma matrix
     * S = V+.G_lead.V
     * Gamma = i.(S-S+) */
    Sigma = conj(V)*G_lead*V;
    Gama = img*(Sigma - conj(Sigma));
    
    // r = trace(Gama.G.Gama.G+)
    P = Gama*G*Gama*conj(G);
    return trace(real(P));
}

// Energy Sweep
vec GreensFunction::EnergySweep(float Emin, float Emax, int N,float B) {
    vec     T(N);
    float   Energy;
    
    for(int i = 0; i < N; i ++) {
        Energy = Emin + (Emax - Emin)*i/N;
        T(i) = Transmission(Energy,B);
    }
    
    return T;
}
        
// Density of States
mat GreensFunction::DOS(float Energy, float B) {
    cx_fmat     V, Gr, Gnn, Sigma, Gama, G, P;
    cx_fcube    M;
    mat         D;
    
    cx_fmat I = eye<cx_fmat>(nr,nr);
    complex<float> img = complex<float>(0,1);
    complex<float> VB;    
    
    // Initializations
    D = zeros(nr,nc);
    M = zeros<cx_fcube>(nr,nr,nc+1);
    
    // Coupling matrix
    VB = img*exrg*B*dx*dy/hbar;
    V = eye<cx_fmat>(nr,nr);
    for(int n = 0; n < nr; n ++)
        V(n,n) = -hex*exp(complex<float>(n+1,0)*VB);    
    
    // Lead Green's function
    Gr = leadGreensFunction(Energy);
    M.slice(nc) = Gr;
    
    // Backward Propagation
    Gnn = Gr;
    for(int j = nc-1; j >= 0; j --) {
        Sigma = V*Gnn*conj(V);
        P = (Energy+img*eta)*I - Hamiltonian(j) - Sigma;
        Gnn = inv(P);
        M.slice(j) = Gnn;
    }
    
    // Forward + Merging
    for(int j = 0; j < nc; j ++) {
        // Forward
        Sigma = conj(V)*Gnn*V;
        P = (Energy+img*eta)*I - Hamiltonian(j) - Sigma;
        Gnn = inv(P);
        
        // Merging
        Sigma = V*M.slice(j+1)*conj(V);
        P = I - Gnn*Sigma;
        G = solve(P,Gnn);
        
        // DOS
        D.col(j) = -conv_to<mat>::from(diagvec(imag(G)));
    }
    
    return D;
}
    
    
