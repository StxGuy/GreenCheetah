using LinearAlgebra
using Plots
#using PyPlot

mutable struct GreenParameters
    nrows   :: Integer
    ncols   :: Integer
    dx      :: Real
    dy      :: Real
    hex     :: Real
    hey     :: Real
    Pot     :: Matrix{Float64}
    iη      :: ComplexF64
end


#-----------------------------#
# Initialize Green's function #
#-----------------------------#
function Green(Pot,height,width,effective_mass)
    ħ = 1.054571817E-34     # J.s
    q = 1.602176634E-19     # C
    mo = 9.1093837015E-31   # kg
    
    rows,cols = size(Pot)
    
    dx = width/cols
    dy = height/rows
    
    hex = ((ħ/dx^2)*(ħ/(effective_mass*mo)))/q
    hey = ((ħ/dy^2)*(ħ/(effective_mass*mo)))/q
    
    G = GreenParameters(rows,cols,dx,dy,hex,hey,Pot,1E-9*im)
    
    #println("Hopping energy x: ",1E3*hex," [meV]")
    #println("Hopping energy y: ",1E3*hey," [meV]")
    
    return G
end 

#-----------------------#
# Lead Green's function #
#-----------------------#
function leadGF(E,par::GreenParameters)
    π = 3.14159265358979323
    
    G = zeros(ComplexF64,par.nrows,par.nrows)
    U = zeros(ComplexF64,par.nrows,par.nrows)
    
    dL = 1.0/(par.nrows+1)
    Vx = -par.hex
    Vy = -par.hey
    q = 2*abs(Vx)
            
    for i in 1:par.nrows
        p = E + 2*(Vx+Vy) + 2*abs(Vy)*cos(Complex(π*i*dL)) + par.iη;
        G[i,i] = (2*p/q^2)*(1.0 - sqrt(Complex(1.0 - (q/p)^2)));
        
        for j in 1:par.nrows
            U[i,j] = sqrt(Complex(2*dL))*sin((j*pi*dL)*i);
        end
    end
        
    return U*G*U';    
end

#-------------#
# Hamiltonian #
#-------------#
function Hamiltonian(n,par::GreenParameters)
    P = par.Pot[:,n]
    v = -par.hey*ones(par.nrows - 1)
    H = Tridiagonal(v, 2*(par.hex + par.hey)*ones(par.nrows) + P, v)
    
    return H
end

#-------------------#
# Density of States #
#-------------------#
function DOS(Energy, B, par::GreenParameters)
    h = 6.62607015E-34      # J.s
    q = 1.602176634E-19     # C

    M = zeros(ComplexF64,par.nrows,par.nrows,par.ncols+1)
    D = zeros(par.nrows, par.ncols)
    
    Φ = B*par.dx*par.dy
    Φo = h/q
    τ = -I(par.nrows)*par.hex*exp(im*Φ/Φo)
    τl = τ'
    
    Gr = leadGF(Energy, par)
    M[:,:,par.ncols+1] = Gr
        
    # Backward propagation
    Gnn = Gr
    for j in par.ncols:-1:1
        Σ = τ*Gnn*τl
        H = Hamiltonian(j,par)
        Gnn = inv((Energy + par.iη)*I(par.nrows) - H - Σ)
        M[:,:,j] = Gnn
    end
   
    # Forward progpagation
    Gnn = Gr
    for j in 1:par.ncols
        # propagate
        Σ = τl*Gnn*τ
        H = Hamiltonian(j,par)
        Gnn = inv((Energy + par.iη)*I(par.nrows) - H - Σ)
        
        #merge
        Σ = τ*M[:,:,j+1]*τl
        G = (I(par.nrows) - Gnn*Σ) \ Gnn
        
        #DOS
        D[:,j] = -imag(diag(G))
    end
    
    return D
end         

#--------------#
# Transmission #
#--------------#
function Transmission(Energy, B, par::GreenParameters)
    h = 6.62607015E-34      # J.s
    q = 1.602176634E-19     # C
    
    Φ = B*par.dx*par.dy
    Φo = h/q
    τ = -I(par.nrows)*par.hex*exp(im*Φ/Φo)
    τl = τ'    
        
    # Begin with left lead
    G_lead = leadGF(Energy, par)
    G_nn = G_lead
    G_Ln = G_lead
    
    # Loop over columns
    for n in 1:par.ncols
        Σ = τl*G_nn*τ
        H = Hamiltonian(n,par)
        G_nn = inv((Energy + par.iη)*I(par.nrows) - H - Σ)
        G_Ln = G_Ln*τ*G_nn
    end
    
    # Add right lead
    Σ = τl*G_nn*τ
    G_nn = (I(par.nrows)-G_lead*Σ) \ G_lead
    G = G_Ln*τ*G_nn
    
    # Compute Γ matrix
    Σ = τl*G_lead*τ
    Γ = im*(Σ - Σ')
    T = Γ*G*Γ*G'
    t = tr(real(T))
    
    return t
end


#==============================================================#
#                           MAIN                               #
#==============================================================#
L = 100
Pot = ones(L,L)
l2 = L÷2 

for j in 1:L
    for i in 1:L
        r = ((i-l2)/l2)^2 + ((j-l2)/l2)^2
        if (r <= 1 || (i > 80 && i < 120))
            Pot[i,j] = 0
        end
    end
end

pars = Green(Pot,675E-9,675E-9,0.041)

t = []
E = []
B = 1
anim = @animate for i in 1:200
#for i in 1:200
    Energy = 50E-4 + i*1E-5

    #push!(t,Transmission(Energy,pars))
    #push!(E,Energy)
    
    D = DOS(Energy,B,pars)
    heatmap(D,c=:redsblues)
end

#plot(E,t)
#show()

gif(anim, "test.gif", fps = 5)
