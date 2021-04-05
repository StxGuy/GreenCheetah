module GreenCheetah
    implicit none
    
    integer,parameter,private   :: dp = kind(1.d0)
    
    real(dp),parameter  :: hbar = 1.054571817E-34   ! J.s
    real(dp),parameter  :: Planck = 6.62607015E-34  ! J.s
    real(dp),parameter  :: mass = 9.1093837015E-31  ! kg
    real(dp),parameter  :: efm = 0.041
    real(dp),parameter  :: exrg = 1.602176634E-19   ! C
    real(dp),parameter  :: kB = 1.380649E-23        ! J/K
    real(dp),parameter  :: eta = 1E-15
           
    private :: Fermi
    
    type,public :: GF
        real,allocatable    :: Pot(:,:)
        integer             :: nr,nc
        real(dp)            :: hex,hey
        real(dp)            :: dx,dy
    
        contains
    
        procedure           :: DOS
        procedure           :: transmission
        procedure           :: SingleT
        procedure           :: SGM
        procedure,private   :: Hamiltonian
        procedure,private   :: leadGreensFunction
        procedure,private   :: GAsm
        
        final                :: GreenDestructor
    end type
           

    interface GF
        procedure   :: GreenConstructor
    end interface
    
    contains
    
    !------------------------------------------------------!
    !             CONSTRUCTOR AND DESTRUCTOR               !
    !------------------------------------------------------!
    ! Constructor
    function GreenConstructor(V,dx,dy) result(self)
        implicit none
        
        real(dp),intent(in) :: V(:,:)
        real(dp),intent(in) :: dx,dy
        type(GType)         :: self
        
        integer             :: nr,nc,s(2)
        
        s = shape(V)
        nr = s(1)
        nc = s(2)
        
        
        allocate(self%Pot(nr,nc))
        self%Pot = V
        self%nr = nr
        self%nc = nc
        self%hex = ((hbar/dx**2)*(hbar/(efm*mass)))/exrg
        self%hey = ((hbar/dy**2)*(hbar/(efm*mass)))/exrg
        self%dx = dx
        self%dy = dy
        
        print *,"Hopping energies in x and y [eV]:"
        print *,self%hex
        print *,self%hey
        print *,"---------------------------------"
        
    end function
    
    ! Destructor
    subroutine GreenDestructor(self)
        implicit none
        
        type(GType)         :: self
        
        deallocate(self%Pot)
    end subroutine
        
    !------------------------------------------------------!
    !                PRIVATE MODULE FUNCTIONS              !
    !------------------------------------------------------!
       
    !--------- Fermi Function ---------!
    function Fermi(Energy,T) result(f)
        implicit none
        
        ! Energy = E - mu
        
        real(dp),intent(in) :: Energy
        real(dp),intent(in) :: T
        real(dp)            :: f

        f = 1.0/(exp(exrg*Energy/(kB*T))+1)
        
        print *,exrg*Energy/(kB*T)
        
    end function
       
     
    !------- Green's function for the leads -------!
    function leadGreensFunction(self,E) result(G)
        use linalg
        use omp_lib
        
        implicit none
        
        class(GType),intent(in) :: self
        real(dp),intent(in)     :: E
        complex(dp)             :: G(self%nr,self%nc)

        real(dp),parameter      :: pi = 3.14159265358979323
        complex(dp)             :: p,q
        integer                 :: i,j,n
        real(dp)                :: Vx,Vy
        
        complex(dp)    :: U(self%nr,self%nc)

        Vx = -self%hex
        Vy = -self%hey
        
        G = 0.0
        q = 2*abs(Vx)
        !-$OMP PARALLEL DO PRIVATE(i,p) SHARED(G)
        do i = 1,self%nr
            p = E + 2*real(Vx+Vy)+2*abs(Vy)*cos(pi*cmplx(i,0)/(self%nr+1)) + cmplx(0,1E-9)
            G(i,i) = (2.0*p/(q*q))*(cmplx(1.0,0)-sqrt(1-q*q/(p*p)))
        end do
        !-$OMP END PARALLEL DO
        
        ! Transformation Unitary Matrix U
        !-$OMP PARALLEL DO PRIVATE(i,j) SHARED(U)
        do j = 1,self%nc
        do i = 1,self%nr
            U(i,j) = sqrt(cmplx(2.0,0)/(self%nr+1))*sin((j*pi/(self%nr+1))*i)
        end do
        end do
        !-$OMP END PARALLEL DO
        
        G = tric(U,G,transpose(U))        
        
    end function
    
    
    !----- Hamiltonian -----!
    function Hamiltonian(self,n) result(H)
        use omp_lib
        implicit none
        
        class(GType),intent(in) :: self
        integer,intent(in)      :: n
        real(dp)                :: H(self%nr,self%nc)
        
        integer :: s,i
        
        H = 0.0
        
        !$OMP PARALLEL PRIVATE(i) SHARED(H) 
        !$OMP DO 
        do i = 1,self%nr
            H(i,i) = 2.0*(self%hex + self%hey) + self%pot(i,n)
            if (i > 1) then
                H(i,i-1) = -self%hey
            end if
            if (i < self%nr) then
                H(i,i+1) = -self%hey
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL
    end function
        

    !------------------------------------------------------!
    !                   PUBLIC FUNCTIONS                   !
    !------------------------------------------------------!
    
    !------ Density of States -----!
    function DOS(self,Energy) result(D)
        use linalg
        implicit none
        
        real(dp),intent(in)     :: Energy
        class(GType),intent(in) :: self
        real(dp)                :: D(self%nr,self%nc)
        
        integer                 :: n,k,j
        
        real(dp)                :: I(self%nr, self%nr)
        real(dp)                :: V(self%nr, self%nr)
        complex(dp)             :: Ho(self%nr, self%nr)
        complex(dp)             :: Gnn(self%nr, self%nr)
        complex(dp)             :: Gl(self%nr, self%nr)
        complex(dp)             :: Gr(self%nr, self%nr)
        complex(dp)             :: G(self%nr, self%nr)
        complex(dp)             :: S(self%nr, self%nr)
        complex(dp)             :: P(self%nr, self%nr)
        complex(dp)             :: M(self%nr, self%nr, self%nc+1)

        
        ! Indentity and perturbation matrices
        I = 0.0
        V = 0.0
        !$OMP PARALLEL DO
        do n = 1,self%nr
            I(n,n) = 1.0
            V(n,n) = -self%hex
        end do
        !$OMP END PARALLEL DO
        
            
        ! Lead Green's function
        Gr = self%leadGreensFunction(Energy)
        M(:,:,self%nc+1) = Gr
        
        ! Sigma = V+.Gnn.V
        ! G = [(E+i.eta).I - H]^(-1)
        !(I-G.Sigma).Gnn = G
        Gnn = Gr
        do j = self%nc,1,-1
            S = d2tri(Gnn,V)
            G = GAssembly1(Energy,self%Hamiltonian(j),eta)
            P = I-matmul(G,S)
            Gnn = solve(P,G)
            M(:,:,j) = Gnn
        end do
        
        ! Sigma = V+.Gnn.V
        ! G = [(E+i.eta).I - H]^(-1)
        ! (I-G.Sigma).Gnn = G
        !
        ! Sigma = V.Gnn.V+
        ! (I-Gnn.Sigma).G = Gnn
        !
        ! D = -image(diag(G))
        Gnn = Gr
        do j = 1,self%nc
            S = d2tri(Gnn,V)
            G = GAssembly1(Energy,self%Hamiltonian(j),eta)
            P = I-matmul(G,S)
            Gnn = solve(P,G)
            
            S = d2tri(M(:,:,j+1),V)
            P = I-matmul(Gnn,S)
            G = solve(P,Gnn)
            
            do k = 1,self%nr
                D(k,j) = -imag(G(k,k))
            end do
        end do
    end function
    
    !------ Single Point Transmission Function ------!
    function SingleT(self,Energy) result(r)
        use linalg
        use omp_lib
        implicit none
        
        class(GType),intent(in) :: self
        real(dp),intent(in)     :: Energy
        real(dp)                :: r
        
        real(dp)        :: I(self%nr,self%nr)
        real(dp)        :: V(self%nr,self%nr)
        complex(dp)     :: G_lead(self%nr,self%nr)
        complex(dp)     :: G_nn(self%nr,self%nr)
        complex(dp)     :: G_Ln(self%nr,self%nr)
        complex(dp)     :: S(self%nr,self%nr)
        real(dp)        :: Ho(self%nr,self%nr)
        complex(dp)     :: G(self%nr,self%nr)
        complex(dp)     :: P(self%nr,self%nr)
        complex(dp)     :: Gama(self%nr,self%nr)
        integer         :: n,k
        
        ! Indentity and perturbation matrices
        I = 0.0
        V = 0.0
        !$OMP PARALLEL DO
        do n = 1,self%nr
            I(n,n) = 1.0
            V(n,n) = -self%hex
        end do
        !$OMP END PARALLEL DO
        
        G_lead = self%leadGreensFunction(Energy)
        G_nn = G_lead
        G_Ln = G_lead
                        
        ! Loop over columns
        ! Sigma = V.G_nn.V+
        ! G_nn = [(E+i.eta).I - H - S]^(-1)
        ! G_Ln = G_Ln.V.G_nn
        do n = 1,self%nc
            G_nn = GAssembly(Energy,self%Hamiltonian(n),d2tri(G_nn,V),eta)
            G_Ln = tri(G_Ln,V,G_nn)
        end do
        
        ! [G_lead.V.G_nn.V).G_nn = G_lead
        ! G = G_Ln.V.G_nn
        P = I - matmul(G_lead,d2tri(G_nn,V))
        G_nn = solve(P,G_lead)
        G = tri(G_Ln,V,G_nn)
        
        ! S = V+.G_lead.V
        ! Gamma = i.(S-S+)
        S = d2tri(G_lead,V)
        Gama = cmplx(0,1)*(S - conjg(S))
                                
        ! Compute trace
        ! r = Trace(Gama.G.Gama.G+)
        r = 0.0
        P = matmul(matmul(Gama,G),matmul(Gama,conjg(G)))
        do k = 1,self%nr
            r = r + P(k,k)
        end do
    end function
        
    
    !------ Transmission Function ------!
    function Transmission(self,sz,Emax) result(r)
        use linalg
        use omp_lib
        implicit none
        
        class(GType),intent(in) :: self
        real(dp),intent(in)     :: Emax
        integer,intent(in)      :: sz
        real(dp)                :: r(sz)
        
        real(dp)        :: I(self%nr,self%nr)
        real(dp)        :: V(self%nr,self%nr)
        complex(dp)     :: G_lead(self%nr,self%nr)
        complex(dp)     :: G_nn(self%nr,self%nr)
        complex(dp)     :: G_Ln(self%nr,self%nr)
        complex(dp)     :: S(self%nr,self%nr)
        real(dp)        :: Ho(self%nr,self%nr)
        complex(dp)     :: G(self%nr,self%nr)
        complex(dp)     :: P(self%nr,self%nr)
        complex(dp)     :: Gama(self%nr,self%nr)
        integer         :: m,n,k
        real(dp)        :: Energy
        
        ! Indentity and perturbation matrices
        I = 0.0
        V = 0.0
        do n = 1,self%nr
            I(n,n) = 1.0
            V(n,n) = -self%hex
        end do
        
        r = 0
        !$OMP PARALLEL DO private(Energy,G_lead,G_nn,G_Ln,P,G,S,Gama,m)
        do m = 1,sz
            Energy = Emax*real(m-1)/(sz-1)
                        
            G_lead = self%leadGreensFunction(Energy)
            G_nn = G_lead
            G_Ln = G_lead
                            
            ! Loop over columns
            ! Sigma = V.G_nn.V+
            ! G_nn = [(E+i.eta).I - H - S]^(-1)
            ! G_Ln = G_Ln.V.G_nn
            do n = 1,self%nc
                G_nn = GAssembly(Energy,self%Hamiltonian(n),d2tri(G_nn,V),eta)
                G_Ln = tri(G_Ln,V,G_nn)
            end do
            
            ! [G_lead.V.G_nn.V).G_nn = G_lead
            ! G = G_Ln.V.G_nn
            P = I - matmul(G_lead,d2tri(G_nn,V))
            G_nn = solve(P,G_lead)
            G = tri(G_Ln,V,G_nn)
            
            ! S = V+.G_lead.V
            ! Gamma = i.(S-S+)
            S = d2tri(G_lead,V)
            Gama = cmplx(0,1)*(S - conjg(S))
                                    
            ! Compute trace
            ! r = Trace(Gama.G.Gama.G+)
            P = matmul(matmul(Gama,G),matmul(Gama,conjg(G)))
            do k = 1,self%nr
                r(m) = r(m) + real(P(k,k))
            end do
            write(*,*) m,Energy,r(m)
        end do
        !$OMP END PARALLEL DO
    end function

    !------ Transmission Function ------!
    function SGM(self,Energy) result(map)
        use linalg
        use omp_lib
        implicit none
        
        class(GType),intent(in)     :: self
        real(dp),intent(in)         :: Energy
        real(dp)                    :: map(self%nr,self%nc)
        
        real(dp)        :: I(self%nr,self%nr)
        real(dp)        :: V(self%nr,self%nr)
        real(dp)        :: SPot(self%nr,self%nc)
        complex(dp)     :: G_lead(self%nr,self%nr)
        complex(dp)     :: G_nn(self%nr,self%nr)
        complex(dp)     :: G_Ln(self%nr,self%nr)
        complex(dp)     :: S(self%nr,self%nr)
        real(dp)        :: Ho(self%nr,self%nr)
        complex(dp)     :: G(self%nr,self%nr)
        complex(dp)     :: P(self%nr,self%nr)
        complex(dp)     :: Gama(self%nr,self%nr)
        integer         :: m,n,k,r,q
        real(dp)        :: tr,tro,s2=0.1
        
        ! Indentity and perturbation matrices
        I = 0.0
        V = 0.0

        do n = 1,self%nr
            I(n,n) = 1.0
            V(n,n) = -self%hex
        end do
        
        !-----------------------------------------------------------------!
        ! COMPUTE SINGLE SHOT TRANSMISSION                                !
        !-----------------------------------------------------------------!
        G_lead = self%leadGreensFunction(Energy)
        G_nn = G_lead
        G_Ln = G_lead
                        
        ! Loop over columns
        ! Sigma = V.G_nn.V+
        ! G_nn = [(E+i.eta).I - H - S]^(-1)
        ! G_Ln = G_Ln.V.G_nn
        do n = 1,self%nc
            G_nn = GAssembly(Energy,self%Hamiltonian(n),d2tri(G_nn,V),eta)
            G_Ln = tri(G_Ln,V,G_nn)
        end do
        
        ! [G_lead.V.G_nn.V).G_nn = G_lead
        ! G = G_Ln.V.G_nn
        P = I - matmul(G_lead,d2tri(G_nn,V))
        G_nn = solve(P,G_lead)
        G = tri(G_Ln,V,G_nn)
        
        ! S = V+.G_lead.V
        ! Gamma = i.(S-S+)
        S = d2tri(G_lead,V)
        Gama = cmplx(0,1)*(S - conjg(S))
                                
        ! Compute trace
        ! r = Trace(Gama.G.Gama.G+)
        P = matmul(matmul(Gama,G),matmul(Gama,conjg(G)))
        tro = 0.0
        do k = 1,self%nr
            tro = tro + real(P(k,k))
        end do
        write(*,*) "Unperturbed transmission:",tro
        !-----------------------------------------------------------------!
        
        !-----------------------------------------------------------------!
        ! SCAN GAUSSIAN POTENTIAL                                         !
        !-----------------------------------------------------------------!
        
        do m = 1,self%nc
        !write(*,*) "SGM Step:",m
        !$OMP PARALLEL DO PRIVATE(n,P,r,k,SPot,G_nn,G_Ln,G,S,Gama)
        do n = 1,self%nr
            SPot = self%Pot
!             do r = max(1,m-10),min(self%nc,m+10)
!             do k = max(1,n-10),min(self%nr,n+10)
!                 SPot(r,k) = SPot(r,k) + 1E-4*exp(-0.5*((r-m)**2+(k-n)**2)/s2)
!             end do
!             end do
            
            SPot(n,m) = SPot(n,m) + 1E-4
            
            G_nn = G_lead
            G_Ln = G_lead
                            
            ! Loop over columns
            ! Sigma = V.G_nn.V+
            ! G_nn = [(E+i.eta).I - H - S]^(-1)
            ! G_Ln = G_Ln.V.G_nn
            do q = 1,self%nc
                G_nn = self%GAsm(Energy,SPot(:,q),d2tri(G_nn,V))
                G_Ln = tri(G_Ln,V,G_nn)
            end do
            
            ! [G_lead.V.G_nn.V).G_nn = G_lead
            ! G = G_Ln.V.G_nn
            P = I - matmul(G_lead,d2tri(G_nn,V))
            G_nn = solve(P,G_lead)
            G = tri(G_Ln,V,G_nn)
            
            ! S = V+.G_lead.V
            ! Gamma = i.(S-S+)
            S = d2tri(G_lead,V)
            Gama = cmplx(0,1)*(S - conjg(S))
                                    
            ! Compute trace
            ! r = Trace(Gama.G.Gama.G+)
            P = matmul(matmul(Gama,G),matmul(Gama,conjg(G)))
            tr = 0.0
            do k = 1,self%nr
                tr = tr + real(P(k,k))
            end do
            
            map(n,m) = tr-tro
            SPot(n,m) = SPot(n,m) - 1E-4
        end do
        !$OMP END PARALLEL DO
        end do

    end function    
    
    ! Assemble a Green's function given the energy, potential and self-energy
    function GAsm(self,E,P,S) result(G)
        use linalg
        implicit none
        
        class(GType),intent(in)     :: self
        real(dp),intent(in)         :: E
        real(dp),intent(in)         :: P(:)
        complex(dp),intent(in)      :: S(:,:)
        complex(dp)                 :: G(self%nr,self%nr)
                
        integer     :: i
        real(dp)    :: H(self%nr,self%nr)
        
        
        H = 0.0
        do i = 1,self%nr
            H(i,i) = E + cmplx(0,eta) - (2.0*(self%hex + self%hey) + P(i))
            if (i > 1) then
                H(i,i-1) = self%hey
            end if
            if (i < self%nr) then
                H(i,i+1) = self%hey
            end if
        end do
        
        G = inv(H-S)
    end function
    
end module


