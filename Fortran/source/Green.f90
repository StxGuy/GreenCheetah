

module GreenCheetah
    implicit none
    
    integer,parameter,private   :: dp = kind(1.d0)
    
    ! Effective mass: Semicond. Sci. Technol. 27 (2012) 035021
    ! InGaAs = 0.041
    
    real(dp),parameter  :: hbar = 1.054571817E-34   ! J.s
    real(dp),parameter  :: Planck = 6.62607015E-34  ! J.s
    real(dp),parameter  :: mass = 9.1093837015E-31  ! kg
    real(dp),parameter  :: exrg = 1.602176634E-19   ! C
    real(dp),parameter  :: kB = 1.380649E-23        ! J/K
    real(dp),parameter  :: eta = 1E-15
    real(dp),parameter  :: pi = 3.14159265358979323
               
    type,public :: Green
        real,allocatable    :: Pot(:,:)
        integer             :: nr,nc
        real(dp)            :: hex,hey
        real(dp)            :: dx,dy
        real                :: efm
        
        real(dp)            :: FermiEnergy
        real(dp)            :: FermiWaveLength
        real(dp)            :: FermiWaveVector
    
        contains
    
        procedure           :: DOS
        procedure           :: SGM
        procedure           :: EnergySweep
        procedure           :: MagneticSweep
        procedure           :: transmission
        procedure           :: getFermi
        procedure,private   :: Hamiltonian
        procedure,private   :: leadGreensFunction
        procedure,private   :: tfun
        
        final                :: GreenDestructor
    end type
           

    interface Green
        procedure   :: GreenConstructor
    end interface
    
    contains
    
    !------------------------------------------------------!
    !             CONSTRUCTOR AND DESTRUCTOR               !
    !------------------------------------------------------!
    function GreenConstructor(V,Height,Width,efm,flag) result(self)
        implicit none
        
        real(dp),intent(in) :: V(:,:)
        real,intent(in)     :: Height, Width
        real,intent(in)     :: efm              ! Effective mass
        logical,intent(in)  :: flag
        type(Green)         :: self
        
        integer             :: nr,nc,s(2)
        
        s = shape(V)
        nr = s(1)
        nc = s(2)
        
        allocate(self%Pot(nr,nc))
        self%Pot = V
        self%nr = nr
        self%nc = nc
        self%dx = Height/nr
        self%dy = Width/nc
        self%efm = efm
        self%hex = ((hbar/self%dx**2)*(hbar/(efm*mass)))/exrg
        self%hey = ((hbar/self%dy**2)*(hbar/(efm*mass)))/exrg
            
            
        if (flag) then
            write(*,*) "-----------------------------------"
            write(*,*) "Lattice constants in y and x [nm]:"
            write(*,*) self%dy/1E-9
            write(*,*) self%dx/1E-9
            write(*,*) ""
            write(*,*) "Hopping energies in y and x [meV]:"
            write(*,*) self%hey/1E-3
            write(*,*) self%hex/1E-3
            write(*,*) "-----------------------------------"
        end if
            
    end function
    
    subroutine GreenDestructor(self)
        implicit none
        
        type(Green)         :: self
        
        deallocate(self%Pot)
    end subroutine
        
    !------------------------------------------------------!
    !                PRIVATE MODULE FUNCTIONS              !
    !------------------------------------------------------!
           
    !----- Lead Green's function -----!
    function leadGreensFunction(self,E) result(G)
        use linalg
        use omp_lib
        
        implicit none
        
        class(Green),intent(in) :: self
        real,intent(in)         :: E
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
        
        G = tri(U,G,transpose(U))        
        
    end function
    
    
    !----- Hamiltonian -----!
    function Hamiltonian(self,n) result(H)
        use omp_lib
        implicit none
        
        class(Green),intent(in) :: self
        integer,intent(in)      :: n
        complex(dp)             :: H(self%nr,self%nc)
        
        integer     :: s,i

        
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
    function DOS(self,Energy,B,flag) result(D)
        use linalg
        implicit none
        
        class(Green),intent(in) :: self
        real,intent(in)         :: Energy
        real,intent(in)         :: B
        logical,intent(in)      :: flag
        real(dp)                :: D(self%nr,self%nc)
        
        integer                 :: n,k,j
        
        real(dp)                :: I(self%nr, self%nr)
        complex(dp)             :: V(self%nr, self%nr)
        complex(dp)             :: Ho(self%nr, self%nr)
        complex(dp)             :: Gnn(self%nr, self%nr)
        complex(dp)             :: Gl(self%nr, self%nr)
        complex(dp)             :: Gr(self%nr, self%nr)
        complex(dp)             :: G(self%nr, self%nr)
        complex(dp)             :: Sigma(self%nr, self%nr)
        complex(dp)             :: P(self%nr, self%nr)
        complex(dp)             :: M(self%nr, self%nr, self%nc+1)
        
        complex(dp)             :: VB

        if (flag .eqv. .True.) then
            write(*,*) "-----------------------------------"
            write(*,*) "Local density of states"
            write(*,*) "  Energy [meV]......:",Energy/1E-3
            write(*,*) "  Magnetic field [T]:",B
            write(*,*) "-----------------------------------"
        end if
        
        ! Indentity and perturbation matrices
        ! Coupling matrix
        V = 0.0
        I = 0.0
        VB = cmplx(0,exrg*B*self%dx*self%dy/hbar)
        do n = 1,self%nr
            V(n,n) = -self%hex*exp(n*VB)
            I(n,n) = 1.0
        end do
        
            
        ! Lead Green's function
        Gr = self%leadGreensFunction(Energy)
        M(:,:,self%nc+1) = Gr
        
        ! Backward Propagation
        Gnn = Gr
        do j = self%nc,1,-1
            Sigma = tri(V,Gnn,conjg(V))
            Gnn = GAssembly(Energy,self%Hamiltonian(j),Sigma,eta)
            M(:,:,j) = Gnn
        end do
        
        ! Forward + Merging
        Gnn = Gr
        do j = 1,self%nc
            ! Forward
            Sigma = tri(conjg(V),Gnn,V)
            Gnn = GAssembly(Energy,self%Hamiltonian(j),Sigma,eta)
            
            ! Merging
            Sigma = tri(V,M(:,:,j+1),conjg(V))
            P = I - matmul(Gnn,Sigma)
            G = solve(P,Gnn)
            
            ! DOS
            do k = 1,self%nr
                D(k,j) = -imag(G(k,k))
            end do
        end do
    end function
      
    
    !------- Energy Sweep ------!
    subroutine EnergySweep(self,Emin,Emax,B,T)
        use omp_lib
        implicit none
        
        class(Green),intent(in) :: self
        real,intent(in)         :: Emin, Emax
        real,intent(in)         :: B
        real,intent(inout)      :: T(:)
        
        integer     :: i, sz
        real        :: Energy
        
        sz = size(T)
        write(*,*) "* Energy sweep [meV] from",Emin/1E-3,"to",Emax/1E-3
        write(*,*) "  Steps:",sz
        write(*,*) "  Magnetic field [T]:",B
        
        !$OMP PARALLEL DO PRIVATE(i,Energy) SHARED(T)
        do i = 1,sz
            Energy = Emin + (Emax-Emin)*real(i-1)/(sz-1)
        
            T(i) = self%Transmission(Energy,B)
        end do
        !$OMP END PARALLEL DO
    end subroutine
    
    !------ Magnetic Sweep -------!
    subroutine MagneticSweep(self,Bmin,Bmax,Energy,T)
        use omp_lib
        implicit none
        
        class(Green),intent(in) :: self
        real,intent(in)         :: Bmin, Bmax
        real,intent(in)         :: Energy
        real,intent(inout)      :: T(:)
        
        integer     :: i,sz
        real        :: B
        
        sz = size(T)
        
        write(*,*) "* Magnetic sweep [T] from",Bmin,"to",Bmax
        write(*,*) "  Steps:",sz
        write(*,*) "  Energy [meV]:",Energy/1E-3
        
        !$OMP PARALLEL DO PRIVATE(i,B) SHARED(T)
        do i = 1,sz
            B = Bmin + (Bmax-Bmin)*real(i-1)/(sz-1)
            
            T(i) = self%Transmission(Energy,B)
        end do
        !$OMP END PARALLEL DO
    end subroutine
        
    !------ Transmission Function ------!
    function Transmission(self,Energy,B) result(T)
        use linalg
        implicit none
        
        class(Green),intent(in) :: self
        real,intent(in)         :: Energy
        real,intent(in)         :: B
        real                    :: T
        
        real                    :: I(self%nr,self%nr)
        complex(dp)             :: V(self%nr,self%nr)
        complex(dp)             :: G_lead(self%nr,self%nr)
        complex(dp)             :: G_nn(self%nr,self%nr)
        complex(dp)             :: G_Ln(self%nr,self%nr)
        complex(dp)             :: S(self%nr,self%nr)
        complex(dp)             :: Sigma(self%nr,self%nr)
        complex(dp)             :: G(self%nr,self%nr)
        complex(dp)             :: P(self%nr,self%nr)
        complex(dp)             :: Gama(self%nr,self%nr)
        integer                 :: m,n,k,sz
                
        complex(dp)             :: VB
        
                
        ! Coupling matrix
        V = 0.0
        I = 0.0
        VB = cmplx(0,exrg*B*self%dx*self%dy/hbar)
        do n = 1,self%nr
            V(n,n) = -self%hex*exp(n*VB)
            I(n,n) = 1.0
        end do
        
        !** Begin with left lead
        G_lead = self%leadGreensFunction(Energy)
        G_nn = G_lead
        G_Ln = G_lead
                        
        ! Loop over columns
        ! Sigma = V+.G_nn.V
        ! G_nn = [(E+i.eta).I - H - S]^(-1)
        ! G_Ln = G_Ln.V.G_nn
        do n = 1,self%nc
            Sigma = tri(conjg(V),G_nn,V)
            G_nn = GAssembly(Energy,self%Hamiltonian(n),Sigma,eta)
            G_Ln = tri(G_Ln,V,G_nn)
        end do
        
        !** Add right lead
        ! [G_lead.V.G_nn.V).G_nn = G_lead
        ! G = G_Ln.V.G_nn
        Sigma = tri(conjg(V),G_nn,V)
        P = I - matmul(G_lead,Sigma)
        G_nn = solve(P,G_lead)
        G = tri(G_Ln,V,G_nn)
        
        !** Compute Gamma matrix
        ! S = V+.G_lead.V
        ! Gamma = i.(S-S+)
        Sigma = tri(conjg(V),G_lead,V)
        Gama = cmplx(0,1)*(Sigma - conjg(Sigma))
                                
        !** Compute trace
        ! r = Trace(Gama.G.Gama.G+)
        P = matmul(matmul(Gama,G),matmul(Gama,conjg(G)))
        T = 0
        do k = 1,self%nr
            T = T + real(P(k,k))
        end do
    end function
        
    !------ Transmission Function Fast ------!
    function tfun(self,Energy,B,G_lead,V,Id,Pot) result(T)
        use linalg
        implicit none
        
        class(Green),intent(in) :: self
        real,intent(in)         :: Energy
        real,intent(in)         :: B
        complex(dp),intent(in)  :: G_lead(self%nr,self%nr)
        complex(dp),intent(in)  :: V(self%nr,self%nr)
        real,intent(in)         :: Id(self%nr,self%nr)
        real,intent(in)         :: Pot(self%nr,self%nc)
        real                    :: T        
        
        complex(dp)             :: H(self%nr,self%nr)
        complex(dp)             :: G_nn(self%nr,self%nr)
        complex(dp)             :: G_Ln(self%nr,self%nr)
        complex(dp)             :: S(self%nr,self%nr)
        complex(dp)             :: Sigma(self%nr,self%nr)
        complex(dp)             :: G(self%nr,self%nr)
        complex(dp)             :: P(self%nr,self%nr)
        complex(dp)             :: Gama(self%nr,self%nr)
        integer                 :: m,n,k,sz,i
                
        complex(dp)             :: VB
                        

        !** Build fundamental Hamiltonian
        H = 0.0
        do i = 1,self%nr
            H(i,i) = 2.0*(self%hex + self%hey)
            
            if (i > 1) then
                H(i,i-1) = -self%hey
            end if
            if (i < self%nr) then
                H(i,i+1) = -self%hey
            end if
        end do        
                        
        !** Begin with left lead
        G_nn = G_lead
        G_Ln = G_lead
                        
        ! Loop over columns
        ! Sigma = V+.G_nn.V
        ! G_nn = [(E+i.eta).I - H - S]^(-1)
        ! G_Ln = G_Ln.V.G_nn
        do n = 1,self%nc
            ! Build Hamiltonian
            do i = 1,self%nr
                H(i,i) = H(i,i) + Pot(i,n)
            end do
        
            Sigma = tri(conjg(V),G_nn,V)
            G_nn = GAssembly(Energy,H,Sigma,eta)
            G_Ln = tri(G_Ln,V,G_nn)
            
            ! Recover Hamiltonian
            do i = 1,self%nr
                H(i,i) = H(i,i) - Pot(i,n)
            end do
        end do
        
        !** Add right lead
        ! [G_lead.V.G_nn.V).G_nn = G_lead
        ! G = G_Ln.V.G_nn
        Sigma = tri(conjg(V),G_nn,V)
        P = Id - matmul(G_lead,Sigma)
        G_nn = solve(P,G_lead)
        G = tri(G_Ln,V,G_nn)
        
        !** Compute Gamma matrix
        ! S = V+.G_lead.V
        ! Gamma = i.(S-S+)
        Sigma = tri(conjg(V),G_lead,V)
        Gama = cmplx(0,1)*(Sigma - conjg(Sigma))
                                
        !** Compute trace
        ! r = Trace(Gama.G.Gama.G+)
        P = matmul(matmul(Gama,G),matmul(Gama,conjg(G)))
        T = 0
        do k = 1,self%nr
            T = T + real(P(k,k))
        end do
    end function    
    
    !--- Get Fermi energy, wavelength, and wavevector ---!
    subroutine getFermi(self,n2D,flag)
        implicit none
        
        class(Green),intent(inout)  :: self
        real,intent(in)             :: n2D  ! [cm^-2]
        logical,intent(in)          :: flag
        
        self%FermiWaveLength = sqrt(2E-4*pi/n2D)
        self%FermiWaveVector = 2*pi/self%FermiWaveLength
        self%FermiEnergy = ((hbar*self%FermiWaveVector)**2/(2*mass*self%efm))
        
        if (flag) then
            write(*,*) "Fermi Energy [meV] = ",self%FermiEnergy/(1E-3*exrg)
            write(*,*) "Fermi Wavelength [nm] = ",self%FermiWaveLength/1E-9
            write(*,*) "Fermi Wavevector [1/nm] = ",1E-9*self%FermiWavevector
        end if
    end subroutine        
        
    !--- Scanning Gate Microscopy (SGM) ---!
    function SGM(self,Energy,B,flag) result(D)
        use linalg
        use omp_lib
        implicit none
        
        class(Green),intent(in) :: self
        real,intent(in)         :: Energy
        real,intent(in)         :: B
        logical,intent(in)      :: flag
        real(dp)                :: D(self%nr,self%nc)
        
        integer                 :: n,k,j,ri,rj,tip_pos
        
        real(dp)                :: I(self%nr, self%nr)
        complex(dp)             :: V(self%nr, self%nr)
        complex(dp)             :: Ho(self%nr, self%nr)
        complex(dp)             :: Gnn(self%nr, self%nr)
        complex(dp)             :: Gnnn(self%nr, self%nr)
        complex(dp)             :: GLn(self%nr, self%nr)
        complex(dp)             :: GnR(self%nr, self%nr)
        complex(dp)             :: GLnn(self%nr, self%nr)
        complex(dp)             :: Gl(self%nr, self%nr)
        complex(dp)             :: Gr(self%nr, self%nr)
        complex(dp)             :: G(self%nr, self%nr)
        complex(dp)             :: Sigma(self%nr, self%nr)
        complex(dp)             :: Sigma1(self%nr, self%nr)
        complex(dp)             :: Sigma2(self%nr, self%nr)
        complex(dp)             :: P(self%nr, self%nr)
        complex(dp)             :: M(self%nr, self%nr, self%nc+1)
        complex(dp)             :: Q(self%nr, self%nr, self%nc+1)
        complex(dp)             :: Gama(self%nr,self%nr)
        
        complex(dp)             :: VB

        if (flag .eqv. .True.) then
            write(*,*) "-----------------------------------"
            write(*,*) "Scanning Gate Microscopy"
            write(*,*) "  Energy [meV]......:",Energy/1E-3
            write(*,*) "  Magnetic field [T]:",B
            write(*,*) "-----------------------------------"
        end if
        
        ! Indentity and perturbation matrices
        ! Coupling matrix
        V = 0.0
        I = 0.0
        VB = cmplx(0,exrg*B*self%dx*self%dy/hbar)
        do n = 1,self%nr
            V(n,n) = -self%hex*exp(n*VB)
            I(n,n) = 1.0
        end do
            
        ! Lead Green's function
        Gr = self%leadGreensFunction(Energy)
        M(:,:,self%nc+1) = Gr
        Q(:,:,self%nc+1) = Gr
        
        ! Backward Propagation
        Gnn = Gr
        GnR = Gr
        do j = self%nc,1,-1
            Sigma = tri(V,Gnn,conjg(V))
            Gnn = GAssembly(Energy,self%Hamiltonian(j),Sigma,eta)
            GnR = tri(Gnn,V,Gnr)
            M(:,:,j) = Gnn
            Q(:,:,j) = Gnr
        end do
        
        ! Forward + Merging
        Gnn = Gr
        GLn = Gr
        do j = 1,self%nc
            write(*,*) j
            ! Merge
            Sigma1 = tri(conjg(V),Gnn,V)
            Sigma2 = tri(V,M(:,:,j+1),conjg(V))
            Gnnn = GAssembly(Energy,self%Hamiltonian(j),Sigma1,Sigma2,eta)
            GLnn = tri(GLn,V,Gnnn)
            GnR = tri(Gnnn,V,Q(:,:,j+1))
        
            ! Forward
            Gnn = GAssembly(Energy,self%Hamiltonian(j),Sigma1,eta)
            GLn = tri(GLn,V,Gnn)            
            
            ! Transmission
            do tip_pos = 1,self%nr
                ! Green's function of Tip 
                do rj = 1,self%nr
                    do ri = 1,self%nr
                        G(ri,rj) = GLnn(ri,tip_pos)*GnR(tip_pos,rj)
                    end do
                end do
                !** Compute Gamma matrix
                ! S = V+.G_lead.V
                ! Gamma = i.(S-S+)
                Sigma = tri(conjg(V),Gr,V)
                Gama = cmplx(0,1)*(Sigma - conjg(Sigma))
                                        
                !** Compute trace
                ! r = Trace(Gama.G.Gama.G+)
                P = matmul(matmul(Gama,G),matmul(Gama,conjg(G)))
                D(tip_pos,j) = 0
                do k = 1,self%nr
                    D(tip_pos,j) = D(tip_pos,j) + real(P(k,k))
                end do
            end do
        end do
    end function    
end module
