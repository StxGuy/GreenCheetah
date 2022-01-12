module linalg
    implicit none
    
    integer,parameter   :: dp = kind(1.d0)    
    
    interface GAssembly
        procedure   GAssembly1
        procedure   GAssembly2
        procedure   GAssembly3
    end interface
    
    contains

    !-------- Matrix inversion --------!
    function inv(A) result(Ainv)
        implicit none
        
        complex(dp),intent(in)       :: A(:,:)
        complex(dp)                  :: Ainv(size(A,1),size(A,2))
        
        integer                      :: s(2), nrows, ncols
        integer                      :: lda, lwork, info
        
        integer                      :: ipiv(size(A,1))
        complex(dp)                  :: work(size(A,1))
        
        external zgetrf
        external zgetri
        
        s = shape(A)
        nrows = s(1)
        ncols = s(2)
        lda = max(1,nrows)
        lwork = max(1,nrows)
                
        
        Ainv = A
        call zgetrf(nrows,ncols,Ainv,lda,ipiv,info)
        
        if (info .ne. 0) then
            print *, "Error in the LU factorization!"
        end if
        
        call zgetri(nrows,Ainv,lda,ipiv,work,lwork,info)
        
        if (info .ne. 0) then
            print *,"Error in the matrix inversion!"
        end if

    end function

    !------- Solve Linear System -------!
    function solve(A,B) result(X)
        implicit none
        
        complex(dp),intent(in)   :: A(:,:), B(:,:)
        complex(dp)              :: X(size(B,1),size(B,2))
        complex(dp)              :: Y(size(A,1),size(A,2))
        
        integer                  :: s(2), nra, nrb, nca, ncb
        integer                  :: lda, ldb, info
        
        integer                  :: ipiv(size(A,1))
        
        external zgesv
        
        s = shape(A)
        nra = s(1)
        nca = s(2)
        s = shape(B)
        nrb = s(1)
        ncb = s(2)
        lda = max(1,nra)
        ldb = max(1,nrb)
        
            
        X = B
        Y = A
        call zgesv(nra,ncb,Y,lda,ipiv,X,ldb,info)
        
        if (info .ne. 0) then
            print *,"Error in the linear system solver!"
        end if
    end function    

    !---------------------------------------------------!
    ! Triple Matrix Product                             !
    !                                                   !
    ! D = A.B.C                                         !
    !                                                   !
    ! A: Complex                                        !
    ! B: Real                                           !
    ! C: Complex                                        !
    !---------------------------------------------------!
    function tri(A,B,C) result(D)
        implicit none
        
        complex(dp),intent(in)     :: A(:,:)
        complex(dp),intent(in)     :: B(:,:)
        complex(dp),intent(in)     :: C(:,:)
        complex(dp)                :: D(size(A,1),size(C,2))
        
        
        ! Hu & Shing matrix chain multiplication
        ! [m] ---B--- [k]
        !  |           |
        !  A           C
        !  |           |
        ! [n]---ABC-- [l]
        
        integer :: n,m,k,l
        integer :: s(2)

        
        s = shape(A)
        n = s(1)
        m = s(2)
        s = shape(C)
        k = s(1)
        l = s(2)
        
        if (n*m*k + n*k*l < m*k*l + n*m*l) then
            D = matmul(matmul(A,B),C)
        else    
            D = matmul(A,matmul(B,C))
        end if
    end function
        
    
    !---------------------------------------------------!
    ! For a diagonal matrix B:                          !
    ! C = B.A.B+ if invert = false                      !
    ! C = B+.A.B if invert = true                       !
    !---------------------------------------------------!
    function d2tri(A,B,invert) result(C)
        use omp_lib
        implicit none
        
        complex(dp),intent(in)  :: A(:,:)
        complex(dp),intent(in)  :: B(:,:)
        logical,intent(in)      :: invert
        complex(dp)             :: C(size(A,1),size(A,2))
        
        integer     :: i,j
        complex(dp) :: p
                
        !$OMP PARALLEL DO
        do j = 1,size(B,1)
            p = B(j,j)
            if (invert .eqv. .false.) then
                p = conjg(p)
            end if
            do i = 1,size(A,1)
                C(i,j) = p*A(i,j)
                if (invert .eqv. .false.) then
                    C(i,j) = C(i,j)*B(i,i)
                else
                    C(i,j) = C(i,j)*conjg(B(i,i))
                end if
            end do
        end do
        !$OMP END PARALLEL DO
    end function    

    !---------------------------------------------------!
    ! G = [(E+i.eta).I - H - S1 - S2]^(-1)                    !
    !---------------------------------------------------!
    function GAssembly3(E,H,S1,S2,eta) result(G)
        implicit none
        
        real,intent(in)         :: E
        complex(dp),intent(in)  :: H(:,:)
        complex(dp),intent(in)  :: S1(:,:)
        complex(dp),intent(in)  :: S2(:,:)
        real(dp),intent(in)     :: eta
        complex(dp)             :: G(size(H,1),size(H,2))
        
        integer             :: i
            
        ! G = (E+i.eta).I - H - S1 - S2
        G = -(H + S1 + S2)
        do i = 1,size(H,2)
            G(i,i) = G(i,i) + E + cmplx(0,1)*eta
        end do
        
        G = inv(G)
    end function     
    
    !---------------------------------------------------!
    ! G = [(E+i.eta).I - H - S]^(-1)                    !
    !---------------------------------------------------!
    function GAssembly2(E,H,S,eta) result(G)
        implicit none
        
        real,intent(in)         :: E
        complex(dp),intent(in)  :: H(:,:)
        complex(dp),intent(in)  :: S(:,:)
        real(dp),intent(in)     :: eta
        complex(dp)             :: G(size(H,1),size(H,2))
        
        integer             :: i
            
        ! G = (E+i.eta).I - H - S
        G = -(H + S)        
        do i = 1,size(H,2)
            G(i,i) = G(i,i) + E + cmplx(0,1)*eta
        end do
        
        G = inv(G)
    end function 
    
    !---------------------------------------------------!
    ! G = [(E+i.eta).I - H]^(-1)                        !
    !---------------------------------------------------!
    function GAssembly1(E,H,eta) result(G)
        implicit none
        
        real,intent(in)         :: E
        complex(dp),intent(in)  :: H(:,:)
        real(dp),intent(in)     :: eta
        complex(dp)             :: G(size(H,1),size(H,2))
        
        integer             :: i
            
        ! G = (E+i.eta).I - H 
        G = -H
        do i = 1,size(H,2)
            G(i,i) = G(i,i) + E + cmplx(0,1)*eta
        end do
        
        G = inv(G)
    end function     
       
    
end module
