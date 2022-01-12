program test
    use GreenCheetah
    implicit none
    
    integer,parameter   :: dp = kind(1.d0)
    type(Green) :: G
    real(dp)        :: r(10,10)

    call random_number(r)
    G = Green(r,1E-6,1E-6,0.041,.false.)
end program
