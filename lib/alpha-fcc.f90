subroutine fcc 
    use variables_definition 
    implicit none 
    integer, parameter ::  NC = 3
    integer, parameter ::  N = 4 * NC ** 3.  
    real, parameter :: RROOT3 = 0.5773503
    real RX(N), RY(N), RZ(N), EX(N), EY(N), EZ(N)
    real CELL, CELL2 
    integer I,IX,IY,IREF,M

!   CALCULATE THE SIDE OF THE UNIT CELL ** 
    CELL = 1.0 / (dble(NC))
    CELL2 = 0.5 * CELL 
!    ** BUILD THE UNIT CELL **

!    ** SUBLATTICE A **

        RX(1) =  0.0
        RY(1) =  0.0
        RZ(1) =  0.0
        EX(1) =  RROOT3
        EY(1) =  RROOT3
        EZ(1) =  RROOT3

!    ** SUBLATTICE B **

        RX(2) =  CELL2
        RY(2) =  CELL2
        RZ(2) =  0.0
        EX(2) =  RROOT3
        EY(2) = -RROOT3
        EZ(2) = -RROOT3

!    ** SUBLATTICE C **

        RX(3) =  0.0
        RY(3) =  CELL2
        RZ(3) =  CELL2
        EX(3) = -RROOT3
        EY(3) =  RROOT3
        EZ(3) = -RROOT3

!    ** SUBLATTICE D **

        RX(4) =  CELL2
        RY(4) =  0.0
        RZ(4) =  CELL2
        EX(4) = -RROOT3
        EY(4) = -RROOT3
        EZ(4) =  RROOT3

!    ** CONSTRUCT THE LATTICE FROM THE UNIT CELL **


    return 
end subroutine fcc 
