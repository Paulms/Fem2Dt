MODULE system
  !
  ! Metodos iterativos
  !
  USE decimal
  USE tipos
  !
  PRIVATE
  !
  PUBLIC system_solution
  !
CONTAINS
  !
  SUBROUTINE system_solution(A,metodo,xx)
    !
    ! meta-subroutina (general)
    !
    TYPE(sparse)              :: A
    INTEGER, INTENT(in)       :: metodo
    REAL(kind=dp)             :: xx(:)
    INTEGER                   :: ierr
    !
    ierr = 0
    !
    SELECT CASE(metodo)
       !
    CASE(1) !metodo directo (usa SuperLU)
       !
       CALL sol_superlu(A,xx)
       !
    CASE(2) ! metodo directo (usa PARDISO)
       !
       CALL sol_pardiso(A,xx)
       !
    CASE(3) ! metodo iterativo (GMRES precondicionado por ILUT)
       !
       CALL sol_gmres(A,xx)
       !
    CASE(4) ! metodo iterativo (BCGSTAB precondicionado por ILUT)
       !
       CALL sol_bcgstab(A,xx)
       !
    CASE(5) ! metodo iterativo (FOM precondicionado por ILUT)
       !
       CALL sol_fom(A,xx)
       !
    CASE default
       PRINT*,'Metdodo aun no implementado ( M: sistemas S:sistemas_solucion)!! ,  metodo =', metodo
    END SELECT
    !
  END SUBROUTINE system_solution
  !
  SUBROUTINE sol_superlu(A,xx)
    !
    ! Subrutina para usar la biblioteca SuperLU (i.e. LU para matrices no simetricas)
    !
    TYPE(sparse)      :: A
    REAL(kind=dp)     :: xx(:)
    INTEGER           :: nrhs,ldb,info,ierr,iopt
    INTEGER(kind=8)   :: factors
    !
    nrhs = 0; ldb = 0; info = 0; ierr = 0; iopt = 0; factors = 0
    !
    ldb  = A%n
    nrhs = 1
    !
    ! Primer paso: Factorizacion de la matriz
    !
    iopt = 1
    !
    CALL c_fortran_dgssv(iopt,A%n,A%nzero,nrhs,A%aa,A%row,A%column,xx,ldb,factors,info)
    !
    IF (info /= 0) THEN
       WRITE (*,*) 'Problemas con la factorizacion LU. INFO = ',info
    ENDIF
    !
    ! Segundo paso: resuelve el sistema usando la factorizacion
    !
    iopt = 2
    !
    CALL c_fortran_dgssv(iopt,A%n,A%nzero,nrhs,A%aa,A%row,A%column,xx,ldb,factors,info)
    !
    IF (info /= 0) THEN
       WRITE (*,*) 'Problemas con la solucion del sistema. INFO = ',info
    ENDIF
    !
    ! Tercer paso: libera la memoria utilizada por SuperLU
    !
    iopt = 3
    !
    CALL c_fortran_dgssv(iopt,A%n,A%nzero,nrhs,A%aa,A%row,A%column,xx,ldb,factors,info)
    !
    IF (info /= 0) THEN
       WRITE (*,*) 'Problemas con la dealocacion de arreglos de SuperLU. INFO = ',info
    ENDIF
    !
  END SUBROUTINE sol_superlu
  !
  SUBROUTINE sol_pardiso(A,xx)
    !
    ! Subrutina para usar la biblioteca PARDISO (i.e. LU para matrices no simetricas)
    !
    TYPE(sparse)      :: A
    REAL(kind=dp)     :: xx(:),bb(SIZE(xx))
    INTEGER           :: maxfct, mnum, mtype, phase, nrhs, error, msglvl
    INTEGER(kind=8)   :: pt(64)
    INTEGER           :: iparm(64)
    REAL(kind=dp)     :: dparm(64) 

    maxfct = 0; mnum = 0; mtype = 0; phase = 0; nrhs = 0; error = 0; msglvl = 0
    pt = 0; iparm = 0; dparm = 0.0_dp; bb = 0.0_dp
    !
    nrhs   = 1
    maxfct = 1
    mnum   = 1
    bb     = xx
    !
    !  Setup Pardiso control parameters und initialize the solvers     
    !  internal adress pointers. This is only necessary for the FIRST   
    !  call of the PARDISO solver.                                     
    !     
    mtype  = 11  ! real and nonsymmetric, complete supernode pivoting
    solver = 0   ! use sparse direct method
    !     
    ! PARDISO license check and initialize solver
    !
    CALL pardisoinit(pt, mtype, solver, iparm, dparm, error)
    !  .. Numbers of Processors ( value of OMP_NUM_THREADS )
    iparm(3)  = 1
    !
    PRINT*,'paso pardisoinit'
    IF (error /= 0) THEN
       IF (error == -10 ) WRITE(*,*) 'No license file found'
       IF (error == -11 ) WRITE(*,*) 'License is expired'
       IF (error == -12 ) WRITE(*,*) 'Wrong username or hostname'
       STOP
    ELSE
       WRITE(*,*) '[PARDISO]: License check was successful ... '
    END IF
    !
    !  .. pardiso_chk_matrix(...)
    !     Checks the consistency of the given matrix.
    !     Use this functionality only for debugging purposes
!!$    CALL pardiso_chkmatrix(mtype, A%n, A%aa, A%row, A%column, error)
!!$    IF (error /= 0) THEN
!!$       WRITE(*,*) 'The following ERROR was detected (1): ', error
!!$       STOP
!!$    ENDIF
    !
    PRINT*,'paso  pardiso_chkmatrix'
    !
    ! ..  pardiso_chkvec(...)
    !     Checks the given vectors for infinite and NaN values
    !     Input parameters (see PARDISO user manual for a description):
    !     Use this functionality only for debugging purposes
!!$    CALL pardiso_chkvec (A%n, nrhs, xx, error)
!!$    IF (error /= 0) THEN
!!$       WRITE(*,*) 'The following ERROR was detected (2): ', error
!!$       STOP
!!$    ENDIF
    !
    PRINT*,'paso  pardiso_chkvec'
    ! ..  pardiso_printstats(...) 
    !     prints information on the matrix to STDOUT.
    !     Use this functionality only for debugging purposes
!!$    CALL pardiso_printstats (mtype, A%n, A%aa, A%row, A%column, nrhs, bb, error)
!!$    IF (error .NE. 0) THEN
!!$       WRITE(*,*) 'The following ERROR was detected (3): ', error
!!$       STOP
!!$    ENDIF
    !
    PRINT*,'paso  pardiso_printstats'
    !..   Reordering and Symbolic Factorization, This step also allocates
    !     all memory that is necessary for the factorization
    !
    phase     = 11     ! only reordering and symbolic factorization
    msglvl    = 1      ! with statistical information
    iparm(33) = 1      ! compute determinant
    !
    CALL pardiso(pt, maxfct, mnum, mtype, phase, A%n, A%aa, A%row, A%column,&
                  idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)
    !
    PRINT*,'paso  pardiso reordenamiento'
    WRITE(*,*) 'Reordering completed ... '
    !
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected (4): ', error
       STOP
    END IF
    !
    WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
    WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)
    !
    !.. Factorization.
    !
    phase     = 22  ! only factorization
    CALL pardiso (pt, maxfct, mnum, mtype, phase, A%n, A%aa, A%row, A%column,& 
                  idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm) 
    !
    IF (iparm(33) == 1)  THEN
       WRITE(*,*) 'Log of determinant is  ',  dparm(33)
    ENDIF
    !
    WRITE(*,*) 'Factorization completed ... '
    IF (error /= 0) THEN
       WRITE(*,*) 'The following ERROR was detected (5): ', error
       STOP
    ENDIF
    !
    !.. Back substitution and iterative refinement
    iparm(8)  = 1   ! max numbers of iterative refinement steps
    phase     = 33  ! only solve
    !
    CALL pardiso (pt, maxfct, mnum, mtype, phase, A%n, A%aa, A%row, A%column,& 
                  idum, nrhs, iparm, msglvl,bb, xx, error, dparm) 
    !
    WRITE(*,*) 'Solve completed ...  '
    !
    !.. Termination and release of memory
    phase     = -1           ! release internal memory
    CALL pardiso (pt, maxfct, mnum, mtype, phase, A%n, ddum, idum, idum,& 
                  idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm) 
    !
  END SUBROUTINE sol_pardiso
  !
  SUBROUTINE sol_gmres(A,xx)
    !
    ! Subrutina para usar GMRES de la libreria SPARSKIT2
    !
    TYPE(sparse)               :: A
    REAL(kind=dp)              :: xx(:),bb(size(xx))
    !
    INTEGER                    :: ipar(16),lfil,nwk,nrow,ierr,maxits,its
    INTEGER, ALLOCATABLE       :: iw(:),ju(:), jlu(:)
    REAL(kind=dp), ALLOCATABLE :: wk(:),alu(:)
    REAL(kind=dp)              :: fpar(16), tol,  res
    !
    ipar = 0; fpar = 0.0_dp; tol = 0.0_dp; lfil = 0; nwk = 0; nrow = 0; ierr = 0; maxits = 0
    tol = 0.0_dp; bb = 0.0_dp

    lfil   = 16
    maxits = A%n
    nwk    = 2*lfil*A%n+A%n+1
    !
    bb = xx
    !
    ALLOCATE(iw(A%n*2),ju(A%n), alu(nwk), jlu(nwk))
    !
    iw = 0; ju = 0; alu = 0.0_dp; jlu = 0
    !
    ipar(1) = 0
    ipar(2) = 2
    ipar(3) = 1
    ipar(5) = 16
    ipar(4) = (A%n + 3)*(ipar(5)+2) + ((ipar(5)+1)*ipar(5))/2
    !
    ALLOCATE(wk(ipar(4)))
    wk = 0.0_dp
    !
    ipar(6) = maxits
    fpar(1) = 1.0e-05_dp
    fpar(2) = 1.0e-10_dp
    !
    !   definiendo el precondicionador ILUT
    !    
    tol = 1.0e-05_dp
    !
    CALL ilut(A%n,A%aa,A%column,A%row,lfil,tol,alu,jlu,ju,nwk,wk,iw,ierr)
    !
    PRINT*,'ierr de ilut = ', ierr
    ipar(1) = 0
    ipar(2) = 2
    !
    its = 0
    res = 0.0_dp
    !
    main: DO
       !
       CALL gmres(A%n,bb,xx,ipar,fpar,wk)
       !
       IF(ipar(7) /= its) its = ipar(7)
       !
       res = fpar(5)
       !
       SELECT CASE(ipar(1))
       CASE(1)
          CALL amux(A%n, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(2)
          CALL atmux(A%n, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(3,5)
          CALL lusol(A%n,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(4,6)
          CALL lutsol(A%n,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(0)
          PRINT *, 'Se alcanzo la convergencia!!'
          EXIT main
       CASE(-1)
          PRINT *, 'El solver itero mas alla de lo permitido. No hay convergencia.'
          STOP
       CASE(-2)
          PRINT *, 'El solver no tiene suficiente espacio para trabajar'
          PRINT *, '(espacio minimo que se requiere es de ', ipar(4),' elementos)'
          STOP
       CASE(-3)
          PRINT *, 'Se dividio por cero en algun momento. No hay convergencia.'
          STOP
       CASE default
          PRINT *, 'Se termino el solver con el codigo de error =', ipar(1)
          STOP
       END SELECT
       !
    END DO main
    !
    PRINT*,'iteraciones = ', ipar(7)
    PRINT*,'residuo     = ', fpar(5)
    !
    DEALLOCATE(iw,wk,ju,alu,jlu, STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas liberando la memoria!! (M:sistemas S: sol_gmres)'
       STOP
    END IF
    !
  END SUBROUTINE sol_gmres
  !
  SUBROUTINE sol_bcgstab(A,xx)
    !
    ! Driver to use Bi Conjugate Gradient stabilized (BCGSTAB) from SPARSKIT2 library
    !
    TYPE(sparse)               :: A
    REAL(kind=dp)              :: xx(:),bb(size(xx))
    !
    INTEGER                    :: ipar(16),lfil,nwk,nrow,ierr,maxits,its
    INTEGER, ALLOCATABLE       :: iw(:),ju(:), jlu(:)
    REAL(kind=dp), ALLOCATABLE :: wk(:),alu(:)
    REAL(kind=dp)              :: fpar(16), tol,  res
    !
    ipar = 0; fpar = 0.0_dp; tol = 0.0_dp; lfil = 0; nwk = 0; nrow = 0; ierr = 0; maxits = 0
    tol = 0.0_dp; bb = 0.0_dp

    lfil   = 16
    maxits = A%n
    nwk    = 2*lfil*A%n+A%n+1
    !
    bb = xx
    !
    ALLOCATE(iw(A%n*2),ju(A%n), alu(nwk), jlu(nwk),wk(8*A%n))
    !
    iw = 0; ju = 0; alu = 0.0_dp; jlu = 0; wk = 0.0_dp
    !
    ipar(1) = 0
    ipar(2) = 2
    ipar(3) = 1
    ipar(4) = 8*A%n
    ipar(5) = 16
    ipar(6) = maxits
    !
    fpar(1) = 1.0e-05_dp
    fpar(2) = 1.0e-10_dp
    !
    !   set-up the preconditioner ILUT(15, 1E-4) ! new definition of lfil
    !    
    tol = 1.0e-05_dp
    !
    CALL ilut(A%n,A%aa,A%column,A%row,lfil,tol,alu,jlu,ju,nwk,wk,iw,ierr)
    !
    PRINT*,'ierr de ilut = ', ierr
    ipar(1) = 0
    ipar(2) = 2
    !
    its = 0
    res = 0.0_dp
    !
    main: DO
       !
       CALL bcgstab(A%n,bb,xx,ipar,fpar,wk)
       !
       IF(ipar(7) /= its) its = ipar(7)
       !
       res = fpar(5)
       !
       SELECT CASE(ipar(1))
       CASE(1)
          CALL amux(A%n, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(2)
          CALL atmux(A%n, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(3,5)
          CALL lusol(A%n,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(4,6)
          CALL lutsol(A%n,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(0)
          PRINT *, 'Iterative solver has satisfied convergence test.'
          EXIT main
       CASE(-1)
          PRINT *, 'Iterative solver has iterated too many times.'
          STOP
       CASE(-2)
          PRINT *, 'Iterative solver was not given enough work space.'
          PRINT *, 'The work space should at least have ', ipar(4),' elements.'
          STOP
       CASE(-3)
          PRINT *, 'Iterative solver is facing a break-down.'
          STOP
       CASE default
          PRINT *, 'Iterative solver terminated. code =', ipar(1)
          STOP
       END SELECT
       !
    END DO main
    !
    PRINT*,'iteraciones = ', ipar(7)
    PRINT*,'residuo     = ', fpar(5)
    !
    DEALLOCATE(iw,wk,ju,alu,jlu, STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas liberando la memoria!! (M:sistemas S: sol_bcgstab)'
       STOP
    END IF
    !
  END SUBROUTINE sol_bcgstab
  !
  SUBROUTINE sol_fom(A,xx)
    !
    ! Driver to use Full Orthogonalization Method (FOM) from SPARSKIT2 library
    !
    TYPE(sparse)               :: A
    REAL(kind=dp)              :: xx(:),bb(size(xx))
    !
    INTEGER                    :: ipar(16),lfil,nwk,nrow,ierr,maxits,its
    INTEGER, ALLOCATABLE       :: iw(:),ju(:), jlu(:)
    REAL(kind=dp), ALLOCATABLE :: wk(:),alu(:)
    REAL(kind=dp)              :: fpar(16), tol,  res
    !
    ipar = 0; fpar = 0.0_dp; tol = 0.0_dp; lfil = 0; nwk = 0; nrow = 0; ierr = 0; maxits = 0
    tol = 0.0_dp; bb = 0.0_dp

    lfil   = 16
    maxits = A%n
    nwk    = 2*lfil*A%n+A%n+1
    !
    bb = xx
    !
    ALLOCATE(iw(A%n*2),ju(A%n), alu(nwk), jlu(nwk))
    !
    iw = 0; ju = 0; alu = 0.0_dp; jlu = 0
    !
    ipar(1) = 0
    ipar(2) = 2
    ipar(3) = 1
    ipar(5) = 16
    ipar(4) = (A%n + 3)*(ipar(5)+2) + ((ipar(5)+1)*ipar(5))/2 
    !
    ALLOCATE(wk(ipar(4)))
    wk = 0.0_dp
    !
    ipar(6) = maxits
    fpar(1) = 1.0e-05_dp
    fpar(2) = 1.0e-10_dp
    !
    !   set-up the preconditioner ILUT(15, 1E-4) ! new definition of lfil
    !    
    tol = 1.0e-05_dp
    !
    CALL ilut(A%n,A%aa,A%column,A%row,lfil,tol,alu,jlu,ju,nwk,wk,iw,ierr)
    !
    PRINT*,'ierr de ilut = ', ierr
    ipar(1) = 0
    ipar(2) = 2
    !
    its = 0
    res = 0.0_dp
    !
    main: DO
       !
       CALL fom(A%n,bb,xx,ipar,fpar,wk)
       !
       IF(ipar(7) /= its) its = ipar(7)
       !
       res = fpar(5)
       !
       SELECT CASE(ipar(1))
       CASE(1)
          CALL amux(A%n, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(2)
          CALL atmux(A%n, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(3,5)
          CALL lusol(A%n,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(4,6)
          CALL lutsol(A%n,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(0)
          PRINT *, 'Iterative solver has satisfied convergence test.'
          EXIT main
       CASE(-1)
          PRINT *, 'Iterative solver has iterated too many times.'
          STOP
       CASE(-2)
          PRINT *, 'Iterative solver was not given enough work space.'
          PRINT *, 'The work space should at least have ', ipar(4),' elements.'
          STOP
       CASE(-3)
          PRINT *, 'Iterative solver is facing a break-down.'
          STOP
       CASE default
          PRINT *, 'Iterative solver terminated. code =', ipar(1)
          STOP
       END SELECT
       !
    END DO main
    !
    PRINT*,'iteraciones = ', ipar(7)
    PRINT*,'residuo     = ', fpar(5)
    !
    DEALLOCATE(iw,wk,ju,alu,jlu, STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas liberando la memoria!! (M:sistemas S: sol_fom)'
       STOP
    END IF
    !
  END SUBROUTINE sol_fom
  !
END MODULE system






