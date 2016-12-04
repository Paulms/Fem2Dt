MODULE util
  !
  ! Modulo que contiene diferentes programas que no pudieron ser 
  ! clasificados en otra parte
  !
  !
  USE decimal
  USE tipos
  !
  IMPLICIT NONE
  !
  REAL(kind=dp),PARAMETER:: unotres=1.0_dp/3.0_dp
  !
CONTAINS
  !
  SUBROUTINE invert(matrix)
    !
    ! Subrutina que invierte una matriz
    !
    ! IN: 
    ! matrix = matriz de 2x2 o 3x3
    !
    ! OUT:
    ! matrix = la inversa de la matriz dada (sobre-escribe)
    !
    IMPLICIT NONE
    REAL (kind=dp),INTENT(inout):: matrix(:,:)
    INTEGER                     :: nn
    REAL (kind=dp)              :: det
    REAL(kind=dp)               :: temp(SIZE(matrix,1),SIZE(matrix,1))
    !
    nn = SIZE(matrix,1)
    !
    det = 0.0_dp; temp = 0.0_dp
    !
    det = determinant(matrix)
    !
    IF(ABS(det) < 1.0e-12_dp) THEN
       PRINT*,"determinante igual a cero!! M: util S: invert"
       STOP
    ENDIF
    !
    SELECT CASE(nn)
    CASE(2)! matrix 2x2
       !
       temp(1,1) = matrix(2,2)
       temp(1,2) =-matrix(1,2)        
       temp(2,1) =-matrix(2,1)
       temp(2,2) = matrix(1,1)
       !
    CASE(3)! matrix 3x3
       !
       temp(1,1) =  matrix(2,2)  * matrix(3,3) - matrix(2,3) * matrix(3,2)
       temp(1,2) = -matrix(1,2)  * matrix(3,3) + matrix(1,3) * matrix(3,2)
       temp(1,3) =  matrix(1,2)  * matrix(2,3) - matrix(1,3) * matrix(2,2)
       temp(2,1) =  matrix(3,1)  * matrix(2,3) - matrix(2,1) * matrix(3,3)
       temp(2,2) = -matrix(3,1)  * matrix(1,3) + matrix(1,1) * matrix(3,3)
       temp(2,3) =  matrix(2,1)  * matrix(1,3) - matrix(1,1) * matrix(2,3)
       temp(3,1) = -matrix(3,1)  * matrix(2,2) + matrix(2,1) * matrix(3,2)
       temp(3,2) =  matrix(3,1)  * matrix(1,2) - matrix(1,1) * matrix(3,2)
       temp(3,3) = -matrix(2,1)  * matrix(1,2) + matrix(1,1) * matrix(2,2)
       !
    CASE default
       PRINT*,'Problema de dimension de la matriz. M: util S: invert'
       STOP
    END SELECT
    !
    matrix = temp/det
    !
  END SUBROUTINE invert
  !
  FUNCTION determinant(matrix) RESULT(det)
    !
    ! Funcion que calcula el determinante de una matriz de 2x2 3x3
    !
    ! IN:
    ! matrix = matriz de   2x2  y/o 3x3.
    !
    ! OUT:
    ! det = determinante de la matriz matrix
    !
    IMPLICIT NONE    
    REAL (kind=dp)           :: det
    REAL (kind=dp),INTENT(in):: matrix(:,:)
    INTEGER                  :: nn 
    !
    det = 0.0_dp
    nn  = SIZE(matrix,1)
    !
    SELECT CASE (nn)
    CASE (2)
       det = matrix(1,1)*matrix(2,2) - matrix(1,2) * matrix(2,1)
    CASE (3)
       !
       det = matrix(1,1)*matrix(2,2)*matrix(3,3)-matrix(1,1)*matrix(2,3)*matrix(3,2)+&
             matrix(2,1)*matrix(3,2)*matrix(1,3)-matrix(2,1)*matrix(1,2)*matrix(3,3)+&
             matrix(3,1)*matrix(1,2)*matrix(2,3)-matrix(3,1)*matrix(2,2)*matrix(1,3)
       !
    CASE default
       PRINT*,' Problema de dimension de la matriz!! M: util F: determinant'
    END SELECT
    !
  END FUNCTION determinant
  !
  REAL(kind=dp) FUNCTION paramametro(nu,sigma,h) RESULT(tau)
    !
    ! Esta funcion calcula el parametro de estabilizacion.
    !
    REAL(kind=dp),INTENT(in):: nu,sigma,h
    REAL(kind=dp)           :: lambda,beta
    !
    tau    = 0.0_dp
    beta   = 0.0_dp
    lambda = 0.0_dp
    !
    IF(ABS(sigma) < 1.0e-15_dp) THEN
       !
       beta = 24.0_dp*nu
       !
    ELSE
       !
       lambda = 12.0_dp*nu/(sigma*h*h)
       !
       IF ( lambda >= 1.0_dp) THEN
          beta = 24.0_dp*nu
       ELSE
          beta = sigma*h*h + 12*nu
       END IF
       !
    END IF
    !
    !    tau  = (h*h)/beta
    tau  = (h*h)/(8.0_dp*nu)
    !
  END FUNCTION paramametro
  !
  REAL(kind=dp) FUNCTION param_adr(nu,vect,sigma,h_T) RESULT(tau)
    !
    REAL(kind=dp),INTENT(in):: nu,vect(:),sigma,h_T
    REAL(kind=dp)           :: Pe(2),alpha,beta
    !
    alpha = 0.0_dp; beta = 0.0_dp; Pe =0.0_dp; tau = 0.0_dp
    !
    alpha = sigma*h_T**2
    !
    IF(ABS(sigma) < 1.0e-15_dp) THEN
       alpha = 0.0_dp
    ELSE
       Pe(1) = (6.0_dp*nu)/(sigma*h_T**2)
       alpha = alpha*MAX(1.0_dp,Pe(1))
    END IF
    !
    beta  = 6.0_dp*nu
    !
    IF(nu < 1.0e-25_dp) THEN
       beta = 0.0_dp
    ELSE
       Pe(2) = (SQRT(DOT_PRODUCT(vect,vect))*h_T)/(3.0_dp*nu)
       beta = beta*MAX(1.0_dp,Pe(2))
    END IF
    !
    tau  = (h_T**2)/(alpha+beta)
    !
  END FUNCTION param_adr
  !
  REAL(kind=dp) FUNCTION param_lps(nu,vect,sigma,h_T) RESULT(tau)
    !
    REAL(kind=dp),INTENT(in):: nu,vect(:),sigma,h_T
    REAL(kind=dp)           :: Pe,alpha,beta,pi
    !
    alpha = 0.0_dp; beta = 0.0_dp; Pe = 0.0_dp; tau = 0.0_dp; pi = 0.0_dp
    !
    pi = 4.0_dp*ATAN(1.0_dp)
    !
    beta = 11.0_dp/(3.0_dp*pi**2)
    !
    Pe = (SQRT(DOT_PRODUCT(vect,vect))*h_T)/(beta*nu)
    !PRINT*,'Pe = ', Pe
    !    
    alpha = MAX(1.0_dp,Pe)
    !
    tau = (3.0_dp*pi**2)/(2.0_dp*nu*alpha)
    !
    !PRINT*,'tau = ', tau,'  nu = ',nu
    !tau = 0.7_dp*( 1.0_dp/(SQRT(DOT_PRODUCT(vect,vect))*h_T))
    !
  END FUNCTION param_lps
  !
  REAL(kind=dp) FUNCTION hvf(x,vect) RESULT(ht)
    !
    REAL(kind=dp),INTENT(in):: x(:,:),vect(:)
    REAL(kind=dp)           :: f1,f2,anorm,a1n,a2n,xi,yi,cf1,cf2
    INTEGER                 :: i,j,k
    !
    f1 = -1.0_dp
    f2 = -1.0_dp
    i = 1
    !
    anorm = SQRT(DOT_PRODUCT(vect,vect))
    !
    IF(anorm < 1.0e-15_dp) THEN
       ht = diametro(x)
       RETURN
    END IF
    !
    a1n = vect(1)/anorm
    a2n = vect(2)/anorm
    !
    DO !WHILE(f1*f2 > 0.0_dp)
       !
       j  = MOD(i+1,4)
       IF (j ==0) j = 1
       !
       f1 = line(x(1,i),x(2,i),x(1,j),x(2,j),a1n,a2n)
       !
       k  = MOD(j+1,4)
       IF (k == 0) k = 1
       !
       f2 = line(x(1,i),x(2,i),x(1,k),x(2,k),a1n,a2n)  
       !
       i  = i+1
       !
       IF(f1*f2 <= 0.0_dp) EXIT
       !
    END DO
    !
    IF (f1 == 0.0_dp) THEN
       ht = SQRT((x(1,j)-x(1,i-1))**2+(x(2,j)-x(2,i-1))**2)
    ELSE
       IF (f2 == 0.0_dp) THEN
          ht = SQRT((x(1,k)-x(1,i-1))**2+(x(2,k)-x(2,i-1))**2)
       ELSE
          IF (a1n == 0.0_dp) THEN
             yi = x(2,i)
             xi = (yi-x(2,j))*(x(1,k)-x(1,j))/(x(2,k)-x(2,j))+x(1,j)
          ELSE
             IF (a2n == 0.0_dp) THEN          
                xi = x(1,i)
                yi = (xi-x(2,j))*(x(2,k)-x(2,j))/(x(1,k)-x(1,j))+x(2,j)
             ELSE 
                IF ((x(1,j)-x(1,k)) == 0.0_dp) THEN
                   xi = x(1,k)        
                   yi = (a2n/a1n)*(xi-x(1,i-1))+x(2,i-1)
                ELSE
                   cf1 = (x(2,j)-x(2,k))/(x(1,j)-x(1,k))
                   cf2 = a2n/a1n
                   xi  = (cf2*x(1,i-1)+(x(2,k)-x(2,i-1))-cf1*x(1,k))/(cf2-cf1)
                   yi  = cf2*(xi-x(1,i-1))+x(2,i-1)
                ENDIF
             ENDIF
          ENDIF
          ht  = SQRT((xi-x(1,i-1))**2+(yi-x(2,i-1))**2)
       ENDIF
    ENDIF
    !
  END FUNCTION hvf
  !
  REAL(kind=dp) FUNCTION line(x1,y1,x2,y2,a1,a2)
    !
    REAL(kind=dp),INTENT(in)::x1,x2,y1,y2,a1,a2
    !
    line = -a2*(x2-x1)+a1*(y2-y1)
    !
  END FUNCTION line
  !
  SUBROUTINE util_input(malla,fisica,frontera,metodo,solver,aprox,visu,name_visu, tiempo)
    !
    INTEGER           :: metodo,solver,aprox,visu,i,iunit
    TYPE(bc)          :: frontera
    CHARACTER(LEN=32) :: name_visu
    TYPE(mesh)        :: malla
    TYPE(param)       :: fisica
    INTEGER           :: nT
    REAL(kind=dp)     :: Tf
    TYPE(time)        :: tiempo
    !    
    iunit = 0
    !    
    ! Archivo desde donde se leeran los datos basicos...
    ! 
    CALL util_get_unit(iunit)
    OPEN(unit=iunit,file='input',status='old',action='read')    
    !
    ! Lectura de la Geometria y de la Interpolacion a usar
    !
    CALL read_geometria(iunit,malla)
    !
    !
    ! Lectura de los parametros fisicos: 
    ! ndominios               = numero de dominios diferentes
    ! nu(i), aa(i,:),sigma(i) = viscosidad, adveccion, reaccion del i-esimo material
    ! ref_dominio(i)          = numero de referencia del i-esimo dominio
    !
    READ(iunit,*) fisica%ndominios
    !
    ALLOCATE(fisica%nu(fisica%ndominios),fisica%aa(fisica%ndominios,malla%ndim),fisica%sigma(fisica%ndominios),&
             fisica%ref_dominio(fisica%ndominios))
    !
    fisica%nu = 0.0_dp; fisica%aa = 0.0_dp; fisica%sigma = 0.0_dp; fisica%ref_dominio = 0
    !
    DO i = 1,fisica%ndominios
       READ(iunit,*) fisica%ref_dominio(i),fisica%nu(i),fisica%aa(i,:),fisica%sigma(i)
    END DO
    !
    ! Lectura del numero total de condiciones de Dirichlet
    !
    frontera%n_dirichlet = 0; frontera%n_neumann = 0
    !
    READ(iunit,*) frontera%n_dirichlet
    !
    IF(frontera%n_dirichlet /= 0) THEN
       !
       ALLOCATE(frontera%DD(frontera%n_dirichlet),frontera%fix(malla%element%dof,frontera%n_dirichlet))
       frontera%DD  = 0
       frontera%fix = 0
       !
       ! Lectura de los numeros de referencias de la cond. Dirichlet y de los bloqueos respectivos
       !
       DO i=1,frontera%n_dirichlet
          READ(iunit,*) frontera%DD(i), frontera%fix(:,i)
       END DO
       !
    END IF
    !
    READ(iunit,*) frontera%n_neumann
    !
    IF(frontera%n_neumann /= 0) THEN
       !
       ALLOCATE(frontera%NN(frontera%n_neumann))
       frontera%NN  = 0
       !
       ! Lectura de los numeros de referencias de la cond. Neumann
       !
       DO i=1,frontera%n_neumann
          READ(iunit,*) frontera%NN(i)
       END DO
       !
    END IF
    !
    ! Metodo de elementos finitos a usar
    !
    READ(iunit,*) metodo
    !
    ! Metodo para resolver el sistema lineal
    !
    READ(iunit,*) solver
    !
    ! Tpo de error a calcular (a priori y/o a posteriori)
    !
    READ(iunit,*) aprox
    !
    ! Metodo para visualizar los resultados
    !
    READ(iunit,*) visu
    !
    ! Lectura del nombre del archivo que contendra los resultados
    !
    READ(iunit,'(a32)')  name_visu
    !
    ! Numero de intervalos en los que se divide el tiempo
    !
    READ(iunit,*) tiempo%nn
    !
    ! Tiempo final
    !
    READ(iunit,*) tiempo%final
    !
    CLOSE(iunit)
    !
  END SUBROUTINE util_input
  !
  SUBROUTINE read_geometria(iunit,malla)
    !
    ! subrutina que lee la malla y la interpolacion asociada al
    ! problema.
    !
    TYPE(mesh)        :: malla
    INTEGER           :: i,iunit,iunita,iunitb,iunitc
    CHARACTER(LEN=32) :: name_node,name_face,name_ele
    !
    malla%element%nod = 0; malla%element%dof = 0; malla%element%npi = 0
    !
    READ(iunit, '(a8)') malla%element%name
    !
    READ(iunit,*) malla%element%npi
    !
    READ(iunit,*) malla%element%dof
    !
    READ(iunit,'(a32)') name_node
    !   
    READ(iunit,'(a32)') name_face
    !
    READ(iunit,'(a32)') name_ele
    !
    malla%element%name = TRIM(malla%element%name)
    !
    SELECT CASE (malla%element%name)
    CASE('TRIAP12D')
       malla%element%tipo     ='triangle            '
       malla%element%nod      = 3
       malla%element%nodf     = 2
       malla%element%face_loc = 3
       malla%ndim             = 2 
    CASE('QUADQ12D')
       malla%element%tipo     ='quadrilateral       '
       malla%element%nod      = 4
       malla%element%nodf     = 2
       malla%element%face_loc = 4
       malla%ndim             = 2 
    CASE('TETRP13D')
       malla%element%tipo     ='tetrahedron         '
       malla%element%nod      = 4
       malla%element%nodf     = 3
       malla%element%face_loc = 4
       malla%ndim             = 3
    CASE('HEXAQ13D')
       malla%element%tipo     ='hexahedron          '
       malla%element%nod      = 8
       malla%element%nodf     = 4
       malla%element%face_loc = 6
       malla%ndim             = 3
    CASE default
       PRINT*,'Elemento aun no implementado!!! M: util, S: util_input'
       STOP
    END SELECT
    !    
    CALL util_get_unit(iunita)
    OPEN(unit=iunita,file=TRIM(name_node),status='old',action='read')
    !
    CALL util_get_unit(iunitb)
    OPEN(unit=iunitb,file=TRIM(name_ele),status='old',action='read')
    !
    CALL util_get_unit(iunitc)
    OPEN(unit=iunitc,file=TRIM(name_face),status='old',action='read')
    !
    ! Leyendo la informacion de los nodos
    !
    READ(iunita,*)
    READ(iunita,*) malla%nnode
    !
    ALLOCATE(malla%coord(malla%ndim,malla%nnode),malla%ref_node(malla%nnode))
    !
    malla%coord    = 0.0_dp
    malla%ref_node = 0
    !
    READ(iunita,*)
    READ(iunita,*) (malla%coord(:,i), malla%ref_node(i), i =1, malla%nnode)
    !
    ! Leyendo la informacion de los elementos
    !
    READ(iunitb,*)
    READ(iunitb,*) malla%nelem        
    !
    ALLOCATE(malla%conec(malla%element%nod,malla%nelem),malla%ref_elem(malla%nelem),&
             malla%ele_face(malla%element%face_loc,malla%nelem),malla%neigh(malla%element%face_loc,malla%nelem))
    !
    malla%conec      = 0
    malla%ref_elem   = 0
    malla%neigh      = 0
    !
    READ(iunitb,*)
    READ(iunitb,*) (malla%conec(:,i), malla%ref_elem(i), i =1, malla%nelem) 
    !
    READ(iunitb,*)
    READ(iunitb,*) (malla%ele_face(:,i), i =1, malla%nelem) 
    !
    READ(iunitb,*)
    READ(iunitb,*) (malla%neigh(:,i), i =1, malla%nelem) 
    !
    ! Leyendo la informacion de los lados (2D) o las caras (3D) de los elementos 
    ! que conforman la malla.
    !
    READ(iunitc,*)
    READ(iunitc,*) malla%nedge
    !
    ALLOCATE(malla%conec_face(malla%element%nodf,malla%nedge),malla%ref_face(malla%nedge))
    !
    malla%ref_face   = 0
    malla%conec_face = 0
    !
    READ(iunitc,*)
    READ(iunitc,*) (malla%conec_face(:,i), malla%ref_face(i),i=1, malla%nedge)
    !
    CLOSE(iunita)
    CLOSE(iunitb)
    CLOSE(iunitc)
    !
    !
  END SUBROUTINE read_geometria
  !
  SUBROUTINE h_malla(malla, h_max)
    !
    TYPE(mesh)                 :: malla
    REAL(kind=dp), INTENT(out) :: h_max
    REAL(kind=dp), ALLOCATABLE :: hh(:),coord(:,:)
    INTEGER, ALLOCATABLE       :: mm(:)
    INTEGER                    :: i,j,ierr
    !
    ALLOCATE(hh(malla%nelem),coord(malla%ndim,malla%element%nod), mm(malla%element%nod))
    !
    hh = 0.0_dp; h_max = 0.0_dp
    !
    DO i = 1, malla%nelem
       !
       mm     = 0
       coord  = 0.0_dp
       !
       mm     = malla%conec(:,i)
       !
       DO j=1,malla%element%nod
          coord(:,j) = malla%coord(:,mm(j))
       END DO
       !
       SELECT CASE (malla%element%name)
       CASE('TRIAP12D')
          !
          hh(i) = diametro(coord)
          !
       CASE('QUADQ12D')
          !
          hh(i) = diametro(coord)
          !
       CASE('TETRP13D')
          !
          hh(i)  = 0.0_dp
          !
       CASE('HEXAQ13D')
          !
          hh(i) = 0.0_dp
          !
       CASE default
          PRINT*,'Elemento aun no implementado!!! S: h_malla  M: util'
          STOP
       END SELECT
       !
    END DO
    !
    h_max = MAXVAL(hh)
    !
    DEALLOCATE( hh, coord,mm,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: h_malla   M:util)!!'
       STOP
    END IF

  END SUBROUTINE h_malla
  !
  FUNCTION diametro(coord) !P1 y Q1 en 2D
    !
    REAL(kind=dp),INTENT(in)  :: coord(:,:)
    REAL(kind=dp),ALLOCATABLE :: hh(:), lado(:,:)
    REAL(kind=dp)             :: diametro
    INTEGER                   :: nod,ndim,i, ierr
    !
    nod = 0; ndim = 0; ierr = 0
    !
    ndim = SIZE(coord,1)
    nod  = SIZE(coord,2)
    !
    diametro = 0.0_dp
    !
    ALLOCATE(lado(ndim,6),hh(6))
    !
    lado = 0.0; hh = 0.0_dp
    !
    SELECT CASE(ndim)
    CASE(2)
       SELECT CASE(nod)
       CASE(3) !P1 dimension 2
          !
          lado(:,1) = coord(:,2)- coord(:,1)
          lado(:,2) = coord(:,3)- coord(:,2)
          lado(:,3) = coord(:,1)- coord(:,3)
          ! 
          DO i=1,3
             hh(i)= SQRT(DOT_PRODUCT(lado(:,i),lado(:,i)))
          END DO
          !
       CASE(4) !Q1 dimension 2
          !
          lado(:,1) = coord(:,1)- coord(:,3)
          lado(:,2) = coord(:,4)- coord(:,2)
          !
          hh(1) = SQRT(DOT_PRODUCT(lado(:,1),lado(:,1)))
          hh(2) = SQRT(DOT_PRODUCT(lado(:,2),lado(:,2)))
          !
       CASE default
          PRINT*,'elemento aun no implementado (funcion diametro)'
          STOP
       END SELECT
    CASE(3) !dimension 3
       SELECT CASE(nod)
       CASE(4) !P1 en dimension 3
          !
          lado(:,1) = coord(:,2) - coord(:,1)
          lado(:,2) = coord(:,3) - coord(:,2)
          lado(:,3) = coord(:,1) - coord(:,3)
          lado(:,4) = coord(:,4) - coord(:,1)
          lado(:,5) = coord(:,4) - coord(:,2)
          lado(:,6) = coord(:,4) - coord(:,3)
          !
          DO i=1,6
             hh(i) = SQRT(DOT_PRODUCT(lado(:,i),lado(:,i)))
          END DO
          !
       CASE(8)! Q1 en dimension 3
          !
          lado = 0.0_dp
          hh   = 0.0_dp
          !
       CASE default
          PRINT*,'elemento aun no implementado (funcion diametro)'
          STOP
       END SELECT
    CASE default
       !
       PRINT*,'mala dimension!!! (F: diametro  M: Util)'
       STOP
    END SELECT
    !
    diametro = MAXVAL(hh)
    !
    DEALLOCATE(lado,hh,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (F: diametro   M: util)!!'
       STOP
    END IF
    !
  END FUNCTION diametro
  !
  SUBROUTINE util_get_unit (iunit)
    !
    !*******************************************************************************
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer IUNIT.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5 and 6).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    IMPLICIT NONE
    !
    INTEGER i, ios, iunit
    LOGICAL lopen
    !  
    iunit = 0
    !
    DO i = 1, 99    
       IF ( i /= 5 .AND. i /= 6 ) THEN        
          INQUIRE ( unit = i, opened = lopen, iostat = ios )        
          IF ( ios == 0 ) THEN
             IF ( .NOT. lopen ) THEN
                iunit = i
                RETURN
             END IF
          END IF
       END IF
    END DO
    !
  END SUBROUTINE util_get_unit
  !
  SUBROUTINE util_medida(malla,medida)
    !
    ! calculo del area de un triangulo o del
    ! volumen de un tetraedro. Las coordenadas
    ! de los vertices son conocidas.
    !
    TYPE(mesh)                 :: malla
    REAL(kind=dp), INTENT(out) :: medida(:)
    REAL(kind=dp), ALLOCATABLE :: coord(:,:), matriz(:,:)
    INTEGER, ALLOCATABLE       :: mm(:)
    REAL(kind=dp)              :: det, factor
    INTEGER                    :: i,j,ndim,nod,nelem
    !
    ndim = 0; nod = 0; nelem = 0; factor = 0.0_dp
    !
    ndim  = malla%ndim
    nod   = malla%element%nod
    nelem = malla%nelem
    !
    ALLOCATE(coord(ndim,nod),mm(nod),matriz(ndim,ndim))
    !
    coord = 0.0_dp; mm = 0; medida = 0.0_dp
    !
    SELECT CASE(ndim)
    CASE(2) 
       factor = 0.5_dp
    CASE(3)
       factor = 1.0_dp/6.0_dp
    CASE default
       PRINT*,' dimension equivocada!! (S: util_medida, M: util)'
       STOP
    END SELECT
    !
    DO i=1,nelem
       !
       det       = 0.0_dp
       matriz    = 0.0_dp
       coord     = 0.0_dp
       mm        = 0
       det       = 0.0_dp
       !
       mm = malla%conec(:,i)
       !
       DO j=1,nod
          coord(:,j) = malla%coord(:,mm(j))
       END DO
       !
       SELECT CASE(ndim)
       CASE(2)
          !
          matriz(1,:) = (/coord(1,2) - coord(1,1), coord(1,3) - coord(1,1)/)
          matriz(2,:) = (/coord(2,2) - coord(2,1), coord(2,3) - coord(2,1)/)
          !
       CASE(3)
          !
          matriz(1,:) = (/coord(1,1)-coord(1,2), coord(1,2) - coord(1,3), coord(1,3) - coord(1,4)/)
          matriz(2,:) = (/coord(2,1)-coord(2,2), coord(2,2) - coord(2,3), coord(2,3) - coord(2,4)/)
          matriz(3,:) = (/coord(3,1)-coord(3,2), coord(3,2) - coord(3,3), coord(3,3) - coord(3,4)/)
          !
       CASE default
          PRINT*,' dimension equivocada!! (S: util_medida, M: util)'
          STOP
       END SELECT
       !
       det = determinant(matriz)
       !
       medida(i) = factor*ABS(det)
       !
    END DO
    !
    DEALLOCATE(coord,mm,matriz)
    !
  END SUBROUTINE util_medida
  !
  SUBROUTINE util_hF_nF(cara,coord,h_face,normal) !P1 2D y 3D
    !
    ! Calculo de la longitud de un lado y su normal exterior (solo P1)
    !
    INTEGER, INTENT(in)        :: cara
    REAL(kind=dp), INTENT (in) :: coord(:,:)
    REAL(kind=dp),INTENT(out)  :: h_face,normal(:)
    REAL(kind=dp),ALLOCATABLE  :: edge(:,:)
    REAL(kind=dp)              :: factor,norma
    INTEGER                    :: ndim,nedge
    !
    ndim   = 0
    nedge  = 0
    h_face = 0.0_dp
    normal = 0.0_dp
    factor = 0.0_dp
    norma  = 0.0_dp
    !
    ndim = SIZE(coord,1)
    !
    SELECT CASE(ndim)
    CASE(2)
       factor = 1.0_dp
       nedge  = 3
    CASE(3)
       factor = 0.5_dp
       nedge  = 6
    CASE default
       PRINT*,'mala dimension (S: util_hF_nF  M: util)'
       STOP
    END SELECT
    !    
    ALLOCATE(edge(ndim,nedge))
    !
    edge = 0.0_dp
    !
    SELECT CASE(ndim)
    CASE(2) !triangulo P1
       !
       edge(:,1) = coord(:,2) - coord(:,1)
       edge(:,2) = coord(:,3) - coord(:,2)
       edge(:,3) = coord(:,1) - coord(:,3)
       !
    CASE(3) !tetraedro P1
       !
       edge(:,1) = coord(:,2) - coord(:,1)
       edge(:,2) = coord(:,3) - coord(:,1)
       edge(:,3) = coord(:,4) - coord(:,1)
       edge(:,4) = coord(:,3) - coord(:,2)
       edge(:,5) = coord(:,4) - coord(:,2)
       edge(:,6) = coord(:,4) - coord(:,3)       
       !
    CASE default
       PRINT*,'elemento aun no implementado (funcion util_hF_nF)'
       STOP
    END SELECT
    !
    SELECT CASE(ndim)
    CASE(2) ! triangulo
       !
       normal = (/edge(2,cara),-edge(1,cara)/)
       !
    CASE(3) !dimension 3
       SELECT CASE(cara) !tetraedro
       CASE(1)
          !
          normal = vect_product(edge(:,1),edge(:,2))
          IF(DOT_PRODUCT(normal,edge(:,3)) > 0.0_dp) normal = -normal
          !
       CASE(2)
          !
          normal = vect_product(edge(:,1),edge(:,3))
          IF(DOT_PRODUCT(normal,edge(:,2)) > 0.0_dp) normal = -normal
          !
       CASE(3)
          !
          normal = vect_product(edge(:,2),edge(:,3))
          IF(DOT_PRODUCT(normal,edge(:,1)) > 0.0_dp) normal = -normal
          !
       CASE(4)
          !
          normal = vect_product(edge(:,4),edge(:,5))
          IF(DOT_PRODUCT(normal,-edge(:,1)) > 0.0_dp) normal = -normal
          !
       CASE default
          PRINT*,'elemento no implementado!! (S: util_hF_nF, M:util)'
          STOP
       END SELECT
       !
    CASE default
       PRINT*,'mala dimension!! (S: util_hF_nF, M: util)'
       STOP
    END SELECT
    !
    norma = DOT_PRODUCT(normal,normal)
    !
    h_face  = factor*norma
    !
    normal = normal/norma 
    !
    DEALLOCATE(edge)
    !
  END SUBROUTINE util_hF_nF
  !
  FUNCTION vect_product(uu,vv)
    !
    ! funcion que calcula el producto vectorial de
    ! dos vectores de R^3
    !
    REAL(kind=dp),INTENT(in) :: uu(:),vv(:)
    REAL(kind=dp)            :: vect_product(SIZE(uu))
    !
    vect_product = 0.0_dp
    !
    vect_product(1) =   uu(2)*vv(3) - uu(3)*vv(2) 
    vect_product(2) = -(uu(1)*vv(3) - uu(3)*vv(1))
    vect_product(3) =   uu(1)*vv(2) - uu(2)*vv(1)
    !
  END FUNCTION vect_product
  !
END MODULE util

