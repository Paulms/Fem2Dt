MODULE storage
  !
  ! Modulo relacionado con la generacion de la matriz del sistema lineal en Morse
  !
  ! (creacion de punteros, ensamble de la matriz)
  !
  USE decimal
  USE tipos
  USE quadrature
  USE funciones
  USE util
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC   storage_A, storage_f, estructura_A
  !
CONTAINS
  SUBROUTINE estructura_A(malla, A)
    TYPE(mesh)          :: malla
    TYPE(sparse)        :: A
    ! sacar assamb
    CALL pointers_crs(malla,A)
    !
  END SUBROUTINE estructura_A
  SUBROUTINE storage_A(malla,fisica,metodo,A, dt)
    TYPE(mesh)          :: malla
    TYPE(sparse)        :: A
    TYPE(param)         :: fisica
    INTEGER, INTENT(in) :: metodo
    REAL(kind=dp)       :: dt
    ! solo malla y A
    A%aa = 0
    CALL assamb_crs(malla,fisica,metodo,A, dt)
  END SUBROUTINE storage_A
  !Add ensamlado_A
  !ensamblado_A(malla, tt, deltat, metodo) ---> assamb_crs
  SUBROUTINE pointers_crs(malla,A)
    ! 
    !    OBJETIVO: esta subrutina crea los tableros row y column
    !              necesarios para el almacenamiento tipo MORSE
    !              de una matriz no simetrica. Los nodos son los
    !              vertices de la malla  
    ! 
    !    IN:
    !      mesh : la malla de elementos finitos
    !
    !    OUT:
    !       A%nzero       ... numero de elementos no nulos de la matriz
    !       A%row(i)      ... ver def. de formato crs
    !       A%column(j)   ... ver def. de formato crs
    !
    TYPE(mesh)             :: malla
    TYPE(sparse)           :: A
    INTEGER,ALLOCATABLE    :: nvecino(:),vecino(:,:),equation(:)
    INTEGER                :: i,j,k,ii,jj,dof,nod,nodof,neq,nnode,nelem,ierr
    !
    nod = 0; dof = 0; nodof = 0; neq = 0; nnode = 0; nelem = 0
    !
    nod     = malla%element%nod
    dof     = malla%element%dof
    nnode   = malla%nnode
    nelem   = malla%nelem
    nodof   = nod*dof                  ! numero de grados de libertad por elemento
    neq     = dof*nnode                ! numero total de ecuaciones y/o incognitas
    ierr    = 0
    !
    A%nzero = 0; A%n = 0
    !
    A%n = neq
    !
    ALLOCATE(A%row(neq+1),nvecino(neq),equation(nodof),vecino(neq,800)) ! 800 = numero maximo de ecuaciones vecinas
    !
    A%row   = 0
    nvecino = 0
    vecino  = 0
    !
    !     creacion de los arreglos:
    !
    !     nvecino(i) --> numero de ecuaciones vecinas de la ecuacion i
    !     vecino(i,j), j =1,...,nvecino(i)  --> numero de la j-esima ecuacion 
    !                                           vecina de la i-esima ecuacion
    !        
    DO  k=1,nelem
       !
       equation = 0
       !
       DO ii=1,nod
          DO jj=1,dof
             equation(dof*(ii-1)+jj)=dof*(malla%conec(ii,k) -1)+jj
          END DO
       END DO
       !
       DO  i=1,nodof
          inner: DO  j=1,nodof
             IF(j==i) CYCLE inner
             IF (nvecino(equation(i))== 0) THEN
                nvecino(equation(i)) = nvecino(equation(i))+1 
                vecino(equation(i),nvecino(equation(i))) = equation(j) 
                CYCLE inner
             ENDIF
             !
             DO jj=1,nvecino(equation(i))
                IF(equation(j)==vecino(equation(i),jj)) CYCLE inner
             END DO
             !
             nvecino(equation(i))                     = nvecino(equation(i)) + 1
             vecino(equation(i),nvecino(equation(i))) = equation(j)
          END DO inner
       END DO
    END DO
    !
    !  se ordena vecino con algoritmo de orden
    !
    DO  i = 1,neq
       vecino(i,nvecino(i)+1) = i
       nvecino(i)             = nvecino(i)+1
       !
       CALL I_inssor(vecino(i,1:nvecino(i)))
       ! 
    END DO
    !
    !  creacion de A%row
    !       
    A%row(1) = 1
    !
    DO  i=1,neq 
       !
       A%row(i+1) = A%row(i) + nvecino(i)
       !
    END DO
    !
    A%nzero = A%row(neq+1) - 1 ! numero de elementos no nulos
    !
    ALLOCATE(A%column(A%nzero),A%aa(A%nzero))
    !
    A%column = 0; A%aa = 0.0_dp
    !
    !  creacion de A%column
    !
    DO i=1,neq 
       !
       DO j = 1,nvecino(i)
          A%column(A%row(i)+j-1) = vecino(i,j)
       END DO
       A%column(A%row(i+1)-1)    = vecino(i,nvecino(i))
    END DO
    !
    DEALLOCATE(nvecino,vecino,equation,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: pointers_crs M:storage)!!'
       STOP
    END IF
    !
  END SUBROUTINE pointers_crs
  !
  SUBROUTINE m_assam(a_local,mm,nn,A)
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                      !
    !                 Ensamblado morse de la matriz elemental              !
    !                                                                      !
    !        VARIABLES:                                                    !
    !                                                                      !
    !          ENTRADA:                                                    !
    !            a_local   :  Matriz elemental                             !
    !            row,column:  Punteros                                     !
    !            mm        :  Numeros de los vertices del elemento         !
    !                                                                      !
    !          ENTRADA/SALIDA:                                             !
    !            aa        :  Vector que contiene la matriz del sistema    !
    !                                                                      !
    !                                                                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE 
    TYPE(sparse)             :: A
    INTEGER,INTENT(in)       :: mm(:),nn(:)
    REAL(kind=dp),INTENT(in) :: a_local(:,:)
    INTEGER                  :: nodof,i,i1,j,j1,l 
    !
    nodof = SIZE(mm)
    !
    DO  i=1,nodof
       !
       i1=mm(i)  ! m(i) = i-esimo numero de ecuacion
       !
       DO  j=1,nodof
          !
          j1=nn(j)  ! m(j) = j-esimo numero de ecuacion
          !
          inner: DO  l=A%row(i1),A%row(i1+1)-1
             IF(A%column(l)==j1) THEN
                A%aa(l) = A%aa(l) + a_local(i,j)
                EXIT inner
             END IF
             !
          END DO inner
          !
       END DO
       !
    END DO
    !                  
  END SUBROUTINE m_assam
  !
  SUBROUTINE m_local(malla,points,weights,coord,nu,velo,sigma,metodo,kt, dt)
    !
    TYPE(mesh)                :: malla
    REAL(kind=dp),INTENT(in)  :: points(:,:),weights(:),coord(:,:),nu,velo(:),sigma
    REAL(kind=dp),INTENT(out) :: kt(:,:)
    INTEGER, INTENT(in)       :: metodo
    !
    INTEGER                   :: i,j,nod,npi,ndim,ierr
    REAL(kind=dp)             :: det, ht, tau,cuociente,norma, dt
    REAL(kind=dp),ALLOCATABLE :: DeP(:,:),Bt(:,:),DPref(:,:), Pref(:,:),difusion(:,:),adveccion(:,:),&
                                 reaccion(:,:),funvel(:,:),gls(:,:),vela(:,:),Masa(:,:),Pv(:,:),x_bar(:,:),&
                                 x_node(:,:),aux(:,:)
    !
    nod = 0; npi = 0; ndim = 0; ierr = 0
    !
    nod   = malla%element%nod
    npi   = malla%element%npi
    ndim  = malla%ndim
    !
    ALLOCATE(DeP(ndim,nod), Bt(ndim,ndim),DPref(ndim,nod),Pref(1,nod),difusion(nod,nod),adveccion(nod,nod),reaccion(nod,nod),&
         gls(nod,nod),vela(1,ndim),funvel(nod,ndim),Masa(nod,nod),Pv(ndim,ndim*nod),x_bar(ndim,1),&
         x_node(ndim*nod,1),aux(ndim,nod))
    !
    kt = 0.0_dp
    !
    difusion  = 0.0_dp
    adveccion = 0.0_dp
    reaccion  = 0.0_dp
    vela      = 0.0_dp
    funvel    = 0.0_dp
    gls       = 0.0_dp
    ht        = 0.0_dp
    tau       = 0.0_dp
    cuociente = 0.0_dp
    x_bar     = 0.0_dp
    x_node    = 0.0_dp
    norma     = 0.0_dp
    !    
    SELECT CASE(ndim)
    CASE(2)
       cuociente = 1.0_dp/3.0_dp
    CASE(3)
       cuociente = 0.25_dp
    CASE default
       PRINT*,'Mala dimension!! (S: m_local, M:storge)'
       STOP
    END SELECT
    !
    x_bar(:,1) = cuociente*SUM(coord,2)
    !
    DO j=1,ndim
       x_node((j-1)*nod+1:j*nod,1) = coord(j,:) 
    END DO
    !
    ! calculo de h local
    !
    !ht  = hvf(coord,velo)
    ht = diametro(coord)
    !
    IF(metodo==2) THEN
       tau = param_adr(nu,velo,sigma,ht)
    ELSEIF(metodo==3) THEN
       tau = param_lps(nu,velo,sigma,ht)
    END IF
    !
    vela(1,:) = velo
    norma     = SQRT(DOT_PRODUCT(velo,velo))
    !
    DO i=1,npi
       !
       Bt    = 0.0_dp
       det   = 0.0_dp
       DeP   = 0.0_dp
       DPref = 0.0_dp
       Pref  = 0.0_dp
       Pv    = 0.0_dp
       Masa  = 0.0_dp
       aux   = 0.0_dp
       !
       gls    = 0.0_dp
       !
       CALL shape_der(DPref, points,i)
       !
       CALL shape_fun(Pref(1,:),points,i)
       !
       DO j=1,ndim
          Pv(j,(j-1)*nod+1:j*nod) = Pref(1,:)
       END DO
       !
       funvel= MATMUL(TRANSPOSE(Pref),vela)
       Masa = MATMUL(TRANSPOSE(Pref),Pref)
       !
       Bt   = MATMUL(DPref,TRANSPOSE(coord))
       det  = determinant(Bt)
       det  = ABS(det)
       !
       CALL invert(Bt)
       DeP = MATMUL(Bt,DPref)
       !
       difusion  =  nu*MATMUL(TRANSPOSE(DeP),DeP)*det*weights(i)
       !
       adveccion =  MATMUL(funvel,DeP)*det*weights(i)
       !
       reaccion  =  (sigma + 1/dt)*Masa*det*weights(i)
       !
       IF(metodo == 2) THEN
          !
          ! Metodo unusual tipo SUPG  (cf, Franca & Valentin 2000)
          !
          gls = gls - sigma**2*tau*Masa*det*weights(i)
          gls = gls + sigma*tau*TRANSPOSE(MATMUL(funvel,DeP))*det*weights(i)
          gls = gls - sigma*tau*MATMUL(funvel,DeP)*det*weights(i)
          gls = gls + tau*MATMUL(TRANSPOSE(MATMUL(vela,DeP)),MATMUL(vela,DeP))*det*weights(i)
          !
          kt = kt + gls
          !
       ELSEIF(metodo==3) THEN
          !
          ! Metodo LPS modificado: P1 2D/3D  (cf. Araya & Mu~noz & Valentin 2011)
          !
          aux = MATMUL(MATMUL(Pv,x_node)-x_bar,MATMUL(vela,DeP))
          !
          gls =  tau*MATMUL(TRANSPOSE(aux),aux)*det*weights(i)
          !
          kt = kt + gls
          !
!!$       ELSEIF(metodo==3) THEN
!!$          !
!!$          ! Metodo LPS modificado: P1 2D/3D  (cf. Araya & Mu~noz & Valentin 2011)
!!$          !
!!$          aux = MATMUL(MATMUL(Pv,x_node)-x_bar,MATMUL(vela,DeP))
!!$          !aux = (1.0_dp/norma)*DOT_PRODUCT(velo,MATMUL(Pv,x_node)-x_bar)*MATMUL(vela,DeP)
!!$          !
!!$          gls =  tau*MATMUL(TRANSPOSE(aux),aux)*det*weights(i)
!!$          !
!!$          kt = kt + gls
          !
       END IF
       !
       kt = kt + difusion + adveccion + reaccion
       !
    END DO
    !
    DEALLOCATE(DeP,Bt,DPref,Pref,difusion,adveccion,reaccion,gls,vela,Pv,x_bar,x_node,aux,&
               funvel, Masa,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: m_local  M:storage)!!'
       STOP
    END IF
      !
  END SUBROUTINE m_local
  !
  SUBROUTINE assamb_crs(malla,fisica,metodo,A, dt)
    !
    TYPE(mesh)                 :: malla
    TYPE(sparse)               :: A
    TYPE(param)                :: fisica
    INTEGER, INTENT(in)        :: metodo
    REAL(kind=dp),ALLOCATABLE  :: points(:,:),coor(:,:),weights(:),kt(:,:)
    INTEGER, ALLOCATABLE       :: mm(:)
    INTEGER                    :: npi,ndim,dof,nod,nodof,ind,i,j,ref_ele,ierr
    REAL(kind=dp)              :: dt
    !
    npi = 0; ndim = 0; dof = 0; nod = 0; nodof = 0; ref_ele = 0; ierr = 0
    !
    npi   = malla%element%npi
    dof   = malla%element%dof
    nod   = malla%element%nod
    ndim  = malla%ndim
    nodof = nod*dof 
    !
    ALLOCATE(points(npi,ndim),weights(npi),coor(ndim,nod),kt(nodof,nodof),mm(nod))
    !
    CALL sample(malla%element%tipo,points,weights)
    !
    DO i=1,malla%nelem
       !
       mm      = 0
       ref_ele = 0
       ind     = 0
       coor    = 0.0_dp
       kt      = 0.0_dp
       !
       mm      = malla%conec(:,i)
       ref_ele = malla%ref_elem(i)
       !
       DO j=1,fisica%ndominios
          IF(fisica%ref_dominio(j) == ref_ele) THEN
             ind = j
             EXIT
          END IF
       END DO
       !
       DO j=1,nod
          coor(:,j) = malla%coord(:,mm(j))
       END DO
       ! ahora depende de dt y tt
       CALL m_local(malla,points,weights,coor,fisica%nu(ind),fisica%aa(ind,:),fisica%sigma(ind),metodo,kt, dt)
       !
       CALL m_assam(kt,mm,mm,A) 
       !
    END DO
    !
     DEALLOCATE(points,weights,coor,kt,mm,STAT=ierr)
      IF(ierr/=0) THEN
         PRINT*,'Problema deallocando variables (S: assamb_crs   M:storage)!!'
         STOP
      END IF
      !
  END SUBROUTINE assamb_crs
  !
  SUBROUTINE storage_f(malla,fisica,metodo,uh, ua, tt, dt)
    !
    TYPE(mesh)                :: malla
    TYPE(param)               :: fisica
    INTEGER, INTENT(in)       :: metodo
    REAL(kind=dp)             :: ua(:), tt, dt
    REAL(kind=dp),ALLOCATABLE :: uh(:)
    !
    ALLOCATE(uh(malla%nnode))
    !
    CALL rhs_pointers(malla,fisica,metodo,uh, ua, tt, dt)
    !
  END SUBROUTINE storage_f
  !
  SUBROUTINE rhs_pointers(malla,fisica,metodo,uh, ua, tt, dt)
    !
    TYPE(mesh)                 :: malla
    TYPE(param)                :: fisica
    INTEGER,INTENT(in)         :: metodo
    REAL(kind=dp)              :: uh(:), ua(:), tt, dt
    REAL(kind=dp),ALLOCATABLE  :: points(:,:),coor(:,:),weights(:),ft(:)
    INTEGER, ALLOCATABLE       :: mm(:)
    INTEGER                    :: npi,ndim,dof,nod,nodof,i,j,ind,ref_ele,ierr
    !
    npi = 0; ndim = 0; dof = 0; nod = 0; nodof = 0; ref_ele = 0; ierr = 0
    !
    npi   = malla%element%npi
    dof   = malla%element%dof
    nod   = malla%element%nod
    ndim  = malla%ndim
    nodof = nod*dof 
    !
    ALLOCATE(points(npi,ndim),weights(npi),coor(ndim,nod),ft(nodof),mm(nod))
    !
    uh = 0.0_dp
    !
    CALL sample(malla%element%tipo,points,weights)
    !
    DO i =1,malla%nelem
       !
       mm      = 0
       coor    = 0.0_dp
       ft      = 0.0_dp
       ref_ele = 0
       ind     = 0
       !
       mm      = malla%conec(:,i)
       ref_ele = malla%ref_elem(i)
       !
       DO j=1,nod
          coor(:,j) = malla%coord(:,mm(j))
       END DO
       !
       DO j=1,fisica%ndominios
          IF(fisica%ref_dominio(j) == ref_ele) THEN
             ind = j
             EXIT
          END IF
       END DO
       !
       CALL rhs_local(malla,points,weights,coor,ref_ele,fisica%nu(ind),fisica%aa(ind,:),fisica%sigma(ind),metodo,ft, ua(mm), tt, dt)
       !
       CALL rhs_assam(ft,mm,uh)
       !
    END DO
    !
    DEALLOCATE(points,weights,coor,ft,mm,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: rhs_pointers   M:storage)!!'
       STOP
    END IF
    !
  END SUBROUTINE rhs_pointers
  !
  SUBROUTINE rhs_local(malla,points,weights,coord,ref_ele,nu,velo,sigma,metodo,ft, ua_loc, tt, dt)
    !
    TYPE(mesh)                :: malla
    REAL(kind=dp),INTENT(in)  :: points(:,:),weights(:),coord(:,:),nu,velo(:),sigma, ua_loc(:), tt, dt
    INTEGER, INTENT(in)       :: metodo,ref_ele
    REAL(kind=dp),INTENT(out) :: ft(:)
    INTEGER                   :: i,j,npi,ndim,nod,dof,nodof,ierr
    REAL(kind=dp),ALLOCATABLE :: Bt(:,:),DPref(:,:),DeP(:,:),Pref(:,:),masa(:,:),vela(:,:),gls_rhs(:),f_node(:),&
                                 Pv(:,:),x_bar(:),x_node(:)
    REAL(kind=dp)             :: det,tau,ht,cuociente,f_bar,aux,norma
    !
    nod = 0; npi = 0; ndim = 0; dof = 0; nodof = 0; ierr = 0
    !
    nod   = malla%element%nod
    npi   = malla%element%npi
    ndim  = malla%ndim
    dof   = malla%element%dof
    nodof = nod*dof
    !
    ALLOCATE(DeP(ndim,nod), Bt(ndim,ndim), DPref(ndim,nod), Pref(1,nod),Masa(nodof,nodof),gls_rhs(nodof),vela(1,ndim),&
             f_node(nodof), Pv(ndim,ndim*nod),x_bar(ndim),x_node(ndim*nod))
    !
    ft = 0.0_dp; ht = 0.0_dp; f_node = 0.0_dp; cuociente = 0.0_dp; x_bar = 0.0_dp; x_node = 0.0_dp; f_bar = 0.0_dp; norma = 0.0_dp
    !
    vela(1,:) = velo
    norma     = SQRT(DOT_PRODUCT(velo,velo))
    !
    SELECT CASE(ndim)
    CASE(2)
       cuociente = 1.0_dp/3.0_dp
    CASE(3)
       cuociente = 0.25_dp
    CASE default
       PRINT*,'Mala dimensio!! (S: m_local, M:storge)'
       STOP
    END SELECT
    !
    ! x_bar =  baricentro del elemento T
    !
    x_bar = cuociente*SUM(coord,2)
    !
    ! x_node = valores en los nodos del vector X
    !
    DO j=1,ndim
       x_node((j-1)*nod+1:j*nod) = coord(j,:) 
    END DO
    !
    ! f_node = valores de f en los nodos:
    !
    DO j=1,nod
       f_node(j) = func(coord(:,j),tt,nu,velo,sigma,ref_ele)
    END DO
    ! Add 1/dt*u^n
    f_node = f_node + 1/dt*ua_loc
    !
    !  f_bar = \frac{1}{|T|} \int_T f,   donde f esta proyectada en P1
    !
    f_bar = cuociente*SUM(f_node)
    !
    !ht  = hvf(coord,velo)
    ht = diametro(coord)
    IF(metodo==2) THEN
       tau = param_adr(nu,velo,sigma,ht)
    ELSEIF(metodo==3) THEN
       tau = param_lps(nu,velo,sigma,ht)
    END IF
    !
    DO i=1,npi
       !
       Bt      = 0.0_dp
       det     = 0.0_dp
       DPref   = 0.0_dp
       DeP     = 0.0_dp
       Pref    = 0.0_dp
       Pv      = 0.0_dp
       Masa    = 0.0_dp
       gls_rhs = 0.0_dp
       aux     = 0.0_dp
       !
       CALL shape_der(DPref, points,i)
       CALL shape_fun(Pref(1,:),points,i)
       !
       DO j=1,ndim
          Pv(j,(j-1)*nod+1:j*nod) = Pref(1,:)
       END DO
       !
       Bt   = MATMUL(DPref,TRANSPOSE(coord))
       det  = determinant(Bt)
       det  = ABS(det)
       !
       CALL invert(Bt)
       DeP = MATMUL(Bt,DPref)
       !
       Masa = MATMUL(TRANSPOSE(Pref),Pref)
       !
       IF(metodo==2) THEN
          gls_rhs =         - tau*sigma*MATMUL(Masa,f_node)*det*weights(i)
          gls_rhs = gls_rhs + tau*MATMUL(TRANSPOSE(MATMUL(vela,DeP)),MATMUL(Pref,f_node))*det*weights(i)
          !
          ft      = ft + gls_rhs
          !
       ELSEIF(metodo==3) THEN
          !
          !aux     = (f_bar/(norma**2))*(DOT_PRODUCT(velo,MATMUL(Pv,x_node)-x_bar))**2
          !
          !gls_rhs =  tau*aux*MATMUL(TRANSPOSE(DeP),velo)*det*weights(i)
          !
          aux     = f_bar*DOT_PRODUCT(MATMUL(Pv,x_node)-x_bar,MATMUL(Pv,x_node)-x_bar)
          !
          gls_rhs =  tau*aux*MATMUL(TRANSPOSE(DeP),velo)*det*weights(i)
          !
          ft = ft + gls_rhs
          !
       END IF
       !
       ft  = ft + MATMUL(Masa,f_node)*det*weights(i)
       !
    END DO
    !
    DEALLOCATE(DeP,Bt,DPref,Pref,Masa,gls_rhs,vela,f_node,x_bar,x_node,Pv,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: rhs_local  M:storage)!!'
       STOP
    END IF
    !
  END SUBROUTINE rhs_local
  !
  SUBROUTINE rhs_assam(ft,mm,uh)
    !
    ! Ensamble del lado derecho del sistema lineal (rhs)
    !
    REAL(kind=dp), INTENT(in)  :: ft(:)
    INTEGER,INTENT(in)         :: mm(:)
    REAL(kind=dp), INTENT(out) :: uh(:)   
    !
    uh(mm) = uh(mm) + ft
    !
  END SUBROUTINE rhs_assam
  !
  SUBROUTINE i_inssor (xdont)
    !
    !  sorts xdont into increasing order (insertion sort)
    ! __________________________________________________________
    !  this subroutine uses insertion sort. it does not use any
    !  work array and is faster when xdont is of very small size
    !  (< 20), or already almost sorted, but worst case behavior
    !  can happen fairly probably (initially inverse sorted).
    !  in many cases, the quicksort or merge sort method is faster.
    !  michel olagnon - apr. 2000
    ! __________________________________________________________
    ! __________________________________________________________
    ! __________________________________________________________
    INTEGER, INTENT (inout)  :: xdont(:)
    ! __________________________________________________________
    INTEGER :: xwrk, xmin
    !
    ! __________________________________________________________
    !
    INTEGER :: icrs, idcr, ndon
    !
    ndon = SIZE (xdont)
    !
    ! we first bring the minimum to the first location in the array.
    ! that way, we will have a "guard", and when looking for the
    ! right place to insert a value, no loop test is necessary.
    !
    IF (xdont (1) < xdont (ndon)) THEN
       xmin = xdont (1)
    ELSE
       xmin = xdont (ndon)
       xdont (ndon) = xdont (1)
    ENDIF
    DO idcr = ndon-1, 2, -1
       xwrk = xdont(idcr)
       IF (xwrk < xmin) THEN
          xdont (idcr) = xmin
          xmin = xwrk
       END IF
    END DO
    xdont (1) = xmin
    !
    ! the first value is now the minimum
    ! loop over the array, and when a value is smaller than
    ! the previous one, loop down to insert it at its right place.
    !
    DO icrs = 3, ndon
       xwrk = xdont (icrs)
       idcr = icrs - 1
       IF (xwrk < xdont(idcr)) THEN
          xdont (icrs) = xdont (idcr)
          idcr = idcr - 1
          DO
             IF (xwrk >= xdont(idcr)) EXIT
             xdont (idcr+1) = xdont (idcr)
             idcr = idcr - 1
          END DO
          xdont (idcr+1) = xwrk
       END IF
    END DO
    !
  END SUBROUTINE i_inssor
  !
END MODULE storage











