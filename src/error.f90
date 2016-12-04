MODULE error
  !
  ! Modulo relacionado con el calculo del error a priori y a posteriori
  !
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
  PUBLIC   estimacion 
  !
CONTAINS
  !  
  SUBROUTINE estimacion(malla,fisica,frontera,uh,aprox,tf)
    TYPE(mesh)                 :: malla 
    TYPE(param)                :: fisica
    TYPE(bc)                   :: frontera
    INTEGER, INTENT(in)        :: aprox
    REAL(kind=dp), INTENT(in)  :: uh(:)
    REAL(kind=dp)              :: h_max, norma_l2, snorma_h1, norma_h1, eta
    INTEGER                    :: iunit
    REAL(kind=dp)              :: tf
    tf = 0.0_dp
    !
    CALL util_get_unit(iunit)
    !
    OPEN(unit=iunit,file='errores.dat',status='replace',action='write')
    !
    h_max = 0.0_dp; norma_l2 = 0.0_dp; snorma_h1 = 0.0_dp; norma_h1 = 0.0_dp; eta = 0.0_dp
    !
    SELECT CASE(aprox)
    CASE(1)          ! Solo arror a priori
       !    
       CALL h_malla(malla,h_max)
       !
       CALL a_priori(malla,fisica,uh,norma_l2,snorma_h1,norma_h1,tf)
       !
       WRITE(iunit,*)' h         ||u-uh||_{0,\Omega}               |u-uh|_{1,\Omega}               ||u-uh||_{\Omega}'
       WRITE(iunit,*) h_max, SQRT(norma_l2), SQRT(snorma_h1), SQRT(norma_h1)
       !
    CASE(2)          ! Solo error a posteriori
       !
       CALL a_posteriori(malla,fisica,frontera,uh,eta,tf)
       !
    CASE(3)          ! A posteriori y a priori
       !
       CALL h_malla(malla,h_max)
       !
       CALL a_priori(malla,fisica,uh,norma_l2,snorma_h1,norma_h1,tf)
       !
       CALL a_posteriori(malla,fisica,frontera,uh,eta,tf)  
       WRITE(iunit,*)' h ||u-uh||_{0,\Omega}  |u-uh|_{1,\Omega}  ||u-uh||_{\Omega}   eta   efectividad   '
       WRITE(iunit,*) h_max, SQRT(norma_l2), SQRT(snorma_h1), SQRT(norma_h1), eta, eta/SQRT(norma_h1)
       !
    CASE default
       PRINT*,'problema de calculo de errores!! (S: estimacion, M: error)'
       STOP
    END SELECT
    !
    CLOSE(iunit)
    !
  END SUBROUTINE estimacion
  !
  SUBROUTINE a_priori(malla,fisica,uh,norma_l2,snorma_h1,norma_h1,tf)
    
    TYPE(mesh)                 :: malla 
    TYPE(param)                :: fisica
    REAL(kind=dp), INTENT(in)  :: uh(:)
    REAL(kind=dp), INTENT(out) :: norma_l2, snorma_h1, norma_h1
    REAL(kind=dp)              :: nl2,snh1,nh1
    REAL(kind=dp),ALLOCATABLE  :: points(:,:),coor(:,:),weights(:),uh_node(:)
    INTEGER, ALLOCATABLE       :: mm(:)
    INTEGER                    :: npi,ndim,dof,nod,nodof,i,j,ind,ref_ele,ierr
    REAL(kind=dp)              :: tf
    !
    npi = 0; ndim = 0; dof = 0; nod = 0; nodof = 0; ref_ele = 0; ierr = 0
    tf = 0.0_dp
    !
    npi   = malla%element%npi
    dof   = malla%element%dof
    nod   = malla%element%nod
    ndim  = malla%ndim
    nodof = nod*dof 
    !
    ALLOCATE(points(npi,ndim),weights(npi),coor(ndim,nod),uh_node(nod),mm(nod))
    !
    points = 0.0_dp; weights = 0.0_dp
    !
    norma_l2 = 0.0_dp; snorma_h1 = 0.0_dp; norma_h1 = 0.0_dp
    !
    CALL sample(malla%element%tipo,points,weights)
    !
    DO i =1,malla%nelem
       !
       mm      = 0
       coor    = 0.0_dp
       uh_node = 0.0_dp
       ref_ele = 0
       ind     = 0
       !
       nl2     = 0.0_dp
       snh1    = 0.0_dp
       nh1     = 0.0_dp
       !
       mm      = malla%conec(:,i)
       ref_ele = malla%ref_elem(i)
       uh_node = uh(mm)
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
       CALL nl2_local(points,weights,coor,uh_node,fisica%nu(ind),ref_ele,nl2,tf)
       !
       norma_l2 = norma_l2 + nl2
       !
       CALL snh1_local(points,weights,coor,uh_node,fisica%nu(ind),ref_ele,snh1,tf)     
       !
       snorma_h1 = snorma_h1 + snh1
       !
       norma_h1 = norma_h1 + nl2 + fisica%nu(ind)*snh1
       !
    END DO
    !
    DEALLOCATE(points,weights,coor,uh_node,mm,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: a_priori   M: error)!!'
       STOP
    END IF
    !
  END SUBROUTINE a_priori
  !
  SUBROUTINE nl2_local(points,weights,coord,uh_node,nu,ref_ele,nl2,tf)
    !
    REAL(kind=dp),INTENT(in)  :: points(:,:),weights(:),coord(:,:),uh_node(:),nu
    INTEGER, INTENT(in)       :: ref_ele
    REAL(kind=dp),INTENT(out) :: nl2
    INTEGER                   :: i,j,ierr,npi,ndim,nod
    REAL(kind=dp),ALLOCATABLE :: Bt(:,:),DPref(:,:),Pref(:),afin(:)
    REAL(kind=dp)             :: det,tf
    !
    nod = 0; npi = 0; ndim = 0; ierr = 0;tf=0.0_dp
    !
    nod   = SIZE(uh_node)
    npi   = SIZE(points,1)
    ndim  = SIZE(points,2)
    !
    ALLOCATE(afin(ndim),Bt(ndim,ndim), DPref(ndim,nod), Pref(nod))
    !
    nl2 = 0.0_dp
    !
    DO i=1,npi
       !
       Bt      = 0.0_dp
       det     = 0.0_dp
       DPref   = 0.0_dp
       Pref    = 0.0_dp
       afin    = 0.0_dp
       !
       CALL shape_der(DPref, points,i)
       CALL shape_fun(Pref,  points,i)
       !
       Bt   = MATMUL(DPref,TRANSPOSE(coord))
       det  = determinant(Bt)
       det = ABS(det)
       !
       afin = MATMUL(coord,Pref)
       !
       nl2  = nl2 + ((u_ex(afin,nu,ref_ele,tf) - DOT_PRODUCT(Pref,uh_node))**2)*det*weights(i)
       !
    END DO
    !
    DEALLOCATE(Bt,DPref,Pref,afin,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: nl2_local  M: error)!!'
       STOP
    END IF
    !
  END SUBROUTINE nl2_local
  !
  SUBROUTINE snh1_local(points,weights,coord,uh_node,nu,ref_ele,snh1,tf) 
    !
    REAL(kind=dp),INTENT(in)  :: points(:,:),weights(:),coord(:,:),uh_node(:),nu
    INTEGER,INTENT(in)        :: ref_ele
    REAL(kind=dp),INTENT(out) :: snh1
    INTEGER                   :: i,j,ierr,npi,ndim,nod
    REAL(kind=dp),ALLOCATABLE :: Bt(:,:),DPref(:,:),Pref(:),DeP(:,:),afin(:)
    REAL(kind=dp)             :: det, tf
    !
    nod = 0; npi = 0; ndim = 0; ierr = 0;tf = 0.0_dp
    !
    nod   = SIZE(uh_node)
    npi   = SIZE(points,1)
    ndim  = SIZE(points,2)
    !
    ALLOCATE(Bt(ndim,ndim), Pref(nod),DPref(ndim,nod), DeP(ndim,nod),afin(ndim))
    !
    snh1 = 0.0_dp
    !
    DO i=1,npi
       !
       Bt     = 0.0_dp
       det    = 0.0_dp
       Pref   = 0.0_dp
       DPref  = 0.0_dp
       DeP    = 0.0_dp
       afin   = 0.0_dp
       !
       CALL shape_fun(Pref,points,i)
       CALL shape_der(DPref, points,i)
       !
       Bt   = MATMUL(DPref,TRANSPOSE(coord))
       det  = determinant(Bt)
       det  = ABS(det)
       !
       afin = MATMUL(coord,Pref)
       !       
       CALL invert(Bt)
       DeP = MATMUL(Bt,DPref)
       !
       snh1  = snh1 + DOT_PRODUCT(grad_ex(afin,nu,ref_ele,tf) - &
       MATMUL(DeP,uh_node),grad_ex(afin,nu,ref_ele,tf) - &
       MATMUL(DeP,uh_node))*det*weights(i)
       !
    END DO
    !
    DEALLOCATE(Bt,Pref,DPref,DeP,afin, STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'Problema deallocando variables (S: snh1_local  M: error)!!'
       STOP
    END IF
    !
  END SUBROUTINE snh1_local
    !
  SUBROUTINE a_posteriori(malla,fisica,frontera,uh,eta,tf)
    !
    TYPE(mesh)                  :: malla
    TYPE(param)                 :: fisica
    TYPE(bc)                    :: frontera
    REAL(kind=dp),INTENT(in)    :: uh(:)
    REAL(kind=dp),ALLOCATABLE   :: medida(:),eta_R(:)
    INTEGER                     :: iunit,i
    REAL(kind=dp)               :: eta, eta_max, theta, factor, tol, tf      
    !
    ALLOCATE(eta_R(malla%nelem),medida(malla%nelem))
    !
    eta = 0.0_dp; eta_R = 0.0_dp; iunit = 0; medida = 0.0_dp; eta_max = 0.0_dp; tol = 0.0_dp
    tf = 0.0_dp
    !
    theta  = 0.25_dp
    !
    factor = 0.1666666666_dp
    !
    tol = 1.0e-09_dp
    !
    CALL util_medida(malla,medida)
    !
    CALL residual(malla,fisica,frontera,uh,eta_R,tf)
    !
    eta = SQRT(SUM(eta_R))
    !
    WRITE(*,*)'estimador residual (eta) = ', eta
    !
    eta_R = SQRT(eta_R)
    !
    eta_max = MAXVAL(eta_R)
    !
    WHERE(eta_R >= theta*eta_max .AND. medida >= tol)
       medida = factor*medida
    ELSEWHERE
       medida = 0.0_dp
    END WHERE
    !
    CALL util_get_unit(iunit)
    !
    OPEN(iunit,file='geometria.area',status='replace', action='write' )
    !
    WRITE(iunit,*) malla%nelem
    DO i=1,malla%nelem
       WRITE(iunit,*) i,' ',medida(i)
    END DO
    !
    CLOSE(iunit)
    !
    DEALLOCATE(medida,eta_R)
    !
  END SUBROUTINE a_posteriori
  !
  SUBROUTINE residual(malla,fisica,frontera,uh,eta_R, tf)
    !
    ! Estimador residual para el pb de transporte   (2D por ahora)
    !
    TYPE(mesh)                :: malla
    TYPE(param)               :: fisica
    TYPE(bc)                  :: frontera
    REAL(kind=dp), INTENT(in) :: uh(:)
    REAL(kind=dp), INTENT(out):: eta_R(:)
    REAL(kind=dp)             :: norma_residuo,norma_salto, tf
    INTEGER                   :: i,nelem
    !
    nelem = 0; eta_R = 0.0_dp; tf = 0.0_dp
    !
    nelem = malla%nelem
    !
    DO i=1,nelem
       !
       norma_residuo = 0.0_dp; norma_salto = 0.0_dp
       !
       CALL n_residuo(malla,fisica,i,uh,norma_residuo,tf)
       !
       CALL n_salto(malla,fisica,frontera,i,uh,norma_salto)
       !
       eta_R(i) = norma_residuo + norma_salto
       !
    END DO
    !
  END SUBROUTINE residual
  !
  SUBROUTINE n_residuo(malla,fisica,ele,uh,norma_residuo,tf)
    !
    TYPE(mesh)                :: malla
    TYPE(param)               :: fisica
    REAL(kind=dp), INTENT(in) :: uh(:)
    INTEGER, INTENT(in)       :: ele
    REAL(kind=dp)             :: norma_residuo,h_elem,theta_T,nu,sigma,det,norma
    INTEGER                   :: i,j,nod,ndim,npi,ref_ele, ind
    INTEGER, ALLOCATABLE      :: mm(:)
    REAL(kind=dp), ALLOCATABLE:: coord(:,:), Bt(:,:), DPref(:,:), DeP(:,:), f_node(:), uh_node(:),&
                                 points(:,:),weights(:), velo(:), Pref(:),afin(:)
    REAL(kind=dp)             :: tf
    !
    norma_residuo = 0.0_dp; nod = 0; ndim = 0; ref_ele = 0; nu = 0.0_dp; sigma = 0.0_dp; theta_T = 0.0_dp
    h_elem = 0.0_dp; npi = 0; tf = 0.0_dp
    !
    nod  = malla%element%nod
    ndim = malla%ndim
    npi  = malla%element%npi
    !
    ALLOCATE(mm(nod),coord(ndim,nod),Bt(ndim,ndim),DPref(ndim,nod),DeP(ndim,nod),f_node(nod),uh_node(nod),&
             Pref(nod),points(npi,ndim),weights(npi),velo(ndim),afin(ndim))
    !
    mm = 0; coord = 0.0_dp; f_node = 0.0_dp; uh_node = 0.0_dp
    velo = 0.0_dp; points = 0.0_dp; weights = 0.0_dp
    !
    mm = malla%conec(:,ele)
    !
    DO j=1,nod
       coord(:,j)  = malla%coord(:,mm(j))
    END DO
    !
    ref_ele = malla%ref_elem(ele)
    !
    DO j=1,fisica%ndominios
       !
       IF(fisica%ref_dominio(j) == ref_ele)  THEN
          ind  = j
          EXIT
       END IF
       !
    END DO
    !
    nu    = fisica%nu(ind)
    velo  = fisica%aa(ind,:)
    sigma = fisica%sigma(ind)
    !
!!$    DO j=1,nod
!!$       f_node(j) = func(coord(:,j),nu,velo,sigma,ref_ele)
!!$    END DO
    !
    h_elem = diametro(coord)
    !
    theta_T = MIN(1.0_dp, h_elem/SQRT(nu))
    !
    uh_node  = uh(mm)
    !
    CALL sample(malla%element%tipo,points,weights)
    !
    DO i=1,npi
       !
       Pref  = 0.0_dp
       DPref = 0.0_dp
       DeP   = 0.0_dp
       Bt    = 0.0_dp
       det   = 0.0_dp
       afin  = 0.0_dp
       !
       CALL shape_fun(Pref,points,i)
       CALL shape_der(DPref, points,i)
       !
       Bt  = MATMUL(DPref,TRANSPOSE(coord))    !Bt = [B^T]
       det = determinant(Bt)
       det = ABS(det)
       !
       CALL invert(Bt)                         !Bt = B^{-T} 
       !
       ! DeP = B^{-T}*DPref derivadas en el elemento actual (ndim x nod)
       !
       DeP = MATMUL(Bt, DPref)
       !
       afin = MATMUL(coord,Pref)
       !
       norma_residuo = norma_residuo + (func(afin,tf,nu,velo,sigma,ref_ele) &
                        - DOT_PRODUCT(velo,MATMUL(DeP,uh_node)))**2*det*weights(i)
       !
    END DO
    !
    norma_residuo = (theta_T**2)*norma_residuo
    !
    DEALLOCATE(mm,coord,Bt,DPref,DeP,f_node,uh_node,Pref,points,weights,velo,afin)
    !
  END SUBROUTINE n_residuo

  SUBROUTINE n_salto(malla,fisica,frontera,ele,uh,norma_salto)
    !
    TYPE(mesh)                :: malla
    TYPE(param)               :: fisica
    TYPE(bc)                  :: frontera
    REAL(kind=dp), INTENT(in) :: uh(:)
    INTEGER, INTENT(in)       :: ele
    REAL(kind=dp)             :: norma_salto,h_face,n_ext(malla%ndim),norma
    INTEGER                   :: j,ele_in,ele_out,face, lado_local,nod,ndim
    INTEGER, ALLOCATABLE      :: mm(:)
    REAL(kind=dp), ALLOCATABLE:: coord(:,:)
    !
    norma_salto = 0.0_dp; lado_local = 0; ele_in = 0; nod = 0; ndim = 0
    !
    nod  = malla%element%nod
    ndim = malla%ndim
    !
    ALLOCATE(mm(nod),coord(ndim,nod))
    mm = 0; coord = 0.0_dp
    !
    mm = malla%conec(:,ele)
    !
    DO j=1,nod
       coord(:,j)  = malla%coord(:,mm(j))
    END DO
    !
    ele_in = ele
    !
    ciclo:DO j=1,malla%element%face_loc
       !
       h_face = 0.0_dp; n_ext = 0.0_dp; face = 0; norma = 0.0_dp; ele_out = 0
       !
       face = malla%ele_face(j,ele)
       !
       IF(malla%neigh(j,ele) == -1 .AND. ALL(malla%ref_face(face)/=frontera%NN)) CYCLE ciclo
       !
       CALL util_hF_nF(j,coord,h_face,n_ext)
       !
       ! SOLO NEUMANN NULA!!!!!
       !
       IF(ANY(malla%ref_face(face) == frontera%NN)) THEN
          !
          CALL error_edge_frontera(malla,fisica,uh,ele_in,n_ext,h_face,norma)
          norma = 2.0_dp*norma
          !
       ELSE
          !
          ele_out = malla%neigh(j,ele)
          CALL error_edge_interno(malla,fisica,uh,ele_in,ele_out,n_ext,h_face,norma)
          !
       END IF
       !
       norma_salto = norma_salto + norma
       !
    END DO ciclo
    !
    norma_salto = 0.5_dp*norma_salto
    !
    DEALLOCATE(mm,coord)
    !
  END SUBROUTINE n_salto
  !
  SUBROUTINE error_edge_interno(malla,fisica,uh,ele_in,ele_out,n_ext,h_face,norma)
    !
    TYPE(mesh)                :: malla
    TYPE(param)               :: fisica
    INTEGER, INTENT(in)       :: ele_in,ele_out
    REAL(kind=dp),INTENT(in)  :: n_ext(:),h_face,uh(:)
    REAL(kind=dp),INTENT(out) :: norma
    INTEGER, ALLOCATABLE      :: mm_in(:),mm_out(:)
    REAL(kind=dp),ALLOCATABLE :: coord_in(:,:),coord_out(:,:),Bt_in(:,:),Bt_out(:,:),DPref(:,:),&
                                 DeP_in(:,:), DeP_out(:,:), uh_in(:),uh_out(:)
    REAL(kind=dp)             :: DePn_in, DePn_out,nu,theta_F
    INTEGER                   :: j,k,nod,ndim,ind_in,ind_out,ref_ele_in,ref_ele_out
    !
    nod = 0; ndim = 0
    !
    ndim  = malla%ndim
    nod   = malla%element%nod
    !
    ALLOCATE(mm_in(nod),mm_out(nod),coord_in(ndim,nod),coord_out(ndim,nod),Bt_in(ndim,ndim), Bt_out(ndim,ndim),&
             DPref(ndim,nod),DeP_in(ndim,nod),DeP_out(ndim,nod),uh_in(nod),uh_out(nod))
    !
    mm_in = 0; mm_out = 0; coord_in =0.0_dp; coord_out = 0.0_dp; Bt_in = 0.0_dp; Bt_out = 0.0_dp
    DPref = 0.0_dp; DeP_in = 0.0_dp; DeP_out = 0.0_dp
    uh_in = 0.0_dp; uh_out = 0.0_dp; norma = 0.0_dp
    ref_ele_in = 0; ref_ele_out = 0; ind_in = 0; ind_out = 0
    DePn_in = 0.0_dp; DePn_out = 0.0_dp
    nu = 0.0_dp
    theta_F = 0.0_dp
    !
    !
    SELECT CASE(ndim)
    CASE(2)
       !
       DPref(1,:) = (/-1.0_dp,1.0_dp,0.0_dp/)
       DPref(2,:) = (/-1.0_dp,0.0_dp,1.0_dp/)
       !
    CASE(3)
       !
       DPref(1,:) = (/-1.0_dp,1.0_dp,0.0_dp,0.0_dp/)
       DPref(2,:) = (/-1.0_dp,0.0_dp,1.0_dp,0.0_dp/)
       DPref(3,:) = (/-1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
       !
    CASE default
       PRINT*,'Problema con normal!! (S: error_edge, M: aposteriori)'
       STOP
    END SELECT
    !
    mm_in  = malla%conec(:,ele_in)
    mm_out = malla%conec(:,ele_out)
    !
    ref_ele_in  = malla%ref_elem(ele_in)
    ref_ele_out = malla%ref_elem(ele_out)
    !
    DO j=1,fisica%ndominios
       !
       IF(fisica%ref_dominio(j) == ref_ele_in)  ind_in  = j
       IF(fisica%ref_dominio(j) == ref_ele_out) ind_out = j
       !
    END DO
    !
    nu = fisica%nu(ind_in)
    !
    theta_F = (1.0_dp/SQRT(nu))*MIN(1.0_dp, h_face/SQRT(nu))
    !
    uh_in  = uh(mm_in)
    uh_out = uh(mm_out)
    !
    DO j=1,nod
       coord_in(:,j)  = malla%coord(:,mm_in(j))
       coord_out(:,j) = malla%coord(:,mm_out(j))
    END DO
    !
    Bt_in  = MATMUL(DPref,TRANSPOSE(coord_in))    !Bt = [B^T]
    Bt_out = MATMUL(DPref,TRANSPOSE(coord_out))   !Bt = [B^T]
    !
    CALL invert(Bt_in)                            ! Bt = B^{-T} 
    CALL invert(Bt_out)                           ! Bt = B^{-T} 
    !
    ! DeP = B^{-T}*DPref derivadas en el elemento actual (ndim x nod)
    !
    DeP_in  = MATMUL(Bt_in, DPref)   
    DeP_out = MATMUL(Bt_out,DPref)
    !
    DePn_in  = DOT_PRODUCT(MATMUL(DeP_in ,uh_in), n_ext)
    DePn_out = DOT_PRODUCT(MATMUL(DeP_out,uh_out),n_ext)
    !
    norma = (nu**2)*h_face*theta_F*((DePn_in - DePn_out)**2)
    !
    DEALLOCATE(mm_in,mm_out,coord_in,coord_out,Bt_in,Bt_out,DPref,DeP_in,DeP_out,uh_in,uh_out)
    !
  END SUBROUTINE error_edge_interno
  !
    SUBROUTINE error_edge_frontera(malla,fisica,uh,ele_in,n_ext,h_face,norma)
    !
    TYPE(mesh)                :: malla
    TYPE(param)               :: fisica
    INTEGER, INTENT(in)       :: ele_in
    REAL(kind=dp),INTENT(in)  :: n_ext(:),h_face,uh(:)
    REAL(kind=dp),INTENT(out) :: norma
    INTEGER, ALLOCATABLE      :: mm_in(:)
    REAL(kind=dp),ALLOCATABLE :: coord_in(:,:),Bt_in(:,:),DPref(:,:),DeP_in(:,:),uh_in(:)
    REAL(kind=dp)             :: DePn_in,nu,theta_F, signo
    INTEGER                   :: j,k,nod,ndim,ind_in,ref_ele_in
    !
    nod = 0; ndim = 0
    !
    ndim  = malla%ndim
    nod   = malla%element%nod
    !
    ALLOCATE(mm_in(nod),coord_in(ndim,nod),Bt_in(ndim,ndim),DPref(ndim,nod),DeP_in(ndim,nod),uh_in(nod))
    !
    mm_in = 0; coord_in =0.0_dp; Bt_in = 0.0_dp
    DPref = 0.0_dp; DeP_in = 0.0_dp
    uh_in = 0.0_dp; norma = 0.0_dp
    ref_ele_in = 0; ind_in = 0
    DePn_in = 0.0_dp
    nu = 0.0_dp
    theta_F = 0.0_dp
    !
    SELECT CASE(ndim)
    CASE(2)
       !
       DPref(1,:) = (/-1.0_dp,1.0_dp,0.0_dp/)
       DPref(2,:) = (/-1.0_dp,0.0_dp,1.0_dp/)
       !
    CASE(3)
       !
       DPref(1,:) = (/-1.0_dp,1.0_dp,0.0_dp,0.0_dp/)
       DPref(2,:) = (/-1.0_dp,0.0_dp,1.0_dp,0.0_dp/)
       DPref(3,:) = (/-1.0_dp,0.0_dp,0.0_dp,1.0_dp/)
       !
    CASE default
       PRINT*,'Problema con normal!! (S: error_edge, M: aposteriori)'
       STOP
    END SELECT
    !
    mm_in  = malla%conec(:,ele_in)
    !
    ref_ele_in  = malla%ref_elem(ele_in)
    !
    DO j=1,fisica%ndominios
       !
       IF(fisica%ref_dominio(j) == ref_ele_in)  THEN
          ind_in  = j
          EXIT
       END IF
       !
    END DO
    !
    nu = fisica%nu(ind_in)
    !
    theta_F = (1.0_dp/SQRT(nu))*MIN(1.0_dp, h_face/SQRT(nu))
    !
    uh_in  = uh(mm_in)
    !
    DO j=1,nod
       coord_in(:,j)  = malla%coord(:,mm_in(j))
    END DO
    !
    Bt_in  = MATMUL(DPref,TRANSPOSE(coord_in))    !Bt = [B^T]
    !
    CALL invert(Bt_in)                            ! Bt = B^{-T} 
    !
    ! DeP = B^{-T}*DPref derivadas en el elemento actual (ndim x nod)
    !
    DeP_in  = MATMUL(Bt_in, DPref)   
    !
    DePn_in  = DOT_PRODUCT(MATMUL(DeP_in ,uh_in), n_ext)
    !
    norma = (nu**2)*h_face*theta_F*((DePn_in)**2)
    !
    DEALLOCATE(mm_in,coord_in,Bt_in,DPref,DeP_in,uh_in)
    !
  END SUBROUTINE error_edge_frontera
  !
END MODULE error
