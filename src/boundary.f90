MODULE boundary
  !
  ! Modulo relacionado con las condiciones de Dirichlet
  !
  !
  USE decimal
  USE funciones
  USE tipos
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC  boundary_condition
  !
CONTAINS
  !
  SUBROUTINE boundary_condition(malla,dirichlet,fisica,A,uh, tt)
    !
    TYPE(mesh)                  :: malla
    TYPE(bc)                    :: dirichlet
    TYPE(param)                 :: fisica
    TYPE(sparse)                :: A
    REAL(KIND=dp),INTENT(inout) :: uh(:)
    REAL(KIND=dp)               :: tt
    !
    CALL boundary_dirichlet_crs(malla,dirichlet,fisica,A,uh, tt)
    !
  END SUBROUTINE boundary_condition
  !
  SUBROUTINE  boundary_dirichlet_crs(malla,dirichlet,fisica,A,uh, tt)
    !
    ! Estrategia para imponer las condiciones de frontera de tipo 
    ! Dirichlet.
    !
    ! IN:
    !      malla: triangulacion del dominio
    !      dirichlet: informacion de la condicion de Dirichlet
    !
    ! OUT:
    !    A: malla del sistema lineal modificada
    !    uh: lado derecho modificado del sistema lineal
    !
    !
    TYPE(mesh)                  :: malla
    TYPE(sparse)                :: A
    TYPE(param)                 :: fisica
    TYPE(bc)                    :: dirichlet
    REAL(KIND=dp),INTENT(inout) :: uh(:)
    REAL(kind=dp)               :: xyz(malla%ndim), u0, nu, velo(malla%ndim), sigma, tt
    INTEGER                     :: i,j,dof,ref_node,nn
    !
    dof = 0
    !
    dof = malla%element%dof
    !
    nu = 0.0_dp;  velo = 0.0_dp; sigma = 0.0_dp  !!FALTA IDENTIFICAR EL DOMINIO SI SON VARIOS!!!!!
    !
    !PRINT*,'ref_node = ', malla%ref_node
    !
    DO i=1,malla%nnode
       !
       nn       = 0
       ref_node = 0
       xyz      = 0.0_dp
       u0       = 0.0_dp
       !
       IF(ALL(malla%ref_node(i) /= dirichlet%DD)) CYCLE ! si el nodo "i" no es Dirichlet, siga iterando
       !
       nn       = i                                 ! numero de ecuaciones asociadas al nodo i 
       ref_node = malla%ref_node(i)                 ! referencia del nodo i
       xyz      = malla%coord(:,i)                  ! coordenadas del nodo i
       !
       DO j = A%row(nn),A%row(nn+1)-1
          A%aa(j) = 0.0_dp
          IF( A%column(j) == nn ) A%aa(j) = 1.0_dp
       END DO
       !
       u0 = u_dirichlet(xyz,nu,velo,sigma,ref_node, tt) ! valor de la cond Dirichlet en el nodo i
       !
       uh(nn) = u0  
       !
    END DO
    !
  END SUBROUTINE boundary_dirichlet_crs
  !
END MODULE boundary
