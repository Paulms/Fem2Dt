MODULE funciones
  !
  USE decimal
  !
  IMPLICIT NONE
  !
CONTAINS
  !
  FUNCTION func(xyz,tt,nu,aa,sigma,ref)
    !
    ! RHS function
    !
    REAL(KIND=dp), INTENT(IN) :: xyz(:),nu,aa(:),sigma, tt
    INTEGER, INTENT(in)       :: ref
    integer                   :: ndim
    REAL(kind=dp)             :: func,x,y,z,x0,y0
    !
    func = 0.0_dp; x = 0.0_dp; y = 0.0_dp; z = 0.0_dp; ndim = 0
    x0 = 0.05_dp; y0 = 0.05_dp
    !
    ndim = SIZE(xyz)
    !
    SELECT CASE(ndim)
    CASE(2)
       x = xyz(1)
       y = xyz(2)
    CASE(3)
       x = xyz(1)
       y = xyz(2)
       z = xyz(3)
    CASE default
       PRINT*,'Problema de dimension!! M: Funciones, F: func'
       STOP
    END SELECT
    !
    func = 0.0_dp! -4.0_dp/(((x-x0)**2+(y-y0)**2)**2)
    !
    !
    ! cuadrado2
!!$    func = -nu*(200.0_dp*x**2*y*(1.0_dp-2.0_dp*y)*(1.0_dp-y)-(800.0_dp*(1.0_dp-x))*x*y*(1.0_dp-2.0_dp*y)*(1.0_dp-y)+&
!!$            200.0_dp*(1.0_dp-x)**2*y*(1.0_dp-2.0_dp*y)*(1.0_dp-y)-400.0_dp*(1.0_dp-x)**2*x**2*(1.0_dp-y)-&
!!$            200.0_dp*(1.0_dp-x)**2*x**2*(1.0_dp-2.0_dp*y)+400.0_dp*(1.0_dp-x)**2*x**2*y)-&
!!$            (400.0_dp*(1.0_dp-x))*x**2*y*(1.0_dp-2.0_dp*y)*(1.0_dp-y)+400.0_dp*(1.0_dp-x)**2*x*y*(1.0_dp-2.0_dp*y)*(1.0_dp-y)+&
!!$            300.0_dp*(1.0_dp-x)**2*x**2*(1.0_dp-2.0_dp*y)*(1.0_dp-y)-600.0_dp*(1.0_dp-x)**2*x**2*y*(1.0_dp-y)-&
!!$            300.0_dp*(1.0_dp-x)**2*x**2*y*(1.0_dp-2.0_dp*y)
!!$    func = 32.0_dp*nu*y-32.0_dp*nu*y**2+32.0_dp*nu*x-32.0_dp*nu*x**2+16.0_dp*y-16.0_dp*y**2-&
!!$           64.0_dp*x*y+32.0_dp*x*y**2+16.0_dp*x-16.0_dp*x**2+32.0_dp*x**2*y
    !
    !
    ! Esta es la funcion del RHS si la solucion del problema esta dada por:
    !
    !    u(x,y) := 16*x*y*(1-x)*(1-y)
    !
!!$    func = 32.0_dp*nu*xyz(2)-32.0_dp*nu*xyz(2)**2+32.0_dp*nu*xyz(1)-32.0_dp*nu*xyz(1)**2+16.0_dp*aa(1)*xyz(2)-16.0_dp*aa(1)*xyz(2)**2&
!!$          -32.0_dp*aa(1)*xyz(1)*xyz(2)+32.0_dp*aa(1)*xyz(1)*xyz(2)**2+16.0_dp*aa(2)*xyz(1)-32.0_dp*aa(2)*xyz(1)*xyz(2)&
!!$          -16.0_dp*aa(2)*xyz(1)**2+32.0_dp*aa(2)*xyz(1)**2*xyz(2)+16.0_dp*sigma*xyz(1)*xyz(2)-16.0_dp*sigma*xyz(1)*xyz(2)**2&
!!$          -16.0_dp*sigma*xyz(1)**2*xyz(2)+16.0_dp*sigma*xyz(1)**2*xyz(2)**2
    !
!!$    func = -0.2393511232e+02_dp * EXP(-0.50e+02_dp * x ** 2 + 0.50e+02_dp * x - 0.25e+02_dp - &
!!$            0.50e+02_dp * y ** 2 + 0.50e+02_dp * y) - 0.4987520800e+02_dp * EXP(-0.50e+02_dp *x ** 2 + &
!!$            0.50e+02_dp * x - 0.25e+02_dp - 0.50e+02_dp * y ** 2 + 0.50e+02_dp * y) * x** 2 +&
!!$            0.4987520800e+02_dp * EXP(-0.50e+02_dp * x ** 2 + 0.50e+02_dp * x - 0.25e+02_dp- &
!!$            0.50e+02_dp * y ** 2 + 0.50e+02_dp * y) * x - 0.4987520800e+02_dp * EXP(-0.50e+02_dp* x** 2 +&
!!$            0.50e+02_dp * x - 0.25e+02_dp - 0.50e+02_dp * y ** 2 + 0.50e+02_dp * y)*y** 2 + &
!!$            0.4987520800e+02_dp * EXP(-0.50e+02_dp * x ** 2 + 0.50e+02_dp * x - 0.25e+02_dp - &
!!$            0.50e+02_dp * y ** 2 + 0.50e+02_dp * y) * y
    !
    !
  END FUNCTION func
  !
  FUNCTION u_dirichlet(xyz,nu,aa,sigma,ref, tt)
    !
    !The Dirichlet BC
    !
    REAL(KIND=dp), INTENT(IN) :: xyz(:),nu,aa(:),sigma
    INTEGER, INTENT(in)       :: ref
    integer                   :: ndim
    REAL(kind=dp)             :: u_dirichlet,x,y,z,x0,y0, tt
    !
    u_dirichlet = 0.0_dp; x = 0.0_dp; y = 0.0_dp; z = 0.0_dp; ndim = 0
    x0 = 0.05_dp; y0 = 0.05_dp
    !
    ndim = SIZE(xyz)
    !
    SELECT CASE(ndim)
    CASE(2)
       x = xyz(1)
       y = xyz(2)
    CASE(3)
       x = xyz(1)
       y = xyz(2)
       z = xyz(3)
    CASE default
       PRINT*,'Problema de dimension!! M: Funciones, F: func'
       STOP
    END SELECT
    !
    SELECT CASE(ref)
    CASE(2)
       !
       u_dirichlet =  5.0_dp
       !
    !CASE (2)
    !    u_dirichlet =  0.0_dp
    CASE default
       PRINT*,'Mala referencia del nodo!! M: Funciones, S: u_dirichlet', ref
       STOP
    END SELECT
    !
  END FUNCTION u_dirichlet
  !
  FUNCTION u_neumann(xyz,ref)
    !
    ! The Neumann  BC
    !
    REAL(KIND=dp), INTENT(IN) :: xyz(:)
    INTEGER, INTENT(in)       :: ref
    REAL(kind=dp)             :: u_neumann,x,y,z
    INTEGER                   :: ndim
    !
    u_neumann = 0.0_dp; ndim = 0; x = 0.0_dp; y = 0.0_dp; z = 0.0_dp
    !
    ndim = SIZE(xyz)
    !
    SELECT CASE(ndim)
    CASE(2)
       !
       x = xyz(1)
       y = xyz(2)
       !
    CASE(3)
       !
       x = xyz(1)
       y = xyz(2)
       z = xyz(3)
       !
    CASE default
       PRINT*,'Problemas con la dimension!!! M: funciones, F: u_newmann'
    END SELECT
    !
    SELECT CASE(ref)
       !
    CASE(30) !caso 3D
       !
       u_neumann = 0.0_dp
       !
    CASE default
       PRINT*,'Mala referencia del lado!! M: Funciones, S: u_neumann'
       STOP
    END SELECT
    !
  END FUNCTION u_neumann
  !
  FUNCTION u_ex(xyz,nu,ref,tf)
    !
    ! Solucion exacta del problema
    !
    REAL(KIND=dp), INTENT(IN) :: xyz(:),nu
    INTEGER, INTENT(in)       :: ref
    REAL(kind=dp)             :: u_ex, x, y,x0,y0,tf
    !
    u_ex = 0.0_dp; x = 0.0_dp; y = 0.0_dp
    x0 = 0.05_dp; y0 = 0.05_dp; tf = 0.0_dp
    !
    x = xyz(1)
    y = xyz(2)
    !
    !cuadrado 2:
    !u_ex = 16.0_dp*x*y*(1.0_dp-x)*(1.0_dp-y)
    !
    !u_ex = 100.0_dp*(1.0_dp-x)**2*x**2*y*(1.0_dp-2.0_dp*y)*(1.0_dp-y)
    !
    u_ex = 1.0_dp/((x-x0)**2 + (y-y0)**2) ! x - (EXP((x-1.0_dp)/nu) - EXP(-1.0_dp/nu))/(1.0_dp - EXP(-1.0_dp/nu)) 
    !
  END FUNCTION u_ex

  FUNCTION grad_ex(xyz,nu,ref,tf)
    !
    ! Gradiente de la solucion exacta del problema
    !
    REAL(KIND=dp), INTENT(IN) :: xyz(:),nu
    INTEGER, INTENT(in)       :: ref
    REAL(kind=dp)             :: grad_ex(SIZE(xyz)), x, y, x0,y0,tf
    !
    grad_ex = 0.0_dp; x = 0.0_dp; y = 0.0_dp
    x0 = 0.05_dp; y0 = 0.05_dp;tf = 0.0_dp
    !
    x = xyz(1)
    y = xyz(2)
    !
    !cuadrado2:
    !grad_ex(1) = 16.0_dp*y*(1.0_dp-x)*(1.0_dp-y)-16.0_dp*x*y*(1.0_dp-y)
    !grad_ex(2) = 16.0_dp*x*(1.0_dp-x)*(1.0_dp-y)-16.0_dp*x*y*(1.0_dp-x)
    !
    grad_ex(1) = (-2.0_dp*(x-x0))/(((x-x0)**2+(y-y0)**2)**2) !1.0_dp - (EXP((x-1.0_dp)/nu))/(nu*(1.0_dp - EXP(-1.0_dp/nu))) 
    grad_ex(2) = (-2.0_dp*(y-y0))/(((x-x0)**2+(y-y0)**2)**2) !0.0_dp
    !
!!$    grad_ex(1) = -200.0_dp*(1.0_dp-x)*x**2*y*(1.0_dp-2.0_dp*y)*(1.0_dp-y)+200.0_dp*(1.0_dp-x)**2*x*y*(1.0-dp-2.0_dp*y)*(1.0_dp-y)
!!$    grad_ex(2) = 100.0_dp*(1.0_dp-x)**2*x**2*(1.0_dp-2.0_dp*y)*(1.0_dp-y)-200.0_dp*(1.0_dp-x)**2*x**2*y*(1.0_dp-y)-&
!!$                 100.0_dp*(1.0_dp-x)**2*x**2*y*(1.0_dp-2.0_dp*y)
    !
  END FUNCTION grad_ex
  !
END MODULE funciones






