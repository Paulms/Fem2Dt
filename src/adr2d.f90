PROGRAM adr2d
  !-----------------------------------------------------------------------------------------
  ! Este programa resuelve el problema de Adveccion-Difusion-Reaccion con                  |
  ! condiciones de Dirichlet/Neumann en la frontera:                                       |
  !                                                                                        |
  ! du/dt -div(\nu*\nabla u) + \mathbf{velocidad}\cdot\nabla u + \sigma u = f,   en \Omega |
  !                                                               u = u_0, en \Gamma_D     |
  !                                        \nu \cdot\nabla u\cdot n = 0,   en \Gamma_N     |
  !                                                                                        |
  ! Utilizamos una aproximacion por elementos finitos P1 o Q1.                             |
  !                                                                                        |
  ! El sistema lineal obtenido es resuelto utilizando una  descomposicion LU de            |
  ! la matriz.                                                                             |
  !                                                                                        |
  ! La matriz esta almacenada en formato CRS o bien CCS                                    | 
  !                                                                                        |
  !                                                                                        |
  !  IN:                                                                                   |
  !  --                                                                                    |
  !      - malla                                                                           |
  !      - segundo miembro f de la ecuacion                                                |
  !      - condicion de frontera u0                                                        |
  !                                                                                        |
  ! OUT:                                                                                   |
  ! ---                                                                                    |
  !      - uh: array que contiene los valores de la solucion en los nodos de la malla.     |
  !                                                                                        |
  !                                                                                        |
  ! OBS:                                                                                   |
  ! ---                                                                                    |
  !                                                                                        |
  !      HECHO POR   : Rodolfo Araya                                                       |
  !                    Departamento de Ingenieria Matematica                               |
  !                    Universidad de Concepcion                                           |
  !                    Casilla 160-C                                                       |
  !                    Concepcion, CHILE                                                   |
  !                                                                                        |
  !                    Tel   : (56-41) 20.31.16                                            |
  !                    Fax   : (56-41) 52.20.55                                            |
  !                    e-mail: rodolfo.araya@udec.cl                                       |
  !     Modificado por: Paul Mendez                                                    |
  !                     e-mail: pmendez@ing-mat.udec.cl                                   |
  !                                                                                        |
  !      Version 1   : Septiembre de 2000 (con almacenamiento banda y Cholesky)            |
  !      Version 2   : Enero de 2001 (con almacenamiento Morse y GMRES de SLATEC)          |
  !      Version 3   : Octubre de 2001 (con tipos malla y crs)                             |
  !      Version 4   : Diciembre 2001 (extension a Q1 usando una malla tipo MESH)          |
  !      Version 5   : Abril 2002 (adaptacion anisotropica de mallas usando BL2D)          |
  !      Version 6   : Mayo 2002 (resolucion del sistema lineal usando SuperLU y ccs)      |
  !      Version 7   : Agosto 2011 (re-escritura general + uso de MKL)                     !
  !      Version 8   : Junio 2015 (revision general + comentarios)                         !
  !      Version 9   : Diciembre 2016 (revision para resolver problemas evolutivos)        !
  !-----------------------------------------------------------------------------------------
  !
  ! Los modulos que usaremos:
  !
  USE decimal       ! declaracion de la precision
  USE tipos         ! definicion de los diferentes tipos a usar
  USE storage       ! lo relacionado con el tipo de almacenamiento de la matriz
  USE boundary      ! lo relativo a las condiciones de frontera
  USE system        ! los algoritmos para resolver el sistema lineal
  USE plot          ! las subrutinas que generan los archivos de visualizacion
  USE util          ! algunas subrutinas utiles (lectura de datos, etc)
  USE error         ! calculo de errores a priori y/o a posteriori
  USE funciones
  !
  IMPLICIT NONE
  !
  INTEGER                    :: visu,solver,aprox,metodo, ierr
  REAL(KIND=dp), ALLOCATABLE :: uh(:)
  REAL(KIND=dp)              :: tiempo1,tiempo2
  CHARACTER(LEN=32)          :: name_visu
  TYPE(mesh)                 :: malla
  TYPE(sparse)               :: A
  TYPE(bc)                   :: frontera
  TYPE(param)                :: fisica
  TYPE(time)                 :: tiempo
  
  INTEGER                    :: nT, i, j
  REAL(KIND=dp)              :: tt
  REAL(KIND=dp), ALLOCATABLE :: ua(:)

  tiempo%final = 0.0_dp
  tiempo%inicial = 0.0_dp
  tiempo%nn = 0
  tiempo%delta = 0.0_dp
  ierr = 0
  tt = 0.0_dp
  !
  CALL cpu_TIME(tiempo1)
  !
  !
  ! Lectura de la malla
  !
  WRITE(*,'(/,a,/)')'Leyendo los datos de entrada: util_input'
  CALL util_input(malla,fisica,frontera,metodo,solver,aprox,visu,name_visu, tiempo)
  !
  ! Creacion de la estructura de la  matriz del sistema lineal (formato CRS)
  !
  CALL estructura_A(malla,A)
  tiempo%delta = tiempo%final/tiempo%nn

  allocate(ua(malla%nnode))
  ua = 0.0_dp
  CALL u_inicial(ua)
!
! Guardamos datos iniciales
!
CALL storage_f(malla,fisica,metodo,uh, ua, 0.0_dp, tiempo%delta)
WRITE(*,'(/,a,/)')'Generando los archivos para la visualizacion: plot_result'
CALL plot_results(malla,visu,name_visu,uh,0)
DEALLOCATE(uh)
  DO i = 1,tiempo%nn
    tt = tt + tiempo%delta
    PRINT *,'Calculos en el paso de tiempo: ',tt
    ! 
    ! Creacion de la matriz del sistema lineal en el formato CRS (Morse)
    !
    WRITE(*,'(/,a,/)')'Creando la matriz del sistema lineal: storage_A'
    CALL storage_A(malla,fisica,metodo,A, tiempo%delta)
    !
    ! Creacion del segundo miembro del sistema lineal
    !
    WRITE(*,'(/,a,/)')'Creando el segundo miembro del sistema lineal: storage_f'
    CALL storage_f(malla,fisica,metodo,uh, ua, tt, tiempo%delta)
    ! Imponiendo la condicion de frontera de tipo Dirichlet
    !
    WRITE(*,'(/,a,/)')'Imponiendo las condiciones de frontera: frontera_condition'
    CALL boundary_condition(malla,frontera,fisica,A,uh, tt)   
    !
    ! Calculo de la solucion del sistema lineal
    !
    WRITE(*,'(/,a,/)')'Resolviendo el sistema lineal: system_solution'
    CALL system_solution(A,solver,uh)
    ua = uh
    !
    ! Generando el archivo para la visualizacion
    !
    WRITE(*,'(/,a,/)')'Generando los archivos para la visualizacion: plot_result'
    CALL plot_results(malla,visu,name_visu,uh,i)
    DEALLOCATE(uh)
  END DO
  !
  CALL cpu_TIME(tiempo2)
  !
  PRINT*,'tiempo de CPU = ',tiempo2-tiempo1
  DEALLOCATE(A%aa,A%row,A%column,STAT=ierr)
  IF(ierr/=0) THEN
      PRINT*,'problemas liberando la memoria!! (M: Sistemas S: solucion)'
      STOP
  END IF
  !
  ! Calculo del error a priori y/o a posteriori
  !
  IF(aprox /=0) THEN
     WRITE(*,'(/,a,/)')'Calculando el error a priori y/o a posteriori: estimacion'
     CALL estimacion(malla,fisica,frontera,uh,aprox,tiempo%final)
  END IF
  !
END PROGRAM adr2d

