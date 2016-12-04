MODULE tipos
  !
  ! Aqui se declaran los nuevos tipos de variables
  !
  ! sparse  = tipo de una matriz almacenada en Morse (CRS) por filas
  ! fem     = tipo de elemento finito que compone la malla
  ! mesh    = tipo de una malla de elementos finitos
  ! param   = tipo que entrega los parametros fisicos asociados al problema
  ! bc      = tipo que define la condicion de frontera
  !
  USE decimal
  !
  IMPLICIT NONE
  !
  TYPE sparse
     !
     ! nzero       = numero de elementos no nulos de la matriz
     ! row, column = arreglos de filas y columnas
     ! aa          = arreglo que contiene los elementos no nulos de la matriz
     !
     INTEGER                    :: n,nzero          
     INTEGER, ALLOCATABLE       :: row(:),column(:) 
     REAL(KIND=dp), ALLOCATABLE :: aa(:)      
     !   
  END TYPE sparse
  !
  TYPE fem
     !
     ! definicion del elemento finito a usar (solo continuos):
     !
     ! name     = nombre TRIAP12D, QUADQ12D, TETRP13D, HEXAQ13D (interpolacion)
     ! nod      = numero de nodos de cada elemento
     ! nodf     = numero de nodos por cara
     ! face_loc = numero de lados (2D) y/o caras (3D) del elemento
     ! dof      = grados de libertad por nodo
     ! npi      = numero de puntos de integracion
     ! tipo     = tipo geometrico del elemento
     !
     CHARACTER(len=8) :: name
     INTEGER          :: nod
     INTEGER          :: nodf
     INTEGER          :: face_loc
     INTEGER          :: dof
     INTEGER          :: npi
     CHARACTER(len=20):: tipo
     !
  END TYPE fem
  !
  TYPE mesh
     !
     ! ndim  = dimension del problema (2 o 3)
     ! nnode = numero total de nodos de la malla
     ! nelem = numero total de elementos de la malla
     ! nedge = numero total de lados (caras) en la frontera de la malla
     ! nod   = numero de nodos de cada elemento
     !
     INTEGER                    :: ndim,nnode,nelem,nedge
     TYPE(fem)                  :: element
     INTEGER, ALLOCATABLE       :: conec(:,:),ref_node(:),ref_elem(:),ele_face(:,:),&
                                   conec_face(:,:),ref_face(:),neigh(:,:)
     REAL(kind=dp), ALLOCATABLE :: coord(:,:)
     !
  END TYPE mesh
  !
  TYPE param
     !
     ! ndominios              = numero total de los diferentes dominios del problema
     ! ref_dominio            = numero de referencias de los diferentes materiales
     ! nu(i), aa(i), sigma(i) = coeficientes fisicos de la ecuacion
     !
     INTEGER                   :: ndominios
     INTEGER, ALLOCATABLE      :: ref_dominio(:)
     REAL(kind=dp), ALLOCATABLE:: nu(:), aa(:,:), sigma(:)
     !
  END TYPE param
  !
  TYPE bc
     ! Tipo donde se definen las condiciones de Frontera
     !
     ! n_dirichlet = numero total de condiciones de Dirichlet
     ! n_neumann   = numero total de condicione de Neumann
     ! DD(i)       = numero de referencia de la i-esima condicion de Dirichlet
     ! fix(i,j)    = Dice si la j-esima variable esta bloqueada (=1) o no (=0) en la i-esima frontera de Dirichlet
     !  NN(i)      = numero de referencia de la i-esima condicion de Neumann
     !     
     INTEGER              :: n_dirichlet,n_neumann
     INTEGER, ALLOCATABLE :: DD(:),fix(:,:), NN(:)
     !
  END TYPE bc
  !
  TYPE time
    REAL (kind=dp)  :: inicial
    REAL (kind=dp)  :: final
    INTEGER         :: nn
    REAL (kind=dp)  :: delta

  END TYPE time
END MODULE tipos
