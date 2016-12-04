MODULE plot
  !
  ! Modulo que contiene las rutinas para generar los archivos graficos de salida 
  !
  ! (formatos leibles por Tecplot, Ensight y Vigie)
  !
  USE decimal
  USE tipos
  USE util
  !
  PRIVATE
  !
  PUBLIC  plot_results
  !
CONTAINS
  !
  SUBROUTINE plot_results(malla,visu,nameb,uh, iteracion)
    !
    ! Subrutina que genera archivos para la visualizacion
    ! de la solucion para diferentes tipos de softwares
    ! (Tecplot, Ensight, Vigie, etc.)
    !
    TYPE(mesh)                   :: malla
    REAL(KIND=dp),INTENT(IN)     :: uh(:)
    INTEGER                      :: visu
    CHARACTER(LEN=*),INTENT(IN)  :: nameb
    CHARACTER(LEN=32)             :: name
    INTEGER                      :: iteracion
    CHARACTER(LEN=8)             :: string
    !
    write(string,'(i8)') iteracion
    name = TRIM(ADJUSTL(nameb))//'.'//TRIM(ADJUSTL(string))
    !

    SELECT CASE(visu)
    CASE(1)
       !
       ! Usando Tecplot para visualizar los resultados
       !
       CALL plot_tecplot(name,malla,uh)
       !
    CASE(2)
       !
       ! Usando Ensight para visualizar los resultados
       !
       CALL plot_ensight(name,malla,uh)
       !
    CASE(3)
       !
       ! Usando Vigie para visualizar los resultados
       !
       CALL plot_vigie(name,malla,uh)
       !
    CASE(4)
       !
       ! Usando Paraview para visualizar los resultados
       !
       CALL plot_paraview(name,malla,uh)
       !
    CASE default
       !
       PRINT*,'Metodo de visualizacion aun no implementado (M: plot, S: plot_results)!!'
       !
    END SELECT
    !
  END SUBROUTINE plot_results
  !
  SUBROUTINE plot_ensight(name,malla,uh)
    !
    ! subrutina que imprime la malla y el resultado usando
    ! el visualizador ENSIGHT
    !
    TYPE(mesh)                   :: malla
    REAL(KIND=dp),INTENT(IN)     :: uh(:)
    CHARACTER(LEN=*),INTENT(IN)  :: name
    CHARACTER(LEN=32)            :: name_case,name_geo,name_concen,name_ele
    INTEGER                      :: i,j
    !
    name_case    = TRIM(ADJUSTL(name))//".case"
    name_geo     = TRIM(ADJUSTL(name))//".geo"
    name_concen  = TRIM(ADJUSTL(name))//".concen"
    !
    OPEN(10,file = name_case,   status='replace', action='write' )
    OPEN(20,file = name_geo,    status='replace', action='write' )
    OPEN(30,file = name_concen, status='replace', action='write' )
    !
    ! creacion del archivo .CASE
    !
    WRITE(10,'(a)')'#'
    WRITE(10,'(a)')'# Using the Ensight visualizator'
    WRITE(10,'(a,a/)') '# Case file: ',name_case
    !
    WRITE(10,'(a,/)')'FORMAT'   
    WRITE(10,'(a,/)') 'type:         ensight gold'
    !
    WRITE(10,'(a,/)')'GEOMETRY' 
    WRITE(10,'(a,a,/)') 'model:      ',name_geo
    !
    WRITE(10,'(a,/)')'VARIABLE' 
    !
    WRITE(10,'(a,a)') 'scalar per node:      concentracion ',name_concen
    !
    ! Genereando el archivo GEO
    !
    WRITE(20,'(a)')'Primer intento de usar Ensight'
    WRITE(20,'(a)')'(ojo con los formatos)'
    !
    WRITE(20,'(a)')'node id assign'
    WRITE(20,'(a)')'element id assign'
    !
    WRITE(20,'(a)')'part'
    WRITE(20,'(i10)') 1
    WRITE(20,'(a)')'Mesh'
    WRITE(20,'(a)')'coordinates'
    WRITE(20,'(i10)') malla%nnode
    !
    DO i=1,malla%ndim
       DO j=1,malla%nnode
          WRITE(20,'(e12.5)') malla%coord(i,j)
       END DO
    END DO
    !
    IF(malla%ndim==2) THEN
       DO i=1,malla%nnode
          WRITE(20,'(e12.5)') 0.0_dp
       END DO
    END IF
    !
    malla%element%name = TRIM(malla%element%name)
    !
    SELECT CASE (malla%element%name)
    CASE('TRIAP12D')
       name_ele = 'tria3'
       !
       WRITE(20,'(a)') name_ele
       WRITE(20,'(i10)') malla%nelem
       !
       DO i=1,malla%nelem
          WRITE(20,'(3(i10))') malla%conec(:,i)
       END DO
       !
    CASE('QUADQ12D')
       name_ele = 'quad4'
       !
       WRITE(20,'(a)') name_ele
       WRITE(20,'(i10)') malla%nelem
       !
       DO i=1,malla%nelem
          WRITE(20,'(4(i10))') malla%conec(:,i)
       END DO

    CASE('TETRP13D')
       name_ele = 'tetra4'
       !
       WRITE(20,'(a)') name_ele
       WRITE(20,'(i10)') malla%nelem
       !
       DO i=1,malla%nelem
          WRITE(20,'(4(i10))') malla%conec(:,i)
       END DO
       !
    CASE('HEXAQ13D')
       name_ele = 'hexa8'
       !
       WRITE(20,'(a)') name_ele
       WRITE(20,'(i10)') malla%nelem
       !
       DO i=1,malla%nelem
          WRITE(20,'(8(i10))') malla%conec(:,i)
       END DO
    CASE default
       PRINT*,'Elemento aun no implementado (M: plot, S: plot_result_ensight)'
       STOP
    END SELECT
    !
    ! Escribiendo el archivo con la concentracion (solucion)
    !    
    WRITE(30,'(a,a)')'Concentracion (nodal), nombre del archivo:',name_concen
    WRITE(30,'(a)')'part'
    WRITE(30,'(i10)') 1
    WRITE(30,'(a)')'coordinates'
    !
    DO j=1,malla%nnode
       WRITE(30,'(e12.5)') uh(j)
    END DO
    !
    CLOSE(10)
    CLOSE(20)
    CLOSE(30)
    !
  END SUBROUTINE plot_ensight
  !
  !
  SUBROUTINE plot_vigie(name,malla,uh)
    !
    ! subrutina que imprime la malla y el resultado usando
    ! el visualizador VIGIE
    !
    TYPE(mesh), INTENT(IN)       :: malla
    REAL(KIND=dp),INTENT(IN)     :: uh(:)
    CHARACTER(LEN=*),INTENT(IN)  :: name
    CHARACTER(LEN=32)            :: name_desc,name_ascii2d
    INTEGER                      :: i,j,iunit1,iunit2
    !
    name_desc    = TRIM(ADJUSTL(name))//".desc"
    name_ascii2d = TRIM(ADJUSTL(name))//".ascii2d"
    !
    CALL util_get_unit(iunit1)
    OPEN(iunit1,file=name_desc,   status='replace', action='write' )
    !
    CALL util_get_unit(iunit2)
    OPEN(iunit2,file=name_ascii2d,status='replace', action='write' )
    !
    ! creacion del archivo .DESC
    !
    WRITE(iunit1,*)" ascii2d"
    WRITE(iunit1,*) "./",name_ascii2d
    !
    ! creacion del archivo .ASCII2D
    !
    WRITE(iunit2,*)" points ", malla%nnode
    !
    DO j=1,malla%nnode
       WRITE(iunit2,*) (malla%coord(i,j),"   ",i=1,2)
    END DO
    !
    SELECT CASE(malla%element%name)
    CASE('TRIAP12D')       
       WRITE(iunit2,*)" triangles ", malla%nelem
    CASE ('QUADQ12D')
       WRITE(iunit2,*)" quadrangles ", malla%nelem
    CASE default
       PRINT*,'wrong element type'
       STOP
    END SELECT
    !
    DO j=1,malla%nelem
       WRITE(iunit2,*) (malla%conec(i,j),"  ", i=1,malla%element%nod)
    END DO
    !
    WRITE(20,*)" scalars  Concentracion "
    DO j=1,malla%nnode
       WRITE(20,*) uh(j)
    END DO
    !
    CLOSE(iunit1)
    CLOSE(iunit2)
    !
  END SUBROUTINE plot_vigie
  !
  SUBROUTINE plot_tecplot(name,malla,uh)
    !
    ! subrutina que imprime la malla y el resultado usando
    ! el visualizador TECPLOT
    ! 
    ! (por el momento solo valido para P1, Q1 en 2D y 3D)
    !
    TYPE(mesh)                   :: malla
    REAL(KIND=dp),INTENT(IN)     :: uh(:)
    CHARACTER(LEN=*),INTENT(IN)  :: name
    CHARACTER(LEN=32)            :: name_dat,name_ele
    INTEGER                      :: i,j,nod,iunit1
    !
    name_dat = TRIM(ADJUSTL(name))//".dat"  
    !
    CALL util_get_unit(iunit1)
    OPEN(iunit1,file=name_dat,   status='replace', action='write' )
    !
    ! creacion del archivo .dat
    !
    SELECT CASE (malla%element%name)
    CASE('TRIAP12D')
       !
       name_ele = 'FETRIANGLE'
       !
    CASE('QUADQ12D')
       !
       name_ele = 'FEQUADRILATERAL'
       !
    CASE('TETRP13D')
       !
       name_ele = 'FETETRAHEDRON'
       !
    CASE('HEXAQ13D')
       !
       name_ele = 'FEBRICK'
       !
    CASE default
       PRINT*,'Elemento aun no implementado (M: plot, S: plot_result_ensight)'
       STOP
    END SELECT
    !
    nod = 0
    nod = malla%element%nod
    !
    WRITE(iunit1,*) 'TITLE="FINITE ELEMENT DISPLACEMENT" '
    IF(malla%ndim==2) THEN
       WRITE(iunit1,*) 'VARIABLES = "X", "Y", "Concentracion"'
    ELSEIF(malla%ndim==3) THEN
       WRITE(iunit1,*) 'VARIABLES = "X", "Y", "Z", "Concentracion"'
    ELSE
       PRINT*,'Problem with dimension!! M: Plot S:plot_solution_tecplot (ndim =',malla%ndim,')'
    END IF
    !
    WRITE(iunit1,*)'ZONE DATAPACKING=POINT, ZONETYPE=',TRIM(name_ele),', NODES=',malla%nnode,', ELEMENTS=',malla%nelem
    !
    DO j=1,malla%nnode
       WRITE(iunit1,*) (malla%coord(i,j),'   ',i=1,malla%ndim), uh(j)
    END DO
    !
    DO j=1,malla%nelem
       WRITE(iunit1,*) (malla%conec(i,j),'  ', i=1,nod)    
    END DO
    !
    CLOSE(iunit1)
    !
  END SUBROUTINE plot_tecplot
  !
  SUBROUTINE plot_paraview(name,malla,uh)
    !
    ! subrutina que imprime la malla y el resultado usando
    ! el visualizador PARAVIEW
    ! 
    ! (por el momento solo valido para P1, Q1 en 2D y 3D)
    !
    TYPE(mesh)                   :: malla
    REAL(KIND=dp),INTENT(IN)     :: uh(:)
    CHARACTER(LEN=*),INTENT(IN)  :: name
    CHARACTER(LEN=32)            :: name_dat
    INTEGER                      :: i,j,nod,iunit1,cell_type
    !
    name_dat = TRIM(ADJUSTL(name))//".vtk"  
    !
    CALL util_get_unit(iunit1)
    OPEN(iunit1,file=name_dat,   status='replace', action='write' )
    !
    ! creacion del archivo .vtk
    !
    SELECT CASE (malla%element%name)
    CASE('TRIAP12D')
       !
       cell_type = 5
       !
    CASE('QUADQ12D')
       !
       cell_type = 9 
       !
    CASE('TETRP13D')
       !
       cell_type = 10
       !
    CASE('HEXAQ13D')
       !
       cell_type = 12
       !
    CASE default
       PRINT*,'Elemento aun no implementado (M: plot, S: plot_result_vtk)'
       STOP
    END SELECT
    !
    nod = 0
    nod = malla%element%nod
    !
    WRITE(iunit1,'(A)') '# vtk DataFile Version 2.0'
    WRITE(iunit1,'(A)') 'Desplazamiento calculado con ADR2D'
    WRITE(iunit1,'(A)') 'ASCII'
    !
    WRITE(iunit1,'(A)')'DATASET UNSTRUCTURED_GRID'
    !
    ! Los nodos:
    !
    WRITE(iunit1,'(A,1x,i12,1x,A)')'POINTS',malla%nnode,'float'
    !
    DO j=1,malla%nnode
       IF(malla%ndim==2) THEN
          WRITE(iunit1,'(3(E30.15,2x))') (malla%coord(i,j),i=1,malla%ndim) ,0.0
       ELSEIF(malla%ndim==3) THEN
          WRITE(iunit1,'(3(E30.15,2x))') (malla%coord(i,j),i=1,malla%ndim)
       ELSE
          PRINT*,'Problem with dimension!! M: Plot S:plot_solution_vtk (ndim =',malla%ndim,')'
       END IF
       !
    END DO
    !
    ! Los elementos
    !
    WRITE(iunit1,'(A)') 
    WRITE(iunit1,'(A,1x,i10,1x,i10)') 'CELLS', malla%nelem, malla%nelem*(1 + nod) 
    !
    DO j=1,malla%nelem
       IF(cell_type == 5) THEN !triangulos P1 en 2D
          WRITE(iunit1,'(i2,1x,3(i12,2x))') nod, (malla%conec(i,j)-1, i=1,nod)
       ELSEIF(cell_type == 10) THEN !tetraedros P1 en 3D
          WRITE(iunit1,'(i2,1x,4(i12,2x))') nod, (malla%conec(i,j)-1, i=1,nod)
       ELSE
          PRINT*,'Problem with cell_type!! M: Plot S:plot_solution_vtk (cell_type =',cell_type,')'
       END IF
    END DO
    !
    WRITE(iunit1,'(A)') 
    WRITE(iunit1,'(A,1x,i12)')'CELL_TYPES', malla%nelem
    DO i=1,malla%nelem
       WRITE(iunit1,'(i2)') cell_type
    END DO
    !
    WRITE(iunit1,'(A)') 
    WRITE(iunit1,'(A,1x,i12)')'CELL_DATA', malla%nelem
    WRITE(iunit1,'(A)')'SCALARS materiales int 1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO i=1,malla%nelem
       WRITE(iunit1,'(i5)') malla%ref_elem(i)
    END DO
    !
    ! Las variables calculadas:
    !
    WRITE(iunit1,'(A)') 
    WRITE(iunit1,'(A,1x,i12)')'POINT_DATA', malla%nnode
    !
    ! Solucion discreta:
    !
    WRITE(iunit1,'(A)')'SCALARS Concentracion float  1'
    WRITE(iunit1,'(A)')'LOOKUP_TABLE default'
    DO i=1,malla%nnode
       WRITE(iunit1,'(E30.15)') uh(i)
    END DO
    !
    CLOSE(iunit1)
    !
  END SUBROUTINE plot_paraview
  !
END MODULE plot



