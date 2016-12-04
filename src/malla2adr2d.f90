PROGRAM malla2adr2d
  !
  ! Programa que genera 
  !
  !
  ! Ultima modificacion: Miercoles 19 de Octubre de 2011
  !
  !
  ! 2D: Triangulo P1 y P2
  !
  ! lado 1:  n1     n2
  ! lado 2:  n2     n3
  ! lado 3:  n3     n1
  !
  !
  ! 3D: Tetraedro P1 y P2
  !
  ! cara 1: n1      n2       n3
  ! cara 2: n1      n2       n4
  ! cara 3: n1      n3       n4
  ! cara 4: n2      n3       n4
  !
  USE decimal
  !
  IMPLICIT NONE
  !
  INTEGER                  :: tipo,i,j,k,ndim,nnode,nelem,nface,face_loc,nod,nod_face,mino,maxo,cara, ss
  REAL(kind=dp),ALLOCATABLE:: coord(:,:) 
  INTEGER, ALLOCATABLE     :: conectivity(:,:),ref_node(:),ref_elem(:)
  INTEGER, ALLOCATABLE     :: conec_face(:,:),ref_face(:), mini(:), maxi(:),link(:), nn(:),&
                              suma(:),ele_face(:,:), vecino(:,:),num_nod_cara(:,:)
  CHARACTER(len=32)        :: name,name_node,name_ele,name_face,name_vecino
  !
  tipo = 0; ndim = 0; nnode = 0; nelem = 0; nface = 0; face_loc = 0; nod_face = 0
  !
  WRITE(*,*) 'Nombre del archivo que contiene la malla? '
  READ(*,*) name
  !
  WRITE(*,*) 'Tipo de malla (1 = TRIANGLE P1, 2 = TRAINGLE P2, 3 = TETGEN P1, 4 = TETGEN P2)'
  READ(*,*) tipo
  !
  OPEN(unit=10,file='nodos.dat'    ,status='replace',action='write')    
  OPEN(unit=20,file='elementos.dat',status='replace',action='write') 
  OPEN(unit=30,file='caras.dat'    ,status='replace',action='write')
  !
  name_node   = TRIM(name)//'.node'
  name_ele    = TRIM(name)//'.ele' 
  name_vecino = TRIM(name)//'.neigh' 
  !
  SELECT CASE(tipo)
  CASE(1) !P1 2D
     name_face   = TRIM(name)//'.edge'
     ndim        = 2
     nod         = 3
     face_loc    = 3
     nod_face    = 2
  CASE(2) !P2 2D
     name_face   = TRIM(name)//'.edge'
     ndim        = 2
     nod         = 6
     face_loc    = 3
     nod_face    = 2
  CASE(3) ! P1 3D
     name_face   = TRIM(name)//'.face'
     ndim        = 3
     nod         = 4
     face_loc    = 4
     nod_face    = 3
  CASE(4) ! P2 3D
     name_face   = TRIM(name)//'.face'
     ndim        = 3
     nod         = 10
     face_loc    = 4
     nod_face    = 3
  CASE default
     PRINT*,'Problema de tipo!'
     STOP
  END SELECT
  !
  OPEN(unit=50,file=TRIM(name_node),  status='old',action='read')    
  OPEN(unit=60,file=TRIM(name_ele),   status='old',action='read')    
  OPEN(unit=70,file=TRIM(name_face),  status='old',action='read')
  OPEN(unit=80,file=TRIM(name_vecino),status='old',action='read')
  !
  PRINT*,'----------------------------------------------------'
  PRINT*,' Malla -----> nodos.dat, caras.dat y elementos.dat'
  PRINT*,'----------------------------------------------------'
  !
  READ(50,*) nnode
  READ(60,*) nelem
  READ(70,*) nface
  !
  ALLOCATE(coord(ndim,nnode),conectivity(nod,nelem),ref_node(nnode),ref_elem(nelem),nn(nod),vecino(face_loc,nelem))
  !
  coord = 0.0_dp; conectivity = 0; ref_node = 0; ref_elem = 0; vecino = 0
  !
  !
  ALLOCATE(conec_face(nod_face,nface), ref_face(nface), mini(nface), maxi(nface),link(MAX(nod_face*nnode,nface)),&
           ele_face(face_loc,nelem), suma(nod_face*nnode), num_nod_cara(face_loc,nod_face))
  !
  conec_face = 0; ref_face = 0; mini = 0; maxi = 0; link = 0; ele_face = 0; suma = 0; num_nod_cara = 0
  !
  SELECT CASE(ndim)
  CASE(2)
     !
     num_nod_cara(1,:) = (/1,2/)
     num_nod_cara(2,:) = (/2,3/)
     num_nod_cara(3,:) = (/3,1/)
     !
  CASE(3)
     !
     num_nod_cara(1,:) = (/1,3,2/)
     num_nod_cara(2,:) = (/1,4,3/)
     num_nod_cara(3,:) = (/1,2,4/)
     num_nod_cara(4,:) = (/2,3,4/)
     !
  CASE default
     PRINT*,'Problema de dimension!!'
     STOP
  END SELECT
  !
  PRINT *,"leyendo ",name_node,"..."
  DO i=1,nnode
     READ(50,*) j,coord(:,i),ref_node(i)
  END DO
  !
  PRINT *,"leyendo ",name_ele,"..."
  DO i=1,nelem
     READ(60,*)j,conectivity(:,i),ref_elem(i)
  END DO
  !
  PRINT *,"leyendo ",name_face,"..."
  DO i=1,nface
     READ(70,*)j,conec_face(:,i),ref_face(i)
  END DO
  !
  DO i=1,nface
     mini(i) = MINVAL(conec_face(:,i))
     maxi(i) = MAXVAL(conec_face(:,i))
  END DO
  !
  !
  PRINT *,"leyendo ",name_vecino,"..."
  READ(80,*)
  !
  ! leyendo los vecinos de los triangulos y/o tetraedros
  !
  SELECT CASE(ndim)
  CASE(2)
     DO i=1,nelem
        READ(80,*)j,vecino(2,i), vecino(3,i),vecino(1,i)
     END DO
  CASE(3)
     DO i=1,nelem
        READ(80,*)j,vecino(4,i),vecino(3,i),vecino(2,i),vecino(1,i)
     END DO
     !
  CASE default
     PRINT*,'Problema de dimension (1)!!'
     STOP
  END SELECT
  !
  CLOSE(50)
  CLOSE(60)
  CLOSE(70)
  CLOSE(80)
  !
  ! Calculo de la conectividad por lados/caras
  !
  PRINT *,"Calculando número de caras..."
  CALL caras_num(nface,conec_face,mini,maxi,link,suma)
  !
  DO i=1,nelem
     !
     nn = 0; ss = 0 
     !
     nn = conectivity(:,i)
     !
     ! Ciclo sobre los lados y/o caras
     !
     DO j=1,face_loc
        mino = 0; ss = 0; cara = 0; maxo = 0
        !
        ss   = SUM(nn(num_nod_cara(j,:)))    !nn(j) + nn(MOD(j,3)+1)
        mino = MINVAL(nn(num_nod_cara(j,:))) !(nn(j),nn(MOD(j,3)+1))
        maxo = MAXVAL(nn(num_nod_cara(j,:))) !(nn(j),nn(MOD(j,3)+1))
        !
        cara = suma(ss)
        !
        DO WHILE (mini(cara) /= mino .OR. maxi(cara) /= maxo)
           cara = link(cara)
        END DO
        !
        ele_face(j,i) = cara
        !
     END DO
     !
  END DO
  !
  DEALLOCATE(mini,maxi,link,suma,num_nod_cara)
  !
  ! generando los archivos de salida
  !
  PRINT *,"Escribiendo datos..."
  WRITE(10,*) '# Numero de nodos de la malla:'
  WRITE(10,'(i7)') nnode
  !
  WRITE(10,*) '# Coordenadas de los nodos y sus referencias:'
  !
  SELECT CASE(ndim)
  CASE(2)
     DO i=1,nnode
        WRITE(10,'(2(F22.10,1x),i4)') coord(:,i),ref_node(i)
     END DO
  CASE(3)
    DO i=1,nnode
        WRITE(10,'(3(F22.10,1x),i4)') coord(:,i),ref_node(i)
     END DO
  CASE default
     PRINT*,'problema ndim en impresion coordenadas!!'
     STOP
  END SELECT
  !
  WRITE(20,*)'# Numero de elementos de la malla:'
  WRITE(20,'(i7)') nelem
  !
  WRITE(20,*)'# Conectividad nodal de los elementos y sus referencias:'
  !
  SELECT CASE(tipo)
     !
  CASE(1) !P1 2D
     DO i=1,nelem
        WRITE(20,'(1x,4(i7,2x))') conectivity(:,i),ref_elem(i)
     END DO
     !
  CASE(2) !P2 2D
     DO i=1,nelem
        WRITE(20,'(1x,7(i8,2x))') conectivity(:,i),ref_elem(i)
     END DO
     !
  CASE(3) ! P1 3D
     DO i=1,nelem
        WRITE(20,'(1x,5(i8,2x))') conectivity(:,i),ref_elem(i)
     END DO
     !
  CASE(4) ! P2 3D
     DO i=1,nelem
        WRITE(20,'(1x,11(i8,2x))') conectivity(:,i),ref_elem(i)
     END DO
  CASE default
    PRINT*,'Problema de tipo de elemento en escritura de conectividad!'
     STOP
  END SELECT
  !
  WRITE(20,*)'# Conectividad de los elementos con respecto a los lados/caras:'
  DO i =1, nelem
     WRITE(20,*)ele_face(:,i)
  END DO
  !
  WRITE(20,*)'# Numero de los elementos vecinos, por lado/cara, a cada elemento:'
  DO i=1,nelem
     WRITE(20,*) vecino(:,i)
  END DO
  !
  WRITE(30,*)'# Numero de lados/caras de la malla'
  WRITE(30,*) nface
  !
  WRITE(30,*) '# Nodos que componen cada lado/cara y referencia del lado/cara:'
  DO i=1,nface
     WRITE(30,*) conec_face(:,i),' ',ref_face(i)
  END DO
  !
  CLOSE(10)
  CLOSE(20)
  CLOSE(30)
  !
  DEALLOCATE(coord,conectivity,ref_node,ref_elem,nn,vecino,conec_face,ref_face, ele_face)
  !
CONTAINS
  !
  SUBROUTINE caras_num(nface,conec_face,mini,maxi,link,suma)
    !
    INTEGER, INTENT(in)    :: nface, conec_face(:,:),mini(:),maxi(:)
    INTEGER, INTENT(inout) :: link(:), suma(:)
    INTEGER                :: i, l, m, ss
    !
    link = 0; suma = 0
    !
    DO i=1,nface
       ss = 0; l = 0; m = 0
       !
       ss = SUM(conec_face(:,i))
       ! ss = conec_face(1,i) + conec_face(2,i)
       !
       IF(suma(ss) == 0) THEN
          suma(ss) = i
          !link(ss) = i
       ELSEIF(mini(suma(ss)) /= mini(i).OR.maxi(suma(ss)) /= maxi(i)) THEN
          l = suma(ss)
          !
          DO
             !IF(link(l) == 0.OR.link(l)==l) THEN
             IF(link(l) == 0) THEN
                link(l) = i
             ELSE
                m = link(l)
                IF(mini(m) /= mini(i).OR.maxi(m) /= maxi(i)) THEN
                   l = m  
                   CYCLE
                END IF
             END IF
             EXIT
          END DO
       END IF
       !
    END DO
    !
  END SUBROUTINE caras_num
  !
END PROGRAM malla2adr2d

