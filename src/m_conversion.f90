MODULE m_conversion
  !
  USE decimal
  USE m_base
  USE m_formats
  USE m_share
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC  side, reference
  !
CONTAINS  
  !
  SUBROUTINE side(nelem,nnode,node,nside,edge,nodes_face)
    !
    ! program to put a number, unique, to each side of the triangulation
    ! it works only for p1 meshes
    !
    INTEGER, DIMENSION(3)                :: turn=(/2,3,1/)
    INTEGER,INTENT(in)                   :: nelem,nnode
    INTEGER, DIMENSION(:,:),INTENT(in)   :: node
    INTEGER, DIMENSION(:,:),INTENT(out)  :: edge
    INTEGER, INTENT(out)                 :: nside
    INTEGER, DIMENSION(:,:),POINTER      :: nodes_face,list
    INTEGER, DIMENSION (:),POINTER       :: head
    INTEGER                              :: nlist,i,j,ierr
    !
    nside = nelem + nnode -1 + 1000      ! Euler's relation con 1000 hoyos
    !
    ALLOCATE(nodes_face(2,nside),stat=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas en la alocacion del array nodes_faces'
       STOP
    ENDIF
    !
    nodes_face = 0
    !
    ALLOCATE(head(nnode),stat=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas en la alocacion del array head'
       STOP
    ENDIF
    !
    head = 0
    !
    ALLOCATE(list(2,2*nside),stat=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas en la alocacion del array list'
       STOP
    ENDIF
    !
    list = 0
    !
    nside    = 0		  ! reset of the variable
    nlist = 0
    !
    DO i = 1, nelem
       DO j = 1, 3
          CALL numbering(node(j,i), node(turn(j),i),edge(j,i))
       END DO
    END DO
    !
  CONTAINS 
    !
    SUBROUTINE numbering(s1, s2,number)
      !
      !numbering of the side j of triangle i
      !
      INTEGER,INTENT(in) :: s1, s2
      INTEGER,INTENT(out)::number
      INTEGER:: ia, smin, smax, ilist
      !    
      IF (s1 < s2) THEN
         smin = s1
         smax = s2
      ELSE
         smin = s2
         smax = s1
      END IF
      
      ! find the face (we know the nodes)
      ilist = head(smin)
      DO
         IF (ilist == 0) EXIT
         ia = list(1,ilist)
         IF (nodes_face(2,ia) == smax) THEN
            number=ia
            RETURN 
         ENDIF
         ilist = list(2,ilist)
      END DO
      !
      !new face
      nside = nside + 1
      number=nside
      ia = nside
      nodes_face(1,ia) = smin
      nodes_face(2,ia) = smax
      ! add to list of faces with smin
      nlist = nlist + 1
      list(1,nlist) = ia
      list(2,nlist) = head(smin)
      head(smin) = nlist
      ! add to list of faces with smax
      nlist = nlist + 1
      list(1,nlist) = ia
      list(2,nlist) = head(smax)
      head(smax) = nlist
    END SUBROUTINE numbering
    
  END SUBROUTINE side
  !
  SUBROUTINE reference(nside,nodes_face,prefix,iterate,ref_side)
    !
    INTEGER,INTENT(in)                 :: nside,iterate
    INTEGER, DIMENSION(:,:),INTENT(in) :: nodes_face
    CHARACTER(len=*),INTENT(in)        :: prefix
    INTEGER, POINTER                   :: ref_side(:)
    TYPE(g_)                           :: x_g
    TYPE(c_)                           :: x_c      
    INTEGER                            :: i, node1,node2
    !
    ALLOCATE(ref_side(nside))
    !
    ref_side = 0
    !
    call ouvrir(stdin, filename3(prefix,0, "g"))
    call g_read(x_g)
    call fermer(stdin)

    call ouvrir(stdin, filename3(prefix,iterate, "c"))
    call c_read(x_c)
    call fermer(stdin)         
    !
    DO i = 1, nside
       node1 = nodes_face(1,i)
       node2 = nodes_face(2,i)
       !
       ref_side(i) = physique_a(node1,node2,x_g,x_c)
       !
    END DO
    !
  END SUBROUTINE reference
  !
END MODULE m_conversion

