MODULE quadrature
  !
  ! Modulo que contiene la info necesaria para hacer las diferentes 
  ! integraciones numericas
  !
  ! (def puntos de integracion y pesos, valores de las derivadas y
  !  las funciones de base en los ptos. de integracion)
  !
  USE tipos
  USE decimal
  !
  IMPLICIT NONE
  !
CONTAINS
  !
  SUBROUTINE sample(element,s,wt)
    ! returns the local coordinates of the integrating points
    IMPLICIT NONE 
    REAL(kind=dp),INTENT(out):: s(:,:),wt(:)  
    CHARACTER(*),INTENT(in)   :: element
    INTEGER                   :: nip 
    REAL(kind=dp)            :: root3, r15,w(3),v(9),b,c
    !
    root3 = 1.0_dp/SQRT(3.0_dp)
    r15 = .2_dp*SQRT(15.0_dp) 
    nip = UBOUND( s , 1 ) 
    w = (/5.0_dp/9.0_dp,8.0_dp/9.0_dp,5.0_dp/9.0_dp/)
    !  v=(/5.0_dp/9.0_dp*w,8.0_dp/9.0_dp*w,5.0_dp/9.0_dp*w/)
    v(1:3)=5.0_dp/9.0_dp*w
    v(4:6)=8.0_dp/9.0_dp*w
    v(7:9)=5.0_dp/9.0_dp*w
    SELECT CASE (element)
    CASE('line')
       SELECT CASE(nip)
       CASE(1)
          s(1,1)=0.0_dp
          wt(1)=2.0_dp
       CASE(2)
          s(1,1)=root3  
          s(2,1)=-s(1,1)  
          wt(1)=1.0_dp  
          wt(2)=1.0_dp
       CASE(3)
          s(1,1)=r15 
          s(2,1)=.0_dp
          s(3,1)=-s(1,1)
          wt = w
       CASE(4)
          s(1,1)=.861136311594053_dp
          s(2,1)=.339981043584856_dp
          s(3,1)=-s(2,1)  
          s(4,1)=-s(1,1)
          wt(1)=.347854845137454_dp
          wt(2)=.652145154862546_dp
          wt(3)=wt(2) 
          wt(4)=wt(1)
       CASE(5)
          s(1,1)=.906179845938664_dp
          s(2,1)=.538469310105683_dp
          s(3,1)=.0_dp
          s(4,1)=-s(2,1) 
          s(5,1)=-s(1,1)
          wt(1)=.236926885056189_dp
          wt(2)=.478628670499366_dp
          wt(3)=.568888888888889_dp
          wt(4)=wt(2) 
          wt(5)=wt(1)
       CASE(6)
          s(1,1)=.932469514203152_dp
          s(2,1)=.661209386466265_dp
          s(3,1)=.238619186083197_dp
          s(4,1)=-s(3,1) 
          s(5,1)=-s(2,1) 
          s(6,1)=-s(1,1)
          wt(1)=.171324492379170_dp
          wt(2)=.360761573048139_dp
          wt(3)=.467913934572691_dp
          wt(4)=wt(3)
          wt(5)=wt(2) 
          wt(6)=wt(1)
       CASE default
          PRINT*,"wrong number of integrating points for a line"
       END SELECT
    CASE('triangle')
       SELECT CASE(nip)
       CASE(1)   ! for triangles weights multiplied by .5
          s(1,1)=1.0_dp/3.0_dp
          s(1,2)=1.0_dp/3.0_dp
          wt(1)= .5_dp
       CASE(3)
          s(1,1) = 0.5_dp
          s(1,2) = 0.5_dp
          s(2,1) = 0.5_dp
          s(2,2) = 0.0_dp
          s(3,1) = 0.0_dp
          s(3,2) = 0.5_dp
          wt(1)  = 1.0_dp/3.0_dp
          wt(2)  = wt(1) 
          wt(3)  = wt(1)   
          wt     = 0.5_dp*wt
       CASE(6)
          s(1,1)=.816847572980459_dp
          s(1,2)=.091576213509771_dp
          s(2,1)=s(1,2)
          s(2,2)=s(1,1) 
          s(3,1)=s(1,2)
          s(3,2)=s(1,2)
          s(4,1)=.108103018168070_dp
          s(4,2)=.445948490915965_dp
          s(5,1)=s(4,2) 
          s(5,2)=s(4,1) 
          s(6,1)=s(4,2)  
          s(6,2)=s(4,2)
          wt(1)=.109951743655322_dp
          wt(2)=wt(1)  
          wt(3)=wt(1)
          wt(4)=.223381589678011_dp
          wt(5)=wt(4)  
          wt(6)=wt(4)    
          wt = .5_dp*wt
       CASE(7)
          s(1,1)=1.0_dp/3.0_dp 
          s(1,2)=1.0_dp/3.0_dp
          s(2,1)=.797426985353087_dp 
          s(2,2)=.101286507323456_dp
          s(3,1)=s(2,2) 
          s(3,2)=s(2,1) 
          s(4,1)=s(2,2) 
          s(4,2)=s(2,2)
          s(5,1)=.470142064105115_dp
          s(5,2)=.059715871789770_dp
          s(6,1)=s(5,2) 
          s(6,2)=s(5,1)
          s(7,1)=s(5,1)
          s(7,2)=s(5,1)
          wt(1)=.225_dp
          wt(2)=.125939180544827_dp
          wt(3)=wt(2)
          wt(4)=wt(2)
          wt(5)=.132394152788506_dp
          wt(6)=wt(5)      
          wt(7)=wt(5)     
          wt = .5_dp*wt
       CASE(12)
          s(1,1)=.873821971016996_dp
          s(1,2)=.063089014491502_dp
          s(2,1)=s(1,2) 
          s(2,2)=s(1,1)
          s(3,1)=s(1,2) 
          s(3,2)=s(1,2)
          s(4,1)=.501426509658179_dp
          s(4,2)=.249286745170910_dp
          s(5,1)=s(4,2)
          s(5,2)=s(4,1)   
          s(6,1)=s(4,2) 
          s(6,2)=s(4,2)
          s(7,1)=.636502499121399_dp
          s(7,2)=.310352451033785_dp
          s(8,1)=s(7,1) 
          s(8,2)=.053145049844816_dp
          s(9,1)=s(7,2) 
          s(9,2) =s(7,1)
          s(10,1)=s(7,2)
          s(10,2)=s(8,2) 
          s(11,1)=s(8,2)
          s(11,2)=s(7,1)
          s(12,1)=s(8,2) 
          s(12,2)=s(7,2)
          wt(1)=.050844906370207_dp
          wt(2)=wt(1)
          wt(3)=wt(1)
          wt(4)=.116786275726379_dp
          wt(5)=wt(4)
          wt(6)=wt(4)
          wt(7)=.082851075618374_dp
          wt(8:12)=wt(7)      
          wt = .5_dp*wt
       CASE(16)
          s(1,1)=1.0_dp/3.0_dp 
          s(1,2)=1.0_dp/3.0_dp
          s(2,1)=.658861384496478_dp
          s(2,2)=.170569307751761_dp
          s(3,1)=s(2,2)   
          s(3,2)=s(2,1)
          s(4,1)=s(2,2)  
          s(4,2)=s(2,2)
          s(5,1)=.898905543365938_dp
          s(5,2)=.050547228317031_dp
          s(6,1)=s(5,2)
          s(6,2)=s(5,1) 
          s(7,1)=s(5,2)  
          s(7,2)=s(5,2)
          s(8,1)=.081414823414554_dp
          s(8,2)=.459292588292723_dp
          s(9,1)=s(8,2)  
          s(9,2)=s(8,1)
          s(10,1)=s(8,2) 
          s(10,2)=s(8,2)
          s(11,1)=.008394777409958_dp
          s(11,2)=.263112829634638_dp
          s(12,1)=s(11,1)    
          s(12,2)=.728492392955404_dp
          s(13,1)=s(11,2) 
          s(13,2)=s(11,1)  
          s(14,1)=s(11,2)
          s(14,2)=s(12,2)
          s(15,1)=s(12,2) 
          s(15,2)=s(11,1) 
          s(16,1)=s(12,2) 
          s(16,2)=s(11,2)
          wt(1)=.144315607677787_dp
          wt(2)=.103217370534718_dp
          wt(3)=wt(2)
          wt(4)=wt(2)
          wt(5)=.032458497623198_dp
          wt(6)=wt(5)   
          wt(7)=wt(5)
          wt(8)=.095091634267284_dp
          wt(9)=wt(8)   
          wt(10)=wt(8)
          wt(11)=.027230314174435_dp
          wt(12:16) = wt(11)  
          wt = .5_dp*wt
       CASE default
          PRINT*,"wrong number of integrating points for a triangle"
       END SELECT
    CASE ('quadrilateral')
       SELECT CASE (nip)
       CASE(1)
          s(1,1) = .0_dp 
          wt(1) = 4.0_dp
       CASE(4)
          s(1,1)=-root3
          s(1,2)= root3
          s(2,1)= root3
          s(2,2)= root3
          s(3,1)=-root3
          s(3,2)=-root3
          s(4,1)= root3
          s(4,2)=-root3
          wt = 1.0_dp
       CASE(9)
          s(1:7:3,1) = -r15
          s(2:8:3,1) = .0_dp
          s(3:9:3,1) =  r15
          s(1:3,2)   = r15
          s(4:6,2)   =  .0_dp
          s(7:9,2)   =-r15
          wt= v
       CASE default
          PRINT*,"wrong number of integrating points for a quadrilateral"
       END SELECT
    CASE('tetrahedron')    
       SELECT CASE(nip)
       CASE(1)          ! for tetrahedra weights multiplied by 1/6
          s(1,1)=.25_dp
          s(1,2)=.25_dp
          s(1,3)=.25_dp
          wt(1)=1.0_dp/6.0_dp
       CASE(4)
          s(1,1)=.58541020_dp
          s(1,2)=.13819660_dp
          s(1,3)=s(1,2)
          s(2,2)=s(1,1) 
          s(2,3)=s(1,2)  
          s(2,1)=s(1,2)
          s(3,3)=s(1,1) 
          s(3,1)=s(1,2)  
          s(3,2)=s(1,2)
          s(4,1)=s(1,2) 
          s(4,2)=s(1,2)  
          s(4,3)=s(1,2) 
          wt(1:4)=.25_dp/6.0_dp
       CASE(5)
          s(1,1)=.25_dp
          s(1,2)=.25_dp
          s(1,3)=.25_dp
          s(2,1)=.5_dp
          s(2,2)=1.0_dp/6.0_dp 
          s(2,3)=s(2,2)
          s(3,2)=.5_dp
          s(3,3)=1.0_dp/6.0_dp  
          s(3,1)=s(3,3)   
          s(4,3)=.5_dp
          s(4,1)=1.0_dp/6.0_dp
          s(4,2)=s(4,1)
          s(5,1)=1.0_dp/6.0_dp
          s(5,2)=s(5,1) 
          s(5,3)=s(5,1) 
          wt(1)=-.8_dp
          wt(2)=9.0_dp/20.0_dp 
          wt(3:5)=wt(2)   
          wt =wt/6.0_dp
       CASE(6)
          wt = 4.0_dp/3.0_dp        
          s(6,3) = 1.0_dp
          s(1,1)=-1.0_dp 
          s(2,1)=1.0_dp 
          s(3,2)=-1.0_dp 
          s(4,2)=1.0_dp 
          s(5,3)=-1.0_dp 
       CASE default
          PRINT*,"wrong number of integrating points for a tetrahedron"
       END SELECT
    CASE('hexahedron')
       SELECT CASE ( nip )
       CASE(1)
          s(1,1) = .0_dp
          wt(1) = 8.0_dp
       CASE(8)   
          s(1,1)= root3
          s(1,2)= root3
          s(1,3)= root3
          s(2,1)= root3
          s(2,2)= root3
          s(2,3)=-root3
          s(3,1)= root3
          s(3,2)=-root3
          s(3,3)= root3
          s(4,1)= root3
          s(4,2)=-root3
          s(4,3)=-root3
          s(5,1)=-root3
          s(5,2)= root3
          s(5,3)= root3
          s(6,1)=-root3
          s(6,2)=-root3
          s(6,3)= root3
          s(7,1)=-root3
          s(7,2)= root3
          s(7,3)=-root3
          s(8,1)=-root3
          s(8,2)=-root3
          s(8,3)=-root3
          wt = 1.0_dp                                               
       CASE(14)
          b=0.795822426_dp   
          c=0.758786911_dp
          wt(1:6)=0.886426593_dp   
          wt(7:) =  0.335180055_dp
          s(1,1)=-b 
          s(2,1)=b  
          s(3,2)=-b 
          s(4,2)=b
          s(5,3)=-b   
          s(6,3)=b
          s(7:,:) = c
          s(7,1)=-c  
          s(7,2)=-c  
          s(7,3)=-c 
          s(8,2)=-c 
          s(8,3)=-c
          s(9,1)=-c  
          s(9,3)=-c  
          s(10,3)=-c
          s(11,1)=-c
          s(11,2)=-c 
          s(12,2)=-c 
          s(13,1)=-c
       CASE(15)
          b=1.0_dp     
          c=0.674199862_dp
          wt(1)=1.564444444_dp 
          wt(2:7)=0.355555556_dp
          wt(8:15)=0.537777778_dp
          s(2,1)=-b  
          s(3,1)=b  
          s(4,2)=-b  
          s(5,2)=b
          s(6,3)=-b  
          s(7,3)=b  
          s(8:,:)=c  
          s(8,1)=-c
          s(8,2)=-c  
          s(8,3)=-c 
          s(9,2)=-c  
          s(9,3)=-c
          s(10,1)=-c 
          s(10,3)=-c  
          s(11,3)=-c 
          s(12,1)=-c
          s(12,2)=-c 
          s(13,2)=-c  
          s(14,1)=-c                          
       CASE(27)
          !         wt = (/5.0_dp/9.0_dp*v,8.0_dp/9.0_dp*v,5.0_dp/9.0_dp*v/)
          wt=0.0_dp !ojo:esto no esta bien!!
          s(1:7:3,1) = -r15
          s(2:8:3,1) = .0_dp
          s(3:9:3,1) =  r15
          s(1:3,3)   = r15
          s(4:6,3)   =  .0_dp
          s(7:9,3)   =-r15
          s(1:9,2)   = -r15
          s(10:16:3,1) = -r15
          s(11:17:3,1) = .0_dp
          s(12:18:3,1) =  r15
          s(10:12,3)   = r15
          s(13:15,3)   =  .0_dp
          s(16:18,3)   =-r15
          s(10:18,2)   = .0_dp
          s(19:25:3,1) = -r15
          s(20:26:3,1) = .0_dp
          s(21:27:3,1) =  r15
          s(19:21,3)   = r15
          s(22:24,3)   =  .0_dp
          s(25:27,3)   =-r15
          s(19:27,2)   =  r15
       CASE default
          PRINT*,"wrong number of integrating points for a hexahedron" 
       END SELECT
    CASE default
       PRINT*,"not a valid element type" 
    END SELECT
    RETURN
  END SUBROUTINE sample
  !
  !
  SUBROUTINE shape_fun(fun,points,i)
    IMPLICIT NONE  
    INTEGER,INTENT(in):: i
    REAL(kind=dp),INTENT(in)::points(:,:)
    REAL(kind=dp),INTENT(out)::fun(:)
    REAL(kind=dp) :: eta,xi,etam,etap,xim,xip,zetam,zetap,c1,c2,c3     !local variables
    REAL(kind=dp) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,x,y,z
    REAL(kind=dp) :: zeta,xi0,eta0,zeta0
    INTEGER::xii(20),etai(20),zetai(20),l,ndim,nod
    ndim = SIZE(points,2)
    nod  = SIZE(fun,1)  
    SELECT CASE (ndim)
    CASE(1) ! one dimensional cases
       xi=points(i,1)
       SELECT CASE(nod)
       CASE(2)
          t1=-1.0_dp-xi 
          t2=1.0_dp-xi
          fun(1)=t2/2.0_dp 
          fun(2)=-t1/2.0_dp
       CASE(3)
          t1=-1.0_dp-xi 
          t2=-xi 
          t3=1.0_dp-xi
          fun(1)=t2*t3/2.0_dp 
          fun(2)=-t1*t3 
          fun(3)=t1*t2/2.0_dp
       CASE(4)
          t1=-1.0_dp-xi 
          t2=-1.0_dp/3.0_dp-xi 
          t3=1.0_dp/3.0_dp-xi 
          t4=1.0_dp-xi
          fun(1)=t2*t3*t4*9.0_dp/16.0_dp  
          fun(2)=-t1*t3*t4*27.0_dp/16.0_dp
          fun(3)=t1*t2*t4*27.0_dp/16.0_dp 
          fun(4)=-t1*t2*t3*9.0_dp/16.0_dp
       CASE(5)
          t1=-1.0_dp-xi 
          t2=-0.5_dp-xi 
          t3=-xi 
          t4=0.5_dp-xi 
          t5=1.0_dp-xi
          fun(1)=t2*t3*t4*t5*2.0_dp/3.0_dp 
          fun(2)=-t1*t3*t4*t5*8.0_dp/3.0_dp
          fun(3)=t1*t2*t4*t5*4.0_dp 
          fun(4)=-t1*t2*t3*t5*8.0_dp/3.0_dp
          fun(5)=t1*t2*t3*t4*2.0_dp/3.0_dp
       CASE default
          PRINT*,"wrong number of nodes in shape_fun"
       END SELECT
    CASE(2) ! two dimensional cases
       c1=points(i,1)
       c2=points(i,2)
       c3=1.0_dp-c1-c2 
       xi=points(i,1)
       eta=points(i,2)
       etam=.25_dp*(1.0_dp-eta)
       etap=.25_dp*(1.0_dp+eta)
       xim=.25_dp*(1.0_dp-xi)
       xip=.25_dp*(1.0_dp+xi)
       SELECT CASE(nod)
       CASE(3)
          fun = (/c3,c1,c2/)  
       CASE(6)
          fun(1)=(2.0_dp*c3-1.0_dp)*c3 
          fun(2)=(2.0_dp*c1-1.0_dp)*c1
          fun(3)=(2.0_dp*c2-1.0_dp)*c2 
          fun(4)=4.0_dp*c3*c1
          fun(5)=4.0_dp*c1*c2      
          fun(6)=4.0_dp*c2*c3
       CASE(15) !no he cambiado el orden
          t1=c1-.25_dp  
          t2=c1-.5_dp 
          t3=c1-.75_dp   
          t4=c2-.25_dp
          t5=c2-.5_dp   
          t6=c2-.75_dp 
          t7=c3-.25_dp  
          t8=c3-.5_dp 
          t9=c3-.75_dp
          fun(1)=32.0_dp/3.0_dp*c1*t1*t2*t3   
          fun(12)=128.0_dp/3.0_dp*c1*c2*t1*t2
          fun(11)=64.0_dp*c1*c2*t1*t4     
          fun(10)=128.0_dp/3.0_dp*c1*c2*t4*t5
          fun(9)=32.0_dp/3.0_dp*c2*t4*t5*t6   
          fun(8)=128.0_dp/3.0_dp*c2*c3*t4*t5
          fun(7)=64.0_dp*c2*c3*t4*t7      
          fun(6)=128.0_dp/3.0_dp*c2*c3*t7*t8
          fun(5)=32.0_dp/3.0_dp*c3*t7*t8*t9   
          fun(4)=128.0_dp/3.0_dp*c3*c1*t7*t8
          fun(3)=64.0_dp*c3*c1*t1*t7      
          fun(2)=128.0_dp/3.0_dp*c3*c1*t1*t2
          fun(13)=128.0_dp*c1*c2*t1*c3    
          fun(15)=128.0_dp*c1*c2*c3*t4
          fun(14)=128.0_dp*c1*c2*c3*t7      
       CASE(4)
          fun=(/4.0_dp*xim*etam,4.0_dp*xip*etap,4.0_dp*xip*etam,4.0_dp*xim*etap/)
       CASE(8)
          fun=(/4.0_dp*etam*xim*(-xi-eta-1.0_dp),4.0_dp*xip*etam*(xi-eta-1.0_dp),&
               4.0_dp*etap*xip*(xi+eta-1.0_dp),4.0_dp*etap*xim*(-xi+eta-1.0_dp),&
               32.0_dp*xim*xip*etam,32.0_dp*etap*xip*etam,              &
               32.0_dp*xim*xip*etap,32.0_dp*etam*xim*etap/)
       CASE(9)
          etam = eta - 1.0_dp
          etap= eta + 1.0_dp
          xim = xi - 1.0_dp
          xip = xi + 1.0_dp
          fun=(/.25_dp*xi*xim*eta*etam,.25_dp*xi*xip*eta*etam,  &
               .25_dp*xi*xip*eta*etap,.25_dp*xi*xim*eta*etap,  &
               -.5_dp*xip*xim*eta*etam,-.5_dp*xi*xip*etap*etam,&
               -.5_dp*xip*xim*eta*etap,-.5_dp*xi*xim*etap*etam,&            
               xip*xim*etap*etam/)
       CASE default
          PRINT*,"wrong number of nodes in shape_fun"
       END SELECT
    CASE(3) ! three dimensional cases
       xi=points(i,1)
       eta=points(i,2)
       zeta=points(i,3)
       etam=1.0_dp-eta 
       xim=1.0_dp-xi  
       zetam=1.0_dp-zeta
       etap=eta+1.0_dp 
       xip=xi+1.0_dp   
       zetap=zeta+1.0_dp
       SELECT CASE(nod)
       CASE(4)
          fun(2) = xi   
          fun(3) = eta 
          fun(4) = zeta 
          fun(1) = 1.0_dp-fun(2)-fun(3)-fun(4)
       CASE(8)
          fun=(/.125_dp*xim*etam*zetam,.125_dp*xim*etam*zetap,.125_dp*xip*etam*zetap,&
               .125_dp*xip*etam*zetam,.125_dp*xim*etap*zetam,.125_dp*xim*etap*zetap,&
               .125_dp*xip*etap*zetap,.125_dp*xip*etap*zetam/)
       CASE(14) !type 6 element
          x = points(i,1)
          y = points(i,2)
          z = points(i,3)
          fun(1)=((x*y+x*z+2.0_dp*x+y*z+2.0_dp*y+2.0_dp*z+2.0_dp)*(x-1.0_dp)*&
               (y-1.0_dp)*(z-1.0_dp))/8.0_dp
          fun(2)=((x*y-x*z-2.0_dp*x+y*z+2.0_dp*y-2.0_dp*z-2.0_dp)*(x-1.0_dp)*&
               (y+1.0_dp)*(z-1.0_dp))/8.0_dp
          fun(3)=((x*y+x*z+2.0_dp*x-y*z-2.0_dp*y-2.0_dp*z-2.0_dp)*(x+1.0_dp)*&
               (y-1.0_dp)*(z-1.0_dp))/8.0_dp
          fun(4)=((x*y-x*z-2.0_dp*x-y*z-2.0_dp*y+2.0_dp*z+2.0_dp)*(x+1.0_dp)*&
               (y+1.0_dp)*(z-1.0_dp))/8.0_dp
          fun(5)=-((x*y-x*z+2.0_dp*x-y*z+2.0_dp*y-2.0_dp*z+2.0_dp)*(x-1.0_dp)*&
               (y-1.0_dp)*(z+1.0_dp))/8.0_dp
          fun(6)=-((x*y+x*z-2.0_dp*x-y*z+2.0_dp*y+2.0_dp*z-2.0_dp)*(x-1.0_dp)*&
               (y+1.0_dp)*(z+1.0_dp))/8.0_dp
          fun(7)=-((x*y-x*z+2.0_dp*x+y*z-2.0_dp*y+2.0_dp*z-2.0_dp)*(x+1.0_dp)*&
               (y-1.0_dp)*(z+1.0_dp))/8.0_dp
          fun(8)=-((x*y+x*z-2.0_dp*x+y*z-2.0_dp*y-2.0_dp*z+2.0_dp)*(x+1.0_dp)*&
               (y+1.0_dp)*(z+1.0_dp))/8.0_dp
          fun(9)=-((x+1.0_dp)*(x-1.0_dp)*(y+1.0_dp)*(y-1.0_dp)*(z-1.0_dp))/2.0_dp
          fun(10)=((x+1.0_dp)*(x-1.0_dp)*(y+1.0_dp)*(y-1.0_dp)*(z+1.0_dp))/2.0_dp
          fun(11)=-((x+1.0_dp)*(x-1.0_dp)*(y-1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp
          fun(12)=((x+1.0_dp)*(x-1.0_dp)*(y+1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp
          fun(13)=-((x-1.0_dp)*(y+1.0_dp)*(y-1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp
          fun(14)=((x+1.0_dp)*(y+1.0_dp)*(y-1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp      
       CASE(20)
          xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
          etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
          zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
          DO l=1,20
             xi0=xi*xii(l)
             eta0=eta*etai(l)
             zeta0=zeta*zetai(l)
             IF(l==4.OR.l==8.OR.l==16.OR.l==20) THEN
                fun(l)=.25_dp*(1.0_dp-xi*xi)*(1.0_dp+eta0)*(1.0_dp+zeta0)
             ELSE IF(l>=9.AND.l<=12)THEN
                fun(l)=.25_dp*(1.0_dp+xi0)*(1.0_dp-eta*eta)*(1.0_dp+zeta0)
             ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
                fun(l)=.25_dp*(1.0_dp+xi0)*(1.0_dp+eta0)*(1.0_dp-zeta*zeta)
             ELSE
                fun(l)=.125_dp*(1.0_dp+xi0)*(1.0_dp+eta0)*(1.0_dp+zeta0)*(xi0+eta0+zeta0-2)
             END IF
          END DO
       CASE default
          PRINT*,"wrong number of nodes in shape_fun"
       END SELECT
    CASE default
       PRINT*,"wrong number of dimensions in shape_fun"  
    END SELECT
    RETURN
  END SUBROUTINE shape_fun
  !
  !
  SUBROUTINE shape_der(der,points,i)
    IMPLICIT NONE
    INTEGER,INTENT(in):: i
    REAL(kind=dp),INTENT(in)::points(:,:)
    REAL(kind=dp),INTENT(out)::der(:,:)
    REAL(kind=dp)::eta,xi,zeta,xi0,eta0,zeta0,etam,etap,xim,xip,c1,c2,c3 ! local variables
    REAL(kind=dp):: t1,t2,t3,t4,t5,t6,t7,t8,t9 ,x2p1,x2m1,e2p1,e2m1,&
         zetam,zetap,x,y,z
    INTEGER :: xii(20), etai(20), zetai(20) ,l,ndim , nod   ! local variables
    ndim = SIZE(der,1)
    nod  = SIZE(der,2)
    SELECT CASE (ndim)
    CASE(1) ! one dimensional case
       xi=points(i,1)
       SELECT CASE (nod)
       CASE(2)
          der(1,1)=-0.5_dp
          der(1,2)=0.5_dp
       CASE(3)
          t1=-1.0_dp-xi 
          t2=-xi  
          t3=1.0_dp-xi
          der(1,1)=-(t3+t2)/2.0_dp
          der(1,2)=(t3+t1)    
          der(1,3)=-(t2+t1)/2.0_dp
       CASE(4)
          t1=-1.0_dp-xi 
          t2=-1.0_dp/3.0_dp-xi 
          t3=1.0_dp/3.0_dp-xi 
          t4=1.0_dp-xi
          der(1,1)=-(t3*t4+t2*t4+t2*t3)*9.0_dp/16.0_dp     
          der(1,2)=(t3*t4+t1*t4+t1*t3)*27.0_dp/16.0_dp 
          der(1,3)=-(t2*t4+t1*t4+t1*t2)*27.0_dp/16.0_dp 
          der(1,4)=(t2*t3+t1*t3+t1*t2)*9.0_dp/16.0_dp
       CASE(5)
          t1=-1.0_dp-xi 
          t2=-0.5_dp-xi 
          t3=-xi 
          t4=0.5_dp-xi 
          t5=1.0_dp-xi
          der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*2.0_dp/3.0_dp
          der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*8.0_dp/3.0_dp
          der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*4.0_dp
          der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*8.0_dp/3.0_dp
          der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*2.0_dp/3.0_dp
       CASE default
          PRINT*,"wrong number of nodes in shape_der"
       END SELECT
    CASE(2)      ! two dimensional elements
       xi=points(i,1)
       eta=points(i,2) 
       c2=xi 
       c3=eta 
       c1=1.0_dp-c2-c3
       etam=.25_dp*(1.0_dp-eta)
       etap=.25_dp*(1.0_dp+eta)
       xim=.25_dp*(1.0_dp-xi)
       xip=.25_dp*(1.0_dp+xi)
       x2p1=2.0_dp*xi+1.0_dp 
       x2m1=2.0_dp*xi-1.0_dp 
       e2p1=2.0_dp*eta+1.0_dp 
       e2m1=2.0_dp*eta-1.0_dp
       SELECT CASE (nod)
       CASE(3)
          der(1,1)=-1.0_dp
          der(1,2)= 1.0_dp
          der(1,3)= 0.0_dp
          der(2,1)=-1.0_dp
          der(2,2)= 0.0_dp
          der(2,3)= 1.0_dp
       CASE(6) 
          der(1,1)=-(4.0_dp*c3-1.0_dp) 
          der(1,2)=4.0_dp*c1-1.0_dp
          der(1,3)=0.0_dp
          der(1,4)=4.0_dp*(c3-c1)  
          der(1,5)=4.0_dp*c2   
          der(1,6)=-4.0_dp*c2
          der(2,1)=-(4.0_dp*c3-1.0_dp) 
          der(2,2)=0.0_dp
          der(2,3)=4.0_dp*c2-1.0_dp
          der(2,4)=-4.0_dp*c1      
          der(2,5)=4.0_dp*c1   
          der(2,6)=4.0_dp*(c3-c2)
       CASE(15)  !no le he cambiado aun 11/03/98                        
          t1=c1-.25_dp
          t2=c1-.5_dp
          t3=c1-.75_dp
          t4=c2-.25_dp
          t5=c2-.5_dp
          t6=c2-.75_dp
          t7=c3-.25_dp
          t8=c3-.5_dp
          t9=c3-.75_dp
          der(1,1)=32.0_dp/3.0_dp*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
          der(1,12)=128.0_dp/3.0_dp*c2*(t2*(t1+c1)+c1*t1) 
          der(1,11)=64.0_dp*c2*t4*(t1+c1)
          der(1,10)=128.0_dp/3.0_dp*c2*t4*t5  
          der(1,9)=0.0_dp
          der(1,8)=-128.0_dp/3.0_dp*c2*t4*t5
          der(1,7)=-64.0_dp*c2*t4*(t7+c3) 
          der(1,6)=-128.0_dp/3.0_dp*c2*(t8*(t7+c3)+c3*t7)
          der(1,5)=-32.0_dp/3.0_dp*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
          der(1,4)=128.0_dp/3.0_dp*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
          der(1,3)=64.0_dp*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
          der(1,2)=128.0_dp/3.0_dp*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
          der(1,13)=128.0_dp*c2*(c3*(t1+c1)-c1*t1) 
          der(1,15)=128.0_dp*c2*t4*(c3-c1)
          der(1,14)=128.0_dp*c2*(c3*t7-c1*(t7+c3))
          der(2,1)=0.0_dp
          der(2,12)=128.0_dp/3.0_dp*c1*t1*t2
          der(2,11)=64.0_dp*c1*t1*(t4+c2)
          der(2,10)=128.0_dp/3.0_dp*c1*(t5*(t4+c2)+c2*t4)
          der(2,9)=32.0_dp/3.0_dp*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
          der(2,8)=128.0_dp/3.0_dp*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
          der(2,7)=64.0_dp*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
          der(2,6)=128.0_dp/3.0_dp*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
          der(2,5)=-32.0_dp/3.0_dp*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
          der(2,4)=-128.0_dp/3.0_dp*c1*(t8*(t7+c3)+c3*t7)
          der(2,3)=-64.0_dp*c1*t1*(t7+c3)  
          der(2,2)=-128.0_dp/3.0_dp*c1*t1*t2
          der(2,13)=128.0_dp*c1*t1*(c3-c2)
          der(2,15)=128.0_dp*c1*(c3*(t4+c2)-c2*t4)
          der(2,14)=128.0_dp*c1*(c3*t7-c2*(c3+t7))        
       CASE (4)                                                              
          der(1,1)=-etam
          der(1,2)=etam 
          der(1,3)=etap
          der(1,4)=-etap
          der(2,1)=-xim 
          der(2,2)=-xip 
          der(2,3)=xip 
          der(2,4)=xim
       CASE(8)
          der(1,1)=etam*(2.0_dp*xi+eta)
          der(1,2)=etam*(2.0_dp*xi-eta)
          der(1,3)=etap*(2.0_dp*xi+eta)
          der(1,4)=etap*(2.0_dp*xi-eta)
          der(1,5)=-4.0_dp*etam*xi     
          der(1,6)=8.0_dp*etap*etam
          der(1,7)=-4.0_dp*etap*xi     
          der(1,8)=-8.0_dp*etam*etap
          der(2,1)=xim*(xi+2.*eta) 
          der(2,2)=xip*(2.0_dp*eta-xi)
          der(2,3)=xip*(xi+2.0_dp*eta) 
          der(2,4)=xim*(2.0_dp*eta-xi)
          der(2,5)=-8.0_dp*xim*xip     
          der(2,6)=-4.0_dp*xip*eta
          der(2,7)=8.0_dp*xim*xip      
          der(2,8)=-4.0_dp*xim*eta
       CASE(9)
          etam = eta - 1.0_dp
          etap = eta + 1.0_dp
          xim = xi - 1.0_dp
          xip = xi + 1.0_dp
          der(1,1)=.25_dp*x2m1*eta*etam  
          der(1,2)=.25*x2p1*eta*etam
          der(1,3)=.25_dp*x2p1*eta*etap  
          der(1,4)=.25_dp*x2m1*eta*etap
          der(1,5)=-xi*eta*etam       
          der(1,6)=-.5_dp*x2p1*etap*etam
          der(1,7)=-xi*eta*etap       
          der(1,8)=-.5_dp*x2m1*etap*etam
          der(1,9)=2.0_dp*xi*etap*etam    
          der(2,1)=.25_dp*xi*xim*e2m1
          der(2,2)=.25_dp*xi*xip*e2m1    
          der(2,3)=.25_dp*xi*xip*e2p1
          der(2,4)=.25_dp*xi*xim*e2p1    
          der(2,5)=-.5_dp*xip*xim*e2m1
          der(2,6)=-xi*xip*eta        
          der(2,7)=-.5_dp*xip*xim*e2p1
          der(2,8)=-xi*xim*eta        
          der(2,9)=2.0_dp*xip*xim*eta
       CASE default
          PRINT*,"wrong number of nodes in shape_der"        
       END SELECT
    CASE(3)  ! three dimensional elements (no la he modificado 11/03/98)
       xi    = points(i,1)
       eta   = points(i,2)
       zeta  = points(i,3)
       etam  = 1.0_dp-eta 
       xim   = 1.0_dp-xi
       zetam = 1.0_dp-zeta
       etap  = eta+1.0_dp 
       xip   = xi+1.0_dp 
       zetap = zeta+1.0_dp
       SELECT CASE (nod)
       CASE(4)
!!$          der(1:3,1:4) = 0.0_dp
!!$          der(1,1)     = 1.0_dp
!!$          der(2,2)     = 1.0_dp  
!!$          der(3,3)     = 1.0_dp
!!$          der(1,4)     =-1.0_dp 
!!$          der(2,4)     =-1.0_dp 
!!$          der(3,4)     =-1.0_dp  
          der          =  0.0_dp
          der(1,1)     = -1.0_dp
          der(1,2)     =  1.0_dp
          der(2,1)     = -1.0_dp
          der(2,3)     =  1.0_dp  
          der(3,1)     = -1.0_dp
          der(3,4)     =  1.0_dp
          !
       CASE(8)
          der(1,1)=-.125_dp*etam*zetam    
          der(1,2)=-.125_dp*etam*zetap
          der(1,3)=.125_dp*etam*zetap     
          der(1,4)=.125_dp*etam*zetam
          der(1,5)=-.125_dp*etap*zetam    
          der(1,6)=-.125_dp*etap*zetap
          der(1,7)=.125_dp*etap*zetap     
          der(1,8)=.125_dp*etap*zetam
          der(2,1)=-.125_dp*xim*zetam     
          der(2,2)=-.125_dp*xim*zetap
          der(2,3)=-.125_dp*xip*zetap     
          der(2,4)=-.125_dp*xip*zetam
          der(2,5)=.125_dp*xim*zetam      
          der(2,6)=.125_dp*xim*zetap
          der(2,7)=.125_dp*xip*zetap      
          der(2,8)=.125_dp*xip*zetam
          der(3,1)=-.125_dp*xim*etam      
          der(3,2)=.125_dp*xim*etam
          der(3,3)=.125_dp*xip*etam       
          der(3,4)=-.125_dp*xip*etam
          der(3,5)=-.125_dp*xim*etap      
          der(3,6)=.125_dp*xim*etap
          der(3,7)=.125_dp*xip*etap       
          der(3,8)=-.125_dp*xip*etap  
       CASE(14) ! type 6 element
          x= points(i,1)    
          y= points(i,2)  
          z= points(i,3) 
          der(1,1)=((2.0_dp*x*y+2.0_dp*x*z+4.0_dp*x+y*z+y+z)*(y-1.0_dp)*(z-1.0_dp))/8.0_dp
          der(1,2)=((2.0_dp*x*y-2.0_dp*x*z-4.0_dp*x+y*z+y-z)*(y+1.0_dp)*(z-1.0_dp))/8.0_dp
          der(1,3)=((2.0_dp*x*y+2.0_dp*x*z+4.0_dp*x-y*z-y-z)*(y-1.0_dp)*(z-1.0_dp))/8.0_dp
          der(1,4)=((2.0_dp*x*y-2.0_dp*x*z-4.0_dp*x-y*z-y+z)*(y+1.0_dp)*(z-1.0_dp))/8.0_dp
          der(1,5)=-((2.0_dp*x*y-2.0_dp*x*z+4.0_dp*x-y*z+y-z)*(y-1.0_dp)*(z+1.0_dp))/8.0_dp
          der(1,6)=-((2.0_dp*x*y+2.0_dp*x*z-4.0_dp*x-y*z+y+z)*(y+1.0_dp)*(z+1.0_dp))/8.0_dp
          der(1,7)=-((2.0_dp*x*y-2.0_dp*x*z+4.0_dp*x+y*z-y+z)*(y-1.0_dp)*(z+1.0_dp))/8.0_dp
          der(1,8)=-((2.0_dp*x*y+2.0_dp*x*z-4.0_dp*x+y*z-y-z)*(y+1.0_dp)*(z+1.0_dp))/8.0_dp
          der(1,9)=-(y+1.0_dp)*(y-1.0_dp)*(z-1.0_dp)*x  
          der(1,10)=(y+1.0_dp)*(y-1.0_dp)*(z+1.0_dp)*x
          der(1,11)=-(y-1.0_dp)*(z+1.0_dp)*(z-1.0_dp)*x 
          der(1,12)=(y+1.0_dp)*(z+1.0_dp)*(z-1.0_dp)*x
          der(1,13)=-((y+1.0_dp)*(y-1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp
          der(1,14)=((y+1.0_dp)*(y-1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp
          der(2,1)=((2.0_dp*x*y+x*z+x+2.0_dp*y*z+4.0_dp*y+z)*(x-1.0_dp)*(z-1.0_dp))/8.0_dp
          der(2,2)=((2.0_dp*x*y-x*z-x+2.0_dp*y*z+4.0_dp*y-z)*(x-1.0_dp)*(z-1.0_dp))/8.0_dp
          der(2,3)=((2.0_dp*x*y+x*z+x-2.0_dp*y*z-4.0_dp*y-z)*(x+1.0_dp)*(z-1.0_dp))/8.0_dp
          der(2,4)=((2.0_dp*x*y-x*z-x-2.0_dp*y*z-4.0_dp*y+z)*(x+1.0_dp)*(z-1.0_dp))/8.0_dp
          der(2,5)=-((2.0_dp*x*y-x*z+x-2.0_dp*y*z+4.0_dp*y-z)*(x-1.0_dp)*(z+1.0_dp))/8.0_dp
          der(2,6)=-((2.0_dp*x*y+x*z-x-2.0_dp*y*z+4.0_dp*y+z)*(x-1.0_dp)*(z+1.0_dp))/8.0_dp
          der(2,7)=-((2.0_dp*x*y-x*z+x+2.0_dp*y*z-4.0_dp*y+z)*(x+1.0_dp)*(z+1.0_dp))/8.0_dp
          der(2,8)=-((2.0_dp*x*y+x*z-x+2.0_dp*y*z-4.0_dp*y-z)*(x+1.0_dp)*(z+1.0_dp))/8.0_dp
          der(2,9)=-(x+1.0_dp)*(x-1.0_dp)*(z-1.0_dp)*y
          der(2,10)=(x+1.0_dp)*(x-1.0_dp)*(z+1.0_dp)*y
          der(2,11)=-((x+1.0_dp)*(x-1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp
          der(2,12)=((x+1.0_dp)*(x-1.0_dp)*(z+1.0_dp)*(z-1.0_dp))/2.0_dp
          der(2,13)=-(x-1.0_dp)*(z+1.0_dp)*(z-1.0_dp)*y
          der(2,14)=(x+1.0_dp)*(z+1.0_dp)*(z-1.0_dp)*y
          der(3,1)=((x*y+2.0_dp*x*z+x+2.0_dp*y*z+y+4.0_dp*z)*(x-1.0_dp)*(y-1.0_dp))/8.0_dp
          der(3,2)=((x*y-2.0_dp*x*z-x+2.0_dp*y*z+y-4.0_dp*z)*(x-1.0_dp)*(y+1.0_dp))/8.0_dp
          der(3,3)=((x*y+2.0_dp*x*z+x-2.0_dp*y*z-y-4.0_dp*z)*(x+1.0_dp)*(y-1.0_dp))/8.0_dp
          der(3,4)=((x*y-2.0_dp*x*z-x-2.0_dp*y*z-y+4.0_dp*z)*(x+1.0_dp)*(y+1.0_dp))/8.0_dp
          der(3,5)=-((x*y-2.0_dp*x*z+x-2.0_dp*y*z+y-4.0_dp*z)*(x-1.0_dp)*(y-1.0_dp))/8.0_dp
          der(3,6)=-((x*y+2.0_dp*x*z-x-2.0_dp*y*z+y+4.0_dp*z)*(x-1.0_dp)*(y+1.0_dp))/8.0_dp
          der(3,7)=-((x*y-2.0_dp*x*z+x+2.0_dp*y*z-y+4.0_dp*z)*(x+1.0_dp)*(y-1.0_dp))/8.0_dp
          der(3,8)=-((x*y+2.0_dp*x*z-x+2.0_dp*y*z-y-4.0_dp*z)*(x+1.0_dp)*(y+1.0_dp))/8.0_dp
          der(3,9)=-((x+1.0_dp)*(x-1.0_dp)*(y+1.0_dp)*(y-1.0_dp))/2.0_dp
          der(3,10)=((x+1.0_dp)*(x-1.0_dp)*(y+1.0_dp)*(y-1.0_dp))/2.0_dp
          der(3,11)=-(x+1.0_dp)*(x-1.0_dp)*(y-1.0_dp)*z  
          der(3,12)=(x+1.0_dp)*(x-1.0_dp)*(y+1.0_dp)*z
          der(3,13)=-(x-1.0_dp)*(y+1.0_dp)*(y-1.0_dp)*z  
          der(3,14)=(x+1.0_dp)*(y+1.0_dp)*(y-1.0_dp)*z
       CASE(20)
          xii=(/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
          etai=(/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1/)
          zetai=(/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
          DO l=1,20
             xi0=xi*xii(l)
             eta0=eta*etai(l)
             zeta0=zeta*zetai(l)
             IF(l==4.OR.l==8.OR.l==16.OR.l==20) THEN
                der(1,l)=-.5_dp*xi*(1.0_dp+eta0)*(1.0_dp+zeta0)
                der(2,l)=.25_dp*etai(l)*(1.0_dp-xi*xi)*(1.0_dp+zeta0)
                der(3,l)=.25_dp*zetai(l)*(1.0_dp-xi*xi)*(1.0_dp+eta0)
             ELSE IF(l>=9.AND.l<=12)THEN
                der(1,l)=.25_dp*xii(l)*(1.0_dp-eta*eta)*(1.0_dp+zeta0)
                der(2,l)=-.5_dp*eta*(1.0_dp+xi0)*(1.0_dp+zeta0)
                der(3,l)=.25_dp*zetai(l)*(1.0_dp+xi0)*(1.0_dp-eta*eta)
             ELSE IF(l==2.OR.l==6.OR.l==14.OR.l==18) THEN
                der(1,l)=.25_dp*xii(l)*(1.0_dp+eta0)*(1.0_dp-zeta*zeta)
                der(2,l)=.25_dp*etai(l)*(1.0_dp+xi0)*(1.0_dp-zeta*zeta)
                der(3,l)=-.5_dp*zeta*(1.0_dp+xi0)*(1.0_dp+eta0)
             ELSE
                der(1,l)=.125_dp*xii(l)*(1.0_dp+eta0)*(1.0_dp+zeta0)*(2.0_dp*xi0+eta0+&
                     zeta0-1.0_dp)
                der(2,l)=.125_dp*etai(l)*(1.0_dp+xi0)*(1.0_dp+zeta0)*(xi0+2.0_dp*eta0+&
                     zeta0-1.0_dp)
                der(3,l)=.125_dp*zetai(l)*(1.0_dp+xi0)*(1.0_dp+eta0)*(xi0+eta0+&
                     2.0_dp*zeta0-1.0_dp)
             END IF
          END DO
       CASE default
          PRINT*,"wrong number of nodes in shape_der"        
       END SELECT
    CASE default
       PRINT*,"wrong number of dimensions in shape_der"
    END SELECT
    RETURN
  END SUBROUTINE shape_der



END MODULE quadrature
