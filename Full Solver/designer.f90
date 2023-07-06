Module Cosmology
  IMPLICIT none
  Real *8, parameter :: aini = 10**(-3.d0)
  Real *8, parameter :: k = 300.d0
  Real *8, parameter :: fR0 = -10**(-4.d0)
  Real *8, parameter :: om = 0.309883d0
  Real *8, parameter :: ndes = 3.386000936329383d0
  Real *8, parameter :: gom= 1.6452934638598856d0
  Real *8, parameter :: eps=10**(-9.d0) ! GR deviation tolerance  
End Module Cosmology

PROGRAM fRmodels
  Use Cosmology
  IMPLICIT none
  DOUBLE PRECISION :: a, ainitial, y(4), c(24), work(4,9), tol=1d-9
  INTEGER :: i,ind
  Real *8 :: H, Hp, F, Fp, Oma
  Real *8 :: dF
  Real *8 :: phiplini, chiini, dmini, Vmini
  EXTERNAL fderiv, fderivGR
 
! The initial conditions.
  phiplini=-1.d0
  chiini=0.d0
  dmini=-2.d0*k**2*phiplini/(3.d0*aini**2*H(aini)**2)
  Vmini=2.d0*k**2*phiplini/(3.d0*aini*H(aini))  
  
  write(*,*)'dmini,Vmini'
  write(*,*)dmini,Vmini
  write(*,*)'phiplini,chiini'
  write(*,*)phiplini,chiini

  write(*,*)'Now solving ODE'
  open(2,file='fDES_k_300.txt')
  ind=1
  do i=1,1000
     a = i*0.001d0
     if (i.eq.1) then 
        ainitial = a
        y(1) = dmini 
        y(2) = Vmini
        y(3) = phiplini
        y(4) = chiini
		write(*,*)'Steps: ',i,'/ 1000 done!'
     else
	    !tol=0.01d0/i**3 !works for fR0=-0.1 & fR0=-10^-4 k=300
	    !tol=0.1d0/i**3.5d0 !works for fR0=-0.1 & fR0=-10^-4 k=300
		!tol=0.001d0/i**3 !works for fR0=-10^-4 k=1500 sort of
		dF=abs(F(a)-1.0d0)
		write(*,*)'dF= ',dF
		if(dF.le.eps) then
		call dverk(4,fderivGR,ainitial,y,a,tol,ind,c,4,work)
		else
		tol=0.001d0/i**3
		call dverk(4,fderiv,ainitial,y,a,tol,ind,c,4,work)
		endif
        ainitial = a
		write(*,*)'Steps: ',i,'/ 1000 done!'
     endif 
	 write(2,*)a,y(1),y(2),y(3),y(4)
	enddo             
  close(2)
  stop
END PROGRAM
!-----------------------------------------------------------------------------
FUNCTION Oma(a)
    Use Cosmology
    Implicit none
    Real *8 :: Oma
	Real *8 :: a
    Oma = om*a**(-3.d0)
    RETURN
 END FUNCTION

FUNCTION H(a)
    Use Cosmology
    Implicit none
    Real *8 :: H
	Real *8 :: a
    H = sqrt(om*a**(-3.d0)+1.d0-om)
    RETURN
 END FUNCTION

FUNCTION Hp(a)
    Use Cosmology
    Implicit none
    Real *8 :: Hp
	Real *8 :: a
    Hp = (-3.d0*om)/(2.d0*a**4*sqrt(1.d0 - om + om/a**3.d0))
    RETURN
 END FUNCTION
 
FUNCTION F(a)
    Use Cosmology
    Implicit none
    Real *8 :: F
	Real *8 :: a
    F = 1.d0+gom*fR0*a**ndes
    RETURN
 END FUNCTION
 
FUNCTION Fp(a)
    Use Cosmology
    Implicit none
    Real *8 :: Fp
	Real *8 :: a
	Fp = ndes*gom*fR0*a**(ndes-1.d0)
    RETURN
 END FUNCTION

SUBROUTINE fderiv(n,a,y,yprime)
  Use Cosmology
  Implicit none 
  INTEGER, intent(IN) :: n ! # of derivatives
  DOUBLE PRECISION :: a,y(n),yprime(n)
  Real *8 :: H, Hp, F, Fp, Oma
  
  yprime(1)= -2.d0*k**2*y(3)*F(a)/(a**4*Fp(a)*H(a)**2)-y(2)/(a**2*H(a))+3.d0*a*y(2)*Hp(a)/k**2-3.d0*y(4)*Hp(a)/(a*Fp(a)*H(a))-3.d0*y(1)*Oma(a)/(a**2*Fp(a)*H(a)**2)
  
  yprime(2)= -y(2)/a+k**2*y(3)/(a**2*H(a))-y(4)*k**2/(2.d0*a**2*F(a)*H(a))
  
  yprime(3)= -y(3)/a+3.d0*y(4)*Fp(a)/(4.d0*F(a)**2)-y(3)*Fp(a)/(2.d0*F(a))+3.d0*y(2)*Oma(a)/(2.d0*k**2*F(a)*H(a))
    
  yprime(4)= y(4)/a+y(3)*Fp(a)-y(4)*Fp(a)/(2.d0*F(a))-4.d0*k**2*y(3)*F(a)**2/(3.d0*a**4*Fp(a)*H(a)**2)-2.d0*y(4)*F(a)*Hp(a)/(a*Fp(a)*H(a))-2.d0*y(1)*F(a)*Oma(a)/(a**2*Fp(a)*H(a)**2)-3.d0*y(2)*Oma(a)/(k**2*H(a))
  return
END SUBROUTINE fderiv

SUBROUTINE fderivGR(n,a,y,yprime)
  Use Cosmology
  Implicit none 
  INTEGER, intent(IN) :: n ! # of derivatives
  DOUBLE PRECISION :: a,y(n),yprime(n)
  Real *8 :: H, Hp, F, Fp, Oma
  
  yprime(1)= -(y(2)/(a**2*H(a)))+(3*a*Hp(a)*y(2))/k**2+(9*Oma(a)*y(2))/(2.*k**2*H(a))
  
  yprime(2)= -(y(2)/a) + (k**2*y(3))/(a**2*H(a))
  
  yprime(3)= (3*Oma(a)*y(2))/(2.*k**2*H(a)) - y(3)/a
    
  yprime(4)= 0.d0
  return
END SUBROUTINE fderivGR
