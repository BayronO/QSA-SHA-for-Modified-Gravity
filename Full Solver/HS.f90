Module Cosmology
  IMPLICIT none
  Real *8, parameter :: aini = 10**(-3.d0)
  Real *8, parameter :: k = 300.d0
  Real *8, parameter :: b = 9.89557132*10**(-4.d0)
  Real *8, parameter :: om = 0.309883d0
  Real *8, parameter :: eps = 10**(-10.d0) ! GR deviation tolerance 
End Module Cosmology

PROGRAM fRmodels
  Use Cosmology
  use, intrinsic :: iso_c_binding
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
  
  write(*,*)'H(aini),Hp(aini)'
  write(*,*)H(aini),Hp(aini)
  write(*,*)'F(aini),Fp(aini)'
  write(*,*)F(aini),Fp(aini)
  write(*,*)'dmini,Vmini'
  write(*,*)dmini,Vmini
  write(*,*)'phiplini,chiini'
  write(*,*)phiplini,chiini

  write(*,*)'Now solving ODE'
  open(2,file='HS_k_300.txt')
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
		dF=abs(F(a)-1.0d0)
		write(*,*)'dF= ',dF
		tol=0.01d0/i**3 !works for b = k = 300
		if(dF.le.eps) then
		call dverk(4,fderivGR,ainitial,y,a,tol,ind,c,4,work)
		else
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
    H = sqrt(1-om+om/a**3+(2*a**2*b*(-1+om)**2*(12*a**7*(-1+om)**2+3*a**4*(-1+om)*om-6*a*om**2))/(4*a**3*(-1+om)-om)**3-(a**5*b**2*(-1+om)**3*(-1024*a**19*(-1+om)**6-9216*a**16*(-1+om)**5*om+22848*a**13*(-1+om)**4*om**2-25408*a**10*(-1+om)**3*om**3+7452*a**7*(-1+om)**2*om**4+4656*a**4*(-1+om)*om**5-37*a*om**6))/(-4*a**3*(-1+om)+om)**8)
    RETURN
 END FUNCTION

FUNCTION Hp(a)
    Use Cosmology
    Implicit none
	Real *8 :: H
    Real *8 :: Hp
	Real *8 :: a
    Hp = ((-3*om)/a**4 + (2*a**2*b*(-1 + om)**2*(84*a**6*(-1 + om)**2 + 12*a**3*(-1 + om)*om - 6*om**2))/(4*a**3*(-1 + om) - om)**3 - (216*a**5*b*(-1 + om)**3*(4*a**6*(-1 + om)**2 + a**3*(-1 + om)*om - 2*om**2))/(-4*a**3*(-1 + om) + om)**4 + (4*a*b*(-1 + om)**2*(12*a**7*(-1 + om)**2 + 3*a**4*(-1 + om)*om - 6*a*om**2))/(4*a**3*(-1 + om) - om)**3 - (a**5*b**2*(-1 + om)**3*(-19456*a**18*(-1 + om)**6 - 147456*a**15*(-1 + om)**5*om + 297024*a**12*(-1 + om)**4*om**2 - 254080*a**9*(-1 + om)**3*om**3 + 52164*a**6*(-1 + om)**2*om**4 + 18624*a**3*(-1 + om)*om**5 - 37*om**6))/(-4*a**3*(-1 + om) + om)**8 - (96*a**7*b**2*(-1 + om)**4*(-1024*a**19*(-1 + om)**6 - 9216*a**16*(-1 + om)**5*om + 22848*a**13*(-1 + om)**4*om**2 - 25408*a**10*(-1 + om)**3*om**3 + 7452*a**7*(-1 + om)**2*om**4 + 4656*a**4*(-1 + om)*om**5 - 37*a*om**6))/(-4*a**3*(-1 + om) + om)**9 - (5*a**4*b**2*(-1 + om)**3*(-1024*a**19*(-1 + om)**6 - 9216*a**16*(-1 + om)**5*om + 22848*a**13*(-1 + om)**4*om**2 - 25408*a**10*(-1 + om)**3*om**3 + 7452*a**7*(-1 + om)**2*om**4 + 4656*a**4*(-1 + om)*om**5 - 37*a*om**6))/(-4*a**3*(-1 + om) + om)**8)/(2.d0*H(a))
    RETURN
 END FUNCTION
 
 FUNCTION Hgr(a)
    Use Cosmology
    Implicit none
    Real *8 :: Hgr
	Real *8 :: a
    Hgr = sqrt(om*a**(-3.d0)+1.d0-om)
    RETURN
 END FUNCTION

FUNCTION Hpgr(a)
    Use Cosmology
    Implicit none
    Real *8 :: Hpgr
	Real *8 :: a
    Hpgr = (-3.d0*om)/(2.d0*a**4*sqrt(1.d0 - om + om/a**3.d0))
    RETURN
 END FUNCTION
 
FUNCTION F(a)
    Use Cosmology
    Implicit none
    Real *8 :: F
	Real *8 :: a
    F = 1 - (2*a**6*b*(-1 + om)**2)/(-4*a**3*(-1 + om) + om)**2 - (4*a**9*b**2*(-1 + om)**3*(128*a**12*(-1 + om)**4 - 32*a**9*(-1 + om)**3*om - 60*a**6*(-1 + om)**2*om**2 + 100*a**3*(-1 + om)*om**3 - om**4))/(4*a**3*(-1 + om) - om)**7
    RETURN
 END FUNCTION
 
FUNCTION Fp(a)
    Use Cosmology
    Implicit none
    Real *8 :: Fp
	Real *8 :: a
	Fp = (12*a**5*b*(-1 + om)**2*om)/(4*a**3*(-1 + om) - om)**3 + (36*a**8*b**2*(-1 + om)**3*om*(256*a**12*(-1 + om)**4 - 224*a**9*(-1 + om)**3*om + 300*a**6*(-1 + om)**2*om**2 + 128*a**3*(-1 + om)*om**3 - om**4))/(-4*a**3*(-1 + om) + om)**8
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
  Real *8 :: Hgr, Hpgr, Oma
  
  yprime(1)= -(y(2)/(a**2*Hgr(a)))+(3*a*Hpgr(a)*y(2))/k**2+(9*Oma(a)*y(2))/(2.*k**2*Hgr(a))
  
  yprime(2)= -(y(2)/a) + (k**2*y(3))/(a**2*Hgr(a))
  
  yprime(3)= (3*Oma(a)*y(2))/(2.*k**2*Hgr(a)) - y(3)/a
    
  yprime(4)= 0.d0
  return
END SUBROUTINE fderivGR
