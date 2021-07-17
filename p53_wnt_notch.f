C	Runge Kutta method for numerical integration.     	

	implicit real*8(a-h,o-z)
	parameter(N=31)
      	dimension x(N),dxdt(N),xout(N)
      	external derivs,rk4
	open(unit=4,file='k42=0.3mdm2.dat',status='unknown')
	
c     initialization

      	write(*,*)"Enter initial time"
     	read(*,*) to
      	write(*,*)"Enter total time to be integrated"
      	read(*,*) tmax
      	write(*,*)"Enter step size"
      	read(*,*) h

      	nMax=int((tmax-to)/h)
c
        istep=100
c     initialization of molecular concentrations
	x(1)=0.01d0!Axin
	x(2)=0.02d0!mRNA_Axin
	x(3)=0.02d0!Unphosphorylated beta catenine
	x(4)=0.01d0!phosphorylated beta catenine
	x(5)=0.008d0!Nuclear beta catenine
	x(6)=2.0d0	!GSK3
	x(7)=10.0d0	!p53
	x(8)=10.0d0	!mdm2
	x(9)=0.0d0	!mdm2_Rna
	x(10)=95.0d0	!mdm2_p53
	x(11)=5.0d0	!nutilin
	x(12)=0.0d0	!mdm2_nutilin
	x(13)=0.0d0	!p53_Gsk3 complex
	x(14)=0.426d0!Dsh_activated
	x(15)=2.0d0!D/A
	x(16)=2.0d0!G/A
	x(17)=0.1d0!Lef
	x(18)=0.01d0!Lef_beta
	x(19)=1.68d0!delta
	x(20)=0.1d0!mRNA Delta cytoplasama
	x(21)=0.1d0!mRNA Delta nuclear
	x(22)=0.33d0!mRNA-Axinc
	x(23)=0.12d0!Notch protein
	x(24)=0.1d0!Notch protein(Cytosolic)
	x(25)=0.15d0!Notch Protein(Nucear)
	x(26)=0.71d0!mRNA_lunatic frienge
	x(27)=0.0071d0!Lunatic Frienge protein
	x(28)=0.2d0!N/D
c
c     Numerical integration of ODEs
	t=to
      	do nc=1,nMax
		t=t+h
c        	Numerical integration
         	call derivs(t,x,dxdt)
         	call rk4(x,dxdt,N,t,h,xout,derivs)
         	do j=1,N
            		x(j)=xout(j)
         	enddo
c		tm=t/60.d0
                if(t.ge.0.and.mod(nc,istep).eq.0)then
		write(*,*)t
		write(4,*)x(8)

        endif
	
                                
      	enddo
      	stop
      	end
c	-------------------------------
      	subroutine derivs(t,x,dxdt)
      	implicit real*8(a-h,o-z)
      	parameter(N=31)
      	dimension x(N),dxdt(N)
	external random

c	Constants
	ck1=1.80d0	
        ck2=0.10d0
	ck3=0.087d0
	ck4=0.70d0
	ck5=1.50d0	
        ck6=5.080d0
	ck7=2.0d0
	ck8=3.0d0									 
	ck9=0.5d0
	ck10=0.28d0
	ck11=1.0d0
	ck12=0.03d0
	ck13=0.0d0
	ck14=7.062d0
	ck15=0.06d0
	ck16=1.64d0
	ck17=0.7d0
	ck18=0.8d0
	ck19=0.48d0
	ck20=0.5d0
	ck21=0.05d0
	ck22=0.02d0
	ck23=0.6d0
        ck24=0.63d0
	ck25=2.0d0
	ck26=1.5d0
	ck27=ck8-x(6)
	ck28=-ck1*x(1)*x(6)+ck2*ck27
	ck29=-ck4*x(3)+ck5*x(5)
	ck30=(ck6*ck9*x(3))/((ck9+ck7)*(ck10+x(3)))
	ck31=(ck11*x(4))/(ck12+x(4))
	ck32=x(5)**ck25
	ck33=ck17**ck25
	ck34=0.000495d0*60.d0		!mdm2 formation
        ck35=0.0001d0*60.d0		!mdm2_Rna
	ck36=0.0001d0*60.d0		!mdm2_Rna
	ck37=0.000433d0*60.d0		!mdm2degradation
	ck38=0.078d0*60.d0		!p53 formation
        ck39=0.00825d0*60.d0		!p53 degradation
	ck40=0.001155d0*60.d0		!p53_mdm2 form
	ck41=0.0001155d0*60.d0		!p53_mdm2 defradation					 
	ck42=0.3d0*60.d0		!nutilin formation0.078d0
	ck43=0.002d0*60.d0		!nutilin_mdm2 formation0.001166d0	
	ck44=0.0005d0*60.d0		!nutilin_mdm2 degradation0.00001166d0
	ck45=0.001d0*60.d0		!nutilin degaradation0.0008d0
	ck46=0.04d0			!p53_Gsk3 complex formation-0.00012,0.000005d0
	ck47=0.12d0			!p53_Gsk3 complex degradation--0.12,0.000012d0
	ck48=0.042d0			!MDM2 Synthesis	
	ck49=5.0d0!vMaDsh
	ck50=1.814d0!wnt
	ck51=1.5d0!kawnt
	ck52=3.33d0!total Dsh
	ck53=ck52-x(14)-x(15)
	ck54=0.95d0!kaDsh
	ck55=1.9d0!VinaDsh
	ck56=1.647d0!Kinadsh
	ck57=0.05d0! k D/A
	ck58=8.0d0!a D/A
	ck59=(ck51+ck50)*(ck54+ck53)	
	ck60=0.9d0!a-B/L
	ck61=2.5d0!d-B/L
	ck62=0.01d0!Kdelta						
	ck63=0.8d0!Vmddelta
	ck64=2.385d0!Kmddelta
	ck65=1.725d0!EXmDeltan						
	ck66=2.0d0!Vmd-mdeltac
	ck67=0.48d0!Kmdeltac	 
	ck68=1.0d0!Vmdelta						
	ck69=1.12d0!VMsmdelta
	ck70=x(18)**ck25
	ck71=ck75**ck25
	ck72=1.0d0!EXmAxin2n
	ck73=1.2d0!VmdmAxin2c
	ck74=1.48d0!KmdmAxin2
	ck75=0.24d0!ksmdelta1
	ck76=0.5d0!ksmaxin2
	ck77=ck76**ck25
	ck78=0.00001d0	
        ck79=2.82d0
	ck80=0.5d0
	ck81=3.45d0
	ck82=0.010d0	
        ck83=0.001d0
	ck84=0.10d0
	ck85=0.001d0									 
	ck86=0.5d0
	ck87=0.1d0
	ck88=0.1d0
	ck89=3.0d0
	ck90=0.05d0
	ck91=1.92d0
	ck92=0.768d0
	ck93=0.3d0
	ck94=0.39d0
	ck95=0.37d0
	ck96=2.5d0
	ck97=2.0d0
	ck98=0.3d0
	ck99=ck87*x(24)-ck88*x(25)
	ck100=ck86**ck25
        ck101=x(27)**ck25
	ck102=x(25)**ck25
	ck103=ck90**ck25
	ck104=(ck79*x(23))/(ck80+x(23))
	ck105=(ck81*x(23)*ck100)/(ck100+ck101)
	ck106=(ck82*x(24))/(ck83+x(24))
	ck107=(ck84*x(25))/(ck85+x(25))
	ck108=(ck89*ck102)/(ck103+ck102)
	ck109=ck91*x(26)/(ck92+x(26))
	ck110=ck94*x(27)/(ck95+x(27))
	ck111=9.0d0!aN/D
	ck112=0.5d0!dN/D
	CK113=x(19)/ck80+x(19)
c										 
c       Deterministic ODEs
        dxdt(1)=ck26*(ck22*x(2)-(ck23*x(1))/(ck24+x(1))+ck28)
c     &   -ck58*x(1)*x(14)
	dxdt(2)=ck26*(ck15+(ck16*ck32)/(ck33+ck32)-
     &  (ck18*x(2))/(ck19+x(2)))+(ck16*ck70)/(ck77+ck70)
	dxdt(3)=ck26*(ck3-(ck30*ck27)/ck8+ck31+ck29-ck13*x(3))
	dxdt(4)=ck26*((ck30*ck27)/ck8-ck31-ck14*x(4))
	dxdt(5)=-ck26*ck29+ck61*x(18)-ck60*x(5)*x(17)
	dxdt(6)=ck26*ck28-ck46*x(7)*x(6)+ck47*x(13)
	dxdt(7)=ck38-ck40*x(7)*x(8)+ck41*x(10)-ck46*x(7)*x(6)+ck47*x(13)
	dxdt(8)=ck34*x(9)-ck37*x(8)+ck39*x(10)-ck40*x(7)*x(8)+ck41*x(10)
     &	-ck43*x(11)*x(8)
     	dxdt(9)=ck35*x(7)-ck36*x(9)+ck48*x(13)
	dxdt(10)=-ck39*x(10)+ck40*x(7)*x(8)-ck41*x(10)
	dxdt(11)=ck42-ck43*x(11)*x(8)+ck44*x(12)-ck45*x(11)
	dxdt(12)=ck43*x(11)*x(8)-ck44*x(12)
	dxdt(13)=ck46*x(7)*x(6)-ck47*x(13)
	dxdt(14)=(ck49*ck50*ck53)/ck59+ck57*x(15)-(ck55*x(14))/
     &   (ck56+x(14))-ck58*x(1)*x(14)-ck111*x(24)*x(14)+ck112*x(28)!new
	dxdt(15)=ck58*x(1)*x(14)-ck9*x(15)
	dxdt(16)=ck1*x(1)*x(6)-ck2*x(16)
	dxdt(17)=ck61*x(18)-ck60*x(5)*x(17)
	dxdt(18)=ck60*x(5)*x(17)-ck61*x(18)
	dxdt(19)=ck62*x(20)-(ck63*x(19))/(ck64+x(19))
	dxdt(20)=ck65*x(21)-ck66*x(20)/(ck67+x(20))
	dxdt(21)=ck68+(ck69*ck70)/(ck71+ck70)-ck65*x(21)
	dxdt(22)=ck72*x(2)-ck73*x(22)/(ck74+x(22))
	dxdt(23)=ck98*ck78-ck98*ck104-ck98*ck105*ck113
	dxdt(24)=ck98*ck105*ck113-ck98*ck106-ck98*ck99
     &  -ck111*x(24)*x(14)+ck112*x(28)!new
	dxdt(25)=ck98*ck99-ck98*ck107
	dxdt(26)=ck98*ck108-ck98*ck109
	dxdt(27)=ck98*ck93*x(26)-ck98*ck110
	dxdt(28)=ck111*x(24)*x(14)-ck112*x(28)!new
c
	return
      	end
c	-----------------------------------------
      	subroutine rk4(y,dydx,n,x,h,yout,derivs)
	implicit real*8(a-h,o-z)
	parameter (NMAX=123)
	dimension dydx(n),y(n),yout(n)
	dimension dym(NMAX),dyt(NMAX),yt(NMAX)
c       integer n,NMAX
c       real h,x,dydx(n),y(n),yout(n)
	external derivs
c       parameter (NMAX=50)
c       integer i
c       real h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
	hh=h*0.5d0
	h6=h/6.0d0
	xh=x+hh
	do 11 i=1,n
	   yt(i)=y(i)+hh*dydx(i)
 11	continue
	call derivs(xh,yt,dyt)
	do 12 i=1,n
	   yt(i)=y(i)+hh*dyt(i)
 12	continue
	call derivs(xh,yt,dym)
	do 13 i=1,n
	   yt(i)=y(i)+h*dym(i)
	   dym(i)=dyt(i)+dym(i)
 13	continue
	call derivs(x+h,yt,dyt)
	do 14 i=1,n
	   yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
 14	continue
	return
	end

