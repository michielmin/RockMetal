	subroutine RockMetal(f_core,f_h2o,Mp,Rp)
	IMPLICIT NONE
	integer nr
	parameter(nr=100)
	real*8 G,Mearth,Rearth,pi,EoS,f_core,f_h2o,Mp,m0
	parameter(G=6.67408e-11)! m3 kg-1 s-2
	parameter(Mearth=5.972e24)! kg
	parameter(Rearth=6371d3)! m
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(m0=1.660539040d-24*6.022136736d23)
	real*8 Mtot,Madd,rho,Rp,Rmin,Rmax,K0(10),K0t(10),rho0(10)
	real*8 Mm(10),M,Mlay,MW,V0,P1,R1,P2,R2
	integer nm,i,j,k,l,nlay
	logical last

!	Ice
c@ARTICLE{1993JChPh..99.5369F,
c       author = {{Fei}, Yingwei and {Mao}, Ho-Kwang and {Hemley}, Russell J.},
c        title = "{Thermal expansivity, bulk modulus, and melting curve of H$_{2}$O-ice VII to 20 GPa}",
c      journal = {\jcp},
c     keywords = {Ice, Melting, Microstructure, Thermal Expansion, Pressure Effects, Temperature Effects, X Ray Diffraction, Thermodynamics and Statistical Physics},
c         year = 1993,
c        month = oct,
c       volume = {99},
c       number = {7},
c        pages = {5369-5373},
c          doi = {10.1063/1.465980},
c       adsurl = {https://ui.adsabs.harvard.edu/abs/1993JChPh..99.5369F},
c      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
c}
	i=1
	K0(i)=24.1d0
	K0t(i)=4.1d0
	V0=12.3
	MW=17.8851
	rho0(i)=m0*MW/V0

c Mantle materials from:
c@ARTICLE{2005GeoJI.162..610S,
c       author = {{Stixrude}, Lars and {Lithgow-Bertelloni}, Carolina},
c        title = "{Thermodynamics of mantle minerals - I. Physical properties}",
c      journal = {Geophysical Journal International},
c     keywords = {bulk modulus, mantle, shear modulus, thermodynamics},
c         year = 2005,
c        month = aug,
c       volume = {162},
c       number = {2},
c        pages = {610-632},
c          doi = {10.1111/j.1365-246X.2005.02642.x},
c       adsurl = {https://ui.adsabs.harvard.edu/abs/2005GeoJI.162..610S},
c      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
c}
!	MgSiO3
	i=2
	K0(i)=211d0
	K0t(i)=4.5d0
	V0=26.35
	MW=100.39
	rho0(i)=m0*MW/V0

!	Fe
c@ARTICLE{2006PhRvL..97u5504D,
c       author = {{Dewaele}, Agn{\`e}s and {Loubeyre}, Paul and {Occelli}, Florent and {Mezouar}, Mohamed and {Dorogokupets}, Peter I. and {Torrent}, Marc},
c        title = "{Quasihydrostatic Equation of State of Iron above 2Mbar}",
c      journal = {\prl},
c     keywords = {62.50.+p, 07.35.+k, 64.30.+t, High-pressure and shock wave effects in solids and liquids, High-pressure apparatus, shock tubes, diamond anvil cells, Equations of state of specific substances},
c         year = 2006,
c        month = nov,
c       volume = {97},
c       number = {21},
c          eid = {215504},
c        pages = {215504},
c          doi = {10.1103/PhysRevLett.97.215504},
c       adsurl = {https://ui.adsabs.harvard.edu/abs/2006PhRvL..97u5504D},
c      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
c}
	i=3
	K0(i)=165d0
	K0t(i)=4.97d0
	V0=11.234*1e-24*6.022136736d23
	MW=55.840
	rho0(i)=m0*MW/V0

	nm=i

!	Case of 3 layers (H2O, MgSiO3, Fe)
	nlay=3

	Mm(1)=f_h2o
	Mm(2)=(1d0-f_h2o)*(1d0-f_core)
	Mm(3)=(1d0-f_h2o)*f_core

	M=Mp*Mearth

	Mm(1:nm)=Mm(1:nm)/sum(Mm(1:nm))

	Rmax=Rearth*(M/Mearth)**(1d0/3d0)
	Rmin=0d0*Rearth
	Rp=(Rmax+Rmin)/2d0

	last=.true.
	do while((Rmax-Rmin).gt.1d-3*Rp.or.last)
		Rp=(Rmax+Rmin)/2d0
		if((Rmax-Rmin).le.1d-3*Rp) then
			Rp=Rmax
			last=.false.
		endif
		R2=Rp
		P2=0d0
		Mtot=0d0
		do j=1,nlay
			Mlay=Mm(j)
			if(Mlay.gt.0d0) then
				do i=1,nr
					rho=EoS(P2,K0(j),K0t(j),rho0(j))*1d3
					Madd=Mlay*M/real(nr)
					R1=(R2**3-Madd*3d0/(4d0*pi*rho))**(1d0/3d0)
					Mtot=Mtot+Madd
					if(.not.R1.gt.0d0) then
						Rmin=Rp
						Rp=Rp*(M/Mtot)**(1d0/3d0)
						if(Rp.gt.Rmax) Rmax=Rp
						goto 10
					endif
					P1=P2+G*(M-Mtot)*rho*
     &			abs(R2-R1)/(0.5d0*(R1+R2))**2
					R2=R1
					P2=P1
				enddo
			endif
		enddo
		Rmax=Rp
10		continue
	enddo

	Rp=Rp/Rearth
	
	end



	real*8 function EoS(P,K0,K0t,rho0)
	IMPLICIT NONE
	real*8 P,K0,K0t,rho0,rhomin,rhomax,eps,P0,f,rho
	rhomin=0.5d0*rho0
	rhomax=5d0*rho0
	
	eps=1d0
	do while(eps.gt.1d-3)
		rho=(rhomin+rhomax)/2d0
		f=rho/rho0
		P0=(3d0*K0/2d0)*(f**(7d0/3d0)-f**(5d0/3d0))*
     &		(1d0-3d0*(4d0-K0t)*(f**(2d0/3d0)-1d0)/4d0)
		P0=P0*1d9
		if(P0.gt.P) then
			rhomax=rho
		else
			rhomin=rho
		endif
		eps=(rhomax-rhomin)/(rhomax+rhomin)
	enddo
	EoS=rho

	return
	end
	
