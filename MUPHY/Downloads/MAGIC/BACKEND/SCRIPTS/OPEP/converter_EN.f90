program converter


! convert OPEP input files into MUPHY input files for MD
! (optionally do Elastic Network)
! for the unit MUPHY will use exactly the same units of OPEP:

! dx = Angstrom = 10^-10 m
! dt = femtosec = 10^-15 s
! dE = kcal/mol 
! dm = kcal/mol * fs^2/Ang^2 = 2390 a.m.u.
! dT = dE / kB = 503.6 Kelvin

IMPLICIT NONE

CHARACTER(len=40) :: input

REAL :: nx,ny,nz

REAL,PARAMETER :: mass_conversion_factor = 2390.0 ! in opep units
REAL,PARAMETER :: timestep=1.d0 ! is fs

! REAL,PARAMETER :: mass_conversion_factor = 1.0 ! in amu
! REAL,PARAMETER :: timestep=0.001d0 ! in amu units....

REAL,PARAMETER :: PI=3.141592653589793d0

REAL :: score(272)

INTEGER,PARAMETER :: natom_max=5000000
INTEGER,PARAMETER :: ntypes_max=100000
INTEGER,PARAMETER :: nres_max=100000
INTEGER :: i,j,k,l,n,kk,m
INTEGER :: natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
           nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
           numbnd,numang,nptra,natyp,nphb,idum,ifbox,nmxrs,ifcap
INTEGER :: nttyp,ntype
CHARACTER(len=4) :: igraph(natom_max),labres(natom_max)
CHARACTER(len=4) :: pdbatm(natom_max)
CHARACTER(len=3) :: pdbres(natom_max)
INTEGER :: ipdbres(natom_max)
REAL :: chrg(natom_max),amass(natom_max)
INTEGER :: iac(natom_max)
INTEGER :: nno(ntypes_max),numex(natom_max)
INTEGER :: ipres(nres_max)
REAL,DIMENSION(natom_max) :: rk,req,tk,teq,pk,pn,phase
REAL,DIMENSION(natom_max) :: solty,cn1,cn2

INTEGER,DIMENSION(natom_max) :: ibh,jbh,icbh,ib,jb,icb
INTEGER,DIMENSION(natom_max) :: ith,jth,kth,icth,it,jt,kt,ict
INTEGER,DIMENSION(natom_max) :: iph,jph,kph,lph,icph,ip,jp,kp,lp,icp
INTEGER :: natex(ntypes_max)

CHARACTER(len=3) :: resp(nres_max)
INTEGER :: ialpha(nres_max), ibeta(nres_max)

INTEGER,PARAMETER :: nb_max=10000000
INTEGER :: nb
INTEGER,DIMENSION(nb_max) :: ni,ivi,nj,ivj,icoeff
REAL,DIMENSION(nb_max) :: rncoe,vamax,ct0lj ! ,ct2lj
REAL,DIMENSION(natom_max) :: x,y,z
REAL :: cm(3)
REAL,DIMENSION(20000) :: rx,ry,rz
REAL,DIMENSION(20) :: foal,fobe,walpha,wbeta,walpha_foal,wbeta_fobe
INTEGER :: n14,i14(natom_max,3)
INTEGER :: ia1,ia2,ibig,isml,ic,ic0,kdiv

REAL,PARAMETER :: scal14 = 1.476    ! taken from ETORS

INTEGER :: ipn(natom_max)
REAL :: gamc(natom_max),gams(natom_max),fmn(natom_max)
REAL :: epshb,epshb_mcmc(natom_max)

CHARACTER(len=100) :: string

CHARACTER(len=6) :: labatm(natom_max)

INTEGER :: nuniq=0,n_ca_sc=0,n_hb=0
CHARACTER(len=6) :: auniq(natom_max)
INTEGER :: iuniq(natom_max),tuniq(natom_max)

INTEGER :: nfrag, lenfrag(1000),ichain(natom_max)

INTEGER :: npar4b,ipar4b(4,4*natom_max)
REAL :: par4b(2,4*natom_max)

LOGICAL :: logic1i,logic21i,logic22,logic42,logic43,logic51i,logic58i,logic1,logic2,logic3, &
           logic6039,logic21,logic23,logic4,logic24,logic25,logic26,logic31,logic32,logic33,logic34, &
           logic5,logic41,logic6,logic51,logic52,logic7,logic8,logic53,logic54,logic55,logic56, &
           logic57,logic58,logic59,logic1333,logic9

INTEGER :: i1,i2,i3,i4,j1,j2,iabi,jcb1,icb1,icb3,iabk,ki31,ki42,kiabi,kiabj,icb2,icb4,ILIM
INTEGER :: i1b,i2b,i3b,i4b
REAL :: parah1,parah2,parah3,parah4,parah=0,parab=0

INTEGER :: nrama, irama(4,natom_max)
REAL :: prama(3,natom_max)
REAL :: c12,c6,xmass
REAL :: scaling_factor
CHARACTER(len=3) :: dums

CHARACTER(len=3),DIMENSION(2) :: crg_plus=(/ 'ARG','LYS' /),crg_minus=(/ 'ASP','GLU' /)
INTEGER,PARAMETER :: npair_max=10000
TYPE ip_infot
      INTEGER :: index_ip(npair_max,3)
      INTEGER :: index_ca(npair_max,2)
      INTEGER :: n_ip
END TYPE
TYPE(ip_infot) :: ip_info
INTEGER,DIMENSION(npair_max) :: index_ip_plus, itype_ip_plus, resindex_ip_plus
INTEGER,DIMENSION(npair_max) :: index_ip_minus, itype_ip_minus, resindex_ip_minus
INTEGER :: ipp,inn,i_cp,i_cn,ii,jj,ires,npair
INTEGER :: idonor,ihydrogen,iacceptor
REAL :: alfa1,alfa2,beta1,beta2
REAL :: e_wall=0.000125,s_wall=4.64
LOGICAL :: lionpair=.false.,lwall=.false.,lnet=.false.,lmass=.false.,lnet_H=.false.


!--- VARAIBLES EXTRACT Sc-Sc INTERMOLECULAR INTERACTIONS ---

!!integer :: i_sc,j_sc,i_c,i_sc_res,j_sc_res
CHARACTER(LEN=6),DIMENSION(17),SAVE :: sc_name = (/ "ARG","ASN","ASP","CYS","GLN","GLU","HIS","ILE","LEU","LYS","MET","PHE", &
     & "SER","THR","TRP","TYR","VAL"/)
INTEGER :: n_scsc,n_scsc_sys,nvdw
INTEGER,DIMENSION(1000) :: pair_sys
REAL,DIMENSION(1000,2) :: pot_data
CHARACTER(LEN=6),DIMENSION(1000) :: pair_pot_type
CHARACTER(LEN=6),DIMENSION(1000,2) ::pair
LOGICAL,DIMENSION(1000,2) :: ipair
CHARACTER(LEN=100) :: pathscsc

!--- Variables for Elastic Network --------------------------
REAL :: net_constant,d0,ddx,ddy,ddz,pair_d0(50000),pair_k(50000),net_rc
INTEGER :: pair_net(50000,2),natom_heavy,net_pair,inet,jnet,heavy(10000),ih
!------------------------------------------------------------

print *,'Insert the prefix to process:'
READ(*,*) input
input = ADJUSTL(input)
print *,input

print *,'Insert nx,ny,nz:'
READ(*,*) nx,ny,nz
print *,nx,ny,nz

print *,'Insert scaling factor (1.6 - 1.8):'
READ(*,*) scaling_factor
print *,scaling_factor

print *,'Insert ion pair flag (.true.,.false.):'
READ(*,*) lionpair
print *,lionpair

print *,'Insert the directory path where the scsc.dat file is found:'
READ(*,'(a)') pathscsc
print *,pathscsc

print *,'Repulsive Wall (.true.,.false.) [for shear....]'
READ(*,*) lwall
print *,lwall


print *, 'Elastic Model (.true.,.false.)'
READ(*,*) lnet
print *,lnet

IF(lnet) THEN
   print *,'Only H? (.true.,.false)'
   READ(*,*) lnet_H
   IF (.not.lnet_H) print *,'The EN will be based on Ca and SC'
ENDIF

IF(.not.lnet) THEN
   print *,'Rescale Mass Backbone (.true.,.flase.)'
   READ(*,*) lmass
   print *,lmass
ENDIF

IF(lnet) THEN
   print *,'Cut-Off for EN (suggested 5-8 A)'
   READ(*,*) net_rc 
   OPEN(unit=80,file='reference.pdb')
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(unit=40,file=TRIM(pathscsc)//'/'//'scsc.dat',STATUS='OLD')
i=0
DO
   i=i+1
   READ(40,*,END=15) pair(i,1),pair(i,2), pair_pot_type(i), pot_data(i,1),pot_data(i,2)
!   print *,'pair',pair(i,1),pair(i,2)
   IF(lionpair) THEN
      IF(pair(i,1)=="ARG".AND.pair(i,2)=="ASP") THEN
         pair_pot_type(i)="ip4b"
         pot_data(i,1)=1
         pot_data(i,2)=0
      ELSEIF(pair(i,1)=="ARG".AND.pair(i,2)=="GLU") THEN
         pair_pot_type(i)="ip4b"
         pot_data(i,1)=3
         pot_data(i,2)=0
      ELSEIF(pair(i,1)=="ASP".AND.pair(i,2)=="LYS") THEN
         pair_pot_type(i)="ip2b"
         pot_data(i,1)=2
         pot_data(i,2)=0
      ELSEIF(pair(i,1)=="GLU".AND.pair(i,2)=="LYS") THEN
         pair_pot_type(i)="ip2b"
         pot_data(i,1)=4
         pot_data(i,2)=0
      ENDIF
   ENDIF
ENDDO
15 CONTINUE
n_scsc=i




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



print *
print *,'................Reading file '//TRIM(input)//'.top.................'
OPEN(10,STATUS='old',FILE=TRIM(input)//'.top')
READ(10,*) 
READ(10,*) natom,ntypes,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
           nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
           numbnd,numang,nptra,natyp,nphb,idum,idum,idum,idum,idum, &
           idum,idum,ifbox,nmxrs,ifcap

nttyp = (ntypes*(ntypes+1))/2
ntype = ntypes**2


READ(10,9108) (igraph(i),i = 1,natom) ! atom name
READ(10,9128) (chrg(i),i = 1,natom) ! atom charge (=0 always)

IF(ALL(ABS(chrg)>1.d-6)) STOP 'some atoms are charged!'

READ(10,9128) (amass(i),i = 1,natom) ! atom mass

! convert from amu to OPEP mass units
! amass = amass * 2390.0 
amass = amass * mass_conversion_factor

READ(10,9118) (iac(i),i = 1,natom) ! atom type 
READ(10,9118) (numex(i),i = 1,natom) ! number of excluded atoms
READ(10,9118) (nno(i),i = 1,ntype) ! non-bond index
READ(10,9108) (labres(i),i=1,nres) ! residue name
READ(10,9118) (ipres(i),i=1,nres) ! id of 1st atom (N) in residue

READ(10,9128) (rk(i),    i = 1,numbnd) ! bonding force constants
READ(10,9128) (req(i),   i = 1,numbnd) ! bonding equilibrium distance
READ(10,9128) (tk(i),    i = 1,numang) ! angle force constant
READ(10,9128) (teq(i),   i = 1,numang) ! angle equilibrium constant
READ(10,9128) (pk(i),    i = 1,nptra) ! dihedral force constant
READ(10,9128) (pn(i),    i = 1,nptra) ! dihedral periodicity
READ(10,9128) (phase(i), i = 1,nptra) ! dihedtral phase
READ(10,9128) (solty(i), i = 1,natyp) ! improper dihedral 
READ(10,9128) (cn1(i),   i = 1,nttyp) ! LJ r12 coefficient
READ(10,9128) (cn2(i),   i = 1,nttyp) ! LJ r6 coefficient

print *,'number of bonds nbonh,nbona:',nbonh,nbona

print *,'number of bond angles ntheth,ntheta:',ntheth,ntheta

print *,'number of bond dihedrals nphih,nphia:',nphih,nphia

! bonds
READ(10,9118) (ibh(i),jbh(i),icbh(i),i = 1,nbonh)
READ(10,9118) (ib(i),jb(i),icb(i),i = 1,nbona)

! angles
READ(10,9118) (ith(i),jth(i),kth(i),icth(i), i=1, ntheth)
READ(10,9118) (it(i),jt(i),kt(i),ict(i), i=1, ntheta)

! dihedrals
READ(10,9118) (iph(i),jph(i),kph(i),lph(i),icph(i),i = 1,nphih)
READ(10,9118) (ip(i),jp(i),kp(i),lp(i),icp(i),i=1,nphia)

READ(10,9118) (natex(i),i=1,NNB)

CLOSE(10)

9108 FORMAT(20a4)
9118 FORMAT(12i6)
9128 FORMAT(5e16.8)

CALL DIHPAR(nptra,pk,pn,phase,gamc,gams,ipn,fmn)

DO i=1,natom
   IF(chrg(i)/=0) print *,'Warning : non zero charge on particle !', i,chrg(i)
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *
print *,'.................Reading file scale.dat.................'
OPEN(10,STATUS='old',FILE='scale.dat')
DO i=1, 272
  READ(10,*,IOSTAT=j) idum, score(i) !PM
  IF (j .NE. 0) EXIT  !PM
ENDDO
CLOSE(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *
print *,'................Reading file ichain.dat................'
OPEN(10,STATUS='old',FILE='ichain.dat')
READ(10,*) nfrag
DO i=1, nfrag
   READ(10,*) j, lenfrag(i)
ENDDO
print *,'numer of chains/fragments:',nfrag

n = 1
DO i=1, nfrag
   DO j=1, lenfrag(i)
      ichain(n) = i
      n = n +1
   ENDDO
ENDDO
CLOSE(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *
print *,'.................Reading file '//TRIM(input)//'.list.................'
OPEN(10,STATUS='old',FILE=TRIM(input)//'.list')
DO i=1, nres
  READ(10,7205) resp(i),ialpha(i),ibeta(i)
ENDDO

nb = 0
DO
   READ(10,'(a)') string
   IF(string=='-999') EXIT
   nb = nb+1
   IF (nb > nb_max) THEN !PM
     PRINT *,"Error: number of hydropatic pairs too high. Increase nb_max and recompile."
     STOP
   END IF
   READ(string,7200) ni(nb), ivi(nb), nj(nb), ivj(nb), rncoe(nb), vamax(nb), icoeff(nb)
   ! READ(string,*) ni(nb), dums, ivi(nb), nj(nb), dums, ivj(nb), rncoe(nb), vamax(nb), icoeff(nb)
   !IF(icoeff(nb)==-2) &
   !print *, ni(nb), ivi(nb), nj(nb), ivj(nb), rncoe(nb), vamax(nb), icoeff(nb)
ENDDO
READ(10,*) ! noms ! UNUSED
READ(10,*) ! nvalr ! UNUSED
CLOSE(10)

print *,'number of hydropatic pairs nb:',nb

! snippet from protein-2006.f - initialise_protein lines ~273

! scalinf_factor = a_scaling_factor = force_scaling_factor = 1 ! SM
      do i=1,266
        score(i) = dble(scaling_factor)*score(i)
      enddo
      score(270) = dble(scaling_factor)*score(270)
! --- NEW SEPT 06
      score(267) = score(267)*3.0d0/1.6d0
      score(268) = score(268)*4.0d0 ! 3.0 before
! -- end scale MD

! snippet from initialize_protein
walpha(1:20) = 1.0d0
wbeta(1:20) = 1.0d0

CALL prop(foal,fobe)

DO i = 1, nb

   ct0lj(i) = rncoe(i)
   ! ct2lj(i) = vamax(i)*vamax(i) ! BASICALLY unneeded

   IF (icoeff(i).eq.-1) THEN

      ct0lj(i) = ct0lj(i) * score(224) !! ponderation CA-CA

   ELSE IF (icoeff(i) /= -2) THEN

      IF (rncoe(i) .gt. 0.) THEN
         ct0lj(i) = 1.5d0*score(icoeff(i))
      ELSE
         ct0lj(i) = score(icoeff(i))
      ENDIF

      ct0lj(i) = rncoe(i) *  ct0lj(i) !! ponderation Sc-Sc

   ENDIF
 

!  first Set-up the right parameters for each H-bond
   IF (ichain(ni(i)) == ichain(nj(i))) THEN !! INTRA VS. INTER-CHAINS

      IF (ivi(i) /= (ivj(i)-4) .and. ivi(i) /= (ivj(i)+4) ) THEN
         epshb = 2.0d0 * score(223) !! this corresponds to j > i+4
      ELSE
         epshb = 1.25d0 * score(222) !! this corresponds to j = i+4
      ENDIF
      IF (ivi(i) == (ivj(i)-5)) THEN
         epshb = 2.25d0 *  score(223) !! 1.75 23-dec to PI-helix
      ENDIF
!     bug resolved 8 dec 05 Ph.D.
      IF (abs(ivi(i)-ivj(i)) >= 15 .AND. abs(ivi(i)-ivj(i)) <= 30)THEN
         epshb = 0.75d0 *  score(223) !! pour 1FSD
      ENDIF
   ELSE
      epshb = 2.0d0 * score(223) !! INTER-CHAIN INTERACTIONS
   ENDIF                   !! INTRA VS INTER-CHAINS

   epshb_mcmc(i) = epshb

ENDDO
print *,'mixmax(icoeff):',MINVAL(icoeff),MAXVAL(icoeff)
print *,'mixmax(score):',MINVAL(score),MAXVAL(score)
print *,'mixmax(rncoe):',MINVAL(rncoe),MAXVAL(rncoe)
print *,'mixmax(ct0lj):',MINVAL(ct0lj),MAXVAL(ct0lj)

DO i=225,244
   walpha(i-224) = score(i)
ENDDO
DO i=245,264
   wbeta(i-244) = score(i)
ENDDO
walpha_foal(1:20) = walpha(1:20)*foal(1:20)
wbeta_fobe(1:20) = wbeta(1:20)*fobe(1:20)


! end snippet from initialize_protein

! 12-6 forces from 1-4 of dihedral lists
n14 = 0
DO n=1, nphih

   i = iph(n)/3+1; j = jph(n)/3+1; k = ABS(kph(n))/3+1; l = ABS(lph(n))/3+1; kk = icph(n)

   kdiv = (2+ISIGN(1,kph(n))+ISIGN(1,lph(n)))/4 ! kdiv=1 iff i>0 and j>0 else kdiv=0
   if(kdiv==0) THEN
      ! print *,'negative kdiv:',iph(n),jph(n),kph(n),lph(n)
      CYCLE
   ENDIF

   n14 = n14 + 1
   i14(n14,1:3)  = (/kk,i,l/)

ENDDO

DO n=1, nphia

   ! i = ABS(ip(n))/3+1; j = ABS(jp(n))/3+1; k = ABS(kp(n))/3+1; l = ABS(lp(n))/3+1; kk = icp(n)
   i = ip(n)/3+1; j = jp(n)/3+1; k = ABS(kp(n))/3+1; l = ABS(lp(n))/3+1; kk = icp(n)

   kdiv = (2+ISIGN(1,kp(n))+ISIGN(1,lp(n)))/4 ! kdiv=1 iff i>0 and j>0 else kdiv=0
   if(kdiv==0) THEN
      ! print *,'negative kdiv:',iph(n),jph(n),kph(n),lph(n)
      CYCLE
   ENDIF

   n14 = n14 + 1
   i14(n14,1:3)  = (/kk,i,l/)
ENDDO
print *
print *,'n14:',n14


n_ca_sc=0
DO n=1, nb

   IF (icoeff(n) /= -2 .and. rncoe(n)<=9.0d0) THEN

      n_ca_sc = n_ca_sc+1

   ENDIF

ENDDO
print *
print *,'n_ca_sc:',n_ca_sc

n_hb=0

print *,'Hbonds computed on the fly based on bioatoms list'


7200 format(i4,4x,2i4,4x,i4,f8.3,f13.3,i9)
7205 format(1x,a,2i7)

! compute 4-bodies hbonds
npar4b = 0

print *,'HB4B computed on-the-fly based on bioatoms list'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *
print *,'...................Reading file '//TRIM(input)//'.pdb...................'
OPEN(10,STATUS='old',FILE=TRIM(input)//'.pdb')
READ(10,*) 
i=0
ipp = 0
inn = 0
ih=0
DO 
  READ(10,'(a)',ERR=1,END=1) string
  i = i+1
  ! READ(string,'(13x,a4,a3,11x,3f8.3)') pdbatm(i),pdbres(i),x(i),y(i),z(i)
  ! READ(string,'(13x,a4,a3,i6,5x,3f8.3)') pdbatm(i),pdbres(i),ires,x(i),y(i),z(i)
  READ(string,'(13x,a4,a3,i6,4x,3f8.3)') pdbatm(i),pdbres(i),ires,x(i),y(i),z(i)
  ipdbres(i) = ires

  
  IF(lnet) THEN
     IF(lnet_H) THEN
        IF(pdbatm(i)(1:1).NE.'H'.OR.(pdbatm(i)(1:3).EQ.'HIS')) THEN   
           ih=ih+1
           WRITE(80,'(a4,2x,i5,2x,a4,a3,i6,4x,3f8.3)') 'ATOM',ih,pdbatm(i),pdbres(i),ires,x(i),y(i),z(i)
        ENDIF
     ELSE
        IF((pdbatm(i)(1:2).EQ.'CA').OR.(pdbatm(i)(1:3).EQ.pdbres(i)(1:3))) THEN
           ih=ih+1
           WRITE(80,'(a4,2x,i5,2x,a4,a3,i6,4x,3f8.3)') 'ATOM',ih,pdbatm(i),pdbres(i),ires,x(i),y(i),z(i)
        ENDIF
     ENDIF
  ENDIF

  
  IF(lionpair) THEN
     ! patch for ion pairs
     do i_cp=1,2
       IF(pdbatm(i) == crg_plus(i_cp)) THEN 
          ipp=ipp+1
          index_ip_plus(ipp)=i
          itype_ip_plus(ipp)=i_cp
          resindex_ip_plus(ipp)=ires
       ENDIF
     enddo

     do i_cn=1,2
       IF(pdbatm(i) == crg_minus(i_cn)) THEN  
          inn=inn+1
          index_ip_minus(inn)=i
          itype_ip_minus(inn)=i_cn
          resindex_ip_minus(inn)=ires
       ENDIF
     enddo

   ENDIF

ENDDO
1 IF(i/=natom) THEN
   print *,i,natom
   STOP 'error in pdb'
ENDIF
CLOSE(10)
IF(lnet) THEN
   CLOSE(80)
ENDIF

!!!-------------------------------------------------------------
!!!-------- Elastic Network Model ------------------------------
!!!-------------------------------------------------------------
ih=0
heavy=0
IF(lnet) THEN
   net_constant=5.0 !kcal mol to test
   k=0
!   net_rc=6
   pair_k(:)=net_constant
   pair_d0=0

   natom_heavy=0
   DO i=1,natom
      IF(lnet_H) THEN
         IF((pdbatm(i)(1:1).NE.'H').OR.(pdbatm(i)(1:3).EQ.'HIS')) THEN
            natom_heavy=natom_heavy+1     
            ih=ih+1
            heavy(ih)=i
         ENDIF
      ELSE
         IF((pdbatm(i)(1:2).EQ.'CA').OR.(pdbatm(i)(1:3).EQ.pdbres(i)(1:3))) THEN
            natom_heavy=natom_heavy+1     
            ih=ih+1
            heavy(ih)=i
         ENDIF
      ENDIF
   ENDDO

   DO i=1,natom_heavy-1
      ii=heavy(i)
      DO j=i+1,natom_heavy
         jj=heavy(j)
         ddx=x(ii)-x(jj)
         ddy=y(ii)-y(jj)
         ddz=z(ii)-z(jj)
         d0=sqrt(ddx**2+ddy**2+ddz**2)
         IF(d0.LT.net_rc) THEN
            k=k+1
            pair_net(k,1)=i
            pair_net(k,2)=j
            pair_d0(k)=d0
         ENDIF
      ENDDO
   ENDDO
   net_pair=k


ENDIF


!--- the 3 value define the type of potential according to the 2*2 matrix ION-PAIR(i_cp,i_cn)
!--- (Asp/Arg,Asp/Lys)
!--- (Glu/Arg,Glu/Lys)
!--- index_ip(3,npair)=ii+jj*(jj-1) pick the type of potential to be read

npair = 0
IF(lionpair) THEN

DO i_cp = 1,ipp
   DO i_cn = 1,inn

     ! print *,'i_cp,i_cn:',i_cp,i_cn,resindex_ip_plus(i_cp),resindex_ip_minus(i_cn)

      IF(resindex_ip_plus(i_cp) == (resindex_ip_minus(i_cn)+1)) CYCLE
      IF(resindex_ip_plus(i_cp) == (resindex_ip_minus(i_cn)-1)) CYCLE
      IF(resindex_ip_plus(i_cp) == (resindex_ip_minus(i_cn)+2)) CYCLE
      IF(resindex_ip_plus(i_cp) == (resindex_ip_minus(i_cn)-2)) CYCLE
!      IF(resindex_ip_plus(i_cp) == (resindex_ip_minus(i_cn)+3)) CYCLE
!      IF(resindex_ip_plus(i_cp) == (resindex_ip_minus(i_cn)-3)) CYCLE

      npair=npair+1

      ip_info%index_ip(npair,1) = index_ip_plus(i_cp)
      ip_info%index_ip(npair,2) = index_ip_minus(i_cn)

      ii = itype_ip_plus(i_cp)
      jj = itype_ip_minus(i_cn)

      ip_info%index_ip(npair,3) = ii+jj*(jj-1)

      ip_info%index_ca(npair,1) = ip_info%index_ip(npair,1)-1
      ip_info%index_ca(npair,2) = ip_info%index_ip(npair,2)-1
   ENDDO
ENDDO
ip_info%n_ip = npair

print *,'---> # Ion pairs :', ip_info%n_ip
do kk=1,ip_info%n_ip
   print *,kk,'Pair: ',pdbres(ip_info%index_ip(kk,1)),ip_info%index_ip(kk,1), &
                       pdbres(ip_info%index_ip(kk,2)),ip_info%index_ip(kk,2), &
              'Type: ',ip_info%index_ip(kk,3)
enddo
ENDIF


! convert from nm to angstrom...no sono gia' in angstrom
! x=x*10; y=y*10; z=z*10

cm(1:3) = (/SUM(x(1:natom)),SUM(y(1:natom)),SUM(z(1:natom))/)/natom

x(1:natom) = x(1:natom) - cm(1) + nx/2
y(1:natom) = y(1:natom) - cm(2) + ny/2
z(1:natom) = z(1:natom) - cm(3) + nz/2


 !!!!!!!!! Ramachandran angles

 nrama = 0

 ! seek PHI angles (C-N-CA-C)
 DO i=5, natom-4

    IF(pdbatm(i)(1:2)/='CA') CYCLE ! select C_alpha  only !

    IF(pdbres(i)=='PRO') CYCLE ! skip Proline

    IF(pdbres(i)=='DPR') CYCLE ! dont know

    IF(ichain(i-4) == ichain(i+2)) THEN ! all atoms must belong to same chain

        nrama = nrama + 1

        IF(pdbres(i)/='GLY') THEN
           irama(1:4,nrama) = (/i-4,i-2,i,i+2/)
        ELSE
           irama(1:4,nrama) = (/i-4,i-2,i,i+1/)
        ENDIF

        IF(pdbres(i)/='GLY' .AND. pdbres(i)/='ASP') THEN
           prama(1,nrama) = 1.1d0 ! force_k / rad2
        ELSE
           prama(1,nrama) = 0.5d0 ! forceg_k / rad2
        ENDIF

        prama(2,nrama) = -160.d0 ! lower angle
        prama(3,nrama) = -60.d0  ! upper angle

    ENDIF
       
 ENDDO

 ! seek PSI angles (N-CA-C-N)
 DO i=5, natom-2

    IF(pdbatm(i)(1:2)/='CA') CYCLE ! select C_alpha  only !

    IF(pdbres(i)=='PRO') CYCLE ! skip Proline

    IF(pdbres(i)=='DPR') CYCLE ! dont know

    IF(ichain(i-2) == ichain(i+4)) THEN ! all atoms must belong to same chain

        nrama = nrama + 1

        IF(pdbres(i)/='GLY') THEN
           irama(1:4,nrama) = (/i-2,i,i+2,i+4/)
        ELSE
           irama(1:4,nrama) = (/i-2,i,i+1,i+3/)
        ENDIF

        IF(pdbres(i)/='GLY' .AND. pdbres(i)/='ASP') THEN
           prama(1,nrama) = 1.1d0 ! force_k / rad2
        ELSE
           prama(1,nrama) = 0.5d0 ! forceg_k / rad2
        ENDIF

        prama(2,nrama) = -60.d0   ! lower angle
        prama(3,nrama) = +160.d0  ! upper angle

    ENDIF
       
 ENDDO
 ! prama(1,nrama) = prama(1,nrama) * (180.d0/PI)**2
 ! prama(1,nrama) = prama(1,nrama) * (180.d0/PI)**2
 ! prama(1,:) = prama(1,:) * (180.d0/PI)**2
 prama(1,:) = prama(1,:) * scaling_factor * 2

 ! prama(1,:) = prama(1,:) * 2.6 ! pompato da fabio
 ! prama(1,:) = prama(1,:) * 1000. ! pompato da fabio

!print *,'ATTTTTTTTNNN nbonh,nbona=0:';nbonh=0; nbona=0
!print *,'ATTTTTTTTNNN ntheth,ntheta=0:';ntheth=0; ntheta=0
!print *,'ATTTTTTTTNNN nphih=0,phia=0:';nphih=0;nphia=0
!print *,'ATTTTTTTTNNN n14=0:'; n14=0          ! dihedral 1-4 forces
!print *,'ATTTTTTTTNNN nrama=0:'; nrama=0          ! ramachandran angles
!print *,'ATTTTTTTTNNN n_ca_sc=0:'; n_ca_sc=0 ! Ca-Ca + Ca-Sc + Sc-Sc Hydropathic forces
!print *,'ATTTTTTTTNNN n_hb=0:'; n_hb=0 ! 2-body hydrogen bonds
!print *,'ATTTTTNNNNNN nb_par4b=0:'; npar4b = 0 ! 4-body hydrogen bonds
!print *,'ATTTTTNNNNNN VDW~score(270)=0:'; score(270) = 0 ! VDW forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!            WRITING atom.inp            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(10,FILE='atom.inp')

WRITE(10,'(a)') 'CONTROL'
WRITE(10,'(2i4,g12.5,a)') 1, 1, timestep,'  ! MD sub-sub-iterations,sub-iterations,smallest timestep'
! WRITE(10,'(2i4,g12.5,a)') 10, 1, timestep/10.0,'  ! MD sub-sub-iterations,sub-iterations,smallest timestep : modified for TITAN'

WRITE(10,'(1000(a,/))') &
'.false.     ! rotational motion',&
'2          ! 1:NGP 2:IB',&
'0. 0. 0.   ! force',&
'.false.     ! linsert',&
'CLOSECONTROL',&
'',&
'MODEL',&
'# title',&
'molecular types    2',&
'_AnyName',&
'nummols    1'

DO i=1, natom
   WRITE(labatm(i),'(a,i0)') igraph(i),iac(i)
   CALL compress_blanks(labatm(i))
ENDDO
print *
print *,'Atom labels : ',labatm(1:4),' ... ... ...   ',labatm(natom-3:natom)

IF(.NOT.lnet) THEN
   WRITE(10,'(a,i8)') 'atoms', natom
ELSE
   WRITE(10,'(a,i8)') 'atoms', natom_heavy
ENDIF
k=0
DO i=1, natom


   xmass=amass(i)
   IF((i.gt.1).AND.(i.lt.natom-2)) THEN
      IF(lmass) THEN
         IF(labatm(i)(1:2).EQ.'C1') xmass=25607.
         IF(labatm(i)(1:1).EQ.'O') xmass=25607.
         IF(labatm(i)(1:1).EQ.'N') xmass=25607.
         IF(labatm(i)(1:2).EQ.'HN') xmass=25607.
      ENDIF
   ENDIF
!   IF(labatm(i)(1:1)=='H') THEN
!      xmass=10*amass(i)
!   ELSE
      
!   ENDIF

   idonor=0
   ihydrogen=0
   iacceptor=0

   IF(.not.lnet) THEN
      IF(i<natom) THEN
         IF(labatm(i)(1:1)=='N' .AND. labatm(i+1)(1:1)=='H') idonor = 1
         IF(labatm(i)(1:1)=='H' .AND. labatm(i-1)(1:1)=='N') ihydrogen = 1
      ENDIF
      IF(labatm(i)(1:1)=='O') iacceptor = 1
   ENDIF

   ! WRITE(10,'(a,2g12.5, a, 3(1x,i0), 2x, 2(1x,i0), a, i0,1x,a,1x,a,2(1x,i0))') &

   IF(.NOT.lnet) THEN
      WRITE(10,'(a,2g12.5, a, a, i0,1x,a,1x,a,2(1x,i0))') &
           labatm(i),xmass,chrg(i), &
           ' 1  0  10. 10. 10. 1 ', &
           ! idonor,ihydrogen,iacceptor, &
           ! ipdbres(i),ichain(i), &
           ' ! name, mass, chge, rept, frzt, inertia x/y/z, frzr, ', &
           ! ' ! name, mass, chge, rept, frzt, inertia x/y/z, frzr, hbond(D,H,A), res,chain ', &
           i,igraph(i),pdbres(i)
   ELSE
      IF(lnet_H) THEN
         IF((labatm(i)(1:1).NE.'H').OR.(labatm(i)(1:3).EQ.'HIS')) THEN
            k=k+1
            WRITE(10,'(a,2g12.5, a, a, i0,1x,a,1x,a,2(1x,i0))') &
                 labatm(i),xmass,chrg(i), &
           ' 1  0  10. 10. 10. 1 ', &
           ! idonor,ihydrogen,iacceptor, &
           ! ipdbres(i),ichain(i), &
           ' ! name, mass, chge, rept, frzt, inertia x/y/z, frzr, ', &
           ! ' ! name, mass, chge, rept, frzt, inertia x/y/z, frzr, hbond(D,H,A), res,chain ', &
           k,igraph(i),pdbres(i)
         ENDIF
      ELSE
         IF((pdbatm(i)(1:2).EQ.'CA').OR.(pdbatm(i)(1:3).EQ.pdbres(i)(1:3))) THEN
            k=k+1
            WRITE(10,'(a,2g12.5, a, a, i0,1x,a,1x,a,2(1x,i0))') &
                 labatm(i),xmass,chrg(i), &
                 ' 1  0  10. 10. 10. 1 ', &
                 ! idonor,ihydrogen,iacceptor, &
                 ! ipdbres(i),ichain(i), &
                 ' ! name, mass, chge, rept, frzt, inertia x/y/z, frzr, ', &
                 ! ' ! name, mass, chge, rept, frzt, inertia x/y/z, frzr, hbond(D,H,A), res,chain ', &
                 k,igraph(i),pdbres(i)
         ENDIF
      ENDIF
   ENDIF

ENDDO

IF(.NOT.lnet) THEN
   WRITE(10,'(a,i8)') 'bioatoms', natom
ELSE
   WRITE(10,'(a,i8)') 'bioatoms', natom_heavy
ENDIF

k=0
DO i=1, natom

!   IF(labatm(i)(1:1)=='H') THEN
!      xmass=10*amass(i)
!   ELSE
      xmass=amass(i)
!   ENDIF

   idonor=0
   ihydrogen=0
   iacceptor=0
   IF(.not.lnet) THEN
      IF(i<natom) THEN
         IF(labatm(i)(1:1)=='N' .AND. labatm(i+1)(1:1)=='H') idonor = 1
         IF(labatm(i)(1:1)=='H' .AND. labatm(i-1)(1:1)=='N') ihydrogen = 1
      ENDIF
      IF(labatm(i)(1:1)=='O'.AND.ipdbres(i).LT.nres) iacceptor = 1
   ENDIF


   IF(.NOT.lnet) THEN
      WRITE(10,'(a, 3(1x,i0), 2x, 2(1x,i0), a, i0,1x,a,1x,a,2(1x,i0))') &
           labatm(i), &
           idonor,ihydrogen,iacceptor, &
           ipdbres(i),ichain(i), &
           ' ! hbond(D,H,A), residue, chain ', &
           i,igraph(i),pdbres(i)
   ELSE
      IF(lnet_H) THEN
         IF(labatm(i)(1:1).NE.'H'.OR.(labatm(i)(1:3).EQ.'HIS')) THEN
            k=k+1
            WRITE(10,'(a, 3(1x,i0), 2x, 2(1x,i0), a, i0,1x,a,1x,a,2(1x,i0))') &
                 labatm(i), &
                 idonor,ihydrogen,iacceptor, &
                 ipdbres(i),ichain(i), &
                 ' ! hbond(D,H,A), residue, chain ', &
                 k,igraph(i),pdbres(i)
         ENDIF
      ELSE
         IF((pdbatm(i)(1:2).EQ.'CA').OR.(pdbatm(i)(1:3).EQ.pdbres(i)(1:3))) THEN
            k=k+1
            WRITE(10,'(a, 3(1x,i0), 2x, 2(1x,i0), a, i0,1x,a,1x,a,2(1x,i0))') &
                 labatm(i), &
                 idonor,ihydrogen,iacceptor, &
                 ipdbres(i),ichain(i), &
                 ' ! hbond(D,H,A), residue, chain ', &
                 k,igraph(i),pdbres(i)
         ENDIF
      ENDIF
   ENDIF  
ENDDO

!!!!!!!!!!!!! bonds
IF(.NOT.lnet) THEN
   WRITE(10,'(a,i8,2f10.5)') 'bonds', nbonh + nbona + n14 + n_ca_sc + ip_info%n_ip, 9., 10.
   
   DO n=1, ip_info%n_ip
      WRITE(10,'(a,2i6,2g16.8,a)') 'harm',ip_info%index_ip(n,1), ip_info%index_ip(n,2), 0., 0.,' ! ion pairs excluded'
   ENDDO
   
   DO n=1, nbonh
      i = ibh(n)/3+1; j = jbh(n)/3+1; kk = icbh(n)
      WRITE(10,'(a,2i6,2g16.8)') 'harm',i,j, 2.0*score(267)*rk(kk), req(kk) ! fix 30 Jan 2013
      ! WRITE(10,'(a,2i6,2g16.8)') 'harm',i,j, score(267)*rk(kk), req(kk)
      ! WRITE(10,'(a,2i6,2g16.8)') 'harm',i,j, 0.0*score(267)*rk(kk), req(kk) ! ATTTN
   ENDDO
   DO n=1, nbona
      i = ib(n)/3+1; j = jb(n)/3+1; kk = icb(n)
      WRITE(10,'(a,2i6,2g16.8)') 'harm',i,j, 2.0*score(267)*rk(kk),req(kk) ! fix 30 Jan 2013
      ! WRITE(10,'(a,2i6,2g16.8)') 'harm',i,j, score(267)*rk(kk),req(kk)
      ! WRITE(10,'(a,2i6,2g16.8)') 'harm',i,j, 0.0*score(267)*rk(kk),req(kk) ! ATTTN
   ENDDO
   
   DO n=1, n14
      
      ! NB: atom ids can be > or < 0
      ic0 = i14(n,1); i = i14(n,2); j = i14(n,3)
      
      ! kdiv = (2+ISIGN(1,i)+ISIGN(1,j))/4
      
      i = ABS(i); j = ABS(j)
      
      ia1 = iac(i); ia2 = iac(j)
      
      ibig = MAX(ia1,ia2)
      isml = MIN(ia1,ia2)
      ic = ibig*(ibig-1)/2+isml
      
      ! print *,'n,ic0,fmn:',n,ic0,fmn(ic0)
      
      ! taken from ETORS, lines 727
      ! NB fmn(ic0)=0 if dihedral periodicity pn < 0
      
      c12 = scal14 * cn1(ic) * fmn(ic0)
      c6 = scal14 * cn2(ic) * fmn(ic0)
      
      WRITE(10,'(a,2i6,2g16.8)') '12-6',i,j, c12, c6
      ! WRITE(10,'(a,2i6,2g16.8)') '12-6',i,j, 0., 0. !ATTT
      
   ENDDO
   
   ! add hydropatic stuff (Ca-Ca, Ca-Sc, Sc-Sc and possibly others)
   IF(n_ca_sc>0) THEN
      
      
      DO n=1, nb
         
         
         IF (icoeff(n) /= -2 .AND. rncoe(n)<=9.0d0) THEN
            
            IF(rncoe(n)>0) THEN
               WRITE(10,'(a,2i6,2g16.8,2i8)') 'opep', ni(n), nj(n), vamax(n), ct0lj(n)
               ! WRITE(10,'(a,2i6,2g16.8)') 'opep', ni(n), nj(n), 0.0*vamax(n), ct0lj(n) ! ATTTN
            ELSE
               WRITE(10,'(a,2i6,2g16.8)') 'opem', ni(n), nj(n), vamax(n), ABS(ct0lj(n))
            ENDIF
            
         ENDIF
         
      ENDDO
      
   ENDIF
!!!---- Elastic Pairs -----
ELSE
   WRITE(10,'(a,i8,2f10.5)') 'bonds', net_pair, 9., 10.
   DO n=1,net_pair
      inet=pair_net(n,1)
      jnet=pair_net(n,2)
      WRITE(10,'(a,2i6,2g16.8)') 'harm',inet,jnet, pair_k(n),pair_d0(n) ! fix 30 Jan 2013
   ENDDO
ENDIF
      

!!!!!!!!!!!!! angles

IF(.NOT.lnet) THEN
   WRITE(10,'(a,i8,2f10.5)') 'angles', ntheth + ntheta + n_hb, 9.0, 10.

   DO n=1, ntheth
      i = ith(n)/3+1; j = jth(n)/3+1; k = kth(n)/3+1; kk = icth(n)
      WRITE(10,'(a,3i6,2f9.3)') 'harm',i,j,k, 2.0*score(268)*tk(kk), teq(kk)*180.0/PI ! fix 30 Jan 3013
      ! WRITE(10,'(a,3i6,2f9.3)') 'harm',i,j,k, score(268)*tk(kk), teq(kk)*180.0/PI
   ENDDO
   DO n=1, ntheta
      i = it(n)/3+1; j = jt(n)/3+1; k = kt(n)/3+1; kk = ict(n)
      WRITE(10,'(a,3i6,2f9.3)') 'harm',i,j,k, 2.0*score(268)*tk(kk), teq(kk)*180.0/PI ! fix 30 Jan 2013
      ! WRITE(10,'(a,3i6,2f9.3)') 'harm',i,j,k, score(268)*tk(kk), teq(kk)*180.0/PI
   ENDDO
   
   IF(n_hb>0) THEN
      
      DO n=1, nb
         
         IF (icoeff(n) /= -2 .and. rncoe(n)<=9.0d0) CYCLE ! complementary to pair hydropatic terms
         
         IF (icoeff(n) /= -1 .and. rncoe(n)>=25.0d0 .and. rncoe(n)<=29.0d0) THEN
            
            WRITE(10,'(a,3i6,2f9.3)') 'hbop', nj(n), nj(n)+1, ni(n), epshb_mcmc(n), 1.8d0
            ! WRITE(10,'(a,3i6,2f9.3)') 'hbop', nj(n), nj(n)+1, ni(n), epshb_mcmc(n), scaling_factor ! ERROR
            
         ENDIF
         
      ENDDO
      
   ENDIF
ENDIF

!!!!!!!!!!!!!! dihedrals
IF(.NOT.lnet) THEN
   WRITE(10,'(a,i8,2f10.5)') 'dihedrals', nphih + nphia + npar4b + nrama + ip_info%n_ip, 9., 10.
   
   DO n=1, nphih
      i = iph(n)/3+1; j = jph(n)/3+1; k = ABS(kph(n))/3+1; l = ABS(lph(n))/3+1; kk = icph(n)
      WRITE(10,'(a4,4i6,4f12.6)') 'cos',i,j,k,l, &
           scaling_factor * 0.232*pk(kk), +phase(kk)*180.0/PI, pn(kk), 0. ! fix 6 Feb
   ENDDO
   DO n=1, nphia
      i = ip(n)/3+1; j = jp(n)/3+1; k = ABS(kp(n))/3+1; l = ABS(lp(n))/3+1; kk = icp(n)
      
      IF(pk(kk)>= 3.05 .AND. pk(kk) <= 3.15) THEN
         ! without LJ exclusion
         WRITE(10,'(a4,4i6,4f12.6)') '-cos',i,j,k,l, &
              scaling_factor * 0.232*pk(kk), +phase(kk)*180.0/PI, pn(kk), 0. ! fix 6 Feb
      ELSE
         ! with LJ exclusion
         WRITE(10,'(a4,4i6,4f12.6)') 'cos',i,j,k,l, &
              scaling_factor * 0.232*pk(kk), +phase(kk)*180.0/PI, pn(kk), 0. ! fix 6 Feb
      ENDIF
   ENDDO
ENDIF
!!!!!    ramachandran restraints

IF(.NOT.lnet) THEN
   DO n=1, nrama
      i = irama(1,n); j = irama(2,n); k = irama(3,n); l = irama(4,n)
      WRITE(10,'(a4,4i6,4f16.6)') 'rama',i,j,k,l, prama(1:3,n), 0.
   ENDDO
   
!!!!    4-body hydrogen bonds
   DO n=1, npar4b
      i = ipar4b(1,n); j = ipar4b(2,n); k = ipar4b(3,n); l = ipar4b(4,n)
      WRITE(10,'(a4,4i6,4f12.6)') 'hb4b',i,j,k,l, par4b(1:2,n), 0., 0.
   ENDDO
   
!!!!    ion pairs
   DO kk=1,ip_info%n_ip
      
      i = ip_info%index_ca(kk,1)
      j = ip_info%index_ip(kk,1)
      k = ip_info%index_ip(kk,2)
      l = ip_info%index_ca(kk,2)
      
      IF    (ip_info%index_ip(kk,3) == 1) THEN
         
         WRITE(10,'(a4,4i6,4f12.6)') 'ip4b',i,j,k,l, 1.,0.,0.,0.
         
      ELSEIF(ip_info%index_ip(kk,3) == 2) THEN
         
         WRITE(10,'(a4,4i6,4f12.6)') 'ip2b',i,j,k,l, 2.,0.,0.,0.
         
      ELSEIF(ip_info%index_ip(kk,3) == 3) THEN
         
         WRITE(10,'(a4,4i6,4f12.6)') 'ip4b',i,j,k,l, 3.,0.,0.,0.
         
      ELSEIF(ip_info%index_ip(kk,3) == 4) THEN
         
         WRITE(10,'(a4,4i6,4f12.6)') 'ip2b',i,j,k,l, 4.,0.,0.,0.
         
      ELSE
         STOP 'unrecognized ion pair potential'
      ENDIF
      
   ENDDO
ENDIF
!!!!!!!!!!!!! other stuff
WRITE(10,'(1000(a,/))') &
'hydro',&
'.false. 1.e-6  5.0  2.0   2.0  0.0  ! passive scalar,tumbling coeff,visc_enhancer,csi,oblate,smooth',&
'0.1   0.0   ! gammas',&
! '1.d-6   0.0   ! gammas',&
'finish',&
'_Wall',&
'nummols    1',&
'atoms    1',&
'W',&
'finish'


nuniq=0
tuniq=0
loop_i: DO i=1, natom


!!!--- INTERMOLECULAR TERMS NEED TO INCLUDE ALSO Sc-Sc interactions similar to INTRAMOLECULAR ONES 

   IF(is_in_alist(labatm(i),auniq,nuniq)) CYCLE
   

   IF(.not.lnet) THEN
      nuniq=nuniq+1
      WRITE(auniq(nuniq),'(a)') labatm(i)
      iuniq(nuniq) = i
   ELSE
      IF(lnet_h) THEN
         IF((labatm(i)(1:1).EQ.'H').AND.(labatm(i)(1:3).NE.'HIS')) CYCLE
         nuniq=nuniq+1
         WRITE(auniq(nuniq),'(a)') labatm(i)
         iuniq(nuniq) = i
      ELSE
         IF((labatm(i)(1:2).NE.'CA').AND.(labatm(i)(1:3).NE.pdbres(i)(1:3))) CYCLE
            nuniq=nuniq+1
            WRITE(auniq(nuniq),'(a)') labatm(i)
            iuniq(nuniq) = i
      ENDIF
   ENDIF

   DO j=1,SIZE(sc_name)                       
      IF(auniq(nuniq)(1:3)==sc_name(j)) tuniq(nuniq)=1
   ENDDO
   

ENDDO loop_i

print *
print *,'Unique list : ',nuniq,' : ',auniq(1:nuniq)
print *

!DO i=1, natom
!   DO j=1, natom
!       ia1=iac(i); ia2=iac(j)
!       ibig = MAX(ia1,ia2); isml = MIN(ia1,ia2); ic = ibig*(ibig-1)/2+isml
!       WRITE(1000,*) igraph(i),igraph(j),cn1(ic)
!   ENDDO
!ENDDO

!!!!!!!!!!!!! write vdw 

!!!!!!!!!!!!! Select Sc-Sc interactions in the systems
!--- Append type to labels first
ipair=.FALSE.
DO i=1,n_scsc

   DO j=1,natom
      IF(TRIM(pair(i,1))==labatm(j)(1:3)) THEN
         pair(i,1)=labatm(j)         
         ipair(i,1)=.TRUE.
      ENDIF
      IF(TRIM(pair(i,2))==labatm(j)(1:3)) THEN
         ipair(i,2)=.TRUE.
         pair(i,2)=labatm(j)
      ENDIF
   ENDDO
ENDDO

n_scsc_sys=0
do i=1,n_scsc
   IF(ipair(i,1).AND.ipair(i,2)) THEN
   n_scsc_sys=n_scsc_sys+1
   pair_sys(n_scsc_sys)=i
   ENDIF
ENDDO

!---- Select no Sc-Sc pair from unique
nvdw=0
IF(.not.lnet) THEN
   do n=1,nuniq
      do m=n,nuniq
         IF(tuniq(n)*tuniq(m)==1) CYCLE
         nvdw=nvdw+1
      enddo
   enddo
ELSE
!--- Elastict Net Exclude H ----
   do n=1,nuniq
      auniq(n)=ADJUSTL(auniq(n))
      !      print *,auniq(n)(1:1)
      IF(lnet_h) THEN
         IF((auniq(n)(1:1).EQ.'H').AND.(auniq(n)(1:3).NE.'HIS')) THEN
            CYCLE
         ENDIF
      ELSE
         IF((auniq(n)(1:2).NE.'CA').AND.(auniq(n)(1:3).NE.pdbres(iuniq(n))(1:3))) THEN
            CYCLE
         ENDIF
      ENDIF

      do m=n,nuniq
         IF(tuniq(n)*tuniq(m)==1) CYCLE
         auniq(m)=ADJUSTL(auniq(m))
         IF(lnet_h) THEN
            IF((auniq(m)(1:1).EQ.'H').AND.(auniq(m)(1:3).NE.'HIS')) CYCLE
            nvdw=nvdw+1
         ELSE
            IF((auniq(m)(1:2).NE.'CA').AND.(auniq(m)(1:3).NE.pdbres(iuniq(m))(1:3))) CYCLE
            nvdw=nvdw+1
         ENDIF
      enddo
   enddo
ENDIF

!----

!--WRITE(10,'(a,i6)') 'vdw    ',((nuniq+1)*(nuniq+2))/2+(n_scsc_sys)
WRITE(10,'(a,i6)') 'vdw    ',nvdw+(n_scsc_sys)+nuniq+1
DO n=1, nuniq
   DO m=n, nuniq
      i = iuniq(n)
      j = iuniq(m)
      ia1=iac(i)
      ia2=iac(j)
      ibig = MAX(ia1,ia2)
      isml = MIN(ia1,ia2)
      ic = ibig*(ibig-1)/2+isml

      ! set the cutoff to 10 !!!!!!!
      IF(tuniq(n)*tuniq(m)==1) CYCLE

      WRITE(10,'(3a,5g13.6)') auniq(n), auniq(m), ' 12-6 ', score(270)*cn1(ic), score(270)*cn2(ic), 0., 0., 10.
      ! WRITE(10,'(3a,5g13.6)') auniq(n), auniq(m), ' 12-6 ', score(270)*cn1(ic), score(270)*cn2(ic), 0., 0., 20.
   ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Side-Chain Side-Chain Interactions

!--- Write table
DO j=1,n_scsc_sys
   i=pair_sys(j)
   write(10,'(3a,5g13.6)') pair(i,1),pair(i,2), pair_pot_type(i), pot_data(i,1),pot_data(i,2)*scaling_factor, 0.,0.,10.
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO n=1, nuniq
   IF(.not.lwall) THEN
      WRITE(10,'(2a)') auniq(n), '  W    lj   0.   0.   0.   0.   10.'
   ELSE
      WRITE(10,'(2a,5g13.6)') auniq(n), '  W    lj ',e_wall,s_wall,0.,0.,10.
   ENDIF
ENDDO

WRITE(10,'(a)') 'W       W    lj   0.   0.   0.   0.   10. '

WRITE(10,'(1000(a,/))') 'CLOSEMODEL'


CALL RANDOM_SEED()

CALL RANDOM_NUMBER(rx(1:natom))
CALL RANDOM_NUMBER(ry(1:natom))
CALL RANDOM_NUMBER(rz(1:natom))

rx=0. ! (rx-0.5)*0.001
ry=0. ! (ry-0.5)*0.001
rz=0. ! (rz-0.5)*0.001

WRITE(10,'(a,i8)') 'CONFIG'
IF(.not.lnet) THEN
   WRITE(10,*) natom
   DO i=1, natom
      ! WRITE(10,*) x(i),y(i),z(i),0.,0.,0.
      WRITE(10,'(2(3f12.5,3x))') x(i) + rx(i), y(i) + ry(i), z(i) + rz(i), 0.,0.,0.
   ENDDO
ELSE
   WRITE(10,*) natom_heavy
   DO i=1, natom
      ! WRITE(10,*) x(i),y(i),z(i),0.,0.,0.
      IF(lnet_h) THEN
         IF((labatm(i)(1:1).EQ.'H').AND.(labatm(i)(1:3).NE.'HIS')) CYCLE
         WRITE(10,'(2(3f12.5,3x))') x(i) + rx(i), y(i) + ry(i), z(i) + rz(i), 0.,0.,0.
      ELSE
         IF((labatm(i)(1:2).EQ.'CA').OR.(labatm(i)(1:3).EQ.pdbres(i)(1:3))) THEN
            WRITE(10,'(2(3f12.5,3x))') x(i) + rx(i), y(i) + ry(i), z(i) + rz(i), 0.,0.,0.
         ENDIF
      ENDIF
   ENDDO
ENDIF
WRITE(10,'(a)') 'CLOSECONFIG'

CALL ip_potential_dat(10,scaling_factor)

CLOSE(10)

CONTAINS

!!!!!!!!!
SUBROUTINE DIHPAR(numphi,pk,pn,phase,gamc,gams,ipn,fmn)

implicit none
!     
!     GET ADDITIONAL PARAMETERS FOR THE VECTOR EPHI 
!     
INTEGER :: numphi
REAL :: pk(:),pn(:),phase(:),gamc(:),gams(:)
INTEGER :: ipn(:)
REAL :: fmn(:)
REAL :: DUM,DUMC,DUMS

REAL,PARAMETER :: ZERO=0.0,ONE=1.0,TENM3=1.e-3,TENM6=1.e-6
REAL,PARAMETER :: FOUR=4.0,PI=3.141592653589793

!     
!     PIM = FOUR*ATAN(ONE)
      DO 100 I = 1,numphi
        DUM = phase(I)
        IF(ABS(DUM-PI).LE.TENM3) DUM = SIGN(PI,DUM)
        DUMC = COS(DUM)
        DUMS = SIN(DUM)
        IF(ABS(DUMC).LE.TENM6) DUMC = ZERO
        IF(ABS(DUMS).LE.TENM6) DUMS = ZERO
        gamc(I) = DUMC*pk(I)
        gams(I) = DUMS*pk(I)
!     print*,' I ', I, gamc(I),gams(I)
        fmn(I) = ONE
        IF(pn(I).LE.ZERO) THEN
           fmn(I) = ZERO
        ENDIF
        pn(I) = ABS(pn(I))
        ipn(I) = INT(pn(I)+TENM3)
 100  CONTINUE

END SUBROUTINE DIHPAR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION is_in_alist(labj,alist,nlist) RESULT(out)
CHARACTER(len=*) :: labj, alist(:)
INTEGER :: nlist
LOGICAL :: out
INTEGER :: i

out = .true.
DO i=1, nlist
   IF(labj==alist(i)) RETURN
ENDDO
out = .false.

END FUNCTION 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE compress_blanks(text)

  ! compress all single blanks

    IMPLICIT NONE

    CHARACTER (len=*) :: text
    INTEGER :: i

    DO
      i = INDEX(TRIM(text), " ")
      IF(i==0) EXIT
      text(i:) = text(i+1:)
    ENDDO

END SUBROUTINE compress_blanks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE prop(foal,fobe)

! ------- energies of residues for alpha and beta
! ------- positive values are required for bad propensities

      implicit double precision (a-h,o-z)

      REAL,DIMENSION(20) :: fobe(:),foal(:)

      i=1  !! CYS
      foal(i) = 0.6
      fobe(i) = 0.2

      i=2  !! LEU 
      foal(i) = 0.19
      fobe(i) = 0.2
!      fobe(i) = 0.0

      i=3  !! VAL
      foal(i) = 0.46
      fobe(i) = -0.3

      i=4  !! ILE 
      foal(i) = 0.35
      fobe(i) = -0.3

      i=5  !! MET 
      foal(i) = 0.21
      fobe(i) = 0.3

      i=6  !! PHE 
      foal(i) = 0.47
      fobe(i) = -0.3

      i=7  !! TYR 
      foal(i) = 0.47
      fobe(i) = -0.3

      i=8  !! LYS 
      foal(i) = 0.15
      fobe(i) = 0.3

      i=9  !! ARG 
      foal(i) = 0.06
      fobe(i) = 0.3

      i=10  !! PRO 
      foal(i) = 1.3
      fobe(i) = 0.3

      i=11  !! GLY 
      foal(i) = 1.10
      fobe(i) = 0.3
!      fobe(i) = 0.0

      i=12  !! ALA 
      foal(i) = 0.1
      fobe(i) = 0.3
!      fobe(i) = 0.0 !! 

      i=13  !! GLN 
      foal(i) = 0.32
      fobe(i) = 0.3

      i=14  !! HIS 
      foal(i) = 0.62
      fobe(i) = 0.3

      i=15  !! ASN 
      foal(i) = 0.60
      fobe(i) = 0.3
!      fobe(i) = 0.0  !!

      i=16  !! ASP 
      foal(i) = 0.59
      fobe(i) = 0.3

      i=17  !! GLU 
      foal(i) = 0.34
      fobe(i) = 0.3

      i=18  !! SER 
      foal(i) = 0.52
      fobe(i) = 0.3

      i=19  !! THR 
      foal(i) = 0.57
      fobe(i) = -0.3

      i=20  !! TRP 
      foal(i) = 0.47
      fobe(i) = 0.3
END SUBROUTINE prop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ip_potential_dat(iout,scaling_factor)

INTEGER,INTENT(in) :: iout
REAL,INTENT(in) :: scaling_factor
REAL,PARAMETER :: conversion = 1.27  ! corresponds to f=1 in Sterpone,JCTC,2013
!REAL,PARAMETER :: conversion = 1.875000 ! f=1.3

REAL :: dum,v1,f1,v2,f2

REAL,PARAMETER :: spacing = 0.10, xmin = 3.00, cutoff = 10.00, xmax = 31.00

! '#-- Asp-Arg   *     ', &
CHARACTER(len=80),PARAMETER :: cc1(71) = [&
' 3.00         2.08381        -6.66965     0.00000         0.00000', &
' 3.10         1.41719        -6.66965     0.00000         0.00000', &
' 3.20         0.74988        -6.69640     0.00000         0.00000', &
' 3.30         0.07791        -6.76935     0.00000         0.00000', &
' 3.40        -0.60399        -6.85480     0.00000         0.00000', &
' 3.50        -1.29305        -6.87350     0.00000         0.00000', &
' 3.60        -1.97869        -6.54125     0.00000         0.00000', &
' 3.70        -2.60130        -5.52000     0.00000         0.00000', &
' 3.80        -3.08269        -4.27885     0.00000         0.00000', &
' 3.90        -3.45707        -2.70590     0.00000         0.00000', &
' 4.00        -3.62387        -0.56315     0.00000         0.00000', &
' 4.10        -3.56970         1.13430     0.00000         0.00000', &
' 4.20        -3.39701         2.29815     0.00000         0.00000', &
' 4.30        -3.11007         3.64240     0.00000         0.00000', &
' 4.40        -2.66853         4.19785     0.00000         0.00000', &
' 4.50        -2.27050         3.04855     0.00000         0.00000', &
' 4.60        -2.05882         2.10882     0.00000         0.00000', &
' 4.70        -1.84874         2.25975     6.00000       -14.69135', &
' 4.80        -1.60687         1.83824     2.04315       -14.69135', &
' 4.90        -1.48109         1.36420     0.06173       -14.16400', &
' 5.00        -1.33403         1.31302    -0.78965        -7.58805', &
' 5.10        -1.21849         0.99790    -1.45588        -6.20445', &
' 5.20        -1.13445         0.71972    -2.03054        -3.60375', &
' 5.30        -1.07454         0.68277    -2.17663         1.29880', &
' 5.40        -0.99790         0.75087    -1.77078         4.27050', &
' 5.50        -0.92437         0.57773    -1.32253         3.59020', &
' 5.60        -0.88235         0.47269    -1.05274         2.17645', &
' 5.70        -0.82983         0.54419    -0.88724         1.19359', &
' 5.80        -0.77351         0.52521    -0.81402         0.75131', &
' 5.90        -0.72479         0.40118    -0.73698         0.58011', &
' 6.00        -0.69328         0.42017    -0.69800         0.44044', &
' 6.10        -0.64076         0.42017    -0.64889         0.52052', &
' 6.20        -0.60924         0.36765    -0.59390         0.60515', &
' 6.30        -0.56723         0.36765    -0.52786         0.49528', &
' 6.40        -0.53571         0.31513    -0.49484         0.38522', &
' 6.50        -0.50420         0.31513    -0.45082         0.38522', &
' 6.60        -0.47269         0.26261    -0.41780         0.38522', &
' 6.70        -0.45168         0.26261    -0.37377         0.33019', &
' 6.80        -0.42017         0.36765    -0.35176         0.22012', &
' 6.90        -0.37815         0.31513    -0.32975         0.33019', &
' 7.00        -0.35714         0.15756    -0.28572         0.38522', &
' 7.10        -0.34664         0.15756    -0.25270         0.16509', &
' 7.20        -0.32563         0.21008    -0.25270         0.11006', &
' 7.30        -0.30462         0.21008    -0.23069         0.27516', &
' 7.40        -0.28361         0.21008    -0.19767         0.22012', &
' 7.50        -0.26261         0.26261    -0.18667         0.05503', &
' 7.60        -0.23109         0.21008    -0.18667         0.11006', &
' 7.70        -0.22059         0.10504    -0.16465         0.16509', &
' 7.80        -0.21008         0.15756    -0.15365         0.05503', &
' 7.90        -0.18908         0.15756    -0.15365         0.05503', &
' 8.00        -0.17857         0.10504    -0.14264         0.05503', &
' 8.10        -0.16807         0.21008    -0.14264         0.11006', &
' 8.20        -0.13655         0.21008    -0.12063         0.11006', &
' 8.30        -0.12605         0.10504    -0.12063        -0.05503', &
' 8.40        -0.11555         0.05252    -0.13163         0.11006', &
' 8.50        -0.11555         0.05252    -0.09862         0.16509', &
' 8.60        -0.10504         0.05252    -0.09862         0.05503', &
' 8.70        -0.10504         0.10504    -0.08761         0.16509', &
' 8.80        -0.08403         0.21008    -0.06560         0.12760', &
' 8.90        -0.06303         0.10504    -0.06209        -0.05503', &
' 9.00        -0.06303         0.00000    -0.07660        -0.07257', &
' 9.10        -0.06303         0.10504    -0.07660         0.05503', &
' 9.20        -0.04202         0.15756    -0.06560         0.05503', &
' 9.30        -0.03151         0.10504    -0.06560         0.11006', &
' 9.40        -0.02101         0.10504    -0.04358         0.11006', &
' 9.50        -0.01050         0.05252    -0.04358         0.00000', &
' 9.60        -0.01050         0.00000    -0.04358         0.11006', &
' 9.70        -0.01050         0.05252    -0.02157         0.05503', &
' 9.80         0.00000         0.10504    -0.03258         0.00000', &
' 9.90         0.01050         0.00000    -0.02157         0.16289', &
'10.00         0.00000         0.00000     0.00000         0.16289']
! '#END']

! '#-- Asp-Lys   *', &
CHARACTER(len=80),PARAMETER :: cc2(71) = [&
' 3.00         7.06628       -10.56985     0.00000         0.00000', &
' 3.10         6.00103       -10.56985     0.00000         0.00000', &
' 3.20         4.95231       -10.31655     0.00000         0.00000', &
' 3.30         3.93771        -9.85755     0.00000         0.00000', &
' 3.40         2.98081        -9.22885     0.00000         0.00000', &
' 3.50         2.09194        -8.73425     0.00000         0.00000', &
' 3.60         1.23395        -8.31460     0.00000         0.00000', &
' 3.70         0.42902        -7.08560     0.00000         0.00000', &
' 3.80        -0.18316        -5.27880     0.00000         0.00000', &
' 3.90        -0.62674        -3.83575     0.00000         0.00000', &
' 4.00        -0.95032        -2.91685     0.00000         0.00000', &
' 4.10        -1.21010        -2.04991     0.00000         0.00000', &
' 4.20        -1.36030        -1.32840     0.00000         0.00000', &
' 4.30        -1.47578        -1.08223     0.00000         0.00000', &
' 4.40        -1.57674        -0.88040     0.00000         0.00000', &
' 4.50        -1.65186        -0.64068     0.00000         0.00000', &
' 4.60        -1.70488        -0.38628     0.00000         0.00000', &
' 4.70        -1.72912        -0.12543     0.00000         0.00000', &
' 4.80        -1.72996         0.09524     0.00000         0.00000', &
' 4.90        -1.71007         0.33605     0.00000         0.00000', &
' 5.00        -1.66275         0.57142     0.00000         0.00000', &
' 5.10        -1.59579         0.74565     0.00000         0.00000', &
' 5.20        -1.51362         0.97600     0.00000         0.00000', &
' 5.30        -1.40059         1.10141     0.00000         0.00000', &
' 5.40        -1.29334         1.14693     0.00000         0.00000', &
' 5.50        -1.17121         1.17214     0.00000         0.00000', &
' 5.60        -1.05892         1.10725     0.00000         0.00000', &
' 5.70        -0.94976         0.95580     0.00000         0.00000', &
' 5.80        -0.86776         0.90370     0.00000         0.00000', &
' 5.90        -0.76901         0.93064     0.00000         0.00000', &
' 6.00        -0.68163         0.85402     0.00000         0.00000', &
' 6.10        -0.59821         0.66081     0.00000         0.00000', &
' 6.20        -0.54946         0.43133     0.00000         0.00000', &
' 6.30        -0.51194         0.24305     0.00000         0.00000', &
' 6.40        -0.50086         0.09635     0.00000         0.00000', &
' 6.50        -0.49267        -0.01765     0.00000         0.00000', &
' 6.60        -0.50438        -0.20991     0.00000         0.00000', &
' 6.70        -0.53466        -0.30450     0.00000         0.00000', &
' 6.80        -0.56529        -0.39642     0.00000         0.00000', &
' 6.90        -0.61394        -0.38120     0.00000         0.00000', &
' 7.00        -0.64152        -0.19539     0.00000         0.00000', &
' 7.10        -0.65302        -0.04597     0.00000         0.00000', &
' 7.20        -0.65072         0.07356     0.00000         0.00000', &
' 7.30        -0.63831         0.16091     0.00000         0.00000', &
' 7.40        -0.61854         0.22527     0.00000         0.00000', &
' 7.50        -0.59325         0.26895     0.00000         0.00000', &
' 7.60        -0.56475         0.29534     0.00000         0.00000', &
' 7.70        -0.53418         0.31143     0.00000         0.00000', &
' 7.80        -0.50246         0.32366     0.00000         0.00000', &
' 7.90        -0.46945         0.32596     0.00000         0.00000', &
' 8.00        -0.43727         0.32366     0.00000         0.00000', &
' 8.10        -0.40472         0.32136     0.00000         0.00000', &
' 8.20        -0.37300         0.30658     0.00000         0.00000', &
' 8.30        -0.34340         0.29865     0.00000         0.00000', &
' 8.40        -0.31327         0.28860     0.00000         0.00000', &
' 8.50        -0.28568         0.28044     0.00000         0.00000', &
' 8.60        -0.25718         0.26895     0.00000         0.00000', &
' 8.70        -0.23189         0.25056     0.00000         0.00000', &
' 8.80        -0.20707         0.24072     0.00000         0.00000', &
' 8.90        -0.18375         0.22693     0.00000         0.00000', &
' 9.00        -0.16168         0.21764     0.00000         0.00000', &
' 9.10        -0.14022         0.20155     0.00000         0.00000', &
' 9.20        -0.12137         0.18849     0.00000         0.00000', &
' 9.30        -0.10252         0.17930     0.00000         0.00000', &
' 9.40        -0.08551         0.16781     0.00000         0.00000', &
' 9.50        -0.06896         0.15980     0.00000         0.00000', &
' 9.60        -0.05355         0.15172     0.00000         0.00000', &
' 9.70        -0.03862         0.14788     0.00000         0.00000', &
' 9.80        -0.02398         0.13333     0.00000         0.00000', &
' 9.90        -0.01195         0.11988     0.00000         0.00000', &
'10.00         0.00000         0.11988     0.00000         0.00000']
! '#END'

! '#-- Glu-Arg   *', &
CHARACTER(len=80),PARAMETER :: cc3(71) = [&
' 3.00        14.30590       -26.95335     0.00000         0.00000', &
' 3.10        11.54973       -26.95335     0.00000         0.00000', &
' 3.20         8.91523       -25.32700     0.00000         0.00000', &
' 3.30         6.48433       -22.93310     0.00000         0.00000', &
' 3.40         4.32861       -19.81700     0.00000         0.00000', &
' 3.50         2.52093       -16.06360     0.00000         0.00000', &
' 3.60         1.11589       -12.27245     0.00000         0.00000', &
' 3.70         0.06644        -9.15565     0.00000         0.00000', &
' 3.80        -0.71524        -6.67975     0.00000         0.00000', &
' 3.90        -1.26951        -4.92315     0.00000         0.00000', &
' 4.00        -1.69987        -3.91665     0.00000         0.00000', &
' 4.10        -2.05284        -2.69050     0.00000         0.00000', &
' 4.20        -2.23797        -1.27630     0.00000         0.00000', &
' 4.30        -2.30810        -0.62055     0.00000         0.00000', &
' 4.40        -2.36208        -0.52950     0.00000         0.00000', &
' 4.50        -2.41400        -0.02555     0.00000         0.00000', &
' 4.60        -2.36719         0.84370     0.00000         0.00000', &
' 4.70        -2.24526         1.75045     0.00000         0.00000', &
' 4.80        -2.01710         2.79655     0.00000         0.00000', &
' 4.90        -1.68595         3.37725     0.00000         0.00000', &
' 5.00        -1.34165         3.77635     0.00000         0.00000', &
' 5.10        -0.93068         3.54360     2.33854         4.57745', &
' 5.20        -0.63293         2.30385     0.91549       -12.27565', &
' 5.30        -0.46991         1.84410    -0.11659        -9.05230', &
' 5.40        -0.26411         1.62059    -0.89497        -6.35840', &
' 5.50        -0.14579         0.89648    -1.38827        -3.46340', &
' 5.60        -0.08482         0.54879    -1.58765        -1.92420', &
' 5.70        -0.03603         0.48782    -1.77311        -1.51935', &
' 5.80         0.01275         0.36586    -1.89152         0.00000', &
' 5.90         0.03714         0.30489    -1.77311         1.90510', &
' 6.00         0.07372         0.24391    -1.51050         3.72900', &
' 6.10         0.08592         0.03914    -1.02731         3.62395', &
' 6.20         0.08155         0.00512    -0.78571         2.23225', &
' 6.30         0.08694         0.03033    -0.58086         1.89075', &
' 6.40         0.08762         0.00000    -0.40756         1.33915', &
' 6.50         0.08694        -0.00849    -0.31303         0.84900', &
' 6.60         0.08592        -0.01685    -0.23776         0.74485', &
' 6.70         0.08357        -0.02858    -0.16406         0.66025', &
' 6.80         0.08020        -0.03707    -0.10571         0.58350', &
' 6.90         0.07616        -0.03707    -0.04736         0.58130', &
' 7.00         0.07279        -0.04381     0.01055         0.47016', &
' 7.10         0.06740        -0.04044     0.04667         0.23918', &
' 7.20         0.06470        -0.02935     0.05839         0.07408', &
' 7.30         0.06153        -0.04718     0.06153        -0.04718', &
' 7.40         0.05527        -0.04479     0.05527        -0.04479', &
' 7.50         0.05257        -0.03033     0.05257        -0.03033', &
' 7.60         0.04920        -0.04718     0.04920        -0.04718', &
' 7.70         0.04314        -0.04718     0.04314        -0.04718', &
' 7.80         0.03977        -0.03033     0.03977        -0.03033', &
' 7.90         0.03707        -0.02696     0.03707        -0.02696', &
' 8.00         0.03437        -0.03033     0.03437        -0.03033', &
' 8.10         0.03100        -0.03033     0.03100        -0.03033', &
' 8.20         0.02831        -0.02696     0.02831        -0.02696', &
' 8.30         0.02561        -0.02359     0.02561        -0.02359', &
' 8.40         0.02359        -0.02022     0.02359        -0.02022', &
' 8.50         0.02157        -0.02359     0.02157        -0.02359', &
' 8.60         0.01887        -0.02359     0.01887        -0.02359', &
' 8.70         0.01685        -0.02022     0.01685        -0.02022', &
' 8.80         0.01483        -0.01685     0.01483        -0.01685', &
' 8.90         0.01348        -0.01685     0.01348        -0.01685', &
' 9.00         0.01146        -0.01685     0.01146        -0.01685', &
' 9.10         0.01011        -0.01348     0.01011        -0.01348', &
' 9.20         0.00876        -0.01685     0.00876        -0.01685', &
' 9.30         0.00674        -0.01685     0.00674        -0.01685', &
' 9.40         0.00539        -0.01348     0.00539        -0.01348', &
' 9.50         0.00404        -0.01348     0.00404        -0.01348', &
' 9.60         0.00270        -0.01348     0.00270        -0.01348', &
' 9.70         0.00135        -0.00674     0.00135        -0.00674', &
' 9.80         0.00135        -0.00674     0.00135        -0.00674', &
' 9.90         0.00000        -0.00674     0.00000        -0.00674', &
'10.00         0.00000        -0.00674     0.00000        -0.00674']
! '#END'

! '#-- Glu-Lys   *     ', &
CHARACTER(len=80),PARAMETER :: cc4(71) = [&
' 3.00         7.72863       -14.15915     0.00000         0.00000', &
' 3.10         6.28710       -14.15915     0.00000         0.00000', &
' 3.20         4.89680       -13.46810     0.00000         0.00000', &
' 3.30         3.59348       -12.44885     0.00000         0.00000', &
' 3.40         2.40703       -11.17900     0.00000         0.00000', &
' 3.50         1.35768        -9.69205     0.00000         0.00000', &
' 3.60         0.46862        -7.87195     0.00000         0.00000', &
' 3.70        -0.21671        -5.82070     0.00000         0.00000', &
' 3.80        -0.69552        -4.19355     0.00000         0.00000', &
' 3.90        -1.05542        -3.26620     0.00000         0.00000', &
' 4.00        -1.34876        -2.37295     0.00000         0.00000', &
' 4.10        -1.53001        -1.36522     0.00000         0.00000', &
' 4.20        -1.62180        -0.78705     0.00000         0.00000', &
' 4.30        -1.68742        -0.65616     0.00000         0.00000', &
' 4.40        -1.75304        -0.54680     0.00000         0.00000', &
' 4.50        -1.79678        -0.39967     0.00000         0.00000', &
' 4.60        -1.83297        -0.21872     0.00000         0.00000', &
' 4.70        -1.84052        -0.09245     0.00000         0.00000', &
' 4.80        -1.85146        -0.02734     0.00000         0.00000', &
' 4.90        -1.84599         0.12265     0.00000         0.00000', &
' 5.00        -1.82693         0.21616     0.00000         0.00000', &
' 5.10        -1.80276         0.39681     0.00000         0.00000', &
' 5.20        -1.74757         0.54936     0.00000         0.00000', &
' 5.30        -1.69289         0.68769     0.00000         0.00000', &
' 5.40        -1.61003         0.68764     0.00000         0.00000', &
' 5.50        -1.55536         0.76132     0.00000         0.00000', &
' 5.60        -1.45777         1.01070     0.00000         0.00000', &
' 5.70        -1.35322         1.16883     0.00000         0.00000', &
' 5.80        -1.22400         1.07475     0.00000         0.00000', &
' 5.90        -1.13827         0.79811     0.00000         0.00000', &
' 6.00        -1.06438         0.78752     0.00000         0.00000', &
' 6.10        -0.98077         0.76011     0.00000         0.00000', &
' 6.20        -0.91236         0.64609     0.00000         0.00000', &
' 6.30        -0.85155         0.57008     0.00000         0.00000', &
' 6.40        -0.79834         0.49407     0.00000         0.00000', &
' 6.50        -0.75273         0.49407     0.00000         0.00000', &
' 6.60        -0.69953         0.41806     0.00000         0.00000', &
' 6.70        -0.66912         0.40818     0.00000         0.00000', &
' 6.80        -0.61789         0.45850     0.00000         0.00000', &
' 6.90        -0.57742         0.34250     0.00000         0.00000', &
' 7.00        -0.54939         0.19695     0.00000         0.00000', &
' 7.10        -0.53803         0.00000     0.00000         0.00000', &
' 7.20        -0.54939        -0.12344     0.00000         0.00000', &
' 7.30        -0.56272        -0.14022     0.00000         0.00000', &
' 7.40        -0.57743        -0.12965     0.00000         0.00000', &
' 7.50        -0.58865        -0.11411     0.00000         0.00000', &
' 7.60        -0.60025        -0.08655     0.00000         0.00000', &
' 7.70        -0.60596        -0.08558     0.00000         0.00000', &
' 7.80        -0.61737        -0.07360     0.00000         0.00000', &
' 7.90        -0.62068         0.06382     0.00000         0.00000', &
' 8.00        -0.60461         0.18541     0.00000         0.00000', &
' 8.10        -0.58360         0.37281     0.00000         0.00000', &
' 8.20        -0.53004         0.52476     0.00000         0.00000', &
' 8.30        -0.47865         0.51399     0.00000         0.00000', &
' 8.40        -0.42725         0.49422     0.00000         0.00000', &
' 8.50        -0.37980         0.43492     0.00000         0.00000', &
' 8.60        -0.34026         0.37561     0.00000         0.00000', &
' 8.70        -0.30468         0.35584     0.00000         0.00000', &
' 8.80        -0.26909         0.33607     0.00000         0.00000', &
' 8.90        -0.23746         0.29653     0.00000         0.00000', &
' 9.00        -0.20979         0.31630     0.00000         0.00000', &
' 9.10        -0.17420         0.29653     0.00000         0.00000', &
' 9.20        -0.15048         0.23723     0.00000         0.00000', &
' 9.30        -0.12676         0.21746     0.00000         0.00000', &
' 9.40        -0.10699         0.19769     0.00000         0.00000', &
' 9.50        -0.08722         0.19769     0.00000         0.00000', &
' 9.60        -0.06745         0.20020     0.00000         0.00000', &
' 9.70        -0.04718         0.21011     0.00000         0.00000', &
' 9.80        -0.02543         0.15564     0.00000         0.00000', &
' 9.90        -0.01605         0.12597     0.00000         0.00000', &
'10.00        -0.00024         0.12597     0.00000         0.00000']
! '#END']

WRITE(iout,'(/a)') 'TABULA'

WRITE(iout,*) 4, '! number of potential tables'

! beware of sign minus for forces

! angle-mediated 1st potential: index 1 ! (Asp-Arg)
WRITE(iout,111) 1, SIZE(cc1), spacing,xmin,cutoff, 2, -2.839, -0.299,&
     & 3.295, -0.265, '   ! id, size, spacing, type, alfa1,beta1,alfa2,beta2'
DO i=1, SIZE(cc1)
   READ(cc1(i),*) dum,v1,f1,v2,f2
   v1 = v1 * conversion * scaling_factor
   f1 = - f1 * conversion * scaling_factor
   v2 = v2 * conversion * scaling_factor
   f2 = - f2 * conversion * scaling_factor
   WRITE(iout,*) v1,f1,v2,f2
ENDDO

! straight 1st potential:       index 2 ! (Asp-Lys)
WRITE(iout,111) 2, SIZE(cc2), spacing,xmin,cutoff, 1, 0., 0., 0., 0.,  '   ! id, size, spacing, type'
DO i=1, SIZE(cc2)
   READ(cc2(i),*) dum,v1,f1
   v1 = v1 * conversion * scaling_factor
   f1 = - f1 * conversion * scaling_factor
   WRITE(iout,*) v1,f1
ENDDO

! angle-mediated 2nd potential: index 3 ! (Glu-Arg)
WRITE(iout,111) 3, SIZE(cc3), spacing,xmin,cutoff, 2, -1.500, -0.387,&
     & 2.034, 0.030, '   ! id, size, spacing, type, alfa1,beta1,alfa2,beta2'
DO i=1, SIZE(cc3)
   READ(cc3(i),*) dum,v1,f1,v2,f2
   v1 = v1 * conversion * scaling_factor
   f1 = - f1 * conversion * scaling_factor
   v2 = v2 * conversion * scaling_factor
   f2 = - f2 * conversion * scaling_factor
   WRITE(iout,*) v1,f1,v2,f2
ENDDO

! straight 2nd potential:       index 4 ! (Glu-Lys)
WRITE(iout,111) 4, SIZE(cc4), spacing,xmin,cutoff, 1, 0., 0., 0., 0., '   ! id, size, spacing, type'
DO i=1, SIZE(cc4)
   READ(cc4(i),*) dum,v1,f1
   v1 = v1 * conversion * scaling_factor
   f1 = - f1 * conversion * scaling_factor
   WRITE(iout,*) v1,f1
ENDDO
111 FORMAT(i3,i6,3f12.5,i6,4f8.3,a)
112 FORMAT(i3,i6,3f12.5,i6,a)

WRITE(iout,'(a/)') 'CLOSETABULA'

END SUBROUTINE ip_potential_dat



END program converter
