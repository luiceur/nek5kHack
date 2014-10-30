# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
c-----------------------------------------------------------------------
c     
c    Stability limits:
c     
c    AB3:    .7236                     w/safety (1.2):   .603
c     
c    RK3:    1.73   (sqrt 3)           w/safety (1.2):   1.44
c     
c    RK4:    2.828                     w/safety (1.2):   2.36
c     
c    SEM Safety factor:  1.52 for N=3
c                     <  1.20 for N=16
c                     ~  1.16 for N=256
c     
c-----------------------------------------------------------------------
      subroutine setup_convect(igeom)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 18 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 18 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 19 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 19 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      logical ifnew
      
      common /cchar/ ct_vx(0:lorder+1) ! time for each slice in c_vx()
      
      if (igeom.eq.1) return
      if (param(99).lt.0) return ! no dealiasing
      
      if (ifchar) then
      
         nelc = nelv
         if (ifmhd) nelc = max(nelv,nelfld(ifldmhd))
         if (ifmhd) call exitti('no characteristics for mhd yet$',istep)
      
         ifnew = .true.
         if (igeom.gt.2) ifnew = .false.
      
         call set_conv_char(ct_vx,c_vx,vx,vy,vz,nelc,time,ifnew)
      
      else
      
         if (.not.ifpert) then
            if (ifcons) then
               call set_convect_cons (vxd,vyd,vzd,vx,vy,vz)
            else
               call set_convect_new  (vxd,vyd,vzd,vx,vy,vz)
            endif
         endif
      
      endif
      
c     write(6,*) istep,' conv',ifnew,igeom,' continu? ',time
c     read(5,*) dum
      
      return
      end
c-----------------------------------------------------------------------
      subroutine char_conv(p0,u,ulag,msk,c,cs,gsl)
c     
c     
c     Convect over last NBD steps using characteristics scheme
c     
c     NOTE:  Here, we assume that ulag is stored by time-slice first,
c            then by field number (this is opposite to prior Nek5000)
c     
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 65 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 65 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 66 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 66 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      real    p0(1),u(1),ulag(1),msk(1),c(1),cs(0:1)
      integer gsl
      
      common /scrns/ ct  (lxd*lyd*lzd*lelv*ldim)
      
      common /scrvh/ bmsk(lx1*ly1*lz1*lelv)
     $             , u1  (lx1*ly1*lz1*lelv)
      
      common /scrmg/ r1  (lx1*ly1*lz1*lelv)
     $             , r2  (lx1*ly1*lz1*lelv)
     $             , r3  (lx1*ly1*lz1*lelv)
     $             , r4  (lx1*ly1*lz1*lelv)
      
c     
      nelc = nelv            ! number of elements in convecting field
      if (ifield.eq.ifldmhd) nelc = nelfld(ifield)
      
      nc  = cs(0)            ! number of stored convecting fields
      
      ln  = lx1*ly1*lz1*lelt
      n   = nx1*ny1*nz1*nelfld(ifield)
      m   = nxd*nyd*nzd*nelc*ndim
      
      if (ifield.eq.ifldmhd) then
         call col3(bmsk,bintm1,msk,n)
      elseif (ifield.eq.1) then
         call col3(bmsk,binvm1,msk,n)
      else ! if (ifield.eq.2) then
         call col3(bmsk,bintm1,msk,n)
      endif
      
      call char_conv1
     $   (p0,bmsk,u,n,ulag,ln,gsl,c,m,cs(1),nc,ct,u1,r1,r2,r3,r4)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine char_conv1
     $   (p0,bmsk,u,n,ulag,ln,gsl,c,m,cs,nc,ct,u1,r1,r2,r3,r4)
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 107 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 107 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 108 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 108 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 109 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 109 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real    p0(n),u(n),ulag(ln,1),bmsk(n),c(m,0:nc),cs(0:nc)
      
      real    ct(m),u1(n),r1(n),r2(n),r3(n),r4(n) ! work arrays
      
      integer gsl
      
      
!     Convect over last NBD steps using characteristics scheme
      
!              n-q                                      n-1
!     Given u(t    , X ),  q = 1,2,...,nbd, compute  phi     ( := p0 )
      
!        n-1       nbd   ~n-q
!     phi     :=  sum    u
!                  q=1
      
!          ~n-q             ~n-q
!     each u     satisfies  u    := v  such that
      
      
!     dv
!     -- + C.grad v = 0  t \in [t^n-q,t^n],   v(t^n-q,X) = u(t^n-q,X)
!     dt
      
      
      tau = time-vlsum(dtlag,nbd)         ! initialize time for u^n-k
      call int_vel (ct,tau,c,m,nc,cs,nid) ! ct(t) = sum w_k c(.,k)
      
      call rzero(p0,n)
      
      do ilag = nbd,1,-1
      
         um = 0
         if (ilag.eq.1) then
            do i=1,n
               p0(i) = p0(i)+bd(ilag+1)*u(i)
               um=max(um,u(i))
            enddo
         else
            do i=1,n
               p0(i) = p0(i)+bd(ilag+1)*ulag(i,ilag-1)
               um=max(um,ulag(i,ilag-1))
            enddo
         endif
      
c        write(6,1) istep,ilag,bd(ilag),bd(ilag+1),um
c 1      format(i5,i4,1p3e14.5,' bdf')
      
         dtau = dtlag(ilag)/ntaubd
         do itau = 1,ntaubd ! ntaubd=number of RK4 substeps (typ. 1 or 2
      
            tau1 = tau + dtau
      
            c1 = 1.
            c2 = -dtau/2.
            c3 = -dtau
            th = tau+dtau/2.
      
            call conv_rhs(r1,p0,ct,bmsk,gsl)         !  STAGE 1
      
            call int_vel (ct,th,c,m,nc,cs,nid)       !  STAGE 2
            call add3s12 (u1,p0,r1,c1,c2,n)
            call conv_rhs(r2,u1,ct,bmsk,gsl)
      
            call add3s12 (u1,p0,r2,c1,c2,n)          !  STAGE 3
            call conv_rhs(r3,u1,ct,bmsk,gsl)
      
            call int_vel (ct,tau1,c,m,nc,cs,nid)     !  STAGE 4
            call add3s12 (u1,p0,r3,c1,c3,n)
            call conv_rhs(r4,u1,ct,bmsk,gsl)
      
            c1 = -dtau/6.
            c2 = -dtau/3.
            do i=1,n
               p0(i) = p0(i)+c1*(r1(i)+r4(i))+c2*(r2(i)+r3(i))
            enddo
            tau = tau1
         enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine int_vel(c_t,t0,c,n,nc,ct,nid)
      
c     Interpolate convecting velocity field c(t_1,...,t_nconv) to
c     time t0 and return result in c_t.
      
c     Ouput:   c_t = sum wt_k * ct_i(k) 
      
c     Here, t0 is the time of interest
      
      real c_t(n),c(n,0:nc),ct(0:nc)
c     
      parameter (lwtmax=10)
      real wt(0:lwtmax)
c     
      if (nc.gt.lwtmax) then
         write(6,*) nid,'ERROR int_vel: lwtmax too small',lwtmax,m0
         call exitt
      endif
c     
      no = nc-1
      call fd_weights_full(t0,ct(0),no,0,wt)  ! interpolation weights
      
      call rzero(c_t,n)
      do j=1,n
      do i=0,no
         c_t(j) = c_t(j) + wt(i)*c(j,i)
      enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine conv_rhs (du,u,c,bmsk,gsl)
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 228 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 228 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 229 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 229 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
c     
c     apply convecting field c(1,ndim) to scalar field u(1)
c     
      real du(1),u(1),c(1),bmsk(1)
      integer gsl
c     
      logical ifconv
c     
c     ifconv = .false.
      ifconv = .true.
c     
      n = nx1*ny1*nz1*nelv
c     
      if (ifconv) then
      
         if (ifcons) then
           if (if3d     ) call convop_cons_3d (du,u,c,nx1,nxd,nelv)
           if (.not.if3d) call convop_cons_2d (du,u,c,nx1,nxd,nelv)
         else
           if (if3d     ) call convop_fst_3d  (du,u,c,nx1,nxd,nelv)
           if (.not.if3d) call convop_fst_2d  (du,u,c,nx1,nxd,nelv)
         endif
      
         call gs_op(gsl,du,1,1,0)  !  +
      
      else
         call rzero   (du,n)
         return
      endif
c     
      do i=1,n
         du(i) = bmsk(i)*du(i)  ! Binv * msk
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine convop_fst_3d(du,u,c,mx,md,nel)
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 269 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 269 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
c     
c     apply convecting field c to scalar field u
c     
      real du(mx*mx*mx,nel)
      real  u(mx*mx*mx,nel)
      real  c(md*md*md,nel,3)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd)
c     
      logical if3d,ifd
      integer e
c     
      if3d = .true.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.
      ifd=.true.
c     
      nrstd = md**3
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convop_fst')
c     
      do e=1,nel
         call grad_rstd(ur,us,ut,u(1,e),mx,md,if3d,ud)
         if (ifd) then    ! dealiased
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DF
               ud(i) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)+c(i,e,3)*ut(i)
            enddo
            call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
         else
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DF
               du(i,e) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)+c(i,e,3)*ut(i)
            enddo
         endif
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine convop_fst_2d(du,u,c,mx,md,nel)
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 311 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 311 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
c     
c     apply convecting field c to scalar field u
c     
      real du(mx*mx,nel)
      real  u(mx*mx,nel)
      real  c(md*md,nel,2)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ud(ldd)
c     
      logical if3d,ifd
      integer e
c     
      if3d = .false.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.
      ifd=.true.
c     
      nrstd = md**2
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convop_fst')
c     
      do e=1,nel
         call grad_rstd(ur,us,ut,u(1,e),mx,md,if3d,ud)
         if (ifd) then    ! dealiased
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DF
               ud(i) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)
            enddo
            call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
         else
            do i=1,nrstd
c              C has the mass matrix factored in per (4.8.5), p. 227, DF
               du(i,e) = c(i,e,1)*ur(i)+c(i,e,2)*us(i)
            enddo
         endif
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_rstd(ur,us,ut,u,mx,md,if3d,ju) ! GLL->GL grad
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 353 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 353 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 354 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 354 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real    ur(1),us(1),ut(1),u(1),ju(1)
      logical if3d
      
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
      
      call intp_rstd(ju,u,mx,md,if3d,0) ! 0 = forward
      
      m0 = md-1
      call get_dgl_ptr (ip,md,md)
      if (if3d) then
         call local_grad3(ur,us,ut,ju,m0,1,dg(ip),dgt(ip))
      else
         call local_grad2(ur,us   ,ju,m0,1,dg(ip),dgt(ip))
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine intp_rstd(ju,u,mx,md,if3d,idir) ! GLL->GL interpolation
      
c     GLL interpolation from mx to md.
      
c     If idir ^= 0, then apply transpose operator  (md to mx)
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 383 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 383 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real    ju(1),u(1)
      logical if3d
      
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
      
      parameter (ld=2*lxd)
      common /ctmp0/ w(ld**ldim,2)
      
      call lim_chk(md,ld,'md   ','ld   ','grad_rstd ')
      call lim_chk(mx,ld,'mx   ','ld   ','grad_rstd ')
      
      ldw = 2*(ld**ldim)
      
      call get_int_ptr (i,mx,md)
c     
      if (idir.eq.0) then
         call specmpn(ju,md,u,mx,jgl(i),jgt(i),if3d,w,ldw)
      else
         call specmpn(ju,mx,u,md,jgt(i),jgl(i),if3d,w,ldw)
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_int(jgl,jgt,mp,np,w)
c     
c     Generate interpolation from np GLL points to mp GL points
c     
c        jgl  = interpolation matrix, mapping from velocity nodes to pre
c        jgt  = transpose of interpolation matrix
c        w    = work array of size (np+mp)
c     
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c     
c     
      real jgl(mp,np),jgt(np*mp),w(1)
c     
      iz = 1
      id = iz + np
c     
      call zwgll (w(iz),jgt,np)
      call zwgl  (w(id),jgt,mp)
c     
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo
c     
      call transpose(jgt,np,jgl,mp)
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_dgl(dgl,dgt,mp,np,w)
c     
c     Generate derivative from np GL points onto mp GL points
c     
c        dgl  = interpolation matrix, mapping from velocity nodes to pre
c        dgt  = transpose of interpolation matrix
c        w    = work array of size (3*np+mp)
c     
c        np   = number of points on GLL grid
c        mp   = number of points on GL  grid
c     
c     
c     
      real dgl(mp,np),dgt(np*mp),w(1)
c     
c     
      iz = 1
      id = iz + np
c     
      call zwgl  (w(iz),dgt,np)  ! GL points
      call zwgl  (w(id),dgt,mp)  ! GL points
c     
      ndgt = 2*np
      ldgt = mp*np
      call lim_chk(ndgt,ldgt,'ldgt ','dgt  ','gen_dgl   ')
c     
      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,1,dgt) ! 1=1st deriv.
         do j=1,np
            dgl(i,j) = dgt(np+j)                       ! Derivative matr
         enddo
      enddo
c     
      call transpose(dgt,np,dgl,mp)
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine lim_chk(n,m,avar5,lvar5,sub_name10)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 485 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 485 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      character*5  avar5,lvar5
      character*10 sub_name10
      
      if (n.gt.m) then
         write(6,1) nid,n,m,avar5,lvar5,sub_name10
    1    format(i8,' ERROR: :',2i12,2(1x,a5),1x,a10)
         call exitti('lim_chk problem. $',n)
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine get_int_ptr (ip,mx,md) ! GLL-->GL pointer
      
c     Get pointer to jgl() for interpolation pair (mx,md)
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 502 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 502 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c     
      parameter (ld=2*lxd)
      common /igrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c     
      ij = md + ld*(mx-1)
      ip = pjgl(ij)
c     
      if (ip.eq.0) then
c     
         nstore   = pjgl(0)
         pjgl(ij) = nstore+1
         nstore   = nstore + md*mx
         pjgl(0)  = nstore
         ip       = pjgl(ij)
c     
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'jgl  ','ldg  ','get_int_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_int_pt')
c     
         call gen_int(jgl(ip),jgt(ip),md,mx,wkd)
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine get_dgl_ptr (ip,mx,md)
c     
c     Get pointer to GL-GL interpolation dgl() for pair (mx,md)
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 540 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 540 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
c     
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
c     
      parameter (ld=2*lxd)
      common /jgrad/ pd    (0:ld*ld)
     $             , pdg   (0:ld*ld)
     $             , pjgl  (0:ld*ld)
      integer pd , pdg , pjgl
c     
      ij = md + ld*(mx-1)
      ip = pdg (ij)
      
      if (ip.eq.0) then
      
         nstore   = pdg (0)
         pdg (ij) = nstore+1
         nstore   = nstore + md*mx
         pdg (0)  = nstore
         ip       = pdg (ij)
c     
         nwrkd = mx + md
         call lim_chk(nstore,ldg ,'dg   ','ldg  ','get_dgl_pt')
         call lim_chk(nwrkd ,lwkd,'wkd  ','lwkd ','get_dgl_pt')
c     
         call gen_dgl(dg (ip),dgt(ip),md,mx,wkd)
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine set_conv_char(ct,c,ux,uy,uz,nelc,tau,ifnew)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 575 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 575 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 576 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 576 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real ct(0:1)               ! time stamps for saved field (0=#flds)
      real c(1)                  ! saved vel. fields, dealiased etc.
      real ux(1),uy(1),uz(1)     ! input vel. field
      integer nelc               ! number of elements in conv. field
      logical ifnew              ! =true if shifting stack of fields
      
      numr      = lxd*lyd*lzd*lelv*ldim*(lorder+1)
      denr      = nxd*nyd*nzd*nelv*ndim
      nconv_max = numr/denr
      if (nconv_max.lt.nbdinp+1) 
     $   call exitti(
     $     'ABORT: not enough memory for characteristics scheme!$',
     $     nconv_max)
      
      nc = ct(0)
      
      m  = nxd*nyd*nzd*nelc*ndim
      
c     write(6,*) nelc,ifnew,' set conv_char',istep,nc,nconv_max
      call set_ct_cvx
     $    (ct,c,m,ux,uy,uz,tau,nc,nconv_max,nelc,ifnew)
      
      nc = min (nc,nbdinp)
      ct(0) = nc  ! store current count
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_ct_cvx(ct,c,m,u,v,w,tau,nc,mc,nelc,ifnew)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 607 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 607 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 608 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 608 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real ct(0:1),c(m,1)
      real u(1),v(1),w(1)
      logical ifnew
      
      if (ifnew) then
      
c        Shift existing convecting fields
c        Note:  "1" entry is most recent
      
         nc = nc+1
         nc = min(nc,mc)
         ct(0) = nc
      
         do i=nc,2,-1
            call copy(c(1,i),c(1,i-1),m)
            ct(i) = ct(i-1)
         enddo
      endif
      
c     Save time and map the current velocity to rst coordinates.
      
      ix = 1
      iy = ix + nxd*nyd*nzd*nelc
      iz = iy + nxd*nyd*nzd*nelc
      
      if (ifcons) then
         call set_convect_cons(c(ix,1),c(iy,1),c(iz,1),u,v,w)
      else
         call set_convect_new (c(ix,1),c(iy,1),c(iz,1),u,v,w)
      endif
      
      ct(1) = tau
      
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_rst(ur,us,ut,u,md,if3d) ! Gauss-->Gauss grad
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 648 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 648 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 649 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 649 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real    ur(1),us(1),ut(1),u(1)
      logical if3d
      
      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt
      
      m0 = md-1
      call get_dgl_ptr (ip,md,md)
      if (if3d) then
         call local_grad3(ur,us,ut,u,m0,1,dg(ip),dgt(ip))
      else
         call local_grad2(ur,us   ,u,m0,1,dg(ip),dgt(ip))
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine convect_new(bdu,u,ifuf,cx,cy,cz,ifcf)
      
C     Compute dealiased form:  J^T Bf *JC .grad Ju w/ correct Jacobians
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 674 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 674 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 675 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 675 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real bdu(1),u(1),cx(1),cy(1),cz(1)
      logical ifuf,ifcf            ! u and/or c already on fine mesh?
      
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)
      
      integer e
      
      call set_dealias_rx
      
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      
      nxyzu = nxyz1
      if (ifuf) nxyzu = nxyzd
      
      nxyzc = nxyz1
      if (ifcf) nxyzc = nxyzd
      
      iu = 1    ! pointer to scalar field u
      ic = 1    ! pointer to vector field C
      ib = 1    ! pointer to scalar field Bdu
      
      
      do e=1,nelv
      
         if (ifcf) then
      
            call copy(tr(1,1),cx(ic),nxyzd)  ! already in rst form
            call copy(tr(1,2),cy(ic),nxyzd)
            if (if3d) call copy(tr(1,3),cz(ic),nxyzd)
      
         else  ! map coarse velocity to fine mesh (C-->F)
      
           call intp_rstd(fx,cx(ic),nx1,nxd,if3d,0) ! 0 --> forward
           call intp_rstd(fy,cy(ic),nx1,nxd,if3d,0) ! 0 --> forward
           if (if3d) call intp_rstd(fz,cz(ic),nx1,nxd,if3d,0) ! 0 --> fo
      
           if (if3d) then  ! Convert convector F to r-s-t coordinates
      
             do i=1,nxyzd
               tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
               tr(i,2)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
               tr(i,3)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
             enddo
      
           else
      
             do i=1,nxyzd
               tr(i,1)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
               tr(i,2)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
             enddo
      
           endif
      
         endif
      
         if (ifuf) then
            call grad_rst(ur,us,ut,u(iu),nxd,if3d)
         else
            call intp_rstd(uf,u(iu),nx1,nxd,if3d,0) ! 0 --> forward
            call grad_rst(ur,us,ut,uf,nxd,if3d)
         endif
      
         if (if3d) then
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)+tr(i,3)*ut(i)
            enddo
         else
            do i=1,nxyzd ! mass matrix included, per DFM (4.8.5)
               uf(i) = tr(i,1)*ur(i)+tr(i,2)*us(i)
            enddo
         endif
         call intp_rstd(bdu(ib),uf,nx1,nxd,if3d,1) ! Project back to coa
      
         ic = ic + nxyzc
         iu = iu + nxyzu
         ib = ib + nxyz1
      
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine convect_cons(bdu,u,ifuf,cx,cy,cz,ifcf)
      
c     Compute dealiased form:  J^T Bf *div. JC Ju w/ correct Jacobians
      
c     conservative form
      
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 770 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 770 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 771 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 771 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real bdu(1),u(1),cx(1),cy(1),cz(1)
      
      logical ifuf,ifcf            ! u and/or c already on fine mesh?
      
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      common /scrcv/ uf(ltd),cf(ltd),cu(ltd)
     $             , cr(ltd),cs(ltd),ct(ltd)
      
      
      integer e
      
      call set_dealias_rx
      
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      
      nxyzu = nxyz1
      if (ifuf) nxyzu = nxyzd
      
      nxyzc = nxyz1
      if (ifcf) nxyzc = nxyzd
      
      iu = 1    ! pointer to scalar field u
      ic = 1    ! pointer to vector field C
      ib = 1    ! pointer to scalar field Bdu
      
      do e=1,nelv
      
        call intp_rstd(uf,u(iu),nx1,nxd,if3d,0) ! 0 --> forward
      
        call rzero(cu,nxyzd)
        do i=1,ndim
      
         if (ifcf) then  ! C is already on fine mesh
      
           call exitt  ! exit for now
      
         else  ! map coarse velocity to fine mesh (C-->F)
      
           if (i.eq.1) call intp_rstd(cf,cx(ic),nx1,nxd,if3d,0) ! 0 --> 
           if (i.eq.2) call intp_rstd(cf,cy(ic),nx1,nxd,if3d,0) ! 0 --> 
           if (i.eq.3) call intp_rstd(cf,cz(ic),nx1,nxd,if3d,0) ! 0 --> 
      
           call col2(cf,uf,nxyzd)   !  collocate C and u on fine mesh
      
           call grad_rst(cr,cs,ct,cf,nxd,if3d)  ! d/dr (C_i*u)
      
           if (if3d) then
      
             do j=1,nxyzd
               cu(j)=cu(j)
     $              +cr(j)*rx(j,i,e)+cs(j)*rx(j,i+3,e)+ct(j)*rx(j,i+6,e)
             enddo
      
           else  ! 2D
      
             do j=1,nxyzd
               cu(j)=cu(j)
     $              +cr(j)*rx(j,i,e)+cs(j)*rx(j,i+2,e)
             enddo
      
           endif
         endif
        enddo
      
        call intp_rstd(bdu(ib),cu,nx1,nxd,if3d,1) ! Project back to coar
      
        ic = ic + nxyzc
        iu = iu + nxyzu
        ib = ib + nxyz1
      
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_convect_cons(cx,cy,cz,ux,uy,uz)
      
c     Put vx,vy,vz on fine mesh (for conservation form)
      
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 854 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 854 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 855 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 855 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      
      real cx(ltd,1),cy(ltd,1),cz(ltd,1)
      real ux(lxy,1),uy(lxy,1),uz(lxy,1)
      
      integer e
      
      call set_dealias_rx
      
      do e=1,nelv    ! Map coarse velocity to fine mesh (C-->F)
      
         call intp_rstd(cx(1,e),ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         call intp_rstd(cy(1,e),uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(cz(1,e),uz(1,e),nx1,nxd,if3d,0) ! 0 --
      
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_convect_new(cr,cs,ct,ux,uy,uz)
C     
C     Put vxd,vyd,vzd into rst form on fine mesh
C     
C     For rst form, see eq. (4.8.5) in Deville, Fischer, Mund (2002).
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 883 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 883 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 884 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 884 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      parameter (lxy=lx1*ly1*lz1,ltd=lxd*lyd*lzd)
      
      real cr(ltd,1),cs(ltd,1),ct(ltd,1)
      real ux(lxy,1),uy(lxy,1),uz(lxy,1)
      
      common /scrcv/ fx(ltd),fy(ltd),fz(ltd)
     $             , ur(ltd),us(ltd),ut(ltd)
     $             , tr(ltd,3),uf(ltd)
      
      integer e
      
      call set_dealias_rx
      
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      
      ic = 1    ! pointer to vector field C
      
      do e=1,nelv 
      
c        Map coarse velocity to fine mesh (C-->F)
      
         call intp_rstd(fx,ux(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         call intp_rstd(fy,uy(1,e),nx1,nxd,if3d,0) ! 0 --> forward
         if (if3d) call intp_rstd(fz,uz(1,e),nx1,nxd,if3d,0) ! 0 --> for
      
c        Convert convector F to r-s-t coordinates
      
         if (if3d) then
      
           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)+rx(i,3,e)*fz(i)
              cs(i,e)=rx(i,4,e)*fx(i)+rx(i,5,e)*fy(i)+rx(i,6,e)*fz(i)
              ct(i,e)=rx(i,7,e)*fx(i)+rx(i,8,e)*fy(i)+rx(i,9,e)*fz(i)
           enddo
      
         else
      
           do i=1,nxyzd
              cr(i,e)=rx(i,1,e)*fx(i)+rx(i,2,e)*fy(i)
              cs(i,e)=rx(i,3,e)*fx(i)+rx(i,4,e)*fy(i)
           enddo
      
         endif
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_char_mask(mask,u,v,w) ! mask for hyperbolic system 
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 937 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 937 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 938 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 938 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 939 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 939 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 940 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 940 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 941 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 941 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      integer msk(0:1)
      character      cb*3
      parameter (lxyz1=lx1*ly1*lz1)
      common /ctmp1/ work(lxyz1,lelt)
      
      real mask(lxyz1,1),u(lxyz1,1),v(lxyz1,1),w(lxyz1,1)
      
      integer e,f
      
      nfaces= 2*ndim
      ntot1 = nx1*ny1*nz1*nelv
      call rzero (work,ntot1)
      call rone  (mask,NTOT1)
      
      ifldv = 1
      do 100 e=1,nelv
      do 100 f=1,nfaces
         cb=cbc(f,e,ifldv)
         if (cb(1:1).eq.'v' .or. cb(1:1).eq.'V') then
      
           call faccl3 (work(1,e),u(1,e),unx(1,1,f,e),f)
           call faddcl3(work(1,e),v(1,e),uny(1,1,f,e),f)
           if (if3d) 
     $     call faddcl3(work(1,e),w(1,e),unz(1,1,f,e),f)
      
           call fcaver (vaver,work,e,f)
      
           if (vaver.lt.0) call facev (mask,e,f,0.0,nx1,ny1,nz1)
         endif
         if (cb(1:2).eq.'ws' .or. cb(1:2).eq.'WS') 
     $   call facev (mask,e,f,0.0,nx1,ny1,nz1)
 100  continue
      
      return
      end
c-----------------------------------------------------------------------
      subroutine advchar
c     
c     Compute convective contribution using 
c     operator-integrator-factor method (characteristics).
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 983 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 983 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 984 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 984 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 985 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 985 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 986 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 986 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 987 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 987 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 988 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 988 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/CTIMER" 1
C     
# 2
      COMMON /CTIMER/ tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
      COMMON /CTIME2/ tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep
     $               ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2
     $               ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn
     $               ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee
     $               ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc
     $               ,twal
C     
      REAL*8          tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
      REAL*8          tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep
     $               ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2
     $               ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn
     $               ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee
     $               ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc
     $               ,twal
C     
      COMMON /ITIMER/ nmxmf,nmxms,ndsum,naxhm,ncopy,ninvc,ninv3
      COMMON /ITIME2/ nsolv,ngsum,ndsnd,ndadd,ncdtp,nmltd,nprep
     $               ,npres,nhmhz,ngop ,ngop1,ndott,nbsol,nbso2
     $               ,nsett,nslvb,nusbc,nddsl,ncrsl,ndsmx,ndsmn
     $               ,ngsmn,ngsmx,neslv,nbbbb,ncccc,ndddd,neeee
     $               ,nvdss,nadvc,nspro,ngop_sync,nsyc,nwal
C     
      COMMON /PTIMER/ pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
     $               ,psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep
     $               ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2
     $               ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn
     $               ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee
     $               ,pvdss,pspro,pgop_sync,psyc,pwal
C     
      REAL*8          pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
      REAL*8          psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep
     $               ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2
     $               ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn
     $               ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee
     $               ,pvdss,pspro,pgop_sync,psyc,pwal
C     
      REAL*8 etime1,etime2,etime0,gtime1,tscrtch
      REAL*8 dnekclock,dnekclock_sync
C     
      COMMON /CTIME3/ etimes,ttotal,tttstp,etims0,ttime
      real*8          etimes,ttotal,tttstp,etims0,ttime
C     
      integer icalld
      save    icalld
      data    icalld /0/
      
      common /ctimel/ ifsync
      logical         ifsync
# 990 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 990 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      common /cchar/ ct_vx(0:lorder) ! time for each slice in c_vx()
      
      common /scruz/ phx  (lx1*ly1*lz1*lelt)
     $ ,             phy  (lx1*ly1*lz1*lelt)
     $ ,             phz  (lx1*ly1*lz1*lelt)
     $ ,             hmsk (lx1*ly1*lz1*lelt)
      

# 999
      if (icalld.eq.0) tadvc=0.0
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()

      
# 1005
      dti = 1./dt
      n   = nx1*ny1*nz1*nelv
      
      call set_char_mask(hmsk,vx,vy,vz) ! mask for hyperbolic system 
      
      call char_conv(phx,vx,vxlag,hmsk,c_vx,ct_vx,gsh_fld(1))
      call char_conv(phy,vy,vylag,hmsk,c_vx,ct_vx,gsh_fld(1))
      
      if (if3d) then
      
        call char_conv(phz,vz,vzlag,hmsk,c_vx,ct_vx,gsh_fld(1))
      
        do i=1,n
           h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
           bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
           bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
           bfz(i,1,1,1) = bfz(i,1,1,1)+phz(i)*h2i
        enddo
      
      else
      
        do i=1,n
           h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
           bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
           bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
        enddo
      
      endif
      

# 1035
      tadvc=tadvc+(dnekclock()-etime1)

# 1037
      return
      end
c-----------------------------------------------------------------------
      subroutine convch
      
c     Compute convective contribution using 
c     operator-integrator-factor method (characteristics).
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 1046 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1046 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 1047 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1047 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 1048 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1048 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 1049 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1049 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 1050 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1050 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 1051 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1051 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/CTIMER" 1
C     
# 2
      COMMON /CTIMER/ tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
      COMMON /CTIME2/ tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep
     $               ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2
     $               ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn
     $               ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee
     $               ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc
     $               ,twal
C     
      REAL*8          tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
      REAL*8          tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep
     $               ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2
     $               ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn
     $               ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee
     $               ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc
     $               ,twal
C     
      COMMON /ITIMER/ nmxmf,nmxms,ndsum,naxhm,ncopy,ninvc,ninv3
      COMMON /ITIME2/ nsolv,ngsum,ndsnd,ndadd,ncdtp,nmltd,nprep
     $               ,npres,nhmhz,ngop ,ngop1,ndott,nbsol,nbso2
     $               ,nsett,nslvb,nusbc,nddsl,ncrsl,ndsmx,ndsmn
     $               ,ngsmn,ngsmx,neslv,nbbbb,ncccc,ndddd,neeee
     $               ,nvdss,nadvc,nspro,ngop_sync,nsyc,nwal
C     
      COMMON /PTIMER/ pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
     $               ,psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep
     $               ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2
     $               ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn
     $               ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee
     $               ,pvdss,pspro,pgop_sync,psyc,pwal
C     
      REAL*8          pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
      REAL*8          psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep
     $               ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2
     $               ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn
     $               ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee
     $               ,pvdss,pspro,pgop_sync,psyc,pwal
C     
      REAL*8 etime1,etime2,etime0,gtime1,tscrtch
      REAL*8 dnekclock,dnekclock_sync
C     
      COMMON /CTIME3/ etimes,ttotal,tttstp,etims0,ttime
      real*8          etimes,ttotal,tttstp,etims0,ttime
C     
      integer icalld
      save    icalld
      data    icalld /0/
      
      common /ctimel/ ifsync
      logical         ifsync
# 1052 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1052 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      common /cchar/ ct_vx(0:lorder) ! time for each slice in c_vx()
      
      common /scruz/ phi  (lx1*ly1*lz1*lelt)
     $ ,             hmsk (lx1*ly1*lz1*lelt)
      
      if (icalld.eq.0) tadvc=0.0
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()
      
      n   = nx1*ny1*nz1*nelv
      dti = 1./dt
      
      if (ifield.eq.2) then  ! set convecting velocity and mask
c        call setup_convect(1)
         call set_char_mask(hmsk,vx,vy,vz) ! mask for hyperbolic system 
      endif
      
      
      call char_conv(phi,t(1,1,1,1,ifield-1),tlag(1,1,1,1,1,ifield-1)
     $        ,hmsk,c_vx,ct_vx,gsh_fld(1))
      
c     pmax = glamax(phi,n)
c     qmax = glamax(vtrans(1,1,1,1,2),n)
c     write(6,*) istep,dti,pmax,' pmax'
      
      do i=1,n
         bq(i,1,1,1,ifield-1) = bq(i,1,1,1,ifield-1)
     $          + phi(i)*bm1(i,1,1,1)*vtrans(i,1,1,1,ifield)*dti
      enddo
      
      tadvc=tadvc+(dnekclock()-etime1)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine convop_cons_3d(du,u,c,mx,md,nel) ! Conservation form
      
c     Apply convecting field c to scalar field u, conservation form d/dx
      
c     Assumes that current convecting field is on dealias mesh, in c()
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 1096 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1096 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 1097 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1097 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 1098 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1098 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      real du(mx*mx*mx,nel)
      real  u(mx*mx*mx,nel)
      real  c(md*md*md,nel,3)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      
      logical if3d,ifd
      integer e
      
      if3d = .true.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.
      ifd=.true.
      
      nrstd = md**3
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convp_cons')
      
      do e=1,nel
      
         call intp_rstd (ju,u(1,e),mx,md,if3d,0) ! 0 = forward; on Gauss
         call rzero     (ud,nrstd)
      
         do j=1,ndim
            do i=1,nrstd
               tu(i)=c(i,e,j)*ju(i)   ! C_j*T
            enddo
            call grad_rst(ur,us,ut,tu,md,if3d)  ! Already on fine (Gauss
      
            j0 = j+0
            j3 = j+3
            j6 = j+6
            do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mes
               ud(i)=ud(i)
     $              +rx(i,j0,e)*ur(i)+rx(i,j3,e)*us(i)+rx(i,j6,e)*ut(i)
            enddo
         enddo
      
         call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
      
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine convop_cons_2d(du,u,c,mx,md,nel) ! Conservation form
      
c     Apply convecting field c to scalar field u, conservation form d/dx
      
c     Assumes that current convecting field is on dealias mesh, in c()
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 1151 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1151 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 1152 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1152 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 1153 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1153 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      
      real du(mx*mx,nel)
      real  u(mx*mx,nel)
      real  c(md*md,nel,2)
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju
      
      logical if3d,ifd
      integer e
      
      if3d = .false.
      ifd  = .false.
      if (md.ne.mx) ifd=.true.
      
      nrstd = md**2
      call lim_chk(nrstd,ldd,'urus ','ldd  ','convp_cons')
      
      if (nid.eq.0.and.istep.lt.3) write(6,*) 'convp_cons',istep
      
      do e=1,nel
      
         call intp_rstd (ju,u(1,e),mx,md,if3d,0) ! 0 = forward; on Gauss
         call rzero     (ud,nrstd)
      
c        call outmat(c(1,e,1),md,md,'fine u',e)
c        call outmat(c(1,e,2),md,md,'fine v',e)
c        call outmat(ju      ,md,md,'fine T',e)
      
         do j=1,ndim
            do i=1,nrstd
               tu(i)=c(i,e,j)*ju(i)   ! C_j*T
            enddo
            call grad_rst(ur,us,ut,tu,md,if3d)  ! Already on fine (Gauss
      
            j0 = j+0
            j2 = j+2
            do i=1,nrstd   ! rx has mass matrix and Jacobian on fine mes
               ud(i)=ud(i)+rx(i,j0,e)*ur(i)+rx(i,j2,e)*us(i)
            enddo
         enddo
      
         call intp_rstd(du(1,e),ud,mx,md,if3d,1) ! 1 --> transpose
      
      enddo
      
c     call exitti('convop_cons_2d$',istep)
      
      return
      end
c-----------------------------------------------------------------------
      

      
c-----------------------------------------------------------------------
# 1209
      subroutine char_conv_acc(p0,u,ulag,msk,c,cs,gsl)
c     
c     
c     Convect over last NBD steps using characteristics scheme
c     
c     NOTE:  Here, we assume that ulag is stored by time-slice first,
c            then by field number (this is opposite to prior Nek5000)
c     
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 1219 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1219 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOTAL" 1

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DXYZ" 1
C     
C     Elemental derivative operators
C     
# 4
      COMMON /DXYZ/ DXM1(LX1,LX1),  DXM12(LX2,LX1)
     $             ,DYM1(LY1,LY1),  DYM12(LY2,LY1)
     $             ,DZM1(LZ1,LZ1),  DZM12(LZ2,LZ1)
     $             ,DXTM1(LX1,LX1), DXTM12(LX1,LX2)
     $             ,DYTM1(LY1,LY1), DYTM12(LY1,LY2)
     $             ,DZTM1(LZ1,LZ1), DZTM12(LZ1,LZ2)
     $             ,DXM3(LX3,LX3),  DXTM3(LX3,LX3)
     $             ,DYM3(LY3,LY3),  DYTM3(LY3,LY3)
     $             ,DZM3(LZ3,LZ3),  DZTM3(LZ3,LZ3)
     $             ,DCM1(LY1,LY1),  DCTM1(LY1,LY1)
     $             ,DCM3(LY3,LY3),  DCTM3(LY3,LY3)
     $             ,DCM12(LY2,LY1), DCTM12(LY1,LY2)
     $             ,DAM1(LY1,LY1),  DATM1(LY1,LY1)
     $             ,DAM12(LY2,LY1), DATM12(LY1,LY2)
     $             ,DAM3(LY3,LY3),  DATM3(LY3,LY3)
      
# 2 "TOTAL" 2
# 2 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DEALIAS" 1
# 1
      common /solnd/  vxd(lxd,lyd,lzd,lelv)
     $             ,  vyd(lxd,lyd,lzd,lelv)
     $             ,  vzd(lxd,lyd,lzd,lelv)
      common /interpd/ imd1(lx1,lxd),imd1t(lxd,lx1)
     $               , im1d(lxd,lx1),im1dt(lx1,lxd)
     $               , pmd1(lx1,lxd),pmd1t(lxd,lx1)
c     common /dedim/ nxd,nyd,nzd
      real imd1,imd1t,im1d,im1dt
c     
# 3 "TOTAL" 2
# 3 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 4 "TOTAL" 2
# 4 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GEOM" 1
C     
C     Geometry arrays
C     
# 4
      COMMON /GXYZ/
     $              XM1   (LX1,LY1,LZ1,LELT)
     $             ,YM1   (LX1,LY1,LZ1,LELT)
     $             ,ZM1   (LX1,LY1,LZ1,LELT)
     $             ,XM2   (LX2,LY2,LZ2,LELV)
     $             ,YM2   (LX2,LY2,LZ2,LELV)
     $             ,ZM2   (LX2,LY2,LZ2,LELV)
C     
      COMMON /GISO1/
     $              RXM1  (LX1,LY1,LZ1,LELT)
     $             ,SXM1  (LX1,LY1,LZ1,LELT)
     $             ,TXM1  (LX1,LY1,LZ1,LELT)
     $             ,RYM1  (LX1,LY1,LZ1,LELT)
     $             ,SYM1  (LX1,LY1,LZ1,LELT)
     $             ,TYM1  (LX1,LY1,LZ1,LELT)
     $             ,RZM1  (LX1,LY1,LZ1,LELT)
     $             ,SZM1  (LX1,LY1,LZ1,LELT)
     $             ,TZM1  (LX1,LY1,LZ1,LELT)
     $             ,JACM1 (LX1,LY1,LZ1,LELT)
     $             ,jacmi (lx1*ly1*lz1,lelt)
      real          jacm1,jacmi
C     
      COMMON /GISO2/
     $              RXM2  (LX2,LY2,LZ2,LELV)
     $             ,SXM2  (LX2,LY2,LZ2,LELV)
     $             ,TXM2  (LX2,LY2,LZ2,LELV)
     $             ,RYM2  (LX2,LY2,LZ2,LELV)
     $             ,SYM2  (LX2,LY2,LZ2,LELV)
     $             ,TYM2  (LX2,LY2,LZ2,LELV)
     $             ,RZM2  (LX2,LY2,LZ2,LELV)
     $             ,SZM2  (LX2,LY2,LZ2,LELV)
     $             ,TZM2  (LX2,LY2,LZ2,LELV)
     $             ,JACM2 (LX2,LY2,LZ2,LELV)
      REAL          JACM2
c     
      common /gisod/ rx(lxd*lyd*lzd,ldim*ldim,lelv)
c     
      COMMON /GMFACT/
     $              G1M1  (LX1,LY1,LZ1,LELT)
     $             ,G2M1  (LX1,LY1,LZ1,LELT)
     $             ,G3M1  (LX1,LY1,LZ1,LELT)
     $             ,G4M1  (LX1,LY1,LZ1,LELT)
     $             ,G5M1  (LX1,LY1,LZ1,LELT)
     $             ,G6M1  (LX1,LY1,LZ1,LELT)
C     
      COMMON /GSURF/
     $              UNX   (LX1,LZ1,6,LELT)
     $             ,UNY   (LX1,LZ1,6,LELT)
     $             ,UNZ   (LX1,LZ1,6,LELT)
     $             ,T1X   (LX1,LZ1,6,LELT)
     $             ,T1Y   (LX1,LZ1,6,LELT)
     $             ,T1Z   (LX1,LZ1,6,LELT)
     $             ,T2X   (LX1,LZ1,6,LELT)
     $             ,T2Y   (LX1,LZ1,6,LELT)
     $             ,T2Z   (LX1,LZ1,6,LELT)
     $             ,AREA  (LX1,LZ1,6,LELT)
     $             ,DLAM
      COMMON /GVOLM/
     $              VNX   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNY   (LX1M,LY1M,LZ1M,LELT)
     $             ,VNZ   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V1Z   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2X   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Y   (LX1M,LY1M,LZ1M,LELT)
     $             ,V2Z   (LX1M,LY1M,LZ1M,LELT)
C     
      COMMON /GLOG/
     $        IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER(LELT),IFQINP(6,LELV),IFEPPM(6,LELV)
     $       ,IFLMSF(0:1),IFLMSE(0:1),IFLMSC(0:1)
     $       ,IFMSFC(6,LELT,0:1)
     $       ,IFMSEG(12,LELT,0:1)
     $       ,IFMSCR(8,LELT,0:1)
     $       ,IFNSKP(8,LELT)
     $       ,IFBCOR
      LOGICAL IFGEOM,IFGMSH3,IFVCOR,IFSURT,IFMELT,IFWCNO
     $       ,IFRZER,IFQINP,IFEPPM
     $       ,IFLMSF,IFLMSE,IFLMSC,IFMSFC
     $       ,IFMSEG,IFMSCR,IFNSKP
     $       ,IFBCOR
C     
# 5 "TOTAL" 2
# 5 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 6 "TOTAL" 2
# 6 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/IXYZ" 1
C     
C     Interpolation operators
C     
# 4
      COMMON /IXYZ/ IXM12 (LX2,LX1),  IXM21 (LX1,LX2)
     $             ,IYM12 (LY2,LY1),  IYM21 (LY1,LY2)
     $             ,IZM12 (LZ2,LZ1),  IZM21 (LZ1,LZ2)
     $             ,IXTM12(LX1,LX2),  IXTM21(LX2,LX1)
     $             ,IYTM12(LY1,LY2),  IYTM21(LY2,LY1)
     $             ,IZTM12(LZ1,LZ2),  IZTM21(LZ2,LZ1)
     $             ,IXM13 (LX3,LX1),  IXM31 (LX1,LX3)
     $             ,IYM13 (LY3,LY1),  IYM31 (LY1,LY3)
     $             ,IZM13 (LZ3,LZ1),  IZM31 (LZ1,LZ3)
     $             ,IXTM13(LX1,LX3),  IXTM31(LX3,LX1)
     $             ,IYTM13(LY1,LY3),  IYTM31(LY3,LY1)
     $             ,IZTM13(LZ1,LZ3),  IZTM31(LZ3,LZ1)
      COMMON /IXYZA/
     $              IAM12 (LY2,LY1),  IAM21 (LY1,LY2)
     $             ,IATM12(LY1,LY2),  IATM21(LY2,LY1)
     $             ,IAM13 (LY3,LY1),  IAM31 (LY1,LY3)
     $             ,IATM13(LY1,LY3),  IATM31(LY3,LY1)
     $             ,ICM12 (LY2,LY1),  ICM21 (LY1,LY2)
     $             ,ICTM12(LY1,LY2),  ICTM21(LY2,LY1)
     $             ,ICM13 (LY3,LY1),  ICM31 (LY1,LY3)
     $             ,ICTM13(LY1,LY3),  ICTM31(LY3,LY1)
     $             ,IAJL1 (LY1,LY1),  IATJL1(LY1,LY1)
     $             ,IAJL2 (LY2,LY2),  IATJL2(LY2,LY2)
     $             ,IALJ3 (LY3,LY3),  IATLJ3(LY3,LY3)
     $             ,IALJ1 (LY1,LY1),  IATLJ1(LY1,LY1)
C     
      REAL   IXM12,IYM12,IZM12,IXM21,IYM21,IZM21
      REAL   IXTM12,IYTM12,IZTM12,IXTM21,IYTM21,IZTM21
      REAL   IXM13,IYM13,IZM13,IXM31,IYM31,IZM31
      REAL   IXTM13,IYTM13,IZTM13,IXTM31,IYTM31,IZTM31
      REAL   IAM12,IAM21,IATM12,IATM21,IAM13,IAM31,IATM13,IATM31
      REAL   ICM12,ICM21,ICTM12,ICTM21,ICM13,ICM31,ICTM13,ICTM31
      REAL   IAJL1,IATJL1,IAJL2,IATJL2,IALJ3,IATLJ3,IALJ1,IATLJ1
# 7 "TOTAL" 2
# 7 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 8 "TOTAL" 2
# 8 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MVGEOM" 1
C     
C     Moving mesh data
C     
# 4
      COMMON /WSOL/   WX    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WY    (LX1M,LY1M,LZ1M,LELT)
     $            ,   WZ    (LX1M,LY1M,LZ1M,LELT)
      COMMON /WLAG/   WXLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WYLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
     $            ,   WZLAG (LX1M,LY1M,LZ1M,LELT,LORDER-1)
      COMMON /WMSU/   W1MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W2MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   W3MASK(LX1M,LY1M,LZ1M,LELT)
     $            ,   WMULT (LX1M,LY1M,LZ1M,LELT)
      COMMON /EIGVEC/ EV1   (LX1M,LY1M,LZ1M,LELV)
     $              , EV2   (LX1M,LY1M,LZ1M,LELV)
     $              , EV3   (LX1M,LY1M,LZ1M,LELV)
# 9 "TOTAL" 2
# 9 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 10 "TOTAL" 2
# 10 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 11 "TOTAL" 2
# 11 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/STEADY" 1
# 1
      COMMON /SSPAR1/ TAUSS(LDIMT1) , TXNEXT(LDIMT1)
      COMMON /SSPAR2/ NSSKIP
      COMMON /SSPAR3/ IFSKIP, IFMODP, IFSSVT, IFSTST(LDIMT1)
     $              ,                 IFEXVT, IFEXTR(LDIMT1)
      LOGICAL         IFSKIP, IFMODP, IFSSVT, IFSTST
     $              ,                 IFEXVT, IFEXTR
      COMMON /SSNORM/ DVNNH1, DVNNSM, DVNNL2, DVNNL8
     $              , DVDFH1, DVDFSM, DVDFL2, DVDFL8
     $              , DVPRH1, DVPRSM, DVPRL2, DVPRL8
# 12 "TOTAL" 2
# 12 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TOPOL" 1
C     
C     Arrays for direct stiffness summation
C     
# 4
      COMMON /CFACES/ NOMLIS(2,3),NMLINV(6),GROUP(6),SKPDAT(6,6)
     $               ,EFACE(6),EFACE1(6)
c     
      INTEGER NOMLIS,NMLINV,GROUP,SKPDAT,EFACE,EFACE1
      COMMON /CEDGES/ ESKIP(-12:12,3),NEDG(3)    ,NCMP
     $               ,IXCN(8),NOFFST(3,0:LDIMT1)
     $               ,MAXMLT,NSPMAX(0:LDIMT1)
     $               ,NGSPCN(0:LDIMT1),NGSPED(3,0:LDIMT1)
     $               ,NUMSCN(LELT,0:LDIMT1),NUMSED(LELT,0:LDIMT1)
     $               ,GCNNUM( 8,LELT,0:LDIMT1),LCNNUM( 8,LELT,0:LDIMT1)
     $               ,GEDNUM(12,LELT,0:LDIMT1),LEDNUM(12,LELT,0:LDIMT1)
     $               ,GEDTYP(12,LELT,0:LDIMT1)
     $               ,NGCOMM(2,0:LDIMT1)
      INTEGER ESKIP,NEDG,IXCN,MAXMLT,NSPMAX,NOFFST,NGSPCN,NGSPED
     $       ,NUMSCN,NUMSED,GCNNUM,LCNNUM,GEDNUM,LEDNUM,GEDTYP
     $       ,NGCOMM
      COMMON /EDGES/ IEDGE(20),IEDGEF(2,4,6,0:1)
     $              ,ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
     $              ,INDX(8),INVEDG(27)
      
# 13 "TOTAL" 2
# 13 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 14 "TOTAL" 2
# 14 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TURBO" 1
C     
C       Common block for turbulence model
C     
# 4
      COMMON /TURBR/ VTURB (LX1M,LY1M,LZ1M,LELV)
     $             , TURBL (LX1M,LY1M,LZ1M,LELV)
     $             , UWALL (LX1M,LZ1M,6,LELV)
     $             , ZWALL (LX1M,LZ1M,6,LELV)
     $             , TWX   (LX1M,LZ1M,6,LELV)
     $             , TWY   (LX1M,LZ1M,6,LELV)
     $             , TWZ   (LX1M,LZ1M,6,LELV)
      COMMON /TURBC/ CMU,CMT,SGK,SGE,CE1,CE2,VKC,BTA,SGT
     $             , BETA1,BETA2
     $             , CMI,SKI,SEI,VKI,BTI,STI
     $             , ZPLDAT,ZPUDAT,ZPVDAT,TLMAX,TLIMUL
      COMMON /TURBI/ IFLDK,IFLDTK,IFLDE,IFLDTE
      COMMON /TURBL/ IFSWALL,IFTWSH(6,LELV),IFCWUZ
C     
      
      LOGICAL IFSWALL,IFTWSH,IFCWUZ
# 15 "TOTAL" 2
# 15 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ESOLV" 1
# 1
      common /econst/ iesolv
      common /efastm/ ifalgn(lelv), ifrsxy(lelv)
      logical         ifalgn, ifrsxy
      common /eouter/ volel(lelv) 
# 16 "TOTAL" 2
# 16 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZ" 1
C     
C     Gauss-Labotto and Gauss points
C     
# 4
      COMMON /GAUSS/ ZGM1(LX1,3), ZGM2(LX2,3), ZGM3(LX3,3)
     $              ,ZAM1(LX1)  , ZAM2(LX2)  , ZAM3(LX3)
C     
C    Weights
C     
      COMMON /WXYZ/ WXM1(LX1), WYM1(LY1), WZM1(LZ1), W3M1(LX1,LY1,LZ1)
     $             ,WXM2(LX2), WYM2(LY2), WZM2(LZ2), W3M2(LX2,LY2,LZ2)
     $             ,WXM3(LX3), WYM3(LY3), WZM3(LZ3), W3M3(LX3,LY3,LZ3)
     $             ,WAM1(LY1), WAM2(LY2), WAM3(LY3)
     $             ,W2AM1(LX1,LY1), W2CM1(LX1,LY1)
     $             ,W2AM2(LX2,LY2), W2CM2(LX2,LY2)
     $             ,W2AM3(LX3,LY3), W2CM3(LX3,LY3)
# 17 "TOTAL" 2
# 17 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/WZF" 1
c     
c Points (z) and weights (w) on velocity, pressure
c     
c     zgl -- velocity points on Gauss-Lobatto points i = 1,...nx
c     zgp -- pressure points on Gauss         points i = 1,...nxp (nxp =
c     
c     parameter (lxm = lx1)
# 8
      parameter (lxq = lx2)
c     
      common /wz1/   zgl(lx1),wgl(lx1)
     $           ,   zgp(lx1),wgp(lxq)
c     
c     Tensor- (outer-) product of 1D weights   (for volumetric integrati
c     
      common /wz2/  wgl1(lx1*lx1),wgl2(lxq*lxq)
     $           ,  wgli(lx1*lx1)
c     
c     
c    Frequently used derivative matrices:
c     
c    D1, D1t   ---  differentiate on mesh 1 (velocity mesh)
c    D2, D2t   ---  differentiate on mesh 2 (pressure mesh)
c     
c    DXd,DXdt  ---  differentiate from velocity mesh ONTO dealiased mesh
c                   (currently the same as D1 and D1t...)
c     
c     
      common /deriv/  d1    (lx1*lx1) , d1t    (lx1*lx1)
     $             ,  d2    (lx1*lx1) , b2p    (lx1*lx1)
     $             ,  B1iA1 (lx1*lx1) , B1iA1t (lx1*lx1)
     $             ,  da    (lx1*lx1) , dat    (lx1*lx1)
     $             ,  iggl  (lx1*lxq) , igglt  (lx1*lxq)
     $             ,  dglg  (lx1*lxq) , dglgt  (lx1*lxq)
     $             ,  wglg  (lx1*lxq) , wglgt  (lx1*lxq)
      real ixd,ixdt,iggl,igglt
c     
# 18 "TOTAL" 2
# 18 "TOTAL"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/DSSUM" 1
# 1
      parameter (lds=lx1*ly1*lz1*lelt)
      common /newdss/ ids_lgl1(-1:lds),ids_lgl2(-1:lds),ids_ptr(lds)
# 1220 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1220 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      real    p0(1),u(1),ulag(1),msk(1),c(1),cs(0:1)
      integer gsl
      
      common /scrns/ ct  (lxd*lyd*lzd*lelv*ldim)
      
      common /scrvh/ bmsk(lx1*ly1*lz1*lelv)
     $             , u1  (lx1*ly1*lz1*lelv)
      
      common /scrmg/ r1  (lx1*ly1*lz1*lelv)
     $             , r2  (lx1*ly1*lz1*lelv)
     $             , r3  (lx1*ly1*lz1*lelv)
     $             , r4  (lx1*ly1*lz1*lelv)
      
c     
      nelc = nelv            ! number of elements in convecting field
      if (ifield.eq.ifldmhd) nelc = nelfld(ifield)
      
      nc  = cs(0)            ! number of stored convecting fields
      
      ln  = lx1*ly1*lz1*lelt
      n   = nx1*ny1*nz1*nelfld(ifield)
      m   = nxd*nyd*nzd*nelc*ndim
      
      if (ifield.eq.ifldmhd) then
         call col3(bmsk,bintm1,msk,n)
      elseif (ifield.eq.1) then
         call col3(bmsk,binvm1,msk,n)
      else ! if (ifield.eq.2) then
         call col3(bmsk,bintm1,msk,n)
      endif
      
      call char_conv1
     $   (p0,bmsk,u,n,ulag,ln,gsl,c,m,cs(1),nc,ct,u1,r1,r2,r3,r4)
      
      return
      end
      
      
c-----------------------------------------------------------------------
      subroutine advchar_acc
c     
c     Compute convective contribution using 
c     operator-integrator-factor method (characteristics).
c     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/pipe_2013_11_27/SIZE" 1
C     Dimension file to be included     
C     
C     HCUBE array dimensions
C     
# 5
      parameter (ldim=3)
      parameter (lx1=10,ly1=lx1,lz1=lx1,lelt=410,lelv=lelt)
      parameter (lxd=10,lyd=lxd,lzd=lxd)
      parameter (lelx=1,lely=1,lelz=1)
      
      parameter (lzl=3 + 2*(ldim-3))
      
c      parameter (lx2=lx1-2)
c      parameter (ly2=ly1-2)
c      parameter (lz2=lz1-2)
      parameter (lx2=lx1)
      parameter (ly2=ly1)
      parameter (lz2=lz1)
      
      parameter (lx3=lx1)
      parameter (ly3=ly1)
      parameter (lz3=lz1)
      
      parameter (lp = 4)
      parameter (lelg = 410)
c     
c     parameter (lpelv=lelv,lpelt=lelt,lpert=3)  ! perturbation
c     parameter (lpx1=lx1,lpy1=ly1,lpz1=lz1)     ! array sizes
c     parameter (lpx2=lx2,lpy2=ly2,lpz2=lz2)
c     
      parameter (lpelv=1,lpelt=1,lpert=1)        ! perturbation
      parameter (lpx1=1,lpy1=1,lpz1=1)           ! array sizes
      parameter (lpx2=1,lpy2=1,lpz2=1)
c     
c     parameter (lbelv=lelv,lbelt=lelt)          ! MHD
c     parameter (lbx1=lx1,lby1=ly1,lbz1=lz1)     ! array sizes
c     parameter (lbx2=lx2,lby2=ly2,lbz2=lz2)
c     
      parameter (lbelv=1,lbelt=1)                ! MHD
      parameter (lbx1=1,lby1=1,lbz1=1)           ! array sizes
      parameter (lbx2=1,lby2=1,lbz2=1)
      
C     LX1M=LX1 when there are moving meshes; =1 otherwise
      parameter (lx1m=1,ly1m=1,lz1m=1)
      parameter (ldimt= 2)                       ! 3 passive scalars + T
      parameter (ldimt1=ldimt+1)
      parameter (ldimt3=ldimt+3)
c     
c     Note:  In the new code, LELGEC should be about sqrt(LELG)
c     
      PARAMETER (LELGEC = 1)
      PARAMETER (LXYZ2  = 1)
      PARAMETER (LXZ21  = 1)
      
      PARAMETER (LMAXV=LX1*LY1*LZ1*LELV)
      PARAMETER (LMAXT=LX1*LY1*LZ1*LELT)
      PARAMETER (LMAXP=LX2*LY2*LZ2*LELV)
      PARAMETER (LXZ=LX1*LZ1)
      PARAMETER (LORDER=3)
      PARAMETER (MAXOBJ=4,MAXMBR=LELT*6)
      PARAMETER (lhis=100)         ! # of pts a proc reads from hpts.in
                                   ! Note: lhis*np > npoints in hpts.in
C     
C     Common Block Dimensions
C     
      PARAMETER (LCTMP0 =2*LX1*LY1*LZ1*LELT)
      PARAMETER (LCTMP1 =4*LX1*LY1*LZ1*LELT)
C     
C     The parameter LVEC controls whether an additional 42 field arrays
C     are required for Steady State Solutions.  If you are not using
C     Steady State, it is recommended that LVEC=1.
C     
      PARAMETER (LVEC=1)
C     
C     Uzawa projection array dimensions
C     
      parameter (mxprev = 20)
      parameter (lgmres = 30)
C     
C     Split projection array dimensions
C     
      parameter(lmvec = 1)
      parameter(lsvec = 1)
      parameter(lstore=lmvec*lsvec)
c     
c     NONCONFORMING STUFF
c     
      parameter (maxmor = lelt)
C     
C     Array dimensions
C     
      COMMON/DIMN/NELV,NELT,NX1,NY1,NZ1,NX2,NY2,NZ2
     $,NX3,NY3,NZ3,NDIM,NFIELD,NPERT,NID
     $,NXD,NYD,NZD
      
c automatically added by makenek
      parameter(lxo   = lx1) ! max output grid size (lxo>=lx1)
      
c automatically added by makenek
      parameter(lpart = 1  ) ! max number of particles
      
c automatically added by makenek
      integer ax1,ay1,az1,ax2,ay2,az2
      parameter (ax1=lx1,ay1=ly1,az1=lz1,ax2=lx2,ay2=ly2,az2=lz2) ! runn
      
c automatically added by makenek
      parameter (lxs=1,lys=lxs,lzs=(lxs-1)*(ldim-2)+1) !New Pressure Pre
      
c automatically added by makenek
      parameter (lfdm=0)  ! == 1 for fast diagonalization method
# 1265 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1265 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 1266 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1266 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/INPUT" 1
C     
C     Input parameters from preprocessors.
C     
C     Note that in parallel implementations, we distinguish between
C     distributed data (LELT) and uniformly distributed data.
C     
C     Input common block structure:
C     
C     INPUT1:  REAL            INPUT5: REAL      with LELT entries
C     INPUT2:  INTEGER         INPUT6: INTEGER   with LELT entries
C     INPUT3:  LOGICAL         INPUT7: LOGICAL   with LELT entries
C     INPUT4:  CHARACTER       INPUT8: CHARACTER with LELT entries
C     
# 14
      COMMON /INPUT1/ PARAM(200)
     $               ,RSTIM,VNEKTON
     $               ,CPFLD(LDIMT1,3)
     $               ,CPGRP(-5:10,LDIMT1,3)
     $               ,QINTEG(LDIMT3,MAXOBJ)
C     
      COMMON /INPUT2/ MATYPE(-5:10,LDIMT1)
     $               ,NKTONV,NHIS,LOCHIS(4,lhis)
     $               ,IPSCAL,NPSCAL,IPSCO, ifldmhd
     $               ,IRSTV,IRSTT,IRSTIM,NMEMBER(MAXOBJ),NOBJ
     $               ,NGEOM
C     
      COMMON /INPUT3/ IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC(LDIMT1),IFTMSH(0:LDIMT1)
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL(LDIMT1)
     $               ,IFVARP(LDIMT1),IFPSCO(LDIMT1),IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO(LDIMT1),IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFSCHCLOB,IFUSERVP,
     $               IFCYCLIC,IFMOAB,IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,
     $               IFXYO_,ifaziv,IFNEKNEK
      LOGICAL         IF3D
     $               ,IFFLOW,IFHEAT,IFTRAN,IFAXIS,IFSTRS,IFSPLIT
     $               ,IFMGRID,IFADVC,IFTMSH
     $               ,IFMVBD,IFNATC,IFCHAR,IFNONL
     $               ,IFVARP,IFPSCO,IFVPS
     $               ,IFMODEL,IFKEPS
     $               ,IFINTQ,IFCONS
     $               ,IFXYO,IFPO,IFVO,IFTO,IFTGO,IFPSO,IFFMTIN
     $               ,IFBO
     $               ,IFANLS,IFANL2,IFMHD,IFESSR,IFPERT,IFBASE
     $               ,IFCVODE,IFLOMACH,IFEXPLVIS,IFUSERVP,IFCYCLIC
     $               ,IFSCHCLOB
     $               ,IFMOAB, IFCOUP, IFVCOUP, IFUSERMV,IFREGUO,IFXYO_
     $               ,ifaziv,IFNEKNEK
      LOGICAL         IFNAV
      EQUIVALENCE    (IFNAV, IFADVC(1))
C     
      COMMON /INPUT4/ HCODE(11,lhis),OCODE(8),RSTV,RSTT,DRIVC(5)
     $               ,INITC(15),TEXTSW(100,2)
      CHARACTER*1     HCODE
      CHARACTER*2     OCODE
      CHARACTER*10    DRIVC
      CHARACTER*14    RSTV,RSTT
      CHARACTER*40    TEXTSW,TURBMOD
      CHARACTER*132    INITC
      EQUIVALENCE    (TURBMOD,TEXTSW(1,1))
C     
      COMMON /CFILES/ REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      CHARACTER*132   REAFLE,FLDFLE,DMPFLE,HISFLE,SCHFLE,OREFLE,NREFLE
      COMMON /CFILE2/ SESSION,PATH,RE2FLE, H5MFLE
      CHARACTER*132   SESSION,PATH,RE2FLE,H5MFLE
C     
C proportional to LELT
C     
      COMMON /INPUT5/ XC(8,LELT),YC(8,LELT),ZC(8,LELT)
     $               ,BC(5,6,LELT,0:LDIMT1)
     $               ,CURVE(6,12,LELT)
     $               ,CERROR(LELT)
C     
      COMMON /INPUT6/ IGROUP(LELT),OBJECT(MAXOBJ,MAXMBR,2)
      INTEGER OBJECT
C     
      COMMON /INPUT8/ CBC(6,LELT,0:LDIMT1),CCURVE(12,LELT)
     $              , CDOF(6,LELT), solver_type
      CHARACTER*1     CCURVE,CDOF
      CHARACTER*3     CBC, solver_type
      COMMON /INPUT9/ IEACT(LELT),NEACT
C     
C material set ids, BC set ids, materials (f=fluid, s=solid), bc types
      PARAMETER (NUMSTS=50)
      COMMON /INPUTMI/ NUMFLU, NUMOTH, NUMBCS 
     $	             , MATINDX(NUMSTS),MATIDS(NUMSTS),IMATIE(LELT)
     $ 	             , IBCSTS (NUMSTS) 
      COMMON /INPUTMR/ BCF    (NUMSTS)
      COMMON /INPUTMC/ BCTYPS (NUMSTS)
      CHARACTER*3 BCTYPS
      integer     bcf
      
# 1267 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1267 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SOLN" 1
# 1
      parameter (lvt1  = lx1*ly1*lz1*lelv)
      parameter (lvt2  = lx2*ly2*lz2*lelv)
      parameter (lbt1  = lbx1*lby1*lbz1*lbelv)
      parameter (lbt2  = lbx2*lby2*lbz2*lbelv)
      
      parameter (lptmsk = lvt1*(5+2*ldimt) + 4*lbt1)
      parameter (lptsol
     $         = lvt1*(12 + 4*ldimt + 2*ldimt1 + (3+ldimt)*(lorder-1))
     $         + lvt2*(lorder-1)
     $         + lbt1*(12 + 3*(lorder-1))
     $         + lbt2*(lorder-1) )
c    $         + lptmsk )
      
      parameter (lorder2 = max(1,lorder-2) )
C     
C     Solution and data
C     
      COMMON /BQCB/    BQ     (LX1,LY1,LZ1,LELT,LDIMT)
      
      COMMON /VPTSOL/  
c     Can be used for post-processing runs (SIZE .gt. 10+3*LDIMT flds)
     $                 VXLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VYLAG  (LX1,LY1,LZ1,LELV,2)
     $               , VZLAG  (LX1,LY1,LZ1,LELV,2)
     $               , TLAG   (LX1,LY1,LZ1,LELT,LORDER-1,LDIMT)
     $               , VGRADT1(LX1,LY1,LZ1,LELT,LDIMT)
     $               , VGRADT2(LX1,LY1,LZ1,LELT,LDIMT)
     $               , ABX1   (LX1,LY1,LZ1,LELV)
     $               , ABY1   (LX1,LY1,LZ1,LELV)
     $               , ABZ1   (LX1,LY1,LZ1,LELV)
     $               , ABX2   (LX1,LY1,LZ1,LELV)
     $               , ABY2   (LX1,LY1,LZ1,LELV)
     $               , ABZ2   (LX1,LY1,LZ1,LELV)
     $               , VDIFF_E(LX1,LY1,LZ1,LELT)
c     Solution data
     $               , VX     (LX1,LY1,LZ1,LELV)
     $               , VY     (LX1,LY1,LZ1,LELV)
     $               , VZ     (LX1,LY1,LZ1,LELV)
     $               , T      (LX1,LY1,LZ1,LELT,LDIMT)
     $               , VTRANS (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , VDIFF  (LX1,LY1,LZ1,LELT,LDIMT1)
     $               , BFX    (LX1,LY1,LZ1,LELV)
     $               , BFY    (LX1,LY1,LZ1,LELV)
     $               , BFZ    (LX1,LY1,LZ1,LELV)
     $               , cflf   (lx1,ly1,lz1,lelv)
     $               , c_vx   (lxd*lyd*lzd*lelv*ldim,lorder+1) ! charact
c     Solution data for magnetic field
     $               , BX     (LBX1,LBY1,LBZ1,LBELV)
     $               , BY     (LBX1,LBY1,LBZ1,LBELV)
     $               , BZ     (LBX1,LBY1,LBZ1,LBELV)
     $               , PM     (LBX2,LBY2,LBZ2,LBELV)
     $               , BMX    (LBX1,LBY1,LBZ1,LBELV)  ! Magnetic field R
     $               , BMY    (LBX1,LBY1,LBZ1,LBELV)
     $               , BMZ    (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX1   (LBX1,LBY1,LBZ1,LBELV) ! Extrapolation ter
     $               , BBY1   (LBX1,LBY1,LBZ1,LBELV) ! magnetic field rh
     $               , BBZ1   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBX2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBY2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BBZ2   (LBX1,LBY1,LBZ1,LBELV)
     $               , BXLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BYLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , BZLAG  (LBX1*LBY1*LBZ1*LBELV,LORDER-1)
     $               , PMLAG  (LBX2*LBY2*LBZ2*LBELV,LORDER2)
      
      common /expvis/  nu_star
      real             nu_star
      
      COMMON /CBM2/  
     $                 PR     (LX2,LY2,LZ2,LELV)
     $               , PRLAG  (LX2,LY2,LZ2,LELV,LORDER2)
c     
      COMMON /DIVERG/  QTL    (LX2,LY2,LZ2,LELT)
     $               , USRDIV (LX2,LY2,LZ2,LELT)
      
      COMMON /VPTMSK/  V1MASK (LX1,LY1,LZ1,LELV)
     $               , V2MASK (LX1,LY1,LZ1,LELV)
     $               , V3MASK (LX1,LY1,LZ1,LELV)
     $               , PMASK  (LX1,LY1,LZ1,LELV)
     $               , TMASK  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , OMASK  (LX1,LY1,LZ1,LELT)
     $               , VMULT  (LX1,LY1,LZ1,LELV)
     $               , TMULT  (LX1,LY1,LZ1,LELT,LDIMT)
     $               , B1MASK (LBX1,LBY1,LBZ1,LBELV)  ! masks for mag. f
     $               , B2MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , B3MASK (LBX1,LBY1,LBZ1,LBELV)
     $               , BPMASK (LBX1,LBY1,LBZ1,LBELV)  ! magnetic pressur
C     
C     Solution and data for perturbation fields
C     
      COMMON /PVPTSL/ VXP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VYP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VZP    (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , PRP    (LPX2*LPY2*LPZ2*LPELV,lpert)
     $              , TP     (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BQP    (LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              , BFXP   (LPX1*LPY1*LPZ1*LPELV,lpert)  ! perturbatio
     $              , BFYP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , BFZP   (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , VXLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VYLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , VZLAGP (LPX1*LPY1*LPZ1*LPELV,LORDER-1,lpert)
     $              , PRLAGP (LPX2*LPY2*LPZ2*LPELV,LORDER2,lpert)
     $              , TLAGP  (LPX1*LPY1*LPZ1*LPELT,LDIMT,LORDER-1,lpert)
     $              , EXX1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! Extrapolatio
     $              , EXY1P  (LPX1*LPY1*LPZ1*LPELV,lpert) ! perturbation
     $              , EXZ1P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXX2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXY2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              , EXZ2P  (LPX1*LPY1*LPZ1*LPELV,lpert)
     $              ,VGRADT1P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
     $              ,VGRADT2P(LPX1*LPY1*LPZ1*LPELT,LDIMT,lpert)
c     
      common /ppointr/ jp
# 1268 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1268 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/TSTEP" 1
# 1
      COMMON /TSTEP1/ TIME,TIMEF,FINTIM,TIMEIO
     $               ,DT,DTLAG(10),DTINIT,DTINVM,COURNO,CTARG
     $               ,AB(10),BD(10),ABMSH(10)
     $               ,AVDIFF(LDIMT1),AVTRAN(LDIMT1),VOLFLD(0:LDIMT1)
     $               ,TOLREL,TOLABS,TOLHDF,TOLPDF,TOLEV,TOLNL,PRELAX
     $               ,TOLPS,TOLHS,TOLHR,TOLHV,TOLHT(LDIMT1),TOLHE
     $               ,VNRMH1,VNRMSM,VNRML2,VNRML8,VMEAN
     $               ,TNRMH1(LDIMT),TNRMSM(LDIMT),TNRML2(LDIMT)
     $               ,TNRML8(LDIMT),TMEAN(LDIMT)
C     
      COMMON /ISTEP2/ IFIELD,IMESH,ISTEP,NSTEPS,IOSTEP,LASTEP,IOCOMM
     $               ,INSTEP
     $               ,NAB,NBD,NBDINP,NTAUBD 
     $               ,NMXH,NMXP,NMXE,NMXNL,NINTER
     $               ,NELFLD(0:LDIMT1)
     $               ,nconv,nconv_max
C     
      COMMON /TSTEP3/ PI,BETAG,GTHETA
      COMMON /TSTEP4/ IFPRNT,if_full_pres
      LOGICAL IFPRNT,if_full_pres
c     
      COMMON /TSTEP5/ lyap(3,lpert)  !  lyapunov simulation history
      real lyap
# 1269 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1269 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/PARALLEL" 1
C     
C     Communication information
C     NOTE: NID is stored in 'SIZE' for greater accessibility
# 4
      COMMON /CUBE1/ NODE,PID,NP,NULLPID,NODE0
      INTEGER        NODE,PID,NP,NULLPID,NODE0
      
C     
C     Maximum number of elements (limited to 2**31/12, at least for now)
      PARAMETER(NELGT_MAX = 178956970)
C     
      COMMON /HCGLB/ NVTOT,NELG(0:LDIMT1)
     $              ,LGLEL(LELT)
     $              ,GLLEL(LELG)
     $              ,GLLNID(LELG)
     $              ,NELGV,NELGT
      
      INTEGER        GLLEL,GLLNID,LGLEL
      INTEGER*8      NVTOT
C     
      COMMON /DIAGL/  IFGPRNT
      LOGICAL IFGPRNT
      COMMON/PRECSN/ WDSIZE,ISIZE,LSIZE,CSIZE
      COMMON/PRECSL/ IFDBLAS
      INTEGER WDSIZE,ISIZE,LSIZE,CSIZE
      LOGICAL IFDBLAS
C     
C     crystal-router, gather-scatter, and xxt handles (xxt=csr grid solv
C     
      common /comm_handles/ cr_h, gsh, gsh_fld(0:ldimt3), xxth(ldimt3)
      integer               cr_h, gsh, gsh_fld          , xxth
# 1270 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1270 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/CTIMER" 1
C     
# 2
      COMMON /CTIMER/ tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
      COMMON /CTIME2/ tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep
     $               ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2
     $               ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn
     $               ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee
     $               ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc
     $               ,twal
C     
      REAL*8          tmxmf,tmxms,tdsum,taxhm,tcopy,tinvc,tinv3
      REAL*8          tsolv,tgsum,tdsnd,tdadd,tcdtp,tmltd,tprep
     $               ,tpres,thmhz,tgop ,tgop1,tdott,tbsol,tbso2
     $               ,tsett,tslvb,tusbc,tddsl,tcrsl,tdsmx,tdsmn
     $               ,tgsmn,tgsmx,teslv,tbbbb,tcccc,tdddd,teeee
     $               ,tvdss,tschw,tadvc,tspro,tgop_sync,tsyc
     $               ,twal
C     
      COMMON /ITIMER/ nmxmf,nmxms,ndsum,naxhm,ncopy,ninvc,ninv3
      COMMON /ITIME2/ nsolv,ngsum,ndsnd,ndadd,ncdtp,nmltd,nprep
     $               ,npres,nhmhz,ngop ,ngop1,ndott,nbsol,nbso2
     $               ,nsett,nslvb,nusbc,nddsl,ncrsl,ndsmx,ndsmn
     $               ,ngsmn,ngsmx,neslv,nbbbb,ncccc,ndddd,neeee
     $               ,nvdss,nadvc,nspro,ngop_sync,nsyc,nwal
C     
      COMMON /PTIMER/ pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
     $               ,psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep
     $               ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2
     $               ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn
     $               ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee
     $               ,pvdss,pspro,pgop_sync,psyc,pwal
C     
      REAL*8          pmxmf,pmxms,pdsum,paxhm,pcopy,pinvc,pinv3
      REAL*8          psolv,pgsum,pdsnd,pdadd,pcdtp,pmltd,pprep
     $               ,ppres,phmhz,pgop ,pgop1,pdott,pbsol,pbso2
     $               ,psett,pslvb,pusbc,pddsl,pcrsl,pdsmx,pdsmn
     $               ,pgsmn,pgsmx,peslv,pbbbb,pcccc,pdddd,peeee
     $               ,pvdss,pspro,pgop_sync,psyc,pwal
C     
      REAL*8 etime1,etime2,etime0,gtime1,tscrtch
      REAL*8 dnekclock,dnekclock_sync
C     
      COMMON /CTIME3/ etimes,ttotal,tttstp,etims0,ttime
      real*8          etimes,ttotal,tttstp,etims0,ttime
C     
      integer icalld
      save    icalld
      data    icalld /0/
      
      common /ctimel/ ifsync
      logical         ifsync
# 1272 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f" 2
# 1272 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/convect.f"
      
      common /cchar/ ct_vx(0:lorder) ! time for each slice in c_vx()
      
      common /scruz/ phx  (lx1*ly1*lz1*lelt)
     $ ,             phy  (lx1*ly1*lz1*lelt)
     $ ,             phz  (lx1*ly1*lz1*lelt)
     $ ,             hmsk (lx1*ly1*lz1*lelt)
      

# 1281
      if (icalld.eq.0) tadvc=0.0
      icalld=icalld+1
      nadvc=icalld
      etime1=dnekclock()

      
# 1287
      dti = 1./dt
      n   = nx1*ny1*nz1*nelv
      
      call set_char_mask(hmsk,vx,vy,vz) ! mask for hyperbolic system 
      
      call char_conv(phx,vx,vxlag,hmsk,c_vx,ct_vx,gsh_fld(1))
      call char_conv(phy,vy,vylag,hmsk,c_vx,ct_vx,gsh_fld(1))
      
      if (if3d) then
      
        call char_conv(phz,vz,vzlag,hmsk,c_vx,ct_vx,gsh_fld(1))
      
        do i=1,n
           h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
           bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
           bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
           bfz(i,1,1,1) = bfz(i,1,1,1)+phz(i)*h2i
        enddo
      
      else
      
        do i=1,n
           h2i = bm1(i,1,1,1)*vtrans(i,1,1,1,1)*dti
           bfx(i,1,1,1) = bfx(i,1,1,1)+phx(i)*h2i
           bfy(i,1,1,1) = bfy(i,1,1,1)+phy(i)*h2i
        enddo
      
      endif
      

# 1317
      tadvc=tadvc+(dnekclock()-etime1)

# 1319
      return
      end
      

      
