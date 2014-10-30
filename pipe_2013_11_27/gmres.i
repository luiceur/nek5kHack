# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
c-----------------------------------------------------------------------
      subroutine uzawa_gmres(res,h1,h2,h2inv,intype,iter)
      
c     Solve the pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)
      

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
# 11 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 11 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 12 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 12 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GMRES" 1
c     parameter (lgmres=50) !max number of gmres iterations between rest
# 2
      common /gmres/  x(lx2*ly2*lz2*lelv)
     $              , r(lx2*ly2*lz2*lelv), w(lx2*ly2*lz2*lelv)
     $              , h(lgmres,lgmres), gamma(lgmres+1)
     $              , c(lgmres), s(lgmres)
      common /gmre1/  v(lx2*ly2*lz2*lelv,lgmres)
      common /gmre2/  z(lx2*ly2*lz2*lelv,lgmres)
      real x,r,w,h,gamma,c,s,v,z
      
      common /spltprec/ ml(lx2*ly2*lz2*lelv), mu(lx2*ly2*lz2*lelv)
      real ml,mu
      
c     w is a work vector
c     c and s store the Givens rotations
c     V stores the orthogonal Krylov subspace basis
c          -1
c     Z = M   V
# 13 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 13 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real             res  (lx2*ly2*lz2*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)
      
      common /scrmg/    wp (lx2,ly2,lz2,lelv)
      
      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      common /cgmres1/ y(lgmres)
      
      real alpha, l, temp
      integer j,m
c     
      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac
c     
      real*8 etime1,dnekclock
c     
      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split0(ml,mu,bm2,bm2inv,nx2*ny2*nz2*nelv)
         norm_fac = 1./sqrt(volvm2)
      endif
c     
      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m = lgmres
c     
      call chktcg2(tolps,res,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
c     if (param(21).lt.0) tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c     
      ntot2  = nx2*ny2*nz2*nelv
c     
      iconv = 0
      call rzero(x,ntot2)
      
      do while(iconv.eq.0.and.iter.lt.100)
      
         if(iter.eq.0) then
                                                  !      -1
            call col3(r,ml,res,ntot2)             ! r = L  res
c           call copy(r,res,ntot2)
         else
            !update residual
            call copy(r,res,ntot2)                ! r = res
            call cdabdtp(w,x,h1,h2,h2inv,intype)  ! w = A x
            call add2s2(r,w,-1.,ntot2)            ! r = r - w
                                                  !      -1
            call col2(r,ml,ntot2)                 ! r = L   r
         endif
                                                  !            ______
         gamma(1) = sqrt(glsc2(r,r,ntot2))        ! gamma  = \/ (r,r) 
                                                  !      1
         if(iter.eq.0) then
            div0 = gamma(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif
      
         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,ntot2)         ! v  = r / gamma
                                                  !  1            1
         do j=1,m
            iter = iter+1
                                                  !       -1
            call col3(w,mu,v(1,j),ntot2)          ! w  = U   v
                                                  !           j
            
            etime2 = dnekclock()
            if(param(43).eq.1) then
               call uzprec(z(1,j),w,h1,h2,intype,wp)
            else                                  !       -1
               call hsmg_solve(z(1,j),w)          ! z  = M   w
c              call copy(z(1,j),w,ntot2)          ! z  = M   w
            endif     
            etime_p = etime_p + dnekclock()-etime2
      
            call cdabdtp(w,z(1,j),                ! w = A z
     $                   h1,h2,h2inv,intype)      !        j
      
                                                  !      -1
            call col2(w,ml,ntot2)                 ! w = L   w
      
c           !modified Gram-Schmidt
c           do i=1,j
c              h(i,j)=glsc2(w,v(1,i),ntot2)       ! h    = (w,v )
c                                                 !  i,j       i
c              call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
c           enddo                                 !          i,j  i
      
      
c           2-PASS GS, 1st pass:
      
            do i=1,j
               h(i,j)=vlsc2(w,v(1,i),ntot2)       ! h    = (w,v )
            enddo                                 !  i,j       i
      
            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs
      
            do i=1,j
               call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
            enddo                                 !          i,j  i
      
      
c           2-PASS GS, 2nd pass:
c     
c           do i=1,j
c              wk1(i)=vlsc2(w,v(1,i),ntot2)       ! h    = (w,v )
c           enddo                                 !  i,j       i
c                                                 !
c           call gop(wk1,wk2,'+  ',j)             ! sum over P procs
c     
c           do i=1,j
c              call add2s2(w,v(1,i),-wk1(i),ntot2)! w = w - h    v
c              h(i,j) = h(i,j) + wk1(i)           !          i,j  i
c           enddo
      
      
            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                  !            ______
            alpha = sqrt(glsc2(w,w,ntot2))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)
      
c            call outmat(h,m,j,' h    ',j)
            
            rnorm = abs(gamma(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nid.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence')
      

# 172
            if (rnorm .lt. tolpss) goto 900  !converged



# 176
            if (j.eq.m) goto 1000 !not converged, restart
      
            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,ntot2)   ! v    = w / alpha
                                                 !  j+1            
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,z(1,i),c(i),ntot2)     ! x = x + c  z
                                                 !          i  i
         enddo
c        if(iconv.eq.1) call dbg_write(x,nx2,ny2,nz2,nelv,'esol',3)
      enddo
 9000 continue
c     
      divex = rnorm
c     iter = iter - 1
c     
c     DIAGNOSTICS
c      call copy   (w,x,ntot2)
       call ortho  (w) ! Orthogonalize wrt null space, if present
c      call copy(r,res,ntot2) !r = res
c      call cdabdtp(r,w,h1,h2,h2inv,intype)  ! r = A w
c      do i=1,ntot2
c         r(i) = res(i) - r(i)               ! r = res - r
c      enddo
c      call uzawa_gmres_temp(r,bm2inv,ntot2)
c                                               !            ______
c      gamma(1) = sqrt(glsc2(r,r,ntot2)/volvm2) ! gamma  = \/ (r,r) 
c                                               !      1
c      print *, 'GMRES end resid:',gamma(1)
c     END DIAGNOSTICS
      call copy(res,x,ntot2)
      
      call ortho (res)  ! Orthogonalize wrt null space, if present
      
      etime1 = dnekclock()-etime1
      if (nid.eq.0) write(6,9999) istep,iter,divex,tolpss,div0,etime_p,
     &                            etime1
c     call flush_hack
 9999 format(i10,' U-PRES gmres:',I7,1p5e12.4)
      
      return
      end
      
c-----------------------------------------------------------------------
      
      subroutine uzawa_gmres_split0(l,u,b,binv,n)
      integer n
      real l(n),u(n),b(n),binv(n)
      integer i
      do i=1,n
         l(i)=sqrt(binv(i))
         u(i)=sqrt(b(i))
         if(abs(u(i)*l(i)-1.0).gt.1e-13) print *, i, u(i)*l(i)
      enddo
      return
      end
      
c-----------------------------------------------------------------------
      subroutine uzawa_gmres_split(l,u,b,binv,n)
      integer n
      real l(n),u(n),b(n),binv(n)
      integer i
      do i=1,n
c        l(i)=sqrt(binv(i))
c        u(i)=sqrt(b(i))
      
c        u(i)=sqrt(b(i))
c        l(i)=1./u(i)
      
c        l(i)=sqrt(binv(i))
         l(i)=1.
         u(i)=1./l(i)
      
      
c        if(abs(u(i)*l(i)-1.0).gt.1e-13)write(6,*) i,u(i)*l(i),' gmr_sp'
      enddo
      return
      end
      
c-----------------------------------------------------------------------
      subroutine uzawa_gmres_temp(a,b,n)
      integer n
      real a(n),b(n)
      integer i
      do i=1,n
         a(i)=sqrt(b(i))*a(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,x,h1,h2,n)

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
# 282 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 282 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 283 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 283 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
c     
c     w = A*x for pressure iteration
c     
      
      integer n
      real w(n),x(n),h1(n),h2(n)
      
      imsh = 1
      isd  = 1
      call axhelm (w,x,h1,h2,imsh,isd)
      call dssum  (w,nx1,ny1,nz1)
      call col2   (w,pmask,n)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine hmh_gmres(res,h1,h2,wt,iter)
      
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      
      

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
# 307 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 307 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 308 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 308 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/FDMH1" 1
c     
c     'FDMH1'
c     
# 4
      common /fdmh1r/ dd(lx1,9),fds(lx1*lx1,9),fdst(lx1*lx1,9)
     $              , elsize(3,lelt)
     $              , bhalf(lx1,ly1,lz1,lelt)
     $              , fds_acc(lx1,lx1,9),fdst_acc(lx1,lx1,9)
c     
      common /fdmh1i/ kfldfdm,ktype(lelt,3,0:4)
c     
      common /fdmh1l/ ifbhalf
      logical         ifbhalf
# 309 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 309 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GMRES" 1
c     parameter (lgmres=50) !max number of gmres iterations between rest
# 2
      common /gmres/  x(lx2*ly2*lz2*lelv)
     $              , r(lx2*ly2*lz2*lelv), w(lx2*ly2*lz2*lelv)
     $              , h(lgmres,lgmres), gamma(lgmres+1)
     $              , c(lgmres), s(lgmres)
      common /gmre1/  v(lx2*ly2*lz2*lelv,lgmres)
      common /gmre2/  z(lx2*ly2*lz2*lelv,lgmres)
      real x,r,w,h,gamma,c,s,v,z
      
      common /spltprec/ ml(lx2*ly2*lz2*lelv), mu(lx2*ly2*lz2*lelv)
      real ml,mu
      
c     w is a work vector
c     c and s store the Givens rotations
c     V stores the orthogonal Krylov subspace basis
c          -1
c     Z = M   V
# 310 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 310 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real             res  (lx1*ly1*lz1*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             wt   (lx1,ly1,lz1,lelv)
      
      common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)
      
      common /cgmres1/ y(lgmres)
      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      real alpha, l, temp
      integer outer
      
      logical iflag,if_hyb
      save    iflag,if_hyb
c     data    iflag,if_hyb  /.false. , .true. /
      data    iflag,if_hyb  /.false. , .false. /
      real    norm_fac
      save    norm_fac
      
      real*8 etime1,dnekclock
      
      
      n = nx1*ny1*nz1*nelv
      
      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres
      
      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml,mu,bm1,binvm1,nx1*ny1*nz1*nelv)
         norm_fac = 1./sqrt(volvm1)
      endif
      
      if (param(100).ne.2) call set_fdm_prec_h1b(d,h1,h2,nelv)
      
      call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c     
      iconv = 0
      call rzero(x,n)
      
      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
      
         if(iter.eq.0) then               !      -1
            call col3(r,ml,res,n)         ! r = L  res
c           call copy(r,res,n)
         else
            !update residual
            call copy  (r,res,n)                  ! r = res
            call ax    (w,x,h1,h2,n)              ! w = A x
            call add2s2(r,w,-1.,n)                ! r = r - w
                                                  !      -1
            call col2(r,ml,n)                     ! r = L   r
         endif
                                                  !            ______
         gamma(1) = sqrt(glsc3(r,r,wt,n))         ! gamma  = \/ (r,r) 
                                                  !      1
         if(iter.eq.0) then
            div0 = gamma(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif
      
         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             ! v  = r / gamma
                                                  !  1            1
         do j=1,m
            iter = iter+1
                                                  !       -1
            call col3(w,mu,v(1,j),n)              ! w  = U   v
                                                  !           j
      
c . . . . . Overlapping Schwarz + coarse-grid . . . . . . .
      
            etime2 = dnekclock()
      
c           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergen
            if (ifmgrid) then
               call h1mg_solve(z(1,j),w,if_hyb)   ! z  = M   w
            else                                  !  j
               kfldfdm = ndim+1
               if (param(100).eq.2) then
                   call h1_overlap_2 (z(1,j),w,pmask)
               else
                   call fdm_h1
     $               (z(1,j),w,d,pmask,vmult,nelv,ktype(1,1,kfldfdm),wk)
               endif
               call crs_solve_h1 (wk,w)           ! z  = M   w
               call add2         (z(1,j),wk,n)    !  j        
            endif
      
      
            call ortho        (z(1,j)) ! Orthogonalize wrt null space, i
            etime_p = etime_p + dnekclock()-etime2
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      
      
            call ax  (w,z(1,j),h1,h2,n)           ! w = A z
                                                  !        j
      
                                                  !      -1
            call col2(w,ml,n)                     ! w = L   w
      
c           !modified Gram-Schmidt
      
c           do i=1,j
c              h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
c                                                 !  i,j       i
      
c              call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
c           enddo                                 !          i,j  i
      
c           2-PASS GS, 1st pass:
      
            do i=1,j
               h(i,j)=vlsc3(w,v(1,i),wt,n)        ! h    = (w,v )
            enddo                                 !  i,j       i
      
            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs
      
            do i=1,j
               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
            enddo                                 !          i,j  i
      
      
c           2-PASS GS, 2nd pass:
c     
c           do i=1,j
c              wk1(i)=vlsc3(w,v(1,i),wt,n)        ! h    = (w,v )
c           enddo                                 !  i,j       i
c                                                 !
c           call gop(wk1,wk2,'+  ',j)             ! sum over P procs
c     
c           do i=1,j
c              call add2s2(w,v(1,i),-wk1(i),n)    ! w = w - h    v
c              h(i,j) = h(i,j) + wk1(i)           !          i,j  i
c           enddo
      
            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = sqrt(glsc3(w,w,wt,n))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)
      
            rnorm = abs(gamma(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nid.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence')
      

# 486
            if (rnorm .lt. tolpss) goto 900  !converged



# 490
            if (j.eq.m) goto 1000 !not converged, restart
      
            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1            
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,z(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
c        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
      enddo
 9000 continue
      
      divex = rnorm
      call copy(res,x,n)
      
      call ortho   (res) ! Orthogonalize wrt null space, if present
      
      etime1 = dnekclock()-etime1
      if (nid.eq.0) write(6,9999) istep,iter,divex,tolpss,div0,etime_p,
     &                            etime1,if_hyb
c     call flush_hack
 9999 format(i9,' PRES gmres:',i5,1p5e12.4,1x,l4)
      
      if (outer.le.2) if_hyb = .false.
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_overlap2
c     
c     Sets up the gather scatter and the SEM operators
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
# 537 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 537 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 538 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 538 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
      common /c_is1/ glo_num(lxs*lys*lzs*lelv)  
      integer*8 glo_num
      common /ivrtx/ vertex ((2**ldim)*lelt)
      common /handle/ gsh_dd
      integer vertex,gsh_dd
      
      mz = ndim-2
      nx = nx1+2
      ny = ny1+2
      nz = nz1+2*mz
      call setupds_no_crn(gsh_dd,nx,ny,nz,nelv,nelgv,vertex,glo_num)
      call swap_lengths ! Set up Matrices for FDM
      call gen_fast_g
      
      return
      end
c-----------------------------------------------------------------------
      subroutine h1_overlap_2(u,v,mask)
c     
c     Local overlapping Schwarz solves with overlap of 2
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
# 561 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 561 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 562 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 562 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
      
      common /cwork1/ v1(lxs,lys,lzs,lelt)
      common /handle/ gsh_dd
      integer gsh_dd
      
      real u(lx1,lx1,lz1,1),v(1),mask(1)
      integer e
      
      n = nx1*ny1*nz1*nelfld(ifield)
      
      call dd_swap_vals(v1,v,gsh_dd)
      
      iz1 = 0
      if (if3d) iz1=1
      do ie=1,nelfld(ifield)
         do iz=1,nz1
         do iy=1,ny1
         do ix=1,nx1
            u(ix,iy,iz,ie) = v1(ix+1,iy+1,iz+iz1,ie)
         enddo
         enddo
         enddo
      enddo
      
      call dssum (u,nx1,ny1,nz1)
      call col2  (u,mask,n)
      
      
      return
      end
c-----------------------------------------------------------------------
      subroutine dd_swap_vals(v1,v0,gsh_dd)
      

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
# 597 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 597 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 598 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 598 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
      common /work1/ w1(lxs,lys,lzs)           ! work arrarys for locals
     $              ,w2(lxs,lys,lzs)
      integer gsh_dd
      
      real v1(lxs,lys,lzs,lelt)
      real v0(lx1,ly1,lz1,lelt)
      
      integer e
      
      mz = ndim-2
      
      nx = nx1+2
      ny = ny1+2
      nz = nz1+2*mz
      
      n  = nx1*ny1*nz1*nelv
      m  = (nx1+2)*(ny1+2)*(nz1+2*mz)*nelv
      
      do e=1,nelv
         call rzero_g        (v1,e,nx,ny,nz)
         call fill_interior_g(v1,v0,e,nx1,nz1,mz,nelv)    ! v0      --> 
         call dface_ext_g    (v1,2,e,nx,nz)               ! v1(int) --> 
      enddo
c     
c                 ~ ~T
c     This is the Q Q  part
c     
      call gs_op(gsh_dd,v1,1,1,0)  ! 1 ==> +         ! swap v1 & add val
      
      do e =1,nelv
         call dface_add1si_g  (v1,-1.,2,e,nx,nz)
         call fastdm1_g       (v1(1,1,1,e),e,w1,w2)
         call s_face_to_int2_g(v1,-1.,2,e,nx,nz)
      enddo
c     
c     Exchange/add elemental solutions
c     
      call gs_op              (gsh_dd,v1,1,1,0)
      do e =1,nelv
         call s_face_to_int2_g(v1,1.,2,e,nx,nz)
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_fast_g
c     
c     Generate fast diagonalization matrices for each element
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
# 649 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 649 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
      
# 650 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 650 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 651 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 651 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 652 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 652 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 653 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 653 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
      parameter (lxss=lxs*lxs)
      common /fastg/  sr(lxss,2,lelv),ss(lxss,2,lelv),st(lxss,2,lelv)
     $             ,  df(lxs*lys*lzs,lelv)
      
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
      
      integer lbr,rbr,lbs,rbs,lbt,rbt
      
      
      call load_semhat_weighted   !   Fills the SEMHAT arrays
      
      
      ierr = 0
      do ie=1,nelv
      
         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,3,ierr)
c     
c        Set up matrices for each element.
c     
         call set_up_fast_1D_sem_g( sr(1,1,ie),lr,nr ,lbr,rbr
     $                      ,llr(ie),lmr(ie),lrr(ie),ie)
         call set_up_fast_1D_sem_g( ss(1,1,ie),ls,ns ,lbs,rbs
     $                      ,lls(ie),lms(ie),lrs(ie),ie)
         if (if3d) then
            call set_up_fast_1D_sem_g( st(1,1,ie),lt,nt ,lbt,rbt
     $                      ,llt(ie),lmt(ie),lrt(ie),ie)
         endif
c     
c        Set up diagonal inverse
c     
         if (if3d) then
            eps = 1.e-5 * (vlmax(lr(2),nr-2)
     $                  +  vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
            l   = 1
            do k=1,nt
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j) + lt(k)
               if (diag.gt.eps) then
                  df(l,ie) = 1.0/diag
               else
c                 write(6,3) ie,'Reset Eig in gen fast:',i,j,k,l
c    $                         ,eps,diag,lr(i),ls(j),lt(k)
c   3             format(i6,1x,a21,4i5,1p5e12.4)
                  df(l,ie) = 0.0
               endif
               l = l+1
            enddo
            enddo
            enddo
         else
            eps = 1.e-5*(vlmax(lr(2),nr-2) + vlmax(ls(2),ns-2))
            l   = 1
            do j=1,ns
            do i=1,nr
               diag = lr(i) + ls(j)
               
               if (diag.gt.eps) then
                  df(l,ie) = 1.0/diag
                  
               else
c                 write(6,2) ie,'Reset Eig in gen fast:',i,j,l
c    $                         ,eps,diag,lr(i),ls(j)
c   2             format(i6,1x,a21,3i5,1p4e12.4)
                  df(l,ie) = 0.0
               endif
               l = l+1
            enddo
            enddo
         endif
c     
c        Next element ....
c     
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem_g(s,lam,n,lbc,rbc,ll,lm,lr,ie)

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
# 741 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 741 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/SEMHAT" 1
# 1
      parameter (lr2=2*lx1*lx1)
      common /ahat/ ah (lr2),bh (lr2),ch (lr2),dh (lr2)
     $             ,dph(lr2),jph(lr2),zh (lr2),wh (lr2)  ! Pressure GLL
     $             ,bgl(lr2),zgl(lr2),dgl(lr2),jgl(lr2)  ! Pressure GL 
      real ah,bh,ch,dh,dph,jph,zh,wh,bgl,zgl,dgl,jgl
c     
      parameter (l3=lx1*(lx1+1)*(lx1+2)/3)
      parameter (l2=lx1*(lx1+1)/2)
      common /hata/ dd(l3),zp(l2),ww(l2)   ! Aggregate arrays
# 742 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 742 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
      parameter (lr3=2*lxs*lxs)
      common /fast1dsem/ g(lr3),w(lr3)
      
      real g,w
      real s(1),lam(1),ll,lm,lr
      integer lbc,rbc
      
      integer bb0,bb1,eb0,eb1,n,n1
      logical l,r
      
      n=nx1-1
      !bcs on E are from normal vel component
      if(lbc.eq.2 .or. lbc.eq.3) then !wall,sym  - Neumann
         eb0=0
         bb0=0
      else !outflow,element                      - Dirichlet
         eb0=1
         bb0=1
      endif
      if(rbc.eq.2 .or. rbc.eq.3) then !wall,sym  - Neumann
         eb1=n
         bb1=n
      else !outflow,element                      - Dirichlet
         eb1=n-1
         bb1=n-1
      endif
      l = (lbc.eq.0)
      r = (rbc.eq.0)
      
c     calculate A tilde operator
      call set_up_fast_1D_sem_op_a(s,eb0,eb1,l,r,ll,lm,lr,ah)
c     calculate B tilde operator
      call set_up_fast_1D_sem_op_b(g,bb0,bb1,l,r,ll,lm,lr,bh)
      
      n=n+3
c     call outmat   (s,n,n,'  A   ',ie)
c     call outmat   (g,n,n,'  B   ',ie)
      call generalev(s,g,lam,n,w)
      if(.not.l) call row_zero(s,n,n,1)
      if(.not.r) call row_zero(s,n,n,n)
      call transpose(s(n*n+1),n,s,n) ! compute the transpose of s
c     call outmat   (s,n,n,'  S   ',ie)
c     call outmat   (s(n*n+1),n,n,'  St  ',1)
c     call exitt
      
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem_op_a(g,b0,b1,l
     $                                   ,r,ll,lm,lr,ah)
c            -1 T
c     G = J B  J
c     
c     gives the inexact restriction of this matrix to
c     an element plus two node on either side
c     
c     g - the output matrix
c     b0, b1 - the range for Bhat indices for the element
c              (enforces boundary conditions)
c     l, r - whether there is a left or right neighbor
c     ll,lm,lr - lengths of left, middle, and right elements
c     ah - hat matrix for A or B
c     
c     result is inexact because:
c        neighbor's boundary condition at far end unknown
c        length of neighbor's neighbor unknown
c        (these contribs should be small for large N and
c         elements of nearly equal size)
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
# 814 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 814 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      real g(0:lx1+1,0:lx1+1)
      real ah(0:lx1-1,0:lx1-1)
      real ll,lm,lr
      integer b0,b1
      logical l,r
      
      real bl,bm,br
      integer n
      
      n =nx1-1
      
c     
c     compute the weight of A hat
c     
      if (l) bl = 2. /ll
             bm = 2. /lm
      if (r) br = 2. /lr
      
      call rzero(g,(n+3)*(n+3))
      do j=1,n+1
         do i=1,n+1
            g(i,j) = ah(i-1,j-1)*bm
         enddo
      enddo
      
      if (l) then
         g(0,0) = g(0,0) + bl*ah(n-1,n-1)
         g(0,1) = g(0,1) + bl*ah(n-1,n  )
         g(1,0) = g(1,0) + bl*ah(n  ,n-1)
         g(1,1) = g(1,1) + bl*ah(n  ,n  )
      elseif (b0.eq.0) then          !Neumann BC
         g(0,0) = 1.
      else                           !Dirichlet BC
         g(0,0)     = 1.
         g(1,1)     = 1.
         do i=2,n+1
            g(i,1)  = 0.
            g(1,i)  = 0.
         enddo
      endif
      
      if (r) then
         g(n+1,n+1) = g(n+1,n+1) + br*ah(0,0)
         g(n+1,n+2) = g(n+1,n+2) + br*ah(0,1)
         g(n+2,n+1) = g(n+2,n+1) + br*ah(0,1)
         g(n+2,n+2) = g(n+2,n+2) + br*ah(1,1)
      elseif (b1.eq.n) then          !Neumann BC
         g(n+2,n+2)=1.
      else                           !Dirichlet BC
         g(n+2,n+2)   = 1.
         g(n+1,n+1)   = 1.
         do i=1,n
            g(i,n+1)  = 0.
            g(n+1,i)  = 0.
         enddo
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem_op_b(g,b0,b1,l
     $                                   ,r,ll,lm,lr,bh)
c            -1 T
c     G = J B  J
c     
c     gives the inexact restriction of this matrix to
c     an element plus two node on either side
c     
c     g - the output matrix
c     b0, b1 - the range for Bhat indices for the element
c              (enforces boundary conditions)
c     l, r - whether there is a left or right neighbor
c     ll,lm,lr - lengths of left, middle, and right elements
c     bh - hat matrix for  B
c     
c     result is inexact because:
c        neighbor's boundary condition at far end unknown
c        length of neighbor's neighbor unknown
c        (these contribs should be small for large N and
c         elements of nearly equal size)
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
# 896 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 896 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      real g(0:lx1+1,0:lx1+1)
      real bh(0:lx1-1)
      real ll,lm,lr
      integer b0,b1
      logical l,r
      
      real bl,bm,br
      integer n
      
      n =nx1-1
c     
c     compute the weight of B hat
c     
      if (l) bl = ll / 2.
             bm = lm / 2.
      if (r) br = lr / 2.
      
      call rzero(g,(n+3)*(n+3))
      do i=1,n+1
            g(i,i) = bh(i-1)*bm
      enddo
      
      if (l) then
         g(0,0) = g(0,0) + bl*bh(n-1)
         g(1,1) = g(1,1) + bl*bh(n  )
      elseif (b0.eq.0) then          !Neumann BC
         g(0,0) = 1.
      else                           !Dirichlet BC
         g(0,0)     = 1.
         g(1,1)     = 1.
         do i=2,n+1
            g(i,1)  = 0.
            g(1,i)  = 0.
         enddo
      endif
      
      if (r) then
         g(n+1,n+1) = g(n+1,n+1) + br*bh(0)
         g(n+2,n+2) = g(n+2,n+2) + br*bh(1)
      elseif (b1.eq.n) then          !Neumann BC
         g(n+2,n+2)=1.
      else                           !Dirichlet BC
         g(n+2,n+2)   = 1.
         g(n+1,n+1)   = 1.
         do i=1,n
            g(i,n+1)  = 0.
            g(n+1,i)  = 0.
         enddo
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine fill_interior_g(v1,v,e,nx,nz,iz1,nel)
      
      real v1(nx+2,nx+2,nz+2*iz1,nel) ! iz1=ndim-2
      real v (nx  ,nx  ,nz      ,nel)
      integer e
      
      ny = nx
      
      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
         v1(ix+1,iy+1,iz+iz1,e) = v(ix,iy,iz,e)
      enddo
      enddo
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine dface_ext_g(x,t,e,nx,nz)
c     Extend interior to face of element
      

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
# 972 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 972 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
      
# 973 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 973 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      real x(nx,nx,nz,1)
      integer e,t
      
      ny = nx
      
      if (if3d) then
        do iz=2,nz-1
        do ix=2,nx-1
            x(ix,1 ,iz,e) = x(ix, 1+t,iz,e)
            x(ix,ny,iz,e) = x(ix,ny-t,iz,e)
        enddo
        enddo
      
        do iz=2,nz-1
        do iy=2,ny-1
            x(1 ,iy,iz,e) = x( 1+t,iy,iz,e)
            x(nx,iy,iz,e) = x(nx-t,iy,iz,e)
        enddo
        enddo
      
        do iy=2,ny-1
        do ix=2,nx-1
            x(ix,iy,1 ,e) = x(ix,iy, 1+t,e)
            x(ix,iy,nz,e) = x(ix,iy,nz-t,e)
        enddo
        enddo
      else
         do ix=2,nx-1
            x(ix,1 ,1,e) = x(ix, 1+t,1,e)
            x(ix,ny,1,e) = x(ix,ny-t,1,e)
         enddo
         do iy=2,ny-1
            x(1 ,iy,1,e) = x( 1+t,iy,1,e)
            x(nx,iy,1,e) = x(nx-t,iy,1,e)
         enddo
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine dface_add1si_g(x,c,t,e,nx,nz)
c     Scale interior and add to face of element
      

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
# 1017 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1017 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
      
# 1018 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1018 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      real x(nx,nx,nz,1)
      
      integer e,t
      
      ny = nx
      
      if (if3d) then
      
         do iz=2,nz-1
         do ix=2,nx-1
            x(ix,1 ,iz,e) = x(ix,1 ,iz,e) + c*x(ix, 1+t,iz,e)
            x(ix,ny,iz,e) = x(ix,ny,iz,e) + c*x(ix,ny-t,iz,e)
         enddo
         enddo
      
         do iz=2,nz-1
         do iy=2,ny-1
            x(1 ,iy,iz,e) = x(1 ,iy,iz,e) + c*x( 1+t,iy,iz,e)
            x(nx,iy,iz,e) = x(nx,iy,iz,e) + c*x(nx-t,iy,iz,e)
         enddo
         enddo
      
         do iy=2,ny-1
         do ix=2,nx-1
            x(ix,iy,1 ,e) = x(ix,iy,1 ,e) + c*x(ix,iy, 1+t,e)
            x(ix,iy,nz,e) = x(ix,iy,nz,e) + c*x(ix,iy,nz-t,e)
         enddo
         enddo
      
      else  ! 2D
      
         do ix=2,nx-1
            x(ix,1 ,1,e) = x(ix,1 ,1,e) + c*x(ix, 1+t,1,e)
            x(ix,ny,1,e) = x(ix,ny,1,e) + c*x(ix,ny-t,1,e)
         enddo
         do iy=2,ny-1
            x(1  ,iy,1,e) = x(1 ,iy,1,e) + c*x( 1+t,iy,1,e)
            x(nx ,iy,1,e) = x(nx,iy,1,e) + c*x(nx-t,iy,1,e)
         enddo
      
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine fastdm1_g(R,ie,w1,w2)
c     
c     Fast diagonalization solver for FEM on mesh 1
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
# 1068 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1068 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
      parameter (lxss=lxs*lxs)
      common /fastg/  sr(lxss,2,lelv),ss(lxss,2,lelv),st(lxss,2,lelv)
     $             ,  df(lxs*lys*lzs,lelv)
      
      parameter (lxyz = lxs*lys*lzs)
      
      real r(1),w1(1),w2(1)
      
      nx = nx1+2
c     
c      T
c     S  r
      call tensr3 (w1,nx,r ,nx,sr(1,2,ie),ss(1,1,ie),st(1,1,ie),w2)
c     
c     
c      -1 T
c     D  S  r
c     
      call col2   (w1,df(1,ie),lxyz)
c     
c     
c        -1 T
c     S D  S  r
c     
      call tensr3 (r ,nx,w1,nx,sr(1,1,ie),ss(1,2,ie),st(1,2,ie),w2)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine s_face_to_int2_g(x,c,t,e,nx,nz)
c     
c     Scale face and add to interior of element
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
# 1103 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1103 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
      
# 1104 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1104 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      real x(nx,nx,nz,1)
      integer t,e
      
      ny=nx
      
      
      if (if3d) then
      
         do iz=2,nz-1
         do ix=2,nx-1
            x(ix, 1+t,iz,e) = c*x(ix,1 ,iz,e) + x(ix, 1+t,iz,e)
            x(ix,ny-t,iz,e) = c*x(ix,ny,iz,e) + x(ix,ny-t,iz,e)
         enddo
         enddo
      
         do iz=2,nz-1
         do iy=2,ny-1
            x( 1+t,iy,iz,e) = c*x(1 ,iy,iz,e) + x( 1+t,iy,iz,e)
            x(nx-t,iy,iz,e) = c*x(nx,iy,iz,e) + x(nx-t,iy,iz,e)
         enddo
         enddo
      
         do iy=2,ny-1
         do ix=2,nx-1
            x(ix,iy, 1+t,e) = c*x(ix,iy,1 ,e) + x(ix,iy, 1+t,e)
            x(ix,iy,nz-t,e) = c*x(ix,iy,nz,e) + x(ix,iy,nz-t,e)
         enddo
         enddo
      
      else
c     2D
         do ix=2,nx-1
            x(ix, 1+t,1,e) = c*x(ix,1 ,1,e) + x(ix, 1+t,1,e)
            x(ix,ny-t,1,e) = c*x(ix,ny,1,e) + x(ix,ny-t,1,e)
         enddo
         do iy=2,ny-1
            x( 1+t,iy,1,e) = c*x(1 ,iy,1,e) + x( 1+t,iy,1,e)
            x(nx-t,iy,1,e) = c*x(nx,iy,1,e) + x(nx-t,iy,1,e)
         enddo
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldr_g(x,txt10,nx,nz,ichk)

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
# 1150 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1150 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 1151 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1151 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      real x(nx,nx,nz,lelt)
      character*10 txt10
      
      integer idum,e
      save idum
      data idum /3/
      
      if (idum.lt.0) return
      
      ny = nx
      
      
      mtot = nx*ny*nz*nelv
      if (nx.gt.8.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
      
      nell = nelt
      rnel = nell
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nell-ne+1
      k = 1
      do ie=1,1
         ne = 0
         write(6,116) txt10,k,ie,xmin,xmax,istep,time,nx
         do l=12,0,-4
            write(6,117)
            do j=ny,1,-1
              if (nx.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.8) write(6,118) ((x(i,j,k,e+l),i=1,nx),e=1,4)
            enddo
         enddo
      enddo
      
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f7.1,5x,6f7.1,5x,6f7.1,5x,6f7.1)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
  118 FORMAT(8f12.9)
      
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5,i2)
  117 FORMAT(' ')
      
      if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldi_g(x,txt10,nx,nz,ichk)

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
# 1213 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1213 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 1214 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1214 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      integer x(nx,nx,nz,lelt)
      character*10 txt10
      
      integer idum,e,xmin,xmax
      save idum
      data idum /3/
      
      if (idum.lt.0) return
      
      ny = nx
      
      
      mtot = nx*ny*nz*nelv
      if (nx.gt.8.or.nelv.gt.16) return
      xmin = iglmin(x,mtot)
      xmax = iglmax(x,mtot)
      
      nell = nelt
      rnel = nell
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nell-ne+1
      k = 1
      do ie=1,1
         ne = 0
         write(6,116) txt10,k,ie,xmin,xmax,istep,time,nx
         do l=12,0,-4
            write(6,117)
            do j=ny,1,-1
              if (nx.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx),e=1,4)
              if (nx.eq.8) write(6,118) ((x(i,j,k,e+l),i=1,nx),e=1,4)
            enddo
         enddo
      enddo
      
  102 FORMAT(4(2i9,2x))
  103 FORMAT(4(3i9,2x))
  104 FORMAT(4(4i7,2x))
  105 FORMAT(5i9,10x,5i9)
  106 FORMAT(6i7,5x,6i7,5x,6i7,5x,6i7)
  107 FORMAT(7i8,5x,7i8)
  108 FORMAT(8i8,4x,8i8)
  118 FORMAT(8i12)
      
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2i12,/,
     $    5X,'       X            ','Step  =',I9,f15.5,i2)
  117 FORMAT(' ')
      
      if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine setupds_no_crn(gs_h,nx,ny,nz,nel,melg,vertex,glo_num)

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
# 1276 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1276 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
      
# 1277 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1277 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 1278 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1278 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/NONCON" 1
c     
c234567
c     
# 4
       common /allr/ 
     $      umult(lx1*ly1*lz1*lelt)
     $   ,  Jmat(lx1,lx1,2,maxmor)
     $   ,  xsp(lx1,lz1),ysp(lx1,lz1),zsp(lx1,lz1)
     $   ,  xch(lx1,lz1),ych(lx1,lz1),zch(lx1,lz1)
     $   ,  rtwid(lx1),stwid(lx1)
     $   ,  dtrk(ldim)
     $   ,  rs( 2 , 2 , 2 )
      real Jmat
c     
       common /alli/ 
     $      noncon_f(maxmor)
     $   ,  noncon_e(maxmor)
     $   ,  noncon_ip(maxmor)
     $   ,  mortar(6,lelt)
     $   ,  imin(3,2)
     $   ,  mort_m
c     
       common /logg/ ifnc,ifhalf,ifJt(maxmor)
       logical ifnc,ifhalf,ifJt
# 1279 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1279 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      integer   gs_h,vertex(1),e
      integer*8 ngv,glo_num(nx,ny,nz,nel)
      
      
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      
c     set up the global numbering
      call set_vert(glo_num,ngv,nx,nel,vertex,.false.)
      
c     zero out corners
      mz1 = max(1,nz-1)
      do e=1,nel
      do k=1,nz,mz1
      do j=1,ny,ny-1
      do i=1,nx,nx-1
         glo_num(i,j,k,e) = 0
      enddo
      enddo
      enddo
      enddo
      
      
      ntot = nx*ny*nz*nel
      
      t0 = dnekclock()
      call gs_setup(gs_h,glo_num,ntot,nekcomm,mp) ! initialize gs code
      t1 = dnekclock()
      
      et = t1-t0
      
c     call gs_chkr(glo_num)
      
      if (nid.eq.0) then
         write(6,1) et,nx,nel,ntot,ngv,gs_h
    1    format('   gs_init time',1pe11.4,' seconds ',i3,4i10)
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine rzero_g(a,e,nx,ny,nz)
      
      real a(nx,ny,nz,1)
      integer e
      
      do  i=1,nx
      do  j=1,ny
      do  k=1,nz
          a(i,j,k,e) = 0.0
      enddo
      enddo
      enddo
      
      return
      END
c-----------------------------------------------------------------------
      
      
      
C====================== OpenACC =================================
      

      
c-----------------------------------------------------------------------
# 1343
      subroutine uzawa_gmres_split_acc(l,u,b,binv,n)
      integer n
      real l(n),u(n),b(n),binv(n)
      integer i
      
!$ACC DATA PRESENT(l,u,b,binv)
!$ACC PARALLEL LOOP
      do i=1,n
c        l(i)=sqrt(binv(i))
c        u(i)=sqrt(b(i))
      
c        u(i)=sqrt(b(i))
c        l(i)=1./u(i)
      
c        l(i)=sqrt(binv(i))
         l(i)=1.
         u(i)=1./l(i)
      
      
c        if(abs(u(i)*l(i)-1.0).gt.1e-13)write(6,*) i,u(i)*l(i),' gmr_sp'
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      end
      
c-----------------------------------------------------------------------
      subroutine ax_acc(w,x,h1,h2,n)

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
# 1373 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1373 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 1374 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1374 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
c     
c     w = A*x for pressure iteration
c     
      
      integer n
      real w(n),x(n),h1(n),h2(n)
      
      imsh = 1
      isd  = 1
      call axhelm_acc (w,x,h1,h2,imsh,isd)
      call gdssum_acc (w,nx1,ny1,nz1)
      call col2_acc   (w,pmask,n)
      
      return
      end
      
c-----------------------------------------------------------------------
      subroutine hmh_gmres_acc(res,h1,h2,wt,iter)
      
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
      
      

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
# 1399 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1399 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

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
# 1400 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1400 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/FDMH1" 1
c     
c     'FDMH1'
c     
# 4
      common /fdmh1r/ dd(lx1,9),fds(lx1*lx1,9),fdst(lx1*lx1,9)
     $              , elsize(3,lelt)
     $              , bhalf(lx1,ly1,lz1,lelt)
     $              , fds_acc(lx1,lx1,9),fdst_acc(lx1,lx1,9)
c     
      common /fdmh1i/ kfldfdm,ktype(lelt,3,0:4)
c     
      common /fdmh1l/ ifbhalf
      logical         ifbhalf
# 1401 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1401 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/GMRES" 1
c     parameter (lgmres=50) !max number of gmres iterations between rest
# 2
      common /gmres/  x(lx2*ly2*lz2*lelv)
     $              , r(lx2*ly2*lz2*lelv), w(lx2*ly2*lz2*lelv)
     $              , h(lgmres,lgmres), gamma(lgmres+1)
     $              , c(lgmres), s(lgmres)
      common /gmre1/  v(lx2*ly2*lz2*lelv,lgmres)
      common /gmre2/  z(lx2*ly2*lz2*lelv,lgmres)
      real x,r,w,h,gamma,c,s,v,z
      
      common /spltprec/ ml(lx2*ly2*lz2*lelv), mu(lx2*ly2*lz2*lelv)
      real ml,mu
      
c     w is a work vector
c     c and s store the Givens rotations
c     V stores the orthogonal Krylov subspace basis
c          -1
c     Z = M   V
# 1402 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1402 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/POPENACC" 1
# 1
      common /spltacc/  ml_acc(lx2*ly2*lz2*lelv)
     $              ,   mu_acc(lx2*ly2*lz2*lelv)
      
      real ml_acc,mu_acc
      
      parameter (lta=lx1*ly1*lz1*lelt)
      common /mgdssum/ mg_lgl1(-1:lta,6)
     $              ,  mg_lgl2(-1:lta,6),mg_ptr(lta,6)
# 1403 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f" 2
# 1403 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/gmres.f"
      
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real             res  (lx1*ly1*lz1*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             wt   (lx1,ly1,lz1,lelv)
      
      common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)
      
      common /gmre3/  v_acc(lx2*ly2*lz2*lelv,lgmres)
     $     ,          z_acc(lx2*ly2*lz2*lelv,lgmres)
      
      common /cgmres1/ y(lgmres)
      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      real alpha, l, temp
      integer outer
      
      logical iflag,if_hyb
      save    iflag,if_hyb
c     data    iflag,if_hyb  /.false. , .true. /
      data    iflag,if_hyb  /.false. , .false. /
      real    norm_fac
      save    norm_fac
      
      real*8 etime1,dnekclock
      
      n = nx1*ny1*nz1*nelv
      
      etime1 = dnekclock()
      etime_p = 0.
      divex = 0.
      iter  = 0
      m     = lgmres
      
      
c === there are different variables between v,z
!$ACC  DATA PRESENT(res,h1,h2,wt,z_acc,v_acc)
!$ACC&      PRESENT(d,pmask,vmult,ktype,bm1,binvm1)
!$ACC&      CREATE(r,x,w,wk)
!$ACC&      PRESENT(ml_acc,mu_acc)
      
      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split_acc
     $        (ml_acc,mu_acc,bm1,binvm1,nx1*ny1*nz1*nelv)
      
         norm_fac = 1./sqrt(volvm1)
      endif
      
      if (param(100).ne.2) call set_fdm_prec_h1b(d,h1,h2,nelv)
      
!$ACC UPDATE DEVICE(d)
      
!$ACC UPDATE HOST(res)
      call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)
      
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c     
      iconv = 0
      
      call rzero_acc(x,n)
      outer = 0
      do while (iconv.eq.0.and.iter.lt.500)
         outer = outer+1
      
         if(iter.eq.0) then               !      -1
            call col3_acc(r,ml_acc,res,n)         ! r = L  res
c           call copy(r,res,n)
         else
            !update residual
            call copy_acc  (r,res,n)                  ! r = res
            call ax_acc    (w,x,h1,h2,n)              ! w = A x
            call add2s2_acc(r,w,-1.,n)                ! r = r - w
                                                  !      -1
            call col2_acc(r,ml_acc,n)                     ! r = L   r
         endif
                                                  !            ______
         gamma(1) = sqrt(glsc3_acc(r,r,wt,n))     ! gamma  = \/ (r,r) 
                                                  !      1
         if(iter.eq.0) then
            div0 = gamma(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif
      
         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
      
         call cmult2_acc(v_acc(1,1),r,temp,n)             ! v  = r / gam
                                                   !  1            1
         do j=1,m
            iter = iter+1
                                                  !       -1
            call col3_acc(w,mu_acc,v_acc(1,j),n)          ! w  = U   v
                                                  !           j
      
c . . . . . Overlapping Schwarz + coarse-grid . . . . . . .
      
            etime2 = dnekclock()
      
c           if (outer.gt.2) if_hyb = .true.       ! Slow outer convergen
            if (ifmgrid) then
               call h1mg_solve_acc(z_acc(1,j),w,if_hyb)   ! z  = M   w
            else                                  !  j
               kfldfdm = ndim+1
               if (param(100).eq.2) then
                   call h1_overlap_2 (z_acc(1,j),w,pmask)
               else
                   call fdm_h1_acc
     $         (z_acc(1,j),w,d,pmask,vmult,nelv,ktype(1,1,kfldfdm),wk)
               endif
               call crs_solve_h1_acc (wk,w)           ! z  = M   w
               call add2_acc         (z_acc(1,j),wk,n)    !  j        
            endif
      
            call ortho        (z_acc(1,j)) ! Orthogonalize wrt null spac
            etime_p = etime_p + dnekclock()-etime2
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
      
      
            call ax_acc  (w,z_acc(1,j),h1,h2,n)           ! w = A z
                                                  !        j
      
                                                  !      -1
            call col2_acc(w,ml_acc,n)                     ! w = L   w
      
c           !modified Gram-Schmidt
      
c           do i=1,j
c              h(i,j)=glsc3(w,v_acc(1,i),wt,n)        ! h    = (w,v )
c                                                 !  i,j       i
      
c              call add2s2(w,v_acc(1,i),-h(i,j),n)    ! w = w - h    v
c           enddo                                 !          i,j  i
      
c           2-PASS GS, 1st pass:
      
      
            do i=1,j
               h(i,j)=vlsc3_acc(w,v_acc(1,i),wt,n)        ! h    = (w,v 
            enddo                                 !  i,j       i
      
            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs
      
            do i=1,j
               call add2s2_acc(w,v_acc(1,i),-h(i,j),n)    ! w = w - h   
            enddo                                 !          i,j  i
      
c           2-PASS GS, 2nd pass:
c     
c           do i=1,j
c              wk1(i)=vlsc3(w,v_acc(1,i),wt,n)        ! h    = (w,v )
c           enddo                                 !  i,j       i
c                                                 !
c           call gop(wk1,wk2,'+  ',j)             ! sum over P procs
c     
c           do i=1,j
c              call add2s2(w,v_acc(1,i),-wk1(i),n)    ! w = w - h    v
c              h(i,j) = h(i,j) + wk1(i)           !          i,j  i
c           enddo
      
            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = sqrt(glsc3_acc(w,w,wt,n))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)
      
            rnorm = abs(gamma(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nid.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence')
      

# 1595
            if (rnorm .lt. tolpss) goto 900  !converged



# 1599
            if (j.eq.m) goto 1000 !not converged, restart
      
            temp = 1./alpha
            call cmult2_acc(v_acc(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1            
      
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2_acc(x,z_acc(1,i),c(i),n)     ! x = x + c  z
         enddo                               !          i  i
      
c        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
      
      enddo
 9000 continue
      
      divex = rnorm
      call copy_acc(res,x,n)
      
      call ortho_acc   (res) ! Orthogonalize wrt null space, if present
      
!$ACC END DATA
      
      etime1 = dnekclock()-etime1
      if (nid.eq.0) write(6,9999) istep,iter,divex,tolpss,div0,etime_p,
     &                            etime1,if_hyb
c     call flush_hack
 9999 format(i9,' PRES gmres:',i5,1p5e12.4,1x,l4)
      
      if (outer.le.2) if_hyb = .false.
      
      return
      end
      
      
      

      
