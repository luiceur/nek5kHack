# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C**********************************************************************
C     
C     ROUTINES FOR ESITMATING AND CALCULATING EIGENVALUES
C     USED IN NEKTON
C     
C**********************************************************************
C     
      SUBROUTINE ESTEIG
C--------------------------------------------------------------
C     
C     Estimate eigenvalues
C     
C-------------------------------------------------------------

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
# 15 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 15 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 16 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 16 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
      
# 17 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 17 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 18 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 18 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 19 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 19 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C     
      NTOT1=NX1*NY1*NZ1*NELFLD(IFIELD)
      XMIN = GLMIN(XM1,NTOT1)
      XMAX = GLMAX(XM1,NTOT1)
      YMIN = GLMIN(YM1,NTOT1)
      YMAX = GLMAX(YM1,NTOT1)
      IF (IF3D) THEN
         ZMIN = GLMIN(ZM1,NTOT1)
         ZMAX = GLMAX(ZM1,NTOT1)
      ELSE
         ZMIN = 0.0
         ZMAX = 0.0
      ENDIF
C     
      XX = XMAX - XMIN
      YY = YMAX - YMIN
      ZZ = ZMAX - ZMIN
      RXY = XX/YY
      RYX = YY/XX
      RMIN = RXY
      IF (RYX .LT. RMIN) RMIN = RYX
      IF (NDIM .EQ. 3) THEN
         RXZ = XX/ZZ
         RZX = ZZ/XX
         RYZ = YY/ZZ
         RZY = ZZ/YY
         IF (RXZ .LT. RMIN) RMIN = RXZ
         IF (RZX .LT. RMIN) RMIN = RZX
         IF (RYZ .LT. RMIN) RMIN = RYZ
         IF (RZY .LT. RMIN) RMIN = RZY
      ENDIF
C     
      XX2    = 1./XX**2
      YY2    = 1./YY**2
      XYZMIN = XX2
      XYZMAX = XX2+YY2
      IF (YY2 .LT. XYZMIN) XYZMIN = YY2
      IF (NDIM .EQ. 3) THEN
         ZZ2 = 1./ZZ**2
         XYZMAX = XYZMAX+ZZ2
         IF (ZZ2 .LT. XYZMIN) XYZMIN = ZZ2
      ENDIF
C     
      one    = 1.
      PI     = 4.*ATAN(one)
      RATIO  = XYZMIN/XYZMAX
      EIGAE  = PI*PI*XYZMIN
      EIGGE  = EIGGA
      IF (NDIM .EQ. 2) EIGAA = PI*PI*(XX2+YY2)/2.
      IF (NDIM .EQ. 3) EIGAA = PI*PI*(XX2+YY2+ZZ2)/3.
      IF (IFAXIS)      EIGAA = .25*PI*PI*YY2
      EIGAS  = 0.25*RATIO
      EIGGS  = 2.0
C     
      IF (NID.EQ.0 .AND. ISTEP.LE.0) THEN
         WRITE (6,*) ' '
         WRITE (6,*) 'Estimated eigenvalues'
         WRITE (6,*) 'EIGAA = ',EIGAA
         WRITE (6,*) 'EIGGA = ',EIGGA
         IF (IFFLOW) THEN
            WRITE (6,*) 'EIGAE = ',EIGAE
            WRITE (6,*) 'EIGAS = ',EIGAS
            WRITE (6,*) 'EIGGE = ',EIGGE
            WRITE (6,*) 'EIGGS = ',EIGGS
         ENDIF
         WRITE (6,*) ' '
      ENDIF
C     
      RETURN
      END
C     
      SUBROUTINE EIGENV
C-----------------------------------------------------------------------
C     
C     Compute the following eigenvalues:
C     EIGAA  = minimum eigenvalue of the matrix A  (=Laplacian)
C     EIGAE  = minimum eigenvalue of the matrix E  (=DB-1DT)
C     EIGAS  = minimum eigenvalue of the matrix S  (=DA-1DT)
C     EIGAST = minimum eigenvalue of the matrix St (=D(A+B/dt)-1DT
C     EIGGA  = maximum eigenvalue of the matrix A
C     EIGGS  = maximum eigenvalue of the matrix S
C     EIGGE  = maximum eigenvalue of the matrix E
C     EIGGST = maximum eigenvalue of the matrix St
C     
C     Method : Power method/Inverse iteration & Rayleigh quotient wo shi
C     
C-----------------------------------------------------------------------

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
# 107 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 107 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 108 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 108 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
      
# 109 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 109 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 110 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 110 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 111 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 111 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C     
      COMMON /SCRVH/ H1 (LX1,LY1,LZ1,LELT)
     $ ,             H2 (LX1,LY1,LZ1,LELT)
      COMMON /SCRHI/ H2INV (LX1,LY1,LZ1,LELV)
C     
      NTOT1  = NX1*NY1*NZ1*NELV
C     
      IF (IFAA) THEN
         NTOT1  = NX1*NY1*NZ1*NELV
         CALL RONE    (H1,NTOT1)
         CALL RZERO   (H2,NTOT1)
         CALL ALPHAM1 (EIGAA1,V1MASK,VMULT,H1,H2,1)
         CALL ALPHAM1 (EIGAA2,V2MASK,VMULT,H1,H2,2)
         EIGAA = MIN  (EIGAA1,EIGAA2)
         IF (NDIM.EQ.3) THEN
            CALL ALPHAM1 (EIGAA3,V3MASK,VMULT,H1,H2,3)
            EIGAA = MIN  (EIGAA,EIGAA3)
         ENDIF
         IF (NID.EQ.0 .AND. ISTEP.LE.0) WRITE (6,*) 'EIGAA = ',EIGAA
      ENDIF
C     
      IF (IFAS) THEN
         INLOC = 0
         CALL RONE    (H1,NTOT1)
         CALL RZERO   (H2,NTOT1)
         CALL RZERO   (H2INV,NTOT1)
         CALL ALPHAM2  (EIGAS,H1,H2,H2INV,INLOC)
         IF (NID.EQ.0 .AND. ISTEP.LE.0) WRITE (6,*) 'EIGAS = ',EIGAS
      ENDIF
C     
      IF (IFAE) THEN
         INLOC = 1
         CALL RZERO   (H1,NTOT1)
         CALL RONE    (H2,NTOT1)
         CALL RONE    (H2INV,NTOT1)
         CALL ALPHAM2  (EIGAE,H1,H2,H2INV,INLOC)
         IF (NID.EQ.0 .AND. ISTEP.LE.0) WRITE (6,*) 'EIGAE = ',EIGAE
      ENDIF
C     
      IF (IFAST) THEN
         INLOC = -1
         CALL SETHLM  (H1,H2,INLOC)
         CALL INVERS2 (H2INV,H2,NTOT1)
         CALL ALPHAM2  (EIGAST,H1,H2,H2INV,INLOC)
         IF (NID.EQ.0 .AND. ISTEP.LE.0) WRITE (6,*) 'EIGAST = ',EIGAST
      ENDIF
C     
      IF (IFGS) THEN
         INLOC = 0
         CALL RONE    (H1,NTOT1)
         CALL RZERO   (H2,NTOT1)
         CALL RZERO   (H2INV,NTOT1)
         CALL GAMMAM2 (EIGGS,H1,H2,H2INV,INLOC)
         IF (NID.EQ.0 .AND. ISTEP.LE.0) WRITE (6,*) 'EIGGS = ',EIGGS
      ENDIF
C     
      IF (IFGE) THEN
         INLOC = 1
         CALL RZERO   (H1,NTOT1)
         CALL RONE    (H2,NTOT1)
         CALL RONE    (H2INV,NTOT1)
         CALL GAMMAM2 (EIGGE,H1,H2,H2INV,INLOC)
         IF (NID.EQ.0 .AND. ISTEP.LE.0) WRITE (6,*) 'EIGGE = ',EIGGE
      ENDIF
C     
      IF (IFGST) THEN
         INLOC = -1
         CALL SETHLM  (H1,H2,INLOC)
         CALL INVERS2 (H2INV,H2,NTOT1)
         CALL GAMMAM2 (EIGGST,H1,H2,H2INV,INLOC)
         IF (NID.EQ.0 .AND. ISTEP.LE.0) WRITE (6,*) 'EIGGST = ',EIGGST
      ENDIF
C     
      IF (IFGA) THEN
         NTOT1  = NX1*NY1*NZ1*NELV
         CALL RONE    (H1,NTOT1)
         CALL RZERO   (H2,NTOT1)
         IF (.NOT.IFSTRS) THEN
            CALL GAMMAM1 (EIGGA1,V1MASK,VMULT,H1,H2,1)
            CALL GAMMAM1 (EIGGA2,V2MASK,VMULT,H1,H2,2)
            EIGGA3 = 0.
            IF (NDIM.EQ.3)
     $      CALL GAMMAM1 (EIGGA3,V3MASK,VMULT,H1,H2,3)
            EIGGA = MAX  (EIGGA1,EIGGA2,EIGGA3)
         ELSE
            CALL GAMMASF (H1,H2)
         ENDIF
      ENDIF
C     
      RETURN
      END
C     
      SUBROUTINE ALPHAM1 (ALPHA,MASK,MULT,H1,H2,ISD)
C-----------------------------------------------------------------------
C     
C     Compute minimum eigenvalue, ALPHA, of the discrete Helmholtz opera
C     
C-----------------------------------------------------------------------

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
# 210 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 210 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 211 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 211 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 212 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 212 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C     
      REAL            MASK (LX1,LY1,LZ1,1)
      REAL            MULT (LX1,LY1,LZ1,1)
      REAL            H1   (LX1,LY1,LZ1,1)
      REAL            H2   (LX1,LY1,LZ1,1)
      COMMON /SCREV/  X1   (LX1,LY1,LZ1,LELT)
     $ ,              Y1   (LX1,LY1,LZ1,LELT)
      CHARACTER NAME*4
C     
      IF (IMESH.EQ.1) NEL  = NELV
      IF (IMESH.EQ.2) NEL  = NELT
      IF (ISD  .EQ.1) NAME = 'EVVX'
      IF (ISD  .EQ.2) NAME = 'EVVX'
      IF (ISD  .EQ.3) NAME = 'EVVX'
C     
      NXYZ1  = NX1*NY1*NZ1
      NTOT1  = NXYZ1*NEL
      EVNEW  = 0.
      CALL STARTX1 (X1,Y1,MASK,MULT,NEL)
C     
      DO 1000 ITER=1,NMXE
         CALL AXHELM (Y1,X1,H1,H2,IMESH,ISD)
         CALL COL2   (Y1,MASK,NTOT1)
         CALL DSSUM  (Y1,NX1,NY1,NZ1)
         RQ     = GLSC3 (X1,Y1,MULT,NTOT1)
         EVOLD  = EVNEW
         EVNEW  = RQ
         write (6,*) 'alphaa = ',rq
         CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
         IF (CRIT.LT.TOLEV)                  GOTO 2000
         CALL COL2    (X1,BM1,NTOT1)
         CALL HMHOLTZ ('NOMG',Y1,X1,H1,H2,MASK,MULT,
     $                              IMESH,TOLHE,NMXH,ISD)
         CALL COL3    (X1,BM1,Y1,NTOT1)
         CALL DSSUM   (X1,NX1,NY1,NZ1)
         YY = GLSC3  (X1,Y1,MULT,NTOT1)
         YNORM = 1./SQRT(YY)
         CALL CMULT   (Y1,YNORM,NTOT1)
         CALL COPY    (X1,Y1,NTOT1)
 1000 CONTINUE
 2000 CONTINUE
C     
      ALPHA = RQ
      RETURN
      END
C     
      SUBROUTINE GAMMAM1 (GAMMA,MASK,MULT,H1,H2,ISD)
C-----------------------------------------------------------------------
C     
C     Compute maximum eigenvalue of the discrete Helmholtz operator
C     
C-----------------------------------------------------------------------

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
# 265 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 265 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 266 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 266 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 267 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 267 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C     
      REAL            MASK (LX1,LY1,LZ1,1)
      REAL            MULT (LX1,LY1,LZ1,1)
      REAL            H1   (LX1,LY1,LZ1,1)
      REAL            H2   (LX1,LY1,LZ1,1)
      COMMON /SCREV/  X1   (LX1,LY1,LZ1,LELT)
     $ ,              Y1   (LX1,LY1,LZ1,LELT)
C     
      IF (IMESH.EQ.1) NEL = NELV
      IF (IMESH.EQ.2) NEL = NELT
      NXYZ1  = NX1*NY1*NZ1
      NTOT1  = NXYZ1*NEL
      EVNEW  = 0.
c     pff (2/15/96)
      if (isd.eq.1) CALL STARTX1 (X1,Y1,MASK,MULT,NEL)
C     
      DO 1000 ITER=1,NMXE
         CALL AXHELM (Y1,X1,H1,H2,IMESH,ISD)
         CALL COL2   (Y1,MASK,NTOT1)
         CALL DSSUM  (Y1,NX1,NY1,NZ1)
         RQ     = GLSC3 (X1,Y1,MULT,NTOT1)
         EVOLD  = EVNEW
         EVNEW  = RQ
         CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
C     
C HMT removed
C     
C         if (nid.eq.0) then
C            write(6,*) iter,' eig_max A:',evnew,crit,tolev
C         endif
         IF (CRIT.LT.TOLEV)                  GOTO 2000
         CALL COL3 (X1,BINVM1,Y1,NTOT1)
         XX     = GLSC3 (X1,Y1,MULT,NTOT1)
         XNORM  = 1./SQRT(XX)
         CALL CMULT (X1,XNORM,NTOT1)
 1000 CONTINUE
 2000 CONTINUE
C     
      GAMMA = RQ
      RETURN
      END
C     
      SUBROUTINE ALPHAM2 (ALPHA,H1,H2,H2INV,INLOC)
C----------------------------------------------------------------------
C     
C     Compute minimum eigenvalue, ALPHA, of one of the matrices
C     defined on the pressure mesh:
C     INLOC =  0  : DA-1DT
C     INLOC =  1  : DB-1DT
C     INLOC = -1  : D(A+B/DT)-1DT
C     
C----------------------------------------------------------------------

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
# 320 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 320 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 321 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 321 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 322 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 322 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C     
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           H2INV(LX1,LY1,LZ1,1)
      COMMON /SCREV/ X2   (LX2,LY2,LZ2,LELV)
     $ ,             Y2   (LX2,LY2,LZ2,LELV)
C     
      NTOT2  = NX2*NY2*NZ2*NELV
      EVNEW  = 0.
      CALL STARTX2 (X2,Y2)
C     
      DO 1000 ITER=1,NMXE
         CALL CDABDTP (Y2,X2,H1,H2,H2INV,INLOC)
         RQ = GLSC2  (X2,Y2,NTOT2)
         EVOLD  = EVNEW
         EVNEW  = RQ
c        write (6,*) 'new eigenvalue ************* eigas = ',evnew
         CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
         IF (CRIT.LT.TOLEV)                GOTO 2000
         CALL COL2   (X2,BM2,NTOT2)
         CALL UZAWA  (X2,H1,H2,H2INV,INLOC,ICG)
         CALL COL3   (Y2,BM2,X2,NTOT2)
         XX = GLSC2 (X2,Y2,NTOT2)
         XNORM = 1./SQRT(XX)
         CALL CMULT  (X2,XNORM,NTOT2)
 1000 CONTINUE
 2000 CONTINUE
C     
      ALPHA = RQ
      RETURN
      END
C     
      SUBROUTINE GAMMAM2 (GAMMA,H1,H2,H2INV,INLOC)
C-------------------------------------------------------------------
C     
C     Compute maximum eigenvalue, GAMMA, of one of the matrices
C     defined on the pressure mesh:
C     INLOC =  0  : DA-1DT
C     INLOC =  1  : DB-1DT
C     INLOC = -1  : D(A+B/DT)-1DT
C     
C-------------------------------------------------------------------

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
# 365 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 365 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 366 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 366 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

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
# 367 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 367 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C     
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      REAL           H2INV (LX1,LY1,LZ1,1)
      COMMON /SCREV/ X2 (LX2,LY2,LZ2,LELV)
     $ ,             Y2 (LX2,LY2,LZ2,LELV)
C     
      NTOT2  = NX2*NY2*NZ2*NELV
      EVNEW  = 0.
      CALL STARTX2 (X2,Y2)
C     
      DO 1000 ITER=1,NMXE
         CALL CDABDTP (Y2,X2,H1,H2,H2INV,INLOC)
         RQ = GLSC2  (X2,Y2,NTOT2)
         EVOLD  = EVNEW
         EVNEW  = RQ
         CRIT   = ABS((EVNEW-EVOLD)/EVNEW)
         IF (CRIT.LT.TOLEV)                GOTO 2000
         CALL INVCOL3 (X2,Y2,BM2,NTOT2)
         XX     = GLSC2  (Y2,X2,NTOT2)
         XNORM  = 1./SQRT(XX)
         CALL CMULT   (X2,XNORM,NTOT2)
 1000 CONTINUE
 2000 CONTINUE
C     
      GAMMA = RQ
      RETURN
      END
C     
c-----------------------------------------------------------------------
      SUBROUTINE STARTX1 (X1,Y1,MASK,MULT,NEL)
      
c     Compute startvector for finding an eigenvalue on mesh 1.
c     Normalization: XT*B*X = 1
      

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
# 403 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 403 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 404 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 404 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
      
      REAL X1   (LX1,LY1,LZ1,1)
      REAL Y1   (LX1,LY1,LZ1,1)
      REAL MASK (LX1,LY1,LZ1,1)
      REAL MULT (LX1,LY1,LZ1,1)
      
      NTOT1 = NX1*NY1*NZ1*NEL
      CALL COPY       (X1,BM1,NTOT1)
      
      
      call rand_fld_h1(y1)            ! pff 3/21/12
      small = 0.001*glamax(x1,ntot1)
      call add2s2(x1,y1,small,ntot1)
      
      
      CALL COL2       (X1,MASK,NTOT1)
      CALL COL3       (Y1,BM1,X1,NTOT1)
      CALL DSSUM      (Y1,NX1,NY1,NZ1)
      XX     = GLSC3 (X1,Y1,MULT,NTOT1)
      XNORM  = 1./SQRT(XX)
      CALL CMULT      (X1,XNORM,NTOT1)
      
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE STARTX2 (X2,Y2)
C------------------------------------------------------------------
C     
C     Compute startvector for finding an eigenvalue on mesh 2.
C     
C------------------------------------------------------------------

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
# 436 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 436 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 437 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f" 2
# 437 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/eigsolv.f"
C     
      REAL X2 (LX2,LY2,LZ2,LELV)
      REAL Y2 (LX2,LY2,LZ2,LELV)
C     
      NXYZ2  = NX2*NY2*NZ2
      NTOT2  = NXYZ2*NELV
      ICONST = 0
      IF ((NDIM .EQ. 2).AND.(NXYZ2 .EQ. 4)) ICONST = 1
      IF ((NDIM .EQ. 3).AND.(NXYZ2 .EQ. 8)) ICONST = 1
C     
      IF (ICONST .EQ. 1) THEN
         DO 1000 IEL=1,NELV
         DO 1000 K=1,NZ2
         DO 1000 J=1,NY2
         DO 1000 I=1,NX2
            X2(I,J,K,IEL) = I*J*K
 1000    CONTINUE
      ELSE
         CALL COPY (X2,BM2,NTOT2)
      ENDIF
C     
      call ortho (x2)
      CALL COL3 (Y2,BM2,X2,NTOT2)
      XX     = GLSC2 (X2,Y2,NTOT2)
      XNORM  = 1./SQRT(XX)
      CALL CMULT (X2,XNORM,NTOT2)
C     
      RETURN
      END
