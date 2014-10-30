# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
c-----------------------------------------------------------------------
      subroutine cbcmesh
C     
C     Generate boundary conditions (CBC arrays) for mesh solver
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
# 7 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 7 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 8 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 8 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 9 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 9 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 10 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 10 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      CHARACTER CBM*1,CBF*3,CBT*3,CB*3
C     
      IFLD   = 0
      NFACE  = 2*NDIM
      IFMELT = .FALSE.
      IF (IFTMSH(IFLD)) IFMELT=.TRUE.
C     
      DO 100 IEL=1,NELT
      DO 100 IFC=1,NFACE
         CBM = CBC(IFC,IEL,0)
         CBF = CBC(IFC,IEL,1)
         CBT = CBC(IFC,IEL,2)
C     
         IF (CBT(1:1).EQ.'M') THEN
             IFLD = 2
             CB   = CBT
             GOTO 200
         ENDIF
         IF (CBF(1:1).EQ.'M' .OR. CBF(1:1).EQ.'m') THEN
             IFLD = 1
             CB   = CBF
C             IF (CBF.EQ.'mv ' .AND. CBT.EQ.'E  ') THEN
C                 IFTMSH(0)=.TRUE.
C             ENDIF
             IF (CBF.EQ.'mv ' .AND. CBM.EQ.'+'  ) THEN
                CB = 'mvn'
                CBC(IFC,IEL,1) = CB
             ENDIF
             GOTO 200
         ENDIF
         IF (CBF(1:1).EQ.'V' .OR. CBF(1:1).EQ.'v' .OR.
     $       CBF(1:1).EQ.'W' ) THEN
             IFLD = 1
             CB   = 'FIX'
             IF (IFMELT .OR. CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IF (CBT.EQ.'T  ' .OR. CBT.EQ.'t  ') THEN
             IFLD = 2
             CB   = 'FIX'
             IF (CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IF (CBF.EQ.'P  ' .OR. CBF.EQ.'E  ') THEN
             IFLD = 1
             CB   = CBF
             IF (CBM.EQ.'-') CB='FIX'
             IF (CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IF (CBT.EQ.'P  ' .OR. CBT.EQ.'E  ') THEN
             IFLD = 2
             CB   = CBT
             IF (CBM.EQ.'-') CB='FIX'
             IF (CBM.EQ.'+') CB='SYM'
             GOTO 200
         ENDIF
         IFLD = 1
         IF (CBF.EQ.'   ') IFLD = 2
         CB   = 'SYM'
         IF (CBM.EQ.'-') CB = 'FIX'
C     
  200    CBC(IFC,IEL,0) = CB
         DO 250 I=1,5
  250    BC(I,IFC,IEL,0)=BC(I,IFC,IEL,IFLD)
  100 CONTINUE
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine admeshv
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
# 84 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 84 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 85 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 85 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 86 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 86 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      COMMON /SCRUZ/ FM1(LX1,LY1,LZ1,LELT)
     $             , FM2(LX1,LY1,LZ1,LELT)
     $             , FM3(LX1,LY1,LZ1,LELT)
     $             , PHI(LX1,LY1,LZ1,LELT)
C     
      NTOT1=NX1*NY1*NZ1*NELV
C     
      CALL RZERO (FM1,NTOT1)
      CALL RZERO (FM2,NTOT1)
      CALL RZERO (FM3,NTOT1)
C     
      CALL DIVWS (FM1,VX,PHI,NELV,1)
      CALL DIVWS (FM2,VY,PHI,NELV,2)
      CALL ADD2  (BFX,FM1,NTOT1)
      CALL ADD2  (BFY,FM2,NTOT1)
      IF (NDIM.EQ.3) THEN
         CALL DIVWS (FM3,VZ,PHI,NELV,3)
         CALL ADD2  (BFZ,FM3,NTOT1)
      ENDIF
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine admesht
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
# 113 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 113 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 114 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 114 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 115 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 115 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      COMMON /SCRUZ/ FMT(LX1,LY1,LZ1,LELT)
     $             , PHI(LX1,LY1,LZ1,LELT)
C     
      IFLD = 0
      NEL  = NELFLD(IFLD)
      NTOT1= NX1*NY1*NZ1*NEL
C     
      CALL RZERO   (FMT,NTOT1)
      CALL DIVWS   (FMT,T(1,1,1,1,IFIELD-1),PHI,NEL,1)
      CALL ADDCOL3 (BQ(1,1,1,1,IFIELD-1),FMT,VTRANS(1,1,1,1,IFIELD),
     &              NTOT1)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine divws (fms,sfv,phi,nel,idir)
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
# 134 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 134 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 135 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 135 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 136 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 136 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 137 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 137 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 138 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 138 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 139 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 139 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      COMMON /SCRSF/ PHR(LX1,LY1,LZ1,LELT)
     $             , PHS(LX1,LY1,LZ1,LELT)
     $             , PHT(LX1,LY1,LZ1,LELT)
      
      DIMENSION FMS(LX1,LY1,LZ1,1)
     $        , SFV(LX1,LY1,LZ1,1)
     $        , PHI(LX1,LY1,LZ1,1)
C     
      NXYZ1 = NX1*NY1*NZ1
      NTOT1 = NXYZ1*NEL
C     
      CALL COL3    (PHI,SFV,WX,NTOT1)
      CALL URST    (PHI,PHR,PHS,PHT,NEL)
      CALL ADDCOL3 (FMS,RXM1,PHR,NTOT1)
      CALL ADDCOL3 (FMS,SXM1,PHS,NTOT1)
      IF (NDIM.EQ.3) CALL ADDCOL3 (FMS,TXM1,PHT,NTOT1)
C     
      CALL COL3    (PHI,SFV,WY,NTOT1)
      CALL URST    (PHI,PHR,PHS,PHT,NEL)
      CALL ADDCOL3 (FMS,RYM1,PHR,NTOT1)
      CALL ADDCOL3 (FMS,SYM1,PHS,NTOT1)
      IF (NDIM.EQ.3) CALL ADDCOL3 (FMS,TYM1,PHT,NTOT1)
C      
      IF (NDIM.EQ.3) THEN
         CALL COL3    (PHI,SFV,WZ,NTOT1)
         CALL URST    (PHI,PHR,PHS,PHT,NEL)
         CALL ADDCOL3 (FMS,RZM1,PHR,NTOT1)
         CALL ADDCOL3 (FMS,SZM1,PHS,NTOT1)
         CALL ADDCOL3 (FMS,TZM1,PHT,NTOT1)
      ENDIF
C     
      CALL COL2    (FMS,BM1,NTOT1)
      CALL INVCOL2 (FMS,JACM1,NTOT1)
C     
      IF (IFAXIS) CALL AXIFMS (FMS,SFV,PHI,NEL,IDIR)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine axifms (fms,sfv,phi,nel,idir)
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
# 182 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 182 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 183 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 183 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 184 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 184 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 185 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 185 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 186 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 186 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 187 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 187 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 188 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 188 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRSF/ PHR(LX1,LY1,LZ1,LELT)
     $             , PHS(LX1,LY1,LZ1,LELT)
     $             , PHT(LX1,LY1,LZ1,LELT)
C     
      DIMENSION FMS(LX1,LY1,LZ1,1)
     $        , PHI(LX1,LY1,LZ1,1)
     $        , SFV(LX1,LY1,LZ1,1)
     $        , WYS(LX1)
      EQUIVALENCE (WYS(1),PHT(1,1,1,1))
C     
      NXYZ1 = NX1*NY1*NZ1
      NTOT1 = NXYZ1*NEL
      CALL COL3 (PHI,SFV,WY,NTOT1)
C     
      DO 100 IEL=1,NEL
      IF ( IFRZER(IEL) ) THEN
         IF (IDIR.EQ.1) THEN
            CALL MXM (WY(1,1,1,IEL),NX1,DATM1,NY1,WYS,1)
            DO 220 IX=1,NX1
            FMS(IX,1,1,IEL)= FMS(IX,1,1,IEL) + WXM1(IX)*WAM1(1)*
     $                       WYS(IX)*SFV(IX,1,1,IEL)*JACM1(IX,1,1,IEL)
  220       CONTINUE
         ENDIF
         DO 320 IX=1,NX1
         DO 320 IY=2,NY1
            FMS(IX,IY,1,IEL)=FMS(IX,IY,1,IEL) + PHI(IX,IY,1,IEL) *
     $                       BM1(IX,IY,1,IEL) / YM1(IX,IY,1,IEL)
  320    CONTINUE
      ELSE
         CALL ADDCOL4 (FMS(1,1,1,IEL),PHI(1,1,1,IEL),JACM1(1,1,1,IEL),
     $                 W2CM1,NXYZ1)
      ENDIF
  100 CONTINUE
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine updcoor
C-----------------------------------------------------------------------
C     
C     Subroutine to update geometry for moving boundary problems
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
# 232 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 232 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 233 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 233 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      IFIELD = 0
      NEL    = NELFLD(IFIELD)
C     
C     Update collocation points coordinates
C     
      CALL UPDXYZ (NEL)
C     
C     Shift lagged mesh velocity
C     
      CALL LAGMSHV (NEL)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine mvbdry (nel)
C     
C     Routine to evaluate mesh velocities at all moving boundaries
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
# 253 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 253 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 254 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 254 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 255 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 255 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 256 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 256 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 257 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 257 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 258 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 258 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRSF/ WVX(LX1,LY1,LZ1,LELT)
     $             , WVY(LX1,LY1,LZ1,LELT)
     $             , WVZ(LX1,LY1,LZ1,LELT)
      COMMON /SCRCH/ WTX(LX1,LY1,LZ1,LELT)
     $             , WTY(LX1,LY1,LZ1,LELT)
      COMMON /SCRMG/ WTZ(LX1,LY1,LZ1,LELT)
     $             , RNX(LX1,LY1,LZ1,LELT)
     $             , RNY(LX1,LY1,LZ1,LELT)
     $             , RNZ(LX1,LY1,LZ1,LELT)
      COMMON /SCRUZ/ DSA(LX1,LY1,LZ1,LELT)
     $             , QNI(LX1,LY1,LZ1,LELT)
     $             , SMT(LX1,LY1,LZ1,LELT)
     $             , TA (LX1,LY1,LZ1,LELT)
C     
      LOGICAL IFALGN,IFNORX,IFNORY,IFNORZ,IFDSMV,IFREGW
      CHARACTER CB*3
C     
      IFIELD = 0
      NXYZ1  = NX1*NY1*NZ1
      NTOT1  = NX1*NY1*NZ1*NEL
      NFACE  = 2*NDIM
      CALL RZERO3  (RNX,RNY,RNZ,NTOT1)
C     
      DO 100 IEL=1,NEL
      DO 100 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MS ' .OR. CB.EQ.'ms ' .OR. 
     $       CB.EQ.'MSI' .OR. CB.EQ.'msi' .OR. 
     $       CB.EQ.'MM ' .OR. CB.EQ.'mm ' .OR. 
     $       CB.EQ.'mv ' .OR. CB.EQ.'mvn' .OR.
     $       CB.EQ.'MLI') THEN
             CALL FACEXV (UNX(1,1,IFC,IEL),UNY(1,1,IFC,IEL),
     $                    UNZ(1,1,IFC,IEL),RNX(1,1,1,IEL),
     $                    RNY(1,1,1,IEL),RNZ(1,1,1,IEL),IFC,1)
         ENDIF
  100 CONTINUE
C     
      CALL DSSUM (RNX,NX1,NY1,NZ1)
      CALL DSSUM (RNY,NX1,NY1,NZ1)
      IF (NDIM.EQ.3) CALL DSSUM (RNZ,NX1,NY1,NZ1)
      CALL UNITVEC (RNX,RNY,RNZ,NTOT1)
C     
      CALL RZERO3 (WVX,WVY,WVZ,NTOT1)
      CALL RZERO3 (WTX,WTY,WTZ,NTOT1)
      DO 1000 ISWEEP=1,2
C     
      IFREGW = .FALSE.
      IFDSMV = .FALSE.
      CALL RZERO  (DSA,NTOT1)
      CALL RZERO  (TA,NTOT1)
C     
      IF (IFFLOW) THEN
C     
      IFIELD = 1
      CALL RZERO  (SMT,NTOT1)
      DO 210 IEL=1,NELV
      DO 210 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'mv ' .OR. CB.EQ.'mvn'  .OR.
     $       CB.EQ.'MM ' .OR. CB.EQ.'mm '  .OR.
     $       CB.EQ.'MS ' .OR. CB.EQ.'ms ') THEN
            IFREGW = .TRUE.
            CALL FACEC3 (WVX(1,1,1,IEL),WVY(1,1,1,IEL),WVZ(1,1,1,IEL),
     $                   VX(1,1,1,IEL),VY(1,1,1,IEL),VZ(1,1,1,IEL),IFC)
            IF (CB.NE.'mv ')
     $      CALL NORCMP (WVX(1,1,1,IEL),WVY(1,1,1,IEL),WVZ(1,1,1,IEL),
     $                   RNX(1,1,1,IEL),RNY(1,1,1,IEL),RNZ(1,1,1,IEL),
     $                 IFC) 
         ENDIF
         IF (CB.EQ.'MSI' .OR. CB.EQ.'msi') THEN
            IFDSMV = .TRUE.
            CALL FACSMT (SMT(1,1,1,IEL),IFC)
         ENDIF
  210 CONTINUE
      
      if (istep.eq.0) call opcopy(wx,wy,wz,wvx,wvy,wvz)
c     call outpost(wvx,wvy,wvz,wtx,wtx,'   ')
c     call outpost(wx,wy,wz,wtx,wtx,'   ')
c     write(6,*) 'quit1'
c     call exitt
      
      iregw = 0                   ! Global handshake on logicals ifregw,
      if (ifregw) iregw = 1       ! pff  9/11/07
      iregw = iglmax(iregw,1)
      if (iregw.eq.1) ifregw = .true.
      
      idsmv = 0
      if (ifdsmv) idsmv = 1
      idsmv = iglmax(idsmv,1)
      if (idsmv.eq.1) ifdsmv = .true.
      
C     
      IF (IFDSMV) THEN
      CALL DSSUM (SMT,NX1,NY1,NZ1)
      DO 215 IEL=1,NELV
      DO 215 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MSI' .OR. CB.EQ.'msi') THEN
         CALL FACEC3 (WTX(1,1,1,IEL),WTY(1,1,1,IEL),WTZ(1,1,1,IEL),
     $                VX(1,1,1,IEL),VY(1,1,1,IEL),VZ(1,1,1,IEL),IFC)
         CALL FACEMV (WTX(1,1,1,IEL),WTY(1,1,1,IEL),WTZ(1,1,1,IEL),
     $                RNX(1,1,1,IEL),RNY(1,1,1,IEL),RNZ(1,1,1,IEL),
     $                SMT(1,1,1,IEL),IFC) 
         ENDIF
  215 CONTINUE
      
         CALL DSSUM (WTX,NX1,NY1,NZ1)
         CALL DSSUM (WTY,NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSSUM (WTZ,NX1,NY1,NZ1)
      
      ENDIF
      
      ENDIF
C     
      IF (IFMELT .AND. ISTEP.GT.0) THEN
         IFIELD = 2
         CALL RZERO (SMT,NTOT1)
         CALL CQNET (QNI,TA,NEL)
      DO 220 IEL=1,NELT
      DO 220 IFC=1,NFACE
         CB   = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MLI') THEN
            CALL FACSMT (SMT(1,1,1,IEL),IFC)
         ENDIF
         IF (CB.EQ.'MLI' .OR. CB.EQ.'MCI') THEN
            CALL FACEXS (AREA(1,1,IFC,IEL),TA,IFC,1)
            CALL ADD2   (DSA(1,1,1,IEL),TA,NXYZ1)
         ENDIF
  220 CONTINUE
         CALL DSSUM (SMT,NX1,NY1,NZ1)
         CALL DSSUM (DSA,NX1,NY1,NZ1)
      DO 280 IEL=1,NELT
      DO 280 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'MLI') THEN
            RHOLA = -0.5 * BC(5,IFC,IEL,IFIELD)
            CALL FACEMT (WTX(1,1,1,IEL),WTY(1,1,1,IEL),WTZ(1,1,1,IEL),
     $                   RNX(1,1,1,IEL),RNY(1,1,1,IEL),RNZ(1,1,1,IEL),
     $                   QNI(1,1,1,IEL),DSA(1,1,1,IEL),SMT(1,1,1,IEL),
     $                   RHOLA,IFC)
         ENDIF
  280 CONTINUE
         CALL DSSUM (WTX,NX1,NY1,NZ1)
         CALL DSSUM (WTY,NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSSUM (WTZ,NX1,NY1,NZ1)
      ENDIF
C     
      IFIELD = 0
      DO 330 IEL=1,NEL
      DO 330 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'SYM') THEN
            CALL CHKNORD (IFALGN,IFNORX,IFNORY,IFNORZ,IFC,IEL)
            IF (IFREGW) THEN
               IF (IFNORX) CALL FACEV (WVX,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORY) CALL FACEV (WVY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORZ) CALL FACEV (WVZ,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (.NOT.IFALGN) 
     $         CALL FACZQN (WVX(1,1,1,IEL),WVY(1,1,1,IEL),
     $                      WVZ(1,1,1,IEL),IFC,IEL)
            ENDIF
            IF (IFDSMV .OR. IFMELT) THEN
               IF (IFNORX) CALL FACEV (WTX,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORY) CALL FACEV (WTY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (IFNORZ) CALL FACEV (WTZ,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (.NOT.IFALGN) 
     $         CALL FACZQN (WTX(1,1,1,IEL),WTY(1,1,1,IEL),
     $                      WTZ(1,1,1,IEL),IFC,IEL)
            ENDIF
         ENDIF
  330 CONTINUE
C     
      DO 350 IEL=1,NEL
      DO 350 IFC=1,NFACE
         CB = CBC(IFC,IEL,IFIELD)
         IF (CB.EQ.'FIX') THEN
            IF (IFREGW) THEN
               CALL FACEV (WVX,IEL,IFC,0.0,NX1,NY1,NZ1)
               CALL FACEV (WVY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (NDIM.EQ.3) CALL FACEV (WVZ,IEL,IFC,0.0,NX1,NY1,NZ1)
            ENDIF
            IF (IFDSMV .OR. IFMELT) THEN
               CALL FACEV (WTX,IEL,IFC,0.0,NX1,NY1,NZ1)
               CALL FACEV (WTY,IEL,IFC,0.0,NX1,NY1,NZ1)
               IF (NDIM.EQ.3) CALL FACEV (WTZ,IEL,IFC,0.0,NX1,NY1,NZ1)
            ENDIF
         ENDIF
  350 CONTINUE 
C     
      IF (ISWEEP.EQ.1) THEN
         IF (IFREGW) THEN
            CALL DSOP (WVX,'MXA',NX1,NY1,NZ1)
            CALL DSOP (WVY,'MXA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WVZ,'MXA',NX1,NY1,NZ1)
         ENDIF
         IF (IFDSMV .OR. IFMELT) THEN
            CALL DSOP (WTX,'MXA',NX1,NY1,NZ1)
            CALL DSOP (WTY,'MXA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WTZ,'MXA',NX1,NY1,NZ1)
         ENDIF
      ELSE
         IF (IFREGW) THEN
            CALL DSOP (WVX,'MNA',NX1,NY1,NZ1)
            CALL DSOP (WVY,'MNA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WVZ,'MNA',NX1,NY1,NZ1)
         ENDIF
         IF (IFDSMV .OR. IFMELT) THEN
            CALL DSOP (WTX,'MNA',NX1,NY1,NZ1)
            CALL DSOP (WTY,'MNA',NX1,NY1,NZ1)
            IF (NDIM.EQ.3) CALL DSOP (WTZ,'MNA',NX1,NY1,NZ1)
         ENDIF
      ENDIF
C     
 1000 CONTINUE
C     
      CALL RMASK (WX,WY,WZ,NEL)
c     call outpost(wx,wy,wz,wtx,wtx,'   ')
c     write(6,*) 'quit2'
c     call exitt
C     
      IF (IFREGW) THEN
         CALL ADD2  (WX,WVX,NTOT1)
         CALL ADD2  (WY,WVY,NTOT1)
         IF (NDIM.EQ.3) CALL ADD2  (WZ,WVZ,NTOT1)
      ENDIF
      IF (IFDSMV .OR. IFMELT) THEN
         CALL ADD2  (WX,WTX,NTOT1)
         CALL ADD2  (WY,WTY,NTOT1)
         IF (NDIM.EQ.3) CALL ADD2  (WZ,WTZ,NTOT1)
      ENDIF
      
c     call outpost(wx,wy,wz,wtx,wtx,'   ')
c     call exitt
      
      return
      end
c-----------------------------------------------------------------------
      subroutine norcmp (wt1,wt2,wt3,rnx,rny,rnz,ifc)
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
# 498 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 498 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)
C     
      DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1)
     $        , RNX(LX1,LY1,LZ1),RNY(LX1,LY1,LZ1),RNZ(LX1,LY1,LZ1)
C     
      NXYZ1 = NX1*NY1*NZ1
C     
      CALL COPY (R1,WT1,NXYZ1)
      CALL COPY (R2,WT2,NXYZ1)
      IF (NDIM.EQ.3) CALL COPY (R3,WT3,NXYZ1)
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C     
      IF (NDIM.EQ.2) THEN
         DO 200 J2=JS2,JF2,JSKIP2
         DO 200 J1=JS1,JF1,JSKIP1
            WN          = R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                    R2(J1,J2,1)*RNY(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
  200    CONTINUE
      ELSE
         DO 300 J2=JS2,JF2,JSKIP2
         DO 300 J1=JS1,JF1,JSKIP1
            WN          = R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                    R2(J1,J2,1)*RNY(J1,J2,1) +
     $                    R3(J1,J2,1)*RNZ(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
            WT3(J1,J2,1) = WN *RNZ(J1,J2,1)
  300    CONTINUE
      ENDIF
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine facemv (wt1,wt2,wt3,rnx,rny,rnz,smt,ifc)
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
# 536 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 536 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)
C     
      DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1)
     $        , RNX(LX1,LY1,LZ1),RNY(LX1,LY1,LZ1),RNZ(LX1,LY1,LZ1)
     $        , SMT(LX1,LY1,LZ1)
C     
      NXYZ1 = NX1*NY1*NZ1
C     
      CALL COPY (R1,WT1,NXYZ1)
      CALL COPY (R2,WT2,NXYZ1)
      IF (NDIM.EQ.3) CALL COPY (R3,WT3,NXYZ1)
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C     
      IF (NDIM.EQ.2) THEN
         DO 200 J2=JS2,JF2,JSKIP2
         DO 200 J1=JS1,JF1,JSKIP1
            WN          = ( R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                      R2(J1,J2,1)*RNY(J1,J2,1) ) / SMT(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
  200    CONTINUE
      ELSE
         DO 300 J2=JS2,JF2,JSKIP2
         DO 300 J1=JS1,JF1,JSKIP1
            WN          = ( R1(J1,J2,1)*RNX(J1,J2,1) + 
     $                      R2(J1,J2,1)*RNY(J1,J2,1) +
     $                      R3(J1,J2,1)*RNZ(J1,J2,1) ) / SMT(J1,J2,1)
            WT1(J1,J2,1) = WN *RNX(J1,J2,1)
            WT2(J1,J2,1) = WN *RNY(J1,J2,1)
            WT3(J1,J2,1) = WN *RNZ(J1,J2,1)
  300    CONTINUE
      ENDIF
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine faczqn (wt1,wt2,wt3,ifc,iel)
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
# 575 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 575 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 576 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 576 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 577 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 577 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRUZ/ R1(LX1,LY1,LZ1),R2(LX1,LY1,LZ1),R3(LX1,LY1,LZ1)
C     
      DIMENSION WT1(LX1,LY1,LZ1),WT2(LX1,LY1,LZ1),WT3(LX1,LY1,LZ1)
C     
      NXYZ1 = NX1*NY1*NZ1
      CALL COPY (R1,WT1,NXYZ1)
      CALL COPY (R2,WT2,NXYZ1)
      IF (NDIM.EQ.3) CALL COPY (R3,WT3,NXYZ1)
C     
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
      I = 0
C     
      IF (NDIM.EQ.2) THEN
      DO 200 J2=JS2,JF2,JSKIP2
      DO 200 J1=JS1,JF1,JSKIP1
         I = I+1
         W1           = R1(J1,J2,1)*T1X(I,1,IFC,IEL) + 
     $                  R2(J1,J2,1)*T1Y(I,1,IFC,IEL)
         WT1(J1,J2,1) = W1 *T1X(I,1,IFC,IEL)
         WT2(J1,J2,1) = W1 *T1Y(I,1,IFC,IEL)
  200 CONTINUE
      ELSE
      DO 300 J2=JS2,JF2,JSKIP2
      DO 300 J1=JS1,JF1,JSKIP1
         I = I+1
         W1           = R1(J1,J2,1)*T1X(I,1,IFC,IEL) + 
     $                  R2(J1,J2,1)*T1Y(I,1,IFC,IEL) +
     $                  R3(J1,J2,1)*T1Z(I,1,IFC,IEL)
         WT1(J1,J2,1) = W1 *T1X(I,1,IFC,IEL)
         WT2(J1,J2,1) = W1 *T1Y(I,1,IFC,IEL)
         WT3(J1,J2,1) = W1 *T1Z(I,1,IFC,IEL)
  300 CONTINUE
      ENDIF
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine facsmt (smt,ifc)
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
# 617 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 617 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      DIMENSION SMT(LX1,LY1,LZ1)
C     
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C     
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         SMT(J1,J2,1)=SMT(J1,J2,1) + 1.0
  100 CONTINUE
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine cqnet (qni,ta,nel)
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
# 632 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 632 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 633 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 633 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 634 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 634 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 635 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 635 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT)
     $             , H2(LX1,LY1,LZ1,LELT)
C     
      DIMENSION QNI(LX1,LY1,LZ1,1)
     $        , TA (LX1,LY1,LZ1,1)
C     
      INTLOC = -1
      IMSHL  =  2
      NTOT1  = NX1*NY1*NZ1*NEL
C     
      CALL SETHLM (H1,H2,INTLOC)
      CALL AXHELM (TA,T(1,1,1,1,IFIELD-1),H1,H2,IMSHL,1)
      CALL SUB3   (QNI,TA,BQ(1,1,1,1,IFIELD-1),NTOT1)
      CALL DSSUM  (QNI,NX1,NY1,NZ1)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine facemt (w1,w2,w3,rnx,rny,rnz,qni,dsa,smt,rhola,ifc)
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
# 656 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 656 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 657 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 657 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      DIMENSION  W1 (LX1,LY1,LZ1)
     $        ,  W2 (LX1,LY1,LZ1)
     $        ,  W3 (LX1,LY1,LZ1)
     $        ,  RNX(LX1,LY1,LZ1)
     $        ,  RNY(LX1,LY1,LZ1)
     $        ,  RNZ(LX1,LY1,LZ1)
     $        ,  QNI(LX1,LY1,LZ1)
     $        ,  DSA(LX1,LY1,LZ1)
     $        ,  SMT(LX1,LY1,LZ1)
C     
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C     
      IF (NDIM.EQ.2) THEN
      DO 200 J2=JS2,JF2,JSKIP2
      DO 200 J1=JS1,JF1,JSKIP1
         AA = QNI(J1,J2,1) / ( DSA(J1,J2,1)*SMT(J1,J2,1)*RHOLA )
         W1(J1,J2,1) = RNX(J1,J2,1) * AA
         W2(J1,J2,1) = RNY(J1,J2,1) * AA
  200 CONTINUE
      ELSE
      DO 300 J2=JS2,JF2,JSKIP2
      DO 300 J1=JS1,JF1,JSKIP1
         AA = QNI(J1,J2,1) / ( DSA(J1,J2,1)*SMT(J1,J2,1)*RHOLA )
         W1(J1,J2,1) = RNX(J1,J2,1) * AA
         W2(J1,J2,1) = RNY(J1,J2,1) * AA
         W3(J1,J2,1) = RNZ(J1,J2,1) * AA
  300 CONTINUE
      ENDIF
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine elasolv (nel)
C     
C     Elastostatic solver for mesh deformation
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
# 695 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 695 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 696 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 696 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 697 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 697 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 698 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 698 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 699 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 699 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 700 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 700 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      COMMON /SCRNS/ DW1  (LX1,LY1,LZ1,LELT)
     $             , DW2  (LX1,LY1,LZ1,LELT)
     $             , DW3  (LX1,LY1,LZ1,LELT)
     $             , AW1  (LX1,LY1,LZ1,LELT)
     $             , AW2  (LX1,LY1,LZ1,LELT)
     $             , AW3  (LX1,LY1,LZ1,LELT)
      COMMON /SCRVH/ H1   (LX1,LY1,LZ1,LELT)
     $             , H2   (LX1,LY1,LZ1,LELT)
      common /scruz/ prt  (lx1,ly1,lz1,lelt)
      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV
C     
C     Set up parameters for mesh solver
C     
c     return
c     if (istep.gt.1) then
c        ifldx = ifield
c        ifield = 1
c        call incomprn (wx,wy,wz,prt) ! project U onto div-free space
c        ifield = ifldx
c        return
c     endif
c     
c     call quickmv
      
      if (ifusermv) return  ! Compute wx,wy,wz in userchk.
c     
      NTOT1  = NX1*NY1*NZ1*NEL
      MAXIT  = 1000
      MATMOD = -1
      IFH2   = .FALSE.
      IFSOLV = .TRUE.
      IMSOLV = 0
      
      vnu    = param(47)
      if (vnu.eq.0) vnu = 0.4 ! Default value
      vnu    = max(0.00,vnu)
      vnu    = min(0.49,vnu)
      
C     Set up elastic material constants
      
      CE = 1./(1. + VNU)
      C2 = VNU * CE / (1. - 2.*VNU)
      C3 = 0.5 * CE
      CALL CFILL (H1,C2,NTOT1)
      CALL CFILL (H2,C3,NTOT1)
C     
C     Solve for interior mesh velocities
C     
      CALL MESHTOL (AW1,TOLMSH,NEL,IMSOLV)
      IF (IMSOLV.EQ.1) return
C     
      CALL AXHMSF  (AW1,AW2,AW3,WX,WY,WZ,H1,H2,MATMOD)
c     call outpost(wx,wy,wz,h1,h2,'   ')
c     call exitt
      CALL CHSIGN  (AW1,NTOT1)
      CALL CHSIGN  (AW2,NTOT1)
      IF (NDIM.EQ.3) CALL CHSIGN (AW3,NTOT1)
      CALL HMHZSF  ('NOMG',DW1,DW2,DW3,AW1,AW2,AW3,H1,H2,
     $               W1MASK,W2MASK,W3MASK,WMULT,TOLMSH,
     $               MAXIT,MATMOD)
C     
C     Update mesh velocities
C     
      CALL ADD2 (WX,DW1,NTOT1)
      CALL ADD2 (WY,DW2,NTOT1)
      IF (NDIM.EQ.3) CALL ADD2 (WZ,DW3,NTOT1)
      
      ifldt = ifield
      ifield=1
      if (ifheat) ifield=2
      call dsavg(wx)
      call dsavg(wy)
      call dsavg(wz)
      ifield = ifldt
      
c     if (istep.gt.1) then
c        ifldx = ifield
c        ifield = 1
c        call incomprn (wx,wy,wz,prt) ! project U onto div-free space
c        ifield = ifldx
c        return
c     endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine meshtol (ta,tolmsh,nel,imsolv)
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
# 791 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 791 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/EIGEN" 1
C     
C     Eigenvalues
C     
# 4
      COMMON /EIGVAL/ EIGAA, EIGAS, EIGAST, EIGAE
     $               ,EIGGA, EIGGS, EIGGST, EIGGE
      COMMON /IFEIG / IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      LOGICAL         IFAA,IFAE,IFAS,IFAST,IFGA,IFGE,IFGS,IFGST
      
# 792 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 792 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 793 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 793 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 794 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 794 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      DIMENSION TA(LX1,LY1,LZ1,1)
C     
      NTOT1 = NX1*NY1*NZ1*NEL
      TOLAB = TOLREL
C     
      DELTA  = 1.0E-9
      X      = 1.0 + DELTA
      Y      = 1.0
      DIFF   = ABS(X - Y)
      IF (DIFF .EQ. 0.0) EPS = 1.0E-05
      IF (DIFF .GT. 0.0) EPS = 1.0E-12
C     
      CALL OPDOT  (TA,WX,WY,WZ,WX,WY,WZ,NTOT1)
C     
      WDOT = GLMAX(TA,NTOT1)
      WMAX = SQRT(WDOT)
      IF (WMAX .LT. EPS) THEN
          IMSOLV = 1
          return
      ENDIF
C     
      TOLMSH = TOLAB * WMAX * SQRT(EIGAA)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine updxyz (nel)
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
# 823 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 823 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 824 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 824 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 825 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 825 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 826 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 826 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRSF/ UX(LX1,LY1,LZ1,LELT)
     $             , UY(LX1,LY1,LZ1,LELT)
     $             , UZ(LX1,LY1,LZ1,LELT)
      DIMENSION ABM(3)
C     
      NTOT1 = NX1*NY1*NZ1*NEL
C     
      DO 10 I=1,NBD
   10 ABM(I) = DT*ABMSH(I)
C     
      IF (ISTEP.EQ.0) THEN
         CALL COPY (UX,WX,NTOT1)
         CALL COPY (UY,WY,NTOT1)
         IF (NDIM.EQ.3) CALL COPY (UZ,WZ,NTOT1)
      ELSE
         CALL CMULT2 (UX,WX,ABM(1),NTOT1)
         CALL CMULT2 (UY,WY,ABM(1),NTOT1)
         IF (NDIM.EQ.3) CALL CMULT2 (UZ,WZ,ABM(1),NTOT1)
         DO 100 ILAG=2,NBD
            CALL ADD2S2 (UX,WXLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
            CALL ADD2S2 (UY,WYLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
            IF (NDIM.EQ.3) 
     $      CALL ADD2S2 (UZ,WZLAG(1,1,1,1,ILAG-1),ABM(ILAG),NTOT1)
  100    CONTINUE
      ENDIF
C     
      CALL ADD2 (XM1,UX,NTOT1)
      CALL ADD2 (YM1,UY,NTOT1)
      IF (NDIM.EQ.3) CALL ADD2 (ZM1,UZ,NTOT1)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine lagmshv (nel)
C-----------------------------------------------------------------------
C     
C     Keep old mesh velocity
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
# 866 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 866 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 867 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 867 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 868 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 868 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 869 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 869 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      NTOT1 = NX1*NY1*NZ1*NEL
C     
      DO 100 ILAG=NBDINP-1,2,-1
         CALL COPY (WXLAG(1,1,1,1,ILAG),WXLAG(1,1,1,1,ILAG-1),NTOT1)
         CALL COPY (WYLAG(1,1,1,1,ILAG),WYLAG(1,1,1,1,ILAG-1),NTOT1)
         IF (NDIM.EQ.3)
     $   CALL COPY (WZLAG(1,1,1,1,ILAG),WZLAG(1,1,1,1,ILAG-1),NTOT1)
 100  CONTINUE
C     
         CALL COPY (WXLAG(1,1,1,1,1),WX,NTOT1)
         CALL COPY (WYLAG(1,1,1,1,1),WY,NTOT1)
         IF (NDIM.EQ.3)
     $   CALL COPY (WZLAG(1,1,1,1,1),WZ,NTOT1)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine facec3 (a1,a2,a3,b1,b2,b3,ifc)
C     
C     Copy the face (IFC) of B1,B2,B3 to A1,A2,A3.
C     IFACE is the input in the pre-processor ordering scheme.
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
# 893 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 893 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      DIMENSION A1(LX1,LY1,LZ1)
     $        , A2(LX1,LY1,LZ1)
     $        , A3(LX1,LY1,LZ1)
     $        , B1(LX1,LY1,LZ1)
     $        , B2(LX1,LY1,LZ1)
     $        , B3(LX1,LY1,LZ1)
C     
      CALL FACIND2 (JS1,JF1,JSKIP1,JS2,JF2,JSKIP2,IFC)
C     
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         A1(J1,J2,1)=B1(J1,J2,1)
         A2(J1,J2,1)=B2(J1,J2,1)
         A3(J1,J2,1)=B3(J1,J2,1)
  100 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine ptbgeom
C-----------------------------------------------------------------------
C     
C     Subroutine to impose perturbation to geometry before solution
C     for moving boundary problems
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
# 919 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 919 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 920 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 920 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 921 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 921 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 922 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 922 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 923 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 923 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 924 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 924 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C     
      IF (ISTEP .EQ. 0) return
      IFIELD = 0
      NEL    = NELFLD(IFIELD)
      NTOT1  = NX1*NY1*NZ1*NEL
      IMESH  = 1
      IF ( IFTMSH(IFIELD) ) IMESH = 2
C     
      CALL IBDGEOM (NEL)
      CALL ELASOLV (NEL)
      CALL UPDXYZ  (NEL)
      CALL GEOM1 (XM3,YM3,ZM3)
      CALL GEOM2
      CALL UPDMSYS (0)
      CALL VOLUME
      CALL SETINVM
      CALL LAGMASS
C     
      CALL RZERO (WX,NTOT1)
      CALL RZERO (WY,NTOT1)
      IF (NDIM.EQ.3) CALL RZERO (WZ,NTOT1)
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine ibdgeom (nel)
C     
C     Routine to evaluate mesh velocities at all moving boundaries
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
# 957 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 957 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 958 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 958 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
      
# 959 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 959 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 960 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 960 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 961 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 961 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 962 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 962 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      CHARACTER CB*1
C     
      NFACE  = 2*NDIM
      NTOT1  = NX1*NY1*NZ1*NEL
C     
      CALL RZERO (WX,NTOT1)
      CALL RZERO (WY,NTOT1)
      CALL RZERO (WZ,NTOT1)
C     
      DO 1000 ISWEEP=1,2
C     
      IFLD = 0
      DO 110 IEL=1,NEL
      DO 110 IFC=1,NFACE
         ieg = lglel(iel)
         CB  = CBC(IFC,IEL,IFLD)
      IF (CB.EQ.'M' .OR. CB.EQ.'m') THEN
         CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX1,NY1,NZ1,IFC)
         DO 140 IZ=KZ1,KZ2
         DO 140 IY=KY1,KY2
         DO 140 IX=KX1,KX2
         CALL INIGEOM (WX(IX,IY,IZ,IEL),WY(IX,IY,IZ,IEL),
     $                 WZ(IX,IY,IZ,IEL),XM1(IX,IY,IZ,IEL),
     $                 YM1(IX,IY,IZ,IEL),ZM1(IX,IY,IZ,IEL),
     $                 IFC,IEG)
  140    CONTINUE
      ENDIF
  110 CONTINUE
C     
      IF (ISWEEP.EQ.1) THEN
         CALL DSOP (WX,'MXA',NX1,NY1,NZ1)
         CALL DSOP (WY,'MXA',NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSOP (WZ,'MXA',NX1,NY1,NZ1)
      ELSE
         CALL DSOP (WX,'MNA',NX1,NY1,NZ1)
         CALL DSOP (WY,'MNA',NX1,NY1,NZ1)
         IF (NDIM.EQ.3) CALL DSOP (WZ,'MNA',NX1,NY1,NZ1)
      ENDIF
C     
 1000 CONTINUE
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine inigeom (ux,uy,uz,x,y,z,iside,iel)
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
# 1010 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1010 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 1011 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1011 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
C     
      UX  =  0.0
      UY  =  0.0
      UZ  =  0.0
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine quickmv

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
# 1021 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1021 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 1022 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1022 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ZPER" 1
      
c     Eigenvalue arrays and pointers for Global Tensor Product 
c     parameter (lelg_sm=2)
c     parameter (ltfdm2 =2)
c     parameter (lelg_sm=lelg)
c     parameter (ltfdm2=2*lx2*ly2*lz2*lelt)
c     parameter (leig2=2*lx2*lx2*(lelx*lelx+lely*lely+lelz*lelz))
c     parameter (leig =2*lx2*(lelx+lely+lelz))
      
# 10
      parameter (lfdm0 = 1-lfdm)
      parameter (lelg_sm=lfdm0+lfdm*lelg)
      parameter (ltfdm2 =lfdm0+lfdm*2*lx2*ly2*lz2*lelt)
      parameter (leig2  =lfdm0+lfdm*2*lx2*lx2
     $                             *(lelx*lelx+lely*lely+lelz*lelz))
      parameter (leig   =lfdm0+lfdm*2*lx2*(lelx+lely+lelz))
      
      
      common /peigi/ neigx ,neigy ,neigz
     $             , pvalx ,pvaly ,pvalz 
     $             , pvecx ,pvecy ,pvecz 
      integer        pvalx ,pvaly ,pvalz 
     $             , pvecx ,pvecy ,pvecz 
      
      common /eigpr/ sp(leig2),spt(leig2),eigp(leig)
     $              ,wavep(ltfdm2)
      common /eigpi/ msp(3,2),mlp(3,2)
      
      
c     Logical, array and geometry data for tensor-product box
      common /perlog/  ifycrv,ifzper,ifgfdm,ifgtp,ifemat
      logical          ifycrv,ifzper,ifgfdm,ifgtp,ifemat
      
      common /gfdmic/  nelx,nely,nelz,nelxy
     $                ,lex2pst(3),pst2lex(3)
     $                ,ngfdm_p(3),ngfdm_v(3,2)
      integer          lex2pst   ,pst2lex
      
c     Complete exchange arrays for pressure
c     common /gfdmcx/  part_in(0:lp),part_out(0:lp)
      parameter (lp_small=256)
      common /gfdmcx/  part_in(0:lp_small),part_out(0:lp_small)
     $                ,msg_id(0:lp_small,2)
     $                ,mcex
      integer          part_in,part_out
      
c     Permutation arrays for gfdm pressure solve
      common /gfdmia/  tpn1(ltfdm2),tpn2(ltfdm2)
     $                ,tpn3(ltfdm2),ind23(ltfdm2)
      integer tpn1,tpn2,tpn3
      
      parameter (lfdx =lfdm0+lfdm*lx2*lelx)
      parameter (lfdy =lfdm0+lfdm*ly2*lely)
      parameter (lfdz =lfdm0+lfdm*lz2*lelz)
      common /pergeo/  xgtp(0:lelx),ygtp(0:lely),zgtp(0:lelz)
     $                ,xmlt(lfdx),ymlt(lfdy),zmlt(lfdz)
      
      
c     Metrics for 2D x tensor-product solver
      common /permet/  rx2  (lx2,ly2,lelv)
     $               , ry2  (lx2,ly2,lelv)
     $               , sx2  (lx2,ly2,lelv)
     $               , sy2  (lx2,ly2,lelv)
     $               , w2d  (lx2,ly2,lelv)
     $               , bxyi (lx1,ly1,lelv)
      
      common /gtpcbc/  gtp_cbc(6,0:ldimt1+1)
      character*3      gtp_cbc
# 1023 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1023 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
c     
      if (if3d) then
         call quickmv3d
      else
         call quickmv2d
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine quickmv2d

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
# 1034 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1034 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 1035 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1035 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ZPER" 1
      
c     Eigenvalue arrays and pointers for Global Tensor Product 
c     parameter (lelg_sm=2)
c     parameter (ltfdm2 =2)
c     parameter (lelg_sm=lelg)
c     parameter (ltfdm2=2*lx2*ly2*lz2*lelt)
c     parameter (leig2=2*lx2*lx2*(lelx*lelx+lely*lely+lelz*lelz))
c     parameter (leig =2*lx2*(lelx+lely+lelz))
      
# 10
      parameter (lfdm0 = 1-lfdm)
      parameter (lelg_sm=lfdm0+lfdm*lelg)
      parameter (ltfdm2 =lfdm0+lfdm*2*lx2*ly2*lz2*lelt)
      parameter (leig2  =lfdm0+lfdm*2*lx2*lx2
     $                             *(lelx*lelx+lely*lely+lelz*lelz))
      parameter (leig   =lfdm0+lfdm*2*lx2*(lelx+lely+lelz))
      
      
      common /peigi/ neigx ,neigy ,neigz
     $             , pvalx ,pvaly ,pvalz 
     $             , pvecx ,pvecy ,pvecz 
      integer        pvalx ,pvaly ,pvalz 
     $             , pvecx ,pvecy ,pvecz 
      
      common /eigpr/ sp(leig2),spt(leig2),eigp(leig)
     $              ,wavep(ltfdm2)
      common /eigpi/ msp(3,2),mlp(3,2)
      
      
c     Logical, array and geometry data for tensor-product box
      common /perlog/  ifycrv,ifzper,ifgfdm,ifgtp,ifemat
      logical          ifycrv,ifzper,ifgfdm,ifgtp,ifemat
      
      common /gfdmic/  nelx,nely,nelz,nelxy
     $                ,lex2pst(3),pst2lex(3)
     $                ,ngfdm_p(3),ngfdm_v(3,2)
      integer          lex2pst   ,pst2lex
      
c     Complete exchange arrays for pressure
c     common /gfdmcx/  part_in(0:lp),part_out(0:lp)
      parameter (lp_small=256)
      common /gfdmcx/  part_in(0:lp_small),part_out(0:lp_small)
     $                ,msg_id(0:lp_small,2)
     $                ,mcex
      integer          part_in,part_out
      
c     Permutation arrays for gfdm pressure solve
      common /gfdmia/  tpn1(ltfdm2),tpn2(ltfdm2)
     $                ,tpn3(ltfdm2),ind23(ltfdm2)
      integer tpn1,tpn2,tpn3
      
      parameter (lfdx =lfdm0+lfdm*lx2*lelx)
      parameter (lfdy =lfdm0+lfdm*ly2*lely)
      parameter (lfdz =lfdm0+lfdm*lz2*lelz)
      common /pergeo/  xgtp(0:lelx),ygtp(0:lely),zgtp(0:lelz)
     $                ,xmlt(lfdx),ymlt(lfdy),zmlt(lfdz)
      
      
c     Metrics for 2D x tensor-product solver
      common /permet/  rx2  (lx2,ly2,lelv)
     $               , ry2  (lx2,ly2,lelv)
     $               , sx2  (lx2,ly2,lelv)
     $               , sy2  (lx2,ly2,lelv)
     $               , w2d  (lx2,ly2,lelv)
     $               , bxyi (lx1,ly1,lelv)
      
      common /gtpcbc/  gtp_cbc(6,0:ldimt1+1)
      character*3      gtp_cbc
# 1036 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1036 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
c     
      integer e,ex,ey,ez,eg
      common /surfa/ zsurf(lx1,lz1,lelx,lely)
     $             , wsurf(lx1,lz1,lelx,lely)
      real nxs,nys,nzs
c     
      icount = 0
      do ex=1,nelx
      do ix=1,nx1
         zsurf(ix,1,ex,1) = -1.e20
         wsurf(ix,1,ex,1) = -1.e20
         ey=nely
         eg  = ex + nelx*(ey-1) 
         mid = gllnid(eg)
         e   = gllel (eg)
         if (mid.eq.nid) then
            zsurf(ix,1,ex,1) = ym1(ix,ny1,1,e)
            vxs              = vx (ix,ny1,1,e)
            vys              = vy (ix,ny1,1,e)
            nxs              = unx(ix,1,3,e)          ! Face 3 is on top
            nys              = uny(ix,1,3,e)
            gamma_s          = (vxs*nxs + vys*nys)/(nys) 
            wsurf(ix,1,ex,1) = gamma_s           ! vertical component of
         endif
         zsurf(ix,1,ex,1) = glmax(zsurf(ix,1,ex,1),1)
         wsurf(ix,1,ex,1) = glmax(wsurf(ix,1,ex,1),1)
         icount = icount+1
      
c        write(6,6) ex,e,ix,xm1(ix,ny1,1,e),ym1(ix,ny1,1,e)
c    $   ,vxs,vys,nxs,nys,gamma_s,wsurf(ix,1,ex,1),zsurf(ix,1,ex,1)
c   6        format(3i3,1p9e12.4,' srf')
      
      enddo
      enddo
      zmin = glmin(ym1,nx1*ny1*nz1*nelv)
c     
      do ex=1,nelx
      do ix=1,nx1
         do ey=1,nely
            eg  = ex + nelx*(ey-1) 
            mid = gllnid(eg)
            e   = gllel (eg)
            if (mid.eq.nid) then
               do iy=1,ny1
                  wy (ix,iy,1,e) = wsurf(ix,1,ex,1)
     $               * (ym1(ix,iy,1,e)-zmin)/(zsurf(ix,1,ex,1)-zmin)
               enddo
            endif
         enddo
      enddo
      enddo
c     
      n = nx1*ny1*nz1*nelv
      call rzero(wx,n)
      call dsavg(wy)
      
c     call opcopy(vx,vy,vz,wx,wy,wz)
c     call outpost(vx,vy,vz,pr,t,'   ')
c     call exitt
      
      return
      end
c-----------------------------------------------------------------------
      subroutine quickmv3d

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
# 1101 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1101 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

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
# 1102 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1102 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/ZPER" 1
      
c     Eigenvalue arrays and pointers for Global Tensor Product 
c     parameter (lelg_sm=2)
c     parameter (ltfdm2 =2)
c     parameter (lelg_sm=lelg)
c     parameter (ltfdm2=2*lx2*ly2*lz2*lelt)
c     parameter (leig2=2*lx2*lx2*(lelx*lelx+lely*lely+lelz*lelz))
c     parameter (leig =2*lx2*(lelx+lely+lelz))
      
# 10
      parameter (lfdm0 = 1-lfdm)
      parameter (lelg_sm=lfdm0+lfdm*lelg)
      parameter (ltfdm2 =lfdm0+lfdm*2*lx2*ly2*lz2*lelt)
      parameter (leig2  =lfdm0+lfdm*2*lx2*lx2
     $                             *(lelx*lelx+lely*lely+lelz*lelz))
      parameter (leig   =lfdm0+lfdm*2*lx2*(lelx+lely+lelz))
      
      
      common /peigi/ neigx ,neigy ,neigz
     $             , pvalx ,pvaly ,pvalz 
     $             , pvecx ,pvecy ,pvecz 
      integer        pvalx ,pvaly ,pvalz 
     $             , pvecx ,pvecy ,pvecz 
      
      common /eigpr/ sp(leig2),spt(leig2),eigp(leig)
     $              ,wavep(ltfdm2)
      common /eigpi/ msp(3,2),mlp(3,2)
      
      
c     Logical, array and geometry data for tensor-product box
      common /perlog/  ifycrv,ifzper,ifgfdm,ifgtp,ifemat
      logical          ifycrv,ifzper,ifgfdm,ifgtp,ifemat
      
      common /gfdmic/  nelx,nely,nelz,nelxy
     $                ,lex2pst(3),pst2lex(3)
     $                ,ngfdm_p(3),ngfdm_v(3,2)
      integer          lex2pst   ,pst2lex
      
c     Complete exchange arrays for pressure
c     common /gfdmcx/  part_in(0:lp),part_out(0:lp)
      parameter (lp_small=256)
      common /gfdmcx/  part_in(0:lp_small),part_out(0:lp_small)
     $                ,msg_id(0:lp_small,2)
     $                ,mcex
      integer          part_in,part_out
      
c     Permutation arrays for gfdm pressure solve
      common /gfdmia/  tpn1(ltfdm2),tpn2(ltfdm2)
     $                ,tpn3(ltfdm2),ind23(ltfdm2)
      integer tpn1,tpn2,tpn3
      
      parameter (lfdx =lfdm0+lfdm*lx2*lelx)
      parameter (lfdy =lfdm0+lfdm*ly2*lely)
      parameter (lfdz =lfdm0+lfdm*lz2*lelz)
      common /pergeo/  xgtp(0:lelx),ygtp(0:lely),zgtp(0:lelz)
     $                ,xmlt(lfdx),ymlt(lfdy),zmlt(lfdz)
      
      
c     Metrics for 2D x tensor-product solver
      common /permet/  rx2  (lx2,ly2,lelv)
     $               , ry2  (lx2,ly2,lelv)
     $               , sx2  (lx2,ly2,lelv)
     $               , sy2  (lx2,ly2,lelv)
     $               , w2d  (lx2,ly2,lelv)
     $               , bxyi (lx1,ly1,lelv)
      
      common /gtpcbc/  gtp_cbc(6,0:ldimt1+1)
      character*3      gtp_cbc
# 1103 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f" 2
# 1103 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/mvmesh.f"
c     
      integer e,ex,ey,ez,eg
      common /surfa/ zsurf(lx1,lz1,lelx,lely)
     $             , wsurf(lx1,lz1,lelx,lely)
      real nxs,nys,nzs
c     
      icount = 0
      do ey=1,nely
      do ex=1,nelx
      do iy=1,ny1
      do ix=1,nx1
         zsurf(ix,iy,ex,ey) = -1.e20
         wsurf(ix,iy,ex,ey) = -1.e20
         ez  = nelz
         eg  = ex + nelx*(ey-1) + nelx*nely*(ez-1)
         mid = gllnid(eg)
         e   = gllel (eg)
         if (mid.eq.nid) then
            zsurf(ix,iy,ex,ey) = zm1(ix,iy,nz1,e)
            vxs                = vx (ix,iy,nz1,e)
            vys                = vy (ix,iy,nz1,e)
            vzs                = vz (ix,iy,nz1,e)
            nxs                = unx(ix,iy,6,e)     ! Face 6 is on top i
            nys                = uny(ix,iy,6,e)
            nzs                = unz(ix,iy,6,e) 
            gamma_s            = (vxs*nxs+vys*nys+vzs*nzs)/(nzs)
            wsurf(ix,iy,ex,ey) = gamma_s  ! vertical component of wsurf
         endif
         zsurf(ix,iy,ex,ey) = glmax(zsurf(ix,iy,ex,ey),1)
         wsurf(ix,iy,ex,ey) = glmax(wsurf(ix,iy,ex,ey),1)
         icount = icount+1
      enddo
      enddo
      enddo
      enddo
      
      n = nx1*ny1*nz1*nelv
      zmin = glmin(zm1,n)
c     
      do ey=1,nely
      do ex=1,nelx
      do iy=1,ny1
      do ix=1,nx1
         do ez=1,nelz
            eg  = ex + nelx*(ey-1) + nelx*nely*(ez-1)
            mid = gllnid(eg)
            e   = gllel (eg)
            if (mid.eq.nid) then
               do iz=1,nz1
                  wz (ix,iy,iz,e) = wsurf(ix,iy,ex,ey)
     $               * (zm1(ix,iy,iz,e)-zmin)/(zsurf(ix,iy,ex,ey)-zmin)
               enddo
            endif
         enddo
      enddo
      enddo
      enddo
      enddo
c     
      n = nx1*ny1*nz1*nelv
      call rzero(wx,n)
      call rzero(wy,n)
      call dsavg(wz)
c     
      return
      end
c-----------------------------------------------------------------------
