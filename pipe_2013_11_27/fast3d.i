# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c-----------------------------------------------------------------------
      subroutine gen_fast(df,sr,ss,st,x,y,z)
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
# 7 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 7 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 8 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 8 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 9 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 9 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 10 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 10 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 11 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 11 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      parameter(lxx=lx1*lx1)
      real df(lx1*ly1*lz1,1),sr(lxx*2,1),ss(lxx*2,1),st(lxx*2,1)
c     
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt 
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
c     
      integer lbr,rbr,lbs,rbs,lbt,rbt,e
c     
      real x(nx1,ny1,nz1,nelv)
      real y(nx1,ny1,nz1,nelv)
      real z(nx1,ny1,nz1,nelv)
      real axwt(lx2)
      
      ierr = 0
      
      do e=1,nelv
c     
         if (param(44).eq.1) then
           call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,2,ierr)
         else
           call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,3,ierr)
         endif
c     
c        Set up matrices for each element.
c     
         if (param(44).eq.1) then
           call set_up_fast_1D_fem( sr(1,e),lr,nr ,lbr,rbr
     $                      ,llr(e),lmr(e),lrr(e),zgm2(1,1),nx2,e)
         else
           call set_up_fast_1D_sem( sr(1,e),lr,nr ,lbr,rbr
     $                      ,llr(e),lmr(e),lrr(e),e)
         endif
         if (ifaxis) then
            xsum = vlsum(wxm2,nx2)
            do i=1,ny2
               yavg = vlsc2(y(1,i,1,e),wxm2,nx2)/xsum
               axwt(i) = yavg
            enddo
            call set_up_fast_1D_fem_ax( ss(1,e),ls,ns ,lbs,rbs
     $                 ,lls(e),lms(e),lrs(e),zgm2(1,2),axwt,ny2,e)
         else
            if (param(44).eq.1) then
               call set_up_fast_1D_fem( ss(1,e),ls,ns ,lbs,rbs
     $                      ,lls(e),lms(e),lrs(e),zgm2(1,2),ny2,e)
            else
               call set_up_fast_1D_sem( ss(1,e),ls,ns ,lbs,rbs
     $                      ,lls(e),lms(e),lrs(e),e)
            endif
         endif
         if (if3d) then
            if (param(44).eq.1) then
               call set_up_fast_1D_fem( st(1,e),lt,nt ,lbt,rbt
     $                      ,llt(e),lmt(e),lrt(e),zgm2(1,3),nz2,e)
            else
               call set_up_fast_1D_sem( st(1,e),lt,nt ,lbt,rbt
     $                      ,llt(e),lmt(e),lrt(e),e)
            endif
         endif
c     
c        DIAGNOSTICS
c     
c        n12 = min(9,nr)
c        write(6,1) e,'1D lr',llr(e),lmr(e),lrr(e),(lr(k),k=1,n12)
c        write(6,1) e,'1D ls',lls(e),lms(e),lrs(e),(ls(k),k=1,n12)
c        if (if3d) 
c    $   write(6,1) e,'1D lt',llt(e),lmt(e),lrt(e),(lt(k),k=1,n12)
c   1    format(i6,1x,a5,1p12e12.4)
c     
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
                  df(l,e) = 1.0/diag
               else
c                 write(6,3) e,'Reset Eig in gen fast:',i,j,k,l
c    $                         ,eps,diag,lr(i),ls(j),lt(k)
c   3             format(i6,1x,a21,4i5,1p5e12.4)
                  df(l,e) = 0.0
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
                  df(l,e) = 1.0/diag
               else
c                 write(6,2) e,'Reset Eig in gen fast:',i,j,l
c    $                         ,eps,diag,lr(i),ls(j)
c   2             format(i6,1x,a21,3i5,1p4e12.4)
                  df(l,e) = 0.0
               endif
               l = l+1
            enddo
            enddo
         endif
c     
c        Next element ....
c     
      enddo
      
      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('E INVALID BC FOUND in genfast$',ierrmx)
      endif
      
      
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_fast_spacing(x,y,z)
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
# 147 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 147 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 148 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 148 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 149 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 149 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 150 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 150 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 151 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 151 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      parameter(lxx=lx1*lx1)
c     
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt 
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
c     
      integer lbr,rbr,lbs,rbs,lbt,rbt,e
c     
      real x(nx1,ny1,nz1,nelv)
      real y(nx1,ny1,nz1,nelv)
      real z(nx1,ny1,nz1,nelv)
      real axwt(lx2)
      
      ierr = 0
      
      if (param(44).eq.1) then
c                                    __ __ __
c        Now, for each element, compute lr,ls,lt between specified plane
c     
         n1 = nx2
         n2 = nx2+1
         nz0 = 1
         nzn = 1
         if (if3d) then
            nz0= 0
            nzn=n2
         endif
         eps = 1.e-7
         if (wdsize.eq.8)  eps = 1.e-14
c     
c        Find mean spacing between "left-most" planes
         call plane_space2(llr,lls,llt, 0,wxm2,x,y,z,n1,n2,nz0,nzn)
c     
c        Find mean spacing between "middle" planes
         call plane_space (lmr,lms,lmt, 1,n1,wxm2,x,y,z,n1,n2,nz0,nzn)
c     
c        Find mean spacing between "right-most" planes
         call plane_space2(lrr,lrs,lrt,n2,wxm2,x,y,z,n1,n2,nz0,nzn)
c     
      else
         call load_semhat_weighted    !   Fills the SEMHAT arrays
      endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine plane_space_std(lr,ls,lt,i1,i2,w,x,y,z,nx,nxn,nz0,nzn)
c     
c     This routine now replaced by "plane_space()"
c     
c     Here, spacing is based on arithmetic mean. 
c     New verision uses harmonic mean.  pff 2/10/07
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
# 211 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 211 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 212 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 212 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c     
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
      j2 = i2
      k2 = i2
c     
      do ie=1,nelv
c     
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
               lr2  = lr2  + ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
     $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
     $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c     
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
               ls2  = ls2  + ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
     $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
     $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c     
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
               lt2  = lt2  + ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
     $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
     $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = sqrt(lt2)
c     
         else
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
               lr2  = lr2  + ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
     $                     +   (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c     
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
               ls2  = ls2  + ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
     $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c           write(6,*) 'lrls',ie,lr(ie),ls(ie)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine plane_space(lr,ls,lt,i1,i2,w,x,y,z,nx,nxn,nz0,nzn)
c     
c     Here, spacing is based on harmonic mean.  pff 2/10/07
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
# 312 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 312 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 313 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 313 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c     
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
      j2 = i2
      k2 = i2
c     
      do ie=1,nelv
c     
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
c              lr2  = lr2  + ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
c    $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
c    $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
c    $                     *   weight
               lr2  = lr2  +   weight /
     $                       ( (x(i2,j,k,ie)-x(i1,j,k,ie))**2
     $                     +   (y(i2,j,k,ie)-y(i1,j,k,ie))**2
     $                     +   (z(i2,j,k,ie)-z(i1,j,k,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = 1./sqrt(lr2)
c     
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
c              ls2  = ls2  + ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
c    $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
c    $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
c    $                     *   weight
               ls2  = ls2  +   weight /
     $                       ( (x(i,j2,k,ie)-x(i,j1,k,ie))**2
     $                     +   (y(i,j2,k,ie)-y(i,j1,k,ie))**2
     $                     +   (z(i,j2,k,ie)-z(i,j1,k,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = 1./sqrt(ls2)
c     
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
c              lt2  = lt2  + ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
c    $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
c    $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
c    $                     *   weight
               lt2  = lt2  +   weight /
     $                       ( (x(i,j,k2,ie)-x(i,j,k1,ie))**2
     $                     +   (y(i,j,k2,ie)-y(i,j,k1,ie))**2
     $                     +   (z(i,j,k2,ie)-z(i,j,k1,ie))**2 )
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = 1./sqrt(lt2)
c     
         else              ! 2D
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
c              lr2  = lr2  + ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
c    $                     +   (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
c    $                     *   weight
               lr2  = lr2  + weight /
     $                       ( (x(i2,j,1,ie)-x(i1,j,1,ie))**2
     $                       + (y(i2,j,1,ie)-y(i1,j,1,ie))**2 )
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = 1./sqrt(lr2)
c     
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
c              ls2  = ls2  + ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
c    $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
c    $                     *   weight
               ls2  = ls2  + weight /
     $                       ( (x(i,j2,1,ie)-x(i,j1,1,ie))**2
     $                     +   (y(i,j2,1,ie)-y(i,j1,1,ie))**2 )
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = 1./sqrt(ls2)
c           write(6,*) 'lrls',ie,lr(ie),ls(ie)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine plane_space2(lr,ls,lt,i1,w,x,y,z,nx,nxn,nz0,nzn)
c     
c     Here, the local spacing is already given in the surface term.
c     This addition made to simplify the periodic bdry treatment.
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
# 432 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 432 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 433 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 433 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      real w(1),lr(1),ls(1),lt(1)
      real x(0:nxn,0:nxn,nz0:nzn,1)
      real y(0:nxn,0:nxn,nz0:nzn,1)
      real z(0:nxn,0:nxn,nz0:nzn,1)
      real lr2,ls2,lt2
c                                    __ __ __
c     Now, for each element, compute lr,ls,lt between specified planes
c     
      ny = nx
      nz = nx
      j1 = i1
      k1 = i1
c     
      do ie=1,nelv
c     
         if (if3d) then
            lr2  = 0.
            wsum = 0.
            do k=1,nz
            do j=1,ny
               weight = w(j)*w(k)
               lr2  = lr2  + ( (x(i1,j,k,ie))**2
     $                     +   (y(i1,j,k,ie))**2
     $                     +   (z(i1,j,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c     
            ls2 = 0.
            wsum = 0.
            do k=1,nz
            do i=1,nx
               weight = w(i)*w(k)
               ls2  = ls2  + ( (x(i,j1,k,ie))**2
     $                     +   (y(i,j1,k,ie))**2
     $                     +   (z(i,j1,k,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c     
            lt2 = 0.
            wsum = 0.
            do j=1,ny
            do i=1,nx
               weight = w(i)*w(j)
               lt2  = lt2  + ( (x(i,j,k1,ie))**2
     $                     +   (y(i,j,k1,ie))**2
     $                     +   (z(i,j,k1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            enddo
            lt2     = lt2/wsum
            lt(ie)  = sqrt(lt2)
c           write(6,1) 'lrlslt',ie,lr(ie),ls(ie),lt(ie)
    1       format(a6,i5,1p3e12.4)
c     
         else
            lr2 = 0.
            wsum = 0.
            do j=1,ny
               weight = w(j)
               lr2  = lr2  + ( (x(i1,j,1,ie))**2
     $                     +   (y(i1,j,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            lr2     = lr2/wsum
            lr(ie)  = sqrt(lr2)
c     
            ls2 = 0.
            wsum = 0.
            do i=1,nx
               weight = w(i)
               ls2  = ls2  + ( (x(i,j1,1,ie))**2
     $                     +   (y(i,j1,1,ie))**2 )
     $                     *   weight
               wsum = wsum + weight
            enddo
            ls2     = ls2/wsum
            ls(ie)  = sqrt(ls2)
c           write(6,*) 'lrls',ie,lr(ie),ls(ie),lt(ie)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_fem(s,lam,n,lbc,rbc,ll,lm,lr,z,nz,e)
      real s(1),lam(1),ll,lm,lr,z(1)
      integer lbc,rbc,e
c     
      parameter (m=100)
      real dx(0:m)
      integer icalld
      save    icalld
      data    icalld/0/
c     
      icalld=icalld+1
c     
      if (nz.gt.m-3) then
         write(6,*) 'ABORT. Error in set_up_fast_1D_fem. Increase m to'
     $             , nz
         call exitt
      endif
c     
c     In the present scheme, each element is viewed as a d-fold
c     tensor of (1+nz+1) arrays, even if funky bc's are applied 
c     on either end of the 1D array.
c     
      n = nz+2
c     
c     Compute spacing, dx()
c     
      call set_up_1D_geom(dx,lbc,rbc,ll,lm,lr,z,nz)
c     
      nn1 = n*n + 1
      call gen_eigs_A_fem(s,s(nn1),lam,n,dx,lbc,rbc)
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_1D_geom(dx,lbc,rbc,ll,lm,lr,z,nz)
c     
      real dx(0:1),ll,lm,lr,z(2)
      integer lbc,rbc
c     
c     
c     Set up the 1D geometry for the tensor-product based overlapping Sc
c     
c     Upon return: 
c     
c       dx() contains the spacing required to set up the stiffness matri
c     
c     
c     Input:
c     
c       lbc (rbc) is 0 if the left (right) BC is Dirichlet, 1 if Neumann
c     
c       ll is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the LEFT element
c     
c       lm is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the MIDDLE element
c     
c       lr is the space between the right-most Gauss point of the middle
c          element and the left-most Gauss point of the RIGHT element
c     
c       --- if ll (lr) is very small (0), it indicates that there is no
c           left (right) spacing, and that the left (right) BC is Neuman
c     
c     
c       z() is the array of nz Gauss points on the interval ]-1,1[.
c     
c     Boundary conditions:
c     
c     bc = 0  -- std. Dirichlet bc applied 2 points away from interior
c     bc = 1  -- Dirichlet bc applied 1 point away from interior (outflo
c     bc = 2  -- Neumann bc applied on interior point (W,v,V,SYM,...)
c     
c     
c     
c     Geometry:
c     
c     
c        dx0       dx1   dx2    dx3   dx5    dx5        dx6
c     
c    bl        |<--ll-->|<------lm------>|<---lr--->|           br
c     0--------x-----|--x---x--------x---x--|-------x-----------0
c     
c       left elem.         middle elem.         right elem.
c                   -1                     +1
c     
c     
c    "bl" = (extrapolated) location of Gauss point NX2-1 in left elem.
c     
c    "br" = (extrapolated) location of Gauss point 2 in right elem.
c     
c    Overlapping Schwarz applied with homogeneous Dirichlet boundary
c    conditions at "bl" and "br", and with a single d.o.f. extending
c    in to each adjacent domain.
c     
      eps = 1.e-5
      call rone(dx,nz+3)
c     
c     Middle
      scale = lm/(z(nz)-z(1))
      do i=1,nz-1
         dx(i+1) = (z(i+1)-z(i))*scale
      enddo
      
c     Left end
      if (lbc.eq.0) then
         dzm0   = z(1) + 1.
         dxm0   = scale*dzm0
         dxln   = ll - dxm0
         scalel = dxln/dzm0
         dx(0)  = scalel*(z(2)-z(1))
         dx(1)  = ll
      elseif (lbc.eq.1) then
         dx(1)  = ll
      endif
c     
c     Right end
      if (rbc.eq.0) then
         dzm0      = z(1) + 1.
         dxm0      = scale*dzm0
         dxr0      = lr - dxm0
         scaler    = dxr0/dzm0
         dx(nz+1)  = lr
         dx(nz+2)  = scaler*(z(2)-z(1))
      elseif (rbc.eq.1) then
         dx(nz+1)  = lr
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_eigs_A_fem(sf,sft,atd,n,l,lbc,rbc)
c     
c     Set up Fast diagonalization solver for FEM on mesh 2
c     
      real sf(n,n),sft(n,n),atd(1),l(0:1)
      integer lbc,rbc
c     
      parameter (m=100)
      real atu(m),ad(m),au(m),c(m),bh(m),li(0:m)
c     
      if (n.gt.m) then
         write(6,*) 'ABORT. Error in gen_eigs_A_fem. Increase m to',n
         call exitt
      endif
c     
c     Get delta x's
c     
      do i=0,n
         li(i) = 1.0/l(i)
      enddo
c                          ^   ^
c     Fill initial arrays, A & B:
c     
      call rzero(ad,n)
      call rzero(au,n-1)
      call rzero(bh,n)
c     
      ie1 = lbc
      ien = n-rbc
      do ie=ie1,ien
c     
c        il,ir are the left and right endpts of element ie.
         il = ie
         ir = ie+1
c     
         if (ie.gt.0) ad(il) = ad(il) +       li(ie)
         if (ie.lt.n) ad(ir) = ad(ir) +       li(ie)
         if (ie.gt.0) au(il) =        -       li(ie)
         if (ie.gt.0) bh(il) = bh(il) + 0.5 * l(ie)
         if (ie.lt.n) bh(ir) = bh(ir) + 0.5 * l(ie)
      enddo
c     
c     Take care of bc's (using blasting)
      bhm = vlmax(bh(2),n-2)/(n-2)
      ahm = vlmax(ad(2),n-2)/(n-2)
c     
      if (lbc.gt.0) then
         au(1) = 0.
         ad(1) = ahm
         bh(1) = bhm
      endif
c     
      if (rbc.gt.0) then
         au(n-1) = 0.
         ad(n  ) = ahm
         bh(n  ) = bhm
      endif
c     
c     
      do i=1,n
         c(i) = sqrt(1.0/bh(i))
      enddo
c                                        ~
c     Scale rows and columns of A by C:  A = CAC
c     
      do i=1,n
         atd(i) = c(i)*ad(i)*c(i)
      enddo
c     
c     Scale upper diagonal
c     
      atu(1) = 0.
      do i=1,n-1
         atu(i) = c(i)*au(i)*c(i+1)
      enddo
c                                             ~
c     Compute eigenvalues and eigenvectors of A
c     
      call calcz(atd,atu,n,dmax,dmin,sf,ierr)
      if (ierr.eq.1) then
         nid = mynode()
         write(6,6) nid,' czfail:',(l(k),k=0,n)
    6    format(i5,a8,1p16e10.2)
         call exitt
      endif
c     
c     Sort eigenvalues and vectors
c     
      call sort(atd,atu,n)
      call transpose(sft,n,sf,n)
      do j=1,n
         call swap(sft(1,j),atu,n,au)
      enddo
      call transpose(sf,n,sft,n)
c     
c     Make "like" eigenvectors of same sign (for ease of diagnostics)
c     
      do j=1,n
         avg = vlsum(sf(1,j),n)
         if (avg.lt.0) call chsign(sf(1,j),n)
      enddo
c     
c     Clean up zero eigenvalue
c     
      eps = 1.0e-6*dmax
      do i=1,n
c        if (atd(i).lt.eps) 
c    $      write(6,*) 'Reset Eig in gen_Afem:',i,n,atd(i)
         if (atd(i).lt.eps) atd(i) = 0.0
      enddo
c     
c     scale eigenvectors by C:
c     
      do i=1,n
         do j=1,n
            sf(i,j) = sf(i,j)*c(i)
         enddo
      enddo
c                                        ^
c     Orthonormalize eigenvectors w.r.t. B inner-product
c     
      do j=1,n
         alpha = vlsc3(bh,sf(1,j),sf(1,j),n)
         alpha = 1.0/sqrt(alpha)
         call cmult(sf(1,j),alpha,n)
      enddo
c     
c     Diagnostics
c     
c     do j=1,n
c        do i=1,n
c           sft(i,j) = vlsc3(bh,sf(1,i),sf(1,j),n)
c        enddo
c     enddo
c     
c     n8 = min(n,8)
c     do i=1,n
c        write(6,2) (sft(i,j),j=1,n8)
c     enddo
c   2 format('Id:',1p8e12.4)
c     
      call transpose(sft,n,sf,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,e,bsym,ierr)
      integer                lbr,rbr,lbs,rbs,lbt,rbt,e,bsym
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
# 806 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 806 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 807 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 807 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 808 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 808 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 809 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 809 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 810 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 810 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      integer fbc(6)
c     
c     ibc = 0  <==>  Dirichlet
c     ibc = 1  <==>  Dirichlet, outflow (no extension)
c     ibc = 2  <==>  Neumann,   
      
      
      do iface=1,2*ndim
         ied = eface(iface)
         ibc = -1
      
         if (ifmhd) call mhd_bc_dn(ibc,iface,e) ! can be overwritten by 
      
         if (cbc(ied,e,ifield).eq.'   ') ibc = 0
         if (cbc(ied,e,ifield).eq.'E  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'msi') ibc = 0
         if (cbc(ied,e,ifield).eq.'MSI') ibc = 0
         if (cbc(ied,e,ifield).eq.'P  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'p  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'O  ') ibc = 1
         if (cbc(ied,e,ifield).eq.'ON ') ibc = 1
         if (cbc(ied,e,ifield).eq.'o  ') ibc = 1
         if (cbc(ied,e,ifield).eq.'on ') ibc = 1
         if (cbc(ied,e,ifield).eq.'MS ') ibc = 1
         if (cbc(ied,e,ifield).eq.'ms ') ibc = 1
         if (cbc(ied,e,ifield).eq.'MM ') ibc = 1
         if (cbc(ied,e,ifield).eq.'mm ') ibc = 1
         if (cbc(ied,e,ifield).eq.'mv ') ibc = 2
         if (cbc(ied,e,ifield).eq.'mvn') ibc = 2
         if (cbc(ied,e,ifield).eq.'v  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'V  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'W  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'SYM') ibc = bsym
         if (cbc(ied,e,ifield).eq.'SL ') ibc = 2
         if (cbc(ied,e,ifield).eq.'sl ') ibc = 2
         if (cbc(ied,e,ifield).eq.'SHL') ibc = 2
         if (cbc(ied,e,ifield).eq.'shl') ibc = 2
         if (cbc(ied,e,ifield).eq.'A  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'S  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'s  ') ibc = 2
         if (cbc(ied,e,ifield).eq.'J  ') ibc = 0
         if (cbc(ied,e,ifield).eq.'SP ') ibc = 0
      
         fbc(iface) = ibc
      
         if (ierr.eq.-1) write(6,1) ibc,ied,e,ifield,cbc(ied,e,ifield)
  1      format(2i3,i8,i3,2x,a3,'  get_fast_bc_error')
      
      enddo
      
      if (ierr.eq.-1) call exitti('Error A get_fast_bc$',e)
      
      lbr = fbc(1)
      rbr = fbc(2)
      lbs = fbc(3)
      rbs = fbc(4)
      lbt = fbc(5)
      rbt = fbc(6)
      
      ierr = 0 
      if (ibc.lt.0) ierr = lglel(e)
      
c     write(6,6) e,lbr,rbr,lbs,rbs,(cbc(k,e,ifield),k=1,4)
c   6 format(i5,2x,4i3,3x,4(1x,a3),'  get_fast_bc')
      
      return
      end
c-----------------------------------------------------------------------
      subroutine outv(x,n,name3)
      character*3 name3
      real x(1)
c     
      nn = min (n,10)
      write(6,6) name3,(x(i),i=1,nn)
    6 format(a3,10f12.6)
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat(a,m,n,name6,ie)
      real a(m,n)
      character*6 name6
c     
      write(6,*) 
      write(6,*) ie,' matrix: ',name6,m,n
      n12 = min(n,12)
      do i=1,m
         write(6,6) ie,name6,(a(i,j),j=1,n12)
      enddo
    6 format(i3,1x,a6,12f9.5)
      write(6,*) 
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_fem_ax
     $   (s,lam,n,lbc,rbc,ll,lm,lr,z,y,nz,ie)
      real s(1),lam(1),ll,lm,lr,z(1),y(1)
      integer lbc,rbc
c     
      parameter (m=100)
      real dx(0:m)
      integer icalld
      save    icalld
      data    icalld/0/
c     
      icalld=icalld+1
c     
      if (nz.gt.m-3) then
         write(6,*) 'ABORT. Error in set_up_fast_1D_fem. Increase m to'
     $             , nz
         call exitt
      endif
c     
c     In the present scheme, each element is viewed as a d-fold
c     tensor of (1+nz+1) arrays, even if funky bc's are applied 
c     on either end of the 1D array.
c     
      n = nz+2
c     
c     Compute spacing, dx()
c     
      call set_up_1D_geom_ax(dx,lbc,rbc,ll,lm,lr,z,y,nz)
c     
      nn1 = n*n + 1
      call gen_eigs_A_fem_ax(s,s(nn1),lam,n,dx,y,lbc,rbc)
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_1D_geom_ax(dx,lbc,rbc,ll,lm,lr,z,y,nz)
c     
      real dx(0:1),ll,lm,lr,z(2),y(1)
      integer lbc,rbc
c     
c     
c     Set up the 1D geometry for the tensor-product based overlapping Sc
c     
c     Upon return: 
c     
c       dx() contains the spacing required to set up the stiffness matri
c     
c     
c     Input:
c     
c       lbc (rbc) is 0 if the left (right) BC is Dirichlet, 1 if Neumann
c     
c       ll is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the LEFT element
c     
c       lm is the space between the left-most Gauss point of the middle
c          element and the right-most Gauss point of the MIDDLE element
c     
c       lr is the space between the right-most Gauss point of the middle
c          element and the left-most Gauss point of the RIGHT element
c     
c       --- if ll (lr) is very small (0), it indicates that there is no
c           left (right) spacing, and that the left (right) BC is Neuman
c     
c     
c       z() is the array of nz Gauss points on the interval ]-1,1[.
c     
c     Boundary conditions:
c     
c     bc = 0  -- std. Dirichlet bc applied 2 points away from interior
c     bc = 1  -- Dirichlet bc applied 1 point away from interior (outflo
c     bc = 2  -- Neumann bc applied on interior point (W,v,V,SYM,...)
c     
c     
c     
c     Geometry:
c     
c     
c        dx0       dx1   dx2    dx3   dx5    dx5        dx6
c     
c    bl        |<--ll-->|<------lm------>|<---lr--->|           br
c     0--------x-----|--x---x--------x---x--|-------x-----------0
c     
c       left elem.         middle elem.         right elem.
c                   -1                     +1
c     
c     
c    "bl" = (extrapolated) location of Gauss point NX2-1 in left elem.
c     
c    "br" = (extrapolated) location of Gauss point 2 in right elem.
c     
c    Overlapping Schwarz applied with homogeneous Dirichlet boundary
c    conditions at "bl" and "br", and with a single d.o.f. extending
c    in to each adjacent domain.
c     
      eps = 1.e-5
      call rone(dx,nz+3)
c     
c     Middle
      scale = lm/(z(nz)-z(1))
      do i=1,nz-1
         dx(i+1) = (z(i+1)-z(i))*scale
      enddo
      
c     Left end
      if (lbc.eq.0) then
         dzm0   = z(1) + 1.
         dxm0   = scale*dzm0
         dxln   = ll - dxm0
         scalel = dxln/dzm0
         dx(0)  = scalel*(z(2)-z(1))
         dx(1)  = ll
      elseif (lbc.eq.1) then
         dx(1)  = ll
      endif
c     
c     Right end
      if (rbc.eq.0) then
         dzm0      = z(1) + 1.
         dxm0      = scale*dzm0
         dxr0      = lr - dxm0
         scaler    = dxr0/dzm0
         dx(nz+1)  = lr
         dx(nz+2)  = scaler*(z(2)-z(1))
      elseif (rbc.eq.1) then
         dx(nz+1)  = lr
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_eigs_A_fem_ax(sf,sft,atd,n,l,y,lbc,rbc)
c     
c     Set up Fast diagonalization solver for FEM on mesh 2
c     
      real sf(n,n),sft(n,n),atd(1),l(0:1),y(1)
      integer lbc,rbc
c     
      parameter (m=100)
      real atu(m),ad(m),au(m),c(m),bh(m),li(0:m)
c     
      if (n.gt.m) then
         write(6,*) 'ABORT. Error in gen_eigs_A_fem. Increase m to',n
         call exitt
      endif
c     
c     Get delta x's
c     
      do i=0,n
         li(i) = 1.0/l(i)
      enddo
c                          ^   ^
c     Fill initial arrays, A & B:
c     
      call rzero(ad,n)
      call rzero(au,n-1)
      call rzero(bh,n)
c     
      ie1 = lbc
      ien = n-rbc
      do ie=ie1,ien
c     
c        il,ir are the left and right endpts of element ie.
         il = ie
         ir = ie+1
c     
         if (ie.gt.0) ad(il) = ad(il) +       li(ie)
         if (ie.lt.n) ad(ir) = ad(ir) +       li(ie)
         if (ie.gt.0) au(il) =        -       li(ie)
         if (ie.gt.0) bh(il) = bh(il) + 0.5 * l(ie)
         if (ie.lt.n) bh(ir) = bh(ir) + 0.5 * l(ie)
      enddo
c     
c     Take care of bc's (using blasting)
      bhm = vlmax(bh(2),n-2)/(n-2)
      ahm = vlmax(ad(2),n-2)/(n-2)
c     
      if (lbc.gt.0) then
         au(1) = 0.
         ad(1) = ahm
         bh(1) = bhm
      endif
c     
      if (rbc.gt.0) then
         au(n-1) = 0.
         ad(n  ) = ahm
         bh(n  ) = bhm
      endif
c     
c     
      do i=1,n
         c(i) = sqrt(1.0/bh(i))
      enddo
c                                        ~
c     Scale rows and columns of A by C:  A = CAC
c     
      do i=1,n
         atd(i) = c(i)*ad(i)*c(i)
      enddo
c     
c     Scale upper diagonal
c     
      atu(1) = 0.
      do i=1,n-1
         atu(i) = c(i)*au(i)*c(i+1)
      enddo
c                                             ~
c     Compute eigenvalues and eigenvectors of A
c     
      call calcz(atd,atu,n,dmax,dmin,sf,ierr)
      if (ierr.eq.1) then
         nid = mynode()
         write(6,6) nid,' czfail2:',(l(k),k=0,n)
    6    format(i5,a8,1p16e10.2)
         call exitt
      endif
c     
c     Sort eigenvalues and vectors
c     
      call sort(atd,atu,n)
      call transpose(sft,n,sf,n)
      do j=1,n
         call swap(sft(1,j),atu,n,au)
      enddo
      call transpose(sf,n,sft,n)
c     
c     Make "like" eigenvectors of same sign (for ease of diagnostics)
c     
      do j=1,n
         avg = vlsum(sf(1,j),n)
         if (avg.lt.0) call chsign(sf(1,j),n)
      enddo
c     
c     Clean up zero eigenvalue
c     
      eps = 1.0e-6*dmax
      do i=1,n
c        if (atd(i).lt.eps) 
c    $      write(6,*) 'Reset Eig in gen_Afm_ax:',i,n,atd(i)
         if (atd(i).lt.eps) atd(i) = 0.0
      enddo
c     
c     scale eigenvectors by C:
c     
      do i=1,n
         do j=1,n
            sf(i,j) = sf(i,j)*c(i)
         enddo
      enddo
c                                        ^
c     Orthonormalize eigenvectors w.r.t. B inner-product
c     
      do j=1,n
         alpha = vlsc3(bh,sf(1,j),sf(1,j),n)
         alpha = 1.0/sqrt(alpha)
         call cmult(sf(1,j),alpha,n)
      enddo
c     
c     Diagnostics
c     
c     do j=1,n
c        do i=1,n
c           sft(i,j) = vlsc3(bh,sf(1,i),sf(1,j),n)
c        enddo
c     enddo
c     
c     n8 = min(n,8)
c     do i=1,n
c        write(6,2) (sft(i,j),j=1,n8)
c     enddo
c   2 format('Id:',1p8e12.4)
c     
      call transpose(sft,n,sf,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine load_semhat_weighted    !   Fills the SEMHAT arrays
c     
c     Note that this routine performs the following matrix multiplies
c     after getting the matrices back from semhat:
c     
c     dgl = bgl dgl
c     jgl = bgl jgl
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
# 1190 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1190 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 1191 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1191 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      nr = nx1-1
      call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,nr,wh)
      call do_semhat_weight(jgl,dgl,bgl,nr)
c     
      return
      end
c----------------------------------------------------------------------
      subroutine do_semhat_weight(jgl,dgl,bgl,n)
      real bgl(1:n-1),jgl(1:n-1,0:n),dgl(1:n-1,0:n)
      
      do j=0,n
         do i=1,n-1
            jgl(i,j)=bgl(i)*jgl(i,j)
         enddo
      enddo
      do j=0,n
         do i=1,n-1
            dgl(i,j)=bgl(i)*dgl(i,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine semhat(a,b,c,d,z,dgll,jgll,bgl,zgl,dgl,jgl,n,w)
c     
c     Generate matrices for single element, 1D operators:
c     
c        a    = Laplacian
c        b    = diagonal mass matrix
c        c    = convection operator b*d
c        d    = derivative matrix
c        dgll = derivative matrix,    mapping from pressure nodes to vel
c        jgll = interpolation matrix, mapping from pressure nodes to vel
c        z    = GLL points
c     
c        zgl  = GL points
c        bgl  = diagonal mass matrix on GL
c        dgl  = derivative matrix,    mapping from velocity nodes to pre
c        jgl  = interpolation matrix, mapping from velocity nodes to pre
c     
c        n    = polynomial degree (velocity space)
c        w    = work array of size 2*n+2
c     
c     Currently, this is set up for pressure nodes on the interior GLL p
c     
c     
      real a(0:n,0:n),b(0:n),c(0:n,0:n),d(0:n,0:n),z(0:n)
      real dgll(0:n,1:n-1),jgll(0:n,1:n-1)
c     
      real bgl(1:n-1),zgl(1:n-1)
      real dgl(1:n-1,0:n),jgl(1:n-1,0:n)
c     
c      real w(0:1) 
      real w(0:2*n+1)
      
      np = n+1
      nm = n-1
      n2 = n-2
c     
      call zwgll (z,b,np)
c     
      do i=0,n
         call fd_weights_full(z(i),z,n,1,w)
         do j=0,n
            d(i,j) = w(j+np)                   !  Derivative matrix
         enddo
      enddo
      
      if (n.eq.1) return                       !  No interpolation for n
      
      do i=0,n
         call fd_weights_full(z(i),z(1),n2,1,w(1))
         do j=1,nm
            jgll(i,j) = w(j   )                  !  Interpolation matrix
            dgll(i,j) = w(j+nm)                  !  Derivative    matrix
         enddo
      enddo
c     
      call rzero(a,np*np)
      do j=0,n
      do i=0,n
         do k=0,n
            a(i,j) = a(i,j) + d(k,i)*b(k)*d(k,j)
         enddo
         c(i,j) = b(i)*d(i,j)
      enddo
      enddo
c     
      call zwgl (zgl,bgl,nm)
c     
      do i=1,n-1
         call fd_weights_full(zgl(i),z,n,1,w)
         do j=0,n
            jgl(i,j) = w(j   )                  !  Interpolation matrix
            dgl(i,j) = w(j+np)                  !  Derivative    matrix
         enddo
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine fd_weights_full(xx,x,n,m,c)
c     
c     This routine evaluates the derivative based on all points
c     in the stencils.  It is more memory efficient than "fd_weights"
c     
c     This set of routines comes from the appendix of 
c     A Practical Guide to Pseudospectral Methods, B. Fornberg
c     Cambridge Univ. Press, 1996.   (pff)
c     
c     Input parameters:
c       xx -- point at wich the approximations are to be accurate
c       x  -- array of x-ordinates:   x(0:n)
c       n  -- polynomial degree of interpolant (# of points := n+1)
c       m  -- highest order of derivative to be approxxmated at xi
c     
c     Output:
c       c  -- set of coefficients c(0:n,0:m).
c             c(j,k) is to be applied at x(j) when
c             the kth derivative is approxxmated by a 
c             stencil extending over x(0),x(1),...x(n).
c     
c     
      real x(0:n),c(0:n,0:m)
c     
      c1       = 1.
      c4       = x(0) - xx
c     
      do k=0,m
      do j=0,n
         c(j,k) = 0.
      enddo
      enddo
      c(0,0) = 1.
c     
      do i=1,n
         mn = min(i,m)
         c2 = 1.
         c5 = c4
         c4 = x(i)-xx
         do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
            do k=mn,1,-1
               c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            enddo
            c(i,0) = -c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo
c     call outmat(c,n+1,m+1,'fdw',n)
      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_fast_1D_sem(s,lam,n,lbc,rbc,ll,lm,lr,ie)

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
# 1352 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1352 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 1353 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1353 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
c     
      common /fast1dsem/ g(lr2),w(lr2)
c     
      real g,w
      real s(1),lam(1),ll,lm,lr
      integer lbc,rbc
      
      integer bb0,bb1,eb0,eb1,n,n1
      logical l,r
      
      n=nx1-1
      !bcs on E are from normal vel component
      if(lbc.eq.2 .or. lbc.eq.3) then !wall,sym - dirichlet velocity
         eb0=1
      else !outflow,element - neumann velocity
         eb0=0
      endif
      if(rbc.eq.2 .or. rbc.eq.3) then !wall,sym - dirichlet velocity
         eb1=n-1
      else !outflow,element - neumann velocity
         eb1=n
      endif
      !bcs on B are from tangent vel component
      if(lbc.eq.2) then !wall - dirichlet velocity
         bb0=1
      else !outflow,element,sym - neumann velocity
         bb0=0
      endif
      if(rbc.eq.2) then !wall - dirichlet velocity
         bb1=n-1
      else !outflow,element,sym - neumann velocity
         bb1=n
      endif
c     
      l = (lbc.eq.0)
      r = (rbc.eq.0)
c     
c     calculate E tilde operator
      call set_up_fast_1D_sem_op(s,eb0,eb1,l,r,ll,lm,lr,bh,dgl,0)
c     call outmat(s,n+1,n+1,'  Et  ',ie)
c     calculate B tilde operator
      call set_up_fast_1D_sem_op(g,bb0,bb1,l,r,ll,lm,lr,bh,jgl,1)
c     call outmat(g,n+1,n+1,'  Bt  ',ie)
      
      n=n+1
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
      subroutine set_up_fast_1D_sem_op(g,b0,b1,l,r,ll,lm,lr,bh,jgl,jscl)
c            -1 T
c     G = J B  J
c     
c     gives the inexact restriction of this matrix to
c     an element plus one node on either side
c     
c     g - the output matrix
c     b0, b1 - the range for Bhat indices for the element
c              (enforces boundary conditions)
c     l, r - whether there is a left or right neighbor
c     ll,lm,lr - lengths of left, middle, and right elements
c     bh - hat matrix for B
c     jgl - hat matrix for J (should map vel to pressure)
c     jscl - how J scales
c            0: J = Jh
c            1: J = (L/2) Jh
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
# 1434 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1434 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
      real g(0:lx1-1,0:lx1-1)
      real bh(0:lx1-1),jgl(1:lx2,0:lx1-1)
      real ll,lm,lr
      integer b0,b1
      logical l,r
      integer jscl
c     
      real bl(0:lx1-1),bm(0:lx1-1),br(0:lx1-1)
      real gl,gm,gr,gll,glm,gmm,gmr,grr
      real fac
      integer n
      n=nx1-1
c     
c     
c     compute the scale factors for J      
      if (jscl.eq.0) then
         gl=1.
         gm=1.
         gr=1.
      elseif (jscl.eq.1) then
         gl=0.5*ll
         gm=0.5*lm
         gr=0.5*lr
      endif
      gll = gl*gl
      glm = gl*gm
      gmm = gm*gm
      gmr = gm*gr
      grr = gr*gr
c     
c     compute the summed inverse mass matrices for
c     the middle, left, and right elements
      do i=1,n-1
         bm(i)=2. /(lm*bh(i))
      enddo
      if (b0.eq.0) then
         bm(0)=0.5*lm*bh(0)
         if(l) bm(0)=bm(0)+0.5*ll*bh(n)
         bm(0)=1. /bm(0)
      endif
      if (b1.eq.n) then
         bm(n)=0.5*lm*bh(n)
         if(r) bm(n)=bm(n)+0.5*lr*bh(0)
         bm(n)=1. /bm(n)
      endif
c     note that in computing bl for the left element,
c     bl(0) is missing the contribution from its left neighbor
      if (l) then
         do i=0,n-1
            bl(i)=2. /(ll*bh(i))
         enddo
         bl(n)=bm(0)
      endif
c     note that in computing br for the right element,
c     br(n) is missing the contribution from its right neighbor
      if (r) then
         do i=1,n
            br(i)=2. /(lr*bh(i))
         enddo
         br(0)=bm(n)
      endif
c      
      call rzero(g,(n+1)*(n+1))
      do j=1,n-1
         do i=1,n-1
            do k=b0,b1
               g(i,j) = g(i,j) + gmm*jgl(i,k)*bm(k)*jgl(j,k)
            enddo
         enddo
      enddo
c      
      if (l) then
         do i=1,n-1
            g(i,0) = glm*jgl(i,0)*bm(0)*jgl(n-1,n)
            g(0,i) = g(i,0)
         enddo
c        the following is inexact
c        the neighbors bc's are ignored, and the contribution
c        from the neighbor's neighbor is left out
c        that is, bl(0) could be off as noted above
c        or maybe i should go from 1 to n
         do i=0,n
            g(0,0) = g(0,0) + gll*jgl(n-1,i)*bl(i)*jgl(n-1,i)
         enddo
      else
         g(0,0)=1.
      endif
c      
      if (r) then
         do i=1,n-1
            g(i,n) = gmr*jgl(i,n)*bm(n)*jgl(1,0)
            g(n,i) = g(i,n)
         enddo
c        the following is inexact
c        the neighbors bc's are ignored, and the contribution
c        from the neighbor's neighbor is left out
c        that is, br(n) could be off as noted above
c        or maybe i should go from 0 to n-1
         do i=0,n
            g(n,n) = g(n,n) + grr*jgl(1,i)*br(i)*jgl(1,i)
         enddo
      else
         g(n,n)=1.
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine swap_lengths
      

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
# 1544 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1544 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 1545 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1545 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 1546 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1546 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 1547 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1547 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
      common /swaplengths/ l(lx1,ly1,lz1,lelv)
      common /ctmpf/  lr(2*lx1+4),ls(2*lx1+4),lt(2*lx1+4)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
      
      real l,l2d
      integer e
      
      n2 = nx1-1
      nz0 = 1
      nzn = 1
      nx  = nx1-2
      if (if3d) then
         nz0 = 0
         nzn = n2
      endif
      call plane_space(lmr,lms,lmt,0,n2,wxm1,xm1,ym1,zm1,nx,n2,nz0,nzn)
      
      n=n2+1
      if (if3d) then
         do e=1,nelv
         do j=2,n2
         do k=2,n2
            l(1,k,j,e) = lmr(e)
            l(n,k,j,e) = lmr(e)
            l(k,1,j,e) = lms(e)
            l(k,n,j,e) = lms(e)
            l(k,j,1,e) = lmt(e)
            l(k,j,n,e) = lmt(e)
         enddo
         enddo
         enddo
         call dssum(l,n,n,n)
         do e=1,nelv
            llr(e) = l(1,2,2,e)-lmr(e)
            lrr(e) = l(n,2,2,e)-lmr(e)
            lls(e) = l(2,1,2,e)-lms(e)
            lrs(e) = l(2,n,2,e)-lms(e)
            llt(e) = l(2,2,1,e)-lmt(e)
            lrt(e) = l(2,2,n,e)-lmt(e)
         enddo
      else
         do e=1,nelv
         do j=2,n2
            l(1,j,1,e) = lmr(e)
            l(n,j,1,e) = lmr(e)
            l(j,1,1,e) = lms(e)
            l(j,n,1,e) = lms(e)
c           call outmat(l(1,1,1,e),n,n,' L    ',e)
         enddo
         enddo
c        call outmat(l(1,1,1,25),n,n,' L    ',25)
         call dssum(l,n,n,1)
c        call outmat(l(1,1,1,25),n,n,' L    ',25)
         do e=1,nelv
c           call outmat(l(1,1,1,e),n,n,' L    ',e)
            llr(e) = l(1,2,1,e)-lmr(e)
            lrr(e) = l(n,2,1,e)-lmr(e)
            lls(e) = l(2,1,1,e)-lms(e)
            lrs(e) = l(2,n,1,e)-lms(e)
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine row_zero(a,m,n,e)
      integer m,n,e
      real a(m,n)
      do j=1,n
         a(e,j)=0.
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine mhd_bc_dn(ibc,face,e)
      integer                  face,e
c     
c     sets Neumann BC on pressure (ibc=2) for face and e(lement) except
c     when ifield normal component has (homogeneous) Neumann
c     boundary condition setting ibc=1 (i.e. Direchlet BC on pressure)
c     
c     Note: 'SYM' on a plane with r,s,t-normal is 'dnn','ndn','nnd'? bsy
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
# 1636 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1636 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 1637 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1637 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
      
# 1638 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1638 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"

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
# 1639 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f" 2
# 1639 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/fast3d.f"
      
      ied = eface(face)	! symmetric -> preprocessor notation
      nfc = face+1
      nfc = nfc/2	! = 1,2,3 for face 1 & 2,3 & 4,5 & 6
      
      if (indx1(cbc(ied,e,ifield),'d',1).gt.0)   ibc=2
      
      if (indx1(cbc(ied,e,ifield),'n',1).gt.nfc) ibc=1 ! 'n' for V_n
      
      return
      end
c-----------------------------------------------------------------------
