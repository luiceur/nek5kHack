# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
c-----------------------------------------------------------------------
      subroutine genbox
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
# 5 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 5 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
      
# 6 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 6 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 7 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 7 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
c     
      character*132 string
      character*1  string1(132)
      equivalence (string,string1)
C     
c     
      call blank (string,132)
c     call gets  (string,132,iend,9) ! read name of "box" file from .rea
c     read (9,132) string       ! *** MESH?
  132 format(a132)
c     
      if (nid.eq.0) then
c     
c        Read in name of previously generated NEKTON data set.
c        (Must have same dimension and number of fields as current run)
c     
         write(6,*) 'opening file:'
c        write(6,*) string
         lfname = ltrunc(reafle,132) - 4
         call blank (string,132)
         call chcopy(string,reafle,lfname)
         call chcopy(string1(lfname+1),'.box',4)
         open (unit=7,file=string,status='old',err=999)
         call gets(string,132,iend,7)  ! This is a dummy read, for compa
c     
c        here is where the 2d/3d determination is made....
c     
         call geti1(ndim,iend,7)
         ndim = iabs(ndim)
         call getrv(rfld,1,iend,7)
         nfld = int(rfld)                  ! Determine number of fields
         if(rfld.ne.nfld) nfld = nfld +1
         call getbox(xgtp,ygtp,zgtp,nfld)  ! Read in the .box file
         close (unit=7)
c     
      endif
c     
      call bcpbox(nfld)
c     
      call mapbox(nelx,nely,nelz)  ! Map elements to processors
      call makebox(xgtp,ygtp,zgtp,nfld) ! Fill up local elements on _thi
c     
c     call outbox_mesh
c     
      ifgfdm = .true.
      param(116) = nelx
      param(117) = nely
      param(118) = nelz
c     
      if3d = .false.
      if (ndim.eq.3) if3d=.true.
c     
      return
      
  999 continue
      if (nid.eq.0) write(6,*) 'ABORT: Could not find box file ',string
      call exitt
      
      end
c-----------------------------------------------------------------------
      subroutine gen_gtp_vertex (vertex,ncrnr)
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
# 70 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 70 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
      
# 71 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 71 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 72 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 72 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 73 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 73 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
c     
c     integer vertex((2**ldim)*lelt)    ! local -- 2D for now !?  long?
      integer vertex(2**ldim,1)         ! local -- 2D for now !?  long?
      integer e,eg,ex,ey,ez
c     
      common /ctmp1/ indx(lelx+1)
     $             , indy(lely+1)
     $             , indz(lelz+1)
c     
      nptsx = nelx + 1
      nptsy = nely + 1
      nptsz = nelz + 1
c     
      call jjnt(indx,nptsx)
      call jjnt(indy,nptsy)
      call jjnt(indz,nptsz)
c     
      ifld = 2                                          ! ifield.gt.2?
      if (ifflow) ifld = 1
c     
      nbc = 2*ndim
      write(6,6) nid,ifld,(gtp_cbc(k,ifld),k=1,nbc)
   6  format(2i4,'  GTP BC:',6(2x,a3))
c     
      if (gtp_cbc(1,ifld).eq.'P  ') indx(nptsx) = 1  ! Set last point to
      if (gtp_cbc(3,ifld).eq.'P  ') indy(nptsy) = 1
      if (gtp_cbc(5,ifld).eq.'P  ') indz(nptsz) = 1
c     
      if (gtp_cbc(1,ifld).eq.'P  ') nptsx = nelx
      if (gtp_cbc(3,ifld).eq.'P  ') nptsy = nely
      if (gtp_cbc(5,ifld).eq.'P  ') nptsz = nelz
c     
c     
      do e=1,nelv                                       ! nelt/nel?
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)
         if (if3d) then
            l = 0
            do k=0,1
            do j=0,1
            do i=0,1
               l = l+1
               i1 = i+ex
               j1 = j+ey
               k1 = k+ez
               vertex(l,e) =  indx(i1)
     $                     + (indy(j1)-1)*nptsx
     $                     + (indz(k1)-1)*nptsx*nptsy
            enddo
            enddo
            enddo
         else
            l = 0
            do j=0,1
            do i=0,1
               l = l+1
               i1 = i+ex
               j1 = j+ey
               vertex(l,e) =  indx(i1)
     $                     + (indy(j1)-1)*nptsx
            enddo
            enddo
         endif
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine getbox(x,y,z,nfld)
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
# 144 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 144 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
      
# 145 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 145 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 146 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 146 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
c     
      real x(0:1),y(0:1),z(0:1)
c     
      logical ifevenx,ifeveny,ifevenz,if_multi_seg
      character*1 boxcirc
c     
      integer nelxyz(3)
c     
c     
      ifevenx = .false.
      ifeveny = .false.
      ifevenz = .false.
c     
      call gets(boxcirc,1,iend,7)
      if (iend.eq.1) goto 99
      if (boxcirc.eq.'c'.or.boxcirc.eq.'C') 
     $                    call getr2(xc,yc,iend,7)
      if_multi_seg = .false.
      if (boxcirc.eq.'m'.or.boxcirc.eq.'M') 
     $      if_multi_seg = .true.
c     
      if (if_multi_seg) then
         maxx = 1
         call get_multi_seg
     $        (nelxyz,x(0),y(0),z(0),maxx,if3d)
         nelx = nelxyz(1)
         nely = nelxyz(2)
         nelz = 1
         if (if3d) nelz = nelxyz(3)
         ifld0 = 2
         if (ifflow) ifld0 = 1
         ifld1 = ifld0-1 + nfld
         nbc = 2*ndim
         do ifld=ifld0,ifld1
            call getcv(gtp_cbc(1,ifld),3,6,iend,7)
            write(6,*) 'CBC1',(gtp_cbc(k,ifld),k=1,nbc),ifld
         enddo
      elseif (if3d) then
         call geti3(nelx,nely,nelz,iend,7)
c     
         if (nelx.lt.0) ifevenx = .true.
         if (nely.lt.0) ifeveny = .true.
         if (nelz.lt.0) ifevenz = .true.
         nelx = abs(nelx)
         nely = abs(nely)
         nelz = abs(nelz)
c     
         if (nelx.gt.lelx.or.nely.gt.lely.or.nelz.gt.lelz) then
            write(6,*) 'Abort, increase lelx,..lelz and recompile',
     $                    nelx,nely,nelz,lelx,lely,lelz
            call exitt
         endif
c     
         nery = nely
         if (boxcirc.eq.'C') nery = 1
         if (iend.eq.1) goto 99
         ninbox = nelx*nely*nelz
         write(6,6) ninbox,nelx,nely,nelz,nfld
c     
         if (ifevenx) then
            call getr3(x0,x1,ratio,iend,7)
            if (ratio.le.0) ratio=1.
            dx = (x1-x0)/nelx
            x(0) = 0.
            do k=1,nelx
               x(k) = x(k-1) + dx
               dx        = dx*ratio
            enddo
            scale = (x1-x0)/x(nelx)
            do k=0,nelx
               x(k) = x0 + scale*x(k)
            enddo
         else
            call getrv(x(0),nelx+1,iend,7)
         endif
c     
         if (ifeveny) then
            call getr3(y0,y1,ratio,iend,7)
            if (ratio.le.0) ratio=1.
            dy = (y1-y0)/nely
            y(0) = 0.
            do k=1,nely
               y(k) = y(k-1) + dy
               dy        = dy*ratio
            enddo
            scale = (y1-y0)/y(nely)
            do k=0,nely
               y(k) = y0 + scale*y(k)
            enddo
         else
            call getrv(y(0),nely+1,iend,7)
         endif
c     
         if (ifevenz) then
            call getr3(z0,z1,ratio,iend,7)
            if (ratio.le.0) ratio=1.
            dz = (z1-z0)/nelz
            z(0) = 0.
            do k=1,nelz
               z(k) = z(k-1) + dz
               dz        = dz*ratio
            enddo
            scale = (z1-z0)/z(nelz)
            do k=0,nelz
               z(k) = z0 + scale*z(k)
            enddo
         else
            call getrv(z(0),nelz+1,iend,7)
         endif
c     
         ifld0 = 2
         if (ifflow) ifld0 = 1
         ifld1 = ifld0-1 + nfld
         nbc   = 2*ndim
         write(6,*) 'this is ifld01:',ifld0,ifld1
         do ifld=ifld0,ifld1
            call getcv(gtp_cbc(1,ifld),3,6,iend,7)
            write(6,*) 'CBC2',(gtp_cbc(k,ifld),k=1,nbc),ifld
         enddo
      else
         nelz = 1
         call geti2(nelx,nely,iend,7)
c     
         if (nelx.lt.0) ifevenx = .true.
         if (nely.lt.0) ifeveny = .true.
         nelx = abs(nelx)
         nely = abs(nely)
         if (nelx.gt.lelx.or.nely.gt.lely) then
            write(6,*) 'Abort, increase lelx,lely and recompile',
     $                    nelx,nely,lelx,lely
            call exitt
         endif
c     
         nery = nely
         if (boxcirc.eq.'C') nery = 1
         if (iend.eq.1) goto 99
         ninbox = nelx*nely*nelz
         write(6,6) ninbox,nelx,nely,nery,nfld
c     
         if (ifevenx) then
            call getr3(x0,x1,ratio,iend,7)
            if (ratio.le.0) ratio=1.
            dx = (x1-x0)/nelx
            x(0) = 0.
            do k=1,nelx
               x(k) = x(k-1) + dx
               dx        = dx*ratio
            enddo
            scale = (x1-x0)/x(nelx)
            do k=0,nelx
               x(k) = x0 + scale*x(k)
            enddo
         else
            call getrv(x(0),nelx+1,iend,7)
         endif
c     
         if (ifeveny) then
            call getr3(y0,y1,ratio,iend,7)
            if (ratio.le.0) ratio=1.
            dy = (y1-y0)/nely
            y(0) = 0.
            do k=1,nely
               y(k) = y(k-1) + dy
               dy        = dy*ratio
            enddo
            scale = (y1-y0)/y(nely)
            do k=0,nely
               y(k) = y0 + scale*y(k)
            enddo
         else
            call getrv(y(0),nely+1,iend,7)
         endif
c     
         ifld0 = 2
         if (ifflow) ifld0 = 1
         ifld1 = ifld0-1 + nfld
         nbc   = 2*ndim
         do ifld=ifld0,ifld1
            call getcv(gtp_cbc(1,ifld),3,4,iend,7)
            write(6,*) 'CBC3',(gtp_cbc(k,ifld),k=1,nbc),ifld
         enddo
      endif
      nel = nelx*nely*nelz
    6 format('Reading',i12,' =',3i6,' elements for box',i3,'.')
   99 continue
c     
      if (nid.eq.0) then
      write(6,*) 'Beginning construction of new mesh.'
     $           ,nel,nelx,nely,nelz
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine mapbox(melx,mely,melz)
c     
c     Given melx  x  mely  x  melz elements, distribute to np processors
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
# 345 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 345 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 346 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 346 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 347 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 347 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
c     
      nelx = melx
      nely = mely
      nelz = melz
c     
      nelgv = nelx*nely*nelz
      nelgt = nelx*nely*nelz
c     
      if (nelgt.gt.lelg.or.nelgv.gt.lelg) then
         write(6,*) nid,' NELGT too large:',nelgt,nelgv,lelg
         call exitt
      endif
c     
      call mapelpr
c     
      if (nelt.gt.lelt.or.nelv.gt.lelv) then
         write(6,*) nid,' NEL too large:',nelt,nelv,lelv,lelt
         call exitt
      endif
c     
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine makebox(x,y,z,nfld)
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
# 374 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 374 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
      
# 375 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 375 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 376 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 376 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 377 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 377 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
c     
      real x(0:1),y(0:1),z(0:1)
      integer e,eg,ex,ey,ez
c     
      character*3 cbc1,   cbc2,   cbc3,   cbc4,   cbc5,   cbc6
      real        rbc1(5),rbc2(5),rbc3(5),rbc4(5),rbc5(5),rbc6(5)
c     
      integer eface(6)
      save    eface
      data    eface  / 4,2,1,3,5,6 /
c     
c     Set up TEMPORARY value for NFIELD - NFLDT
c     
      nfldt = 1
      if (ifheat) nfldt=2+npscal
      if (ifmhd ) nfldt=2+npscal+1
      nbcs      = nfldt
      ibcs      = 2
      if (ifflow) ibcs = 1
c     
      nbc = 5*6*lelt*(ldimt1+1)
      call rzero(bc,nbc)
c     
      eg = 0
      do ez=1,nelz
      do ey=1,nely
      do ex=1,nelx
         eg   = eg+1
         jnid = gllnid(eg)
         if (jnid.eq.nid) then
            e = gllel(eg)
c     
            xc(1,e) = x(ex-1)
            xc(2,e) = x(ex  )
            xc(3,e) = x(ex  )
            xc(4,e) = x(ex-1)
c     
            yc(1,e) = y(ey-1)
            yc(2,e) = y(ey-1)
            yc(3,e) = y(ey  )
            yc(4,e) = y(ey  )
c     
            if (if3d) then
c     
               xc(5,e) = x(ex-1)
               xc(6,e) = x(ex  )
               xc(7,e) = x(ex  )
               xc(8,e) = x(ex-1)
c     
               yc(5,e) = y(ey-1)
               yc(6,e) = y(ey-1)
               yc(7,e) = y(ey  )
               yc(8,e) = y(ey  )
c     
               zc(1,e) = z(ez-1)
               zc(2,e) = z(ez-1)
               zc(3,e) = z(ez-1)
               zc(4,e) = z(ez-1)
c     
               zc(5,e) = z(ez  )
               zc(6,e) = z(ez  )
               zc(7,e) = z(ez  )
               zc(8,e) = z(ez  )
c     
            endif
c     
            do ifld=ibcs,nbcs
c     
               ilftz = mod1(ez-1+nelz,nelz)
               irgtz = mod1(ez+1+nelz,nelz)
               ilfty = mod1(ey-1+nely,nely)
               irgty = mod1(ey+1+nely,nely)
               ilftx = mod1(ex-1+nelx,nelx)
               irgtx = mod1(ex+1+nelx,nelx)
c     
               rbc1(1) = nely*nelx*(ez-1)+nelx*(ey-1) + ilftx
               rbc1(2) = eface(2)
               rbc2(1) = nely*nelx*(ez-1)+nelx*(ey-1) + irgtx
               rbc2(2) = eface(1)
c     
               rbc3(1) = nely*nelx*(ez-1)+nelx*(ilfty-1) + ex
               rbc3(2) = eface(4)
               rbc4(1) = nely*nelx*(ez-1)+nelx*(irgty-1) + ex
               rbc4(2) = eface(3)
c     
               rbc5(1) = nely*nelx*(ilftz-1)+nelx*(ey-1) + ex
               rbc5(2) = eface(6)
               rbc6(1) = nely*nelx*(irgtz-1)+nelx*(ey-1) + ex
               rbc6(2) = eface(5)
c     
               if (ex.eq.1) then
                  cbc1=gtp_cbc(1,ifld)
                  if (cbc1.ne.'P  ') call rzero(rbc1,5)
               else
                  cbc1='E  '
               endif
c     
               if (ex.eq.nelx) then
                  cbc2=gtp_cbc(2,ifld)
                  if (cbc2.ne.'P  ') call rzero(rbc2,5)
               else
                  cbc2='E  '
               endif
c     
               if (ey.eq.1) then
                  cbc3=gtp_cbc(3,ifld)
                  if (cbc3.ne.'P  ') call rzero(rbc3,5)
               else
                  cbc3='E  '
               endif
c     
               if (ey.eq.nely) then
                  cbc4=gtp_cbc(4,ifld)
                  if (cbc4.ne.'P  ') call rzero(rbc4,5)
               else
                  cbc4='E  '
               endif
c     
               if (ez.eq.1) then
                  cbc5=gtp_cbc(5,ifld)
                  if (cbc5.ne.'P  ') call rzero(rbc5,5)
               else
                  cbc5='E  '
               endif
c     
               if (ez.eq.nelz) then
                  cbc6=gtp_cbc(6,ifld)
                  if (cbc6.ne.'P  ') call rzero(rbc6,5)
               else
                  cbc6='E  '
               endif
c     
               cbc(1,e,ifld) = cbc3 
               cbc(2,e,ifld) = cbc2
               cbc(3,e,ifld) = cbc4
               cbc(4,e,ifld) = cbc1
               cbc(5,e,ifld) = cbc5
               cbc(6,e,ifld) = cbc6
c     
               call copy(bc(1,1,e,ifld),rbc3,5)
               call copy(bc(1,2,e,ifld),rbc2,5)
               call copy(bc(1,3,e,ifld),rbc4,5)
               call copy(bc(1,4,e,ifld),rbc1,5)
               call copy(bc(1,5,e,ifld),rbc5,5)
               call copy(bc(1,6,e,ifld),rbc6,5)
c     
            enddo
         endif
c     
      enddo
      enddo
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine geti1(i,iend,io)
c     
c     Get an integer from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) i
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine geti2(i1,i2,iend,io)
c     
c     Get an integer from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) i1,i2
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine geti3(i1,i2,i3,iend,io)
c     
c     Get an integer from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) i1,i2,i3
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine getr2(r1,r2,iend,io)
c     
c     Get two reals from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) r1,r2
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine getr3(r1,r2,r3,iend,io)
c     
c     Get two reals from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) r1,r2,r3
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine getrv(r,n,iend,io)
      real r(n)
      parameter (big=1.e29)
c     
c     Get reals from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
c     
c     serious hack for f90
      call cfill(r,big,n)
c     
      open(unit=99,file='box.tmp')
      read(99,*,end=1) (r(k),k=1,n)
      close(unit=99)
      return
    1 continue
      close(unit=99)
      do k=1,n
         if (r(k).eq.big) goto 15
      enddo
   15 continue
      write(6,*) 'this is k:',k,n
      read( 7,*,end=2) (r(j),j=k,n)
      return
    2 continue
      iend=1
      return
      end
c-----------------------------------------------------------------------
      subroutine getiv(r,n,iend,io)
      integer r(n)
      parameter (ibig=99999999)
c     
c     Get integers from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
c     
c     serious hack for f90
      call ifill(r,ibig,n)
c     
      open(unit=99,file='box.tmp')
      read(99,*,end=1) (r(k),k=1,n)
      close(unit=99)
      return
    1 continue
      close(unit=99)
      do k=1,n
         if (r(k).eq.ibig) goto 15
      enddo
   15 continue
      write(6,*) 'this is k:',k,n
      read( 7,*,end=2) (r(j),j=k,n)
      return
    2 continue
      iend=1
      return
      end
c-----------------------------------------------------------------------
      subroutine getcv0(c,m,n,iend,io)
      character*1 c(m,n)
      character*1 adum(132)
c     
c     Get character strings, with no seperator, from first uncommented l
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,1,end=2) (adum(k),k=1,132)
    1 format(132a1)
    2 continue
c     
      i = 0
      do l=1,n
         do k=1,m
            i = i+1
            c(k,l) = adum(i)
         enddo
      enddo
      do j=1,n
         write(6,3) m,n,(c(i,j),i=1,m)
    3    format(2i4,'getcv0:',132a1)
      enddo
c     
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine getcv(c,m,n,iend,io)
      character*1 c(m,n)
      character*1 adum(132)
c     
c     Get character strings, with single seperator, from first uncomment
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,1,end=2) (adum(k),k=1,132)
    1 format(132a1)
    2 continue
c     
      i = 0
      do l=1,n
         do k=1,m
            i = i+1
            c(k,l) = adum(i)
         enddo
c     
c        bump pointer for "," in input string
         i = i+1
      enddo
      do j=1,n
         write(6,3) m,n,(c(i,j),i=1,m)
    3    format(2i4,'getcv:',132a1)
      enddo
c     
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine gets(c,n,iend,io)
      character*1 c(n)
c     
c     Get character string from first uncommented line
c     
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,1) (c(k),k=1,n)
    1 format(132a1)
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine get_multi_seg(nelxyz,x,y,z,m,if3d)
c     
      integer nelxyz(3)
      real x(0:m),y(0:m),z(0:m)
      logical if3d
c     
      integer nels (1000)
      real    gains(1000)
      real    xs   (0:1000)
c     
c     This routine allows the user to input a complex sequence of
c     segments for each of the x, y and z directions, so that
c     genbox will generate the desired tensor product array.
c     
c     Input format for each coordinat is
c     
c     nsegs
c     nel_1    nel_2 ...  nel_nsegs
c     gain_1  gain_2 ... gain_nsegs
c     x_0, x_1, x_2, ...    x_nsegs
c     
c     The output will be:
c     
c     x(0)...x(nel_1)...x_(nel_2)...x(nel_segs),
c     
c     where the segment between x(e-1) and x(e) is filled with
c     a distribution determined by gain_e.   Uniform spacing
c     corresponds to gain_e = 1, otherwise, a geometric sequence
c     is generated, with dx_i+1 = dx_i * gain_e
c     
      ndim = 2
      if (if3d) ndim=3
c     
      do id=1,ndim
c     
         call geti1(nseg,iend,7)
         call getiv(nels ,nseg  ,iend,7)
         call getrv(xs   ,nseg+1,iend,7)
         call getrv(gains,nseg  ,iend,7)
         write(6,*) 'NSEG: ',id,nseg,(nels(k),k=1,nseg)
c     
         nelxyz(id) = ilsum(nels,nseg)
         if (nelxyz(id).gt.m) then
            write(6,*) 'NEL too large:',nelxyz(id),m,id
            write(6,*) 'Abort in get_multi_seg.'
            call exitt
         endif
c     
         k = 0
         do is=1,nseg
            nel = nels(is)
            x0  = xs(is-1)
            x1  = xs(is)
            if (id.eq.1) call geometric_x(x(k),nel,x0,x1,gains(is))
            if (id.eq.2) call geometric_x(y(k),nel,x0,x1,gains(is))
            if (id.eq.3) call geometric_x(z(k),nel,x0,x1,gains(is))
            k = k + nel
         enddo
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine geometric_x(x,n,x0,x1,gain)
c     
      real x(0:n)
c     
      x(0) = 0.
      dx = (x1-x0)/n
      do i=1,n
         x(i) = x(i-1) + dx
         dx   = dx*gain
      enddo
c     
      scale = (x1-x0)/(x(n)-x(0))
      do i=0,n
         x(i) = x0 + scale*x(i)
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine get_xyz_distribution (x,nelx)
      real x(0:2)
      if (nelx.ge.0) return ! else, uniform or geometric spacing
c     
      x0 = x(0)
      x1 = x(1)
      r  = x(2)
      write(6,*) 'x0:',x0,x1,r,nelx
c     
      nelx = abs(nelx)
      if (r.eq.1.0) then
         dx = (x1-x0)/nelx
         do i=1,nelx
            x(i) = x0 + i*dx
         enddo
         x(nelx) = x1
      else
         dx = (x1-x0)/nelx
         x(0) = 0.
         do i=1,nelx
            x(i) = x(i-1) + dx
            dx = r*dx
         enddo
         scale = (x1-x0)/x(nelx)
         call cmult(x,scale,nelx+1)
         call cadd (x,x0   ,nelx+1)
         x(nelx) = x1
      endif
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine scannocom(iend,infile)
c     
c     scan through infile until "no comment" is found
c     
      character*132 string
      character*1  string1(132)
      equivalence (string1,string)
c     
      character*1 comment
      save        comment
      data        comment /'#'/
c     
      iend = 0
      do line=1,100000
         call blank(string,132)
         read (infile ,80,end=100,err=100) string
c     
         write(6,80) string
         if   (indx1(string,comment,1).ne.1) then
              open(unit=99,file='box.tmp')
              len = ltrunc(string,132)
              write(99,81) (string1(k),k=1,len)
              close(unit=99)
              return
         else
              i1 = indx1(string,comment,1)
c             write(6,*) i1,' i1', comment
         endif
c     
      enddo
   80 format(a132)
   81 format(132a1)
c     
  100 continue
      iend = 1
      return
      end
c-----------------------------------------------------------------------
      function ilsum(x,n)
      integer x(1)
      ilsum = 0
      do i=1,n
         ilsum = ilsum + x(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine jjnt(x,n)
      integer x(1)
      do i=1,n
         x(i) = i
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine bcpbox(nfld)
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
# 900 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 900 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
      
# 901 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 901 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 902 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 902 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 903 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 903 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
c     
      call bcast(ndim,isize)
      call bcast(nfld,isize)
      call bcast(nelx,isize)
      call bcast(nely,isize)
      call bcast(nelz,isize)
c     
      lnx = wdsize*(nelx+1)
      lny = wdsize*(nely+1)
      lnz = wdsize*(nelz+1)
      call bcast(xgtp,lnx)
      call bcast(ygtp,lny)
      call bcast(zgtp,lnz)
c     
      ifld0 = 2
      if (ifflow) ifld0 = 1
      ifld1 = ifld0-1 + nfld
      do ifld=ifld0,ifld1
         call bcast(gtp_cbc(1,ifld),18)
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine outbox_mesh

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
# 929 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 929 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
      
# 930 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 930 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 931 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 931 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"

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
# 932 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f" 2
# 932 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/genbox.f"
      integer e
c     
      if (np.gt.1) return
      if (ndim.eq.2) return
c     
      do e=1,nelv
         write(88,1) e
         write(88,2) (xc(k,e),k=1,4)
         write(88,2) (yc(k,e),k=1,4)
         write(88,2) (zc(k,e),k=1,4)
         write(88,2) (xc(k,e),k=5,8)
         write(88,2) (yc(k,e),k=5,8)
         write(88,2) (zc(k,e),k=5,8)
      enddo
    1 format('            ELEMENT',i12,' [    1a]    GROUP     0')
    2 format(1p4e15.7)
c     
      return
      end
c-----------------------------------------------------------------------
