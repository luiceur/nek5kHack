c-----------------------------------------------------------------------
      subroutine nek_init(intracomm)

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'ZPER'
c
      include 'OPCTR'
      include 'CTIMER'

C     used scratch arrays
C     NOTE: no initial declaration needed. Linker will take 
c           care about the size of the CBs automatically
c
c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)
  
      real kwave2
      real*8 t0, tpp

      logical ifemati,ifsync_


      call get_session_info(intracomm)
      
C     Initalize Nek (MPI stuff, word sizes, ...)
c     call iniproc (initalized in get_session_info)

      etimes = dnekclock()
      istep  = 0
      tpp    = 0.0

      call opcount(1)

C     Initialize and set default values.
      call initdim
      call initdat
      call files

C     Read .rea +map file
      etime1 = dnekclock()
      call readat

      ifsync_ = ifsync
      ifsync = .true.

C     Initialize some variables
      call setvar  

c     Map BCs
      if (ifmoab) then
#ifdef MOAB
        call nekMOAB_bcs
#endif
      endif

c     Check for zero steps
      instep=1
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

C     Setup domain topology  
      igeom = 2
      call setup_topo

C     Compute GLL stuff (points, weights, derivate op, ...)
      call genwz

C     Initalize io unit
      call io_init

C     Set size for CVODE solver
      if(ifcvode .and. nsteps.gt.0) call cv_setsize(0,nfield)

C     USRDAT
      if(nid.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrdat' 

C     generate geometry (called after usrdat in case something changed)
      call gengeom (igeom)

      if (ifmvbd) call setup_mesh_dssum ! Set up dssum for mesh (needs geom)

C     USRDAT2
      if(nid.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrdat2' 

      call geom_reset(1)    ! recompute Jacobians, etc.
      call vrdsmsh          ! verify mesh topology

      call echopar ! echo back the parameter stack
      call setlog  ! Initalize logical flags

C     Zero out masks corresponding to Dirichlet boundary points.

      call bcmask

C     Need eigenvalues to set tolerances in presolve (SETICS)
      if (fintim.ne.0.0.or.nsteps.ne.0) call geneig (igeom)

C     Verify mesh topology
      call vrdsmsh

C     Pressure solver initialization (uses "SOLN" space as scratch)
      if (ifflow.and.(fintime.ne.0.or.nsteps.ne.0)) then
         call estrat
         if (iftran.and.solver_type.eq.'itr') then
            call set_overlap
         elseif (solver_type.eq.'fdm'.or.solver_type.eq.'pdm')then
            ifemati = .true.
            kwave2  = 0.0
            if (ifsplit) ifemati = .false.
            call gfdm_init(nx2,ny2,nz2,ifemati,kwave2)
         elseif (solver_type.eq.'25D') then
            call g25d_init
         endif
      endif

C     Initialize optional plugin
      call init_plugin

C     USRDAT3
      if(nid.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nid.eq.0) write(6,'(A,/)') ' done :: usrdat3'

C     Set initial conditions + compute field properties
      call setics
      call setprop

C     USRCHK
      if(instep.ne.0) then
        if(nid.eq.0) write(6,*) 'call userchk'
        call userchk
        if(nid.eq.0) write(6,'(A,/)') ' done :: userchk' 
      endif

C     Initialize CVODE
      if(ifcvode .and. nsteps.gt.0) call cv_init

      call comment
      call sstest (isss) 

      call dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow

      call in_situ_init()

C     Initalize timers to ZERO
      call time00
      call opcount(2)

      etims0 = dnekclock_sync()
      IF (NID.EQ.0) THEN
        WRITE (6,*) ' '
        IF (TIME.NE.0.0) WRITE (6,'(A,E14.7)') ' Initial time:',TIME
        WRITE (6,'(A,g13.5,A)') 
     &              ' Initialization successfully completed ',
     &              etims0-etimes, ' sec'
      ENDIF

      ifsync = ifsync_ ! restore initial value

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'

      real*4 papi_mflops
      integer*8 papi_flops

#ifdef _OPENACC
      include 'GEOM'
      include 'DXYZ'
      include 'IXYZ'
      include 'MASS'
      include 'WZ'
      include 'DSSUM'
      include 'SOLN'
      include 'FDMH1'
      include 'DOMAIN'

      include 'HSMG'

      include 'POPENACC'

      COMMON /SCRUZ/ TA1   (LX1,LY1,LZ1,LELT)
     $ ,             TA2   (LX1,LY1,LZ1,LELT)
     $ ,             TA3   (LX1,LY1,LZ1,LELT)
     $ ,             TA4   (LX1,LY1,LZ1,LELT)

      common /ctmp0/ work (lx2,ly2,lz2,lelt)
     $ ,             work1(lx1,ly1,lz1,lelt)
    
      common /ctmp1/ ur(lx1,ly1,lz1,lelt)
     $     ,         us(lx1,ly1,lz1,lelt)
     $     ,         ut(lx1,ly1,lz1,lelt)

      common /ctmp2/ ug(lx1*ly1*lz1*lelt)

      COMMON /CTMP3/ TMP1  (LX1,LY1,LZ1,LELT)
     $ ,             TMP2  (LX1,LY1,LZ1,LELT)
     $ ,             TMP3  (LX1,LY1,LZ1,LELT)
     $ ,             DUDR  (LX1,LY1,LZ1,LELT)
     $ ,             DUDS  (LX1,LY1,LZ1,LELT)
     $ ,             DUDT  (LX1,LY1,LZ1,LELT)

      COMMON /SCRNS/ RES1  (LX1,LY1,LZ1,LELT)
     $ ,             RES2  (LX1,LY1,LZ1,LELT)
     $ ,             RES3  (LX1,LY1,LZ1,LELT)
     $ ,             DV1   (LX1,LY1,LZ1,LELT)
     $ ,             DV2   (LX1,LY1,LZ1,LELT)
     $ ,             DV3   (LX1,LY1,LZ1,LELT)
     $ ,             RESPR (LX2,LY2,LZ2,LELT)
     $ ,             DFC   (LX1,LY1,LZ1,LELT)
     $ ,             DVC   (LX1,LY1,LZ1,LELT)

      COMMON /SCRSF/ VEXT1  (LX1*LY1*LZ1*LELT)
     $     ,         VEXT2  (LX1*LY1*LZ1*LELT)
     $     ,         VEXT3  (LX1*LY1*LZ1*LELT)

! Conflict between SCRMG and GMRES !!
      common /hsmgw2/ work2(lx1+2,ly1+2,lz1+2,lelv)
     $     ,          work3(lx1+2,ly1+2,lz1+2,lelv)

      common /SCRMG2/  r (LX1*LY1*LZ1*LELT) 
     $              , w (LX1*LY1*LZ1*LELT) 
     $              , p (LX1*LY1*LZ1*LELT) 
     $              , z (LX1*LY1*LZ1*LELT)

      common /SCRMG3/ e (2*LX1*LY1*LZ1*LELT) 

      COMMON /SCRCG/ d (LX1*LY1*LZ1*LELT)

      common /scrpre/ uc(lcr*lelt)
      common /scrpr2/ vc(lcr*lelt)


      common /scrvh/ h1    (lx1,ly1,lz1,lelt)
     $ ,             h2    (lx1,ly1,lz1,lelt)

c      parameter(lgmres=50)
      common /gmre3/ v_acc(lx2*ly2*lz2*lelv,lgmres)
     $     ,         z_acc(lx2*ly2*lz2*lelv,lgmres)

!$ACC DATA CREATE(TMP1,TMP2,TMP3,DUDR,DUDS,DUDT)
!$ACC&     CREATE(G1M1,G2M1,G3M1,G4M1,G5M1,G6M1)
!$ACC&     CREATE(DXM1,DYM1,DZM1,DXTM1,DYTM1,DZTM1)
!$ACC&     CREATE(BM1,BINVM1,JACM1,JACMI)
!$ACC&     CREATE(TA1,TA2,TA3,TA4)
!$ACC&     CREATE(UR,US,UT,WORK,WORK1)
!$ACC&     CREATE(W3M1,DXM12,DYM12,DZM12)
!$ACC&     CREATE(DXTM12,DYTM12,DZTM12)
!$ACC&     CREATE(RXM1,SXM1,TXM1,RXM2,SXM2,TXM2)
!$ACC&     CREATE(RYM1,SYM1,TYM1,RYM2,SYM2,TYM2)
!$ACC&     CREATE(RZM1,SZM1,TZM1,RZM2,SZM2,TZM2)
!$ACC&     CREATE(IDS_LGL1,IDS_LGL2,IDS_PTR,UG)
!$ACC&     CREATE(RES1,RES2,RES3,DV1,DV2,DV3,DVC)
!$ACC&     CREATE(VEXT1,VEXT2,VEXT3)
!$ACC&     CREATE(V1MASK,V2MASK,V3MASK,PMASK)
!$ACC&     CREATE(r,w,p,z,d,e)
!$ACC&     CREATE(uc,vc,ktype,fds_acc,fdst_acc)
!$ACC&     CREATE(h1_basis_acc,h1_basist_acc)
!$ACC&     CREATE(v_acc,z_acc)
!$ACC&     CREATE(ml_acc,mu_acc)

!$ACC&     CREATE(mg_lgl1,mg_ptr)
!$ACC&     CREATE(mg_imask,mg_fast_s,mg_fast_d)
!$ACC&     CREATE(mg_schwarz_wt,mg_rstr_wt)
!$ACC&     CREATE(mg_jh,mg_jht,mg_work,work2,work3)

!$ACC&     CREATE(h1,h2)
#endif 

      call nekgsync()

      if (instep.eq.0) then
        if(nid.eq.0) write(6,'(/,A,/,A,/)') 
     &     ' nsteps=0 -> skip time loop',
     &     ' running solver in post processing mode'
      else
        if(nid.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
      endif

#ifdef PAPI
      call nek_flops(papi_flops,papi_mflops)
#endif
#ifdef HPCT
      call summary_start()
      call hpm_start("nek_advance")
#endif
      isyc  = 0
      itime = 0
      if(ifsync) isyc=1
#ifndef NOTIMER
      itime = 1
#endif
      call nek_comm_settings(isyc,itime)

      call nek_comm_startstat()

#ifdef _OPENACC

!$ACC UPDATE DEVICE(G1M1,G2M1,G3M1,G4M1,G5M1,G6M1)
!$ACC UPDATE DEVICE(DXM1,DYM1,DZM1,DXTM1,DYTM1,DZTM1)
!$ACC UPDATE DEVICE(BM1,BINVM1,JACM1,JACMI,W3M1)
!$ACC UPDATE DEVICE(DXM12,DYM12,DZM12)
!$ACC UPDATE DEVICE(DXTM12,DYTM12,DZTM12)
!$ACC UPDATE DEVICE(RXM1,SXM1,TXM1,RXM2,SXM2,TXM2)
!$ACC UPDATE DEVICE(RYM1,SYM1,TYM1,RYM2,SYM2,TYM2)
!$ACC UPDATE DEVICE(RZM1,SZM1,TZM1,RZM2,SZM2,TZM2)
!$ACC UPDATE DEVICE(V1MASK,V2MASK,V3MASK,PMASK)

!$ACC UPDATE DEVICE(ids_lgl1,ids_lgl2,ids_ptr)
!$ACC UPDATE DEVICE(mg_lgl1,mg_ptr)
!$ACC UPDATE DEVICE(mg_imask,mg_fast_s,mg_fast_d)
!$ACC UPDATE DEVICE(mg_schwarz_wt,mg_rstr_wt,mg_jh,mg_jht)

#endif

      DO ISTEP=1,NSTEPS
         call nek_advance
         call userchk
         call prepost (.false.,'his')
         call in_situ_check()
         if (lastep .eq. 1) goto 1001
      ENDDO
 1001 lastep=1

      call nek_comm_settings(isyc,0)

      call comment

c     check for post-processing mode
      if (instep.eq.0) then
         nsteps=0
         istep=0
         if(nid.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nid.eq.0) write(6,*) 'done :: userchk'
         call prepost (.true.,'his')
      else
         if (nid.eq.0) write(6,'(/,A,/)') 
     $      'end of time-step loop' 
      endif

#ifdef HPCT
      call hpm_stop("nek_advance")
      call summary_stop()
#endif

#ifdef _OPENACC
!$ACC END DATA
#endif

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine nek_advance

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /cgeom/ igeom

      call nekgsync
      IF (IFTRAN) CALL SETTIME
      if (ifmhd ) call cfl_check
      CALL SETSOLV
      CALL COMMENT

      if (ifsplit) then   ! PN/PN formulation

         igeom = 1
         if (ifheat)          call heat     (igeom)
         call setprop
         call qthermal
         igeom = 1
         if (ifflow)          call fluid    (igeom)
         if (param(103).gt.0) call q_filter(param(103))
         call setup_convect (2) ! Save convective velocity _after_ filter

      else                ! PN-2/PN-2 formulation

         call setprop
         do igeom=1,ngeom

            if (igeom.gt.2) call userchk_set_xfer

            if (ifgeom) then
               call gengeom (igeom)
               call geneig  (igeom)
            endif

            if (ifmhd) then
               if (ifheat)      call heat     (igeom)
                                call induct   (igeom)
            elseif (ifpert) then
               if (ifbase.and.ifheat)  call heat          (igeom)
               if (ifbase.and.ifflow)  call fluid         (igeom)
               if (ifflow)             call fluidp        (igeom)
               if (ifheat)             call heatp         (igeom)
            else  ! std. nek case
               if (ifheat)             call heat          (igeom)
               if (ifflow)             call fluid         (igeom)
               if (ifmvbd)             call meshv         (igeom)
            endif

            if (igeom.eq.ngeom.and.param(103).gt.0) 
     $          call q_filter(param(103))

            call setup_convect (igeom) ! Save convective velocity _after_ filter

         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'OPCTR'

      if(instep.ne.0)  call runstat
      if(xxth(1).gt.0) call crs_stats(xxth(1))

   
      call in_situ_end()
      return
      end
c-----------------------------------------------------------------------
