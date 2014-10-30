# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"
      subroutine byte_sync_mpi(mpi_fh)
      




      
# 8
      return
      end
C-----------------------------------------------------------------------
      subroutine byte_open_mpi(fname,mpi_fh,ierr)
      

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
# 14 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 14 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/RESTART" 1
      
c     parameter (lelr=max(lelt,lelg/16)) ! THIS IS THE MEMORY conservati
# 3
      parameter (lelr=lelg)              ! THIS IS THE MEMORY INTENSIVE 
      
      common /crst_i/ max_rst            ! for full restart
      
      common /cmfi_i/ nxr,nyr,nzr,nelr,nelgr,istpr,ifiler,nfiler
     $              , nxo,nyo,nzo,nrg
     $              , wdsizr,wdsizo
     $              , nfileo,nproc_o,nfldr
     $              , er(lelr),nelB,nelBr
      integer wdsizr,wdsizo,er
      
      parameter(iHeaderSize=132)
      
      common /cmfi_r/ timer
      
      common /cmfi_c/ ihdr,rdcode,mfi_fname
      character*3  ihdr
      character*10 rdcode
      character*80 mfi_fname
      
      character*1  rdcode1(10)
      equivalence (rdcode,rdcode1)
      
      common /cmfi_l/ 
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps (ldimt1),ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr(ldimt1),ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      logical
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps         ,ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr        ,ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      
      common /cmfi_p/ fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00 
      integer         fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00
      
      integer          nekcomm_io,ifh_mbyte
      common /i4mpiio/ nekcomm_io,ifh_mbyte
# 15 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 15 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"
      
















# 32
      write(6,*) 'byte_open_mpi: No MPI-IO support!'
      ierr=1
      return

# 36
      ierr=0
      return
      end
C-----------------------------------------------------------------------
      subroutine byte_read_mpi(buf,icount,iorank,mpi_fh,ierr)
      

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
# 43 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 43 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/RESTART" 1
      
c     parameter (lelr=max(lelt,lelg/16)) ! THIS IS THE MEMORY conservati
# 3
      parameter (lelr=lelg)              ! THIS IS THE MEMORY INTENSIVE 
      
      common /crst_i/ max_rst            ! for full restart
      
      common /cmfi_i/ nxr,nyr,nzr,nelr,nelgr,istpr,ifiler,nfiler
     $              , nxo,nyo,nzo,nrg
     $              , wdsizr,wdsizo
     $              , nfileo,nproc_o,nfldr
     $              , er(lelr),nelB,nelBr
      integer wdsizr,wdsizo,er
      
      parameter(iHeaderSize=132)
      
      common /cmfi_r/ timer
      
      common /cmfi_c/ ihdr,rdcode,mfi_fname
      character*3  ihdr
      character*10 rdcode
      character*80 mfi_fname
      
      character*1  rdcode1(10)
      equivalence (rdcode,rdcode1)
      
      common /cmfi_l/ 
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps (ldimt1),ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr(ldimt1),ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      logical
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps         ,ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr        ,ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      
      common /cmfi_p/ fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00 
      integer         fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00
      
      integer          nekcomm_io,ifh_mbyte
      common /i4mpiio/ nekcomm_io,ifh_mbyte
# 44 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 44 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"
      






















# 67
      write(6,*) 'byte_read_mpi: No MPI-IO support!'
      ierr=1
      return

      
# 72
      ierr=0
      
      return
      end
C-----------------------------------------------------------------------
      subroutine byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)
      

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
# 80 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 80 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/RESTART" 1
      
c     parameter (lelr=max(lelt,lelg/16)) ! THIS IS THE MEMORY conservati
# 3
      parameter (lelr=lelg)              ! THIS IS THE MEMORY INTENSIVE 
      
      common /crst_i/ max_rst            ! for full restart
      
      common /cmfi_i/ nxr,nyr,nzr,nelr,nelgr,istpr,ifiler,nfiler
     $              , nxo,nyo,nzo,nrg
     $              , wdsizr,wdsizo
     $              , nfileo,nproc_o,nfldr
     $              , er(lelr),nelB,nelBr
      integer wdsizr,wdsizo,er
      
      parameter(iHeaderSize=132)
      
      common /cmfi_r/ timer
      
      common /cmfi_c/ ihdr,rdcode,mfi_fname
      character*3  ihdr
      character*10 rdcode
      character*80 mfi_fname
      
      character*1  rdcode1(10)
      equivalence (rdcode,rdcode1)
      
      common /cmfi_l/ 
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps (ldimt1),ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr(ldimt1),ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      logical
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps         ,ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr        ,ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      
      common /cmfi_p/ fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00 
      integer         fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00
      
      integer          nekcomm_io,ifh_mbyte
      common /i4mpiio/ nekcomm_io,ifh_mbyte
# 81 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 81 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"
      






















# 104
      write(6,*) 'byte_write_mpi: No MPI-IO support!'
      ierr=1
      return

# 108
      ierr=0
      return
      end
C-----------------------------------------------------------------------
      subroutine byte_close_mpi(mpi_fh,ierr)
      

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
# 115 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 115 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/RESTART" 1
      
c     parameter (lelr=max(lelt,lelg/16)) ! THIS IS THE MEMORY conservati
# 3
      parameter (lelr=lelg)              ! THIS IS THE MEMORY INTENSIVE 
      
      common /crst_i/ max_rst            ! for full restart
      
      common /cmfi_i/ nxr,nyr,nzr,nelr,nelgr,istpr,ifiler,nfiler
     $              , nxo,nyo,nzo,nrg
     $              , wdsizr,wdsizo
     $              , nfileo,nproc_o,nfldr
     $              , er(lelr),nelB,nelBr
      integer wdsizr,wdsizo,er
      
      parameter(iHeaderSize=132)
      
      common /cmfi_r/ timer
      
      common /cmfi_c/ ihdr,rdcode,mfi_fname
      character*3  ihdr
      character*10 rdcode
      character*80 mfi_fname
      
      character*1  rdcode1(10)
      equivalence (rdcode,rdcode1)
      
      common /cmfi_l/ 
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps (ldimt1),ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr(ldimt1),ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      logical
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps         ,ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr        ,ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      
      common /cmfi_p/ fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00 
      integer         fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00
      
      integer          nekcomm_io,ifh_mbyte
      common /i4mpiio/ nekcomm_io,ifh_mbyte
# 116 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 116 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"
      










# 127
      if(nid.eq.0) write(6,*) 'byte_close_mpi: No MPI-IO support!'
      ierr=1
      return

      
# 132
      return
      end
C-----------------------------------------------------------------------
      subroutine byte_set_view(ioff_in,mpi_fh)
      

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
# 138 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 138 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/RESTART" 1
      
c     parameter (lelr=max(lelt,lelg/16)) ! THIS IS THE MEMORY conservati
# 3
      parameter (lelr=lelg)              ! THIS IS THE MEMORY INTENSIVE 
      
      common /crst_i/ max_rst            ! for full restart
      
      common /cmfi_i/ nxr,nyr,nzr,nelr,nelgr,istpr,ifiler,nfiler
     $              , nxo,nyo,nzo,nrg
     $              , wdsizr,wdsizo
     $              , nfileo,nproc_o,nfldr
     $              , er(lelr),nelB,nelBr
      integer wdsizr,wdsizo,er
      
      parameter(iHeaderSize=132)
      
      common /cmfi_r/ timer
      
      common /cmfi_c/ ihdr,rdcode,mfi_fname
      character*3  ihdr
      character*10 rdcode
      character*80 mfi_fname
      
      character*1  rdcode1(10)
      equivalence (rdcode,rdcode1)
      
      common /cmfi_l/ 
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps (ldimt1),ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr(ldimt1),ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      logical
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps         ,ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr        ,ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      
      common /cmfi_p/ fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00 
      integer         fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00
      
      integer          nekcomm_io,ifh_mbyte
      common /i4mpiio/ nekcomm_io,ifh_mbyte
# 139 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 139 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"
      


















      
# 159
      return
      end
C-----------------------------------------------------------------------
      subroutine nek_comm_io(nn)
      

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
# 165 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 165 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/RESTART" 1
      
c     parameter (lelr=max(lelt,lelg/16)) ! THIS IS THE MEMORY conservati
# 3
      parameter (lelr=lelg)              ! THIS IS THE MEMORY INTENSIVE 
      
      common /crst_i/ max_rst            ! for full restart
      
      common /cmfi_i/ nxr,nyr,nzr,nelr,nelgr,istpr,ifiler,nfiler
     $              , nxo,nyo,nzo,nrg
     $              , wdsizr,wdsizo
     $              , nfileo,nproc_o,nfldr
     $              , er(lelr),nelB,nelBr
      integer wdsizr,wdsizo,er
      
      parameter(iHeaderSize=132)
      
      common /cmfi_r/ timer
      
      common /cmfi_c/ ihdr,rdcode,mfi_fname
      character*3  ihdr
      character*10 rdcode
      character*80 mfi_fname
      
      character*1  rdcode1(10)
      equivalence (rdcode,rdcode1)
      
      common /cmfi_l/ 
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps (ldimt1),ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr(ldimt1),ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      logical
     $        ifgetx ,ifgetu ,ifgetp ,ifgett ,ifgtps         ,ifgtim
     $       ,ifgetxr,ifgetur,ifgetpr,ifgettr,ifgtpsr        ,ifgtimr
     $       ,if_byte_sw
     $       ,ifgetz,ifgetw
     $       ,ifdiro
      
      common /cmfi_p/ fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00 
      integer         fid0,fid0r,pid0,pid1,pid0r,pid1r,pid00
      
      integer          nekcomm_io,ifh_mbyte
      common /i4mpiio/ nekcomm_io,ifh_mbyte
# 166 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 166 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"

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
# 167 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f" 2
# 167 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/byte_mpi.f"
      


















































      
# 219
      return
      end
