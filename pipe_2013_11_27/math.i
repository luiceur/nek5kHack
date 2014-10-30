# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
c-----------------------------------------------------------------------
      SUBROUTINE BLANK(A,N)
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
C     
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE VSQ (A,N)
      DIMENSION  A(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 18 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 18 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 20
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'vsq   '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 32
      DO 100 I = 1, N
 100     A(I) = A(I)**2
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE VSQRT(A,N)
      DIMENSION  A(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 41 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 41 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 43
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'vsqrt '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 55
      DO 100 I = 1, N
 100     A(I) = SQRT(A(I))
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine invers2(a,b,n)
      REAL A(1),B(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 64 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 64 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 66
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'inver2'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 78
      DO 100 I=1,N
         A(I)=1./B(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine invcol1(a,n)
      REAL A(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 88 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 88 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 90
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl1'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 102
      DO 100 I=1,N
         A(I)=1./A(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine invcol2(a,b,n)
C     
      REAL A(1),B(1)

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
# 112 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 112 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 113 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 113 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 115
      if (icalld.eq.0) tinvc=0.0
      icalld=icalld+1
      ninvc=icalld
      etime1=dnekclock()
C     
C     
C     
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl2'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 134
      DO 100 I=1,N
         A(I)=A(I)/B(I)
 100  CONTINUE

# 138
      tinvc=tinvc+(dnekclock()-etime1)

# 140
      return
      END
c-----------------------------------------------------------------------
      subroutine invcol3(a,b,c,n)
      REAL A(1),B(1),C(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 147 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 147 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

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
# 148 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 148 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 150
      if (icalld.eq.0) tinv3=0.0
      icalld=icalld+1
      ninv3=icalld
      etime1=dnekclock()
C     
C     
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl3'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 168
      DO 100 I=1,N
         A(I)=B(I)/C(I)
 100  CONTINUE

# 172
      tinv3=tinv3+(dnekclock()-etime1)

# 174
      return
      END
c-----------------------------------------------------------------------
      subroutine col4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 181 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 181 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 183
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col4  '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 195
      DO 100 I=1,N
         A(I)=B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine Xaddcol3(a,b,c,n)
      REAL A(1),B(1),C(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 205 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 205 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 207
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl3'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 219
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine addcol4(a,b,c,d,n)
      REAL A(n),B(n),C(n),D(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 229 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 229 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 231
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl4'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 243
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ascol5 (a,b,c,d,e,n)
      REAL A(1),B(1),C(1),D(1),E(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 253 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 253 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 255
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ascol5'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 267
      DO 100 I=1,N
         A(I) = B(I)*C(I)-D(I)*E(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine sub2(a,b,n)
      REAL A(1),B(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 277 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 277 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 279
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'sub2  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 291
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine sub3(a,b,c,n)
      REAL A(1),B(1),C(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 301 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 301 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 303
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'sub3  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 315
      DO 100 I=1,N
         A(I)=B(I)-C(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine subcol3(a,b,c,n)
      REAL A(1),B(1),C(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 325 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 325 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 327
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'subcl3'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 339
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine subcol4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 349 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 349 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 351
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'subcl4'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 363
      DO 100 I=1,N
         A(I)=A(I)-B(I)*C(I)*D(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine rzero(a,n)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END
c-----------------------------------------------------------------------
      subroutine izero(a,n)
      INTEGER A(1)
C     
      DO 100 I = 1, N
 100     A(I ) = 0
      return
      END
c-----------------------------------------------------------------------
      subroutine ione(a,n)
      INTEGER   A(1)
      DO 100 I = 1, N
 100     A(I ) = 1
      return
      END
c-----------------------------------------------------------------------
      subroutine rone(a,n)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 1.0
      return
      END
c-----------------------------------------------------------------------
      subroutine cfill(a,b,n)
      DIMENSION  A(1)
C     
      DO 100 I = 1, N
 100     A(I) = B
      return
      END
c-----------------------------------------------------------------------
      subroutine ifill(ia,ib,n)
      DIMENSION IA(1)
C     
      DO 100 I = 1, N
 100     IA(I) = IB
      return
      END
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)
      
      do i=1,n
         a(i)=b(i)
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)
      CHARACTER*1 A(1), B(1)
C     
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
C     
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
C     
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine i8copy(a,b,n)
      INTEGER*8 A(1), B(1)
C     
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine chsign(a,n)
      REAL A(1)
C     
      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
      return
      END
C     
c-----------------------------------------------------------------------
      subroutine cmult(a,const,n)
      REAL A(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 462 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 462 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 464
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cmult '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 476
      DO 100 I=1,N
         A(I)=A(I)*CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine cadd(a,const,n)
      REAL A(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 486 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 486 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 488
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cadd  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 500
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine iadd(i1,iscal,n)
      DIMENSION I1(1)
C     
      DO 10 I=1,N
         I1(I)=I1(I)+ISCAL
   10 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine cadd2(a,b,const,n)
      REAL A(1),B(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 519 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 519 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 521
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cadd2 '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 533
      DO 100 I=1,N
         A(I)=B(I)+CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      real function vlmin(vec,n)
      REAL VEC(1)
      TMIN = 99.0E20
C     
      DO 100 I=1,N
         TMIN = MIN(TMIN,VEC(I))
 100  CONTINUE
      VLMIN = TMIN
      return
      END
c-----------------------------------------------------------------------
      integer function ivlmin(vec,n)
      integer vec(1),tmin
      if (n.eq.0) then
         ivlmin=0
         return
      endif
      tmin = 8888888
      do i=1,n
         tmin = min(tmin,vec(i))
      enddo
      ivlmin = tmin
      return
      end
c-----------------------------------------------------------------------
      integer function ivlmax(vec,n)
      integer vec(1),tmax
      if (n.eq.0) then
         ivlmax=0
         return
      endif
      TMAX =-8888888
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      Ivlmax = tmax
      return
      end
c-----------------------------------------------------------------------
      real function vlmax(vec,n)
      REAL VEC(1)
      TMAX =-99.0E20
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      VLMAX = TMAX
      return
      END
c-----------------------------------------------------------------------
      real function vlamax(vec,n)
      REAL VEC(1)
      TAMAX = 0.0
C     
      DO 100 I=1,N
         TAMAX = MAX(TAMAX,ABS(VEC(I)))
 100  CONTINUE
      VLAMAX = TAMAX
      return
      END
c-----------------------------------------------------------------------
      real function vlsum(vec,n)
      REAL VEC(1)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 602 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 602 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 604
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'vlsum '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 616
      SUM = 0.
C     
      DO 100 I=1,N
         SUM=SUM+VEC(I)
 100  CONTINUE
      VLSUM = SUM
      return
      END
c-----------------------------------------------------------------------
      subroutine vcross (u1,u2,u3,v1,v2,v3,w1,w2,w3,n)
C     
C     Compute a Cartesian vector cross product.
C     
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
      DIMENSION W1(1),W2(1),W3(1)
C     
C     
      DO 100 I=1,N
         U1(I) = V2(I)*W3(I) - V3(I)*W2(I)
         U2(I) = V3(I)*W1(I) - V1(I)*W3(I)
         U3(I) = V1(I)*W2(I) - V2(I)*W1(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vdot2 (dot,u1,u2,v1,v2,n)
C     
C     Compute a Cartesian vector dot product. 2-d version
C     
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1)
      DIMENSION V1(1),V2(1)
C     
C     
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) 
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vdot3 (dot,u1,u2,u3,v1,v2,v3,n)
C     
C     Compute a Cartesian vector dot product. 3-d version
C     
      DIMENSION DOT(1)
      DIMENSION U1(1),U2(1),U3(1)
      DIMENSION V1(1),V2(1),V3(1)
C     
C     
      DO 100 I=1,N
         DOT(I) = U1(I)*V1(I) + U2(I)*V2(I) + U3(I)*V3(I)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine addtnsr(s,h1,h2,h3,nx,ny,nz)
C     
C     Map and add to S a tensor product form of the three functions H1,H
C     This is a single element routine used for deforming geometry.
C     
      DIMENSION H1(1),H2(1),H3(1)
      DIMENSION S(NX,NY,NZ)
C     
      DO 200 IZ=1,NZ
      DO 200 IY=1,NY
         HH = H2(IY)*H3(IZ)
         DO 100 IX=1,NX
            S(IX,IY,IZ)=S(IX,IY,IZ)+HH*H1(IX)
  100    CONTINUE
  200 CONTINUE
      return
      END
      function ltrunc(string,l)
      CHARACTER*1 STRING(L)
      CHARACTER*1   BLNK
      DATA BLNK/' '/
C     
      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      return
      END
c-----------------------------------------------------------------------
      function mod1(i,n)
C     
C     Yields MOD(I,N) with the exception that if I=K*N, result is N.
C     
      MOD1=0
      IF (I.EQ.0) THEN
         return
      ENDIF
      IF (N.EQ.0) THEN
         WRITE(6,*) 
     $  'WARNING:  Attempt to take MOD(I,0) in function mod1.'
         return
      ENDIF
      II = I+N-1
      MOD1 = MOD(II,N)+1
      return
      END
c-----------------------------------------------------------------------
      integer function log2(k)
      RK=(K)
      RLOG=LOG10(RK)
      RLOG2=LOG10(2.0)
      RLOG=RLOG/RLOG2+0.5
      LOG2=INT(RLOG)
      return
      END
c-----------------------------------------------------------------------
      subroutine iflip(i1,n)
      DIMENSION I1(1)
      N1=N+1
      N2=N/2
      DO 10 I=1,N2
         ILAST=N1-I
         ITMP=I1(ILAST)
         I1(ILAST)=I1(I)
         I1(I)=ITMP
   10 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine iswap(b,ind,n,temp)
      INTEGER B(1),IND(1),TEMP(1)
C***  
C***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
C***  INTO ITEM(I), WHERE JJ=IND(I).
C***  
      DO 20 I=1,N
         JJ=IND(I)
         TEMP(I)=B(JJ)
   20 CONTINUE
      DO 30 I=1,N
   30 B(I)=TEMP(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine col2(a,b,n)
      real a(1),b(1)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 762 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 762 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 764
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 777
      do i=1,n
         a(i)=a(i)*b(i)
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine col2c(a,b,c,n)
      real a(1),b(1),c
      
      do i=1,n
         a(i)=a(i)*b(i)*c
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      subroutine col3(a,b,c,n)
      real a(1),b(1),c(1)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 797 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 797 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 799
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col3  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 812
      do i=1,n
         a(i)=b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2(a,b,n)
      real a(1),b(1)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 821 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 821 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 823
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ADD2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 836
      do i=1,n
         a(i)=a(i)+b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add3(a,b,c,n)
      real a(1),b(1),c(1)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 845 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 845 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 847
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ADD3  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 860
      do i=1,n
         a(i)=b(i)+c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine addcol3(a,b,c,n)
      real a(1),b(1),c(1)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 869 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 869 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 871
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl3'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 884
      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2s1(a,b,c1,n)
      real a(1),b(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 894 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 894 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 896
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s1'
      endif
      isbcnt = 2*N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 908
      DO 100 I=1,N
        A(I)=C1*A(I)+B(I)
  100 CONTINUE
      return
      END
C     
c-----------------------------------------------------------------------
      subroutine add2s2(a,b,c1,n)
      real a(1),b(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 919 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 919 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 921
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s2'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 933
      DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
  100 CONTINUE
      return
      END
C     
c-----------------------------------------------------------------------
      subroutine add3s2(a,b,c,c1,c2,n)
      real a(1),b(1),c(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 944 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 944 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 946
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add3s2'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 958
      DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)
  100 CONTINUE
      return
      END
C     
c-----------------------------------------------------------------------
      subroutine add4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 969 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 969 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 971
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add4  '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 983
      DO 100 I=1,N
         A(I)=B(I)+C(I)+D(I)
 100  CONTINUE
      return
      END
      real function vlsc2(x,y,n)
      REAL X(1),Y(1)

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
# 991 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 991 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 992 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 992 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

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
# 993 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 993 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 995
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 1007
      s = 0.
      do i=1,n
         s = s + x(i)*y(i)
      enddo
      vlsc2=s
      return
      end
c-----------------------------------------------------------------------
      real function vlsc21(x,y,n)
      real x(1),y(1)

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
# 1018 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1018 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 1019 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1019 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

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
# 1020 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1020 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 1022
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC21'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 1034
      s = 0.
      do i=1,n
         s = s + x(i)*x(i)*y(i)
      enddo
      vlsc21=s
      return
      end
      
      
C-----------------------------------------------------------------------
C     
C     Vector reduction routines which require communication 
C     on a parallel machine. These routines must be substituted with
C     appropriate routines which take into account the specific architec
C     
C-----------------------------------------------------------------------
      
      
      function glsc3(a,b,mult,n)
C     
C     Perform inner-product in double precision
C     
      REAL A(1),B(1),MULT(1)
      REAL TMP,WORK(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 1060 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1060 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 1062
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc3 '
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 1074
      TMP = 0.0
      DO 10 I=1,N
         TMP = TMP + A(I)*B(I)*MULT(I)
 10   CONTINUE
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC3 = TMP
      return
      END
c-----------------------------------------------------------------------
      function glsc2(x,y,n)
C     
C     Perform inner-product in double precision
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 1088 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1088 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
c     
      real x(1), y(1)
      real tmp,work(1)
C     

# 1093
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
# 1105
      tmp=0.0
      do 10 i=1,n
         tmp = tmp+ x(i)*y(i)
   10 continue
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC2 = TMP
      return
      END
c-----------------------------------------------------------------------
      function glsc23(x,y,z,n)
c     
C     Perform inner-product  x*x*y*z
c     
      real x(1), y(1),z(1)
      real tmp,work(1)
      
      ds = 0.0
      do 10 i=1,n
         ds=ds+x(i)*x(i)*y(i)*z(i)
   10 continue
      tmp=ds
      call gop(tmp,work,'+  ',1)
      glsc23 = tmp
      return
      end
c-----------------------------------------------------------------------
      real function gl2norm(a,n)
      

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
# 1134 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1134 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/MASS" 1
# 1
      COMMON /MASS/
     $       BM1   (LX1,LY1,LZ1,LELT),  BM2   (LX2,LY2,LZ2,LELV)
     $      ,BINVM1(LX1,LY1,LZ1,LELV),  BINTM1(LX1,LY1,LZ1,LELT)
     $      ,BM2INV(LX2,LY2,LZ2,LELT),  BAXM1 (LX1,LY1,LZ1,LELT)
     $      ,BM1LAG(LX1,LY1,LZ1,LELT,LORDER-1)
     $      ,VOLVM1,VOLVM2,VOLTM1,VOLTM2
     $      ,YINVM1(LX1,LY1,LZ1,LELT)
# 1135 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1135 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      
      real a(1)
      
      common /scrsf/ w1 (lx1,ly1,lz1,lelt)
      
      call col3 (w1,a,a,n)
      call col2 (w1,bm1,n)
      gl2norm = sqrt(glsum (w1,n)/volvm1)
      
      return
      end
c-----------------------------------------------------------------------
      function glsum (x,n)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      TSUM = 0.
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
      TMP(1)=TSUM
      CALL GOP(TMP,WORK,'+  ',1)
      GLSUM = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      real function glamax(a,n)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMAX = 0.0
      DO 100 I=1,N
         TMAX = MAX(TMAX,ABS(A(I)))
 100  CONTINUE
      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      GLAMAX=ABS(TMP(1))
      return
      END
c-----------------------------------------------------------------------
      real function glamin(a,n)
      real a(1)
      dimension tmp(1),work(1)
      tmin = 9.e28
      do 100 i=1,n
         tmin = min(tmin,abs(a(i)))
 100  continue
      tmp(1)=tmin
      call gop(tmp,work,'m  ',1)
      glamin=abs(tmp(1))
      return
      end
c-----------------------------------------------------------------------
      function iglmin(a,n)
      integer a(1),tmin
      integer tmp(1),work(1)
      tmin=  999999999
      do i=1,n
         tmin=min(tmin,a(i))
      enddo
      tmp(1)=tmin
      call igop(tmp,work,'m  ',1)
      iglmin=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      function iglmax(a,n)
      integer a(1),tmax
      integer tmp(1),work(1)
      tmax= -999999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call igop(tmp,work,'M  ',1)
      iglmax=tmp(1)
      return
      end
c-----------------------------------------------------------------------
      function iglsum(a,n)
      integer a(1),tsum
      integer tmp(1),work(1)
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call igop(tmp,work,'+  ',1)
      iglsum=tmp(1)
      return
      end
C-----------------------------------------------------------------------
      integer*8 function i8glsum(a,n)
      integer*8 a(1),tsum
      integer*8 tmp(1),work(1)
      tsum= 0
      do i=1,n
         tsum=tsum+a(i)
      enddo
      tmp(1)=tsum
      call i8gop(tmp,work,'+  ',1)
      i8glsum=tmp(1)
      return
      end
C-----------------------------------------------------------------------
      function glmax(a,n)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMAX=-99.0e20
      DO 100 I=1,N
         TMAX=MAX(TMAX,A(I))
  100 CONTINUE
      TMP(1)=TMAX
      CALL GOP(TMP,WORK,'M  ',1)
      GLMAX=TMP(1)
      return
      END
c-----------------------------------------------------------------------
      function glmin(a,n)
      REAL A(1)
      DIMENSION TMP(1),WORK(1)
      TMIN=99.0e20
      DO 100 I=1,N
         TMIN=MIN(TMIN,A(I))
  100 CONTINUE
      TMP(1)=TMIN
      CALL GOP(TMP,WORK,'m  ',1)
      GLMIN = TMP(1)
      return
      END
c-----------------------------------------------------------------------
      subroutine gllog(la,lb)
C     
C     If ANY LA=LB, then ALL LA=LB.
C     
      LOGICAL LA,LB
      DIMENSION TMP(1),WORK(1)
C     
      TMP(1)=1.0
      IF (LB) THEN
         IF (LA) TMP(1)=0.0
      ELSE
         IF (.NOT.LA) TMP(1)=0.0
      ENDIF
      CALL GOP(TMP,WORK,'*  ',1)
      IF (TMP(1).EQ.0.0) LA=LB
      return
      END
c-----------------------------------------------------------------------
      function fmdian(a,n,ifok)
C     find the Median of the (global) set A

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
# 1285 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1285 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      DIMENSION A(1)
      DIMENSION WORK1(5),WORK2(5)
      DIMENSION GUES(100)
      LOGICAL IFOK
C     
      AMP  =1.5
      AFAC =1.5
      GMIN =GLMIN(A,N)
      GMAX =GLMAX(A,N)
      GMIN0=GLMIN(A,N)
      GMAX0=GLMAX(A,N)
      GUESS=(GMAX+GMIN)/2.0
      EPS  =(GMAX-GMIN)
      IF (EPS.EQ.0.0) THEN
         FMDIAN=GMAX
         return
      ENDIF
      WORK1(1)=N
      CALL GOP(WORK1,WORK2,'+  ',1)
      NTOT=WORK1(1)
      N2 = (NTOT+1)/2
      IF (.NOT.IFOK) THEN
        WRITE(6,8) NID,N,(A(I),I=1,N)
        WRITE(6,9) NID,NTOT,N2,N,GMIN,GMAX
    8   FORMAT(I5,'N,A:',I5,10(6F10.5,/)) 
    9   FORMAT(I5,'mnx:',3I6,2F10.5)
      ENDIF
C     
C     This is the trial loop
C     
      ITRY=-1
   10 CONTINUE
      ITRY=ITRY+1
      II=ITRY+1
      IF (II.LE.100) GUES(II)=GUESS
C     error check for infinite loop
      IF (ITRY.GT.2*NTOT) GOTO 9000
      CALL RZERO(WORK1,5)
      NLT=0
      NGT=0
      CLT=GMIN0
      CGT=GMAX0
      DO 100 I=1,N
         AA=A(I)
         IF (AA.NE.GUESS) THEN
            IF (AA.LT.GUESS) THEN
               NLT=NLT+1
C              CLT - closest value to GUESS Less Than GUESS
               IF (AA.GT.CLT) CLT=AA
            ENDIF
            IF (AA.GT.GUESS) THEN
               NGT=NGT+1
C              CGT - closest value to GUESS Greater Than GUESS
               IF (AA.LT.CGT) CGT=AA
            ENDIF
            DUM=1./(EPS+ABS(AA-GUESS))
            WORK1(1)=WORK1(1)+DUM
            WORK1(2)=WORK1(2)+DUM*AA
         ELSE
C           detected values equaling the guess.
            WORK1(5)=WORK1(5)+1.0
         ENDIF
  100 CONTINUE
C     Invoke vector reduction across processors:
      WORK2(1)=CLT
      CLT=GLMAX(WORK2,1)
      WORK2(1)=CGT
      CGT=GLMIN(WORK2,1)
      WORK1(3)=NLT
      WORK1(4)=NGT
      CALL GOP(WORK1,WORK2,'+  ',5)
      NLT=WORK1(3)
      NGT=WORK1(4)
      IF (.NOT.IFOK) THEN
         WRITE(6,101) NID,GUESS,CLT,CGT
         WRITE(6,102) NID,(WORK1(I),I=1,5)
  101    FORMAT(I5,'Glg:',3F12.5)
  102    FORMAT(I5,'WORK1:',5F12.5)
      ENDIF
C     
C     Done?
C     
      IF (NLT.GT.N2.OR.NGT.GT.N2) THEN
C        we're not done.....
         IF (NGT.GT.NLT) THEN
C           guess is too low
            GMIN=CGT
            G2=CGT+MAX(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2.GT.GMAX) G2=0.5*(GUESS+GMAX)
            EPS=AFAC*ABS(G2-GUESS)
C           see that we move at least as far as the next closest value.
            GUESS=MAX(G2,CGT)
            GOTO 10
         ELSE IF (NLT.GT.NGT) THEN
C           guess is too high
            GMAX=CLT
            G2=CLT+MIN(0.,WORK1(2)/WORK1(1)-GUESS)*AMP
            IF (G2.LT.GMIN) G2=0.5*(GUESS+GMIN)
            EPS=AFAC*ABS(G2-GUESS)
C           see that we move at least as far as the next closest value.
            GUESS=MIN(G2,CLT)
            GOTO 10
         ENDIF
      ELSE
C     
C        we're done....
         IF (WORK1(5).NE.0) THEN
C           the median is (usually) one of the values 
            FMDIAN=GUESS
            IF (WORK1(5).EQ.1.0) THEN
               IF (MOD(NTOT,2).EQ.0) THEN
                  IF (NGT.GT.NLT) THEN
                     FMDIAN=0.5*(GUESS+CGT)
                  ELSE
                     FMDIAN=0.5*(GUESS+CLT)
                  ENDIF
               ELSE
                  IF (NGT.EQ.NLT) THEN
                     FMDIAN=GUESS
                  ELSE IF(NGT.GT.NLT) THEN
                     FMDIAN=CGT
                  ELSE
                     FMDIAN=CLT
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            IF (MOD(NTOT,2).EQ.0) THEN
               IF (NGT.EQ.NLT) THEN
                  FMDIAN=0.5*(CLT+CGT)
               ELSE IF(NGT.GT.NLT) THEN
                  FMDIAN=0.5*(GUESS+CGT)
               ELSE
                  FMDIAN=0.5*(GUESS+CLT)
               ENDIF
            ELSE
               IF (NGT.EQ.NLT) THEN
                  FMDIAN=GUESS
               ELSE IF(NGT.GT.NLT) THEN
                  FMDIAN=CGT
               ELSE
                  FMDIAN=CLT
               ENDIF
           ENDIF
         ENDIF
C     
      ENDIF
       IF (.NOT.IFOK) WRITE(6,*) NID,'FMDIAN2',FMDIAN,(A(I),I=1,N)
      return
C     
C     Error handling
C     
 9000 CONTINUE
      WRITE(6,11) NTOT,GMIN0,GMAX0,GUESS
   11 FORMAT('ABORTING IN FMDIAN: N,AMIN,AMAX:',I6,3G14.6)
      DO 13 I1=1,N,5
        IN=I1+5 
        IN=MIN(IN,N)
        WRITE(6,12) NID,(A(I),I=I1,IN)
   12   FORMAT(I4,' FMA:',5G14.6)
   13 CONTINUE
      DO 15 I1=1,ITRY,5
        IN=I1+5
        IN=MIN(IN,ITRY)
        WRITE(6,14) NID,(GUES(I),I=I1,IN)
   14   FORMAT(I4,' FMG:',5G14.6)
   15 CONTINUE
      call exitt
      END
      
C=======================================================================
C     Double precision matrix and vector routines
C=======================================================================
      
c-----------------------------------------------------------------------
      subroutine dcadd(a,const,n)
      real*8 A(1),CONST
C     
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine dsub2(a,b,n)
      real*8 A(1), B(1)
C     
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
      return
      END
C     
c-----------------------------------------------------------------------
      subroutine dadd2(a,b,n)
      real*8 A(1), B(1)
C     
      DO 100 I=1,N
         A(I)=A(I)+B(I)
 100  CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine chswapr(b,L,ind,n,temp)
      INTEGER IND(1)
      CHARACTER*6 B(1),TEMP(1)
C***  
C***  SORT ASSOCIATED ELEMENTS BY PUTTING ITEM(JJ)
C***  INTO ITEM(I), WHERE JJ=IND(I).
C***  
      DO 20 I=1,N
         JJ=IND(I)
         TEMP(I)=B(JJ)
   20 CONTINUE
      DO 30 I=1,N
   30 B(I)=TEMP(I)
      return
      END
c-----------------------------------------------------------------------
      subroutine drcopy(r,d,N)
      real*8    d(1)
      dimension r(1)
      do 10 i=1,n
         r(i)=d(i)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine rrcopy(r,d,N)
      real*4 d(1)
      real*4 r(1)
      do 10 i=1,n
         r(i)=d(i)
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine sorts(xout,xin,work,n)
      real xout(1),xin(1),work(1)
      call copy(xout,xin,n)
      call sort(xout,work,n)
      return
      end
C     
c-----------------------------------------------------------------------
      function ivlsum(a,n)
      INTEGER A(1)
      INTEGER TSUM
      if (n.eq.0) then
         ivlsum = 0
         return
      endif
      TSUM=A(1)
      DO 100 I=2,N
         TSUM=TSUM+A(I)
  100 CONTINUE
      IVLSUM=TSUM
      return
      END
c-----------------------------------------------------------------------
      subroutine icadd(a,c,n)
      INTEGER A(1),C
      DO 100 I = 1, N
 100     A(I) = A(I) + C
      return
      END
      subroutine isort(a,ind,n)
C     
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C     
      integer a(1),ind(1)
      integer aa
C     
      dO 10 j=1,n
         ind(j)=j
   10 continue
C     
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
      subroutine sort(a,ind,n)
C     
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C     
      real a(1),aa
      integer ind(1)
C     
      dO 10 j=1,n
         ind(j)=j
   10 continue
C     
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      subroutine iswap_ip(x,p,n)
      integer x(1),xstart
      integer p(1)
c     
c     In-place permutation: x' = x(p)
c     
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! iswap_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  x(last) = x(next)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
c     
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine iswapt_ip(x,p,n)
      integer x(1),t1,t2
      integer p(1)
c     
c     In-place permutation: x'(p) = x
c     
      
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
               else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c     
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine swap_ip(x,p,n)
      real    x(1),xstart
      integer p(1)
c     
c     In-place permutation: x' = x(p)
c     
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! swap_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  x(last) = x(next)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
c     
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine swapt_ip(x,p,n)
      real    x(1),t1,t2
      integer p(1)
c     
c     In-place permutation: x'(p) = x
c     
      
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! swapt_ip problem.',j,k,n,next
                  call exitt
               elseif (next.eq.loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
               else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c     
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine glvadd(x,w,n)
      real x(1),w(1)
      call gop(x,w,'+  ',n)
      return
      end
c-----------------------------------------------------------------------
      subroutine add3s12(x,y,z,c1,c2,n)
      real x(1),y(1),z(1),c1,c2
      do i=1,n
         x(i) = c1*y(i)+c2*z(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      integer*8 function i8glmax(a,n)
      integer*8 a(1),tmax
      integer*8 tmp(1),work(1)
      tmax= -999999
      do i=1,n
         tmax=max(tmax,a(i))
      enddo
      tmp(1)=tmax
      call i8gop(tmp,work,'M  ',1)
      i8glmax=tmp(1)
      if (i8glmax .eq. -999999) i8glmax=0
      return
      end
c-----------------------------------------------------------------------
      subroutine admcol3(a,b,c,d,n)
      REAL A(1),B(1),C(1),D
C     
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D
  100 CONTINUE
      return
      END
      
c-----------------------------------------------------------------------
      subroutine add2col2(a,b,c,n)
      real a(1),b(1),c(1)
c     
      do i=1,n
         a(i) = a(i) + b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2sxy(x,a,y,b,n)
      real x(1),y(1)
c     
      do i=1,n
         x(i) = a*x(i) + b*y(i)
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine col2s2(x,y,s,n)
      real x(n),y(n)
c     
      do i=1,n
         x(i)=s*x(i)*y(i)
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      
      
      
      
C====================== OPENACC ============
      

      
# 1881
      subroutine i8copy_acc(a,b,n)
      INTEGER*8 A(n), B(n)
C     
!$ACC DATA PRESENT(a,b)
!$ACC PARALLEL LOOP
      DO 100 I = 1, N
 100          A(I) = B(I)
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      subroutine rzero_acc(a,n)
      DIMENSION  A(n)
      
!$ACC DATA PRESENT(a(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(64)
      DO 100 I = 1, N
 100     A(I ) = 0.0
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
      
c-----------------------------------------------------------------------
      subroutine rone_acc(a,n)
      DIMENSION  A(n)
      
!$ACC DATA PRESENT(a(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(64)
      DO 100 I = 1, N
 100     A(I ) = 1.0
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
      
c-----------------------------------------------------------------------
      subroutine copy_acc(a,b,n)
      real a(n),b(n)
      
!$ACC DATA PRESENT(a(1:n),b(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(64)
      do i=1,n
         a(i)=b(i)
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      end
      
      
c-----------------------------------------------------------------------
      subroutine invers2_acc(a,b,n)
      REAL A(n),B(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 1945 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1945 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 1947
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'inver2'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
      
# 1960
!$ACC DATA PRESENT(a(1:n),b(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      DO 100 I=1,N
         A(I)=1./B(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      subroutine add2_acc(a,b,n)
      real a(n),b(n)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 1975 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 1975 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 1977
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ADD2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 1990
!$ACC DATA PRESENT(a(1:n),b(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      do i=1,n
         a(i)=a(i)+b(i)
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      end
      
c-----------------------------------------------------------------------
      subroutine col2_acc(a,b,n)
      real a(n),b(n)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2005 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2005 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 2007
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 2020
!$ACC DATA PRESENT(a(1:n),b(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
       do i=1,n
         a(i)=a(i)*b(i)
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA
      return
      end
      
c-----------------------------------------------------------------------
      subroutine invcol2_acc(a,b,n)
C     
      REAL A(n),B(n)

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
# 2035 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2035 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2036 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2036 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2038
      if (icalld.eq.0) tinvc=0.0
      icalld=icalld+1
      ninvc=icalld
      etime1=dnekclock()
C     
C     
C     
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl2'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2057
!$ACC DATA PRESENT(a(1:n),b(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      DO 100 I=1,N
         A(I)=A(I)/B(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      

# 2066
      tinvc=tinvc+(dnekclock()-etime1)

# 2068
      return
      END
c-----------------------------------------------------------------------
      
      
c-----------------------------------------------------------------------
      subroutine addcol3_acc(a,b,c,n)
      real a(n),b(n),c(n)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2077 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2077 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 2079
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl3'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
# 2092
!$ACC DATA PRESENT(a(1:n),b(1:n),c(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo
!$ACC END DATA
      return
      end
      
      
c-----------------------------------------------------------------------
      subroutine cmult_acc(a,const,n)
      REAL A(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2107 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2107 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2109
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cmult '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2121
!$ACC DATA PRESENT(a(1:n)) 
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      DO 100 I=1,N
         A(I)=A(I)*CONST
 100  CONTINUE
!$ACC END PARALLEL LOOP 
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      subroutine add2col2_acc(a,b,c,n)
      real a(n),b(n),c(n)
c     
!$ACC DATA PRESENT(a(1:n),b(1:n),c(1:n)) 
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      do i=1,n
         a(i) = a(i) + b(i)*c(i)
      enddo
!$ACC END PARALLEL LOOP 
!$ACC END DATA
      
      return
      end
      
c-----------------------------------------------------------------------
      subroutine cadd_acc(a,const,n)
      REAL A(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2152 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2152 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2154
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'cadd  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2166
!$ACC DATA PRESENT(a(1:n)) 
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      DO 100 I=1,N
         A(I)=A(I)+CONST
 100  CONTINUE
!$ACC END PARALLEL LOOP 
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      subroutine col3_acc(a,b,c,n)
      real a(n),b(n),c(n)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2181 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2181 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 2183
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col3  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
!xbm* unroll (10)
C     
# 2197
!$ACC DATA PRESENT(a(1:n),b(1:n),c(1:n)) 
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      do i=1,n
         a(i)=b(i)*c(i)
      enddo
!$ACC END PARALLEL LOOP 
!$ACC END DATA
      
      return
      end
      
c-----------------------------------------------------------------------
      subroutine add2s2_acc(a,b,c1,n)
      real a(n),b(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2213 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2213 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2215
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s2'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
C     
# 2228
!$ACC DATA PRESENT(a(1:n),b(1:n)) 
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
  100 CONTINUE
!$ACC END PARALLEL LOOP 
!$ACC END DATA
      return
      END
C     
      
c-----------------------------------------------------------------------
      real function vlsc2_acc(x,y,n)
      REAL X(n),Y(n)

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
# 2243 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2243 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2244 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2244 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

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
# 2245 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2245 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2247
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2259
      s = 0.
      
!$ACC DATA PRESENT(x(1:n),y(1:n))
!$ACC PARALLEL LOOP REDUCTION(+:s) 
      do i=1,n
         s = s + x(i)*y(i)
      enddo
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      vlsc2_acc=s
      return
      end
c-----------------------------------------------------------------------
      
c-----------------------------------------------------------------------
      function glsum_acc(x,n)
      DIMENSION X(1)
      DIMENSION TMP(1),WORK(1)
      TSUM = 0.
      
!$ACC DATA PRESENT(x(1:n))
!$ACC PARALLEL LOOP REDUCTION(+:TSUM) 
      DO 100 I=1,N
         TSUM = TSUM+X(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      TMP(1)=TSUM
      CALL GOP(TMP,WORK,'+  ',1)
      GLSUM_ACC = TMP(1)
      return
      END
      
      
c-----------------------------------------------------------------------
      function glsc2_acc(x,y,n)
C     
C     Perform inner-product in double precision
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2301 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2301 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
c     
      real x(n), y(n)
      real tmp,work(1)
C     

# 2306
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

      
# 2318
      tmp=0.0
      
!$ACC DATA PRESENT(x(1:n),y(1:n))
!$ACC PARALLEL LOOP REDUCTION(+:tmp) 
      do 10 i=1,n
         tmp = tmp+ x(i)*y(i)
   10 continue
!$ACC END PARALLEL LOOP
!$ACC END DATA
      CALL GOP(TMP,WORK,'+  ',1)
      
      GLSC2_ACC = TMP
      
      return
      END
      
      
C-----------------------------------------------------------------------
C     
C     Vector reduction routines which require communication 
C     on a parallel machine. These routines must be substituted with
C     appropriate routines which take into account the specific architec
C     
C-----------------------------------------------------------------------
      
      function glsc3_acc(a,b,mult,n)
C     
C     Perform inner-product in double precision
C     
      REAL A(n),B(n),MULT(n)
      REAL TMP,WORK(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2351 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2351 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2353
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc3 '
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2365
      TMP = 0.0
!$ACC DATA PRESENT(A(1:n),B(1:n),MULT(1:n))
!$ACC PARALLEL LOOP REDUCTION(+:TMP) 
      DO 10 I=1,N
         TMP = TMP + A(I)*B(I)*MULT(I)
 10   CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC3_ACC = TMP
      return
      END
      
c-----------------------------------------------------------------------
      subroutine sub2_acc(a,b,n)
      REAL A(n),B(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2384 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2384 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2386
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'sub2  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2398
!$ACC DATA PRESENT(A(1:n),B(1:n))
!$ACC PARALLEL LOOP VECTOR_LENGTH(128)
      DO 100 I=1,N
         A(I)=A(I)-B(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      return
      END
      
c-----------------------------------------------------------------------
      subroutine add3s2_acc(a,b,c,c1,c2,n)
      real a(n),b(n),c(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2413 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2413 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2415
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add3s2'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2427
!$ACC DATA PRESENT(A,B,C)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)
  100 CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      subroutine addcol4_acc(a,b,c,d,n)
      REAL A(n),B(n),C(n),D(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2443 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2443 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2445
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl4'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2457
!$ACC DATA PRESENT(A,B,C,D)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D(I)
  100 CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      subroutine chsign_acc(a,n)
      REAL A(N)
C     
!$ACC DATA PRESENT(A)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I) = -A(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
C     
      
c-----------------------------------------------------------------------
      subroutine add2s1_acc(a,b,c1,n)
      real a(n),b(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2489 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2489 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2491
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s1'
      endif
      isbcnt = 2*N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2503
!$ACC DATA PRESENT(a,b)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
        A(I)=C1*A(I)+B(I)
  100 CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      subroutine sub3_acc(a,b,c,n)
      REAL A(n),B(n),C(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2519 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2519 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2521
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'sub3  '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
      
# 2534
!$ACC DATA PRESENT(a,b,c)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I)=B(I)-C(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
      
c-----------------------------------------------------------------------
      subroutine invcol3_acc(a,b,c,n)
      REAL A(n),B(n),C(n)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2551 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2551 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"

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
# 2552 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2552 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
      

# 2554
      if (icalld.eq.0) tinv3=0.0
      icalld=icalld+1
      ninv3=icalld
      etime1=dnekclock()
C     
C     
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl3'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
      
# 2573
!$ACC DATA PRESENT(a,b,c)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I)=B(I)/C(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      

# 2582
      tinv3=tinv3+(dnekclock()-etime1)

# 2584
      return
      END
      
c-----------------------------------------------------------------------
      subroutine admcol3_acc(a,b,c,d,n)
      REAL A(n),B(n),C(n),D
C     
!$ACC DATA PRESENT(a,b,c)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D
  100 CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
      
c-----------------------------------------------------------------------
      subroutine invcol1_acc(a,n)
      REAL A(1)
C     

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2608 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2608 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2610
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'invcl1'
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2622
!$ACC DATA PRESENT(a)
!$ACC PARALLEL LOOP
      DO 100 I=1,N
         A(I)=1./A(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      return
      END
      
c-----------------------------------------------------------------------
      real function vlsum_acc(vec,n)
      REAL VEC(n)

# 1 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/OPCTR" 1
C     OPCTR is a set of arrays for tracking the number of operations,
C     and number of calls for a particular subroutine
      
# 4
      PARAMETER (MAXRTS=1000)
      COMMON /OPCTRC/ rname(MAXRTS)
      character*6     rname
C     
      COMMON /OPCTRD/ dct(MAXRTS),rct(MAXRTS),dcount
      real*8 dcount,dct,rct
C     
      COMMON /OPCTRI/ ncall(MAXRTS),nrout
C     
      integer myrout,isclld
      save    myrout,isclld
      data    myrout /0/
      data    isclld /0/
# 2637 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f" 2
# 2637 "/lustre/atlas/scratch/csep22/trn001/nek5_acc/nek_2013_11_27/math.f"
C     

# 2639
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'vlsum '
      endif
      isbcnt = n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)

C     
# 2651
      SUM = 0.
C     
!$ACC DATA PRESENT(vec)
!$ACC PARALLEL LOOP REDUCTION(+:SUM)
      DO 100 I=1,N
         SUM=SUM+VEC(I)
 100  CONTINUE
!$ACC END PARALLEL LOOP
!$ACC END DATA
      
      VLSUM_ACC = SUM
      return
      END
      
      

      
      
      
      
      
      
