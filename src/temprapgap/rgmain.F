      PROGRAM RGMAIN
	Implicit None
#include "rgfull.inc"
#include "rgdisdif.inc"
#include "rgluco.inc"
#include "rgludat1.inc"
#include "rgpara.inc"
#include "rgpypars.inc"
#include "rgdiffr.inc"
#include "rgscq2.inc"
#include "rglq2.inc"
#include "rghsunts.inc"
#include "rghsoptn.inc"
#include "rghscuts.inc"
#include "rghsvglp.inc"
	Integer Nevent
	Common/steer1/Nevent
	Integer K1,K2,LUX

	real timeleft
      Integer Minuts
      External Minuts
	External pydata

C definitions for outputfile and so on
      integer asciiLun, asciiOut
      parameter (asciiLun=29)
      parameter (asciiOut=1)
      character outputfilename*80
      parameter (outputfilename='rapgap.txt')

      Integer I,ISEED
C---initialise ARIADNE parameters, now done via block data
C---initialise PYTHIA 6 parameters, via pythia block data 
C     initialize random number generator
      ISEED = 213123
      ISEED = Iabs(MINUTS())
      ISEED = 44788029
      LUX = 4
      K1=0
      K2=0
      CALL RLUXGO(LUX,ISEED,K1,K2)
C---initialise RAPGAP parameters
      CALL GRAINI
C-- read in parameters from file 
      Call Steer
C-- change standard parameters of RAPGAP 	
	Call rapcha
C-- change standard parameters of HERACLES	
	Call hercha
C-- change standard parameters of JETSET
	Call pytcha
	IF(IFPS.EQ.10) then
C Initialize ARIADNE
         CALL ARINIT('RAPGAP')
C-- change standard parameters of ARIADNE
	   Call aricha
      endif

C ---------------------------------------------------------------------
c     Open ascii output file
c ---------------------------------------------------------------------
       if (asciiOut.eq.1) then
          open(asciiLun, file=outputfilename)
       endif

        write(29,*)' RAPGAP EVENT FILE '
        write(29,*)'============================================'
        write(29,40)
40      format('I, ievent, genevent, subprocess, idir, idisdif, 
     &  c.s., sigma(c.s.), s, q2, y, xgam, xpr, Pt_hat, pt2_hat, 
     &  s_hat, t, x_pom, s_hat, z, x, phi, nrTracks')
        write(29,*)'============================================'
        write(29,*)' I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)
     &  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V9I,2)  V(I,3)'
        write(29,*)'============================================'

C--- CALCULATE X SECTION
      CALL PTIME(' rapgap  ',1,0)
      CALL RAPGAP
C--- print x section
      CALL RAEND(1)
C--- event generation
      DO 10 I=1,Nevent
         CALL PTIME(' event  ',1,0)
c	   call timel(timeleft)
c	   if(timeleft.le.10.) then
c	      write(6,*) ' *** time limit reached. time left = ',timeleft
c            goto 123
c	   endif
         CALL EVENT
C--- user analysis routine
          CALL ANALYS
C--- user analysis routine (hztool)
cc         CALL hzraana
      CALL PTIME(' event  ',2,0)
C---
   10 CONTINUE
123   Continue
C---PRINT NR OF GENERATED EVENTS
      CALL RAEND(20)	
C--- user terminate routine (hztool)
cc	CALL hzraend
	
      CALL PTIME(' rapgap  ',2,0)
      CALL PTIME('        ',2,99)
      STOP
      END
C#include "pythia-bkdata/pydata.f"
