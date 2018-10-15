      SUBROUTINE EVENT
      IMPLICIT NONE
#include "rglujets.inc"
#include "rgluco.inc"
#include "rgpart.inc"     
#include "rgpara1.inc"     
#include "rgherup.inc"
      logical first
	data first/.true./
      if(ILHA.ge.1) then 
	   if(first) then 
            IF(IABS(KE).NE.2212) call eventep
            IF(IABS(KE).EQ.2212) call eventpp
	      call rgupinit
		first=.false.
         endif
      endif
      IF(IABS(KE).NE.2212) call eventep
      IF(IABS(KE).EQ.2212) call eventpp
             
      if(ILHA.ge.1) then 
         CALL PYROBO(0,N,0.d0,0.d0,-CM(1)/CM(4),-CM(2)/CM(4),
     &      -CM(3)/CM(4))
c         call pylist(1)
	   call rgupevnt
         CALL PYROBO(0,N,0.d0,0.d0,CM(1)/CM(4),CM(2)/CM(4),CM(3)/CM(4))
	endif
      return
	end
