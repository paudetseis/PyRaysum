! Copyright 2020 Andrew Frederiksen and Pascal Audet

! This file is part of PyRaysum.

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modIFy, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!===========================================================================
!
! MODULE conf
!
! Configuration module that contains global variables used in rmat and plane 
! modules to interface with the Python codes.
!
!===========================================================================

      MODULE conf

      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793d0
      INTEGER, PARAMETER :: nlaymx = 30
!
! Model parameters
!
      DOUBLE PRECISION :: a(3,3,3,3,nlaymx), thickn(nlaymx)
      DOUBLE PRECISION :: rho(nlaymx)
      INTEGER :: isoflg(nlaymx)
!
! Wavefield parameters
!
      DOUBLE PRECISION :: dt, slow, baz
!
! Eigen values and eigen vectors
!
      DOUBLE COMPLEX :: evals(6,nlaymx), evecs(6,6,nlaymx)
!
! R/T matrices
!
      DOUBLE COMPLEX :: tui(3,3,nlaymx), rui(3,3,nlaymx), &
                        tdi(3,3,nlaymx), rdi(3,3,nlaymx)

      END MODULE conf


!===========================================================================
!
! MODULE raysum
!
! Contains subroutines get_arrivals and raysum that compute R/T matrices
! for stacks of generally anisotropic layers.
!
!===========================================================================

      MODULE raysum

      CONTAINS

!---------------------------------------------------------------------------
! Subroutine get_arrivals
!
! Andrew Frederiksen
!---------------------------------------------------------------------------

      SUBROUTINE get_arrivals()

      USE conf

      IMPLICIT NONE

      DOUBLE PRECISION :: 
      INTEGER :: ntr, numph

        DO  itr = 1, ntr
          
          bailout = .true.

          DO iph = 1, numph

            CALL calc_raysum(travel_time(iph,itr),amplitude(1,iph,itr),
     &                 thick_shift,ar_list,rot,
     &                 baz(itr),slow(itr),phaselist(1,1,iph),nseg(iph),
     &                 nlay,amp_in,fchg,bailout)
          END DO 
      
        END DO

      END SUBROUTINE get_arrivals


!---------------------------------------------------------------------------
! Subroutine calc_raysum
!
! Andrew Frederiksen
!
! Main arrival-calculation routine. See raysum2.m
! Update -- keeps complex matrices as long as possible, for
! accuracy. p and amplitude are kept real -- complex values
! make no sense in ray theory.
!---------------------------------------------------------------------------

      SUBROUTINE calc_raysum()
     
      USE conf

      IMPLICIT NONE
        
!
! Interface variables:
!
      DOUBLE PRECISION tt,amp(3)
      DOUBLE PRECISION ar_list(3,3,3,3,maxlay,2),rot(3,3,maxlay)
      DOUBLE PRECISION amp_in
      INTEGER phase(maxseg,2),nseg,nlay,fchg
      logical isoflag(maxlay),bailout
!
! Scratch variables:
!
      INTEGER seg, lay1, lay2, phase1, phase2, laytop, laybot, rnum, mult
      DOUBLE PRECISION p1r(3), p2r(3), Rt(3,3)
      DOUBLE COMPLEX evalbot(6), evecbot(6,6), evaltop(6), evectop(6,6)
      DOUBLE COMPLEX eval1(6), evec1(6,6), eval2(6)
      DOUBLE COMPLEX evec2(6,6), evecin_r(3)
      DOUBLE COMPLEX QQ(6,6),MM(3,3),Nu(3,3),Nd(3,3),wrk3(3,3),wrk6(6,6)
      LOGICAL upflag,rflag,fsflag,errflag
!       
! Iteration-to-iteration storage variables 
!
      DOUBLE PRECISION amp_list(maxseg),p(3,maxseg),tt_list(maxseg)
      DOUBLE COMPLEX evecin_list(3,maxseg)
      save amp_list,p,tt_list,evecin_list
      
      bailout = .false.

!
! Initialize variables, get incident slowness vector
! Assume half-space is isotropic -- save some grief.
!
        IF (fchg .eq. 0) THEN      
          tt_list(1) = 0
          amp_list(1) = amp_in
          p(1,1) = -slow*cos(baz)
          p(2,1) = -slow*sin(baz)
          CALL isotroc(aa(1,1,1,1,nlay),rho(nlay),p(1,1),
     &                 p(2,1),evalbot,evecbot)
          phase1 = mod(phase(1,2)+2,6)+1
          p(3,1) = real(evalbot(phase1))
!        
! Store active eigenvector for later sign check.
!
          CALL cextractvec(evecbot,evecin_list(1,1),1,3,phase1,6,6)
          
        END IF
       
        DO seg = max((fchg-1), 1), (nseg-1)
        
!       
! 1 is incident, 2 is transmitted/reflected 
!
          lay1 = phase(seg, 1)
          lay2 = phase(seg+1, 1)
!
! First and second phase, converted to eigenvector index
!
          phase1 = mod(phase(seg, 2) + 2, 6)+1
          phase2 = mod(phase(seg+1, 2) + 2, 6)+1
!
! Distinguish upper and lower layers.
!
          laytop = min(lay1,lay2)
          laybot = max(lay1,lay2)
!
! Upflag: true if the incident wave is upgoing.    
!
          upflag = phase1 .gt. 3
!
! Rflag: true for a reflection.
!
          rflag = lay1 .eq. lay2
          
          IF (rflag) THEN
            IF (upflag) THEN
              laytop = laybot-1
            ELSE
              laybot = laytop+1
            END IF
          END IF
!
! Fsflag: free-surface reflection 
!
          fsflag = laytop .eq. 0
          IF (fsflag .and. (.not. rflag)) THEN
            write (*,*) 'ERROR: can''t have free-surface transmission.'
          END IF 
                   
!
! laytop and laybot now contain correct layers for R/T calculation.
! Find correct rotator, determining active interface:
!
          IF (upflag) THEN
            rnum = lay1
          ELSE
            rnum = lay1+1
          END IF
          CALL rtransp(rot(1,1,rnum), Rt, 3, 3)
!
! Project slowness onto interface
!
          CALL rmatvec3(Rt,p(1,seg),p1r)
!
! IF the slowness vector is incorrectly oriented wrt the interface,
! we have a trapped phase.
!
          IF ((upflag .and. (p1r(3) .gt. 0)) .or.
     &        ((.not. upflag) .and. (p1r(3) .lt. 0))) THEN
            amp(1) = 0
            amp(2) = 0
            amp(3) = 0
            tt = 0
            tt_list(seg+1) = 0
            amp_list(seg+1) = 0
            bailout = .true.
            RETURN
          END IF          
!     
! Find evals/evecs. This is where it gets hairy...
!
! Lower eigenvalues/vectors
!
          IF (isoflag(laybot)) THEN
            CALL isotroc(aa(1,1,1,1,laybot), rho(laybot),
     &              p1r(1), p1r(2), evalbot, evecbot)
          ELSE
!
! Upper interface of laybot
!
            CALL anisotroc(ar_list(1,1,1,1,laybot,2), rho(laybot),
     &              p1r(1), p1r(2), evalbot, evecbot)
          END IF
!        
! Handle degeneracy (lower)
!
          IF (abs((evalbot(2)-evalbot(3))/evalbot(2)) .le. ztol) THEN
            CALL rot_evec(evecbot, rot(1,1,rnum))
          END IF
!
! Upper eigenvalues/vectors. Don't bother if we're at the
! free surface
!
          IF (.not.fsflag) THEN
            IF (isoflag(laytop)) THEN
              CALL isotroc(aa(1,1,1,1,laytop), rho(laytop),
     &           p1r(1), p1r(2), evaltop, evectop)
            ELSE
!
! Lower interface of laytop
!
              CALL anisotroc(ar_list(1,1,1,1,laytop,1), rho(laytop),
     &              p1r(1), p1r(2), evaltop, evectop)
            END IF
!
! Handle degeneracy (upper)             
!
            IF ((evaltop(2) .ne. 0.) .and. 
     &         (abs((evaltop(2)-evaltop(3))/evaltop(2)) .le. ztol)) THEN
              CALL rot_evec(evectop, rot(1,1,rnum))
            END IF
          END IF                  
!         
! Distribute evals/evecs between first and second propagation legs:
!
          IF (rflag) THEN
!
! Reflection
!
            IF (upflag) THEN
!
! Upgoing wave reflects off top layer -- incl. free-surf.
!
              CALL cveccopy(evalbot,eval1,6)
              CALL cmatcopy(evecbot,evec1,6,6)
              CALL cveccopy(eval1,eval2,6)
              CALL cmatcopy(evec1,evec2,6,6)
            ELSE
!
! Downgoing wave reflects off bottom layer.
!
              CALL cveccopy(evaltop,eval1,6)
              CALL cmatcopy(evectop,evec1,6,6)
              CALL cveccopy(eval1,eval2,6)
              CALL cmatcopy(evec1,evec2,6,6)
            END IF
          ELSE
!
! Transmission
!
            IF (upflag) THEN
!
! Upgoing wave enters top layer
!
              CALL cveccopy(evalbot,eval1,6)
              CALL cmatcopy(evecbot,evec1,6,6)
              CALL cveccopy(evaltop,eval2,6)
              CALL cmatcopy(evectop,evec2,6,6)
            ELSE
!
! Downgoing wave enters bottom layer
!
              CALL cveccopy(evaltop,eval1,6)
              CALL cmatcopy(evectop,evec1,6,6)
              CALL cveccopy(evalbot,eval2,6)
              CALL cmatcopy(evecbot,evec2,6,6)
            END IF
          END IF
!
!  If the eigenvalue is zero, we didn't get a transmission -- bail out
!
          IF (abs(real(eval2(phase2))) .lt. ztol**2) THEN
            amp(1) = 0
            amp(2) = 0
            amp(3) = 0
            tt = 0
            tt_list(seg+1) = 0
            amp_list(seg+1) = 0
            bailout = .true.
            RETURN
          END IF
!
! Consistency check: are the incident phase and sign correct?
! Consistency may regained by altering phase1 or using a sign
! multiplier (mult)
!
          CALL rcmatvec3(Rt, evecin_list(1,seg), evecin_r)
          CALL evec_check(evec1, evecin_r, phase1, mult, errflag)
          IF (errflag) THEN
            write (*,*) 'evaltop:',evaltop
            write (*,*) 'evalbot:',evalbot
          END IF
!       
! Find new slowness vector. By Snell's law, interface-parallel
! components carry over. Eigenvalue provides remaining component.
!
          p2r(1) = p1r(1)
          p2r(2) = p1r(2)
          p2r(3) = real(eval2(phase2))
          CALL rmatvec3(rot(1,1,rnum), p2r, p(1,seg+1))
!
! Now we're ready to calculate the amplitude. Two cases to consider.
!
          IF (fsflag) THEN   
!       
! Free-surface multiple
! Partition N: Nu=evecbot(4:6,4:6), Nd=evecbot(4:6,1:3)
!
            CALL cextractblock(evecbot,Nu,4,6,4,6,6,6)
            CALL cextractblock(evecbot,Nd,4,6,1,3,6,6)
!
! Free-surface reflection matrix: MM=-Nd^(-1)*Nu
!
            CALL cmatinv3(Nd,wrk3)
            CALL cmatmul3(wrk3,Nu,MM)
            CALL cmatconst3(MM,-1.)           
          ELSE
!
! Layer interaction
! scattering matrix QQ=evecbot^(-1)*evectop
!
            CALL eiginv(evecbot,wrk6)
            CALL cmatmul(wrk6,evectop,QQ,6,6,6)
!
!            Pull out reflection or transmission matrix MM
!
            IF (rflag) THEN
              IF (upflag) THEN
                CALL cextractblock(QQ,Nd,1,3,4,6,6,6)
                CALL cextractblock(QQ,Nu,4,6,4,6,6,6)
                CALL cmatinv3(Nu,wrk3)
                CALL cmatmul3(Nd,wrk3,MM)
              ELSE
                CALL cextractblock(QQ,Nu,4,6,4,6,6,6)
                CALL cextractblock(QQ,Nd,4,6,1,3,6,6)
                CALL cmatinv3(Nu,wrk3)
                CALL cmatmul3(wrk3,Nd,MM)
                CALL cmatconst3(MM,-1.)
              END IF
            ELSE
              IF (upflag) THEN
                CALL cextractblock(QQ,wrk3,4,6,4,6,6,6)
                CALL cmatinv3(wrk3,MM)
              ELSE
                CALL cextractblock(QQ,Nu,4,6,4,6,6,6)
                CALL cmatinv3(Nu,wrk3)
                CALL cextractblock(QQ,Nu,4,6,1,3,6,6)
                CALL cmatmul3(wrk3,Nu,Nd)
                CALL cextractblock(QQ,Nu,1,3,4,6,6,6)
                CALL cmatmul3(Nu,Nd,wrk3)
                CALL cextractblock(QQ,Nu,1,3,1,3,6,6)
                CALL cmatlincomb3(cmplx(1.),Nu,cmplx(-1.),wrk3,MM)
              END IF
            END IF            
          END IF
!          
! Carry amplitude across interface
!
          amp_list(seg+1) = amp_list(seg)*real(mult)*
     &       real(MM(mod(phase2-1,3)+1,mod(phase1-1,3)+1))
!     
!       Add post-interface segment to travel time
!
          IF (thick(lay2) .lt. 0.) THEN
            write (*,*) 'ERROR: negative thickness!'
          END IF
          tt_list(seg+1) = abs(p(3,seg+1))*thick(lay2)
!
! Pull out relevant eigenvector for next-segment consistency
! check. evecin=evec_out(:,seg+1)=R*evec2(1:3,phase2)
!
          CALL cextractvec(evec2,evecin_r,1,3,phase2,6,6)
          CALL rcmatvec3(rot(1,1,rnum),evecin_r,evecin_list(1,seg+1))
                              
        END do
!               
! Propagation complete. Now just need to convert amplitude to
! displacement, using the free-surface transfer matrix.
! MAY BE BUGGY FOR ANISOTROPIC TOPMOST LAYER.
!
! Check IF last segment really is upgoing. Otherwise, bail out.
!
        IF (p(3,nseg) .gt. 0.) THEN
          amp(1) = 0
          amp(2) = 0
          amp(3) = 0
          tt = 0
          RETURN          
        END IF
!
! Eigenvectors, again:
!
        IF (isoflag(1)) THEN
          CALL isotroc(aa(1,1,1,1,1),rho(1),p(1,nseg),
     &           p(2,nseg),evaltop,evectop)
        ELSE
          CALL anisotroc(aa(1,1,1,1,2),rho(1),p(1,nseg),
     &           p(2,nseg),evaltop,evectop)
        END IF
        phase1 = mod(phase(nseg,2)+2,6)+1
        CALL evec_check(evectop,evecin_list(1,nseg),phase1,mult,errflag)
        IF (errflag) THEN
          write (*,*) 'surface evals:',evaltop
        END IF
!
!  Get amplitude, from free-surface transfer, as follows:
!  Md=evec1(1:3,1:3),Mu=evec1(1:3,4:6),Nd=evec1(4:6,1:3),Nu=evec1(4:6,4:6)
!
        p1r(1) = 0
        p1r(2) = 0
        p1r(3) = 0
        p1r(mod(phase1-1,3)+1) = amp_list(nseg)*mult
        CALL cextractblock(evectop,Nu,4,6,4,6,6,6)
        CALL cextractblock(evectop,Nd,4,6,1,3,6,6)
        CALL cmatinv3(Nd,wrk3)
        CALL cmatmul3(wrk3,Nu,Nd)
        CALL cextractblock(evectop,Nu,1,3,1,3,6,6)
        CALL cmatmul3(Nu,Nd,wrk3)
        CALL cextractblock(evectop,Nu,1,3,4,6,6,6)
        CALL cmatlincomb3(cmplx(-1.),Nu,cmplx(1.),wrk3,MM)
        CALL crmatvec3(MM,p1r,amp)
!        
! Assemble travel-time:
!
        tt = 0
        do seg = 1, nseg
          tt = tt + tt_list(seg)
        END do
                        
      END SUBROUTINE calc_raysum 


!---------------------------------------------------------------------------
! Subroutine rot_evec
!
! Andrew Frederiksen
!
! Rotate degenerate post-rotation eigenvectors into a consistent
! coordinate system
!---------------------------------------------------------------------------
      SUBROUTINE rot_evec(evec,R)
      
      USE conf

      IMPLICIT NONE

      real R(3,3),theta
      complex evec(6,6),prod(3,6),evec2(6,2),x3(3),A(2,2)
      complex cdot3
      integer i,j
      data x3(1)/(0.,0.)/,x3(2)/(0.,0.)/,x3(3)/(1.,0.)/
        
        DO i = 1, 6
          CALL rcmatvec3(R,evec(1,i),prod(1,i))
        END DO
!
! Check if SH is off-horizontal in rotated frame.
! If (R*u3).x3 > 0 or (R*u6).x3 > 0...        
!
        IF ((abs(cdot3(prod(1,3),x3)) .gt. ztol) .or.
     &      (abs(cdot3(prod(1,6),x3)) .gt. ztol)) THEN
!                   
!         Downgoing set:
!
          theta = atan2(real(cdot3(prod(1,3),x3)),
     &                real(cdot3(prod(1,2),x3)))
          
          A(1,1) = cmplx(cos(theta))
          A(1,2) = cmplx(-sin(theta))
          A(2,1) = cmplx(sin(theta))
          A(2,2) = cmplx(cos(theta))
          CALL cmatmul(evec(1,2), A, evec2, 6, 2, 2)
          DO i = 1, 6
            DO j = 1, 2
              evec(i,j+1) = evec2(i,j)
            END DO
          END DO
!
!         Upgoing set:
!
          theta = atan2(real(cdot3(prod(1,6),x3)),
     &                real(cdot3(prod(1,5),x3)))
          A(1,1) = cmplx(cos(theta))
          A(1,2) = cmplx(-sin(theta))
          A(2,1) = cmplx(sin(theta))
          A(2,2) = cmplx(cos(theta))
          CALL cmatmul(evec(1,5),A,evec2,6,2,2)
          DO i = 1, 6
            DO j = 1, 2
              evec(i,j+4) = evec2(i,j)
            END DO
          END DO

        END if
        
      END SUBROUTINE rot_evec


!---------------------------------------------------------------------------
! Subroutine evec_check
!
! Andrew Frederiksen
!
! Consistency check for eigenvectors. Checks that evec (current
! eigenvectors) is consistent with invec (eigenvector for current
! phase calculated at previous step). phase1 (indicating which
! column in evec corresponds to invec) may be modified if a 90-degree
! rotation has occurred; if the sign is wrong, mult is set to -1.
!---------------------------------------------------------------------------

      SUBROUTINE evec_check(evec,invec,phase1,mult,errflag)
      
      IMPLICIT NONE

      complex evec(6,6),invec(3),echk(3,3),comp_list(3)
      real maxval
      integer phase1,mult,index,maxindex
      logical errflag
        
        mult = 1
        errflag = .false.
!
! index: 1 if downgoing, 4 if up
!
        index=int((sign(1.,real(phase1)-3.5)/2. + 0.5)*3.+1.)
!
! echk = evec(1:3,index:(index+2))  [relevant eigenvectors]
!
        CALL cextractblock(evec,echk,1,3,index,index+2,6,6)
!
! comp_list = real(invec'*echk) [check which evec is closest]
!
        CALL cvecmat3(invec,echk,comp_list)
        CALL cmax_vec_abs(comp_list,maxindex,maxval,3)
!
! check if the phase matches
!
        IF (maxindex.ne.(phase1-index+1)) THEN
          phase1 = maxindex+index-1
        END IF
!
! check if the sign matches
!
        IF (real(comp_list(maxindex)) .lt. 0) THEN
          mult = -1
        END IF
!
! check if the evecs don't line up. This shouldn't happen.
!
        IF (maxval .lt. 0.99) THEN
          WRITE (*,*) 'ERROR in evec_check, maxval: ',maxval,
     &                ' (should be 1)'
          WRITE (*,*) 'evec_check -- comp_list is',comp_list
          errflag = .true.
          WRITE (*,*) 'invec is',invec
        END IF
        
      END SUBROUTINE evec_check
      
      
!---------------------------------------------------------------------------
! Subroutine eig_inv
!
! Andrew Frederiksen
!
! Invert eigenvectors. Cribbed from Thomson's xeveci routine,
! dependent on characteristics of eigenvector matrix.
!---------------------------------------------------------------------------
      SUBROUTINE eiginv(eig,eig_i)
      
      IMPLICIT NONE

      complex eig(6,6),eig_i(6,6),wrk(6,6)
      integer i,j
        
        DO i = 1,3
          DO j = 1,3
            eig_i(i,j) = eig(j+3,i)
            eig_i(i,j+3) = eig(j,i)
            eig_i(i+3,j) = eig(j+3,i+3)
            eig_i(i+3,j+3) = eig(j,i+3)
          END DO
        END DO
        
        CALL cmatmul(eig_i,eig,wrk,6,6,6)
        
        DO i = 1, 6
          DO j = 1, 6
            eig_i(i,j) = eig_i(i,j)/wrk(i,i)
          END DO
        END DO
        
      END SUBROUTINE eig_inv

