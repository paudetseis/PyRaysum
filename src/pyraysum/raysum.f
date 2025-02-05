c####&
c New version with more efficient algorithm


      subroutine get_arrivals(travel_time,amplitude,thick,rho,isoflag,
     &           strike,dip,aa,ar_list,rot,baz,slow,sta_dx,sta_dy,
     &           phaselist,ntr,nseg,numph,nlay,amp_in)
c Get a list of arrivals, given model, geometry, and desired phases.
      
        implicit none
        include 'params.h'
        
c Interface variables:
        real travel_time(maxph,maxtr),amplitude(3,maxph,maxtr)
        real rho(maxlay),aa(3,3,3,3,maxlay),ar_list(3,3,3,3,maxlay,2)
        real strike(maxlay),dip(maxlay)
        real rot(3,3,maxlay),baz(maxtr),slow(maxtr),sta_dx(maxtr)
        real sta_dy(maxtr),thick(maxlay),amp_in
        integer phaselist(maxseg,2,maxph),ntr,nseg(maxph),numph,nlay
        logical isoflag(maxlay),bailout

c Function        
c        real rdot3
        
c Scratch variables:
        integer itr,iph,ilay,jlay,fchg
        real dzdx(maxlay),dzdy(maxlay),dthdx(maxlay),dthdy(maxlay)
        real thick_shift(maxlay)
c        real incvec(3),normvec(3)
        
c Calculate thickness changes caused by shifts from baseline
        do ilay=1,nlay          
          dzdx(ilay)=tan(dip(ilay))*sin(-strike(ilay))
          dzdy(ilay)=tan(dip(ilay))*sin(pi/2. -strike(ilay))
        end do
        do ilay=1,(nlay-1)
c          dthdx(ilay)=dzdx(ilay+1)
c          dthdy(ilay)=dzdy(ilay+1)
c          do jlay=1,ilay
c            dthdx(ilay)=dthdx(ilay)-dzdx(jlay)
c            dthdy(ilay)=dthdy(ilay)-dzdy(jlay)
c          end do
          dthdx(ilay)=dzdx(ilay+1)-dzdx(ilay)
          dthdy(ilay)=dzdy(ilay+1)-dzdy(ilay)
        end do
c        write(*,*) dthdy
        
c        write(*,'(X,5A11)')
c     &      'Trace','baz (deg)','slow (s/m)','dx (m)','dy (m)'
        do itr=1,ntr
c          write(*,'(X,I10,X,F10.3,1X,F10.8,1X,2F10.0)') 
c     &        itr,baz(itr)/pi*180,slow(itr),sta_dx(itr),sta_dy(itr)
          do ilay=1,nlay           
            thick_shift(ilay)=thick(ilay)+
     &            sta_dx(itr)*dthdx(ilay)+sta_dy(itr)*dthdy(ilay)
          end do

c          write (*,*) thick_shift

c          Check if incident phase is valid.
c          incvec(1)=-slow(itr)*cos(baz(itr))
c          incvec(2)=-slow(itr)*sin(baz(itr))
c          if (phaselist(1,1,1) .eq. 1) then
c            incvec(3)=-sqrt(1./aa(3,3,3,3,nlay)
c     &                     -incvec(1)**2-incvec(2)**2)
c          else
c            incvec(3)=-sqrt(1./aa(2,3,2,3,nlay)
c     &                     -incvec(1)**2-incvec(2)**2)
c          end if
c          normvec(1)=-sin(strike(nlay))*sin(dip(nlay))
c          normvec(2)=cos(strike(nlay))*sin(dip(nlay))
c          normvec(3)=-cos(dip(nlay))
c          write (*,*) incvec,normvec,rdot3(incvec,normvec)
c          if (rdot3(incvec,normvec) .le. 0) then
c            write (*,*) 'Incident ray arrives above bottom interface,',
c     &                  ' trace ',itr,' skipped.'
c          end if

          bailout = .true.
          do iph=1,numph
            if (bailout) then
              fchg=0
            else
              fchg=1
              if (nseg(iph) .eq. nseg(iph-1)) then                
                do while 
     &            ((phaselist(fchg,2,iph) .eq. phaselist(fchg,2,iph-1))
     &       .and. (phaselist(fchg,1,iph) .eq. phaselist(fchg,1,iph-1)))
                  fchg=fchg+1
                end do
              end if
            end if
c            write (*,*) 'iph=',iph,' fchg=',fchg
            call raysum(travel_time(iph,itr),amplitude(1,iph,itr),
     &                  thick_shift,rho,isoflag,aa,ar_list,rot,
     &                  baz(itr),slow(itr),phaselist(1,1,iph),nseg(iph),
     &                  nlay,amp_in,fchg,bailout)
          end do
        end do 
      
      end



c --------------------------------------------------------

      subroutine raysum(tt,amp,thick,rho,isoflag,aa,ar_list,rot,
     &                  baz,slow,phase,nseg,nlay,amp_in,fchg,bailout)
     
        implicit none
        include 'params.h'
c Main arrival-calculation routine. See raysum2.m
c Update -- keeps complex matrices as long as possible, for
c accuracy. p and amplitude are kept real -- complex values
c make no sense in ray theory.
        
c Interface variables:
        real tt,amp(3),thick(maxlay),rho(maxlay),aa(3,3,3,3,maxlay)
        real ar_list(3,3,3,3,maxlay,2),rot(3,3,maxlay),baz,slow
        real amp_in
        integer phase(maxseg,2),nseg,nlay,fchg
        logical isoflag(maxlay),bailout
 
c Scratch variables:
        integer seg,lay1,lay2,phase1,phase2,laytop,laybot,rnum,mult
        real p1r(3),p2r(3),Rt(3,3)
        complex evalbot(6),evecbot(6,6),evaltop(6)
        complex evectop(6,6),eval1(6),evec1(6,6),eval2(6)
        complex evec2(6,6),evecin_r(3)
        complex QQ(6,6),MM(3,3),Nu(3,3),Nd(3,3),wrk3(3,3),wrk6(6,6)
        logical upflag,rflag,fsflag,errflag
        
c Iteration-to-iteration storage variables 
        real amp_list(maxseg),p(3,maxseg),tt_list(maxseg)
        complex evecin_list(3,maxseg)
        save amp_list,p,tt_list,evecin_list
        
        bailout = .false.

        if (fchg .eq. 0) then      
c Initialize variables, get incident slowness vector
c Assume half-space is isotropic -- save some grief.
          tt_list(1)=0
          amp_list(1)=amp_in
c          write (*,*) 'amp_in=',amp_in
          p(1,1)=-slow*cos(baz)
          p(2,1)=-slow*sin(baz)
c          write (*,*) 'p1=',p1
          call isotroc(aa(1,1,1,1,nlay),rho(nlay),p(1,1),
     &                 p(2,1),evalbot,evecbot)
          phase1=mod(phase(1,2)+2,6)+1
          p(3,1)=real(evalbot(phase1))
        
c          print *,phase
        
c Store active eigenvector for later sign check.
          call cextractvec(evecbot,evecin_list(1,1),1,3,phase1,6,6)
          
        end if
       
        do seg=max((fchg-1),1),(nseg-1)
        
c          write (*,*) 'raysum: leg',seg
        
c 1 is incident, 2 is transmitted/refelcted 
          lay1=phase(seg,1)
          lay2=phase(seg+1,1)
c First and second phase, converted to eigenvector index
          phase1=mod(phase(seg,2)+2,6)+1
          phase2=mod(phase(seg+1,2)+2,6)+1
c Distinguish upper and lower layers.
          laytop=min(lay1,lay2)
          laybot=max(lay1,lay2)
c Upflag: true if the incident wave is upgoing.    
          upflag=phase1 .gt. 3
c Rflag: true for a reflection.
          rflag=lay1 .eq. lay2
          
          if (rflag) then
            if (upflag) then
              laytop=laybot-1
            else
              laybot=laytop+1
            end if
          end if

c Fsflag: free-surface reflection 
          fsflag=laytop .eq. 0
          if (fsflag .and. (.not. rflag)) then
            write (*,*) 'ERROR: can''t have free-surface transmission.'
          end if 
          
c          write (*,*) 'lay1,lay2,laytop,laybot,upflag,rflag,fsflag'
c          write (*,*) lay1,lay2,laytop,laybot,upflag,rflag,fsflag        
          

c laytop and laybot now contain correct layers for R/T calculation.
c Find correct rotator, determining active interface:
          if (upflag) then
            rnum=lay1
          else
            rnum=lay1+1
          end if
          call rtransp(rot(1,1,rnum),Rt,3,3)
          
c Project slowness onto interface
          call rmatvec3(Rt,p(1,seg),p1r)
c          write(*,*) 'p1r: ',p1r

c If the slowness vector is incorrectly oriented wrt the interface,
c we have a trapped phase.
          if ((upflag .and. (p1r(3) .gt. 0)) .or.
     &        ((.not. upflag) .and. (p1r(3) .lt. 0))) then
c            write (*,*) 'Ray does not intersect interface'
            amp(1)=0
            amp(2)=0
            amp(3)=0
            tt=0
            tt_list(seg+1)=0
            amp_list(seg+1)=0
            bailout = .true.
            return
          end if          
          
c Find evals/evecs. This is where it gets hairy...

c    Lower eigenvalues/vectors
          if (isoflag(laybot)) then
            call isotroc(aa(1,1,1,1,laybot),rho(laybot),
     &              p1r(1),p1r(2),evalbot,evecbot)
          else
c           Upper interface of laybot
            call anisotroc(ar_list(1,1,1,1,laybot,2),rho(laybot),
     &              p1r(1),p1r(2),evalbot,evecbot)
          end if
c          write (*,*) 'raysum: lower evals:',evalbot
          
c    Handle degeneracy (lower)
          if (abs((evalbot(2)-evalbot(3))/evalbot(2)) .le. ztol) then
c            write (*,*) 'Rotating evecbot',evalbot(2),evalbot(3)
            call rot_evec(evecbot,rot(1,1,rnum))
          end if
            
c     Upper eigenvalues/vectors. Don't bother if we're at the
c     free surface
c          print *,laytop,isoflag(laytop)
          if (.not.fsflag) then
            if (isoflag(laytop)) then
              call isotroc(aa(1,1,1,1,laytop),rho(laytop),
     &           p1r(1),p1r(2),evaltop,evectop)
            else
c           Lower interface of laytop
              call anisotroc(ar_list(1,1,1,1,laytop,1),rho(laytop),
     &              p1r(1),p1r(2),evaltop,evectop)
            end if
c            write (*,*) 'raysum: upper evals:',evaltop
c              Handle degeneracy (upper)             
            if ((evaltop(2) .ne. 0.) .and. 
     &         (abs((evaltop(2)-evaltop(3))/evaltop(2)) .le. ztol)) then
c              write (*,*) 'Rotating evaltop',evaltop(2),evaltop(3)
              call rot_evec(evectop,rot(1,1,rnum))
            end if
          end if                  
          
c Distribute evals/evecs between first and second propagation legs:

          if (rflag) then
c          Reflection
            if (upflag) then
c            Upgoing wave reflects off top layer -- incl. free-surf.
c            evec1=evecbot; eval1=evalbot; evec2=evec1; eval2=eval1;
              call cveccopy(evalbot,eval1,6)
              call cmatcopy(evecbot,evec1,6,6)
              call cveccopy(eval1,eval2,6)
              call cmatcopy(evec1,evec2,6,6)
            else
c            Downgoing wave reflects off bottom layer.
c            evec1=evectop; eval1=evaltop; evec2=evec1; eval2=eval1;
              call cveccopy(evaltop,eval1,6)
              call cmatcopy(evectop,evec1,6,6)
              call cveccopy(eval1,eval2,6)
              call cmatcopy(evec1,evec2,6,6)
            end if
          else
c          Transmission
            if (upflag) then
c            Upgoing wave enters top layer
c            evec1=evecbot; eval1=evalbot; evec2=evectop; eval2=evaltop;
              call cveccopy(evalbot,eval1,6)
              call cmatcopy(evecbot,evec1,6,6)
              call cveccopy(evaltop,eval2,6)
              call cmatcopy(evectop,evec2,6,6)
            else
c            Downgoing wave enters bottom layer
c            evec1=evectop; eval1=evaltop; evec2=evecbot; eval2=evalbot;
              call cveccopy(evaltop,eval1,6)
              call cmatcopy(evectop,evec1,6,6)
              call cveccopy(evalbot,eval2,6)
              call cmatcopy(evecbot,evec2,6,6)
            end if
          end if
          
c  If the eigenvalue is zero, we didn't get a transmission -- bail out
          if (abs(real(eval2(phase2))) .lt. ztol**2) then
c            write (*,*) 'Phase not transmitted'
            amp(1)=0
            amp(2)=0
            amp(3)=0
            tt=0
            tt_list(seg+1)=0
            amp_list(seg+1)=0
            bailout = .true.
            return
          end if

          
c Consistency check: are the incident phase and sign correct?
c Consistency may regained by altering phase1 or using a sign
c multiplier (mult)
          call rcmatvec3(Rt,evecin_list(1,seg),evecin_r)
          call evec_check(evec1,evecin_r,phase1,mult,errflag)
          if (errflag) then
            amp(1)=0
            amp(2)=0
            amp(3)=0
            tt=0
            tt_list(seg+1)=0
            amp_list(seg+1)=0
            bailout = .true.
            return
c           write (*,*) 'evaltop:',evaltop
c           write (*,*) 'evalbot:',evalbot
          end if
          
c Find new slowness vector. By Snell's law, interface-parallel
c components carry over. Eigenvalue provides remaining component.
          p2r(1)=p1r(1)
          p2r(2)=p1r(2)
          p2r(3)=real(eval2(phase2))
          call rmatvec3(rot(1,1,rnum),p2r,p(1,seg+1))
          
c Now we're ready to calculate the amplitude. Two cases to consider.
          if (fsflag) then          
c          Free-surface multiple
c            Partition N: Nu=evecbot(4:6,4:6), Nd=evecbot(4:6,1:3)
            call cextractblock(evecbot,Nu,4,6,4,6,6,6)
            call cextractblock(evecbot,Nd,4,6,1,3,6,6)
c Free-surface reflection matrix: MM=-Nd^(-1)*Nu
            call cmatinv3(Nd,wrk3)
            call cmatmul3(wrk3,Nu,MM)
            call cmatconst3(MM,(-1.,0.))           
          else
c          Layer interaction
c            scattering matrix QQ=evecbot^(-1)*evectop
            call eiginv(evecbot,wrk6)
            call cmatmul(wrk6,evectop,QQ,6,6,6)

c            Pull out reflection or transmission matrix MM
            if (rflag) then
              if (upflag) then
c                Ru = Q(1:3,4:6)*inv(Q(4:6,4:6))
                call cextractblock(QQ,Nd,1,3,4,6,6,6)
                call cextractblock(QQ,Nu,4,6,4,6,6,6)
                call cmatinv3(Nu,wrk3)
                call cmatmul3(Nd,wrk3,MM)
              else
c                Rd = -inv(Q(4:6,4:6))*Q(4:6,1:3)
                call cextractblock(QQ,Nu,4,6,4,6,6,6)
                call cextractblock(QQ,Nd,4,6,1,3,6,6)
                call cmatinv3(Nu,wrk3)
                call cmatmul3(wrk3,Nd,MM)
                call cmatconst3(MM,(-1.,0.))
              end if
            else
              if (upflag) then
c                Tu=inv(Q(4:6,4:6))
                call cextractblock(QQ,wrk3,4,6,4,6,6,6)
                call cmatinv3(wrk3,MM)
              else
c                Td=Q(1:3,1:3)-Q(1:3,4:6)*inv(Q(4:6,4:6))*Q(4:6,1:3)
                call cextractblock(QQ,Nu,4,6,4,6,6,6)
                call cmatinv3(Nu,wrk3)
                call cextractblock(QQ,Nu,4,6,1,3,6,6)
                call cmatmul3(wrk3,Nu,Nd)
                call cextractblock(QQ,Nu,1,3,4,6,6,6)
                call cmatmul3(Nu,Nd,wrk3)
                call cextractblock(QQ,Nu,1,3,1,3,6,6)
                call cmatlincomb3(cmplx(1.),Nu,cmplx(-1.),wrk3,MM)
              end if
            end if            
          end if
          
c       Carry amplitude across interface
          amp_list(seg+1)=amp_list(seg)*real(mult)*
     &       real(MM(mod(phase2-1,3)+1,mod(phase1-1,3)+1))
     
c          if (abs(amp_list(seg+1)/amp_in) .lt. ztol) then
cc            write (*,*) 'No significant amplitude'
c            amp(1)=0
c            amp(2)=0
c            amp(3)=0
c            tt=0
c            tt_list(seg+1)=abs(p(3,seg+1))*thick(lay2)
c            amp_list(seg+1)=0
c            bailout = .true.
c            return
c          end if
     
c       Add post-interface segment to travel time
          if (thick(lay2) .lt. 0.) then
            write (*,*) 'ERROR: negative thickness!'
          end if
          tt_list(seg+1)=abs(p(3,seg+1))*thick(lay2)

c       Pull out relevant eigenvector for next-segment consistency
c       check. evecin=evec_out(:,seg+1)=R*evec2(1:3,phase2)
          call cextractvec(evec2,evecin_r,1,3,phase2,6,6)
          call rcmatvec3(rot(1,1,rnum),evecin_r,evecin_list(1,seg+1))
                              
        end do
                
c      Propagation complete. Now just need to convert amplitude to
c      displacement, using the free-surface transfer matrix.
c      MAY BE BUGGY FOR ANISOTROPIC TOPMOST LAYER.

c         Check if last segment really is upgoing. Otherwise, bail out.
        if (p(3,nseg) .gt. 0.) then
c          write (*,*) 'Ray does not reach surface'
          amp(1)=0
          amp(2)=0
          amp(3)=0
          tt=0
          return          
        end if

c    Eigenvectors, again:
        if (isoflag(1)) then
          call isotroc(aa(1,1,1,1,1),rho(1),p(1,nseg),
     &           p(2,nseg),evaltop,evectop)
        else
          call anisotroc(aa(1,1,1,1,1),rho(1),p(1,nseg),
     &           p(2,nseg),evaltop,evectop)
        end if
        phase1=mod(phase(nseg,2)+2,6)+1
        call evec_check(evectop,evecin_list(1,nseg),phase1,mult,errflag)
        if (errflag) then
          amp(1)=0
          amp(2)=0
          amp(3)=0
          tt=0
          tt_list(seg+1)=0
          amp_list(seg+1)=0
          bailout = .true.
          return
c         write (*,*) 'surface evals:',evaltop
        end if

c  Get amplitude, from free-surface transfer, as follows:
c  Md=evec1(1:3,1:3),Mu=evec1(1:3,4:6),Nd=evec1(4:6,1:3),Nu=evec1(4:6,4:6)
c      cu(mod(phase1-1,3)+1)=amp_out(nseg)*mult;
c      amp=-(Mu-Md*inv(Nd)*Nu)*cu; p1r stands for cu (wavevector)
        p1r(1)=0
        p1r(2)=0
        p1r(3)=0
        p1r(mod(phase1-1,3)+1)=amp_list(nseg)*mult
        call cextractblock(evectop,Nu,4,6,4,6,6,6)
        call cextractblock(evectop,Nd,4,6,1,3,6,6)
        call cmatinv3(Nd,wrk3)
        call cmatmul3(wrk3,Nu,Nd)
        call cextractblock(evectop,Nu,1,3,1,3,6,6)
        call cmatmul3(Nu,Nd,wrk3)
        call cextractblock(evectop,Nu,1,3,4,6,6,6)
        call cmatlincomb3(cmplx(-1.),Nu,cmplx(1.),wrk3,MM)
        call crmatvec3(MM,p1r,amp)
        
c Assemble travel-time:
        tt=0
        do seg=1,nseg
          tt=tt+tt_list(seg)
        end do
        
c        write (*,'(A10,X,$)') 'amp_list='
c        do seg=1,nseg
c          write (*,'(X,F10.7,$)') amp_list(seg)
c        end do
c        write (*,*)
c        write (*,*) 'p_list:'
c        do seg=1,nseg
c          write (*,'(3(X,F9.7))') p(1,seg),p(2,seg),p(3,seg)
c        end do
                
      end
      
      
c ---------------------------------------

      subroutine rot_evec(evec,R)
c       Rotate degenerate post-rotation eigenvectors into a consistent
c       coordinate system
      
        implicit none
        include 'params.h'
        real R(3,3),theta
        complex evec(6,6),prod(3,6),evec2(6,2),x3(3),A(2,2)
        complex cdot3
        integer i,j
        data x3(1)/(0.,0.)/,x3(2)/(0.,0.)/,x3(3)/(1.,0.)/
        
        do i=1,6
c             prod_i = R*evec(1:3,i)
          call rcmatvec3(R,evec(1,i),prod(1,i))
        end do

c       Check if SH is off-horizontal in rotated frame.
c       If (R*u3).x3 > 0 or (R*u6).x3 > 0...        
        if ((abs(cdot3(prod(1,3),x3)) .gt. ztol) .or.
     &      (abs(cdot3(prod(1,6),x3)) .gt. ztol)) then
                    
c         Downgoing set:
          theta=atan2(real(cdot3(prod(1,3),x3)),
     &                real(cdot3(prod(1,2),x3)))
          
          A(1,1)=cmplx(cos(theta))
          A(1,2)=cmplx(-sin(theta))
          A(2,1)=cmplx(sin(theta))
          A(2,2)=cmplx(cos(theta))
c          evec(:,2:3)*A
          call cmatmul(evec(1,2),A,evec2,6,2,2)
          do i=1,6
            do j=1,2
              evec(i,j+1)=evec2(i,j)
            end do
          end do
          
c         Upgoing set:
          theta=atan2(real(cdot3(prod(1,6),x3)),
     &                real(cdot3(prod(1,5),x3)))
          A(1,1)=cmplx(cos(theta))
          A(1,2)=cmplx(-sin(theta))
          A(2,1)=cmplx(sin(theta))
          A(2,2)=cmplx(cos(theta))
c          evec(:,5:6)*A
          call cmatmul(evec(1,5),A,evec2,6,2,2)
          do i=1,6
            do j=1,2
              evec(i,j+4)=evec2(i,j)
            end do
          end do
          
        
        end if
        
      end


c -----------------------------------------

      subroutine evec_check(evec,invec,phase1,mult,errflag)
c Consistency check for eigenvectors. Checks that evec (current
c eigenvectors) is consistent with invec (eigenvector for current
c phase calculated at previous step). phase1 (indicating which
c column in evec corresponds to invec) may be modified if a 90-degree
c rotation has occurred; if the sign is wrong, mult is set to -1.
      
        implicit none
        complex evec(6,6),invec(3),echk(3,3),comp_list(3)
        real maxval
        integer phase1,mult,index,maxindex
        logical errflag
        
        mult=1
        errflag=.false.
        
c       index: 1 if downgoing, 4 if up
        index=int((sign(1.,real(phase1)-3.5)/2. + 0.5)*3.+1.)
c          echk=evec(1:3,index:(index+2))  [relevant eigenvectors]
        call cextractblock(evec,echk,1,3,index,index+2,6,6)
c          comp_list=real(invec'*echk) [check which evec is closest]
        call cvecmat3(invec,echk,comp_list)
        call cmax_vec_abs(comp_list,maxindex,maxval,3)
c       check if the phase matches
        if (maxindex.ne.(phase1-index+1)) then
c          write (*,*) 'evec_check: Phase mismatch ',maxindex,
c     &                (phase1-index+1)
          phase1=maxindex+index-1
        end if
c       check if the sign matches
        if (real(comp_list(maxindex)) .lt. 0) then
c          write (*,*) 'evec_check: Sign mismatch'
          mult=-1
        end if
c       check if the evecs don't line up. This shouldn't happen.
        if (maxval .lt. 0.99) then
          write (*,*) 'WARNING in evec_check, maxval: ',maxval,
     &                ' (should be 1). Ignoring phase.'
c         write (*,*) 'evec_check -- comp_list is',comp_list
          errflag=.true.
c         write (*,*) 'invec is',invec
        end if
        
      end
      
      
c ------------------------------------------

      subroutine eiginv(eig,eig_i)
c      Invert eigenvectors. Cribbed from Thomson's xeveci routine,
c      dependant on characteristics of eigenvector matrix.
      
        implicit none
        complex eig(6,6),eig_i(6,6),wrk(6,6)
        integer i,j
        
        do i=1,3
          do j=1,3
            eig_i(i,j)=eig(j+3,i)
            eig_i(i,j+3)=eig(j,i)
            eig_i(i+3,j)=eig(j+3,i+3)
            eig_i(i+3,j+3)=eig(j,i+3)
          end do
        end do
        
        call cmatmul(eig_i,eig,wrk,6,6,6)
        
        do i=1,6
          do j=1,6
            eig_i(i,j)=eig_i(i,j)/wrk(i,i)
          end do
        end do
        
      end
           
      
      

        
        
        
