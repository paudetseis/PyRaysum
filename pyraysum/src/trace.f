c####&

c Build traces from arrivals. tt and amp are the travel-time and
c amplitude lists, ntr is the number of traces, nph is the number
c of phases, dt is the sample spacing (in seconds), align is the
c index of the phase to align on (0 for no alignment), shift is the
c trace time for tt=0 (seconds), and Tr_cart is the resulting
c traces in Cartesian coodinates, Tr_cart(:,1) being N-S, 2 E-W, and
c 3 vertical.
      subroutine make_traces(tt,amp,ntr,nph,nsamp,dt,align,shift,
     &                       verbose,Tr_cart)
        
        implicit none
        include 'params.h'
        
c        Interface variables
        real tt(maxph,maxtr),amp(3,maxph,maxtr),dt,shift
        real Tr_cart(3,maxsamp,maxtr)
        integer ntr,nph,align,nsamp
        
c        Scratch variables
        integer itr,iph,isamp,icomp
        real delta,max_t,tt_s,curamp(3)
        logical verbose
c        Zero traces
        do isamp=1,nsamp
          do itr=1,ntr
            Tr_cart(1,isamp,itr)=0.
            Tr_cart(2,isamp,itr)=0.
            Tr_cart(3,isamp,itr)=0.
          end do
        end do
        
c        Max. time that fits. The first sample is at t=0.
        max_t=dt*real(nsamp-1)
        
        do itr=1,ntr
          if (align .eq. 0) then
            delta=-shift
          else
            delta=tt(align,itr)-shift
          end if
          do iph=1,nph
            tt_s=tt(iph,itr)-delta
c           shift output arrivals as well
            tt(iph,itr)=tt_s
            if ((tt_s .lt. 0.) .or. (tt_s .gt. max_t)) then
              if (verbose) then
                write (*,*) 'Phase ',iph,' cropped out of trace ',itr
              end if
            else
              do icomp=1,3
                curamp(icomp)=amp(icomp,iph,itr)
                if (curamp(icomp) .ne. curamp(icomp)) then
                  write (*,*) '!!! ERROR -- amplitude is not a number'
                  write (*,*) 'Trace ',itr,' phase ',iph,' comp ',icomp
                  write (*,*) 'Setting to zero, continuing.'
                  curamp(icomp)=0
                end if
              end do
              call putspike(Tr_cart(1,1,itr),tt_s/dt,curamp,nsamp) 
            end if
          end do
        end do
        
      end
      

c -------------------------

c Insert a spike (delta fct) into a trace      
      subroutine putspike(traces,mu,amp,nsamp)
      
        implicit none
        include 'params.h'
        
        real traces(3,maxsamp),mu,amp(3),gauss
        integer nsamp
        
        integer isamp

c        write (*,*) 'putgauss:',mu,amp
				isamp = int(mu)
        if ((isamp .gt. 0) .and. (isamp .le. nsamp)) then
          traces(1,isamp)=traces(1,isamp)+1.*amp(1)
          traces(2,isamp)=traces(2,isamp)+1.*amp(2)
          traces(3,isamp)=traces(3,isamp)+1.*amp(3)
        end if
        
      end
      
c ---------------------------

c  Copy traces from Tr_cart to Tr_ph
      subroutine copy_traces(Tr_cart,ntr,nsamp,Tr_ph)
      
        implicit none 
        include 'params.h'
        
        real Tr_cart(3,maxsamp,maxtr),Tr_ph(3,maxsamp,maxtr)
        integer ntr,itr,isamp,nsamp
        
        do itr=1,ntr
          do isamp=1,nsamp
            Tr_ph(1,isamp,itr)=Tr_cart(1,isamp,itr)
            Tr_ph(2,isamp,itr)=Tr_cart(2,isamp,itr)
            Tr_ph(3,isamp,itr)=Tr_cart(3,isamp,itr)
          end do
        end do
        
      end


c ---------------------------

c  Copy amplitudes from amp_cart to amp_ph
      subroutine copy_amplitudes(amp_cart,ntr,npha,amp_ph)
      
        implicit none 
        include 'params.h'
        
        real amp_cart(3,maxph,maxtr),amp_ph(3,maxph,maxtr)
        integer ntr,itr,ipha,npha
        
        do itr=1,ntr
          do ipha=1,npha
            amp_ph(1,ipha,itr)=amp_cart(1,ipha,itr)
            amp_ph(2,ipha,itr)=amp_cart(2,ipha,itr)
            amp_ph(3,ipha,itr)=amp_cart(3,ipha,itr)
          end do
        end do
        
      end


c ---------------------------

c  Rotate traces into R-T-Z system. Tr_cart and Tr_ph can be the same
c  variable.
      subroutine rot_traces(Tr_cart,baz,ntr,nsamp,Tr_ph)
      
        implicit none 
        include 'params.h'
        
        real Tr_cart(3,maxsamp,maxtr),baz(maxtr),Tr_ph(3,maxsamp,maxtr)
        real ct,st,dr,dt,dz
        integer ntr,itr,isamp,nsamp
        
        do itr=1,ntr
          ct=cos(baz(itr))
          st=sin(baz(itr))
          do isamp=1,nsamp
            dr=ct*Tr_cart(1,isamp,itr)+st*Tr_cart(2,isamp,itr)
            dt=-st*Tr_cart(1,isamp,itr)+ct*Tr_cart(2,isamp,itr)
            dz=Tr_cart(3,isamp,itr)
            Tr_ph(1,isamp,itr)=dr
            Tr_ph(2,isamp,itr)=dt
            Tr_ph(3,isamp,itr)=dz
          end do
        end do
        
      end

c ---------------------------

c  Rotate amplitudes into R-T-Z system.
      subroutine rot_amplitudes(amp_cart,baz,ntr,npha,amp_ph)
      
        implicit none 
        include 'params.h'
        
        real amp_cart(3,maxph,maxtr),baz(maxtr),amp_ph(3,maxph,maxtr)
        real ct,st,dr,dt,dz
        integer ntr,itr,ipha,npha
        
        do itr=1,ntr
          ct=cos(baz(itr))
          st=sin(baz(itr))
          do ipha=1,npha
            dr=ct*amp_cart(1,ipha,itr)+st*amp_cart(2,ipha,itr)
            dt=-st*amp_cart(1,ipha,itr)+ct*amp_cart(2,ipha,itr)
            dz=amp_cart(3,ipha,itr)
            amp_ph(1,ipha,itr)=dr
            amp_ph(2,ipha,itr)=dt
            amp_ph(3,ipha,itr)=dz
          end do
        end do
        
      end



c ---------------------------

c Rotate traces from NS-EW-Z coordinates to P-SV-SH system.
c Tr_cart and Tr_ph can be the same variable.
      subroutine fs_traces(Tr_cart,baz,slow,alpha,beta,rho,ntr,
     &                     nsamp,Tr_ph)
      
        implicit none
        include 'params.h'
        
        real Tr_cart(3,maxsamp,maxtr),baz(maxtr),slow(maxtr),alpha,beta
        real rho,Tr_ph(3,maxsamp,maxtr)
        integer ntr,nsamp
        
        integer itr,isamp,j
        real p1,p2,a(3,3,3,3)
        complex eval(6),evec(6,6),Md(3,3),Mu(3,3),Nd(3,3),Nu(3,3)
        complex wrk(3,3),invop(3,3),u(3),v(3)
        
        a(3,3,3,3)=alpha**2
        a(2,3,2,3)=beta**2
                
        do itr = 1,ntr
        
          p1=-slow(itr)*cos(baz(itr))
          p2=-slow(itr)*sin(baz(itr))
          call isotroc(a,rho,p1,p2,eval,evec)
          call cextractblock(evec,Md,1,3,1,3,6,6)
          call cextractblock(evec,Mu,1,3,4,6,6,6)
          call cextractblock(evec,Nd,4,6,1,3,6,6)
          call cextractblock(evec,Nu,4,6,4,6,6,6)
          call cmatinv3(Nd,invop)        
          call cmatmul3(invop,Nu,wrk)
          call cmatmul3(Md,wrk,invop)
          call cmatlincomb3((1.,0),Mu,(-1.,0),invop,wrk)
          call cmatinv3(wrk,invop)
          
          do isamp=1,nsamp
            do j=1,3
              u(j)=Tr_cart(j,isamp,itr)
            end do
            call cmatvec3(invop,u,v)
            Tr_ph(1,isamp,itr)=-real(v(1))
            Tr_ph(2,isamp,itr)=-real(v(2))
            Tr_ph(3,isamp,itr)=-real(v(3))
          end do
          
        end do
        
      end
      
c ---------------------------

c Rotate amplitudes from NS-EW-Z coordinates to P-SV-SH system.
      subroutine fs_amplitudes(amp_cart,baz,slow,alpha,beta,rho,ntr,
     &                     npha,amp_ph)
      
        implicit none
        include 'params.h'
        
        real amp_cart(3,maxph,maxtr),baz(maxtr),slow(maxtr),alpha,beta
        real rho,amp_ph(3,maxph,maxtr)
        integer ntr,npha
        
        integer itr,ipha,j
        real p1,p2,a(3,3,3,3)
        complex eval(6),evec(6,6),Md(3,3),Mu(3,3),Nd(3,3),Nu(3,3)
        complex wrk(3,3),invop(3,3),u(3),v(3)
        
        a(3,3,3,3)=alpha**2
        a(2,3,2,3)=beta**2
                
        do itr = 1,ntr
        
          p1=-slow(itr)*cos(baz(itr))
          p2=-slow(itr)*sin(baz(itr))
          call isotroc(a,rho,p1,p2,eval,evec)
          call cextractblock(evec,Md,1,3,1,3,6,6)
          call cextractblock(evec,Mu,1,3,4,6,6,6)
          call cextractblock(evec,Nd,4,6,1,3,6,6)
          call cextractblock(evec,Nu,4,6,4,6,6,6)
          call cmatinv3(Nd,invop)        
          call cmatmul3(invop,Nu,wrk)
          call cmatmul3(Md,wrk,invop)
          call cmatlincomb3((1.,0),Mu,(-1.,0),invop,wrk)
          call cmatinv3(wrk,invop)
          
          do ipha=1,npha
            do j=1,3
              u(j)=amp_cart(j,ipha,itr)
            end do
            call cmatvec3(invop,u,v)
            amp_ph(1,ipha,itr)=-real(v(1))
            amp_ph(2,ipha,itr)=-real(v(2))
            amp_ph(3,ipha,itr)=-real(v(3))
          end do
          
        end do
        
      end
      
      
      
c ----------------------------------------------------



      subroutine norm_arrivals(amp,baz,slow,alpha,beta,rho,ntr,numph,
     &                         arr,comp)
c      Normalize arrival amplitudes amp by the amplitude of arrival arr,
c      component (in P-SV-SH system) comp. The arrivals are assumed to
c      be in NS-EW-Z coordinates.

        implicit none
        include 'params.h'
        
        real amp(3,maxph,maxtr),baz(maxtr),slow(maxtr),alpha,beta,rho
        integer ntr,numph,arr,comp
        
        integer itr,iph,j
        real normamp,p1,p2,a(3,3,3,3)
        complex eval(6),evec(6,6),Md(3,3),Mu(3,3),Nd(3,3),Nu(3,3)
        complex wrk(3,3),invop(3,3),u(3),v(3)
        
        a(3,3,3,3)=alpha**2
        a(2,3,2,3)=beta**2
        
        do itr=1,ntr
        
          p1=-slow(itr)*cos(baz(itr))
          p2=-slow(itr)*sin(baz(itr))
          call isotroc(a,rho,p1,p2,eval,evec)
          call cextractblock(evec,Md,1,3,1,3,6,6)
          call cextractblock(evec,Mu,1,3,4,6,6,6)
          call cextractblock(evec,Nd,4,6,1,3,6,6)
          call cextractblock(evec,Nu,4,6,4,6,6,6)
          call cmatinv3(Nd,invop)        
          call cmatmul3(invop,Nu,wrk)
          call cmatmul3(Md,wrk,invop)
          call cmatlincomb3((1.,0),Mu,(-1.,0),invop,wrk)
          call cmatinv3(wrk,invop)
          
          do j=1,3
            u(j)=amp(j,arr,itr)
          end do
          call cmatvec3(invop,u,v)
          normamp=-real(v(comp))
          
          do iph=1,numph
            do j=1,3
              if (normamp .gt. 0.) then
                amp(j,iph,itr)=amp(j,iph,itr)/normamp
              else
                amp(j,iph,itr)=0.
              end if
            end do
          end do
          
        end do
        
      end
