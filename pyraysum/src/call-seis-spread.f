c Calculates a spread of seismograms for a given model.
c soubroutine provides python itnerface

c####&

      subroutine call_seis_spread(
     &                            thick, rho, alpha, beta, isoflag, pct,
     &                            trend, plunge, strike, dip, nlay,
     &                            baz, slow, sta_dx, sta_dy, ntr,
     &                            iphname, mults, nsamp, dt, align,
     &                            shift, out_rot, verb,
     &                            Tr_ph, Tr_cart)

c Anyone who fails to start a Fortran program with this line
c should be severely beaten:
        implicit none
    
c Include constants:
        include 'params.h'

c ========================================
c Input parameters !! intent(in) them !!

c Incoming wave type
        character iphname*(namelen)
Cf2py intent(in) :: iphname

c Multiple flag, number of samples, trace alignment
c coordinate system of output
        integer mults, nsamp, align
        integer out_rot
Cf2py   integer intent(in) :: mults, nsamp, align
Cf2py   integer intent(in) :: out_rot

c Time step, trace shift 
        real dt, shift
Cf2py   intent(in) :: dt, shift

c model parameters
        integer nlay
        real thick(maxlay), rho(maxlay), alpha(maxlay), beta(maxlay)
        real pct(maxlay), trend(maxlay), plunge(maxlay)
        real strike(maxlay), dip(maxlay)
        logical isoflag(maxlay)
Cf2py   integer intent(in) :: nlay
Cf2py   intent(in) :: thick, rho, alpha, beta
Cf2py   intent(in) :: pct, trend, plunge
Cf2py   intent(in) :: strike, dip
Cf2py   intent(in) :: isoflag

c Geometry parameters
        real baz(maxtr), slow(maxtr), sta_dx(maxtr), sta_dy(maxtr)
        integer ntr
Cf2py   integer intent(in) :: ntr
Cf2py   intent(in) :: baz, slow, sta_dx, sta_dy
        
c ========================================
c Output parameters !! intent(out) them !!
c Traces

        real Tr_cart(3,maxsamp,maxtr)
        real Tr_ph(3,maxsamp,maxtr)
Cf2py intent(out) :: Tr_cart, Tr_ph

c ==================
c Internal variables
c Scratch variables:
        integer i, j, iargc, verb
        real amp_in
        logical verbose

c Phase parameters
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph,iphase

c Arrivals
        real travel_time(maxph,maxtr),amplitude(3,maxph,maxtr)
        
c   aa is a list of rank-4 tensors (a_ijkl = c_ijkl/rho)
c   rot is a list of rotator matrices, used to rotate into the local
c   coordinate system of each interface.
        real aa(3,3,3,3,maxlay), rot(3,3,maxlay)
c   ar_list is aa, pre-rotated into reference frames for the interfaces
c   above (last index=2) and below (1) each respective layer.
        real ar_list(3,3,3,3,maxlay,2)
        
        
c Read in parameters from file 'raysum-params'

        verbose=.false.
        if (verb .eq. 1) then
          verbose=.true.
        end if

c Determine initial phase
        if (iphname .eq. 'P') then
          iphase=1
        else if (iphname .eq. 'SV') then
          iphase=2
        else if (iphname .eq. 'SH') then
          iphase=3
        else
          write (*,*) iphname,' is not a valid phase type.'
          write (*,*) 'Valid phases are P, SV and SH.'
          stop
        end if
        if (verbose) then
          write (*,*) 'Initial phase is ',iphname
        end if

c Write out model for testing      
        if (verbose) then
          call writemodel(6,thick,rho,alpha,beta,isoflag,
     &                  pct,trend,plunge,strike,dip,nlay)
        end if
          
c Set up model for calculation, make rotators
        call buildmodel(aa,ar_list,rot,thick,rho,alpha,beta,isoflag,
     &                  pct,trend,plunge,strike,dip,nlay)
     
c Return geometry for testing 
        if (verbose) then
          call writegeom(6,baz,slow,sta_dx,sta_dy,ntr)
        end if
        
c Generate phase list
        numph=0
        if (mults .ne. 3) then
          call ph_direct(phaselist,nseg,numph,nlay,iphase)
        end if
        if (mults .eq. 1) then
          call ph_fsmults(phaselist,nseg,numph,nlay,1,iphase)
        end if
        if (mults .eq. 2) then
          do j=1,nlay-1
            call ph_fsmults(phaselist,nseg,numph,nlay,j,iphase)
          end do
        end if
        if (verbose) then
          call printphases(phaselist,nseg,numph)
        end if
        
c Perform calculation                   
        amp_in=1.
        call get_arrivals(travel_time,amplitude,thick,rho,isoflag,
     &       strike,dip,aa,ar_list,rot,baz,slow,sta_dx,sta_dy,
     &       phaselist,ntr,nseg,numph,nlay,amp_in)
     
c Normalize arrivals
        if (iphase .eq. 1) then 
          call norm_arrivals(amplitude,baz,slow,alpha(1),beta(1),
     &                       rho(1),ntr,numph,1,1)
        end if
                 
c Assemble traces
        call make_traces(travel_time,amplitude,ntr,numph,nsamp,
     &                   dt,align,shift,Tr_cart)
     
        if (out_rot .eq. 0) then
c          call writetraces(iounit2,Tr_cart,ntr,nsamp,dt,align,shift)
        else
          if (out_rot .eq. 1) then
c Rotate to RTZ
            call rot_traces(Tr_cart,baz,ntr,nsamp,Tr_ph)
          else
c Rotate to wavevector coordinates
            call fs_traces(Tr_cart,baz,slow,alpha(1),beta(1),
     &                     rho(1),ntr,nsamp,Tr_ph)
          end if
c Write results
c          call writetraces(iounit2,Tr_ph,ntr,nsamp,dt,align,shift)
        end if
                
      end