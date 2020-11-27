c####&


      subroutine readmodel(filename,thick,rho,alpha,beta,isoflag,
     &                     pct_a,pct_b,trend,plunge,strike,dip,nlay)
     
c Read in a model file from disk.
          
        implicit none
        include 'params.h'
        
        character filename*(namelen), buffer*(buffsize)
        integer eof,nlay,i,isoint
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct_a(maxlay),pct_b(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay)
        logical isoflag(maxlay)
        
        open(unit=iounit1,file=filename,status='old')
        
c Read in lines from file, skipping lines that start with '#' (comments)
	eof=0
        nlay=0
        call getline(iounit1,buffer,eof)
        do while ((eof .eq. 0) .and. (nlay .lt. maxlay))
          nlay=nlay+1
c          write(*,*) buffer
          read(buffer,*) thick(nlay),rho(nlay),alpha(nlay),beta(nlay),
     &                   isoint,pct_a(nlay),pct_b(nlay),
     &                   trend(nlay),plunge(nlay),strike(nlay),dip(nlay)
          if (isoint .eq. 0) then
            isoflag(nlay) = .false.
          else
            isoflag(nlay) = .true.
          end if
          call getline(1,buffer,eof)
        end do 
c        write(*,*) 'Got ',nlay,' layers.'   
               
        close(iounit1)
        
c Convert radians to degrees
        do i=1,nlay
          strike(i)=strike(i)/180. * pi
          dip(i)=dip(i)/180. * pi
          trend(i)=trend(i)/180. * pi
          plunge(i)=plunge(i)/180. * pi
        end do
        
        
      end
      
      
c ------------------------------------
      
c Read geometry file
      subroutine readgeom(filename,baz,slow,sta_dx,sta_dy,ntr)
      
        implicit none
        include 'params.h'
        
        character filename*(namelen),buffer*(buffsize)
        real baz(maxtr),slow(maxtr),sta_dx(maxtr),sta_dy(maxtr)
        integer eof,ntr        

        open(unit=iounit1,file=filename,status='old')
        
c       Read in lines from file, skipping comments
	eof=0
        ntr=0
        call getline(iounit1,buffer,eof)
        do while ((eof .eq. 0) .and. (ntr .lt. maxtr))
          ntr=ntr+1
          read(buffer,*) baz(ntr),slow(ntr),sta_dx(ntr),sta_dy(ntr)
          call getline(1,buffer,eof)
c         Convert to radians
          baz(ntr)=baz(ntr)/180. * pi
        end do 
c        write(*,*) 'Got ',ntr,' trace geometries.'   
        close(iounit1)
        
      end

c ----------------------------------

      subroutine readtraces(filename,traces,ntr,nsamp,dt,align,shift)
      
        implicit none
        include 'params.h'
        
        integer ntr,align,eof,itr,isamp,nsamp
        real traces(3,maxsamp,maxtr),dt,shift
        character filename*(namelen),buffer*(buffsize)
        
        open(unit=iounit1,file=filename,status='old')
        
        eof=0
        call getline(iounit1,buffer,eof)
        read(buffer,*) ntr,nsamp,dt,align,shift
        if (nsamp .gt. maxsamp) then
          write (*,*) 'ERROR -- trace too long'
          stop
        end if
        
        do itr=1,ntr
          do isamp=1,nsamp
            call getline(iounit1,buffer,eof)
            read (buffer,*) traces(1,isamp,itr),traces(2,isamp,itr),
     &                      traces(3,isamp,itr)
          end do
        end do
        
      end
        

c ----------------------------------

c Get the next line from a file that isn't a comment
      subroutine getline(unit,buffer,eof)
      
        implicit none
        include 'params.h'
      
        character buffer*(buffsize),firstnonblank
        integer unit, eof
                
        read(unit,'(A)',iostat=eof) buffer
        do while (((firstnonblank(buffer) .eq. '#') .or.
     &           (firstnonblank(buffer) .eq. ' ')) .and. (eof .eq. 0))
          read(unit,'(A)',iostat=eof) buffer
        end do
                    
      end
      
      
c --------------------------------

c Obtain the first non-blank character in a string.
      character function firstnonblank(buffer)

        implicit none
        include 'params.h'        
        character buffer*(buffsize)
        integer i
        
        i=1
        do while ((buffer(i:i) .eq. ' ') .and. (i .lt. buffsize))
          i=i+1
        end do
        firstnonblank=buffer(i:i)
        
      end

c -----------------------------------
      
c Print out a model on the given unit:
      subroutine writemodel(unit,thick,rho,alpha,beta,isoflag,
     &                      pct_a,pct_b,trend,plunge,strike,dip,nlay)
     
        implicit none
        include 'params.h'
        
        integer unit, nlay, i, isonum
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct_a(maxlay),pct_b(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay)
        logical isoflag(maxlay)
        
        write(unit,'(A1,4A7,A4,6A7)') '#','thick','rho','alpha','beta',
     &        'iso','%P','%S','trend','plunge','strike','dip'
        
        do i=1,nlay
          if (isoflag(i)) then
            isonum=1
          else
            isonum=0
          end if
          write(unit,'(1X,4F7.0,I4,2F7.1,4F7.0)')
     &          thick(i),rho(i),alpha(i),beta(i),isonum,
     &          pct_a(i),pct_b(i),trend(i)*180/pi,plunge(i)*180/pi,
     &          strike(i)*180./pi,dip(i)*180./pi
        end do
             
      end
      
c -----------------------------------
      
      subroutine writegeom(unit,baz,slow,sta_dx,sta_dy,ntr)
      
        implicit none
        include 'params.h'
        
        real baz(maxtr),slow(maxtr),sta_dx(maxtr),sta_dy(maxtr)
        integer unit,ntr,i
        
c        write(unit,'(A1,2A7,3A10)') '#','trace','baz','slow','dx','dy'
        write(unit,'(A1,A7,3A10)') '#','baz','slow','dx','dy'
        
        do i=1,ntr
c          write(unit,'(1X,1I4,3X,F7.2,1X,F9.8,1X,2F10.0)')
c     &       i,baz(i)*180./pi,slow(i),sta_dx(i),sta_dy(i)        
          write(unit,'(1X,F7.2,1X,F9.8,1X,2F10.0)')
     &       baz(i)*180./pi,slow(i),sta_dx(i),sta_dy(i)        
        end do
        
      end
      
c ----------------------------------

      subroutine printphases(phaselist,nseg,numph)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer iseg,iph
        
        do iph=1,numph
          do iseg=1,nseg(iph)
            write(*,'(1X,I2,1X,I1)') phaselist(iseg,1,iph),
     &                                  phaselist(iseg,2,iph)
          end do
          write(*,*)
        end do
       end

c --------------------------------
      subroutine readphases(filename,phaselist,nseg,numph)
      
        implicit none
        include 'params.h'
        
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph,eof
        integer dummy1,dummy2,iseg
        character filename*(namelen),buffer*(buffsize)
        
        open(unit=iounit1,file=filename,status='old')
        
        eof=0
        numph=0
        call getline(iounit1,buffer,eof)
        do while ((eof .eq. 0) .and. (numph .lt. maxph))
          numph=numph+1
          read(buffer,*) dummy1,nseg(numph)
          read(buffer,*) dummy1,dummy2,(phaselist(iseg,1,numph),
     &                   phaselist(iseg,2,numph),iseg=1,nseg(numph))
          call getline(1,buffer,eof)
        end do
        
       close(iounit1)        
      
      end
       
c --------------------------------

      subroutine writephases(unit,phaselist,nseg,numph)
      
        implicit none
        include 'params.h'
        
        integer unit,phaselist(maxseg,2,maxph),nseg(maxph),numph
        integer iph,iseg
        
        write(unit,'(A1,2A7,A30)') '#','phase','nseg',
     &                  'layer,phase,layer,phase,...'
        
        do iph=1,numph
          write(unit,'(X,I6,X,I6,A5,$)') iph,nseg(iph),' '
          do iseg=1,nseg(iph)
            write(unit,'(2(X,I2),$)')
     &              phaselist(iseg,1,iph),phaselist(iseg,2,iph)
          end do
          write(unit,*)
        end do
      
      end
       
c ---------------------------------

      subroutine writearrivals(unit,tt,amp,ntr,nph)
      
        implicit none
        include 'params.h'
        
        integer unit,ntr,nph,itr,iph
        real tt(maxph,maxtr),amp(3,maxph,maxtr)
        
        write (unit,'(A1,2A6,4A16)') '#','trace','phase','tt (s)',
     &         'amp n/s','amp e/w','amp z'
        
        do itr=1,ntr
          do iph=1,nph
            write (unit,'(1X,2I6,4(X,G16.9))') itr,iph,tt(iph,itr),
     &             amp(1,iph,itr),amp(2,iph,itr),amp(3,iph,itr)
          end do
        end do
        
      end
      
c -------------------------------------------

      subroutine readweight(filename,weight,ntr)
        
        implicit none
        include 'params.h'
        
        character filename*(namelen),buffer*(buffsize)
        real weight(maxtr)
        integer ntr,ios,itr,eof,dum
        
        open(unit=iounit1,status='old',file=filename,iostat=ios)
        
        if (ios .eq. 0) then
          do itr=1,ntr
            call getline(iounit1,buffer,eof)
            read (buffer,*) dum,weight(itr)
c            write (*,*) weight(itr)
          end do
        else
          write (*,*) 'Weights file nonexistent -- creating.'
          close(unit=iounit1)
          open(unit=iounit1,status='unknown',file=filename)
          write(iounit1,*) '#  Trace, weight'
          do itr=1,ntr
            write(iounit1,*) itr,1.
            weight(itr)=1.
          end do
        end if
        close(unit=iounit1)
        
      end
      
            
c ---------------------------------

      subroutine writetraces(unit,traces,ntr,nsamp,dt,align,shift)
      
        implicit none
        include 'params.h'
c        use csv_file

        integer unit,ntr,align,itr,isamp,nsamp
        real traces(3,maxsamp,maxtr),dt,shift
        
        write(unit,'(A1,5A11)') '#','#traces','#samples',
     &                          'dt (s)','align','shift'
        write(unit,'(X,2I11,F11.3,I11,F11.3)') ntr,nsamp,dt,align,shift
       
        do itr=1,ntr
          write(unit,'(A20)') '#-------------------'
          write(unit,'(A14,1X,I5)') '# Trace number',itr
          write(unit,'(A20)') '#-------------------'
          do isamp=1,nsamp
            a = (traces(1, isamp, itr), traces(2, isamp, itr), traces(3, isamp, itr))
c            call csv_write(1,a,.true.)
c             write (unit,'(3(X,G15.7))') traces(1,isamp,itr),
c     &              traces(2,isamp,itr),traces(3,isamp,itr)
          end do
        end do
        
      end

      
