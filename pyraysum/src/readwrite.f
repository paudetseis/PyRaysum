c####&


      subroutine readmodel(filename,thick,rho,alpha,beta,isoflag,
     &                     pct,trend,plunge,strike,dip,nlay)
     
c Read in a model file from disk.
          
        implicit none
        include 'params.h'
        
        character filename*(namelen), buffer*(buffsize)
        integer eof,nlay,i,isoint
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay)
        logical isoflag(maxlay)
        
        open(unit=iounit1,file=filename,status='old')
        
c Read in lines from file, skipping lines that start with '#' (comments)
	      eof=0
        nlay=0
        call getline(iounit1,buffer,eof)
        do while ((eof .eq. 0) .and. (nlay .lt. maxlay))
          nlay=nlay+1
          read(buffer,*) thick(nlay),rho(nlay),alpha(nlay),beta(nlay),
     &                   isoint,pct(nlay),trend(nlay),plunge(nlay),
     &                   strike(nlay),dip(nlay)
          if (isoint .eq. 0) then
            isoflag(nlay) = .false.
          else
            isoflag(nlay) = .true.
          end if
          call getline(1,buffer,eof)
        end do 
               
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
        
c Read in lines from file, skipping comments
	eof=0
        ntr=0
        call getline(iounit1,buffer,eof)
        do while ((eof .eq. 0) .and. (ntr .lt. maxtr))
          ntr=ntr+1
          read(buffer,*) baz(ntr),slow(ntr),sta_dx(ntr),sta_dy(ntr)
          call getline(1,buffer,eof)
c Convert to radians
          baz(ntr)=baz(ntr)/180. * pi
        end do 
        close(iounit1)
        
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
     &                      pct,trend,plunge,strike,dip,nlay)
     
        implicit none
        include 'params.h'
        
        integer unit, nlay, i, isonum
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay)
        logical isoflag(maxlay)
        
        write(unit,'(A1,4A7,A4,5A7)') '#','thick','rho','alpha','beta',
     &        'iso','%','trend','plunge','strike','dip'
        
        do i=1,nlay
          if (isoflag(i)) then
            isonum=1
          else
            isonum=0
          end if
          write(unit,'(1X,4F7.0,I4,F7.1,4F7.0)')
     &          thick(i),rho(i),alpha(i),beta(i),isonum,
     &          pct(i),trend(i)*180/pi,plunge(i)*180/pi,
     &          strike(i)*180./pi,dip(i)*180./pi
        end do
             
      end
      
c -----------------------------------
      
      subroutine writegeom(unit,baz,slow,sta_dx,sta_dy,ntr)
      
        implicit none
        include 'params.h'
        
        real baz(maxtr),slow(maxtr),sta_dx(maxtr),sta_dy(maxtr)
        integer unit,ntr,i
        
        write(unit,'(A1,A7,3A10)') '#','baz','slow','dx','dy'
        
        do i=1,ntr
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

            
c ---------------------------------

      subroutine writetraces(unit,traces,ntr,nsamp,dt,align,shift)
      
        implicit none
        include 'params.h'

        integer unit,ntr,align,itr,isamp,nsamp
        real traces(3,maxsamp,maxtr),dt,shift
        
        write(unit, '(7A)') 'itr',',','trace1',',','trace2',
     &            ',','trace3'

        do itr=1,ntr
          do isamp=1,nsamp
            write (unit,'(I5,A,G15.7,A,G15.7,A,G15.7)') itr-1,',',
     &              traces(1,isamp,itr),
     &              ',',traces(2,isamp,itr),',',traces(3,isamp,itr)
          end do
        end do
        
      end


      
