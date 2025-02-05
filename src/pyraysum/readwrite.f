c####&


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
          write(*,*) 'phase: ', iph
          write(*,*) 'layer: type:'
          do iseg=1,nseg(iph)
            write(*,'(5X,I2,5X,I1)') phaselist(iseg,1,iph),
     &                                  phaselist(iseg,2,iph)
          end do
          write(*,*)
        end do
       end

            


      
