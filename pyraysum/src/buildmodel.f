c####&

      subroutine buildmodel(a,ar,rot,thick,rho,alpha,beta,isoflag,
     &                      pct_a,pct_b,trend,plunge,strike,dip,nlay)
c Set up model for use.
     
        implicit none
        include 'params.h'
        
c Model variables:
        integer nlay        
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct_a(maxlay),pct_b(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay)
        logical isoflag(maxlay)
        
c Outputs:
        real a(3,3,3,3,maxlay),ar(3,3,3,3,maxlay,2),rot(3,3,maxlay)
        
c Scratch variables:
        integer ilay
        real rot_axis(3,3),d_a,d_b,AA,CC,LL,NN,FF,a_temp(3,3,3,3)
        
c Parameter eta, combined with percentages of P and S anisotropy,
c is sufficient to obtain all coefficients for a hexagonally-symmetric
c medium. Fixed from Farra et al., 1991
        real eta
        parameter (eta=1.03)
                
        do ilay=1,nlay
        
          if (isoflag(ilay)) then
c Store isotropic coefficients (only 2 are used)
            a(3,3,3,3,ilay) = alpha(ilay)**2.
            a(2,3,2,3,ilay) = beta(ilay)**2.
          else
c Build anisotropic coefficients for a hexagonally-symmetric medium,
c after Farra et al, 1991
            d_a=alpha(ilay)*pct_a(ilay)/100.
            d_b=beta(ilay)*pct_b(ilay)/100.
            AA=rho(ilay)*(alpha(ilay) - d_a/2.)**2.
            CC=rho(ilay)*(alpha(ilay) + d_a/2.)**2.
            LL=rho(ilay)*(beta(ilay) + d_b/2.)**2.
            NN=rho(ilay)*(beta(ilay) - d_b/2.)**2.
            FF=eta*(AA-2.*LL)
c            write (*,*) 'Anis. params are:',AA,CC,FF,LL,NN
c Get tensor with unrotated axes
            call tritensr(a_temp,AA,CC,FF,LL,NN,rho(ilay))
c Rotate axes:
            rot_axis(1,1)=cos(trend(ilay))*cos(plunge(ilay))
            rot_axis(2,1)=-sin(trend(ilay))
            rot_axis(3,1)=-cos(trend(ilay))*sin(plunge(ilay))
            rot_axis(1,2)=sin(trend(ilay))*cos(plunge(ilay))
            rot_axis(2,2)=cos(trend(ilay))
            rot_axis(3,2)=-sin(trend(ilay))*sin(plunge(ilay))
            rot_axis(1,3)=sin(plunge(ilay))
            rot_axis(2,3)=0.
            rot_axis(3,3)=cos(plunge(ilay))
c            write (*,*) 'rot_axis:',rot_axis
            call rot_tensor(a(1,1,1,1,ilay),a_temp,rot_axis)
            
          end if                      
        end do
        
c        write (*,*) 'build-model: a(:,:,2,3,2)'
c        do ilay=1,3
c          write (*,*) a(ilay,1,2,3,2),a(ilay,2,2,3,2),a(ilay,3,2,3,2)
c        end do
               
c Make rotator-matrix list:
        call make_rotator(rot,strike,dip,nlay)
c        write (*,*) 'build-model: rot(:,:,3)'
c        do ilay=1,3
c          write (*,*) rot(ilay,1,3),rot(ilay,2,3),rot(ilay,3,3)
c        end do        
        
        
c Pre-rotate tensors (skip isotropic ones)
        do ilay=1,nlay
          if (.not.(isoflag(ilay))) then
c            Upper interface:
            call rot_tensor(ar(1,1,1,1,ilay,2),a(1,1,1,1,ilay),
     &           rot(1,1,ilay))
c            Lower interface. Bottom layer is a half-space.
            if (ilay .lt. nlay) then
              call rot_tensor(ar(1,1,1,1,ilay,1),a(1,1,1,1,ilay),
     &             rot(1,1,ilay+1))
            end if
          end if
        end do
        
c        write (*,*) 'build-model: ar(:,:,2,3,2,1)'
c        do ilay=1,3
c          write (*,*) ar(ilay,1,2,3,2,1),ar(ilay,2,2,3,2,1),
c     &                ar(ilay,3,2,3,2,1)
c        end do       
      
      end
      
      
c -----------------------------
      
      
      subroutine tritensr(a,AA,CC,FF,LL,NN,rho)
      
c Sets up transversely isotropic tensor with horizontal
c symmetry axis.
      
        implicit none
        
        integer i,j,k,l
        real a(3,3,3,3),AA,CC,FF,LL,NN,rho
        
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                a(i,j,k,l)=0
              end do
            end do
          end do
        end do
        
        
        a(3,3,3,3)=AA/rho
        a(2,2,2,2)=AA/rho
        a(1,1,1,1)=CC/rho

        a(3,3,2,2)=(AA-2*NN)/rho
        a(2,2,3,3)=(AA-2*NN)/rho

        a(3,3,1,1)=FF/rho
        a(1,1,3,3)=FF/rho

        a(2,2,1,1)=FF/rho
        a(1,1,2,2)=FF/rho

        a(2,1,2,1)=LL/rho
        a(1,2,1,2)=LL/rho
        a(1,2,2,1)=LL/rho
        a(2,1,1,2)=LL/rho

        a(1,3,1,3)=LL/rho
        a(3,1,3,1)=LL/rho
        a(1,3,3,1)=LL/rho
        a(3,1,1,3)=LL/rho
 
        a(3,2,3,2)=NN/rho
        a(2,3,2,3)=NN/rho
        a(3,2,2,3)=NN/rho
        a(2,3,3,2)=NN/rho
      
      end
    

c ---------


      subroutine rot_tensor(CR,CC,R)
c Rotate tensor CC according to rotator matrix R. Result is output
c as tensor CR.
      
        implicit none
        
        real CC(3,3,3,3),CR(3,3,3,3),R(3,3)
        integer i,j,k,l,a,b,c,d
        
        do i=1,3
        do j=1,3
        do k=1,3
        do l=1,3
          
          CR(i,j,k,l)=0
          
          do a=1,3
          do b=1,3
          do c=1,3
          do d=1,3
          
            CR(i,j,k,l)=CR(i,j,k,l) + R(a,i)*R(b,j)*R(c,k)*R(d,l)
     &                  * CC(a,b,c,d)
     
          end do
          end do
          end do
          end do

        end do
        end do
        end do
        end do
        
      end
      
      
c ---------------------------

      subroutine make_rotator(R,strike,dip,nlay)
c Build list of rotator matrices from a list of interface strikes
c and dips. See make_rotator.m for more info.

      
        implicit none
        include 'params.h'
        
        real R(3,3,maxlay),strike(maxlay),dip(maxlay)
        integer nlay,layer,i,j
        
c Topmost rotator is the identity matrix

        do i=1,3
          do j=1,3
            if (i.eq.j) then
              R(i,j,1)=1
            else
              R(i,j,1)=0
            end if
          end do
        end do
        
c Build other rotators
        
        do layer=2,nlay
          
          R(1,1,layer)=cos(strike(layer))
          R(2,1,layer)=sin(strike(layer))
          R(3,1,layer)=0
          
          R(1,2,layer)=-cos(dip(layer))*sin(strike(layer))
          R(2,2,layer)=cos(dip(layer))*cos(strike(layer))
          R(3,2,layer)=sin(dip(layer))
          
          R(1,3,layer)=sin(dip(layer))*sin(strike(layer))
          R(2,3,layer)=-sin(dip(layer))*cos(strike(layer))
          R(3,3,layer)=cos(dip(layer))
          
        end do
        
        
c        do layer=1,nlay
c          write (*,*) 'Rotator for',layer,strike(layer),dip(layer)
c          write (*,*) R(1,1,layer),R(1,2,layer),R(1,3,layer)
c          write (*,*) R(2,1,layer),R(2,2,layer),R(2,3,layer)
c          write (*,*) R(3,1,layer),R(3,2,layer),R(3,3,layer)
c        end do
      end
          
      
      
      
      
      
      
