      real function getmisfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3,mtype)
     
        implicit none
        include 'params.h'
       
        real data_tr(3,maxsamp,maxtr),syn_tr(3,maxsamp,maxtr)
        real weight(maxtr),m1,m2,m3
        integer ntr,nsamp,mtype

        real l1_misfit,l1sc_misfit,l2_misfit,l2sc_misfit,corr_misfit
        
        if (mtype .eq. 0) then
          getmisfit=corr_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
        end if
        if (mtype .eq. 1) then
          getmisfit=l1sc_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
        end if
        if (mtype .eq. 2) then
          getmisfit=l2sc_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
        end if
        if (mtype .eq. 3) then
          getmisfit=l1_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
        end if
        if (mtype .eq. 4) then
          getmisfit=l2_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
        end if
        
      end

c ----------------------------------------------

      real function l1_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
c L1-norm misfit between two sets of traces
      
        implicit none
        include 'params.h'
        
        real data_tr(3,maxsamp,maxtr),syn_tr(3,maxsamp,maxtr)
        real weight(maxtr),m1,m2,m3
        integer ntr,nsamp,itr,isamp
        real diff,trdiff,wtot
        
        
        diff=0.
        wtot=0.
        do itr=1,ntr
          trdiff=0.
          do isamp=1,nsamp
            trdiff=trdiff+m1*
     &            abs(data_tr(1,isamp,itr)-syn_tr(1,isamp,itr))
            trdiff=trdiff+m2*
     &            abs(data_tr(2,isamp,itr)-syn_tr(2,isamp,itr))
            trdiff=trdiff+m3*
     &            abs(data_tr(3,isamp,itr)-syn_tr(3,isamp,itr))
          end do
          diff=diff+weight(itr)*trdiff
          wtot=wtot+weight(itr)
        end do
        
        l1_misfit=diff/(wtot*(m1+m2+m3)*real(nsamp))
        
      end
      
      
c ------------------------------------------------------

      real function l1sc_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
c Scaled L1-norm misfit between two sets of traces
      
        implicit none
        include 'params.h'
        
        real data_tr(3,maxsamp,maxtr),syn_tr(3,maxsamp,maxtr)
        real weight(maxtr),m1,m2,m3,m(3)
        integer ntr,nsamp,itr,isamp,j
        real diff,trdiff,wtot
        
        m(1)=m1
        m(2)=m2
        m(3)=m3
        
        call scaletr(data_tr,ntr,nsamp,m,1)
        call scaletr(syn_tr,ntr,nsamp,m,1)
c        call scaletr(data_tr,ntr,nsamp,m,2)
c        call scaletr(syn_tr,ntr,nsamp,m,2)
        
        diff=0.
        wtot=0.
        do itr=1,ntr
          trdiff=0.
          do isamp=1,nsamp
            do j=1,3
              trdiff=trdiff+m(j)*
     &          abs(data_tr(j,isamp,itr)-syn_tr(j,isamp,itr))
            end do
          end do
          diff=diff+weight(itr)*trdiff
          wtot=wtot+weight(itr)
        end do
        
        l1sc_misfit=diff/(wtot*(m1+m2+m3)*real(nsamp))
        
      end
      
      
c ------------------------------------------------------
      real function l2_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
c L2-norm misfit between two sets of traces
      
        implicit none
        include 'params.h'
        
        real data_tr(3,maxsamp,maxtr),syn_tr(3,maxsamp,maxtr)
        real weight(maxtr),m1,m2,m3
        integer ntr,nsamp,itr,isamp
        real diff,trdiff,wtot
        
        diff=0.
        wtot=0.
        do itr=1,ntr
          trdiff=0.
          do isamp=1,nsamp
            trdiff=trdiff+m1*
     &            (data_tr(1,isamp,itr)-syn_tr(1,isamp,itr))**2
            trdiff=trdiff+m2*
     &            (data_tr(2,isamp,itr)-syn_tr(2,isamp,itr))**2
            trdiff=trdiff+m3*
     &            (data_tr(3,isamp,itr)-syn_tr(3,isamp,itr))**2
          end do
          diff=diff+weight(itr)*sqrt(trdiff)
          wtot=wtot+weight(itr)
        end do
        
        l2_misfit=diff/(wtot*sqrt(m1+m2+m3)*real(nsamp))
        
      end
      
c --------------------------------
      
      real function l2sc_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
c Scaled L2-norm misfit between two sets of traces
      
        implicit none
        include 'params.h'
        
        real data_tr(3,maxsamp,maxtr),syn_tr(3,maxsamp,maxtr)
        real weight(maxtr),m1,m2,m3,m(3)
        integer ntr,nsamp,itr,isamp,j
        real diff,trdiff,wtot 
               
        m(1)=m1
        m(2)=m2
        m(3)=m3
        call scaletr(data_tr,ntr,nsamp,m,2)
        call scaletr(syn_tr,ntr,nsamp,m,2)
        
        diff=0.
        wtot=0.
        do itr=1,ntr
          trdiff=0.
          do isamp=1,nsamp
            do j=1,3
              trdiff=trdiff+m(j)*
     &          (data_tr(j,isamp,itr)-syn_tr(j,isamp,itr))**2
            end do
          end do
          diff=diff+weight(itr)*sqrt(trdiff)
          wtot=wtot+weight(itr)
        end do
        
        l2sc_misfit=diff/(wtot*sqrt(m1+m2+m3)*real(nsamp))
        
      end
      
      
c -------------------------------------------------------

      real function corr_misfit(data_tr,syn_tr,weight,ntr,
     &                          nsamp,m1,m2,m3)
     
c A correlation-based misfit function
        
        implicit none
        include 'params.h'
        
        real data_tr(3,maxsamp,maxtr),syn_tr(3,maxsamp,maxtr)
        real weight(maxtr),m1,m2,m3
        integer ntr,nsamp,itr,isamp
        real prod,d1,d2,d3,s1,s2,s3
        real sprod,dprod,corr
        
        prod=0.
        sprod=0.
        dprod=0.
        do itr=1,ntr
          do isamp=1,nsamp
            d1=data_tr(1,isamp,itr)
            s1=syn_tr(1,isamp,itr)
            d2=data_tr(2,isamp,itr)
            s2=syn_tr(2,isamp,itr)
            d3=data_tr(3,isamp,itr)
            s3=syn_tr(3,isamp,itr)
            prod=prod + (m1*d1*s1 + m2*d2*s2 + m3*d3*s3)*weight(itr)
            sprod=sprod + (m1*s1*s1 + m2*s2*s2 + m3*s3*s3)*weight(itr)
            dprod=dprod + (m1*d1*d1 + m2*d2*d2 + m3*d3*d3)*weight(itr)
          end do
        end do
        
c        corr=1 for perfect correlation, -1 for perfect
c        anticorrelation
        corr=prod/(sqrt(sprod)*sqrt(dprod))
        
        corr_misfit=1-corr
        
      end
      

c ----------------------------

      subroutine scaletr(tr,ntr,nsamp,m,l)
      
        implicit none
        include 'params.h'
        
        real tr(3,maxsamp,maxtr),m(3),norm
        integer ntr,nsamp,l,isamp,itr,j
        
        do itr=1,ntr
          norm=0.
          do isamp=1,nsamp
            do j=1,3
              norm=norm+abs(tr(j,isamp,itr))**l
            end do
          end do
          norm=norm**(1./real(l))
          do isamp=1,nsamp
            do j=1,3
              tr(j,isamp,itr)=tr(j,isamp,itr)/norm
            end do
          end do
        end do
        
      end
