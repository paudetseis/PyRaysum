c####&

c Eigenvector/eigenvalue problem for isotropic and anisotropic layers.
c Real values only.

      subroutine isotroc(a,rho,p1,p2,eval,evec)
c   Obtain isotropic eigenvectors/values. This is a little obscure,
c   but it's just the Zoeppritz equations.
      
        implicit none
        real a(3,3,3,3),rho,p1,p2
        complex eval(6),evec(6,6)
        complex vp2,vs2,mu,pp,p,qdp,qds,qup,qus
        real xnorm
        integer i,j
        
c        write (*,*) 'isotroc, (p1,p2):',p1,p2
        
c   Load up some convenient values
        vp2=a(3,3,3,3)
        vs2=a(2,3,2,3)
        mu=rho*vs2
        pp=p1*p1+p2*p2
        p=sqrt(pp)

c    Set eigenvalues
c        write (*,*) 1./vp2,1./vs2,pp       
        qdp=sqrt(1/vp2-pp)
        qds=sqrt(1/vs2-pp)
        qup=-qdp
        qus=-qds
        eval(1)=qdp
        eval(2)=qds
        eval(3)=qds
        eval(4)=qup
        eval(5)=qus
        eval(6)=qus
        
c    Set eigenvector matrix, column by column,
c      N(1:6,1)=[p1,p2,qdp,2*mu*p1*qdp,2*mu*p2*qdp,(rho-2*mu*pp)]'
        evec(1,1)=p1
        evec(2,1)=p2
        evec(3,1)=qdp
        evec(4,1)=2.*mu*p1*qdp
        evec(5,1)=2.*mu*p2*qdp
        evec(6,1)=rho-2.*mu*pp
c      N(1:6,2)=[p1,p2,-pp/qds,p1*N(6,1)/qds,p2*N(6,1)/qds,-2*mu*pp]';
        evec(1,2)=p1
        evec(2,2)=p2
        evec(3,2)=-pp/qds
        evec(4,2)=p1*evec(6,1)/qds
        evec(5,2)=p2*evec(6,1)/qds
        evec(6,2)=-2.*mu*pp
c      N(1:6,3)=[-p2,p1,0,-p2*qds*mu,p1*qds*mu,0]'
        evec(1,3)=-p2
        evec(2,3)=p1
        evec(3,3)=0.
        evec(4,3)=-p2*qds*mu
        evec(5,3)=p1*qds*mu
        evec(6,3)=0.
c      N(1:6,4)=[p1,p2,qup,2*mu*p1*qup,2*mu*p2*qup,N(6,1)]';
        evec(1,4)=p1
        evec(2,4)=p2
        evec(3,4)=qup
        evec(4,4)=2.*mu*p1*qup
        evec(5,4)=2.*mu*p2*qup
        evec(6,4)=evec(6,1)
c      N(1:6,5)=[p1,p2,-pp/qus,p1*N(6,1)/qus,p2*N(6,1)/qus,-2*mu*pp]';
        evec(1,5)=p1
        evec(2,5)=p2
        evec(3,5)=-pp/qus
        evec(4,5)=p1*evec(6,1)/qus
        evec(5,5)=p2*evec(6,1)/qus
        evec(6,5)=-2.*mu*pp
c      N(1:6,6)=[-p2,p1,0,-p2*qus*mu,p1*qus*mu,0]'
        evec(1,6)=-p2
        evec(2,6)=p1
        evec(3,6)=0.
        evec(4,6)=-p2*qus*mu
        evec(5,6)=p1*qus*mu
        evec(6,6)=0.
        
c      Normalize wrt displacement magnitude:
        do j=1,6         
          xnorm=sqrt(real(evec(1,j))**2+real(evec(2,j))**2+
     &               real(evec(3,j))**2)
          do i=1,6
            evec(i,j)=evec(i,j)/xnorm
          end do
        end do
             
c        write(*,*) 'Isotroc: eigens'
c        write(*,*) eval
c        do i=1,6
c          write (*,*) real(evec(i,1)),real(evec(i,2)),real(evec(i,3)),
c     &                real(evec(i,4)),real(evec(i,5)),real(evec(i,6))
c        end do
        
      end
      
      
c ---------------------------------------

      subroutine anisotroc(a,rho,p1,p2,eval,evec)
c   Obtain anisotropic eigenvectors/values. Rather hairy.
      
        implicit none
        real a(3,3,3,3),rho,p1,p2
        complex eval(6),evec(6,6)
                
c        integer worksize
c        parameter (worksize=100)
        real AA(6,6),CC(3,3,3,3),iC33(3,3),T(3,3),S(3,3),p(3)
        real wrk1(3,3),wrk2(3,3),wrk3(3,3),wrk4(3,3),eye(3,3)
        real zero6(6,6),evalr(6),evecr(6,6),evali(6),eveci(6,6)
c        real dummy(1,6),evec_la(6,6),work_la(worksize)
        integer i,j,k,l,ierr
        real xnorm
        
c       Identity matrix
        data eye(1,1)/1./,eye(1,2)/0./,eye(1,3)/0./
        data eye(2,1)/0./,eye(2,2)/1./,eye(2,3)/0./
        data eye(3,1)/0./,eye(3,2)/0./,eye(3,3)/1./
        
c        write (*,*) 'anisotroc, (p1,p2):',p1,p2
       
c       Build partion matrices CIJ where CIJ(k,l)=rho*a(k,i,l,j)
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                CC(k,l,i,j)=rho*a(k,i,l,j)
              end do
            end do
          end do
        end do
        
c        write (*,*) 'anisotroc -- c(:,:,1,1)'
c        do i=1,3
c          write (*,*) a(1,i,1,1)*rho,a(2,i,1,1)*rho,a(3,i,1,1)
c        end do
       
        call rmatinv3(CC(1,1,3,3),iC33)
        
c       T=(-p1*C13-p2*C23)*iC33
        call rmatlincomb3(-p1,CC(1,1,1,3),-p2,CC(1,1,2,3),wrk1)
        call rmatmul3(wrk1,iC33,T)
        
c       S=rho*eye(3,3)-sum(i=1:2,j=1:2,p(i)*p(j)*(Cij-Ci3*C33^(-1)*C3j))
        p(1)=p1
        p(2)=p2
        do i=1,3
          do j=1,3
            wrk1(i,j)=0
          end do
        end do
        do j=1,2
          call rmatmul3(iC33,CC(1,1,3,j),wrk2)
          do i=1,2
            call rmatmul3(CC(1,1,i,3),wrk2,wrk3)
            call rmatlincomb3(1.,CC(1,1,i,j),-1.,wrk3,wrk4)
            call rmatlincomb3(1.,wrk1,p(i)*p(j),wrk4,wrk3)
            call rmatcopy3(wrk3,wrk1)
          end do
        end do
        call rmatlincomb3(rho,eye,-1.,wrk1,S)
        
c       Build system matrix AA=[T',iC33;S,T]
        do i=1,3
          do j=1,3
            AA(i,j)=T(j,i)
            AA(i,j+3)=iC33(i,j)
            AA(i+3,j)=S(i,j)
            AA(i+3,j+3)=T(i,j)
          end do
        end do

c        write (*,*) 'anisotroc: system sub-matrices'
c        write (*,*) 'T='
c        do i=1,3
c          write (*,*) T(i,1),T(i,2),T(i,3)
c        end do
c        write (*,*) 'iC33='
c        do i=1,3
c          write (*,*) iC33(i,1),iC33(i,2),iC33(i,3)
c        end do
c        write (*,*) 'S='
c        do i=1,3
c          write (*,*) S(i,1),S(i,2),S(i,3)
c        end do
        
c        write (*,*) 'anisotroc: system matrix'
c        do i=1,6
c          write (*,*) AA(i,1),AA(i,2),AA(i,3),AA(i,4),AA(i,5),AA(i,6)
c        end do
        
c       Obtain eigenvalues/vectors using EISPACK.
        do i=1,6
          do j=1,6
            zero6(i,j)=0.
          end do
        end do
        call cg(6,6,AA,zero6,evalr,evali,1,evecr,eveci,
     &          wrk1,wrk2,wrk3,ierr)
        if (ierr .ne. 0) then
          write (*,*) '*** Error from EISPACK cg : ',ierr
        end if
        

c       Obtain eigenvalues/vectors using LAPACK. Hopefully will
c       give better results than EISPACK,although the way that it
c       stores eigenvectors is supremely annoying.
c        call sgeev('N','V',6,AA,6,evalr,evali,dummy,1,evec_la,6,
c     &       work_la,worksize,ierr)
c        if (ierr .ne. 0) then
c          write (*,*) '*** Error from LAPACK sgeev : ',ierr
c        end if
c        j=1
cc       Rearrange evecs in a more reasonable fashion
c        do while (j .le. 6)
c          if (evali(j) .eq. 0.) then
c            do i=1,6
c              evecr(i,j)=evec_la(i,j)
c              eveci(i,j)=0.
c            end do
c            j=j+1
c          else
cc          Complex eigens are in conjugate pairs -- rearrange and skip
c            do i=1,6
c              evecr(i,j)=evec_la(i,j)
c              evecr(i,j+1)=evec_la(i,j)
c              eveci(i,j)=evec_la(i,j+1)
c              eveci(i,j+1)=-evec_la(i,j+1)
c            end do
c            j=j+2
c          end if
c        end do
        
c        write (*,*) 'Anisotroc: unsorted eigenvalues:',evalr
c        write (*,*) 'Anisotroc: unsorted eigenvectors'
c        do i=1,6
c          write (*,*) evecr(i,1),evecr(i,2),evecr(i,3),evecr(i,4),
c     &                evecr(i,5),evecr(i,6)
c        end do        
                
c        Sort evecs/evals
        call sort_evec(evalr,evali,evecr,eveci,eval,evec)
            
c        Normalize evecs to unit displacement
        do j=1,6
          xnorm=sqrt(real(evec(1,j))**2+real(evec(2,j))**2+
     &               real(evec(3,j))**2)
          do i=1,6
            evec(i,j)=evec(i,j)/xnorm
          end do
        end do
        
c        write (*,*) 'Anisotroc: sorted, normalized eigens'
c        write (*,*) 'evals:',eval
c        write (*,*) 'evec:'
c        do i=1,6
c          write (*,*) real(evec(i,1)),real(evec(i,2)),real(evec(i,3)),
c     &                real(evec(i,4)),real(evec(i,5)),real(evec(i,6))
c        end do        
        
      end
      
      
c --------------------

c Sort eigenvectors in the same manner as the MATLAB anisotroc.m
c routine. If all evals are real, the order will be like [1 2 3
c -1 -2 3]; otherwise, evanescent waves go first, e.g. [2i,1+2i,
c 1 -2i, -1-2i, -1].
      subroutine sort_evec(evalr,evali,evecr,eveci,eval,evec)
      
        implicit none
        include 'params.h'
        
        real evalr(6),evali(6),evecr(6,6),eveci(6,6)
        complex eval(6),evec(6,6)
        integer i,index(6),imagpos(6),imagneg(6),realpos(6),realneg(6)
        integer nrp,nrn,nip,nin
        integer j
        
c Got burned by using a DATA statement here...
        nrp=0
        nrn=0
        nip=0
        nin=0        
        
c Divide eigenvalues up into real positive, real negative, complex
c positive, complex negative.
        do i=1,6
c         Is it real?
          if (abs(evali(i)/evalr(i)) .lt. ztol*100) then
            evali(i)=0
            do j=1,6
              eveci(j,i)=0
            end do
            if (evalr(i) .ge. 0.) then
              nrp=nrp+1
              realpos(nrp)=i
            else
              nrn=nrn+1
              realneg(nrn)=i
            end if
          else
            if (evali(i) .ge. 0.) then
              nip=nip+1
              imagpos(nip)=i
            else
              nin=nin+1
              imagneg(nin)=i
            end if
          end if
        end do
        
c        print *,nip,nrp,nin,nrn

c        Sort sub-groups
        call sort(evalr,realpos,6,nrp)
        call sort(evalr,realneg,6,nrn)
        call sort(evali,imagpos,6,nip)
        call sort(evali,imagneg,6,nin)
        
c        Assemble sub-ranges
        do i=1,nip
          index(i)=imagpos(i)
        end do
        do i=1,nrp
          index(i+nip)=realpos(i)
        end do
        do i=1,nin
          index(i+nip+nrp)=imagneg(nin-i+1)
        end do
        do i=1,nrn
          index(i+nip+nrp+nin)=realneg(nrn-i+1)
        end do

c        print *,'sort-evec:',evalr
c        print *,'sort-evec:',evali
c        print *,'sort-evec:',index

c        write (*,*) 'sort-evec: index is',index
        
        call reorder_evec(evalr,evali,evecr,eveci,index,eval,evec)
          
      end

c ----------------------------------------------------------

      subroutine reorder_evec(evalr,evali,evecr,eveci,index,eval,evec)
        
        implicit none
        real evalr(6),evali(6),evecr(6,6),eveci(6,6)
        integer index(6),j,k
        complex eval(6),evec(6,6)
        
        do j=1,6
          eval(j)=cmplx(evalr(index(j)),evali(index(j)))
          do k=1,6
            evec(k,j)=cmplx(evecr(k,index(j)),eveci(k,index(j)))
          end do
        end do
        
      end        
        
c ----------------------------------------------------------

c Sort an array of n real-valued elements, from smallest (most negative)
c to largest (most positive). Actual elements aren't moved -- index
c returns new positions. Uses a basic algorithm (insertion sort)
c since we're sorting few elements. m is the number of elements in a;
c n is the number of elements in index (n <= m).
c
c Index should be pre-initialized with initial locations.
        
      subroutine sort(a,index,m,n)
      
        implicit none
        integer m,n
        integer index(n),j,pos,i
        real a(m),value
        
        if (n .le. 1) then
          return
        end if
        
        do j=2,n
          i=index(j)
          value=a(i)          
          pos=j
          do while ((pos .gt. 1) .and. (value .lt. a(index(pos-1))))
            index(pos)=index(pos-1)
            pos=pos-1
          end do
          index(pos)=i
        end do
        
      end
 





