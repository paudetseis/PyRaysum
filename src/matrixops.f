c####&

c Matrix operation subroutines



c --- generic real matrix linear combination. xA+yB=C, all mXn
      subroutine rmatlincomb(x,A,y,B,C,m,n)
      
        implicit none
        integer i,j,m,n
        real x,A(m,n),y,B(m,n),C(m,n)
        
        do i=1,m
          do j=1,n
            C(i,j)=x*A(i,j)+y*B(i,j)
          end do
        end do
        
      end
      
c --- generic compled matrix linear combination. xA+yB=C, all mXn
      subroutine cmatlincomb(x,A,y,B,C,m,n)
      
        implicit none
        integer i,j,m,n
        complex x,A(m,n),y,B(m,n),C(m,n)
        
        do i=1,m
          do j=1,n
            C(i,j)=x*A(i,j)+y*B(i,j)
          end do
        end do
        
      end

c --- Generic real matrix-constant multiply. A=c*A, A:mxn
      subroutine rmatconst(A,c,m,n)
      
        implicit none
        integer m,n,i,j
        real A(m,n),c
        
        do i=1,m
          do j=1,n
            A(i,j)=c*A(i,j)
          end do
        end do
        
      end
      
      
c --- Generic complex matrix-constant multiply. A=c*A, A:mxn
      subroutine cmatconst(A,c,m,n)
      
        implicit none
        integer m,n,i,j
        complex A(m,n),c
        
        do i=1,m
          do j=1,n
            A(i,j)=c*A(i,j)
          end do
        end do
        
      end
      
      
c --- Generic real matrix multiplication, A*B=C. A:mXn, B:nXo, C:mXo
      subroutine rmatmul(a,b,c,m,n,o)
      
        implicit none
        integer n,m,o,i,j,k
        real a(m,n),b(n,o),c(m,o)

        do i=1,m
          do j=1,o
            c(i,j)=0
            do k=1,n
              c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
          end do
        end do
                
      end
      
      
c --- Generic complex matrix multiplication, A*B=C. A:mXn, B:nXo, C:mXo
      subroutine cmatmul(a,b,c,m,n,o)
      
        implicit none
        integer n,m,o,i,j,k
        complex a(m,n),b(n,o),c(m,o)

        do i=1,m
          do j=1,o
            c(i,j)=0
            do k=1,n
              c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
          end do
        end do
                
      end      
      
      
c --- Generic real matrix-vector multiply, Ax=y. A:mXn,x:n,y:m
      subroutine rmatvec(a,x,y,m,n)
      
        implicit none
        integer m,n,i,j
        real a(m,n),x(n),y(m)
        
        do i=1,m
          y(i)=0
          do j=1,n
            y(i)=y(i)+a(i,j)*x(j)
          end do
        end do
        
      end
      
      
c --- Generic complex matrix-vector multiply, Ax=y. A:mXn,x:n,y:m
      subroutine cmatvec(a,x,y,m,n)
      
        implicit none
        integer m,n,i,j
        complex a(m,n),x(n),y(m)
        
        do i=1,m
          y(i)=0
          do j=1,n
            y(i)=y(i)+a(i,j)*x(j)
          end do
        end do
        
      end
      
        
c --- 3x3 real matrix linear combination. xA+yB=C
      subroutine rmatlincomb3(x,A,y,B,C)
      
        implicit none
        integer i,j
        real x,A(3,3),y,B(3,3),C(3,3)
        
        do i=1,3
          do j=1,3
            C(i,j)=x*A(i,j)+y*B(i,j)
          end do
        end do
        
      end
      
      
c --- 3x3 complex matrix linear combination. xA+yB=C
      subroutine cmatlincomb3(x,A,y,B,C)
      
        implicit none
        integer i,j
        complex x,A(3,3),y,B(3,3),C(3,3)
        
        do i=1,3
          do j=1,3
            C(i,j)=x*A(i,j)+y*B(i,j)
          end do
        end do
        
      end
      
      
c --- 3x3 matrix-constant multiply. A=c*A
      subroutine rmatconst3(A,c)
      
        implicit none
        real A(3,3),c
        integer i,j
        
        do i=1,3
          do j=1,3
            A(i,j)=c*A(i,j)
          end do
        end do
        
      end

      
c --- 3x3 complex matrix-constant multiply. A=c*A
      subroutine cmatconst3(A,c)
      
        implicit none
        complex A(3,3),c
        integer i,j
        
        do i=1,3
          do j=1,3
            A(i,j)=c*A(i,j)
          end do
        end do
        
      end
      
           
c --- 3x3 real matrix multiplication, A*B=C 
      subroutine rmatmul3(a,b,c)
      
        implicit none
        real a(3,3),b(3,3),c(3,3)
        integer i,j,k
        
        do i=1,3
          do j=1,3
            c(i,j)=0
            do k=1,3
              c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
          end do
        end do
                
      end

c --- 3x3 complex matrix multiplication, A*B=C 
      subroutine cmatmul3(a,b,c)
      
        implicit none
        complex a(3,3),b(3,3),c(3,3)
        integer i,j,k
        
        do i=1,3
          do j=1,3
            c(i,j)=0
            do k=1,3
              c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
          end do
        end do
                
      end
   
      
c --- 3x3 real matrix-vector multiplication,Ax=y
      subroutine rmatvec3(a,x,y)
      
        implicit none
        real a(3,3),x(3),y(3)
        integer i,j
        
        do i=1,3
          y(i)=0
          do j=1,3
            y(i)=y(i)+a(i,j)*x(j)
          end do
        end do
        
      end

      
c --- 3x3 complex matrix-vector multiplication,Ax=y
      subroutine cmatvec3(a,x,y)
      
        implicit none
        complex a(3,3),x(3),y(3)
        integer i,j
        
        do i=1,3
          y(i)=0
          do j=1,3
            y(i)=y(i)+a(i,j)*x(j)
          end do
        end do
        
      end
      
c --- 3x3 real/complex matrix-vector multiplication, complex output
      subroutine rcmatvec3(a,x,y)
      
        implicit none
        real a(3,3)
        complex x(3),y(3)
        integer i,j
        
        do i=1,3
          y(i)=0
          do j=1,3
            y(i)=y(i)+cmplx(a(i,j))*x(j)
          end do
        end do
        
      end
      
c --- 3x3 complex/real matrix-vector multiplication, real output
      subroutine crmatvec3(a,x,y)
      
        implicit none
        complex a(3,3)
        real x(3),y(3)
        integer i,j
        
        do i=1,3
          y(i)=0
          do j=1,3
            y(i)=y(i)+real(a(i,j))*x(j)
          end do
        end do
        
      end

c --- 3x3 real vector-matrix multiplication, x'A=y
      subroutine rvecmat3(x,a,y)
      
        implicit none
        real x(3),a(3,3),y(3)
        integer i,j
        
        do i=1,3
          y(i)=0
          do j=1,3
            y(i)=y(i)+x(j)*a(j,i)
          end do
        end do
      
      end
      
      
c --- 3x3 complex vector-matrix multiplication, x'A=y
      subroutine cvecmat3(x,a,y)
      
        implicit none
        complex x(3),a(3,3),y(3)
        integer i,j
        
        do i=1,3
          y(i)=0
          do j=1,3
            y(i)=y(i)+x(j)*a(j,i)
          end do
        end do
      
      end
      
c --- Real dot product, 3-vectors, z=x.y
       real function rdot3(x,y)
      
        implicit none
        real x(3),y(3)
        
        rdot3=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
        
      end
              
c --- Complex dot product, 3-vectors, z=x.y     
      complex function cdot3(x,y)
      
        implicit none
        complex x(3),y(3)
        
        cdot3=conjg(x(1))*y(1)+conjg(x(2))*y(2)+conjg(x(3))*y(3)
        
      end
      

      
c --- Real 3x3 matrix inverse. ai=a^(-1)
c     Adapted from Thompson's routine; based on cofactors.
      subroutine rmatinv3(a,ai)
      
        implicit none
        real a(3,3),ai(3,3),co(3,3),det
        integer i,j
        
        co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
        co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
        co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))
        co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
        co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
        co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))
        co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
        co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
        co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))
        det=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
        do i=1,3
          do j=1,3
            ai(i,j)=co(j,i)/det
          end do
        end do
        
      end

c --- Complex 3x3 matrix inverse. ai=a^(-1)
c     Adapted from Thompson's routine; based on cofactors.
      subroutine cmatinv3(a,ai)
      
        implicit none
        complex a(3,3),ai(3,3),co(3,3),det
        integer i,j
        
        co(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))
        co(1,2)=-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
        co(1,3)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))
        co(2,1)=-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
        co(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))
        co(2,3)=-(a(1,1)*a(3,2)-a(1,2)*a(3,1))
        co(3,1)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))
        co(3,2)=-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
        co(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))
        det=a(1,1)*co(1,1)+a(1,2)*co(1,2)+a(1,3)*co(1,3)
        do i=1,3
          do j=1,3
            ai(i,j)=co(j,i)/det
          end do
        end do
        
      end        
            

c --- Generic real matrix transpose, B=transpose of A. A: mxn,B:nxm
      subroutine rtransp(A,B,m,n)
      
        implicit none
        integer m,n,i,j
        real a(m,n),b(n,m)
        
        do i=1,m
          do j=1,n
            b(j,i)=a(i,j)
          end do
        end do
        
      end
      
c --- Generic complex matrix transpose, B=transpose of A. A: mxn,B:nxm
      subroutine ctransp(A,B,m,n)
      
        implicit none
        integer m,n,i,j
        complex a(m,n),b(n,m)
        
        do i=1,m
          do j=1,n
            b(j,i)=a(i,j)
          end do
        end do
        
      end
      
c --- 3x3 real matrix transpose, B=transpose of A.
      subroutine rtransp3(A,B)
      
        implicit none
        integer i,j
        real a(3,3),b(3,3)
        
        do i=1,3
          do j=1,3
            b(j,i)=a(i,j)
          end do
        end do
        
      end      
   
c --- 3x3 complex matrix transpose, B=transpose of A.
      subroutine ctransp3(A,B)
      
        implicit none
        integer i,j
        complex a(3,3),b(3,3)
        
        do i=1,3
          do j=1,3
            b(j,i)=a(i,j)
          end do
        end do
        
      end      
   


c --- Real matrix copy, arbitrary dimensions. A is copied to B, both mXn
      subroutine rmatcopy(a,b,m,n)
      
        implicit none
        integer m,n,i,j
        real a(m,n),b(m,n)
        
        do i=1,m
          do j=1,n
            b(i,j)=a(i,j)
          end do
        end do
        
      end
      
c --- Complex matrix copy, arbitrary dimensions. A is copied to B, both mXn
      subroutine cmatcopy(a,b,m,n)
      
        implicit none
        integer m,n,i,j
        complex a(m,n),b(m,n)
        
        do i=1,m
          do j=1,n
            b(i,j)=a(i,j)
          end do
        end do
        
      end

c --- 3x3 real matrix copy. A is copied to B.
      subroutine rmatcopy3(a,b)
      
        implicit none
        integer i,j
        real a(3,3),b(3,3)
        
        do i=1,3
          do j=1,3
            b(i,j)=a(i,j)
          end do
        end do
        
      end
      
c --- 3x3 complex matrix copy. A is copied to B.
      subroutine cmatcopy3(a,b)
      
        implicit none
        integer i,j
        complex a(3,3),b(3,3)
        
        do i=1,3
          do j=1,3
            b(i,j)=a(i,j)
          end do
        end do
        
      end
      
     
c --- Real vector copy, arbitrary dimention. x is copied to y, dimension n
      subroutine rveccopy(x,y,n)
      
        implicit none
        integer n,i
        real x(n),y(n)
        
        do i=1,n
          y(i)=x(i)
        end do
        
      end
      
c --- Complex vector copy, arbitrary dimention. x is copied to y, dimension n
      subroutine cveccopy(x,y,n)
      
        implicit none
        integer n,i
        complex x(n),y(n)
        
        do i=1,n
          y(i)=x(i)
        end do
        
      end      
      
c --- Extract column vector v from real matrix A.
c     v=A(i1:i2,k), A:mXn, v:i2-i1+1
      subroutine rextractvec(a,v,i1,i2,k,m,n)
      
        implicit none
        integer i1,i2,k,m,n,i
        real a(m,n),v(i2-i1+1)
        
        do i=i1,i2
          v(i-i1+1)=a(i,k)
        end do
        
      end
      
c --- Extract column vector v from complex matrix A.
c     v=A(i1:i2,k), A:mXn, v:i2-i1+1
      subroutine cextractvec(a,v,i1,i2,k,m,n)
      
        implicit none
        integer i1,i2,k,m,n,i
        complex a(m,n),v(i2-i1+1)
        
        do i=i1,i2
          v(i-i1+1)=a(i,k)
        end do
        
      end      
      
c --- Extract submatrix B from real matrix A. B=A(i1:i2,j1:j2). A is mxn
      subroutine rextractblock(a,b,i1,i2,j1,j2,m,n)
      
        implicit none
        integer i1,i2,j1,j2,m,n,ii,jj
        real a(m,n),b(i2-i1+1,j2-j1+1)
        
        do ii=i1,i2
          do jj=j1,j2
            b(ii-i1+1,jj-j1+1)=a(ii,jj)
          end do
        end do
        
      end
      
c --- Extract submatrix B from complex matrix A. B=A(i1:i2,j1:j2). A is mxn
      subroutine cextractblock(a,b,i1,i2,j1,j2,m,n)
      
        implicit none
        integer i1,i2,j1,j2,m,n,ii,jj
        complex a(m,n),b(i2-i1+1,j2-j1+1)
        
        do ii=i1,i2
          do jj=j1,j2
            b(ii-i1+1,jj-j1+1)=a(ii,jj)
          end do
        end do
        
      end      
      
c --- Find index and value of highest-magnitude (absolute-value max)
c     element in real vector v. v:m elements. In case of equality,
c     finds first element.
      subroutine rmax_vec_abs(v,index,value,m)
      
        implicit none
        integer index,m,j
        real v(m),value
        
        value=abs(v(1))
        index=1
        do j=2,m
          if (abs(v(j)) .gt. value) then
            value=abs(v(j))
            index=j
          end if
        end do
        
      end
        
c --- Find index and value of highest-magnitude (absolute-value max)
c     element in complex vector v. v:m elements. In case of equality,
c     finds first element.
      subroutine cmax_vec_abs(v,index,value,m)
      
        implicit none
        integer index,m,j
        complex v(m)
        real value
        
        value=abs(v(1))
        index=1
        do j=2,m
          if (abs(v(j)) .gt. value) then
            value=abs(v(j))
            index=j
          end if
        end do
        
      end
                
        
c --- Normalize a real 3-vector to 1. y=x/||x||
      subroutine rnorm3(x,y)
      
        implicit none
        real x(3),y(3),norm
        
        norm=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
        if (norm .ne. 0) then
          y(1)=x(1)/norm
          y(2)=y(2)/norm
          y(3)=y(3)/norm
c        else
c          write (*,*) 'rnorm3: Norm is zero'
        end if
        
      end
      
c --- Normalize a complex 3-vector to 1. y=x/||x||
      subroutine cnorm3(x,y)
      
        implicit none
        complex x(3),y(3),norm
        
        norm=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
        if (norm .ne. 0) then
          y(1)=x(1)/norm
          y(2)=y(2)/norm
          y(3)=y(3)/norm
        end if
        
      end       
        
c --- Print a generic real MxN matrix.
      subroutine rprintmat(A,m,n)
      
        implicit none
        integer m,n,i,j
        real A(m,n)
        character*9 fmt
        parameter (fmt='(G14.6,$)')
        
        do i=1,m
          do j=1,n
            write(*,fmt) A(i,j)
          end do
          write(*,*)
        end do
      
      end
        
c --- Print a generic complex MxN matrix.
      subroutine cprintmat(A,m,n)
      
        implicit none
        integer m,n,i,j
        complex A(m,n)
        
        do i=1,m
          do j=1,n
            write(*,'(X,''('',G14.6,'','',G14.6,'')'',$)')
     &             real(A(i,j)),aimag(A(i,j))
          end do
          write(*,*)
        end do
      
      end
