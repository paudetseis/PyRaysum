! Copyright 2020 Andrew Frederiksen and Pascal Audet

! This file is part of PyRaysum.

! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modIFy, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:

! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!===========================================================================
!
! MODULE conf
!
! Configuration module that contains global variables used in rmat and plane 
! modules to interface with the Python codes.
!
!===========================================================================

      MODULE conf

      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: pi = 3.141592653589793d0
      INTEGER, PARAMETER :: nlaymx = 30
!
! Model parameters
!
      DOUBLE PRECISION :: a(3,3,3,3,nlaymx), thickn(nlaymx)
      DOUBLE PRECISION :: rho(nlaymx)
      INTEGER :: isoflg(nlaymx)
!
! Wavefield parameters
!
      DOUBLE PRECISION :: dt, slow, baz

      END MODULE conf


