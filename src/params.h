c####&
c namelen is the length of filenames (in characters)
c maxlay is the maximum number of layers allowed in the model
c maxtr is the maximum number of traces allowed
c maxseg: maximum # of segments (should be 3*maxlay for 1st-order
c         multiples
c maxph: maximum number of phases per trace
c buffsize is the max. line length assumed for reading files.
      integer namelen, maxlay, maxtr, maxseg, maxph, buffsize
      parameter (namelen=40,maxlay=15,maxtr=200,maxseg=45)
      parameter (maxph=40000,buffsize=120)
      
c Units for reading and writing
      integer iounit1,iounit2
      parameter (iounit1=1,iounit2=2)
      
c pi: duh. ztol: tolerance considered to be equiv. to zero.
      real pi,ztol
      parameter (pi=3.141592653589793,ztol=1.e-7)
      
c nsamp is the number of samples per trace.
      integer maxsamp
      parameter (maxsamp=2000)
