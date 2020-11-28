c####&

      program test
      
        implicit none
        include 'params.h'
        
        character filename*(namelen)
        integer phaselist(maxseg,2,maxph),nseg(maxph),numph
        
        filename='../Forward/Sample/sample.ph'
        call readphases(filename,phaselist,nseg,numph)  
        call writephases(6,phaselist,nseg,numph)
        
      end
