c####&

      program seis_misfit
      
        implicit none
        include 'params.h'
        
        integer nargs,mtype
        character cur_arg*(namelen),fname_1*(namelen),fname_2*(namelen)
        character weight_name*(namelen)
        real m1,m2,m3
        
        integer nsamp1,ntr1,nsamp2,ntr2,itr
        real data_1(3,maxsamp,maxtr),data_2(3,maxsamp,maxtr)
        real dt,align,shift,trweight(maxtr)
        
        real getmisfit
        integer iargc
        
        nargs=iargc()
        if (nargs .lt. 6) then
          write (*,*) 'Usage: seis-misfit traces_1 traces_1 type',
     &                ' w1 w2 w3 [trace_weights]'
          stop
        end if
        
        call getarg(1,fname_1)
        call readtraces(fname_1,data_1,ntr1,nsamp1,dt,align,shift)
        call getarg(2,fname_2)
        call readtraces(fname_2,data_2,ntr2,nsamp2,dt,align,shift)
        if ((ntr1 .ne. ntr2) .or. (nsamp1 .ne. nsamp2)) then
          write (*,*) 'Mismatch in number of traces or samples!'
          stop
        end if
        
        call getarg(3,cur_arg)
        read (cur_arg,*) mtype
        call getarg(4,cur_arg)
        read (cur_arg,*) m1
        call getarg(5,cur_arg)
        read (cur_arg,*) m2
        call getarg(6,cur_arg)
        read (cur_arg,*) m3
        
        if (nargs .ge. 7) then
          call getarg(7,weight_name)
          call readweight(weight_name,trweight,ntr1)
        else
          do itr=1,ntr1
            trweight(itr) = 1.
          end do
        end if
        
        write (*,*) m1,m2,m3
        
        write (*,*) 'Misfit type ',mtype,' = ',
     &    getmisfit(data_1,data_2,trweight,ntr1,nsamp1,m1,m2,m3,mtype)
        
      end
        
        
