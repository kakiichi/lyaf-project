program lya_visibility
  use param
  implicit none

  character(len=128) :: infile, outfile

   call read_arg

  open(10,file='param.input')
  read(10) r0
  read(10) slope
  read(10) local_phr
  read(10) req
  read(10) betaN

                                          
contains                                           
  subroutine read_arg
    implicit none
    integer            :: i,n 
    integer            :: iargc                          
    character(len=4)   :: opt                       
    character(len=128) :: arg                       

    n = iargc()
    if (n < 3) then
       print *, 'usage: RSD  -in  param.in'
       print *, '            -out RSD.dat'
       print *, 'e.g: skewers -in param.in -out RSD.dat'
       print *, ' '
       stop
    end if
    
    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case('-in')
          infile=trim(arg)
       case ('-out')
          outfile=trim(arg)
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
   
    return    
  end subroutine read_arg

end program lya_visibility
