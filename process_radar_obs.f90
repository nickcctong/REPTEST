program read
  implicit none

  character(len=500) :: obsfile
  character(len=4) :: pe_name
  character(len=3) :: obtype
  integer  :: iunit, nchar, nreal, ii, mype, ios, idate, i, ipe, ioff0, k
  integer,dimension(2) :: nn,nobst, nobsps, nobsq, nobsuv, nobsgps, &
                          nobstcp,nobstcx,nobstcy,nobstcz,nobssst, nobsspd, nobsdw, nobsrw, nobspw, &
                          nobsdbz, nobssrw ! CAPS added nobssrw 
  character(8),allocatable,dimension(:):: cdiagbuf
  real(4),allocatable,dimension(:,:)::rdiagbuf
  real(4) :: errorlimit,errorlimit2,error,pres,obmax
  real(4) :: errorlimit2_obs,errorlimit2_bnd
  logical :: fexist, init_pass

  obsfile="diag_conv_ges.ensmean"
  open(88,form="unformatted",file=obsfile,iostat=ios)
  open(11,form="formatted",file="dbz_diag_conv_ges.ensmean",iostat=ios)
  open(22,form="formatted",file="rw_diag_conv_ges.ensmean",iostat=ios)
  read(88) idate
10 continue

  read(88,err=20,end=30) obtype,nchar,nreal,ii,mype,ioff0
  allocate(cdiagbuf(ii),rdiagbuf(nreal,ii))
  read(88) cdiagbuf(1:ii),rdiagbuf(:,1:ii)

  if (obtype .eq. "dbz") then
     do k=1,ii ! if you want to extract other data, please refer to 'contents_binary_CAPS_diag' subroutine in 'setupdbz.f90'
                                  ! lat            ! lon          ! hgt          ! obs            ! air_pressure(Pascal)
           write (11,'(5(f12.5,1x))') rdiagbuf(3,k), rdiagbuf(4,k), rdiagbuf(7,k), rdiagbuf(6,k)*100., rdiagbuf(17,k)
     end do
  end if
  
  if (obtype .eq. " rw") then
     do k=1,ii ! if you want to extract other data, please refer to 'contents_binary_CAPS_diag' subroutine in 'setupdbz.f90'
                                  ! lat            ! lon          ! hgt  
                                  ! ! obs            ! air_pressure(Pascal)
           write (22,'(5(f12.5,1x),2(f12.7,1x))') rdiagbuf(3,k), rdiagbuf(4,k), rdiagbuf(7,k), rdiagbuf(6,k)*100., & 
                                                 rdiagbuf(17,k), rdiagbuf(20,k), rdiagbuf(21,k)
     end do
  end if

  deallocate(cdiagbuf)
  deallocate(rdiagbuf)
  go to 10
20     continue
       print *,'error reading diag_conv file',obtype
30     continue
  
  close(11)
  close(22)
  close(88)
  print *,"finished "
end program read
