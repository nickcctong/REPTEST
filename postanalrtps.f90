PROGRAM POSTANALRTPS

!  Module usage
   USE netcdf
   IMPLICIT NONE
!  Directory location of input data
   CHARACTER*10  :: varname
   CHARACTER*100 :: bkgdir, anldir
   CHARACTER*100 :: dynvfile,tracfile
   CHARACTER*100 :: prior_inp, postr_inp
   CHARACTER*100 :: gesfilenm, anlfilenm

   REAL, ALLOCATABLE :: VAR_PRIOR(:,:,:,:)
   REAL, ALLOCATABLE :: VAR2D_PRIOR(:,:,:)
   REAL, ALLOCATABLE :: VAR_POSTR(:,:,:,:)
   REAL, ALLOCATABLE :: VAR2D_POSTR(:,:,:)
   REAL, ALLOCATABLE :: ENS_PRIOR(:,:,:,:)
   REAL, ALLOCATABLE :: ENS_POSTR(:,:,:,:)
   REAL, ALLOCATABLE :: PTB_PRIOR(:,:,:,:)
   REAL, ALLOCATABLE :: PTB_POSTR(:,:,:,:)
   REAL, ALLOCATABLE :: INFPTB_POSTR(:,:,:,:)
   REAL, ALLOCATABLE :: INFENS_POSTR(:,:,:,:)
   INTEGER, ALLOCATABLE :: ncidm(:)
   REAL :: fsprd, asprd, tmp_chunk2, clip

   CHARACTER (len=3) :: ce
   INTEGER :: err_stat, ncid, varid 
   INTEGER :: NX, NY, NZ, NT
   INTEGER :: xdim, ydim, zdim
   integer :: istart(1:4),iend(1:4)   
   INTEGER :: nn,mm,ii,jj,kk  
   INTEGER :: ensmem, nvar
   REAL :: analpertwt
   CHARACTER (len=10) :: vars(1:20)

!  Namelist file parameters
   namelist/iosetup/ ensmem,bkgdir,anldir,dynvfile,tracfile
   namelist/rtpscoeff/ analpertwt
   namelist/varinfo/ nvar, vars

!---------------------------------------------------------------------------------------------
   write(6,'(/a)')' [1] Initialize information.'
!---------------------------------------------------------------------------------------------
   ensmem = 30
   bkgdir = './'
   anldir = './'
   dynvfile = 'fv_core.res.tile1.nc'
   tracfile = 'fv_tracer.res.tile1.nc'
   analpertwt = 0.95 
   nvar = 1
   vars = 'ua'

   open(unit=666, file='namelist.postrtps', &
        form='formatted', status='old', action='read')
   read(666, iosetup)
   read(666, rtpscoeff)
   read(666, varinfo)   
   close(666)

   write(6,'(2a)')'   Reading ens.Background from ',trim(bkgdir)
   write(6,'(2a)')'       and ens.Analysis from ',trim(anldir)
   write(6,'(a,f10.2)')'   Performing RTPS inflation with analpertwt: ',analpertwt
   write(6,'(a,i4)')'  Number of ensemble members = ',ensmem
   write(6,'(a,i4)')'  Number of variables to inflate = ',nvar  
   write(6,'(20a)')'   List of variables = ', vars(1:nvar)

!  PROGRAM BEGINS HERE -------------------------------
   CALL griddims_fv3(trim(bkgdir)//'/001/'//trim(tracfile),NX,NY,NZ,NT)

   DO nn = 1, nvar
     SELECT CASE (vars(nn))
       CASE ('delp','ua','va','W','T')
         gesfilenm=dynvfile
         anlfilenm=dynvfile
         varname=vars(nn)
         xdim=NX
         ydim=NY
         zdim=NZ
       CASE ('sphum','liq_wat','rainwat','ice_wat','snowwat','graupel')
         gesfilenm=tracfile
         anlfilenm=tracfile
         varname=vars(nn)
         xdim=NX
         ydim=NY
         zdim=NZ
     END SELECT
     istart=1
     iend(1)=xdim
     iend(2)=ydim
     iend(3)=zdim
     iend(4)=NT

     WRITE(*,*)'Processing variable:',varname

!  Loop thru ens firstguess and ens analysis for variables
     ALLOCATE(ncidm(ensmem))    
     ALLOCATE(ENS_PRIOR(xdim,ydim,zdim,ensmem),ENS_POSTR(xdim,ydim,zdim,ensmem))
     ALLOCATE(PTB_PRIOR(xdim,ydim,zdim,ensmem),PTB_POSTR(xdim,ydim,zdim,ensmem))
     ALLOCATE(INFPTB_POSTR(xdim,ydim,zdim,ensmem),INFENS_POSTR(xdim,ydim,zdim,1))
     ALLOCATE(VAR_PRIOR(xdim,ydim,zdim,NT))
     ALLOCATE(VAR_POSTR(xdim,ydim,zdim,NT))
     DO mm = 1, ensmem
        WRITE(UNIT=ce,FMT='(i3.3)')mm
        prior_inp = trim(bkgdir)//'/'//trim(ce)//'/'//gesfilenm
        err_stat  = NF90_OPEN(prior_inp,NF90_NOWRITE,ncid)
        err_stat  = NF90_INQ_VARID(ncid,varname,varid)
        err_stat  = NF90_GET_VAR(ncid,varid,VAR_PRIOR)
        ENS_PRIOR(:,:,:,mm)=VAR_PRIOR(:,:,:,1)
        err_stat = NF90_CLOSE(ncid)

        postr_inp = trim(anldir)//'/'//trim(ce)//'/'//anlfilenm
        err_stat  = NF90_OPEN(postr_inp,NF90_WRITE,ncid)
        ncidm(mm)=ncid
        err_stat  = NF90_INQ_VARID(ncid,varname,varid)
        err_stat  = NF90_GET_VAR(ncid,varid,VAR_POSTR)
        ENS_POSTR(:,:,:,mm)=VAR_POSTR(:,:,:,1) 
     END DO   
!  Calculate ens. prior and posterior perturbations
     DO ii = 1, xdim
       DO jj = 1, ydim
         DO kk = 1, zdim
           PTB_PRIOR(ii,jj,kk,:)=ENS_PRIOR(ii,jj,kk,:)-sum(ENS_PRIOR(ii,jj,kk,:))/ensmem
           fsprd = sum(PTB_PRIOR(ii,jj,kk,:)**2)/39
           fsprd = max(fsprd,tiny(fsprd))
           fsprd = sqrt(fsprd)
           fsprd = max(fsprd,tiny(fsprd))
           PTB_POSTR(ii,jj,kk,:)=ENS_POSTR(ii,jj,kk,:)-sum(ENS_POSTR(ii,jj,kk,:))/ensmem
           asprd = sum(PTB_POSTR(ii,jj,kk,:)**2)/39
           asprd = max(asprd,tiny(asprd)) 
           asprd = sqrt(asprd)
           asprd = max(asprd,tiny(asprd))
!  Apply inflation to posterior perturbations
           tmp_chunk2 = analpertwt*((fsprd-asprd)/asprd) + 1.0
           INFPTB_POSTR(ii,jj,kk,:) = tmp_chunk2*PTB_POSTR(ii,jj,kk,:)
!  Add inflated perturnations back to posterior mean for ens. members
           ENS_POSTR(ii,jj,kk,:) = sum(ENS_POSTR(ii,jj,kk,:))/ensmem + INFPTB_POSTR(ii,jj,kk,:) 
           if(isnan(ENS_POSTR(ii,jj,kk,1))) then
             WRITE(*,*)'@ i,j,k=',ii,jj,kk
             WRITE(*,*)'  fsprd=',fsprd,'asprd=',asprd,'tmp_chunk2=',tmp_chunk2
             WRITE(*,*)'INFPTB_POSTR=',INFPTB_POSTR(ii,jj,kk,:)
             WRITE(*,*)'   PTB_POSTR=',PTB_POSTR(ii,jj,kk,:)
           end if
         END DO
       END DO
     END DO
!  Clip Q variables to no smaller than 0 before write-out
     IF (anlfilenm .EQ. tracfile) THEN
       clip = tiny(ENS_POSTR(1,1,1,1))
       where (ENS_POSTR < clip) ENS_POSTR = clip
     END IF
!  Write out the inflated ens. analysis back to corresponding member files  
     DO mm = 1, ensmem
        INFENS_POSTR(:,:,:,1)=ENS_POSTR(:,:,:,mm)
!        err_stat  = NF90_INQ_VARID(ncidm(mm),varname,varid)
        call ncvpt( ncidm(mm), varid, istart, iend, INFENS_POSTR, err_stat) 
        err_stat = NF90_CLOSE(ncidm(mm))
     END DO
     DEALLOCATE(ENS_PRIOR,ENS_POSTR,PTB_PRIOR,PTB_POSTR)
     DEALLOCATE(INFPTB_POSTR,INFENS_POSTR)
     DEALLOCATE(VAR_PRIOR,VAR_POSTR)
     DEALLOCATE(ncidm)

   END DO
      

END PROGRAM POSTANALRTPS

!Get dimensions of a fv3_tracer data
  SUBROUTINE griddims_fv3(FV3_INPUT,NX,NY,NZ,NT)
  use netcdf
  INTEGER(KIND=4) :: ncid, status, DIM_ID
  INTEGER, INTENT(OUT) :: NX,NY,NZ,NT
  CHARACTER(*),INTENT(IN) :: FV3_INPUT

  status = nf90_open(FV3_INPUT, nf90_nowrite, ncid)
  status = nf90_inq_dimid(ncid, "xaxis_1", DIM_ID)
  status = nf90_inquire_dimension(ncid, DIM_ID, len = NX)
  status = nf90_inq_dimid(ncid, "yaxis_1", DIM_ID)
  status = nf90_inquire_dimension(ncid, DIM_ID, len = NY)
  status = nf90_inq_dimid(ncid, "zaxis_1", DIM_ID)
  status = nf90_inquire_dimension(ncid, DIM_ID, len = NZ)
  status = nf90_inq_dimid(ncid, "Time", DIM_ID)
  status = nf90_inquire_dimension(ncid, DIM_ID, len = NT)
  status = nf90_close(ncid)

  END SUBROUTINE griddims_fv3
