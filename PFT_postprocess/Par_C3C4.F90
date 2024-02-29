PROGRAM pftlai_c3c4

!=======================================================
! USAGE:
!
! - to compile:
! gfortran -g -mcmodel=large -fbounds-check -o mkmod mkmod.F90 -I/usr/include -lnetcdf -lnetcdff
!
! - to run (3 ways):
! ./mkmod                !default region file reg_5x5 and year 2005
! ./mkmod reg_5x5 2005   !input region file name and year
! ./mkmod.sh 2005        !bash script in paralell mode (details see mkmod.sh)
!
! ________________
! History:
!   2019/06: Hua Yuan, Initial R code version
!   2022/02: Wenzong Dong, Rewrite R code to Fortran version
!   2023/07: Jiahao Shi, Make code cuts to fit individualized needs. 
!
!=======================================================================
    
     USE netcdf
    
    
     IMPLICIT NONE
    
     INTEGER , parameter :: r8 = selected_real_kind(12)
     INTEGER , parameter :: pfts = 16, mon = 12
    
     CHARACTER(len=*), parameter :: Title   = "Land surface model input vagetation data"
     CHARACTER(len=*), parameter :: Authors = "Dai YJ group at Sun Yat-sen University"
     CHARACTER(len=*), parameter :: Address = "School of Atmospheric Sciences, Sun Yat-sen University, Zhuhai, China"
     CHARACTER(len=*), parameter :: Email   = "yuanh25@mail.sysu.edu.cn"
   
     INTEGER , dimension(12)    :: mons, dom
     INTEGER , dimension(16)    :: pftnum
     REAL(r8), dimension(12)    :: phimin  !=0.5
     REAL(r8), dimension(16)    :: laimax, sgdd, &
                                   minfr, saimin, sairtn, saimin1, saires, &
                                   pctpft, laiinimax
     REAL(r8), dimension(16,12) :: phi, laiini, laiini1, saiini2
    
     ! output data
     REAL(r8) :: lat, lon
     REAL(r8), dimension(1,1,12)    :: lclai, lcsai
     REAL(r8), dimension(1,1,16)    :: ppft
     REAL(r8), dimension(1,1,16,12) :: pftlai, pftsai
     REAL(r8), dimension(1,1,12) :: laitot_out
     CHARACTER (len=255) :: filename
    
     ! input vars id
     INTEGER :: ncid
     ! output vars id
     INTEGER :: csai_id, lat_dimid, lat_vid, lon_dimid, lon_vid, mon_dimid, &
                mon_vid, pft_dimid, pft_vid, plai_id, ppft_id, psai_id, mlai_id
    
     ! vars
     INTEGER         :: loc1(1), loc2(1)
     INTEGER         :: i, j, imonth, itmin, nextmonth, prevmonth, &
                        imonth_prev, laimaxloc, iloop, inx, ipft, mcnt
     INTEGER , dimension(5)     :: indx1=(/2,3,5,6,10/)
     INTEGER , dimension(10)    :: indx2=(/4,7,8,9,11,12,13,14,15,16/)
     REAL                       :: fillvalue
     REAL(r8), dimension(46)    :: lai
     REAL(r8), dimension(12)    :: laitot, dd2, dd5, gdd2, gdd5, &
                                   rgdd2, rgdd5, tmax, tmin, tavg
     REAL(r8), dimension(16,12) :: saiini, saiini1, laidiff
    
     ! short data
     REAL(r8), dimension(12) :: prec
     REAL(r8) :: laitot_nonevg, laiup
     REAL(r8) :: herb_cover, summ, davg, sumwgt, &
                 sumevg, frac_c4, sumnon, x1, x2, sum_judg
     INTEGER  :: XY2D(2), XY3D(3), XY4D(4), XY3F(3)    
    
     DO i=1,12,1
        mons(i) = i
     ENDDO
    
     DO i=1,16,1
        pftnum(i) = i-1
     ENDDO
    
     dom = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    

    
     ! values in the table below are from Lawrence et al., 2007
     ! with modifications according to Sitch et al., 2003
     laimax = (/0, 5, 5, 5, 7, 7, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4/)
     sgdd   = (/0, 0, 0, 100, 0, 0, 0, 200, 200, 0, 100, 100, 100, 100, 100, 100/)
     minfr  = (/0., 0.7, 0.7, 0., 0.8, 0.8, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./)
     saimin = (/0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.1/)
     sairtn = (/0., 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, &
                0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0./)
     phimin(:)   = 0.5
     phi   (:,:) = 1.
     laiini(:,:) = 0.

    
    !-----------Input：Proportion of vegetation cover, mean monthly precipitation, air temperature, mean monthly maximum and minimum air temperature, mean monthly LAI----------------

    ppft(1,1,:) = (/40.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,60.0,0.0,0.0/) 
    pctpft(:) = (/40.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,60.0,0.0,0.0/)/100. 

    prec(:) = (/14.3,10.2,9.8,2.3,4.2,7.0,87.4,80.6,46.5,8.6,9.8,12.0/)                    
    tavg(:) = (/8.9,9.8,13.0,16.5,21.0,26.1,24.7,23.9,22.6,18.5,13.7,8.8/)              
    tmax(:) = (/15.8,17.0,21.3,24.2,30.3,32.8,32.8,29.4,29.8,27.1,23.8,17.1/) 
    tmin(:) = (/3.7,3.4,5.2,9.5,12.7,19.0,19.5,19.3,17.4,11.9,5.9,2.1/) 
        
    laitot = (/0.2,0.2,0.2,0.2,0.3,0.2,0.3,0.6,0.6,0.2,0.2,0.1/)                        

    
     ! calculate GDD and phi
     phi(:,:) = 1.
    
     IF (any(tavg<-3.4e+2)) THEN  ! 修改缺省值bug (all(tavg==-3.4e+38)-->any(tavg==-3.4e+38))
        ! when T and P are missing
        PRINT*, 'NA temperature! Error!'
     ELSE
        loc1    = minloc(tavg) ! check difference between min and which.min(R)
        itmin   = loc1(1)
        mcnt    = 0
        dd2  (:)= 0.
        dd5  (:)= 0.
        gdd2 (:)= 0.
        gdd5 (:)= 0.
        rgdd2(:)= 0.
        rgdd5(:)= 0.
     ENDIF
    
     DO WHILE (mcnt < 12)
    
        imonth = mod((itmin+mcnt-1),12) + 1
    
        ! calculate gdd according to tmin, tmax, tavg, 2, 5degree
        ! ------------------------------
    
        ! calcualte day index of tavg
        IF ((tmax(imonth)-tmin(imonth)) > 0.) THEN
           davg = (tmax(imonth)-tavg(imonth))/ &
                 (tmax(imonth)-tmin(imonth))*dom(imonth)    !dom为常数
        ELSE
           davg = 0
        ENDIF
    
        ! gdd2 in different cases
        IF (2. < tmin(imonth)) THEN
           dd2(imonth) = (tavg(imonth)-2.)* dom(imonth)
        ELSE IF (2. >= tmax(imonth)) THEN
           dd2(imonth) = 0.
        ELSE IF (2.>=tmin(imonth) .and. 2.<tavg(imonth) .and. tavg(imonth)/=tmin(imonth)) THEN
           dd2(imonth) = (tavg(imonth)-2.)/2.*(tavg(imonth)-2.)/ &
                       (tavg(imonth)-tmin(imonth))*davg + &
                       (tmax(imonth)+tavg(imonth)-4.)/2.* &
                       (dom(imonth)-davg)
        ELSE IF (2.>=tavg(imonth) .and. 2.<tmax(imonth) &
              .and. tmax(imonth)/=tavg(imonth)) THEN
           dd2(imonth) = (tmax(imonth)-2.)/2.*(tmax(imonth)-2.)/ &
                       (tmax(imonth)-tavg(imonth))*(dom(imonth)-davg)
        ELSE
           dd2(imonth) = max(0., ((tavg(imonth)-2.)*dom(imonth)))
        ENDIF
    
        ! gdd5 in different cases
        IF (5. < tmin(imonth)) THEN
           dd5(imonth) = (tavg(imonth)-5.)* dom(imonth)
        ELSE IF (5. >= tmax(imonth)) THEN
           dd5(imonth) = 0.
        ELSE IF (5.>=tmin(imonth) .and. 5.<tavg(imonth) &
              .and. tavg(imonth)/=tmin(imonth)) THEN
           dd5(imonth) = (tavg(imonth)-5.)/2.*(tavg(imonth)-5.)/ &
                       (tavg(imonth)-tmin(imonth))*davg + &
                       (tmax(imonth)+tavg(imonth)-10.)/2.* &
                       (dom(imonth)-davg)
        ELSE IF (5.>=tavg(imonth) .and. 5.<tmax(imonth) &
              .and. tmax(imonth)/=tavg(imonth)) THEN
           dd5(imonth) = (tmax(imonth)-5.)/2.*(tmax(imonth)-5.)/ &
                       (tmax(imonth)-tavg(imonth))*(dom(imonth)-davg)
        ELSE
           dd5(imonth) = max(0., ((tavg(imonth)-5.)*dom(imonth)))
        ENDIF
    
        mcnt = mcnt + 1
     ENDDO
    
     ! split the tmin month gdd into 2 parts
     ! gdd, rgdd
    
     IF (dd2(itmin) < 1e-6) THEN
        gdd2 (itmin) = 0.
        rgdd2(itmin) = 0.
     ELSE
        nextmonth = mod(itmin,12) + 1
        prevmonth = mod((10+itmin),12) + 1
        summ      = max(0., (tmin(nextmonth)-2.)) + &
                    max(0., (tmin(prevmonth)-2.))
    
        IF (summ > 0.) THEN
           gdd2 (itmin) = dd2(itmin)*max(0., (tmin(nextmonth)-2.))/summ
           rgdd2(itmin) = dd2(itmin)*max(0., (tmin(prevmonth)-2.))/summ
        ELSE
           IF ((tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin)) > 0) THEN
              gdd2 (itmin) = dd2(itmin)*(tavg(nextmonth)-tavg(itmin))/ &
                             (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
              rgdd2(itmin) = dd2(itmin)*(tavg(prevmonth)-tavg(itmin))/ &
                             (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
           ELSE
              gdd2 (itmin) = dd2(itmin) * 0.5
              rgdd2(itmin) = dd2(itmin) * 0.5
           ENDIF
        ENDIF
     ENDIF
    
     IF (dd5(itmin) < 1e-6) THEN
        gdd5 (itmin) = 0.
        rgdd5(itmin) = 0.
     ELSE
        nextmonth = mod(itmin,12) + 1
        prevmonth = mod((itmin+10),12) + 1
        summ      = max(0., (tmin(nextmonth)-5.)) + &
                    max(0., (tmin(prevmonth)-5.))
    
        IF (summ > 0.) THEN
           gdd5 (itmin) = dd5(itmin)*max(0., tmin(nextmonth)-5.)/summ
           rgdd5(itmin) = dd5(itmin)*max(0., tmin(prevmonth)-5.)/summ
        ELSE
           IF ((tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin)) > 0) THEN
              gdd5 (itmin) = dd5(itmin)*(tavg(nextmonth)-tavg(itmin))/ &
                             (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
              rgdd5(itmin) = dd5(itmin)*(tavg(prevmonth)-tavg(itmin))/ &
                             (tavg(nextmonth)+tavg(prevmonth)-2*tavg(itmin))
           ELSE
              gdd5 (itmin) = dd5(itmin) * 0.5
              rgdd5(itmin) = dd5(itmin) * 0.5
           ENDIF
        ENDIF
     ENDIF
    
     imonth_prev = itmin
     mcnt        = 1
     DO WHILE (mcnt < 12)
    
        imonth = mod((itmin+mcnt-1),12) + 1
    
        gdd2(imonth) = gdd2(imonth_prev) + dd2(imonth)
        gdd5(imonth) = gdd5(imonth_prev) + dd5(imonth)
    
        imonth_prev  = imonth
        mcnt         = mcnt + 1
     ENDDO
    
     imonth_prev = itmin
     mcnt        = 1
     DO WHILE (mcnt < 12)
    
        imonth = mod((12+itmin-mcnt-1),12) + 1
    
        rgdd2(imonth) = rgdd2(imonth_prev) + dd2(imonth)
        rgdd5(imonth) = rgdd5(imonth_prev) + dd5(imonth)
    
        imonth_prev = imonth
        mcnt        = mcnt + 1
     ENDDO
    
     gdd2 (itmin) = dd2(itmin)
     rgdd2(itmin) = dd2(itmin)
     gdd5 (itmin) = dd5(itmin)
     rgdd5(itmin) = dd5(itmin)
    
     ! calculate phi using gdd/sgdd
     !!!!!!
     DO imonth=1,12,1
        gdd2(imonth)   = min(gdd2(imonth), rgdd2(imonth))
        gdd5(imonth)   = min(gdd5(imonth), rgdd5(imonth))
    
        phi (4,imonth) = max(phimin(imonth), &
                             min(1., gdd2(imonth)/sgdd(4)))
        phi (8,imonth) = max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(8)))
        phi (9,imonth) = max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(9)))
        phi (11,imonth)= max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(11)))
        phi (12,imonth)= max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(12)))
        phi (13,imonth)= max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(13)))
        phi (14,imonth)= max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(14)))
        phi (15,imonth)= max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(15)))
        phi (16,imonth)= max(phimin(imonth), &
                             min(1., gdd5(imonth)/sgdd(16)))
     ENDDO
    
     ! distribution LAI
     ! sum(lai*PFT%) = MONTHLY_LAI
    
     ! calculate initial LAI
     laiini(:,:) = 0.
     
    
     DO imonth=1,12,1
        sumwgt = sum(phi(:,imonth)*laimax(:)*pctpft(:))
    
        IF (sumwgt > 0.) THEN
           laiini(:,imonth) = phi(:,imonth)*laimax(:)/sumwgt*laitot(imonth)
        ELSE
           laiini(:,imonth) = 0.
        ENDIF
     ENDDO
    
     ! adjust LAI
     ! -----------------------------
     !laiinimax = apply(laiini, 1, max) !check
     loc2      = maxloc(laitot)
     laimaxloc = loc2(1)
     laiinimax(:) = laiini(:,laimaxloc)
    
     DO imonth=1,12,1
    
        ! evergreen tree phenology (Zeng et al., 2002)
        ! minimum fraction of maximum initial PFT LAI
        ! laiini[,imonth] = pmax(laiini[,imonth], laiinimax*minfr)
    
        ! max value limited
        DO iloop=1,16
           IF (laiini(iloop,imonth) > laiinimax(iloop)) THEN
              laiini(iloop,imonth) = laiinimax(iloop)
           ENDIF
        ENDDO
    
        ! calculate for nonevergreen PFT
        !indx1  = (/2,3,5,6,10/)
        sumevg = sum(laiini(indx1(:),imonth)*pctpft(indx1(:)))
        laitot_nonevg = max(laitot(imonth)-sumevg, 0.)
    
        !indx2  = (/4,7,8,9,11,12,13,14,15,16/)
    
        sumnon = sum(phi(indx2(:),imonth)*laimax(indx2(:))*pctpft(indx2(:)))
        IF (sumnon > 0.) THEN
           DO inx=1,10
              laiini(indx2(inx),imonth) = phi(indx2(inx),imonth) * &
                                  laimax(indx2(inx)) / sumnon * laitot_nonevg
           ENDDO
        ELSE
           sumnon = sum(laimax(indx2(:))*pctpft(indx2(:)))
           IF (sumnon > 0.) THEN
              ! do not consider phenology
              DO inx=1,10
                 laiini(indx2(inx),imonth) = laimax(indx2(inx)) / sumnon * laitot_nonevg
              ENDDO
           ELSE
              ! no percentage cover
              DO inx=1,10
                 laiini(indx2(inx),imonth) = 0._r8
              ENDDO
           ENDIF
        ENDIF
     ENDDO
    
     ! max LAI value limited (need more test and check)
     IF (count(laiini>10.)>0) THEN
        WHERE (laiini>10.) laiini=10.
     ENDIF
    
     ! split C3/C4 nature grass (Still et al., 2003)
     ! ------------------------------
    
     herb_cover = 0
     herb_cover = sum(ppft(1,1,13:16))
    ! print *, '*******************************'
     print *, herb_cover
     ! C4 case
     IF (count(tavg>22.)==12 .and. count(prec>25.)==12) THEN
        ppft(1,1,14+1)  = sum(ppft(1,1,13:14)) ! index check
        ppft(1,1,13:14) = 0.                   ! index check
     ELSE IF (herb_cover>75 .and. count(laiini(14,:)>0)>0) THEN
        ! minxed C3/C4 case
        laiup = 0.                                            
        DO inx=1,12
           IF (tavg(inx)>22. .and. prec(inx)>25.) THEN
              laiup = laiup + laiini(13+1,inx)
           ENDIF
        ENDDO
        frac_c4 = laiup/sum(laiini(14,:))
    
        ppft(1,1,14+1) = ppft(1,1,13+1) * frac_c4
        ppft(1,1,13+1) = ppft(1,1,13+1) * (1.-frac_c4)
     ENDIF
    
    
    
    
     
     ! calculate SAI
     ! --------------------
     DO ipft=1,pfts,1
        saimin1(ipft) = saimin(ipft)*maxval(laiini(ipft,:))
     ENDDO
     saimin1(2:16)= saimin1(2:16)/laimax(2:16)
     saiini (:,:) = 0.
     saiini1(:,:) = saiini(:,:)
    
     DO iloop=1,12
        saiini1(:,iloop) = saimin(:)
     ENDDO
    
     ! PFT LAI diff
     laidiff(:,1) = laiini(:,12) - laiini(:,1)
     DO iloop=1,11,1
        inx = iloop+1
        laidiff(:,inx) = laiini(:,iloop) - laiini(:,inx)
     ENDDO
    
     WHERE(laidiff<0.) laidiff=0.
    
     iloop = 1
     DO WHILE (iloop<13)
        saires(:) = sairtn*saiini1(:,mod((iloop+10),12)+1)
        DO inx=1,16
           x1 = (saires(inx)+laidiff(inx,iloop)*0.5)
           x2 = saimin1(inx)
           saiini(inx,iloop) = max(x1,x2) !maxval((saires(:)+laidiff(:,imonth)*0.5), saimin1())
        ENDDO
        iloop = iloop + 1
    
        IF (iloop == 13) THEN
           sum_judg = abs(sum(saiini-saiini1))
           IF (sum_judg > 1e-6) THEN
              iloop  = 1
              saiini1(:,:) = saiini(:,:)
           ENDIF
        ENDIF
     ENDDO
    

     IF (count(laiini>20) > 0) THEN
        PRINT*, 'check for laiini!'
        STOP
     ENDIF
    
     ! max SAI value limited
     IF (count(saiini>3) > 0) THEN
        WHERE(saiini>3.) saiini = 3.
     ENDIF
    
    
     ! output data
     pftlai(1,1,:,:) = laiini(:,:)
     pftsai(1,1,:,:) = saiini(:,:)
     laitot_out(1,1,:) = laitot(:)
    
    filename = '/stu01/shijh21/data/parlai/US-Wkg_pftlai_c3c4.nc' 
    lat = 31.7364 
    lon = -109.9419 
    
   DO iloop=1,12
      laiini1(:,iloop) = pctpft(:)*laiini(:,iloop)
      saiini2(:,iloop) = pctpft(:)*saiini(:,iloop)
   ENDDO

   DO iloop=1,12
      lclai(:,:,iloop) = sum(laiini1(:,iloop))
      lcsai(:,:,iloop) = sum(saiini2(:,iloop))
   ENDDO
      

    CALL check( nf90_create(filename, NF90_NETCDF4, ncid) )
        ! define dimensions
    CALL check( nf90_def_dim(ncid, "lat", 1, lat_dimid) )
    CALL check( nf90_def_dim(ncid, "lon", 1, lon_dimid) )
    CALL check( nf90_def_dim(ncid, "pft", pfts , pft_dimid) )
    CALL check( nf90_def_dim(ncid, "mon", mon  , mon_dimid) )
    
    
        ! define variables
        ! ---------------------------------------
    CALL check( nf90_def_var(ncid, "lat", NF90_FLOAT, lat_dimid, lat_vid, deflate_level=6) )
    CALL check( nf90_def_var(ncid, "lon", NF90_FLOAT, lon_dimid, lon_vid, deflate_level=6) )
    CALL check( nf90_def_var(ncid, "pft", NF90_INT  , pft_dimid, pft_vid, deflate_level=6) )
    CALL check( nf90_def_var(ncid, "mon", NF90_INT  , mon_dimid, mon_vid, deflate_level=6) )
    
    CALL check( nf90_put_att(ncid, lat_vid, "long_name", "Latitude"     ) )
    CALL check( nf90_put_att(ncid, lat_vid, "units"    , "degrees_north") )
    CALL check( nf90_put_att(ncid, lon_vid, "long_name", "Longitude"    ) )
    CALL check( nf90_put_att(ncid, lon_vid, "units"    , "degrees_east" ) )
    CALL check( nf90_put_att(ncid, pft_vid, "long_name", "Index of PFT" ) )
    CALL check( nf90_put_att(ncid, mon_vid, "long_name", "Month of year") )
    CALL checK( nf90_put_att(ncid, mon_vid, "units"    , "month"        ) )
    
    
    CALL check( nf90_put_var(ncid, lat_vid, lat  ) )
    CALL check( nf90_put_var(ncid, lon_vid, lon  ) )
    CALL check( nf90_put_var(ncid, pft_vid, pftnum) )
    CALL check( nf90_put_var(ncid, mon_vid, mons  ) )
    
    
    XY2D = (/lon_dimid, lat_dimid/)
    XY3D = (/lon_dimid, lat_dimid, mon_dimid/)
    XY3F = (/lon_dimid, lat_dimid, pft_dimid/)
    XY4D = (/lon_dimid, lat_dimid, pft_dimid, mon_dimid/)
    fillvalue = -999.
    CALL check( nf90_def_var(ncid, "MONTHLY_PFT_LAI", NF90_FLOAT  , XY4D, plai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, plai_id          , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, plai_id          , "long_name" , "Monthly PFT LAI values") )
    CALL check( nf90_put_att(ncid, plai_id          , "_FillValue", fillvalue    ) )
    
    CALL check( nf90_def_var(ncid, "MONTHLY_PFT_SAI", NF90_FLOAT  , XY4D, psai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, psai_id          , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, psai_id          , "long_name" , "Monthly PFT SAI values") )
    CALL check( nf90_put_att(ncid, psai_id          , "_FillValue", fillvalue    ) )
    
    CALL check( nf90_def_var(ncid, "PCT_PFT"        , NF90_FLOAT  , XY3F, ppft_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, ppft_id          , "units"     , "%"          ) )
    CALL check( nf90_put_att(ncid, ppft_id          , "long_name" , "Percent PFT cover") )
    CALL checK( nf90_put_att(ncid, ppft_id          , "_FillValue", fillvalue    ) )
    
    CALL check( nf90_def_var(ncid, "MONTHLY_LAI", NF90_FLOAT  , XY3D, mlai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, mlai_id          , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, mlai_id          , "long_name" , "Monthly LAI values") )
    CALL check( nf90_put_att(ncid, mlai_id          , "_FillValue", fillvalue    ) )
    

    CALL check( nf90_def_var(ncid, "MONTHLY_LC_SAI" , NF90_FLOAT  , XY3D, csai_id, deflate_level=6) )
    CALL check( nf90_put_att(ncid, csai_id          , "units"     , "m^2/m^2"    ) )
    CALL check( nf90_put_att(ncid, csai_id          , "long_name" , "Monthly landcover SAI values") )
    CALL check( nf90_put_att(ncid, csai_id          , "_FillValue", fillvalue    ) )

    
    CALL check( nf90_inq_varid(ncid, "MONTHLY_PFT_LAI", plai_id ) )
    CALL check( nf90_put_var  (ncid, plai_id          , pftlai  ) )
    CALL check( nf90_inq_varid(ncid, "MONTHLY_PFT_SAI", psai_id ) )
    CALL check( nf90_put_var  (ncid, psai_id          , pftsai  ) )
    CALL check( nf90_inq_varid(ncid, "PCT_PFT"        , ppft_id ) )
    CALL check( nf90_put_var  (ncid, ppft_id          , ppft    ) )
    CALL check( nf90_inq_varid(ncid, "MONTHLY_LAI"    , mlai_id ) )
    CALL check( nf90_put_var  (ncid, mlai_id          , laitot_out ) )

    CALL check( nf90_inq_varid(ncid, "MONTHLY_LC_SAI" , csai_id ) )
    CALL check( nf90_put_var  (ncid, csai_id          , lcsai   ) )
    
    CALL check( nf90_close(ncid) )
    print *, laitot
    CONTAINS
    
       SUBROUTINE check(status)
          INTEGER, intent(in) :: status
    
          IF (status /= nf90_noerr) THEN
             print *, trim( nf90_strerror(status))
             stop 2
          ENDIF
       END SUBROUTINE check
    
    END PROGRAM pftlai_c3c4
    
    
    