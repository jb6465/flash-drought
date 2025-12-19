      PROGRAM GENERATE_7DAY_TS
C     ************************
C     Creates time series of 7-day running mean drivers from daily data,
C     for use in attribution analysis of changes in ETrs into contributions
C     from its drivers. Converts daily maximum and minimum 2-m air
C     temperatures (Tmax and Tmin) to daily mean 2-m air temperature (T_2m),
C     and 10-m windspeed (U10m) to 2-m windspeed (U_2m). Does not estimate 
C     data for February 29 in non-leap years.

C     Input data, from NLDAS datafiles:
C       Varin = daily reference ET drivers

C     Output:
C       Varmean = 7-day running mean reference ET drivers

C     Author: Mike Hobbins
C     Date: May 21, 2024

      IMPLICIT  NONE
      INTEGER   mlen(12),nrows,ncols,bnd_s,bnd_n,bnd_w,bnd_e,bnd_lo,
     +          bnd_hi,yrstart,yrend,row,col,colin,rowin,year,mon,lpdy,
     +          day,var,winday,wday,wmon,wyear,i
      REAL      Lat,Lon,Elev
      CHARACTER varnamein(5)*4,varnameout(5)*4,root*99,maskfile*99,
     +          outfile*99,inpath*99,infile*99
      DATA mlen/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA varnamein/'ETrs','SWdn','SpHm','U10m','Tmax'/
      DATA varnameout/'ETrs','SWdn','SpHm','U_2m','T_2m'/

C     Define nrows and ncols, geographic bounds of analysis, and input
C     and output path.
      PARAMETER (nrows   = 224,
     +           ncols   = 464,
     +           bnd_s   = 25,
     +           bnd_n   = 53,
     +           bnd_w   = -125,
     +           bnd_e   = -67,
     +           bnd_lo  = -90,
     +           bnd_hi  = 4500,
     +           yrstart = 1980,
     +           yrend   = 2021,
     +           root    = '/Volumes/mhobbins_drobo/')

      INTEGER Mask(nrows,ncols)
      REAL    Varin(8,nrows,ncols),Tmin(ncols),Varmean(ncols)


C     Read in geographic data to set mask for analysis
      Mask = 0
      WRITE(maskfile,'(A35,A31)') 
     +      '/System/Volumes/Data/data/mhobbins/',
     +      'mac159/scripts/gtopomean15k.asc'
      OPEN(UNIT=10,FILE=maskfile,STATUS='OLD',ERR=99)
      DO 20 row = nrows, 1, -1
        DO 10 col = 1, ncols
C         - read input data: Col; Row; Lat [deg], latitude; Lon [deg],
C         longitude; Elev [m], elevation
          READ(10,*) colin,rowin,Lat,Lon,Elev
C         Generate geographic analysis mask
          IF((Lat.GT.bnd_s).AND.(Lat.LT.bnd_n).AND.
     +       (Lon.GT.bnd_w).AND.(Lon.LT.bnd_e).AND.
     +       (Elev.GT.bnd_lo).AND.(Elev.LT.bnd_hi))
     +       Mask(row,col) = 1
10      ENDDO                                    ! end col loop
20    ENDDO                                      ! end row loop
      CLOSE(10)


C     Loop through years
      DO 90 year = yrstart, yrend

C       Loop through months
        DO 80 mon = 1, 12

          lpdy = 0
          IF((mon.EQ.2).AND.(MOD(year,4).EQ.0)) lpdy = 1

C         Loop through days
          DO 70 day = 1, mlen(mon) + lpdy

            PRINT*,'Date:',year * 10000 + mon * 100 + day

C           Initialize variable array
            Varin = -9999.

C           Loop through input and output variables
            DO 60 var = 1, 5

C             Open daily output file for weekly mean variable
              WRITE(outfile,'(A24,A32,A4,A1,I4,2I2.2,A4)')
     +            root,'ETrc_attribution/CONUS/7_day_ts/',
     +            varnameout(var),'_',year,mon,day,'.asc'
              OPEN(UNIT=13,FILE=outfile,STATUS='UNKNOWN')
              WRITE(13,'(A9)') 'ncols 464'
              WRITE(13,'(A9)') 'nrows 224'
              WRITE(13,'(A14)') 'xllcorner -125'
              WRITE(13,'(A12)') 'yllcorner 25'
              WRITE(13,'(A14)') 'cellsize 0.125'
              WRITE(13,'(A19)') 'nodata_value -9999.'

C             Loop through days in the week window
              DO 30 winday = -3, 3

C               Determine date of start of week
                wday  = day + winday
                wmon  = mon
                wyear = year
                lpdy = 0
                IF((wmon.EQ.2).AND.(MOD(wyear,4).EQ.0)) lpdy = 1
                IF(wday.GT.mlen(mon)+lpdy) THEN
                  wmon = mon + 1
                  IF(wmon.GT.12) THEN
                    wmon = 1
                    wyear = year + 1
                  ENDIF
                  wday = wday - (mlen(mon) + lpdy)
                ELSEIF(wday.LT.1) THEN
                  wmon = mon - 1
                  lpdy = 0
                  IF((wmon.EQ.2).AND.(MOD(wyear,4).EQ.0)) lpdy = 1
                  IF(wmon.LE.0) THEN
                    wmon = 12
                    wyear = wyear - 1
                  ENDIF
                  wday = wday + mlen(wmon) + lpdy
                ENDIF

C               Open daily input file for variable
                inpath = 'ETrc_data/daily/'
                WRITE(infile,'(A24,A16,A4,A1,I4,2I2.2,A4)') root,
     +                  inpath,varnamein(var),'_',wyear,wmon,wday,'.asc'
                IF(var.GT.1) THEN
                  inpath = 'NLDAS_drivers/daily/'
                  WRITE(infile,'(A24,A20,A4,A1,I4,2I2.2,A4)') root,
     +                  inpath,varnamein(var),'_',wyear,wmon,wday,'.asc'
                ENDIF
                OPEN(UNIT=11,FILE=infile,STATUS='OLD',ERR=99)

C               Read time-series grids
                DO i = 1, 6                      ! read headers
                  READ(11,*)
                ENDDO
                DO row = 1, nrows                ! read data
                  READ(11,*) (Varin(winday+4,row,col),col=1,ncols)
                ENDDO
                CLOSE(11)
C               Read time-series grids for Tmin to combine with Tmax for T_2m
                IF(var.EQ.5) THEN
                  WRITE(infile,'(A24,A20,A5,I4,2I2.2,A4)')
     +                  root,inpath,'Tmin_',wyear,wmon,wday,'.asc'
                  OPEN(UNIT=12,FILE=infile,STATUS='OLD',ERR=99)
                  DO i = 1, 6                      ! read headers
                    READ(12,*)
                  ENDDO
                  DO row = 1, nrows                ! read data
                    READ(12,*) (Tmin(col),col=1,ncols)
                    DO col = 1, ncols
                      Varin(winday+4,row,col) = (Varin(winday+4,row,col)
     +                                           + Tmin(col)) / 2.
                    ENDDO
                  ENDDO
                  CLOSE(12)
                ENDIF
30            ENDDO                                ! end winday loop

C             Initialise output array
              Varmean = -9999.

C             Estimate 7-day means
              DO 50 row = 1,nrows

                DO 40 col = 1,ncols

                  IF(Mask(row,col).EQ.1) THEN
                    Varmean(col) = (Varin(1,row,col) + Varin(2,row,col)
     +                            + Varin(3,row,col) + Varin(4,row,col)
     +                            + Varin(5,row,col) + Varin(6,row,col)
     +                            + Varin(7,row,col)) / 7.
                    IF(var.EQ.4) Varmean(col) = Varmean(col)
     +                                  * 4.87 / log(67.8 * 10. - 5.42)
                  ELSE
                    Varmean(col) = -9999.
                  ENDIF

40              ENDDO                            ! end col loop

C               Write output data
                WRITE(13,*) (Varmean(col),col=1,ncols)

50            ENDDO                              ! end row loop

C             Close output file
              CLOSE(13)

60          ENDDO                                ! end var loop

70        ENDDO                                  ! end day loop

80      ENDDO                                    ! end mon loop

90    ENDDO                                      ! end year loop

99    END