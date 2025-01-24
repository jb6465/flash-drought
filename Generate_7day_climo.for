      PROGRAM GENERATE_7DAY_CLIMO
C     ***************************
C     Creates climatological, multi-year means of 7-day ETrs and 7-day drivers,
C     for use in attribution analysis of changes in ETrs into its drivers.

C     Assumes all Februarys have 29 days, using 7-day timeseries of February 29
C     in leap years and March 1 for non-leap years.

C     Input data, from NLDAS datafiles:
C       Varin = 7-day means of driver variables and ETrs

C     Output:
C       Varout = multi-year means of 7-day driver variables and ETrs

C     Author: Mike Hobbins
C     Date: May 21, 2024

      IMPLICIT  NONE
      INTEGER   mlen(12),nrows,ncols,bnd_s,bnd_n,bnd_w,bnd_e,bnd_lo,
     +          bnd_hi,yrstart,yrend,row,col,colin,rowin,mon,day,var,
     +          year,i,count
      REAL      Lat,Lon,Elev,tot
      CHARACTER varname(5)*4,root*99,maskfile*99,tracker*17,outfile*99,
     +          infile*99

      DATA mlen/31,29,31,30,31,30,31,31,30,31,30,31/
      DATA varname/'ETrs','T_2m','SWdn','SpHm','U_2m'/

C     Define nrows and ncols, geographic bounds of analysis, and root path.
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
      REAL    Varin(yrend-yrstart+1,nrows,ncols),Varout(ncols)

 
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
          IF((rowin.NE.nrows-row+1).OR.(colin.NE.col)) THEN
            PRINT*,'Static grid input error: col ',colin,'; row ',rowin
            GOTO 99
          ENDIF
C         Generate geographic analysis mask
          IF((Lat.GT.bnd_s).AND.(Lat.LT.bnd_n).AND.
     +       (Lon.GT.bnd_w).AND.(Lon.LT.bnd_e).AND.
     +       (Elev.GT.bnd_lo).AND.(Elev.LT.bnd_hi))
     +       Mask(row,col) = 1
10      ENDDO                                    ! end col loop
20    ENDDO                                      ! end row loop
      CLOSE(10)


C     Loop through months
      DO 90 mon = 1, 12

C       Loop through days
        DO 80 day = 1, mlen(mon)

          WRITE(tracker,'(A5,I2.2,A8,I2.2)') 'Mon: ',mon,',  Day: ',day
          PRINT*,tracker

C         Loop through variables
          DO 70 var = 1, 5

C           Initialize input array
            Varin = -9999.

C           Open daily output file for climatological mean variable
            WRITE(outfile,'(A24,A35,A4,A1,2I2.2,A4)')
     +          root,'ETrc_attribution/CONUS/7_day_climo/',
     +          varname(var),'_',mon,day,'.asc'
            OPEN(UNIT=11,FILE=outfile,STATUS='UNKNOWN')
            WRITE(11,'(A9)') 'ncols 464'
            WRITE(11,'(A9)') 'nrows 224'
            WRITE(11,'(A14)') 'xllcorner -125'
            WRITE(11,'(A12)') 'yllcorner 25'
            WRITE(11,'(A14)') 'cellsize 0.125'
            WRITE(11,'(A19)') 'nodata_value -9999.'

C           Loop through years
            DO 30 year = yrstart, yrend

C             Open daily input file for variable
              WRITE(infile,'(A24,A32,A4,A1,I4,2I2.2,A4)')
     +              root,'ETrc_attribution/CONUS/7_day_ts/',
     +              varname(var),'_',year,mon,day,'.asc'
!             use March 1 for February 29 in non-leap years
              IF((mon.EQ.2).AND.(day.EQ.29)
     +           .AND.(MOD(year,4).NE.0)) THEN
                WRITE(infile,'(A24,A32,A4,A1,I4,A8)')
     +                root,'ETrc_attribution/CONUS/7_day_ts/',
     +                varname(var),'_',year,'0301.asc'
              ENDIF
              OPEN(UNIT=12,FILE=infile,STATUS='OLD',ERR=99)
C             Read time-series grids
              DO i = 1, 6                        ! read headers
                READ(12,*)
              ENDDO
              DO row = 1, nrows                  ! read data
                READ(12,*) (Varin(year-yrstart+1,row,col),col=1,ncols)
              ENDDO
              CLOSE(12)

30          ENDDO                                ! end year loop


C           Estimate daily long-term mean
            DO 60 row = 1,nrows

C             Initialise output array
              Varout = -9999.

              DO 50 col = 1,ncols

                IF(Mask(row,col).EQ.1) THEN
                  tot = 0.
                  count = 0
                  DO 40 year = yrstart, yrend
                    tot = tot + Varin(year-yrstart+1,row,col)
                    count = count + 1
40                ENDDO
                  Varout(col) = tot / count
                ELSE
                  Varout(col) = -9999.
                ENDIF

50            ENDDO                            ! end col loop

C             Write output data
              WRITE(11,*) (Varout(col),col=1,ncols)

60          ENDDO                              ! end row loop

C           Close output file
            CLOSE(11)

70        ENDDO                                ! end var loop

80      ENDDO                                  ! end day loop

90    ENDDO                                    ! end mon loop

99    END