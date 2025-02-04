      PROGRAM NET_CONTRIBUTION
C     ************************
C     For a series of defined intervals, estimates net contribution of
C     each driver to the anomaly in ETrs.
C
C     Add cases to the "Select case" code to add extra output periods. 
C
C     Author: Mike Hobbins
C     Date: May 21, 2024 

      IMPLICIT  NONE
      INTEGER   monlen(12),nrows,ncols,bnd_s,bnd_n,bnd_w,bnd_e,bnd_lo,
     +          bnd_hi,year,mon,day,row,col,colin,rowin,daysback,date,
     +          var,i,output,lpdy
      REAL      Lat,Lon,Elev
      CHARACTER varname(6)*4,root*99,specs*8,maskfile*99,infile*99,
     +          tsname*6,outsumfile*99,outavgfile*99
      DATA      monlen/31,28,31,30,31,30,31,31,30,31,30,31/
      DATA      varname/'T_2m','SWdn','SpHm','U_2m','ClEr','ETrs'/

C     Define nrows and ncols, and input and output paths.
      PARAMETER (nrows  = 224,
     +           ncols  = 464,
     +           bnd_s  = 25,
     +           bnd_n  = 53,
     +           bnd_w  = -125,
     +           bnd_e  = -67,
     +           bnd_lo = -90,
     +           bnd_hi = 4500,
     +           root='/Volumes/mhobbins_drobo/ETrc_attribution/CONUS/')

      INTEGER   Mask(nrows,ncols),count(nrows,ncols)
      REAL      sumdat(6,nrows,ncols),dat(ncols),avgdat(ncols)


C     Read in date from argument
      CALL get_command_argument(1,specs)
      READ(specs,'(I4,2I2)') year,mon,day

C     Read in geographic data to set mask for analysis
      Mask = 0
      WRITE(maskfile,'(A35,A31)') 
     +      '/System/Volumes/Data/data/mhobbins/',
     +      'mac159/scripts/gtopomean15k.asc'
      OPEN(UNIT=9,FILE=maskfile,STATUS='OLD',ERR=99)
      DO 10 row = nrows, 1, -1
        DO col = 1, ncols
C         - read input data: Col; Row; Lat [deg], latitude; Lon [deg],
C         longitude; Elev [m], elevation
          READ(9,*) colin,rowin,Lat,Lon,Elev
C         Generate geographic analysis mask
          IF((Lat.GT.bnd_s).AND.(Lat.LT.bnd_n).AND.
     +       (Lon.GT.bnd_w).AND.(Lon.LT.bnd_e).AND.
     +       (Elev.GT.bnd_lo).AND.(Elev.LT.bnd_hi))
     +       Mask(row,col) = 1
        ENDDO                                    ! end col loop
10    ENDDO                                      ! end row loop
      CLOSE(9)


C     Initialize arrays
      count  = 0
      sumdat = 0

      DO 80 daysback = 1, 360

        date = year * 10000 + mon * 100 + day

        PRINT*,'Date: ',date

C       Read in date's contribution file
        DO 40 var = 1, 6
          WRITE(infile,'(A47,A19,A4,A1,I8,A4)')
     +          root,'contributions/Cont_',varname(var),'_',date,'.asc'
          IF(var.EQ.6) WRITE(infile,'(A47,A24,I8,A4)')
     +          root,'contributions/ETrs_anom_',date,'.asc'
          OPEN(UNIT=11,FILE=infile,STATUS='OLD',ERR=99)

C         Skip header of input file
          DO i = 1, 6
            READ(11,*)
          ENDDO

C         Loop through rows
          DO 30 row = 1, nrows

C           Read in data
            READ(11,*) (dat(col),col=1,ncols)

C           Loop through cols
            DO 20 col = 1, ncols

              IF(dat(col).NE.-9999.) THEN
C               Add 1 to count of data
                count(row,col) = count(row,col) + 1
C               Add data to sum of data
                sumdat(var,row,col) = sumdat(var,row,col) + dat(col)
              ENDIF

20          ENDDO                            ! end col loop

30        ENDDO                              ! end row loop

          CLOSE(11)

40      ENDDO                                ! end variable loop

C       Check if date is required for output
        output = 0
        SELECT CASE (daysback)
          CASE (7)
            output = 1
            tsname = '_01wk_'
          CASE (14)
            output = 1
            tsname = '_02wk_'
          CASE (30)
            output = 1
            tsname = '_01mn_'
          CASE (60)
            output = 1
            tsname = '_02mn_'
          CASE (90)
            output = 1
            tsname = '_03mn_'
          CASE (180)
            output = 1
            tsname = '_06mn_'
          CASE (360)
            output = 1
            tsname = '_12mn_'
          CASE DEFAULT
            output = 0
        END SELECT

        IF(output.EQ.1) THEN

C         Write net contribution for each variable
          DO 70 var = 1, 6

C           Open output file and write header
            WRITE(outsumfile,'(A47,A16,A4,A6,I8,A4)') 
     +           root,'results/Netcont_',varname(var),tsname,date,'.asc'
            WRITE(outavgfile,'(A47,A16,A4,A6,I8,A4)') 
     +          root,'results/Avgcont_',varname(var),tsname,date,'.asc'
            IF(var.EQ.6) THEN
              WRITE(outsumfile,'(A47,A23,A6,I8,A4)')
     +             root,'results/Net_ETo_anomaly',tsname,date,'.asc'
              WRITE(outavgfile,'(A47,A23,A6,I8,A4)')
     +             root,'results/Avg_ETo_anomaly',tsname,date,'.asc'
            ENDIF
            OPEN(UNIT=12,FILE=outsumfile,STATUS='UNKNOWN')
            OPEN(UNIT=13,FILE=outavgfile,STATUS='UNKNOWN')
            DO i = 12, 13
              WRITE(i,'(A9)') 'ncols 224'
              WRITE(i,'(A9)') 'nrows 464'
              WRITE(i,'(A14)') 'xllcorner -125'
              WRITE(i,'(A12)') 'yllcorner 25'
              WRITE(i,'(A14)') 'cellsize 0.125'
              WRITE(i,'(A19)') 'NODATA_value -9999.'
            ENDDO

C           Loop through rows
            DO 60 row = 1, nrows

C             Initialize array
              avgdat = -9999.

C             Loop through cols
              DO 50 col = 1, ncols

                IF(count(row,col).GT.0)
     +              avgdat(col) = sumdat(var,row,col) / count(row,col)
                IF(count(row,col).LE.0)
     +              sumdat(var,row,col) = -9999.

50            ENDDO                                ! end col loop

C             Write to output files
              WRITE(12,*) (sumdat(var,row,col),col=1,ncols)
              WRITE(13,*) (avgdat(col),col=1,ncols)

60          ENDDO                                  ! end row loop

            CLOSE(12)
            CLOSE(13)

70        ENDDO                                    ! end variable loop

        ENDIF

C       Generate date for previous day for next loop
        day = day - 1
        IF(day.eq.0) THEN
          mon = mon - 1
          IF(mon.EQ.0) THEN
            year = year - 1
            mon = 12
          ENDIF
          lpdy = 0
          IF((mon.EQ.2).AND.(MOD(year,4).EQ.0)) lpdy = 1
          day = monlen(mon) + lpdy
        ENDIF

80    ENDDO

99    END