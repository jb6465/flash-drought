      PROGRAM GENERATE_T2M_U2M_HISTORY
C     ********************************
C     Generates a multi-year timeseries of daily T_2m as the average of
C     daily Tmin and Tmax and of daily U_2m as a function of daily U10m
C     and sensor height.

C     Author: Mike Hobbins
C     Date: May 21, 2024

      IMPLICIT  NONE
      INTEGER   mlen(12),nrows,ncols,bnd_s,bnd_n,bnd_w,bnd_e,bnd_lo,
     +          bnd_hi,yrstart,yrend,row,col,colin,rowin,year,mon,lpday,
     +          day,date,i
      REAL      Lat,Lon,Elev
      CHARACTER root*99,maskfile*99,outfile*99,infile*99
      DATA      mlen/31,28,31,30,31,30,31,31,30,31,30,31/

C     Define nrows and ncols, geographic bounds of analysis, and path.
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
     +           root  = '/Volumes/mhobbins_drobo/NLDAS_drivers/daily/')

      INTEGER Mask(nrows,ncols)
      REAL    Tmin(ncols),Tmax(ncols),T_2m(ncols),U10m(ncols),
     +        U_2m(ncols)

 
C     Read in geographic data to set mask for analysis
      Mask = 0
      WRITE(maskfile,'(A35,A31)') 
     +      '/System/Volumes/Data/data/mhobbins/',
     +      'mac159/scripts/gtopomean15k.asc'
      OPEN(UNIT=10,FILE=maskfile,STATUS='OLD',ERR=99)
      DO 20 row = nrows, 1, -1
        DO 10 col = 1, ncols
C         Read input data: Col; Row; Lat [deg], latitude; Lon [deg],
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


C     Loop through years
      DO 50 year = yrstart, yrend

C       Loop through months
        DO 40 mon = 1, 12

          lpday = 0
          IF((mon.EQ.2).AND.(MOD(year,4).EQ.0.)) lpday = 1

C         Loop through days
          DO 30 day = 1, mlen(mon) + lpday

            date = year * 10000 + mon * 100 + day
            PRINT*,'Date:',date

C           Open daily T_2m and U_2m output files
            WRITE(outfile,'(A44,A5,I8,A4)') root,'T_2m_',date,'.asc'
            OPEN(UNIT=11,FILE=outfile,STATUS='UNKNOWN')
            WRITE(outfile,'(A44,A5,I8,A4)') root,'U_2m_',date,'.asc'
            OPEN(UNIT=12,FILE=outfile,STATUS='UNKNOWN')
            DO i = 11,12
              WRITE(i,'(A9)') 'ncols 464'
              WRITE(i,'(A9)') 'nrows 224'
              WRITE(i,'(A14)') 'xllcorner -125'
              WRITE(i,'(A12)') 'yllcorner 25'
              WRITE(i,'(A14)') 'cellsize 0.125'
              WRITE(i,'(A19)') 'nodata_value -9999.'
            ENDDO

C           Open daily Tmin, Tmax, and U10m input files
            WRITE(infile,'(A44,A5,I8,A4)') root,'Tmin_',date,'.asc'
            OPEN(UNIT=13,FILE=infile,STATUS='OLD',ERR=99)
            WRITE(infile,'(A44,A5,I8,A4)') root,'Tmax_',date,'.asc'
            OPEN(UNIT=14,FILE=infile,STATUS='OLD',ERR=99)
            WRITE(infile,'(A44,A5,I8,A4)') root,'U10m_',date,'.asc'
            OPEN(UNIT=15,FILE=infile,STATUS='OLD',ERR=99)
            DO i = 1, 6                          ! read headers
              READ(13,*)
              READ(14,*)
              READ(15,*)
            ENDDO

C           Read daily Tmin, Tmax, and U10m grids
            DO row = 1, nrows
              Tmin = -9999.
              Tmax = -9999.
              T_2m = -9999.
              U10m = -9999.
              U_2m = -9999.
              READ(13,*) (Tmin(col),col=1,ncols)
              READ(14,*) (Tmax(col),col=1,ncols)
              READ(15,*) (U10m(col),col=1,ncols)

C             Estimate T_2m and U_2m and write to output files
              DO col = 1, ncols
                IF(Mask(row,col).EQ.1) THEN
                   T_2m(col) = (Tmin(col) + Tmax(col)) / 2.
                   U_2m(col) = U10m(col) * 4.87 / log(67.8 * 10. - 5.42)
                ENDIF
              ENDDO
              WRITE(11,*) (T_2m(col),col=1,ncols)
              WRITE(12,*) (U_2m(col),col=1,ncols)
            ENDDO
            DO i = 11,15
              CLOSE(i)
            ENDDO

30        ENDDO                                  ! end day loop

40      ENDDO                                    ! end mon loop

50    ENDDO                                      ! end year loop

99    END