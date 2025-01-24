      PROGRAM GENERATE_T2M_U2M
C     ************************
C     Generates daily T_2m as the average of daily Tmin and Tmax and
C     daily U_2m as a function of daily U10m.

C     Author: Mike Hobbins
C     Date: May 21, 2024

      IMPLICIT  NONE
      INTEGER   nrows,ncols,bnd_s,bnd_n,bnd_w,bnd_e,bnd_lo,bnd_hi,row,
     +          col,colin,rowin,i
      REAL      Lat,Lon,Elev
      CHARACTER root*99,specs*8,date*8,maskfile*99,outfile*99,infile*99

C     Define nrows and ncols, geographic bounds of analysis, and path.
      PARAMETER (nrows  = 224,
     +           ncols  = 464,
     +           bnd_s  = 25,
     +           bnd_n  = 53,
     +           bnd_w  = -125,
     +           bnd_e  = -67,
     +           bnd_lo = -90,
     +           bnd_hi = 4500,
     +           root  = '/Volumes/mhobbins_drobo/NLDAS_drivers/daily/')

      INTEGER Mask(nrows,ncols)
      REAL    Tmin(ncols),Tmax(ncols),U10m(ncols),T_2m(ncols),
     +        U_2m(ncols)

 
C     Read date from arguments
      CALL get_command_argument(1,specs)
      READ(specs,'(A8)') date

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

C     Open daily T_2m and U_2m output files
      WRITE(outfile,'(A44,A5,A8,A4)') root,'T_2m_',date,'.asc'
      OPEN(UNIT=11,FILE=outfile,STATUS='UNKNOWN')
      WRITE(outfile,'(A44,A5,A8,A4)') root,'U_2m_',date,'.asc'
      OPEN(UNIT=12,FILE=outfile,STATUS='UNKNOWN')
      DO 30 i = 11,12
        WRITE(i,'(A9)') 'ncols 464'
        WRITE(i,'(A9)') 'nrows 224'
        WRITE(i,'(A14)') 'xllcorner -125'
        WRITE(i,'(A12)') 'yllcorner 25'
        WRITE(i,'(A14)') 'cellsize 0.125'
        WRITE(i,'(A19)') 'nodata_value -9999.'
30    ENDDO

C     Open daily Tmin, Tmax, and U10m input files
      WRITE(infile,'(A44,A5,A8,A4)') root,'Tmin_',date,'.asc'
      OPEN(UNIT=13,FILE=infile,STATUS='OLD',ERR=99)
      WRITE(infile,'(A44,A5,A8,A4)') root,'Tmax_',date,'.asc'
      OPEN(UNIT=14,FILE=infile,STATUS='OLD',ERR=99)
      WRITE(infile,'(A44,A5,A8,A4)') root,'U10m_',date,'.asc'
      OPEN(UNIT=15,FILE=infile,STATUS='OLD',ERR=99)
      DO 40 i = 1, 6                      ! read headers
        READ(13,*)
        READ(14,*)
        READ(15,*)
40    ENDDO

C     Read daily Tmin, Tmax, and U10m grids
      DO 50 row = 1, nrows
        T_2m = -9999.
        Tmin = -9999.
        Tmax = -9999.
        U_2m = -9999.
        U10m = -9999.
        READ(13,*) (Tmin(col),col=1,ncols)
        READ(14,*) (Tmax(col),col=1,ncols)
        READ(15,*) (U10m(col),col=1,ncols)

C       Estimate T_2m and U_2m and write to output files
        DO col = 1, ncols
          IF(Mask(row,col).EQ.1) THEN
            T_2m(col) = (Tmin(col) + Tmax(col)) / 2.
            U_2m(col) = U10m(col) * 4.87 / log(67.8 * 10. - 5.42)
          ENDIF
        ENDDO
        WRITE(11,*) (T_2m(col),col=1,ncols)
        WRITE(12,*) (U_2m(col),col=1,ncols)
50    ENDDO

      DO i = 11,15
        CLOSE(i)
      ENDDO

99    END