      PROGRAM ATTRIBUTION
C     *******************
C     Uses data derived from NLDAS-2 to calculate the depth contribution
C     from each of the four ETrs driving variables to daily ETrs anomalies.

C     Input data, from NLDAS datafiles:
C       T_2m = 2-m daily air temperature [K]
C       SWdn = surface downward shortwave radiation flux [W/m2]
C       SpHm = 2-m specific humidity [kg/kg]
C       U_2m = 2-m wind speed [m/sec]
C       ETrs = daily reference ET [mm]

C     Output:
C       cont_[var] = mm contribution of each variable [var] to the
C                    daily ETrs anomaly

C     Author: Mike Hobbins
C     Date: May 21, 2024

      IMPLICIT  NONE
      INTEGER   nrows,ncols,bnd_s,bnd_n,bnd_w,bnd_e,bnd_lo,bnd_hi,year,
     +          mon,day,row,col,colin,rowin,i,var
      REAL      Lat,Lon,Elev
      CHARACTER varname(5)*4,root1*99,root2*99,specs*10,maskfile*99,
     +          infile*99,outfile*99
      DATA      varname/'T_2m','SWdn','SpHm','U_2m','ClEr'/

C     Define nrows and ncols and root path.
      PARAMETER (nrows  = 224,
     +           ncols  = 464,
     +           bnd_s  = 25,
     +           bnd_n  = 53,
     +           bnd_w  = -125,
     +           bnd_e  = -67,
     +           bnd_lo = -90,
     +           bnd_hi = 4500,
     +           root1  = '/Volumes/mhobbins_drobo/',
     +           root2  = 'ETrc_attribution/CONUS/')

      INTEGER    Mask(nrows,ncols)
      REAL       ETrs_anom(nrows,ncols),sens(4,ncols),mean(4,ncols),
     +           ts(4,ncols),cont(5,nrows,ncols),ETrs_mean(ncols),
     +           ETrs(ncols)

C     Read date from arguments
      CALL get_command_argument(1,specs)
      READ(specs,'(I4,2I2.2)') year,mon,day

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
      ETrs_anom = -9999.
      sens      = -9999.
      mean      = -9999.
      ts        = -9999.
      cont      = -9999.

C     Estimate ETrs anomaly

C     Open mean grid for ETrs
      WRITE(infile,'(A24,A23,A17,2I2.2,A4)')
     +      root1,root2,'7_day_climo/ETrs_',mon,day,'.asc'
      OPEN(UNIT=10,FILE=infile,STATUS='OLD',ERR=99)
C     Open time-series grid for ETrs
      WRITE(infile,'(A24,A21,I4,2I2.2,A4)')
     +      root1,'ETrc_data/daily/ETrs_',year,mon,day,'.asc'
      OPEN(UNIT=11,FILE=infile,STATUS='OLD',ERR=99)
C     Read mean and time-series grids
      DO i = 1, 6                                ! skip file headers
        READ(10,*)
        READ(11,*)
      ENDDO
      DO 20 row = 1, nrows
        READ(10,*) (ETrs_mean(col),col=1,ncols)
        READ(11,*) (ETrs(col),col=1,ncols)
        DO col = 1,ncols
          IF((ETrs(col).GT.-9000.).AND.(ETrs_mean(col).GT.-9000.))
     +      ETrs_anom(row,col) = ETrs(col) - ETrs_mean(col) 
        ENDDO
20    ENDDO
      CLOSE(10)
      CLOSE(11)

C     Loop through all four ETrs-input variables {T_2m, SWdn, SpHm, U_2m}
      DO 50 var = 1, 4

C       Open variable sensitivity grids
        WRITE(infile,'(A24,A23,A23,A4,A1,2I2.2,A4)')
     +        root1,root2,'7_day_sensitivity/Sens_',
     +        varname(var),'_',mon,day,'.asc'
        OPEN(UNIT=12,FILE=infile,STATUS='OLD',ERR=99)
C       Open variable mean grids
        WRITE(infile,'(A24,A23,A12,A4,A1,2I2.2,A4)')
     +        root1,root2,'7_day_climo/',
     +        varname(var),'_',mon,day,'.asc'
        OPEN(UNIT=13,FILE=infile,STATUS='OLD',ERR=99)
C       Open variable time-series grids
        WRITE(infile,'(A24,A20,A4,A1,I4,2I2.2,A4)')
     +        root1,'NLDAS_drivers/daily/',
     +        varname(var),'_',year,mon,day,'.asc'
        OPEN(UNIT=14,FILE=infile,STATUS='OLD',ERR=99)

C       Skip file headers
        DO i = 1, 6
          READ(12,*)
          READ(13,*)
          READ(14,*)
        ENDDO

C       Read sensitivity, mean, and time-series grids
        DO 40 row = 1, nrows

          READ(12,*) (sens(var,col),col=1,ncols)
          READ(13,*) (mean(var,col),col=1,ncols)
          READ(14,*) (ts(var,col),col=1,ncols)

C         Estimate daily contribution to ETrs anomaly
          DO 30 col = 1, ncols
            IF((ts(var,col).GT.-9000.).AND.(mean(var,col).GT.-9000.)
     +         .AND.(ETrs_anom(row,col).GT.-9000.)) THEN
              cont(var,row,col) = (ts(var,col) - mean(var,col))
     +                            * sens(var,col)
              IF(var.eq.4) THEN                  ! calculate closure error
                cont(5,row,col) = ETrs_anom(row,col) - (cont(1,row,col)
     +                            + cont(2,row,col) + cont(3,row,col)
     +                            + cont(4,row,col))
              ENDIF
            ENDIF

30        ENDDO                                  ! end col loop

40      ENDDO                                    ! end row loop

        CLOSE(12)
        CLOSE(13)
        CLOSE(14)

50    ENDDO                                      ! end var loop


C     Output contribution depths for each driver and closure error
      DO 60 var = 1, 5

        WRITE(outfile,'(A24,A23,A19,A4,A1,I4,2I2.2,A4)')
     +        root1,root2,'contributions/Cont_',
     +        varname(var),'_',year,mon,day,'.asc'
        OPEN(UNIT=15,FILE=outfile,STATUS='UNKNOWN',ERR=99)
        WRITE(15,'(A6,I3)') 'ncols ',ncols
        WRITE(15,'(A6,I3)') 'nrows ',nrows
        WRITE(15,'(A14)') 'xllcorner -125'
        WRITE(15,'(A12)') 'yllcorner 25'
        WRITE(15,'(A14)') 'cellsize 0.125'
        WRITE(15,'(A19)') 'nodata_value -9999.'
        DO row = 1, nrows
          WRITE(15,*) (cont(var,row,col),col=1,ncols)
        ENDDO
        CLOSE(15)

60    ENDDO

C     Output ETrs anomaly depth
      WRITE(outfile,'(A24,A23,A24,I4,2I2.2,A4)')
     +      root1,root2,'contributions/ETrs_anom_',year,mon,day,'.asc'
      OPEN(UNIT=16,FILE=outfile,STATUS='UNKNOWN',ERR=99)
      WRITE(16,'(A6,I3)') 'ncols ',ncols
      WRITE(16,'(A6,I3)') 'nrows ',nrows
      WRITE(16,'(A14)') 'xllcorner -125'
      WRITE(16,'(A12)') 'yllcorner 25'
      WRITE(16,'(A14)') 'cellsize 0.125'
      WRITE(16,'(A19)') 'nodata_value -9999.'
      DO row = 1, nrows
        WRITE(16,*) (ETrs_anom(row,col),col=1,ncols)
      ENDDO
      CLOSE(16)

99    END
