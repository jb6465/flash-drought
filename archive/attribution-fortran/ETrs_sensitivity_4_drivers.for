      PROGRAM ETRS_SENSITIVITY_4_DRIVERS
C     **********************************
C     Uses 42-year climatological means of 7-day surfaces of North American
C     Land Data Assimilation System phase-2 (NLDAS-2) drivers to generate
C     surfaces of daily sensitivities of ETrs (simplified relative to the
C     strict Penman-Monteith equation recommended by the ASCE-EWRI) to each
C     of its drivers.

C     ASCE-EWRI (2005), The ASCE Standardized Reference Evapotranspiration
C       Equation, Task Committee on Standardization of Reference
C       Evapotranspiration, Ed. Richard G. Allen.

C     Input data:
C       Elev = elevation [m]
C       Lat  = latitude [degs]
C     and 30-year means of the following input variables:
C       T_2m = 2-m daily mean air temperature [K]
C       SWdn = surface downward shortwave radiation flux [W/m2]
C       SpHm = 2-m specific humidity [kg/kg]
C       U_2m = 2-m wind speed [m/sec]

C     Output data are the sensitivities of ETrs to the input variables:
C       sens_T = to T_2m [(mm/day)/K]
C       sens_R = to SWdn [(mm/day)/(W/m2)]
C       sens_q = to SpHm [(mm/day)/(kg/kg)]
C       sens_U = to U_2m [(mm/day)/(m/sec)]

C     Notes:
C     * all input variables in the main program are capitalized (not so in
C         subroutine).
C     * most parameters are in SI units; exceptions are Cn and ETrs [mm/day].
C     * monlen array specifies a 29-day February, as sensitivities are needed
C         for all days including February 29 in leap years.
C     * jday array specifies a 28-day February as solar calculations need an
C         accurate, integer Julian Day count, and most years have only 28 days
C         in February.

C     Author: Mike Hobbins
C     Date: May 21, 2024


      IMPLICIT  NONE
      INTEGER   monlen(12),jday(12),nrows,ncols,bnd_s,bnd_n,bnd_w,bnd_e,
     +          row,col,Colin,Rowin,mon,day,jdaysub,var,i,year,offset
      REAL      bnd_lo,bnd_hi,Lon,pi
      CHARACTER varname(4)*4,root*99,maskfile*99,outpath*99,outfile*99,
     +          infile*99

      DATA varname/'T_2m','SWdn','SpHm','U_2m'/
      DATA monlen/31,29,31,30,31,30,31,31,30,31,30,31/
      DATA jday/0,31,59,90,120,151,181,212,243,273,304,334/

C     Define nrows [-] and ncols [-], geographic bounds of analysis, pi,
C     and output path.
      PARAMETER (nrows  = 224,
     +           ncols  = 464,
     +           bnd_s  = 25,
     +           bnd_n  = 53,
     +           bnd_w  = -125,
     +           bnd_e  = -67,
     +           bnd_lo = -90.,
     +           bnd_hi = 4500.,
     +           pi     = 3.1415926536,
     +           root='/Volumes/mhobbins_drobo/ETrc_attribution/CONUS/')

      INTEGER   mask(nrows,ncols)
      REAL      Latdeg(nrows),Elev(nrows,ncols),avg(4,nrows,ncols),
     +          Datavar(ncols)
      REAL*8    sens_T(ncols),sens_R(ncols),sens_q(ncols),sens_U(ncols)


C     Create mask to limit analysis
      Mask = 0
      WRITE(maskfile,'(A35,A31)') 
     +      '/System/Volumes/Data/data/mhobbins/',
     +      'mac159/scripts/gtopomean15k.asc'
      OPEN(UNIT=10,FILE=maskfile,STATUS='OLD',ERR=99)
      DO row = nrows, 1, -1
        DO col = 1, ncols
C         Read input data: Colin [-], column; Rowin [-], row;
C         Lat [deg], latitude; Lon [deg], longitude; Elev [m], elevation
          READ(10,*) Colin,Rowin,Latdeg(row),Lon,Elev(row,col)
          IF((Rowin.NE.nrows-row+1).OR.(Colin.NE.col)) THEN
            PRINT*,'Grid input error: colin ',Colin,'; rowin ',Rowin
            GOTO 99
          ENDIF
          IF((Latdeg(row).GE.bnd_s).AND.(Latdeg(row).LE.bnd_n).AND.
     +       (Lon.GE.bnd_w).AND.(Lon.LE.bnd_e).AND.
     +       (Elev(row,col).GE.bnd_lo).AND.(Elev(row,col).LE.bnd_hi))
     +       mask(row,col) = 1
        ENDDO
      ENDDO
      CLOSE(10)


      DO 60 mon = 1, 12                          ! loop through months

        DO 50 day = 1, monlen(mon)               ! loop through days

          PRINT*,'Month: ',mon,'; Day: ',day

          avg = 0.

C         Set the jday variable for the subroutine
          jdaysub = jday(mon) + day
          IF((mon.EQ.2).AND.(day.EQ.29)) jdaysub = 60

          DO 20 var = 1, 4                       ! loop through variables

C           Open variable sensitivity output files, write headers
            WRITE(outfile,'(A47,A23,A4,A1,2I2.2,A4)')
     +            root,'7_day_sensitivity/Sens_',
     +            varname(var),'_',mon,day,'.asc'
            OPEN(UNIT=10+var,FILE=outfile,STATUS='UNKNOWN')
            WRITE(10+var,'(A6,I3)') 'ncols ',ncols
            WRITE(10+var,'(A6,I3)') 'nrows ',nrows
            WRITE(10+var,'(A10,I4)') 'xllcorner ',bnd_w
            WRITE(10+var,'(A10,I2)') 'yllcorner ',bnd_s
            WRITE(10+var,'(A14)') 'cellsize 0.125'
            WRITE(10+var,'(A19)') 'nodata_value -9999.'

C           Read climatologic input data for the day
C             1 = T2m, daily mean T at 2 m [K]
C             2 = SWd, downward shortwave radiation at surface [W/m2]
C             3 = SH, specific humidity at 2 m [kg/kg]
C             4 = U2m, wind speed at 10 m [m/sec]
            WRITE(infile,'(A47,A12,A4,A1,2I2.2,A4)')
     +            root,'7_day_climo/',varname(var),'_',mon,day,'.asc'
            OPEN(UNIT=14+var,FILE=infile,STATUS='OLD')
            DO i = 1, 6
              READ(14+var,*)
            ENDDO
C           Read variable climatological means, filter for no data 
            DO row = 1, nrows                    ! loop through rows
              READ(14+var,*) (Datavar(col),col=1,ncols)
              DO col = 1, ncols                  ! loop through cols
                IF(mask(row,col).NE.0) THEN
                  avg(var,row,col) = Datavar(col)
                ELSE
                  avg(var,row,col) = -9999.
                ENDIF
              ENDDO                              ! end col loop
            ENDDO                                ! end row loop
            CLOSE(14+var)

20        ENDDO                                  ! end variable loop


C         Calculate sensitivities

          DO 40 row = 1, nrows                   ! loop through rows
            sens_T = -9999.
            sens_R = -9999.
            sens_q = -9999.
            sens_U = -9999.
            DO 30 col = 1, ncols                 ! loop through cols
              IF((mask(row,col).NE.0).AND.(avg(1,row,col).GE.0.).AND.
     +           (avg(2,row,col).GE.0.).AND.(avg(3,row,col).GE.0.).AND.
     +           (avg(4,row,col).GE.0.))
     +           CALL ETrc_sensitivity(avg(1,row,col),avg(2,row,col),
     +                  avg(3,row,col),avg(4,row,col),Elev(row,col),
     +                  Latdeg(row)*pi/180.,jdaysub,sens_T(col),
     +                  sens_R(col),sens_q(col),sens_U(col))
30          ENDDO                                ! end col loop


C           Write output surfaces to files

C           Write sensitivities to output files
            WRITE(11,*) (sens_T(col),col=1,ncols)
            WRITE(12,*) (sens_R(col),col=1,ncols)
            WRITE(13,*) (sens_q(col),col=1,ncols)
            WRITE(14,*) (sens_U(col),col=1,ncols)
40        ENDDO                                  ! end row loop

C         Close output files
          DO i = 11,14
            CLOSE(i)
          ENDDO

50      ENDDO                                    ! end day loop

60    ENDDO                                      ! end month loop

99    END


      SUBROUTINE ETrc_sensitivity(T2m,SWdn,SpHm,U2m,Elev,Lat,jday,
     +                            sens_T,sens_R,sens_q,sens_U)
C     ******************************************************************
C     Estimates ETrs at the pixel in question. Inputs are scalars drawn
C     from vectors in the main program. Output is a scalar that is
C     inserted into a vector of ETrs in the main program.

      IMPLICIT NONE
      INTEGER  jday
      REAL     T2m,SWdn,SpHm,U2m,Elev,Lat,Cn,Cd,alb
      REAL*8   sens_T,sens_R,sens_q,sens_U,pi,Gsc,sigma,Rs,T,Pa,gamma,
     +         delta,esat,w,ea,Rns,G,decl,ws,dr,Ra,Rso,RsRso,fcd,Rnl,Rn,
     +         ETrs,A,B,C,ATa,ATb,dATadT,dATbdT,dAdT,dAdq,dAdR,dBdT,
     +         dBdq,dBdU,dCdT,dCdU


C     Define physical parameter constants:

C     Numerator constants, Cn [mm.K.sec/(day.m.Pa)], referred to in
C     ASCE-EWRI as units of [K.mm.sec3/(Kg.day)], which are the same.
C= this text may not be true any more.
      PARAMETER (Cn = 1.6)  !for tall reference crop (alfalfa) on a daily basis

C     Denominator constants, Cd [sec/m], equivalent to a wind function parameter
      PARAMETER (Cd = 0.38) !for tall reference crop (alfalfa) on a daily basis

C     Pi, pi [-]
      PARAMETER (pi = 3.1415926536)

C     Albedo of reference crop, alb [-]
      PARAMETER (alb = 0.23)

C     Solar constant, Gsc [J/m2/hour]
      PARAMETER (Gsc = 4.92e6)

C     Stefan-Boltzmann constant, sigma [J/K4/m2/day]
      PARAMETER (sigma = 4.901e-3)


C     Downwelling SW radiation Rs [J/m2/day] = f(SWdn [W/m2])
      Rs = SWdn * 86400.

C     Daily mean temperature at 2 m, T [degC] = f(T2m [K])
      T = T2m - 273.15

C     Pressure, Pa [Pa] = f(elev [m])
      Pa = 101300. * ((293. - 0.0065 * Elev) / 293.)**5.26

C     Psychrometric constant, gamma [Pa/K] = f(Pa [Pa])
      gamma = 0.000665 * Pa

C     Slope of the saturated vapor pressure curve at T, delta [Pa/K] = f(T [degC])
      delta = 2.503e6 * exp(17.27 * T / (T + 237.3)) / (T + 237.3)**2.

C     Saturated vapor pressure, esat [Pa] = f(T [degC])
      esat = 610.8 * exp(17.27 * T / (T + 237.3))

C     Mixing ratio, w [kg/kg] = f(SpHm [kg/kg])
      w = SpHm / (1. - SpHm)

C     Actual vapor pressure, ea [Pa] = f(Patm [Pa], w [kg/kg])
      ea = Pa * w / (0.622 + w)

C     Limit ea to below or equal to esat
      IF(ea.GT.esat) ea = esat

C     Net SW radiation downward, Rns [J/m2/day] = f(alb [-], Rs [J/m2/day])
      Rns = (1. - alb) * Rs

C     Ground heat flux, G [J/m2/day], assumed to be zero at the daily time-step
      G = 0.

C     Declination, decl [rad] = f(pi [-], jday [-])
      decl = 0.409 * sin((2. * pi * jday / 365) - 1.39)

C     Sunset hour angle, ws [rad] = f(Lat [rad], decl [rad])
      ws = acos (-1. * tan(Lat) * tan(decl))

C     Inverse relative distance of earth from sun, dr [-] = f(pi [-], jday [-])
      dr = 1. + 0.033 * cos(2. * pi * jday / 365)

C     Extra-terrestrial (TOA) SW radiation, Ra [J/m2/day] = f(Gsc [J/m2/hour],
C     dr [-], ws [rad], Lat [rad], decl [rad])
      Ra = (24. / pi) * Gsc * dr * (ws * sin(Lat) * sin(decl)
     +     + cos(Lat) * cos(decl) * sin(ws))

C     Clear-sky SW radiation at surface, Rso [J/m2/day] = f(Ra [J/m2/day],
C     Elev [m])
      Rso = Ra * (0.75 + 2.e-5 * Elev)

C     Relative shortwave radiation, RsRso = f(Rs [J/m2/day], Rso [J/m2/day])
      RsRso = Rs / Rso
      IF(RsRso.LT.0.3) RsRso = 0.3
      IF(RsRso.GT.1.)  RsRso = 1.

C     Cloudiness function, fcd [-] = f(swratio [-])
      fcd = 1.35 * RsRso - 0.35
      IF(fcd.LT.0.05) fcd = 0.05
      IF(fcd.GT.1.)   fcd = 1.

C     Net LW radiation upward, Rnl [J/m2/day] = f(sigma [J/K4/m2/day],
C     fcd [-], ea [Pa], Tmx [degC], Tmn [degC])
      Rnl = sigma * fcd * (0.34 - 0.14 * (0.001 * ea)**0.5)
     +      * (T + 273.15)**4.

C     Net radiation, Rn [J/m2/day] = f(Rns [J/m2/day], Rnl [J/m2/day])
      Rn = Rns - Rnl

C     ASCE-EWRI ETrs, ETrs [mm/day] = f(0.408e-6 [mm.m2/J], delta [Pa/K],
C     gamma [Pa/K], u2 [m/sec], Rn [J/m2/day], G [J/m2/day], T [degC],
C     esat [Pa], ea [Pa], Cn [K.mm.sec3/kg/day], Cd [sec/m])
      ETrs = (0.408e-6 * delta * (Rn - G)
     +       + gamma * (Cn / (T + 273.)) * U2m * (esat - ea))
     +       / (delta + gamma * (1. + Cd * U2m))

      A = 0.408e-6 * delta * (Rn - G)
      B = gamma * (Cn / (T + 273.)) * U2m * (esat - ea)
      C = delta + gamma * (1. + Cd * U2m)

      ATa = 0.408e-6 * (2.503e6 * exp(17.27 * T / (T + 237.3))
     +    / (T + 237.3)**2.) 
      ATb = (Rns - (sigma * fcd * (0.34 - 0.14 * (0.001 * ea)**0.5)
     +    * (T + 273.15)**4.))
      dATadT = -(2.04245 * exp((17.27 * T) / (T + 237.3))
     +       * (T**2. * (T + 237.3)**2. - 2049.09 *(T + 237.3)**3.
     +       + 474.6 * T * (T + 237.3)**2. + 56311.3 * (T + 237.3)**2.))
     +       / (T + 237.3)**7.
      dATbdT = 0.0177088 * fcd * (-76.7982 + ea**0.5) * sigma
     +       * (T + 273.15)**3.
      dAdT = ATa * dATbdT + ATb * dATadT
      dAdq = (2.41719e-9 * delta * fcd * Pa * (SpHm - 1.)
     +     * sigma * (273.15 + T)**4.)
     +     / ((SpHm - 1.) * ((Pa * SpHm)
     +     / (1.6455 + SpHm))**0.5 * (1.6455 + SpHm)**2.)
      dAdR = -0.0352512 * (alb - 1.) * delta
     +     + (0.000210686 * delta * (-76.7982 + ea**0.5)
     +     * sigma * (273.15 + T)**4) / Rso

      B = (0.000665 * 101300. * ((293. - 0.0065 * Elev) / 293.)**5.26)
     +  * (Cn / (T + 273.15)) * U2m
     +  * ((610.8 * exp(17.27 * T / (T + 237.3)))
     +  - (Pa * (SpHm / (1. - SpHm)) / (0.622 + (SpHm / (1. -SpHm)))))
      dBdT = -1.*(67.3431 * Cn * (610.8 * exp((17.27 * T) / (237.3 + T))
     +     - (Pa * SpHm)
     +     / ((1 - SpHm) * (0.622 + SpHm / (1 - SpHm)))) * U2m)
     +     / (273.15 + T)**2. + (41133.2 * Cn * exp((17.27 * T)
     +     / (237.3 + T)) * (-(17.27 * T) / (237.3 + T)**2. + 17.27
     +     / (237.3 + T)) * U2m) / (273.15 + T)
      dBdq = (Cn * Pa * 293.157 * (1 - SpHm) * U2m)
     +     / ((SpHm - 1.) * (1.6455 + SpHm)**2. * (273.15 + T))
      dBdU = (0.000665 * Cn * Pa * (610.8 * exp((17.27 * T)/(237.3 + T))
     +     - (Pa * w)/(0.622 + w)))/(273.15 + T)

      C = 2.503e6 * exp(17.27 * T / (T + 237.3)) / (T + 237.3)**2.
     +   + 0.000665 * Pa * (1. + Cd * U2m)
      dCdT = -(5006000. * exp((17.27 * T)/(237.3 + T)) * (-1811.79 + T))
     +     / (237.3 + T)**4.
      dCdU = 0.000665 * Cd * Pa

C     Compiling sensitivity components into driver sensitivities
      sens_T = (dAdT + dBdT - (A + B) * dCdT / C) / C
      sens_R  = dAdR / C
      sens_q  = (dAdq + dBdq) / C
      sens_U  = ((C * dBdU) - (A + B) * dCdU) / C**2.

      RETURN
      END