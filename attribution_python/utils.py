def netRadiation(rsds, elev, tasmax, tasmin, ea, lat_name='lat'):
    """
    Calculates net radation according to Allen et al 2005.
    Input:
        - rsds: (array) downwelling shortwave radiation at the surface [W/m2]
        - elev: (array) elevation [m]
        - tasmax: (array) maximum temperature [˚C]
        - tasmin: (array) minimum temperature [˚C]
        - ea: (array) actual vapour pressure (can be calculated with actVapourPressure()) [kPa]
        - lat_name: name of latutide coordinate in rsds dataset. Default 'lat'.
    Output:
        - Rn: (array) net radiation [MJ/m2/day]
    """
    # Net SW radiation downward, Rns [MJ/m2/day] = f(alb [-], Rs [MJ/m2/day]), Allen et al 2005 equaton 43
    # Albedo of reference crop, alb [-]
    # Allen et al 205 equation 43 - constant value for daily and hourly.
    alb = 0.23

    # Stefan-Boltzmann constant, sigma [MJ/K4/m2/day]
    # SB constant = 5.67*10^-8 J s^-1 m^-2 K^-4 *10^-6 J->MJ * 86400 s day^-1 = 4.8988 MJ K^-4 m^-2 day^-1
    # Allen et al 2005 equation 44 sigma = 2.42 * 10^-10 MJ K^-4 m^-2 h^-1
    sigma = 4.901e-9

    _jday = rsds.time.dt.dayofyear.values
    jday = xr.DataArray(_jday, coords=[('day', _jday)], dims=["day"])
    
    # Downwelling SW radiation Rs [MJ/m2/day] = f(rsds [W/m2])
    # 24 hrs/day * 60 mins/hr * 60 sec/min = 86400 sec/day
    # J -> MJ => * 10^-6
    # rsds  = W m^-2 = J m^-2 s^-1 * 86400 s day^-1 * 10^-6 MJ J^-1 = MJ m^-2 day^-1
    Rs = rsds * 86400. / 1.e6
    Rns = (1. - alb) * Rs
                
    # Declination, decl [rad] = f(pi [-], jday [-])
    # Allen et al 2005 equation 51 solar declination
    decl = 0.409 * np.sin((2. * np.pi * jday / 365) - 1.39)
        
    # Sunset hour angle, ws [rad] = f(Lat [rad], decl [rad])
    # Allen et al 2005 equation 59 - is this correct? should it be tan of lat???
    Lat  = Rns[lat_name]
    ws = np.arccos(-1 * np.sin(np.deg2rad(Lat)) * np.tan(decl))
        
    # Inverse relative distance of earth from sun, dr [-] = f(pi [-], jday [-]), Allen et al 2005 equation 50
    dr = 1. + 0.033 * np.cos(2. * np.pi * jday / 365)

    # Solar constant, Gsc [MJ/m2/hour]
    # 1.366 kW m^-2 * 1000 = 1366 W m^-2 = 1366 J s^-1 m^-2 * 3600 s hr^-1 *10^-6 = 4.9176 MJ m^-2 hr^-1
    Gsc = 4.92

    # Extra-terrestrial (TOA) SW radiation, Ra [MJ/m2/day] = f(Gsc [MJ/m2/hour],
    # dr [-], ws [rad], Lat [rad], decl [rad]), Allen et al 2005 equation 48. SW solar radiation in the absence of an atmosphere. Well-behaved function of the day of the year, time of day, lat and lon.
    Ra = (24. / np.pi) * Gsc * dr * (ws * np.sin(Lat * np.pi/180) * np.sin(decl) + np.cos(Lat * np.pi/180) * np.cos(decl) * np.sin(ws))
    
    # Clear-sky SW radiation at surface, Rso [MJ/m2/day] = f(Ra [MJ/m2/day], Elev [m]), Allen et al 2005 equation 47
    Rso = Ra * (0.75 + 2.e-5 * elev)
    
    # Relative shortwave radiation, RsRso = f(Rs [MJ/m2/day], Rso [MJ/m2/day])
    # Allen et al 2005 limited to between 0.3 and 1.0
    # Resolve misaligned units 'day' and 'time'. Align Rso using the 'time' dimension from Rs.
######
    # Rso_aligned = Rso.reindex(time=Rs['time'])
    Rso = Rso.rename({'day':'time'})
    Rso['time']=Rs['time']
    
    # Perform the computation of RsRso
    RsRso = Rs / Rso
#######    
    # Clip values between 0.3 and 1.0
    RsRso = RsRso.clip(min=0.3, max=1.0)

    # Cloudiness function, fcd [-] = f(swratio [-]), See Allen et al 2005 equation 44/45 - limited to between 0.05 and 1.0, dimensionless
    fcd = 1.35 * RsRso - 0.35
    # Clip values between 0.05 and 1.0
    fcd = fcd.clip(min=0.05, max=1.0)
        
    # Net LW radiation upward, Rnl [MJ/m2/day] = f(sigma [MJ/K4/m2/day],
    # fcd [-], ea [kPa], Tmx [degC], Tmn [degC])
    # Net LW radiation is the difference between the LW radiation radiated upward from the standardised surface and LW radiation radiated downward from the atmospere. This method was introduced by Brunt (1932, 1952) and uses near surface vapour pressure to predict net surface emissivity. Allen et al 2005 equation 44
    Rnl = sigma * fcd * (0.34 - 0.14 * ea**0.5) * 0.5 * ((tasmax+273.16)**4. + (tasmin+273.16)**4.)
    
    # Net radiation, Rn [MJ/m2/day] = f(Rns [MJ/m2/day], Rnl [MJ/m2/day])
    # Allen et al 2005 equation 42
    Rn = Rns - Rnl
    return Rn
    
def wind2m(sfcWind,scaling_factor=None):
    """
    Estimates 2m surface wind from near surface wind (usually 10m).
    Input:
        - sfcWind: (array) 10m near surface winds
        - scaling_factor: (array) gridded roughness map for Australia based on observational station data (Zhang et al., 2022) to convert 10m wind speed datasets to 2m
    Output:
        - u2: (array) 2m surface wind estimate
    """
    if scaling_factor:
        # Use scaling factors from gridded roughness map for Australia based on observational station data (Zhang et al., 2022) to convert 10m wind speed datasets to 2m
        u2 = sfcWind * scaling_factor
    else:
        # Wind speed at 2 m, u2 [m/sec] = f(sfcWind [m/sec])
        # Log wind profile to 2 m from 10 m uses Allen et al 2005 equation 33, which is for wind measurements taken above a short grass or similar surface. Eqn B.14 in Appendix B is given for wind measured above alfalfa or similar vegetation having about 0.5 m height. Wind speed data collected above 2m are acceptable for use following adjustment to 2 m, and may be preferred if vegetation commonly exceeds 0.5 m. Measurement at a greater height reduces the influence of the taller vegetation.
        u2 = sfcWind * 4.87 / np.log(67.8 * 10. - 5.42) 
    return u2


def meanTemperature(tasmax,tasmin):        
    # Air temperature at 2 m, T [degC] = f(Tmx [degC], Tmn [degC])
    # This gives approximation of daily mean 2 m T. As stated by Allen et al 2005, this is the preferred method to calculate the mean, rather than an average of hourly T measurements, to provide consistency across all data sets.
    tas = (tasmax + tasmin) / 2.
    
    return tas

def satVapourPressure(tasmax, tasmin):
    # Saturated vapor pressure, esat [kPa] = f(Tmx [degC], Tmn [degC])
    # See Chapter 3 - Meteorological Data for formulae, Allen 2005, Equations 6 and 7
    # esat = 0.6108 exp{(17.27T/(T+237.3)}
    esat = (0.6108 * np.exp(17.27 * tasmax / (tasmin + 237.3)) + 0.6108 * np.exp(17.27 * tasmin / (tasmin + 237.3))) / 2.
    return esat

def actVapourPressure(esat, psl=None, huss=None, hurs=None):
    # Check that either 'hurs' is provided or both 'psl' and 'huss' are provided
    if (hurs is None and (psl is None or huss is None)):
        raise ValueError("Provide either 'hurs', or both 'psl' and 'huss'.")

    if hurs is not None:
        # Calculate ea from hurs (Relative Humidity in %) and esat
        ea = hurs * esat / 100  # Convert hurs from percentage to fraction and calculate ea
        
    elif psl is not None and huss is not None:
        # Specific humidity at 2 m (huss)
        w = huss / (1. - huss)
        
        # Actual vapor pressure (ea) = f(psl [Pa], w [kg/kg])
        ea = 0.001 * psl * w / (0.622 + w)  # Convert psl from Pa to kPa
        
    # Clip ea values to below or equal to esat
    ea = ea.clip(max=esat)
    
    return ea 

def calculate_ASCE_pmpet(Rn,t,u2,esat,ea,ps):
        
    # Define physical parameter constants:
    # See Allen at el 2005 for details of the parameters needed for this calculation
    print("Files read. Start calculation...")
    # Numerator constants, Cn_alf and Cn_grs [mm.K.sec/(day.m.kPa)], referred to in ASCE-EWRI as units of [K.mm.sec3/(Mg.day)], which are the same.
    Cn_alf = 1600.    # for tall reference crop (alfalfa) on a daily basis
    Cn_grs = 900.     # for short reference crop (grass) on a daily basis
        
    # Denominator constants, Cd_alf and Cd_grs [sec/m], equivalent to a wind function parameter
    Cd_alf = 0.38     # for tall reference crop (alfalfa) on a daily basis
    Cd_grs = 0.34     # for short reference crop (grass) on a daily basis
        
    # Psychometric constant, gamma [kPa/degC] = f(Pa [kPa])
    gamma = 0.000665 * ps * 10**(-3)

    # Ground heat flux, G [MJ/m2/day], assumed to be zero at the daily time-step
    G = 0.
        
    # Slope of the saturated vapor pressure curve at T,
    # delta [kPa/degC] = f(T [degC])
    # See Chapter 3 - Meteorological Data for formulae, Allen 2005, Equation 5    
    delta = 2503. * np.exp(17.27 * t / (t + 237.3)) / (t + 237.3)**2.

    # ASCE-EWRI ETrs, ETrs [mm/day] = f(0.408 [mm.m2/MJ], delta [kPa/K],
    # gamma [kPa/K], u2 [m/sec], Rn [MJ/m2/day], G [MJ/m2/day], T [degC],
    # esat [kPa], ea [kPa], Cn_alf [K.mm.sec3/Mg/day], Cd_alf [sec/m])
    ETrs = (0.408 * delta * (Rn - G) + gamma * (Cn_alf / (t + 273.)) * u2 * (esat - ea)) / (delta + gamma * (1. + Cd_alf * u2))
        
    # ASCE-EWRI ETos, ETos [mm/day] = f(0.408 [mm.m2/MJ], delta [kPa/K],
    # gamma [kPa/K], u2 [m/sec], Rn [MJ/m2/day], G [MJ/m2/day], T [degC],
    # esat [kPa], ea [kPa], Cn_grs [K.mm.sec3/Mg/day], Cd_grs [sec/m])
    ETos = (0.408 * delta * (Rn - G) + gamma * (Cn_grs / (t + 273.)) * u2 * (esat - ea)) / (delta + gamma * (1. + Cd_grs * u2))

    # Attributes
    ETrs.attrs['units'] = 'mm'
    ETrs.attrs['long_name'] = 'ASCE-EWRI tall reference crop ET for 0.5 m alfalfa [mm/day]'
    ETos.attrs['units'] = 'mm'
    ETos.attrs['long_name'] = 'ASCE-EWRI short reference crop ET for 12 cm grass [mm/day]'
    
    # List of coordinates to drop    
    coords_to_drop = ['height', 'level_height', 'model_level_number', 'sigma']
    # Drop coordinates only if they are present in the dataset
    ETrs = ETrs.drop([coord for coord in coords_to_drop if coord in ETrs.coords]).rename('ETrs').to_dataset()
    ETos = ETos.drop([coord for coord in coords_to_drop if coord in ETos.coords]).rename('ETos').to_dataset()

    ET = xr.merge([ETrs, ETos])
    # if 'rlat' in ET.dims:
    #     ET = ET.drop(['lat','lon']).rename_dims({'rlat': 'lat','rlon': 'lon'})
        
    ET.attrs['description'] = 'Uses daily or monthly data from downscaled CMIP6 outputs to calculate mean daily ETrs and ETos, using the version of the Penman-Monteith expression recommended by the American Society of Civil Engineers Environmental Water Resources Institute (ACSCE-EWRI), and known as the ASCE Standardized Reference ET equation. Reference: ASCE-EWRI (2005): The ASCE Standardized Reference Evapotranspiration Equation, Task Committee on Standardization of Reference Evapotranspiration, Ed. Richard G. Allen'
    return ET
    
def convertTemperature(temp):
    if temp.attrs['units'] not in ('C', 'K'):
        raise ValueError("The 'unit' argument for tasmax and tasmin is neither 'C' (Celsius) or 'K' (Kelvin).")
    if temp.attrs['units'] == 'K':
        temp = temp - 273.16
        temp.attrs['units'] = 'C'
    return temp
        