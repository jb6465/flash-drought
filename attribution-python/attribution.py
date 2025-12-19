'''
###########################################################################################
## BUREAU OF METEOROLOGY
## AUSTRALIAN CLIMATE SERVICE
## DROUGHT AND CHANGES IN ARIDITY HAZARD TEAM
##
## DATE:             Nov-2025
## SCRIPT:           decomposition.py
## AUTHOR:           jessica.bhardwaj@bom.gov.au
##
## DESCRIPTION:      Script for decomposing the sensivity of daily ETrs to its drivers and 
##                   then deriving spatial and temporal contributions for event analysis.
##
############################################################################################
'''

def compute_ET0(Rs, T, Tmax, Tmin, Tdew, Patm, U, q, RHmean, RHmax, RHmin, elev, **kwargs):
    """
    Inputs:
    Arrays (surface daily aggregated from hourly) 
    *ensure input units are consistent with the below  |  Input units   |
    - Rs: downwave shortwave radiation flux               [J/(m2.day)]
    - T: mean temperature                                 [°C]
    - Tmax: max temperature                               [°C]
    - Tmin: min temperature                               [°C]
    - Tdew: mean dewpoint temperature                     [°C]
    - Patm: mean pressure                                 [Pa] 
    - U: mean windspeed                                   [m/s]
    - q: specific humidity                                [kg/kg]
    - RHmean: mean relative humidity                      [%]
    - RHmax: max relative humidity                        [%]
    - RHmin: min temperature                              [%]
    - elev: model derived elevation                       [m]
    - lat: input lat, used for solar calculations         [degrees]
    - jday: day of year for each time step, used for solar calculations
    
    Key word arguments (with string inputs listed in order of preference)
    - derive_2m_windspeed: switch to derive 2m windspeed from 10m inputs [True, False]
    - short_tall_crop_switch: string to switch between tall or short ET0 ['tall', 'short']
    - T2M_model_or_derived: string to switch between T2M input ['derived', 'model']
    - Patm_model_or_derived: string to switch between Patm input ['model', 'derived']
    - esat_method: string to select esat calculation method ['TmaxTmin', 'Tmean']
    - eact_method: string to select eact calculation method ['Tdew', 'Patm_q', 'RHmaxRHmean', 'RH']
    - Rnl_method: string to select Rnl calculation method ['TmaxTmin', 'Tmean']
    
    Additional notes
    * Using mean air temperature instead of maximum and minimum air temperatures yields a lower saturation vapour pressure, and hence a lower vapour pressure difference (es - ea), and a lower ET0 estimate.
    * Table 3 in ASCE and Box 7 in FAO56 both stipulate preferenced order of calculating actual vapor pressure as (1)  derived from dewpoint temperature (2) derived from maximum and minimum relative humidity (3) derived from maximum relative humidity (4)  derived from mean relative humidity [less recommended]
    
    Output:
    Xarrays of:
    - daily reference evapotranspiration computed using ASCE’s standardized reference evapotranspiration equations
    """
    import numpy as np
    import xarray as xr

    ##### Access key word arguments and set defaults in case of no input #####
    derive_2m_windspeed = kwargs.get('derive_2m_windspeed', False)
    short_tall_crop_switch = kwargs.get('short_tall_crop_switch', 'tall')
    T2M_model_or_derived = kwargs.get('T2M_model_or_derived', 'derived')
    Patm_model_or_derived = kwargs.get('Patm_model_or_derived', 'model')
    esat_method = kwargs.get('esat_method', 'TmaxTmin')
    eact_method = kwargs.get('eact_method', 'Tdew')
    Rnl_method = kwargs.get('Rnl_method', 'TmaxTmin')
    
    ##### Constants #####
    Cn = 1.6 if short_tall_crop_switch == 'tall' else 0.9         # Cn: numerator crop constant [K.m.s3/(kg.day)]
    Cd = 0.38 if short_tall_crop_switch == 'tall' else 0.34       # Cd: denominator crop constant [s/m]
    Gsc = 4.92e6                                                  # Gsc: solar constant [J/(m2.hr)]
    sigma = 4.901e-3                                              # sigma: Stefan-Boltzmann constant [J/(K4.m2.day)]
    albedo = 0.23                                                 # albedo [unitless]

    ##### ET0 computation prelude #####
    # T: daily mean temperature [°C]  {ASCE, Eq.2 / FAO56 Eq. 9}
    T = 0.5*(Tmax + Tmin) if T2M_model_or_derived == 'derived' else T

    # U: daily 2-m windspeed [m/s]  {ASCE, Eq.33 / FAO56 Eq. 47}
    U = U * (4.87 / np.log(67.8 * 10.0 - 5.42)) if derive_2m_windspeed == True else U
    
    # Patm: surface pressure [Pa], preference to use direct model outputs otherwise derive as f(elev[m])  {ASCE, Eq.3 / FAO56 Eq. 7}
    Patm = 101300 * ((293. - 0.0065 * elev) / 293.)**5.26 if Patm_model_or_derived == 'derived' else Patm
    
    # gamma: psychrometric constant [Pa/°C]  {ASCE, Eq.4 / FAO56 Eq. 8}
    gamma = 0.665e-3 * Patm

    # esat: saturated vapor pressure [Pa], preference to use Tmax and Tmin otherwise derive from Tmean f(T[°C])  {ASCE, Eq.3 / FAO56 Eq. 11}
    esat = 0.5* (610.8 * np.exp((17.27*Tmax/(Tmax+237.3))) + 610.8 * np.exp((17.27*Tmin/(Tmin+237.3))) ) if esat_method == 'TmaxTmin' else 610.8 * np.exp((17.27*T/(T+237.3))) 

    # eact: actual vapor pressure [Pa] = f(Patm[Pa], w[kg/kg])
    w = q / (1. - q)
    if eact_method == 'Tdew':
        eact = 610.8 * np.exp((17.27*Tdew/(Tdew+237.3)))
    elif eact_method == 'RHmaxRHmean':
        esatTmax = 610.8 * np.exp((17.27*Tmax/(Tmax+237.3)))
        esatTmin = 610.8 * np.exp((17.27*Tmin/(Tmin+237.3)))
        eact = 0.5 * (RHmean * esatTmax / 100 + RHmean * esatTmin / 100)
    elif eact_method == 'RHmean':
        eact = RHmean * esat / 100
    else:
        eact = (Patm * w) / (0.622 + w) #note this method is the same as Tdew ~ just a rearranged equation. Not provided in ASCE or FAO but derived physically by Mike.
    eact = eact.clip(max=esat) #clip eact ≤ esat

    # delta: slope of the saturated vapor pressure curve at T [Pa/K] = f(T[°C])  {ASCE, Eq.5 / FAO56 Eq. 13}
    delta = (4098. * 610.8 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3)**2    

    # G: ground heat flux  [MJ/(m2.day)], on daily time scales, G ≈ 0 {ASCE Eq. 30 / FAO56 Eq. 42]
    G = 0       

    # decl: declination [radians] = f(pi[unitless], jday[unitless]) {ASCE Eq. 24 / FAO56 Eq. 24}
    jday = xr.DataArray(Rs.time.dt.dayofyear.values, coords=[('day', Rs.time.dt.dayofyear.values)], dims=["day"])
    decl = 0.409 * np.sin((2. * np.pi * jday / 365) - 1.39)                                                                   
    # ws: sunset hour angle ws[radians] = f(latitude[rad], dec[rad]) {ASCE Eq. 27 / FAO56 Eq. 25}
    ws = np.arccos(-1 * np.tan(np.deg2rad(Rs.lat)) * np.tan(decl))

    # dr: inverse relative distance of earth from sun [unitless]  {ASCE Eq. 23 / FAO56 Eq. 23}
    dr = 1. + 0.033 * np.cos(2. * np.pi * jday / 365)                           
    
    # Ra: extra-terrestrial (TOA) SW radiation [J/(m2.day)] = f(Gsc [J/m2/hour], dr [unitless], ws [rad], Lat [rad], decl [rad])  {ASCE Eq. 21 / FAO56 Eq. 21}
    Ra = (24. / np.pi) * Gsc * dr * (ws * np.sin(np.deg2rad(Rs.lat)) * np.sin(decl) + np.cos(np.deg2rad(Rs.lat)) * np.cos(decl) * np.sin(ws))
    
    # Rso: clear-sky SW radiation at surface [J/m2/day] = f(Ra[J/(m2.day)], elev[m])  {ASCE Eq. 19 / FAO56 Eq. 37}
    Rso = Ra * (0.75 + 2e-5 * elev)

    # fcd: cloudiness function [unitless]  {ASCE Eq. 45 / FAO56 pg. 52}
    RsRso = xr.DataArray((Rs.values)/(Rso.values), coords=[('time', Rs.time.values), ('lat', Rs.lat.values), ('lon', Rs.lon.values)], dims=["time", "lat", "lon"])
    RsRso = RsRso.clip(min=0.3, max=1.0)
    fcd = 1.35 * RsRso - 0.35 
    fcd = fcd.clip(min=0.05, max=1.0) #clip values between 0.05 and 1.0

    # Rns: Net shortwave radiation  [J/(m2.day)]  {ASCE Eq. 16 / FAO56 Eq. 38}
    Rns = (1 - albedo) * Rs

    # Rnl: net LW radiation upward[J/(m2.day)] = f(sigma [J/(K4.m2.day)], eact[Pa], T[°C])  {ASCE Eq. 44 / FAO56 Eq. 39}
    Rnl = sigma * fcd * (0.34 - 0.14 * (0.001 * eact)**0.5) * (0.5 *((Tmax + 273.15)**4. + (Tmin + 273.15)**4.)) if Rnl_method == 'TmaxTmin' else sigma * fcd * (0.34 - 0.14 * (0.001 * eact)**0.5) * (T + 273.15)**4.

    # Rn: net radiation [J/(m2.day)]  {ASCE Eq. 15 / FAO56 Eq. 40}
    Rn = Rns - Rnl

    ##### Compute ET0 #####
    A = (0.408*1e-6) * (delta) * (Rn - G)            # A left numerator
    B = (gamma) * (Cn/(T+273)) * U * ((esat-eact))   # B right numerator
    C = ((delta) + (gamma) * (1+Cd*U))               # C denominator
    
    ET0 = (A+B)/C
    
    ET0 = ET0.rename(f"ET0_{short_tall_crop_switch}_crop")
    ET0.attrs["units"] = "mm/day"
    ET0.attrs["long_name"] = f"Reference evapotranspiration for {short_tall_crop_switch} crops"
    ET0.attrs["description"] = "Computed using ASCE standardized ET0 equation"
    ET0.attrs["function_kwargs"] = str({k: v for k, v in kwargs.items()})
    
    return ET0
    
def compute_driver_sensivities(Rs, T, U, q, Patm, elev, short_tall_crop_switch):
    """
    Inputs:
    !!! ensure inputs are in the same format as for ET0 computation and DOY means for an appropriate reference period
    Output:
    Xarrays of:
    - sensitvity of daily reference evapotranspiration to its input drivers (T, Rs, q, U)

    Additional notes
    * These sensitvity equations are derived for ET0 calculated with a 
    """
    import numpy as np
    import xarray as xr
    
    #constants
    albedo = 0.23
    sigma = 4.901e-3 
    Cn = 1.6 if short_tall_crop_switch == 'tall' else 0.9
    Cd = 0.38 if short_tall_crop_switch == 'tall' else 0.34
    Gsc = 4.92e6

    #equations ~ same as for compute ET0 function but need to be replicated here for the computation with doy means required for sensitivity
    Rns = (1 - albedo) * Rs
    sigma = 4.901e-3
    gamma = 0.665e-3 * Patm
    esat = 610.8 * np.exp((17.27*T/(T+237.3))) 

    w = q / (1. - q)
    eact = (Patm * w) / (0.622 + w)
    eact = eact.clip(max=esat)

    delta = (4098. * 610.8 * np.exp(17.27 * T / (T + 237.3))) / (T + 237.3)**2
    G = 0       
    jday = xr.DataArray(Rs.doy.values, coords=[('day', Rs.doy.values)], dims=["day"])
    decl = 0.409 * np.sin((2. * np.pi * jday / 365) - 1.39)                                                                   
    ws = np.arccos(-1 * np.tan(np.deg2rad(Rs.lat)) * np.tan(decl))
    dr = 1. + 0.033 * np.cos(2. * np.pi * jday / 365)                           
    Ra = (24. / np.pi) * Gsc * dr * (ws * np.sin(np.deg2rad(Rs.lat)) * np.sin(decl) + np.cos(np.deg2rad(Rs.lat)) * np.cos(decl) * np.sin(ws))
    Rso = Ra * (0.75 + 2e-5 * elev)
    RsRso = xr.DataArray((Rs.values)/(Rso.values), coords=[('doy', Rs.doy.values), ('lat', Rs.lat.values), ('lon', Rs.lon.values)], dims=["doy", "lat", "lon"])
    RsRso = RsRso.clip(min=0.3, max=1.0)
    fcd = 1.35 * RsRso - 0.35 
    fcd = fcd.clip(min=0.05, max=1.0) #clip values between 0.05 and 1.0
    
    Rns = (1 - albedo) * Rs
    Rnl = sigma * fcd * (0.34 - 0.14 * (0.001 * eact)**0.5) * (T + 273.15)**4.
    Rn = Rns - Rnl

    A = (0.408*1e-6) * (delta) * (Rn - G)
    ATa = 0.408e-6 * (2.503e6 * np.exp(17.27 * T / (T + 237.3))/ (T + 237.3)**2.) 
    ATb = (Rns - (sigma * fcd * (0.34 - 0.14 * (0.001 * (eact))**0.5)* (T + 273.15)**4.))
    dATadT = -(2.04245 * np.exp((17.27 * T) / (T + 237.3))* (T**2. * (T + 237.3)**2. - 2049.09 *(T + 237.3)**3 + 474.6 * T * (T + 237.3)**2. + 56311.3 * (T + 237.3)**2.))/ (T + 237.3)**7.
    dATbdT = 0.0177088 * fcd * (-76.7982 + (eact)**0.5) * sigma * (T + 273.15)**3.
    dAdT = ATa * dATbdT + ATb * dATadT
    dAdq = (2.41719e-9 * delta * fcd * Patm * (q - 1.)* sigma * (273.15 + T)**4.)/ ((q - 1.) * ((Patm * q)/ (1.6455 + q))**0.5 * (1.6455 + q)**2.)
    dAdR = -0.0352512 * (albedo - 1.) * (delta) + (0.000210686 * (delta) * (-76.7982 + (eact)**0.5) * sigma * (273.15 + T)**4) / Rso.values

    B = (gamma) * (Cn/(T+273)) * U * ((esat-eact))
    B = B.transpose('doy', 'lat', 'lon')
    dBdT = -1.*(67.3431 * Cn * (610.8 * np.exp((17.27 * T) / (237.3 + T)) - (Patm * q) / ((1 - q) * (0.622 + q / (1 - q)))) * U) / (273.15 + T)**2. + (41133.2 * Cn * np.exp((17.27 * T) / (237.3 + T)) * (-(17.27 * T) / (237.3 + T)**2. + 17.27 / (237.3 + T)) * U) / (273.15 + T)
    dBdq = (Cn * Patm * 293.157 * (1 - q) * U)/ ((q - 1.) * (1.6455 + q)**2. * (273.15 + T))
    dBdU = (0.000665 * Cn * Patm * (610.8 * np.exp((17.27 * T)/(237.3 + T)) - (Patm * w)/(0.622 + w)))/(273.15 + T)

    C = ((delta) + (gamma) * (1+Cd*U))   
    dCdT = -(5006000. * np.exp((17.27 * T)/(237.3 + T)) * (-1811.79 + T))/ (237.3 + T)**4.
    dCdU = 0.000665 * Cd * Patm

    sens_T = (dAdT + dBdT - (A + B) * dCdT / C) / C
    sens_R  = dAdR / C
    sens_q  = (dAdq + dBdq) / C
    sens_U  = ((C * dBdU) - (A + B) * dCdU) / C**2
    return sens_T.rename("T_sens"), sens_R.rename("R_sens"), sens_q.rename("q_sens"), sens_U.rename("U_sens")

def compute_rolling_xday(daily, n_rolling_days, centre_window=True):
    x_day = daily.rolling(time=n_rolling_days, center=centre_window).mean()
    return x_day
    
def compute_doy_climo(input_xr): 
    climo = input_xr.assign_coords(doy=input_xr['time'].dt.dayofyear).groupby('doy').mean('time')
    return climo
    
def compute_daily_anom(input_xr): 
    climo = compute_doy_climo(input_xr)
    x_day_anom = input_xr.groupby('time.dayofyear') - climo.rename({'doy': 'dayofyear'})
    return x_day_anom