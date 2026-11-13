"""
Python version of WASP model implemented in SURFEX/wasp_flux.F90 (version from MNH-V5-7-2)
"""

# pylint: disable=missing-function-docstring,invalid-name

import sys
import math
from . import Cst


def QSAT_SEAWATER(PT,PP,CQSAT):
    cst = Cst()
    ZFOES  = PSAT(PT,CQSAT)
    ZFOES  = 0.98*ZFOES
    # vapor pressure reduction of 2% over saline seawater could have a significant
    # impact on the computation of surface latent heat flux under strong wind
    # conditions (Zeng et al, 1998).
    #
    ZWORK1 = ZFOES/PP
    ZWORK2    = cst.Rd/cst.Rv
    PQSAT = ZWORK2*ZWORK1 / (1.+(ZWORK2-1.)*ZWORK1)
    return PQSAT

def QSAT(PT,PP,CQSAT):
    cst = Cst()
    ZFOES  = PSAT(PT,CQSAT)
    ZWORK1 = ZFOES/PP
    ZWORK2 = cst.Rd/cst.Rv
    PQSAT = ZWORK2*ZWORK1 / (1.+(ZWORK2-1.)*ZWORK1)
    return PQSAT

def PSAT(PT,CQSAT):
    cst = Cst()
    ZALP  = cst.alpw
    ZBETA = cst.betaw
    ZGAM  = cst.gamw
    if(CQSAT=='NEW' and PT<=cst.Tt):
        ZALP  = cst.alpi
        ZBETA = cst.betai
        ZGAM  = cst.gami

    PPSAT = math.exp( ZALP - ZBETA/PT - ZGAM*math.log(PT) )
    return PPSAT

def CHARNOCK_WAGE(PWIND,PWAGE):
    ZCOEFU = [ 0.70, -2.52 ]
    ZCOEFA2 = [ 2.27, -6.67E-02 ]
    ZCOEFB2 = [ -2.41, 4.30E-02 ]
    ZCOEFA = [ -9.202, 2.265, -0.134, 2.35E-03 ]
    ZCOEFB =[ -0.4124, -0.2225, 0.01178, -1.616E-04 ]
    ZPOLYU = [ 0.0981, -4.13E-03, 4.34E-5, 1.16E-08 ]
    ZLIMCHAR = 0.018
    ZLIMCHAR2 = 0.002
    ZLIMCHAR1 = 0.1
    PCHARNWA = ZCOEFU[1-1]*(PWIND**ZCOEFU[2-1])
    if PWIND >= 7.0:
        ZAA = ZCOEFA[1-1] + ZCOEFA[2-1]*PWIND    \
                        + ZCOEFA[3-1]*PWIND**2 \
                        + ZCOEFA[4-1]*PWIND**3
        ZBB = ZCOEFB[1-1] + ZCOEFB[2-1]*PWIND    \
                        + ZCOEFB[3-1]*PWIND**2 \
                        + ZCOEFB[4-1]*PWIND**3
        PCHARNWA = ZAA * (PWAGE**ZBB)

    if PWIND >= 23.0:
        ZAA = ZCOEFA2[1-1] + ZCOEFA2[2-1]*PWIND
        ZBB = ZCOEFB2[1-1] + ZCOEFB2[2-1]*PWIND
        PCHARNWA= ZAA * (PWAGE**ZBB)
        if PCHARNWA < ZLIMCHAR:
            PCHARNWA = ZLIMCHAR

    if PWIND >= 25.0:
        PCHARNWA = ZPOLYU[1-1] + ZPOLYU[2-1]*PWIND    \
                           + ZPOLYU[3-1]*PWIND**2 \
                           + ZPOLYU[4-1]*PWIND**3
        if PCHARNWA < ZLIMCHAR2:
            PCHARNWA = ZLIMCHAR2

    # Final check, to avoid too large CHAR
    if PCHARNWA > ZLIMCHAR1:
        PCHARNWA = ZLIMCHAR1
    return PCHARNWA

def  WIND_THRESHOLD(PWIND,PUREF,XCISMIN, XVMODMIN, LALDTHRES):
    if not LALDTHRES:
        #  minimum value for exchange coefficients computations : 1m/s / 10m
        PWIND_NEW = max(PWIND , 0.1 * min(10.,PUREF) )
    else:
        #  minimum value for exchange coefficients computations : 1m/s / 10m
        PWIND_NEW = max( XVMODMIN, math.sqrt(math.pow(PWIND,2) + math.pow(XCISMIN*PUREF),2) )
    return PWIND_NEW

def SURFACE_RI(XRIMAX, PTG, PQS, PEXNS, PEXNA, PTA, PQA,PZREF, PUREF,
               PDIRCOSZW, PVMOD,XCISMIN, XVMODMIN, LALDTHRES ):
    cst = Cst()
    ZTHVA=PTA/PEXNA*( 1.+(cst.Rv/cst.Rd-1.)*PQA )
    ZTHVS=PTG/PEXNS*( 1.+(cst.Rv/cst.Rd-1.)*PQS )
    ZVMOD = WIND_THRESHOLD(PVMOD,PUREF,XCISMIN, XVMODMIN, LALDTHRES)
    PRI = cst.g * PDIRCOSZW * PUREF * PUREF              \
              * (ZTHVA-ZTHVS) / (0.5 * (ZTHVA+ZTHVS) )  \
              / (ZVMOD*ZVMOD) /PZREF
    PRI = min(PRI,XRIMAX)
    return PRI


def PSifCTUW(PZL):
    if PZL<0.:
        ZX   = math.pow(1.0 - 15. * PZL, 0.25)         # Kansas unstable
        ZPSIK= 2.0 * math.log((1.0+ZX       )/2.0) \
                 +       math.log((1.0+ZX*ZX)/2.0) \
                 - 2.0 * math.atan(ZX) \
                 + 2.0 * math.atan(1.0)
        #
        ZY= (1.0 - 10.15 * PZL)**0.3333     # Convective
        ZPSIC= 1.5 * math.log((ZY*ZY+ZY+1.)/3.) \
                 - (3.0**0.5) * math.atan((2.0*ZY+1.0)/(3.0**0.5)) \
                 + 4.0        * math.atan(1.0)/(3.0**0.5)
        #
        ZF=PZL * PZL / (1.0+PZL*PZL)
        #
        PSIFCTUW=(1.-ZF) * ZPSIK + ZF * ZPSIC
    else:
        ZC=min(50.,0.35*PZL)           # Stable
        PSIFCTUW=-((1.+1.*PZL) + 0.6667*(PZL-14.28)/math.exp(ZC) + 8.525)
    return PSIFCTUW

def PSifCTTW(PZL):
    if PZL<0.:
        ZX   = (1. - 15. * PZL)**.5         # Kansas unstable
        ZPSIK= 2.0 * math.log((1.0+ZX       )/2.0)
        #
        ZY   = math.pow(1.0 - 34.15 * PZL,0.3333)  # Convective
        ZPSIC= 1.5 * math.log((ZY*ZY+ZY+1.0)/3.) \
                 - (3.0**0.5) * math.atan((2.0*ZY+1.0)/(math.pow(3.0,0.5))) \
                 + 4.0        * math.atan(1.0)/(3.0**0.5)
        #
        ZF   = PZL * PZL / (1.0+PZL*PZL)
        #
        PSIFCTTW= (1.-ZF) * ZPSIK + ZF * ZPSIC
    else:
        ZC=min(50.,0.35*PZL)           # Stable
        PSIFCTTW=-(math.pow(1.+2.*PZL/3.,1.5) + 0.6667*(PZL-14.28)/math.exp(ZC) + 8.525)
    return PSIFCTTW

def WASP_FLUX(PTA,PQA, PEXNA, PRHOA, PVMOD, PZREF, PUREF, PSST, PEXNS, PPS, PRAIN, PHS, PTP,
              XVZ0CM=0.,XRIMAX=0.2,XCISMIN=6.7E-5, XVMODMIN=0., LALDTHRES=False,CQSAT="NEW"):
    cst = Cst()

    ZRVSRDM1  = cst.Rv/cst.Rd-1. # 0.607766
    ZRDSRV    = cst.Rd/cst.Rv    # 0.62198
    ZBETAGUST = 1.2        # value based on TOGA-COARE experiment
    ZZBL      = 600.       # Set a default value for boundary layer depth
    ZS        = 10.0        # Standard heigth =10m
    ZCH10     = 0.00115
    #
    #       1.2   Array initialization by undefined values
    #
    PSFTH =0.0
    PSFTQ =0.0
    PUSTAR=0.0
    #
    PCD = 0.0
    PCDN = 0.0
    PCH = 0.0
    PCE =0.0
    PRI = 0.0
    #
    PRESA=0.0

    XSURF_EPSILON=sys.float_info.epsilon
    LPWG=False
    LCPL_WAVE=False
    LWAVEWIND=False
    LPRECIP=False
    #
    #-------------------------------------------------------------------------------
    #       2. INITIAL GUESS FOR THE ITERATIVE METHOD
    #          -------------------------------------
    #
    #       2.0     Temperature
    #
    # Set a non-zero value for the temperature gradient
    #
    if (PTA*PEXNS/PEXNA-PSST)==0.:
        ZTA=PTA-1E-3
    else:
        ZTA=PTA

    #       2.1     Wind and humidity
    #
    # Sea surface specific humidity
    #
    PQSAT=QSAT_SEAWATER(PSST,PPS,CQSAT)
    #
    # Set a minimum value to wind
    #
    ZVMOD = WIND_THRESHOLD(PVMOD,PUREF,XCISMIN, XVMODMIN, LALDTHRES)
    #
    # Specific humidity at saturation at the atm. level
    #
    ZPA = cst.P0* (math.pow(PEXNA,(cst.Cpd/cst.Rd)))
    ZQASAT = QSAT(ZTA,ZPA,CQSAT)
    #
    #
    ZO  = 0.0001
    ZWG = 0.
    if LPWG: ZWG = 0.5
    #
    ZCHARN = 0.011
    #
    #      2.2       initial guess
    #
    ZDU = ZVMOD   #wind speed difference with surface current(=0) (m/s)
                  #initial guess for gustiness factor
    ZDT = -(ZTA/PEXNA) + (PSST/PEXNS) #potential temperature difference
    ZDQ = PQSAT-PQA                         #specific humidity difference
    #
    ZDUWG = math.sqrt(math.pow(ZDU,2)+math.pow(ZWG,2))     #wind speed difference including gustiness ZWG
    #
    #      2.3   initialization of neutral coefficients
    #
    ZU10  = ZDUWG*math.log(ZS/ZO)/(PUREF/ZO)
    ZUSR  = 0.035*ZU10
    ZVISA = 1.326E-5*(1.+6.542E-3*(ZTA-cst.Tt)+8.301E-6*(ZTA-cst.Tt)**2-4.84E-9*(ZTA-cst.Tt)**3) #Andrea (1989) CRREL Rep. 89-11
    #
    ZO10 = ZCHARN*ZUSR*ZUSR/cst.g+0.11*ZVISA/ZUSR
    ZCD  = (cst.Karman/math.log(PUREF/ZO10))**2  #drag coefficient
    ZCD10= math.pow((cst.Karman/math.log(ZS/ZO10)),2)
    ZCT10= ZCH10/math.sqrt(ZCD10)
    ZOT10= ZS/math.exp(cst.Karman/ZCT10)
    #
    #-------------------------------------------------------------------------------
    #             Grachev and Fairall (JAM, 1997)
    ZCT = cst.Karman/math.log(PZREF/ZOT10)      #temperature transfer coefficient
    ZCC = cst.Karman*ZCT/ZCD               #z/L vs Rib linear coef.
    #
    ZRIBCU = -PUREF/(ZZBL*0.004*ZBETAGUST**3) #saturation or plateau Rib
    ZRIBU  = -cst.g*PUREF*(ZDT+ZRVSRDM1*ZTA*ZDQ)/(ZTA*ZDUWG**2)
    #
    if ZRIBU<0.:
        ZETU = ZCC*ZRIBU/(1.+ZRIBU/ZRIBCU)    #Unstable G and F
    else:
        ZETU = ZCC*ZRIBU/(1.+27./9.*ZRIBU/ZCC)#Stable
    #
    ZL10 = PUREF/ZETU #MO length
    #

    #
    #  First guess M-O stability dependent scaling params. (u*,T*,q*) to estimate ZO and z/L (ZZL)
    ZUSR = ZDUWG*cst.Karman/(math.log(PUREF/ZO10)-PSifCTUW(PUREF/ZL10))
    ZTSR = -ZDT*cst.Karman/(math.log(PZREF/ZOT10)-PSifCTTW(PZREF/ZL10))
    ZQSR = -ZDQ*cst.Karman/(math.log(PZREF/ZOT10)-PSifCTTW(PZREF/ZL10))
    #
    ZZL = 0.0
    #
    #
    if ZETU>50.:
        ITERMAX = 1
    else:
        ITERMAX = 3 #number of iterations
    #
    #                3.  ITERATIVE LOOP TO COMPUTE USR, TSR, QSR
    #                -------------------------------------------
    #
    if (LWAVEWIND and not LCPL_WAVE):
        ZTWAVE = 0.5*PVMOD
    else:
        ZTWAVE = PTP
    # to avoid the nullity of HS and TP
    if ZTWAVE == 0.0: ZTWAVE = 0.5*PVMOD
    if ZTWAVE > 30.0: ZTWAVE = 0.5*PVMOD
    ZCWAVE = cst.g*ZTWAVE/(2.*math.pi)
    ZWAGE = ZCWAVE/ZUSR
    #
    #
    #
    # Boucle principale
    for JLOOP in range(1, ITERMAX + 1):
        if JLOOP > ITERMAX:
            continue

        ZCHARN = CHARNOCK_WAGE(ZVMOD, ZWAGE)

        ZO = ZCHARN * math.pow(ZUSR,2) / cst.g + 0.11 * ZVISA / ZUSR  # Smith 1988

        ZOT = PZREF * math.exp(-(cst.Karman**2) / (ZCH10 * math.log(PUREF / ZO)))
        ZOQ = ZOT

        ZZL = cst.Karman * cst.g * PUREF * (ZTSR * (1. + ZRVSRDM1 * PQA) + ZRVSRDM1 * ZTA * ZQSR) / (ZTA * ZUSR**2 * (1. + ZRVSRDM1 * PQA))

        ZZTL = ZZL * PZREF / PUREF  # for T

        ZPUZ = PSifCTUW(ZZL)
        ZPTZ = PSifCTTW(ZZTL)
        ZPQZ = ZPTZ  # simplification

        # 3.1 Scale parameters
        ZUSR = ZDUWG * cst.Karman / (math.log(PUREF / ZO) - ZPUZ)
        ZTSR = -ZDT * cst.Karman / (math.log(PZREF / ZOT) - ZPTZ)
        ZQSR = -ZDQ * cst.Karman / (math.log(PZREF / ZOQ) - ZPQZ)

        # 3.2 Gustiness factor (ZWG)
        if LPWG:
            ZBF = -cst.g / ZTA * ZUSR * (ZTSR + ZRVSRDM1 * ZTA * ZQSR)
            if ZBF > 0.0:
                ZWG = ZBETAGUST * (ZBF * ZZBL)**(1. / 3.)
            else:
                ZWG = 0.2
        ZDUWG = math.sqrt(ZVMOD**2 + ZWG**2)

    # Initialisation des flux
    ZTAU = 0.0
    ZHF  = 0.0
    ZEF  = 0.0
    ZTAUR = 0.0
    ZRF   = 0.0

    # 4.1 Coefficients de transfert
    PCD = (ZUSR / ZDUWG)**2
    PCH = ZUSR * ZTSR / (ZDUWG * (ZTA * PEXNS / PEXNA - PSST))
    PCE = ZUSR * ZQSR / (ZDUWG * (PQA - PQSAT))

    PCDN = (cst.Karman / math.log(ZS / ZO))**2

    PZ0SEA = ZCHARN * ZUSR**2 / cst.g + XVZ0CM * PCD / PCDN

    ZLV = cst.LvTt + (cst.Cpv - cst.Cl) * (PSST - cst.Tt)

    # 4.2 Flux de surface
    ZTSR = -ZTSR
    ZQSR = -ZQSR
    ZTAU = -PRHOA * ZUSR**2 * ZVMOD / ZDUWG
    ZHF  = PRHOA * cst.Cpd * ZUSR * ZTSR
    ZEF  = PRHOA * ZLV  * ZUSR * ZQSR

    # 4.3 Flux de pluie
    if LPRECIP:
        ZTAC = ZTA - cst.Tt
        ZXLR = cst.LvTt + (cst.Cpv - cst.Cl) * ZTAC
        ZDQSDT = ZQASAT * ZXLR / (cst.Rd * ZTA**2)
        ZDTMP = (1.0 + 3.309e-3 * ZTAC - 1.44e-6 * ZTAC**2) * 0.02411 / (PRHOA * cst.Cpd)
        ZDWAT = 2.11e-5 * (cst.P0 / ZPA) * (ZTA / cst.Tt)**1.94

        ZALFAC = 1.0 / (1.0 + ZRDSRV * ZDQSDT * ZXLR * ZDWAT / (ZDTMP * cst.Cpd))

        ZCPLW = 4224.8482 + ZTAC * (-4.707 + ZTAC * (0.08499 + ZTAC * (1.2826e-3 + ZTAC * (4.7884e-5 - 2.0027e-6 * ZTAC))))

        ZRF = PRAIN * ZCPLW * ZALFAC * (PSST - ZTA + (PQSAT - PQA) * ZXLR / cst.Cpd)
        ZTAUR = -0.85 * (PRAIN * ZVMOD)


    # 4.5 Friction velocity (avec pluie)
    ZUSTAR2 = - (ZTAU + ZTAUR) / PRHOA
    PUSTAR = math.sqrt(ZUSTAR2)

    # 4.6 Flux totaux de surface
    PSFTH = ZHF + ZRF
    PSFTQ = ZEF / ZLV

    # 5.1 Nombre de Richardson
    ZDIRCOSZW = 1.0
    PRI = SURFACE_RI(XRIMAX, PSST, PQSAT, PEXNS, PEXNA, ZTA, ZQASAT, PZREF, PUREF, ZDIRCOSZW, PVMOD,XCISMIN, XVMODMIN, LALDTHRES)

    # 5.2 Conductance et résistance aérodynamiques
    ZAC = PCH * ZVMOD
    PRESA = 1. / max(ZAC, XSURF_EPSILON)

    # 5.3 Hauteurs de rugosité sur mer
    PZ0SEA = ZCHARN * ZUSTAR2 / cst.g + XVZ0CM * PCD / PCDN
    PZ0HSEA = PZ0SEA

    # PSFTH W m-2
    # PSFTQ kg m-2 s-1
    return PSFTH, PSFTQ, PUSTAR#, PQSAT, PCD, PCDN, PCH, PCE, PRI, PRESA, PZ0HSEA, PZ0SEA
