######################################################### SCATTERING FUNCTIONS#######################

function Scatter(Ene, Wvk, ScattTable, ScatterBin, TauMax, m_eeBin)
    # this function identifies the relevant scattering mechanism and calls functions to update their energy and momentum
    ParEn = Ene/e_c
    EnInx = floor(Int, ParEn/EnergyStep) 
    frac = (ParEn - EnInx*EnergyStep)/EnergyStep
    if EnInx == 0
       EnInx = 1
       frac = 0
    elseif EnInx >= NumEnergyLevel
       EnInx = NumEnergyLevel-1
       frac = 1
    end
    GammaInx = randopen()*(1/TauMax)
    SelectMech = 0
    RateLo = 0.0
    RateHi = 0.0
    for j in 1:NumScatMech
        RateLo = RateHi
        RateHi = (1-frac)*ScattTable[ScatterBin, EnInx, j] + frac*ScattTable[ScatterBin, EnInx+1,j] + RateHi
        if (GammaInx >= RateLo) && (GammaInx < RateHi)
            SelectMech = j
            break
        end
    end
    if SelectMech == 1 || SelectMech == 2
        Ene, Wvk = PolarOpticalUpdate(SelectMech, ParEn, Wvk, ScatterBin, m_eeBin)
    elseif (SelectMech == 3) || (SelectMech == 4) || (SelectMech == 5)
        Wvk = AcouPhononUpdate(Wvk) 
    elseif SelectMech == 6
        Wvk = IonImpurityUpdate(Wvk, ParEn, ScatterBin, m_eeBin)
    else #SelectMech == 0
        # self-scatter
        Ene = Ene
        Wvk = Wvk
    end
    return Ene, Wvk, SelectMech
end


##########################################
function PolarOpticalUpdate(SelectMech, ParEn, Wvk, ScatterBin, m_eeBin)
    # add energy to the particle
    if SelectMech == 1
        NewEn = (ParEn + EpO)*e_c
    else
        NewEn = (ParEn - EpO)*e_c
    end
    # adjust particle momentum using Vasileska's approach
    Wvk_xy = sqrt(Wvk[1]^2 + Wvk[2]^2)
    Wvk_xy_z = sqrt(Wvk_xy^2 + Wvk[3]^2)
    costh0 = Wvk[3]/Wvk_xy_z
    sinth0 = Wvk_xy/Wvk_xy_z
    cosph0 = Wvk[1]/Wvk_xy
    sinph0 = Wvk[2]/Wvk_xy
    
    Wvk_p = sqrt(((2 * (m_eeBin[ScatterBin]*m0) * NewEn))/hbar^2)
    alpha = 0.0
    ge = ParEn*e_c*(1 + alpha*ParEn*e_c)
    gnew = NewEn*(1 + alpha*NewEn)
    zeta = 2*sqrt(ge*gnew)/(ge+gnew-2*sqrt(ge*gnew))
    RandomNum = randopen()
    costh = ((zeta+1)-(2*zeta+1)^RandomNum)/zeta
    sinth = sqrt(1-costh^2)
    fai = 2*pi*randopen()
    cosfai = cos(fai)
    sinfai = sin(fai)
    kxp = Wvk_p*sinth*cosfai
    kyp = Wvk_p*sinth*sinfai
    kzp = Wvk_p*costh
    
    
    Wvk[1] = kxp*cosph0*costh0-kyp*sinph0+kzp*cosph0*sinth0
    Wvk[2] = kxp*sinph0*costh0+kyp*cosph0+kzp*sinph0*sinth0
    Wvk[3] = -kxp*sinth0+kzp*costh0
    Ene = NewEn
    return Ene, Wvk
end

function AcouPhononUpdate(Wvk)
    # completely trandomized process
    Wvk_p = norm(Wvk)
    fi = 2*pi*randopen()
    costh = 1-2*randopen()
    sinth = sqrt(1 - costh^2)
    Wvk[1] = Wvk_p*sinth*cos(fi)
    Wvk[2] = Wvk_p*sinth*sin(fi)
    Wvk[3] = Wvk_p*costh
    return Wvk
end

function IonImpurityUpdate(Wvk, ParEn, ScatterBin, m_eeBin)
    # strong forward angle preference, based on Seungha/Vasileska code
    kxy = sqrt(Wvk[1]^2 + Wvk[2]^2)
    kp = norm(Wvk)
    cth0 = Wvk[3]/kp
    sth0 = kxy/kp
    cfi0 = Wvk[1]/kxy
    sfi0 = Wvk[2]/kxy
   # kp = sqrt(((2 * (m_eeBin[ScatterBin]*m0) * ParEn*e_c))/hbar^2)
    RandNum = randopen()
    cth = 1 - 2*RandNum/(1 + 4*(1-RandNum)*kp^2*LenDebye^2)
    sth = sqrt(1 - cth^2)
    fi = 2*pi*randopen()
    cfi = cos(fi)
    sfi = sin(fi)
    kxp = kp*sth*cfi
    kyp = kp*sth*sfi
    kzp = kp*cth
    Wvk[1] = kxp*cfi0*cth0-kyp*sfi0+kzp*cfi0*sth0
    Wvk[2] = kxp*sfi0*cth0+kyp*cfi0+kzp*sfi0*sth0
    Wvk[3] = -kxp*sth0+kzp*cth0
    return Wvk
end