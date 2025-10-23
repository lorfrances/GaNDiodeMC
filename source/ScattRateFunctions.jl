function PolarRate(E_k,m_ee) # polar optical phonon scattering rate (absorption + emission)
    f_p = 1/(exp((EpO*e_c)/(kB*T_e)) - 1)
    Const = ((e_c^2*Omega_pO)/(2*pi*hbar*(2*E_k/(m_ee*m0))^0.5))*((1/(eps0*eps_inf))-(1/(eps0*eps_s)))
    AbsTerm = f_p*asinh((E_k/(EpO*e_c))^0.5)
    EmisTerm = 0.0
    if (E_k/e_c) >= EpO
        EmisTerm =  (f_p+1)*asinh(((E_k/(EpO*e_c))-1)^0.5)
    end
    return Const*AbsTerm, Const*EmisTerm
end

function AcousticRate(E_k, m_ee) # aoustic phonon scattering rate
    W_ac = (2^0.5 * (m_ee*m0)^1.5 * kB * T_e * (DefPot*e_c)^2 * (E_k)^0.5)/(pi * MatDens * VelSound^2 * hbar^4)
    return W_ac
end

function AlloyRate(E_k, m_ee, Al_x, Al) # alloy attom scattering rate (only in AlGaN)
    if Al == 1
        Pot = AlGaNPot
    else
        Pot = InGaNPot
    end
    DOS = (2^0.5*(m_ee*m0)^1.5 * E_k^0.5)/(pi^2*hbar^3)
    V0 = sqrt(3)*Lat_a^2*Lat_c/2
    W_all = (3*pi^3/(8*hbar))* V0 * (Pot*e_c)^2 *DOS * Al_x * (1 - Al_x)
    return W_all
end

function IonImpurityRate(E_k, m_ee) # ionized impurity scattering rate (Brooks-Herring model with Ridley modification)
    lambda = sqrt((DopDens_n*e_c^2)/(eps0*eps_s*kB*T_e))
    kWvk = sqrt((2* m_ee * m0 * E_k)/hbar^2)
    prefactor = (sqrt(2) * DopDens_n * e_c^4 * (m_ee*m0)^0.5)/(8 * (eps_s*eps0)^2 * hbar^2 * lambda^2)
    W_ii_BH = prefactor*(1/E_k^0.5) * (1/(1 + (lambda/(2*kWvk))^2))
    
    #  Ridley modification
    kVel = hbar*kWvk/(m_ee*m0)
    a = (2*pi*DopDens_n)^(-1/3)
    W_ii_R = (kVel/a)*(1-exp(-a*W_ii_BH/kVel))
    return W_ii_R
end

function ScatteringTable(m_eeBin, AlContent, InContent) # construct the complete scattering taple
    ScattTable = zeros(NumBin, NumEnergyLevel, NumScatMech)
    for i in 1:NumBin
        m_ee = m_eeBin[i]
        Al_x = AlContent[i]
        In_y = InContent[i]
        for j in 1:NumEnergyLevel
            E_k = j*EnergyStep*e_c
            POPRates = PolarRate(E_k, m_ee)
            AcouRate = AcousticRate(E_k, m_ee)
            AllRate = AlloyRate(E_k, m_ee, Al_x, 1)
            InlRate = AlloyRate(E_k, m_ee, In_y, 0)
            IIRate = IonImpurityRate(E_k, m_ee)
            ScattTable[i,j,1] = POPRates[1]
            ScattTable[i,j,2] = POPRates[2]
            ScattTable[i,j,3] = AcouRate
            ScattTable[i,j,4] = AllRate
            ScattTable[i,j,5] = InlRate
            ScattTable[i,j,6] = IIRate
        end
    end
    TotMaxBin = maximum(ScattTable[:,:,1]) + maximum(ScattTable[:,:,2]) + maximum(ScattTable[:,:,3]) + maximum(ScattTable[:,:,4]) + maximum(ScattTable[:,:,5]) + maximum(ScattTable[:,:,6])
    TauMax = 1/TotMaxBin
    return ScattTable, TauMax
end
