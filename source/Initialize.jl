function InitializeMaterial()
    LocNode = zeros(NumBin)
    for i in 1:NumBin
        BinCenter = (i*BinLen) - 0.5*BinLen   
        LocNode[i] = round(BinCenter,digits=15)
    end
    ParnTypeBin = zeros(NumBin) .+ NumPar_n/NumBin
    ParnTypeBin[end] = NumPar_n/NumBin/2
    ParCharge = DopDens_n./ParnTypeBin
    # material properties, we consider the Al heterobarrier if that feature is turned on
    m_eeBin = zeros(NumBin) .+ m_ee0
    AlContent = zeros(NumBin)
    InContent = zeros(NumBin)
    PotHeight = zeros(NumBin)
    PolarCharge = zeros(NumBin) .+ PolarChargeBulk
    PolarDeriv = zeros(NumBin)
    
    if HB == 1 && AlloyFractionMethod == 0
        NewEg = Eg + BarHeight*(1/RatEcEv)
        println(NewEg)
        AlFraction = FindAl(NewEg)
    elseif HB == 1 && AlloyFractionMethod == 1
        AlFraction = AlFractionConst
    end

    for i in 1:NumBin
        if (ConsiderAl == 1 && HB == 1) #simulate HB and aluminum
            if (LocNode[i] >= HBPosL) && (LocNode[i] <= HBPosR)
                #AlFraction =  ((Eg + ((BarHeight/RatEcEv) + BarHeight)) - Eg)/(EgAlN-Eg)
                #println(AlFraction)# 10 percent aluminum
                BarBin =  BinID(LocNode[i]) - BinID(HBPosL+1e-18)
                PotentialHeight = BarHeight - BarGrad*BarBin # potential on conduction band
                NewEg = Eg + PotentialHeight*(1/RatEcEv) # first term is correposnding cahnge in valence band so tot bandgap change is added to above line
               # println(NewEg)
                InFraction = FindIn(NewEg, AlFraction)
                if InFraction <= 0.0
                    InFraction = 0.0
                end
                y = 1-AlFraction-InFraction
                x = AlFraction
                m_eeBin[i] = x*m_ee0AlN + y*m_ee0 + (1-x-y)*m_ee0InN 
                AlContent[i] = AlFraction
                InContent[i] = InFraction
                PolarCharge[i] = RatioPolar * (-0.090*AlFraction - 0.034*(1-AlFraction) + 0.038*AlFraction*(1-AlFraction))
                PotHeight[i] = PotentialHeight
            end
        end
    end
    BinHBEnd = BinID(HBPosR+1.0e-18)
    m_eeBin[BinHBEnd:end] .= m_eeBin[BinHBEnd-1]
    AlContent[BinHBEnd:end] .= AlContent[BinHBEnd-1]
    InContent[BinHBEnd:end] .= InContent[BinHBEnd-1]
    PolarCharge[BinHBEnd:end] .= PolarCharge[BinHBEnd-1]
    println(m_eeBin)
    println(InContent)
    println(AlContent)
    #println(AlContent)
    AlloyChange = m0*(m_eeBin[BinID(HBPosL+1.0e-18)] - m_eeBin[BinID(HBPosR-1.0e-18)])/(HBPosR-HBPosL)
    if ConsiderAlPolar == 1
        for i in 2:NumBin
            PolarDeriv[i] = (PolarCharge[i]-PolarCharge[i-1])/(BinLen)
        end
    end
    # now return the bin-dependent scattering table
    ScattTable, TauMax = ScatteringTable(m_eeBin, AlContent, InContent)
    writedlm("Scattering_Table.dat", ScattTable[64,:,:])
    return LocNode, ParnTypeBin, ParCharge, m_eeBin, ScattTable, TauMax, PolarCharge, PolarDeriv, AlContent, AlloyChange
end

function InitializeParticles(TauMax)
    dn_dE = zeros(NumEnergyLevel)
    CumEnDist = zeros(NumEnergyLevel)

    # calculate electron occupation
    for i in 1:NumEnergyLevel
        E_k = i*EnergyStep*e_c
        DOS = (m_ee0*m0*(2*m_ee0*m0*E_k)^0.5)/(pi^2*hbar^3)
        f_e = 1/(exp((E_k - FermiEnergy*e_c)/(kB*T_e)) + 1)
        dn_dE[i] = DOS*f_e  
    end

    CumEn = 0.0
    # cumulative energy distribution
    for i in 1:NumEnergyLevel 
        CumEn = dn_dE[i] + CumEn
        CumEnDist[i] = CumEn/sum(dn_dE)
    end

    # Assign initial particle energy and positions and free flight
    ParEnergy = zeros(NumPar)
    ParStat = zeros(NumPar)
    ParTau = zeros(NumPar)
    ParPos = zeros(NumPar,3)
    ParWvk = zeros(NumPar,3)
    TotPar = 0
    for i in 1:NumBin
        for j in 1:ParPerBin
            inx = (i-1)*ParPerBin + j
            RandomNum = randopen()
            NearestEn = argmin(abs.(CumEnDist .- RandomNum))
            ParEnergy[inx] = NearestEn*EnergyStep*e_c
            ParWvkk2 = (2 * (m_eeBin[i]*m0) * (EnergyStep*NearestEn*e_c))/hbar^2
            # initial wavevector, from Seungha (Vasileska) randomized
            fai = 2*pi*RandomNum
            RandomNum = randopen()
            ct = 1 - 2*RandomNum
            st = sqrt(1-ct^2)
            kx = sqrt(ParWvkk2)*st*cos(fai)
            ky = sqrt(ParWvkk2)*st*sin(fai)
            kz = sqrt(ParWvkk2)*ct
            ParWvk[inx,1] = kx
            ParWvk[inx,2] = ky
            ParWvk[inx,3] = kz
            ParPos[inx,1] = randopen()*LenX
            ParTau[inx] = -TauMax*log(RandomNum)
            TotPar = TotPar + 1
        end
    end
    return ParWvk, ParTau, ParPos, ParStat, ParEnergy
end