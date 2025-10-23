function InstantaneousProperties(ParPos, ParWvk, LocNode, m_eeBin, ParCharge, ParStat, WeightedCharge, BinDist, VelDist, EneDist, ParEneDist, ParVelDist, CurrentSpectrum)
    WeightedCharge .= 0.0
    BinDist .= 0
    VelDist .= 0.0
    EneDist .= 0.0
    ParEneDist .= 0.0
    ParVelDist .= 0.0
    CurrentSpectrum .= 0.0
    for i in 1:NumPar
        if ParStat[i] == 0
            # fetch individual partilce properties
            Posx = ParPos[i,1]
            Bin = BinID(Posx)
            AlloyBandTerm_i, m_eePos = ReturnAlloyBandTerm(ParPos[i,:], m_eeBin)                 # alloy band term for governing equation and effective mass
            Velx = (hbar*ParWvk[i,1])/(m0*m_eePos)
            Ene = (hbar^2*norm(ParWvk[i,:])^2)/(2*m0*m_eePos*e_c)
            BinDist[Bin] = BinDist[Bin] + 1
            VelDist[Bin] = VelDist[Bin] + Velx
            EneDist[Bin] = EneDist[Bin] + Ene

            # statistical property update
            EnInx = Int.(round((Ene/EnergyStep)))
            if EnInx <= 1
                EnInx = 1
            elseif EnInx > NumEnergyLevel
                EnInx = NumEnergyLevel
            end
            VelInx = Int.(round((Velx/VelStep))) + NumVelLevel + 1
            if Velx > MaxVel
                VelInx =  NumVelLevel*2
            elseif Velx < -MaxVel
                VelInx = 1
            end
            ParEneDist[EnInx, Bin] = ParEneDist[EnInx, Bin] + 1
            ParVelDist[VelInx, Bin] = ParVelDist[VelInx, Bin] + 1
            CurrentSpectrum[EnInx, Bin] = CurrentSpectrum[EnInx, Bin] + Velx
            # assign charge to the bins
            if PoiScheme == 2
                if Posx > LocNode[Bin] # right
                    NextNode = Bin + 1 # the next nearest node is the next one on the right
                    if NextNode > NumBin
                        NextNode = Bin # if we are at the last section then there is no next nearest bin
                    else
                        NextNode = NextNode
                    end
                elseif Posx < LocNode[Bin]
                    NextNode = Bin - 1
                    if NextNode < 1
                        NextNode = 1
                    else
                        NextNode = NextNode
                    end
                else
                    NextNode = Bin
                end
                WeightN = 1 - (abs(Posx - LocNode[Bin])/BinLen)
                WeightNN = 1 - WeightN
                WeightedCharge[Bin] = WeightedCharge[Bin] + WeightN*ParCharge[Bin]
                WeightedCharge[NextNode] = WeightedCharge[NextNode] + WeightNN*ParCharge[NextNode]
            elseif PoiScheme == 1
                WeightedCharge[Bin] = WeightedCharge[Bin] + 1*ParCharge[Bin]
            end
        end
    end
    # current calculation
    CurrentDist = ParCharge .* VelDist .* e_c
    CurrentStep = mean(CurrentDist[85:end])
    CurrentSpectrum = CurrentSpectrum .* e_c .* ParCharge'
    return CurrentDist, CurrentStep, BinDist, WeightedCharge, VelDist, EneDist, ParEneDist, ParVelDist, CurrentSpectrum
end

function UpdateAverages(CurrentDist, CurrentStep, BinDist, VelDist, EneDist, TotCurrentBin, TotCurrentStep, TotBinCount, TotVelDist, TotEneDist, ParEneDist, ParVelDist, CurrentSpectrum, TotParEneDist, TotParVelDist, TotCurrentSpectra, NodeField, NodePotential, TotNodeField, TotNodePotential, step, StepCross, TotCrossFreq)
    TotCrossFreq = TotCrossFreq + StepCross
    TotCurrentStep = TotCurrentStep + CurrentStep
    TotCurrentBin = TotCurrentBin .+ CurrentDist
    TotBinCount = TotBinCount .+ BinDist
    TotVelDist = TotVelDist .+ VelDist
    TotEneDist = TotEneDist .+ EneDist
    TotNodeField = TotNodeField .+ NodeField
    TotNodePotential = TotNodePotential .+ NodePotential

    if step >= stepStatisticsBegin
        TotParEneDist = TotParEneDist .+ ParEneDist
        TotParVelDist = TotParVelDist .+ ParVelDist
        TotCurrentSpectra = TotCurrentSpectra .+ CurrentSpectrum
    end
    return TotCurrentStep, TotCurrentBin, TotBinCount, TotVelDist, TotEneDist, TotParEneDist, TotParVelDist, TotCurrentSpectra, TotNodeField, TotNodePotential, TotCrossFreq
end

function WriteToFileMain(TotCurrentStep, TotCurrentBin, TotBinCount, TotVelDist, TotEneDist, ScattTableData, PhononCountData, AllTimeScattCount, PhononScattCount, FieldTimeData, PotentialTimeData, CurrentDataAve, CurrentDataBinAve, BinCountAve, VelDistAve, EnDistAve, TotNodeField, TotNodePotential, LocNode, ParEneDistData, ParVelDistData, TotParEneDist, TotParVelDist, TotCurrentSpectra, CurrentSpectra, step, TotCrossFreq, CurrentCross)
    DeviceCurrent = TotCurrentStep/FreqWriteOut
    DeviceCurrentBin = TotCurrentBin ./ FreqWriteOut
    AveBinDist = TotBinCount ./ FreqWriteOut
    AveVelDist = TotVelDist ./ FreqWriteOut ./ AveBinDist
    AveEneDist = TotEneDist ./ FreqWriteOut ./ AveBinDist
    AveNodeField = TotNodeField ./ FreqWriteOut
    AveNodePotential = TotNodePotential ./ FreqWriteOut
    CurrentByCross = e_c*TotCrossFreq/(FreqWriteOut*Dt)
    CurrentDensityByCross = CurrentByCross*1e-4/A_eff
    PrintCrossData = [TotCrossFreq CurrentByCross CurrentDensityByCross]

    # scattering outputs
    NetPhononEmissionBin = (DopDens_n*(PhononScattCount[:,2] - PhononScattCount[:,1])*e_c*EpO)./(ParnTypeBin*WriteOutTime)
    

    # write all averages to file
    open(CurrentDataAve,"a") do io
        writedlm(io, (DeviceCurrent/1e4)', ' ')
    end

    open(CurrentDataBinAve,"a") do io
        writedlm(io, (DeviceCurrentBin./1e4)', ' ')
    end

    open(BinCountAve,"a") do io
        writedlm(io, (AveBinDist)', ' ')
    end

    open(VelDistAve,"a") do io
        writedlm(io, (AveVelDist)', ' ')
    end   
    
    open(EnDistAve,"a") do io
        writedlm(io, (AveEneDist)', ' ')
    end   

    open(FieldTimeData,"a") do io
        writedlm(io, (AveNodeField)', ' ')
    end   

    open(PotentialTimeData,"a") do io
        writedlm(io, (AveNodePotential)', ' ')
    end   

    open(PhononCountData, "a") do io
        writedlm(io, NetPhononEmissionBin')
    end

    open(CurrentCross, "a") do io
        writedlm(io, PrintCrossData)
    end

    writedlm(ScattTableData, AllTimeScattCount)

    # write the band diagram (a bit slow...can be turned off)
    BandProfile(AveNodeField, AveNodePotential, LocNode)

    # write statistical datas
    if step >= stepStatisticsBegin
        AveParEneDist = TotParEneDist ./ (step - stepStatisticsBegin)
        AveParVelDist = TotParVelDist ./ (step - stepStatisticsBegin)
        AveCurrentSpectra = TotCurrentSpectra ./ (step-stepStatisticsBegin) ./ 1e4
        writedlm(ParEneDistData, AveParEneDist)
        writedlm(ParVelDistData, AveParVelDist)
        writedlm(CurrentSpectra, AveCurrentSpectra)
    end
    # since we wrote to file, we reset all averages
    TotCurrentStep = 0.0
    TotCurrentBin = zeros(NumBin)
    TotBinCount = zeros(NumBin)
    TotVelDist = zeros(NumBin)
    TotEneDist = zeros(NumBin)
    TotNodeField = zeros(NumBin)
    TotNodePotential = zeros(NumBin)
    PhononScattCount = zeros(NumBin, 2)
    TotCrossFreq = 0
    return TotCurrentStep, TotCurrentBin, TotBinCount, TotVelDist, TotEneDist, PhononScattCount, TotNodeField, TotNodePotential, TotCrossFreq
end


function BinID(Posx)
    # assigns particle to its bin 1D x-direction
    Bin = 1
    for i in 1:NumBin
        BinL = BinLen*(i-1)
        BinR = i*BinLen
        if (Posx >= BinL) && (Posx < BinR)
            Bin = i
            break
        end
    end
    return Bin
end

function PosWeightedGrid(Posx, LocNode, GridVector, Bin)
    if Posx > LocNode[Bin] # right
        NextNode = Bin + 1 # the next nearest node is the next one on the right
        if NextNode > NumBin
            NextNode = Bin # if we are at the last section then there is no next nearest bin
        else
            NextNode = NextNode
        end
    elseif Posx < LocNode[Bin]
        NextNode = Bin - 1
        if NextNode < 1
            NextNode = 1
        else
            NextNode = NextNode
        end
    else
        NextNode = Bin
    end
    WeightN = 1 - (abs(Posx - LocNode[Bin])/BinLen)
    WeightNN = 1 - WeightN
    WeightedVector = WeightN*GridVector[Bin] + WeightNN*GridVector[NextNode]
    return WeightedVector
end

function BandProfile(AveNodeField, AveNodePotential, LocNode)
    BandsWithPar = "BandsWithPar" * suff * ".dat"
    BandDiag = "BandDiag" * suff * ".dat"
    
    # RETURNS VISUALIZATION OF BANDPROFILE, HAS ITS OWN FREQUENCY if desired 
    Bands = zeros(NumBin)
    for i in 1:NumBin
        Pot_i = 0.0
        if LocNode[i] >= HBPosL && LocNode[i] < HBPosR && HB == 1
            BarBin =  BinID(LocNode[i]) - BinID(HBPosL+1e-18)
            Pot_i = Pot_i  + BarHeight - BarGrad*BarBin
        end

        SlfCons = 0.0
        if LocNode[i] <= LenX
            SlfCons =  AveNodePotential[i]
        end 
	Bands[i] = Pot_i + SlfCons
    end
    BandDiagDat = zeros(NumBin,2)
    BandDat = zeros(NumPar,2)
    for i in 1:NumPar
        Posx = ParPos[i,1]
        Bin = BinID(Posx)
        Enek = ParEnergy[i]/e_c
        BandDat[i,1] = Posx
        BarBin =  BinID(Posx) - BinID(HBPosL+1e-18)
        PotentialHeight =  # potential on conduction band
        if Posx >= HBPosL && Posx < HBPosR && HB == 1
            Pot_i = PosWeightedGrid(Posx, LocNode, AveNodePotential, Bin) + BarHeight - BarGrad*BarBin
        else
            Pot_i = PosWeightedGrid(Posx, LocNode, AveNodePotential, Bin)
        end
        BandDat[i,2] = Enek + Pot_i
    end
    BandDiagDat[:,1] = LocNode
    BandDiagDat[:,2] = Bands
    writedlm(BandsWithPar, BandDat)
    writedlm(BandDiag, BandDiagDat)
end   