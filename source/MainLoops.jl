function MCMain(ParTau, ParWvk, ParPos, ParEnergy, ParStat,  ParCharge, ScattTable, LocNode, m_eeBin, TauMax, PolarCharge, ParnTypeBin, PolarDeriv)
    # start simulation counters
    time = 0.0
    CountWrite = 1
    CountWriteSub = 1
    CountPoi = 1
    step = 1

    # initialize all file systems for writing output
    ScattTableData = output_path * "ScattBin" * suff * ".dat"
    PhononCountData = output_path*"PhononCountData" * suff * ".dat"
    FieldTimeData = output_path*"FieldTimeData" * suff * ".dat"
    PotentialTimeData = output_path*"PotentialTimeData" * suff * ".dat"
    CurrentDataAve = output_path*"CurrentDataAve" * suff * ".dat"
    CurrentDataBinAve = output_path*"CurrentDataBinAve" * suff * ".dat"
    BinCountAve = output_path*"BinCountAve" * suff * ".dat"
    VelDistAve = output_path*"VelDistAve" * suff * ".dat"
    EnDistAve = output_path*"EnDistAve" * suff * ".dat"
    ParEneDistData = output_path*"ParEneDistData" * suff * ".dat"
    ParVelDistData = output_path*"ParVelDistData" * suff * ".dat"
    CurrentSpectra = output_path*"CurrentSpectra" * suff * ".dat"
    CurrentCross = output_path*"CurrentCross" * suff * ".dat"
    # initialize some system variables for averaging
    PhononScattCount = zeros(NumBin, 2)
    AllTimeScattCount = zeros(NumBin, NumScatMech)
    NodeField = zeros(NumBin)
    NodePotential = zeros(NumBin)
    TotNodeField = zeros(NumBin)
    TotNodePotential = zeros(NumBin)
    TotBinCount = zeros(NumBin)
    TotCurrentBin = zeros(NumBin)
    TotVelDist = zeros(NumBin)
    TotEneDist = zeros(NumBin)
    TotParEneDist = zeros(NumEnergyLevel, NumBin)
    TotParVelDist = zeros(NumVelLevel*2+1, NumBin)
    TotCurrentSpectra = zeros(NumEnergyLevel, NumBin)
    TotCurrentStep = 0.0
    TotCrossFreq = 0
    # recycled vectors

    WeightedCharge = zeros(NumBin)
    BinDist = zeros(NumBin)
    VelDist = zeros(NumBin)
    EneDist = zeros(NumBin)
    ParEneDist = zeros(NumEnergyLevel, NumBin)
    ParVelDist = zeros(NumVelLevel*2+1, NumBin)
    CurrentSpectrum = zeros(NumEnergyLevel, NumBin)

    # write intial file message
    open(ScattTableData, "w") do io
        println(io, "NET SCATTERING EVENTS AND TYPE in bins during write-out timesteps (rows)")
    end

    open(PhononCountData, "w") do io
        println(io,"NET PHONON EMISSION POWER (W/m^3) in bins during write-out timesteps (columns)")
    end

    open(EnDistAve, "w") do io
        println(io,"AVERAGE PARTICLE ENERGY (eV) in bins (columns)")
    end

    open(FieldTimeData, "w") do io
        println(io,"AVERAGE ELECTRIC FIELD (V/m) in bins (columns)")
    end

    open(PotentialTimeData, "w") do io
        println(io,"AVERAGE INTERNAL POTENTIAL (V) in bins (columns)")
    end

    open(CurrentDataAve, "w") do io
        println(io,"AVERAGE CURRENT (A/cm^2)")
    end

    open(CurrentDataBinAve, "w") do io
        println(io,"AVERAGE CURRENT (A/cm^2) in bins (columns)")
    end

    open(BinCountAve, "w") do io
        println(io,"AVERAGE PARTICLE COUNT (N) in bins (columns)")
    end

    open(VelDistAve, "w") do io
        println(io,"AVERAGE VELOCITY (m/s) in bins (columns)")
    end

    open(CurrentCross,"w") do io
        println(io, "CURRENT BY PARTICLES CROSSING RIGHT BOUNDARY")
    end

    while time < TotTime # begin simulation
        ParTau, ParWvk, ParPos, ParEnergy, AllTimeScattCount, PhononScattCount, InActPar, ParStat, StepCross = LoopParticles(ParTau, ParWvk, ParPos, ParEnergy, PhononScattCount, AllTimeScattCount, NodeField, ScattTable, ParStat, m_eeBin, TauMax)
        CurrentDist, CurrentStep, BinDist, WeightedCharge, VelDist, EneDist, ParEneDist, ParVelDist, CurrentSpectrum = InstantaneousProperties(ParPos, ParWvk, LocNode, m_eeBin, ParCharge, ParStat, WeightedCharge, BinDist, VelDist, EneDist, ParEneDist, ParVelDist, CurrentSpectrum)

        TotCurrentStep, TotCurrentBin, TotBinCount, TotVelDist, TotEneDist, TotParEneDist, TotParVelDist, TotCurrentSpectra, TotNodeField, TotNodePotential, TotCrossFreq = UpdateAverages(CurrentDist, CurrentStep, BinDist, VelDist, EneDist, TotCurrentBin, TotCurrentStep, TotBinCount, TotVelDist, TotEneDist, ParEneDist, ParVelDist, CurrentSpectrum, TotParEneDist, TotParVelDist, TotCurrentSpectra, NodeField, NodePotential, TotNodeField, TotNodePotential, step, StepCross, TotCrossFreq)

        ParPos, ParWvk, ParStat = InjectOhmic(ParPos, ParWvk, ParStat, ParnTypeBin, BinDist)

        if (Poi == 1 && time > timeBeginPoisson && CountPoi >= PoiUpdateFreq)
            ElCharge, ElDop = ChargeDopDist(LocNode, WeightedCharge)
            NodeField, NodePotential = PoissonSolver(ElCharge, ElDop, PolarDeriv)
            CountPoi = 0
        elseif (Poi == 0 || time <= timeBeginPoisson)
            NodeField, NodePotential = FixedPotentialFieldSchottky(LocNode)
        end

        if CountWrite >= FreqWriteOut
            println("TIME: ", time, "  INACTIVE PARTICLES: ", InActPar)
            TotCurrentStep, TotCurrentBin, TotBinCount, TotVelDist, TotEneDist, PhononScattCount, TotNodeField, TotNodePotential, TotCrossFreq =  WriteToFileMain(TotCurrentStep, TotCurrentBin, TotBinCount, TotVelDist, TotEneDist, ScattTableData, PhononCountData, AllTimeScattCount, PhononScattCount, FieldTimeData, PotentialTimeData, CurrentDataAve, CurrentDataBinAve, BinCountAve, VelDistAve, EnDistAve, TotNodeField, TotNodePotential, LocNode, ParEneDistData, ParVelDistData, TotParEneDist, TotParVelDist, TotCurrentSpectra, CurrentSpectra, step, TotCrossFreq, CurrentCross)
            CountWrite = 0
        end


        CountWrite = CountWrite + 1
        CountWriteSub = CountWriteSub + 1
        CountPoi = CountPoi + 1
        step = step + 1
        time = time + Dt
    end

end

function LoopParticles(ParTau, ParWvk, ParPos, ParEnergy, PhononScattCount, AllTimeScattCount, NodeField, ScattTable, ParStat, m_eeBin, TauMax) # primary MC loop where particles are drifted and scattered every timestep
    InActPar = 0
    Wvk = zeros(3)
    Pos = zeros(3)
    Vel_i = zeros(3)
    Pos_i = zeros(3)
    Wvk_i = zeros(3)
    Vel_f = zeros(3)
    EField_i = zeros(3)
    EField_f = zeros(3)
    AddAlloyBandTerm_i = zeros(3)
    AddAlloyBandTerm_f = zeros(3)
    initWvk = zeros(3)
    initPos = zeros(3)
    Ene = 0.0
    StepCross = 0
    for i in 1:NumPar
        CarryTime = 0.0
        if ParStat[i] == 0 # particle is active (not out of bounds)
            RemainTime = Dt # initialize the remaining drift time as the timestep
            CrossCurr = 0
            kk = 1
            while RemainTime > 0.0 # iterate until remaining drift time reaches 0.0
                # extract particle
                initWvk .= ParWvk[i,:] # particle properties at beiginning of loop
                initPos .= ParPos[i,:]
                initBin = BinID(initPos[1])

                if kk == 1 # first loop
                    #println("CHECKED")
                    NewTau = ParTau[i]
                else
                    NewTau = -TauMax*log(randopen())
                end

                if NewTau <= RemainTime
                    DriftTime = NewTau
                    RemainTime = RemainTime - NewTau
                    Wvk, Pos, Ene, LastBin, CrossCurr  = Drift(DriftTime, initWvk, initPos, NodeField, Wvk, Pos, Ene, LocNode, initBin, m_eeBin, Vel_i, Pos_i, Wvk_i, Vel_f, EField_i, EField_f, AddAlloyBandTerm_i, AddAlloyBandTerm_f)
                    Ene, Wvk, SelectMech = Scatter(Ene, Wvk, ScattTable, LastBin, TauMax, m_eeBin)
                    if SelectMech != 0
                        AllTimeScattCount[LastBin, SelectMech] = AllTimeScattCount[LastBin, SelectMech] + 1
                        if SelectMech == 1 || SelectMech == 2
                            PhononScattCount[LastBin, SelectMech] = PhononScattCount[LastBin, SelectMech] + 1
                        end
                    end
                else
                    DriftTime = RemainTime
                    CarryTime = NewTau - RemainTime
                    RemainTime = 0.0
                    Wvk, Pos, Ene, LastBin, CrossCurr = Drift(DriftTime, initWvk, initPos, NodeField, Wvk, Pos, Ene, LocNode, initBin, m_eeBin, Vel_i, Pos_i, Wvk_i, Vel_f, EField_i, EField_f, AddAlloyBandTerm_i, AddAlloyBandTerm_f)
                end
                ParWvk[i,:] = Wvk
                ParPos[i,:] = Pos
                ParEnergy[i] = Ene
                if Pos[1] == 1
                    ParStat[i] = 1
                    CarryTime = -TauMax*log(randopen())
                    break
                end
                kk = kk + 1
            end
            ParTau[i] = CarryTime
            StepCross = StepCross + CrossCurr  
        else
            InActPar =  InActPar + 1
        end 
    end
    #println("crossed: ",StepCross)
    return ParTau, ParWvk, ParPos, ParEnergy,  AllTimeScattCount, PhononScattCount, InActPar, ParStat, StepCross
end