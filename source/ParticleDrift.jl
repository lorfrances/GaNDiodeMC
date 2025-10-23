function Drift(DriftTime, initWvk, initPos, NodeField, Wvk, Pos, Ene, LocNode, initBin, m_eeBin, Vel_i, Pos_i, Wvk_i, Vel_f, EField_i, EField_f, AddAlloyBandTerm_i, AddAlloyBandTerm_f)
    # the ith electron drifts for drift time
    RemainDriftTime = DriftTime
    checkInt = 0
    Bin_f = 0
    CrossCurr = 0
    while RemainDriftTime > 0.0 # the drift is split into parts if interface is crossed
        AddField = IdentifyFieldHB(initPos) # HB local field in graded region
        
        # RK2 first approximation
        AddFieldPoi_i = IdentifyFieldPoi(initBin, initPos, NodeField, LocNode) # Poisson field at initial position
        EField_i = [AddField + AddFieldPoi_i, 0.0, 0.0] # total field at initial possition
        AlloyBandTerm_i, m_eePos = ReturnAlloyBandTerm(initPos, m_eeBin)  # alloy band term for governing equation and effective mass
        AddAlloyBandTerm_i = [AlloyBandTerm_i, 0.0, 0.0] # initial possiton alloy band term
        
        
        initK2 = initWvk[1]^2 + initWvk[2]^2 + initWvk[3]^2 # k^2 wavevector
        Vel_i .= (hbar.*initWvk)./(m0*m_eePos)                                   # intiial velocity approximation 
        Pos_i .= initPos .+ Vel_i.*RemainDriftTime                             # positiona after drift - first apporximation
        Wvk_i .= initWvk .+ RemainDriftTime.*((-e_c.*EField_i./hbar) .- (0.5.*hbar.*initK2 .*AddAlloyBandTerm_i)) # first wavevector approximation
        
        ParK2_i = Wvk_i[1]^2 + Wvk_i[2]^2 + Wvk_i[3]^2 # k^2 after 1st approximation
        finBin = BinID(Pos_i[1])                                                # bin based on first approcimation of position after the drift 
        AddFieldPoi_f = IdentifyFieldPoi(finBin, Pos_i, NodeField, LocNode)     # Poisson field at the first position approximation after first drift
        AlloyBandTerm_f, m_eePosf = ReturnAlloyBandTerm(Pos_i, m_eeBin)
        AddAlloyBandTerm_f = [AlloyBandTerm_f, 0.0, 0.0]
        
        EField_f = [EField_i[1] - AddFieldPoi_i + AddFieldPoi_f, 0.0, 0.0]                 # total electric field at the first position approximation after drift   # total electric field at the first position approximation after d
        Vel_f .= (hbar.*Wvk_i)./(m0*m_eePosf)                                       # velocity at final position based on first approximation of wavevector

        if (m_eePos == m_ee0 && m_eePosf != m_ee0) || (m_eePos != m_ee0 && m_eePosf == m_ee0)
            Vel_f .= (hbar.*Wvk_i)./(m0*m_eePos)
        end  

        # second approximation
        Pos .= initPos .+ 0.5.*RemainDriftTime.*(Vel_i .+ Vel_f)
        Wvk .= initWvk .+ RemainDriftTime.*((-e_c.*0.5.*(EField_i .+ EField_f)./hbar) .- (0.5.*hbar.*0.5*(ParK2_i + initK2) .* 0.5.*(AddAlloyBandTerm_i .+ AddAlloyBandTerm_f)))
        

        Bin_f = BinID(Pos[1])
        AlloyBandTerm_ff, m_eePos_ff = ReturnAlloyBandTerm(Pos, m_eeBin)
        ParK2 = Wvk[1]^2 + Wvk[2]^2 + Wvk[3]^2
        Ene = (ParK2 * hbar^2)/(2*m0*m_eePos_ff)
        if (m_eePos == m_ee0 && m_eePos_ff != m_ee0) || (m_eePos != m_ee0 && m_eePos_ff == m_ee0)
            Ene = (ParK2 * hbar^2)/(2*m0*m_eePos)
        end 

        #################################################################
        Pos, Wvk, HBCrossCondition, CrossCurr = CrossID(initPos, initBin, Pos, Wvk, Bin_f, LocNode) # did the particle cross some interface in the process of drifting?
        if (sum(HBCrossCondition)) > 0 # if it has passed the interface we need to consider that the properties change during the crossing, so we split up the free flight to consider both sections (should be two parts as long as the timestep is small enought)            
            initWvk, initPos, initBin, TimeToInt = DriftToInterface(HBCrossCondition, Bin_f, initPos, initWvk, Pos, Wvk, Wvk_i,Pos_i, Vel_i, Vel_f, RemainDriftTime, NodeField, LocNode, m_eeBin, initBin, EField_i, EField_f, AddAlloyBandTerm_i, AddAlloyBandTerm_f)
            UsedTime = TimeToInt
            if TimeToInt > RemainDriftTime
                println("ERROR: The time to reach an interface has surpassed the drift time!")
            end
            checkInt = checkInt + 1 
        else
            UsedTime = RemainDriftTime
        end
        RemainDriftTime = RemainDriftTime - UsedTime
    end
    if checkInt > 20
        println("WARNING: Particle crossed more than 20 interface crossed in a single timestep!  ", checkInt)
    end
    return Wvk, Pos, Ene, Bin_f, CrossCurr
end


function IdentifyFieldHB(initPos)
    if (HB == 1) && (initPos[1] >= HBPosL) && (initPos[1] <= HBPosR)
        AddField = HBField
    else
        AddField = 0.0
    end
    return AddField
end

function ReturnAlloyBandTerm(initPos, m_eeBin)
    if (HB == 1) && (initPos[1] >= HBPosL) && (initPos[1] <= HBPosR)
        m_eePos = m_eeBin[BinID(initPos[1])]
        AlloyBandTerm = dm_dz/(m0*m_eePos)^2
    else
        AlloyBandTerm = 0.0
        m_eePos = m_ee0
    end

    return AlloyBandTerm, m_eePos

end

function IdentifyFieldPoi(initBin, initPos, NodeField, LocNode)
        Bin = initBin
        if PoiScheme == 2 
            Posx = initPos[1]
            WeightedField = PosWeightedGrid(Posx, LocNode, NodeField, Bin)
            AddFieldPoi = WeightedField 
        elseif PoiScheme == 1
            AddFieldPoi = NodeField[Bin]
        end         
    return AddFieldPoi
end

function CrossID(initPos, initBin, Pos, Wvk, Bin_f, LocNode)
    # cross boundary 
    CrossEnd = 0
    CrossCurr = 0
    if Pos[1] > LenX
        Pos[1] = 1.0
        Wvk[1] = Wvk[1]
        CrossEnd = 1
        CrossCurr = 1
    elseif Pos[1] < 0.0 # particle leaves through the left boundary
        Pos[1] = 1.0 # particle has left the x-domain
        Wvk[1] = Wvk[1]
        CrossEnd = 1
    end
    
    ############################################ passing HB barrier part ########################################################################
    HBCrossL_in = 0
    HBCrossL_out = 0
    HBCrossR_out = 0
    HBCrossR_in = 0

    if HB == 1 && CrossEnd != 1
        if (initPos[1] < HBPosL) && (Pos[1] > HBPosL) # crossed barrier movin riht
            HBCrossL_in = 1
        elseif (initPos[1] > HBPosL) && (Pos[1] < HBPosL) # crossed barrier moving left
            HBCrossL_out = 1
        end

        if (initPos[1] < HBPosR) && (Pos[1] > HBPosR) # crossed barrier movin riht
            HBCrossR_out = 1
        elseif (initPos[1] > HBPosR) && (Pos[1] < HBPosR) # crossed barrier moving left
            HBCrossR_in = 1
        end

    end

    HBCrossCondition = (HBCrossL_in, HBCrossL_out, HBCrossR_out, HBCrossR_in)

    ######################################### let's see if Al bin is changed ##################################################################


    return Pos, Wvk, HBCrossCondition, CrossCurr
end

function UpdateInt(HBCrossCondition, IntPos, Ene, Wvk, Pos, m_eeBin)
    if HBCrossCondition[1] == 1 # particle is going to hit the barrier but may transmit - conservation of lateral momentum
        NewE = Ene - BarHeight*e_c
        NewParK2 = (2*m_eeBin[BinID(HBPosL+1.0e-18)]*m0*NewE)/hbar^2
        #println(m_eeBin[BinID(HBPosL+1.0e-18)])
        OldParK2_yz = Wvk[2]^2 + Wvk[3]^2
        FinWvk_x2 = NewParK2 - OldParK2_yz
        if (NewE >= 0) && (FinWvk_x2 >= 0.0) 
            #particle is  transmitted over barrier
            if Wvk[1] < 0.0
                println("ERROR: Particle approaching barrier from left has negative wavevector!")
            end

            Wvk[1] = sqrt(FinWvk_x2)
            Pos[1] = HBPosL + 1.0e-18
        else
            # particle is reflected by barrier
            if Wvk[1] < 0.0
                println("ERROR: Particle approaching barrier from left has negative wavevector!")
            end

            Pos[1] = HBPosL-1.0e-18 # should not feel the field again
            Wvk[1] = -Wvk[1] # reverse the wavevector
        end
    elseif HBCrossCondition[2] == 1 # particle falls down the barrier in the other direction
        NewE = Ene + BarHeight*e_c
        NewParK2 = (2*m_ee0*m0*NewE)/hbar^2
        OldParK2_yz = Wvk[2]^2 + Wvk[3]^2
        FinWvk_x2 = NewParK2 - OldParK2_yz
        if (NewE >= 0) && (FinWvk_x2 >= 0.0) 
            #particle is  transmitted over barrier
            if Wvk[1] > 0.0
                println("ERROR: Particle approaching barrier from left has negative wavevector!")
            end

            Wvk[1] = -sqrt(FinWvk_x2)
            Pos[1] = HBPosL - 1.0e-18
        else
            println("Condition met!   ", 1000*Ene/e_c)
            # particle is reflected by effective mass difference!!!!!!
            if Wvk[1] > 0.0
                println("ERROR: Particle approaching barrier from left has negative wavevector!")
            end

            Pos[1] = HBPosL + 1.0e-18 # should not feel the field again
            Wvk[1] = -Wvk[1] # reverse the wavevector
        end
    elseif HBCrossCondition[3] == 1 
        Wvk[1] = Wvk[1]
        Pos[1] = HBPosR + 1.0e-18
    elseif HBCrossCondition[4] == 1  # we need to ensure that the position is not modified again if modified in QW code section since the interface is shared
        Wvk[1] = Wvk[1]
        Pos[1] = HBPosR - 1.0e-18
    end


    ########################################### handle condition #########################################################################
    return Wvk, Pos
end


function DriftToInterface(HBCrossCondition, Bin_f, initPos, initWvk, Pos, Wvk, Wvk_i, Pos_i, Vel_i, Vel_f, RemainDriftTime, NodeField, LocNode, m_eeBin, initBin, EField_i, EField_f, AddAlloyBandTerm_i, AddAlloyBandTerm_f)
# drift the particle up to the inerface, # first find the time until it reaches it

    
    if (HBCrossCondition[1] == 1) || (HBCrossCondition[2] == 1)
        IntPos = HBPosL
    elseif (HBCrossCondition[3] == 1) || (HBCrossCondition[4] == 1)
        IntPos = HBPosR
    end
    
    # it is possible (due to short bin size) that more than one Al interface is crossed, so we need to use the first one in the direction fo travel

    DistToInt = IntPos - initPos[1]
    if DistToInt >= 0.0 
        MovingDir = 1
    else
        MovingDir = -1
    end

    initVelo = Vel_i[1]
    finalVelo = Vel_f[1]
    Accel = (finalVelo - initVelo)/RemainDriftTime

    if Accel != 0.0
        TimeSol1 = (-sqrt(2*Accel*DistToInt + initVelo^2) - initVelo)/Accel
        TimeSol2 = (sqrt(2*Accel*DistToInt + initVelo^2) - initVelo)/Accel
        if initVelo*DistToInt > 0
            if Accel*DistToInt > 0
                if (TimeSol2 > TimeSol1) && (TimeSol2 > 0.0)
                    TimeToInt = TimeSol2
                elseif (TimeSol1 > TimeSol2) && (TimeSol1 > 0.0)
                    TimeToInt = TimeSol1
                else
                    TimeToInt = 0.0
                    println("Time error!   ",IntPos)
                end
            else
                if (TimeSol2 < TimeSol1) && (TimeSol2 > 0.0)
                    TimeToInt = TimeSol2
                elseif (TimeSol1 < TimeSol2) && (TimeSol1 > 0.0)
                    TimeToInt = TimeSol1
                else
                    TimeToInt = 0.0
                    println("Time error!   ",IntPos)
                end
            end
        else
            if Accel*DistToInt < 0.0
                TimeToInt = 0.0
                println("Time error!   ",IntPos)
            else
                if (TimeSol2 > TimeSol1) && (TimeSol2 > 0.0)
                    TimeToInt = TimeSol2
                elseif (TimeSol1 > TimeSol2) && (TimeSol1 > 0.0)
                    TimeToInt = TimeSol1
                else
                    TimeToInt = 0.0
                    println("Time error   ",IntPos)
                end
            end
        end
    else
        TimeToInt = DistToInt/initVelo
    end


    if TimeToInt < 0
        println("ERROR: The time to reach an interface is negative!")
    end

    
    AddField = IdentifyFieldHB(initPos) # HB local field

    # first approximation
    AddFieldPoi_i = IdentifyFieldPoi(initBin, initPos, NodeField, LocNode) # Poisson field at initial position
    EField_i = [AddField + AddFieldPoi_i, 0.0, 0.0]
    AlloyBandTerm_i, m_eePos = ReturnAlloyBandTerm(initPos, m_eeBin)                 # alloy band term for governing equation and effective mass
    AddAlloyBandTerm_i = [AlloyBandTerm_i, 0.0, 0.0]

    Vel_i .= (hbar.*initWvk)./(m0*m_eePos)                                   # intiial velocity approximation 
    Pos_i .= initPos .+ Vel_i.*TimeToInt
    Pos_i[1] = IntPos - MovingDir*1.0e-18 
    initK2 = initWvk[1]^2 + initWvk[2]^2 + initWvk[3]^2    # positiona after drift - first apporximation                                               # positiona after drift - first apporximation
    Wvk_i .= initWvk .+ TimeToInt.*((-e_c.*EField_i./hbar) .- (0.5.*hbar.*initK2.*AddAlloyBandTerm_i))              # wavevector after drift - first approximation
    ParK2_i = Wvk_i[1]^2 + Wvk_i[2]^2 + Wvk_i[3]^2

    finBin = BinID(Pos_i[1])                                               # bin based on first approcimation of position after the drift 
    AddFieldPoi_f = IdentifyFieldPoi(finBin, Pos_i, NodeField, LocNode)     # Poisson field at the first position approximation after first drift
    AlloyBandTerm_f, m_eePosf = ReturnAlloyBandTerm(Pos_i, m_eeBin)
    AddAlloyBandTerm_f = [AlloyBandTerm_f, 0.0, 0.0]
    EField_f = [EField_i[1] - AddFieldPoi_i + AddFieldPoi_f, 0.0, 0.0]                 # total electric field at the first position approximation after drift   # total electric field at the first position approximation after d
    Vel_f .= (hbar.*Wvk_i)./(m0*m_eePosf)                                       # velocity at final position based on first approximation of wavevector
 
    

    # second approximation
    Pos .= Pos_i
    Wvk .= initWvk .+ TimeToInt.*((-e_c.*0.5.*(EField_i .+ EField_f)./hbar) .- (0.5.*hbar.*0.5*(ParK2_i + initK2) .* 0.5.*(AddAlloyBandTerm_i .+ AddAlloyBandTerm_f)))
    ParK2 = Wvk[1]^2 + Wvk[2]^2 + Wvk[3]^2
    Ene = (ParK2 * hbar^2)/(2*m0*m_eePosf)

    # now, we consider exactly what the particle does at the interface and update the properties for continued drift
    Wvk, Pos = UpdateInt(HBCrossCondition, IntPos, Ene, Wvk, Pos, m_eeBin)
    initWvk .= Wvk # initial values for the continued drift 
    initPos .= Pos
    initBin = BinID(initPos[1])
    return initWvk, initPos, initBin, TimeToInt
end

function InjectOhmic(ParPos, ParWvk, ParStat, ParnTypeBin, BinDist)
    # function injects particle to the simulation based on maintaining a constant carrier concentration in first bin
    TargetPop = ParnTypeBin[1]
    CurrentPop = BinDist[1]
    MissingPar = TargetPop - CurrentPop
    if MissingPar > 0
        CountInj = 0
        for i in 1:NumPar
            # update the particle status 
            if ParStat[i] != 0 
                ParPos[i,1] = 0.0
                RandomNum = randopen()
                ParWvk[i,1] = 1*(-2*m_ee0*m0*kB*T_e*log(RandomNum))^0.5/hbar
                ParStat[i] = 0
                CountInj = CountInj + 1
            end
            if CountInj >= MissingPar
                break
            end
        end
    end
    return ParPos, ParWvk, ParStat
end