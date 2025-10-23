function FixedPotentialFieldSchottky(LocNode)
    NodeField = zeros(NumBin)
    NodePotential = zeros(NumBin)
    for i in 1:NumBin
        if (LocNode[i] >= DeplPos_n)
            NodeField[i] = e_c*DopDens_n*(LocNode[i]-DeplPos_n)/(eps0*eps_s)
        else
            NodeField[i] = 0.0
        end
    end
    return NodeField, NodePotential
end

function ChargeDopDist(LocNode, WeightedCharge)
    ElDop = zeros(NumBin)
    ElCharge = zeros(NumBin)
    for i in 1:NumBin
        ElDop[i] = DopDens_n
        ElCharge[i] = WeightedCharge[i]
    end
    return ElCharge, ElDop
end

function PoissonSolver(ElCharge, ElDop, PolarDeriv)
    NodeField = zeros(NumBin)
    b = -(e_c*(ElDop .- ElCharge).-PolarDeriv).*(BinLen^2)./(eps0*eps_s)
    b[1] = 0.0 # zero potential on left side
    b[NumBin] = NodePotential_n # total potential on right side
    MainDiag = 2 * ones(NumBin)
    OffDiag = -1 * ones(NumBin-1)
    A = Matrix(Tridiagonal(OffDiag, MainDiag, OffDiag))
    A[1,1] = 1
    A[1,2] = 0
    A[NumBin,NumBin] = 1
    A[NumBin,NumBin-1] = 0
    NodePotential = A\b   
    # evaluate the derivative to determine the electric field. Use a central differnec for internal points, and forward or backward differences for the edges
    for i in 2:NumBin-1
        NodeField[i] = (NodePotential[i+1] - NodePotential[i-1])/(2*BinLen)
    end
    NodeField[1] = (NodePotential[2] - NodePotential[1])/BinLen
    NodeField[NumBin] = (NodePotential[NumBin] - NodePotential[NumBin-1])/BinLen
    return NodeField, NodePotential
end