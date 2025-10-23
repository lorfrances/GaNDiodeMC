include("Initialize.jl")
include("ScattRateFunctions.jl")
include("MainLoops.jl")
include("randopen.jl")
include("ParticleDrift.jl")
include("ParticleScatter.jl")
include("Properties.jl")
include("PoissonSolver.jl")
include("FindIn.jl")

using LinearAlgebra
using DelimitedFiles
using Statistics
using Printf
using Base.Threads
#using Plots

# physical Constants
const kB =  1.380649e-23
const hbar = 1.054571817e-34
const m0 = 9.1093837e-31
const e_c = 1.60217663e-19
const eps0 = 8.8541878188e-12

# read-in input files
input_path = ARGS[1]
output_path = ARGS[2]
InputData = readdlm(input_path, comments=true, comment_char='#')

# Electrical calculation input
const ApplVoltage = InputData[1]
const BarHeightSchottky = 0.5

# Domain length and bins
const LenX = InputData[2]
const NumBin = round.(Int,InputData[3])
const BinLen = 1.0e-9

# HB calculation conditions
const HB = round.(Int,InputData[5])
const HBPosL = InputData[7]
const HBPosR = InputData[8]
const BarHeight = InputData[6]
const BarHeightEnd = 0.0 # no potential gain
const CLMSetting = 0 # barrier transition conserved lateral momentum mechanism (0 - full (lax), 1 - little more restrictive)
# general calculation input
const Dt = 0.75e-15
const TotTime = 15.0e-9
const NumPar = round.(Int,InputData[4])
const FreqWriteOut = 1200 # property write out frequnecy
const PoiUpdateFreq = 2
const T_e = InputData[10]
const NumScatMech = 6
const stepStatisticsBegin = 12000

# consider aluminum contant
const ConsiderAl = 1
const ConsiderAlPolar = 1
const RatioPolar = 0.5

# Doping
const DopDens_n = InputData[9]*1e6 # dopant density electron
const DopDens_p = DopDens_n

# system discretization 
const NumEnergyLevel = 800
const MaxEnergy = 0.8
const NumVelLevel = 70
const MaxVel = 8e5

# GaN material properties
const m_he = 1.428
const eps_s = 10.0
const eps_inf = 5.5
const DefPot = 10.1 #eV
const VelSound = 6560 #m/s
const MatDens = 6150 #kg/m3
const EpO = 0.091 #eV
const Lat_a = 3.189e-10 #m
const Lat_c = 5.186e-10
const RatEcEv = 0.70 #

const m_ee0 = 0.193
const m_ee0AlN = 0.3271
const m_ee0InN = 0.046
const EgAlN = 6.30
const Eg = 3.47
const EgInN = 0.78
const AlGaNPot = 0.60 #eV # for alloy scattering
const InGaNPot = 0.60 #eV

#bandgap parameters
const bAlGaN = 1.0
const bInGaN = 1.348
const bAlInN = 3.678
const Beg = 0


const PolarChargeBulk = RatioPolar * -0.034
const Omega_pO = EpO*e_c/hbar

const AlloyFractionMethod = 0 # 0 - aluminum 100% is used for barrier height, indium to drop initial, 1 - alumnium fraction is constant, indium used to calibrate (may be insufficient for the barrier height)
const AlFractionConst = 0.10 # may be insufficient for barrier height
# Poisson calculation conditions
const Poi = 1
const PoiScheme = 2 # 1 = nearest grid, 2 = weighted grid
const timeParEneDist = 5e-12
const timeBeginPoisson = 0.0 # do not solve Poisson (use approximation) until now

# electrical initialization calculations
const DopDens_i = 2*(kB*T_e/(2*pi*hbar^2))^1.5 * (m_ee0*m0*m_he*m0)^0.75 * exp(-Eg*e_c/(2*kB*T_e))
const Nc = 2*((kB*T_e*m_ee0*m0*2*pi)/(hbar*2*pi)^2)^1.5
const FermiEnergy = (kB*T_e/e_c)*log(DopDens_n/Nc)
const LenDebye = sqrt((eps0*eps_s*kB*T_e)/(e_c^2*DopDens_n))
const BIVoltage = BarHeightSchottky + FermiEnergy# diode built-in voltage (approximation)

# solve for the depletion width in Schottky diode based on 

const DeplWidth_n = (2*eps0*eps_s*(BIVoltage-ApplVoltage)/(e_c*DopDens_n))^0.5
const DeplPos_n = LenX - DeplWidth_n
const NodePotential_n = BIVoltage - ApplVoltage

# barrier initialization
const BarLen = HBPosR - HBPosL
const HBField = -(BarHeight - BarHeightEnd)/BarLen
const BarGrad = round(BarHeight/BarLen/1e9,digits=15)
println(BarGrad)
# particle properties initialization
const NumPar_n = NumPar
const EnergyStep = MaxEnergy/NumEnergyLevel
const ParPerBin = Int(NumPar/NumBin)
const VelStep = (MaxVel)/NumVelLevel
const rRMS = 1.0e-9
const A_eff = ParPerBin/(DopDens_n*BinLen)

# call primary initialization function
LocNode, ParnTypeBin, ParCharge, m_eeBin, ScattTable, TauMax, PolarCharge, PolarDeriv, AlContent, AlloyChange = InitializeMaterial()
ParWvk, ParTau, ParPos, ParStat, ParEnergy = InitializeParticles(TauMax)
const WriteOutTime = FreqWriteOut * Dt
const dm_dz = AlloyChange
println(dm_dz)
# simulation start message
if HB == 0
    const suff = "_NB"   ############ suffix for ouput files (to distinguish whtehr barrier is included)
else
    const suff = ""
end

println("SIMULATION START", " Threads used: ", Threads.nthreads())
println("Applied voltage: ", ApplVoltage,  "  Built-in voltage: ", BIVoltage, " Approximate depletion length (n) (nm):  ", DeplWidth_n*1e9, "  Right-side BC (accounting for polarizations) ", NodePotential_n )
println("HB:  ", HB, "  HB field:  ", HBField, "   Maximum Al_x fraction: ", maximum(AlContent), "  Al Polar considered?  ", ConsiderAlPolar)
println("HB start pos (m): ", HBPosL, " HB end pos (m): ", HBPosR, " HB height (eV): ", BarHeight)
println("Dopant Density:  ", DopDens_n, "  Intrinsic doping:  ", DopDens_i)
println("Poisson: ", Poi)
println("Particle free flight time (maximum): ", TauMax)
println("initial particles - n type: ", NumPar_n, "  particles per n-type bin: ", ParnTypeBin)

# enter transient loop
MCMain(ParTau, ParWvk, ParPos, ParEnergy, ParStat, ParCharge, ScattTable, LocNode, m_eeBin, TauMax, PolarCharge, ParnTypeBin, PolarDeriv)