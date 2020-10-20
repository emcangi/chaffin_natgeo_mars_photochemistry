################################################################################
# PARAMETERS.jl
# TYPE: (1) Model files - required
# DESCRIPTION: Just some standard global constants that need to get used 
# EVERYWHERE. Also the chemical reaction network.
# 
# Eryn Cangi
# Created December 2019
# Last edited: 21 July 2020
# Currently tested for Julia: 1.4.1
################################################################################

research_dir = "/home/emc/GDrive-CU/Research-Modeling/FractionationFactor/Code/"
results_dir = research_dir*"Results/"
main_cases_dir = "MainCases/"
det_cases_dir = "DetailedCases/"


# fundamental constants ========================================================
const boltzmannK = 1.38e-23;    # J/K
const bigG = 6.67e-11;          # N m^2/kg^2
const mH = 1.67e-27;            # kg
const marsM = 0.1075*5.972e24;  # kg
const radiusM = 3396e5;         # cm
DH = 5.5 * 1.6e-4               # SMOW value from Yung 1988

# Altitude grid discretization =================================================
const alt = convert(Array, (0:2e5:250e5))
const intaltgrid = round.(Int64, alt/1e5)[2:end-1];
const zmin = alt[1]
const zmax = alt[end];
const dz = alt[2]-alt[1];
n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])

hygropause_alt = 40e5

# Temperatures and water stuff =================================================
global meanTs = 216.0
global meanTt = 130.0
global meanTe = 205.0
global meantemps = [meanTs, meanTt, meanTe]

global meanTsint = 216
global meanTtint = 130
global meanTeint = 205

# These are the low and high values for the "standard atmosphere and reasonable 
# climate variations" cases. NOT the full range of temperatures used to make the
# detailed cases stuff, because those include temps up to 350 K.
global lowTs = 160.0
global hiTs = 270.0
global lowTt = 100.0
global hiTt = 160.0
global lowTe = 150.0
global hiTe = 250.0

global highestTe = 350.0  # This is the highest exobase temperature considered.
                          # AFAIK this parameter is only used in making the 3-panel
                          # temperature profile plot.

MR_mean_water = 1.38e-4

# Lists ========================================================================
const fullspecieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                         :O1D, :H, :N2, :Ar, :CO2pl, :HOCO,
                         # species for deuterium chemistry:
                         :HDO, :OD, :HDO2, :D, :DO2, :HD, :DOCO];

specieslist=fullspecieslist;  

# array of species for which photolysis is important. All rates should
# start with J and contain a species in specieslist above, which is used to 
# compute photolysis. 
const Jratelist=[:JCO2ion,:JCO2toCOpO,:JCO2toCOpO1D,:JO2toOpO,:JO2toOpO1D,
                 :JO3toO2pO,:JO3toO2pO1D,:JO3toOpOpO,:JH2toHpH,:JOHtoOpH,
                 :JOHtoO1DpH,:JHO2toOHpO,:JH2OtoHpOH,:JH2OtoH2pO1D,:JH2OtoHpHpO,
                 :JH2O2to2OH,:JH2O2toHO2pH,:JH2O2toH2OpO1D,
                 # deuterated species J rates:
                 :JHDOtoHpOD, :JHDOtoDpOH, :JHDO2toOHpOD,
                 # new March 2018
                 :JHDOtoHDpO1D, :JHDOtoHpDpO, :JODtoOpD, :JHDtoHpD, :JDO2toODpO,
                 :JHDO2toDO2pH, :JHDO2toHO2pD, :JHDO2toHDOpO1D, :JODtoO1DpD];

const nochemspecies = [:N2, :Ar, :CO2pl, :H2O, :HDO];
const chemspecies = setdiff(specieslist, nochemspecies);
const notransportspecies = [:CO2pl, :H2O, :HDO];
const transportspecies = setdiff(specieslist, notransportspecies);
const speciesmolmasslist = Dict(:CO2=>44, :O2=>32, :O3=>48, :H2=>2, :OH=>17,
                                :HO2=>33, :H2O=>18, :H2O2=>34, :O=>16, :CO=>28,
                                :O1D=>16, :H=>1, :N2=>28, :Ar=>40, :CO2pl=>44,
                                :HOCO=>45,
                                # deuterium chemistry:
                                :HDO=>19, :OD=>18, :HDO2=>35, :D=>2, :DO2=>34,
                                :HD=>3, :DOCO=>46);


# Plotty plot plot stuff =======================================================

rcParams = PyDict(matplotlib."rcParams")

speciescolor = Dict( # H group
                :H => "#ff0000", :D => "#ff0000", # red
                :H2 => "#e526d7", :HD =>  "#e526d7", # dark pink/magenta

                # hydroxides
                :OH => "#7700d5", :OD => "#7700d5", # purple

                # water group (roughly, I ain't a chemist)
                :H2O => "#0083dc", :HDO => "#0083dc", # cornflower blue
                :H2O2 => "#0000ff", :HDO2 => "#0000ff", # true blue
                :HO2 => "#046868", :DO2 => "#046868",  # dark teal

                # O group
                :O1D => "#808000", # olive
                :O => "#1a6115",   # forest green
                :O2 => "#15da09",  # kelly/grass green
                :O3 => "#269e56",  # light green

                # CO group
                :CO2 => "#d18564",   # dark peach
                :CO2pl => "#614215", # brown
                :CO => "#ff6600",    # orange
                :HOCO => "#e8ba8c", :DOCO => "#e8ba8c",  #tannish

                # nonreactants
                :Ar => "#808080", :N2 => "#cccccc",);  # grays

speciesstyle = Dict( # H group
                :H => "-", :D => "--",
                :H2 => "-",  :HD => "--",
                # hydroxides
                :OH => "-", :OD => "--",
                # "water group" (roughly, I ain't a chemist)
                :H2O => "-", :HDO => "--",
                :H2O2 => "-", :HDO2 => "--",
                :HO2 => "-", :DO2 => "--",
                # O group
                :O1D => "-", :O => "-", :O2 => "-", :O3 => "-",
                # CO group
                :CO => "-", :CO2 => "-", :CO2pl => "-",
                :HOCO => "-", :DOCO => "--",
                # nonreactants
                :Ar => "-", :N2 => "-",);
                
medgray = "#444444"

# Chemistry ====================================================================
# function to replace three body rates with the recommended expression
threebody(k0, kinf) = :($k0*M/(1+$k0*M/$kinf)*0.6^((1+(log10($k0*M/$kinf))^2)^-1))
threebodyca(k0, kinf) = :($k0/(1+$k0/($kinf/M))*0.6^((1+(log10($k0/($kinf*M)))^2)^-1))

################################################################################
############################### REACTION NETWORK ###############################
################################################################################

# reactions and multipliers on base rates for deuterium reactions from Yung
# 1988; base rates from this work or Chaffin+ 2017. Note: below, H-ana means 
# the same reaction but with only H-bearing species.
reactionnet = [   #Photodissociation
             [[:CO2], [:CO, :O], :JCO2toCOpO],
             [[:CO2], [:CO, :O1D], :JCO2toCOpO1D],
             [[:O2], [:O, :O], :JO2toOpO],
             [[:O2], [:O, :O1D], :JO2toOpO1D],
             [[:O3], [:O2, :O], :JO3toO2pO],
             [[:O3], [:O2, :O1D], :JO3toO2pO1D],
             [[:O3], [:O, :O, :O], :JO3toOpOpO],
             [[:H2], [:H, :H], :JH2toHpH],
             [[:HD], [:H, :D], :JHDtoHpD],
             [[:OH], [:O, :H], :JOHtoOpH],
             [[:OH], [:O1D, :H], :JOHtoO1DpH],
             [[:OD], [:O, :D], :JODtoOpD],
             [[:OD], [:O1D, :D], :JODtoO1DpD],
             [[:HO2], [:OH, :O], :JHO2toOHpO], # other branches should be here, but
                                               # have not been measured
             [[:DO2], [:OD, :O], :JDO2toODpO],
             [[:H2O], [:H, :OH], :JH2OtoHpOH],
             [[:H2O], [:H2, :O1D], :JH2OtoH2pO1D],
             [[:H2O], [:H, :H, :O], :JH2OtoHpHpO],
             [[:HDO], [:H, :OD], :JHDOtoHpOD], 
             [[:HDO], [:D, :OH], :JHDOtoDpOH], 
             [[:HDO], [:HD, :O1D], :JHDOtoHDpO1D], # inspiration from Yung89
             [[:HDO], [:H, :D, :O], :JHDOtoHpDpO], # inspiration from Yung89
             [[:H2O2], [:OH, :OH], :JH2O2to2OH],
             [[:H2O2], [:HO2, :H], :JH2O2toHO2pH],
             [[:H2O2], [:H2O, :O1D], :JH2O2toH2OpO1D],
             [[:HDO2], [:OH, :OD], :JHDO2toOHpOD], # Yung89
             [[:HDO2], [:DO2, :H], :JHDO2toDO2pH],
             [[:HDO2], [:HO2, :D], :JHDO2toHO2pD],
             [[:HDO2], [:HDO, :O1D], :JHDO2toHDOpO1D],

             # recombination of O
             [[:O, :O, :M], [:O2, :M], :(1.8*3.0e-33*(300 ./ T)^3.25)],
             [[:O, :O2, :N2], [:O3, :N2], :(5e-35*exp(724 ./ T))],
             [[:O, :O2, :CO2], [:O3, :CO2], :(2.5*6.0e-34*(300 ./ T)^2.4)],
             [[:O, :O3], [:O2, :O2], :(8.0e-12*exp(-2060 ./ T))],  # Sander 2011
             [[:O, :CO, :M], [:CO2, :M], :(2.2e-33*exp(-1780 ./ T))],

             # O1D attack
             [[:O1D, :O2], [:O, :O2], :(3.2e-11*exp(70 ./ T))], # verified NIST 4/3/18
             [[:O1D, :O3], [:O2, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :O3], [:O, :O, :O2], :(1.2e-10)], # verified NIST 4/3/18
             [[:O1D, :CO2], [:O, :CO2], :(7.5e-11*exp(115 ./ T))], # Sander2011. NIST: 7.41e-11*exp(120/T)
             ## O1D + H2
             [[:O1D, :H2], [:H, :OH], :(1.2e-10)],  # Sander2011. Yung89: 1e-10; NIST 1.1e-10
             [[:O1D, :HD], [:H, :OD], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1.3e-10 @298K
             [[:O1D, :HD], [:D, :OH], :(0.41*1.2e-10)], # Yung88: rate 0.41*H-ana (assumed). NIST 1e-10 @298K
             ## O1D + H2O
             [[:O1D, :H2O], [:OH, :OH], :(1.63e-10*exp(60 ./ T))], # Sander2011. Yung89: 2.2e-10; NIST: 1.62e-10*exp(65/T)
             [[:O1D, :HDO], [:OD, :OH], :(1.63e-10*exp(60 ./ T))], # Yung88: rate same as H-ana.

             # loss of H2
             [[:H2, :O], [:OH, :H], :(6.34e-12*exp(-4000 ./ T))], # KIDA <-- Baulch, D. L. 2005
             [[:HD, :O], [:OH, :D], :(4.40e-12*exp(-4390 ./ T))], # NIST
             [[:HD, :O], [:OD, :H], :(1.68e-12*exp(-4400 ./ T))], # NIST
             # HD and H2 exchange
             [[:H, :HD], [:H2, :D], :(6.31e-11*exp(-4038 ./ T))], # rate: Yung89. NIST rate is from 1959 for 200-1200K.
             [[:D, :H2], [:HD, :H], :(6.31e-11*exp(-3821 ./ T))], # NIST (1986, 200-300K): 8.19e-13*exp(-2700/T)

             ## OH + H2
             [[:OH, :H2], [:H2O, :H], :(2.8e-12*exp(-1800 ./ T))], # Sander2011. Yung89: 5.5e-12*exp(-2000/T). KIDA: 7.7E-12*exp(-2100/T). old rate from Mike: 9.01e-13*exp(-1526/T)
             [[:OH, :HD], [:HDO, :H], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # Yung88: rate (3/20)*H-ana. Sander2011: 5e-12*exp(-2130 ./ T)
             [[:OH, :HD], [:H2O, :D], :((3 ./ 20.)*2.8e-12*exp(-1800 ./ T))], # see prev line
             [[:OD, :H2], [:HDO, :H], :(2.8e-12*exp(-1800 ./ T))], # Yung88: rate same as H-ana (assumed)
             [[:OD, :H2], [:H2O, :D], :(0)], # Yung88 (assumed)
             ### [[:OD, :HD], [:HDO, :D], :(???)],  # possibilities for which I 
             ### [[:OD, :HD], [:D2O, :H], :(???)],  # can't find a rate...?

             # recombination of H. Use EITHER the first line OR the 2nd and 3rd.
             #[[:H, :H, :CO2], [:H2, :CO2],:(1.6e-32*(298 ./ T)^2.27)],
             [[:H, :H, :M], [:H2, :M], :(1.6e-32*(298 ./ T)^2.27)], # general version of H+H+CO2, rate: Justin Deighan.
             [[:H, :D, :M], [:HD, :M], :(1.6e-32*(298 ./ T)^2.27)], # Yung88: rate same as H-ana.

             [[:H, :OH, :CO2], [:H2O, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # Can't find in databases. Mike's rate.
             [[:H, :OD, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # not in Yung88. assumed rate
             [[:D, :OH, :CO2], [:HDO, :CO2], :(1.9*6.8e-31*(300 ./ T)^2)], # not in Yung88. assumed rate

             ## H + HO2
             [[:H, :HO2], [:OH, :OH], :(7.2e-11)], # Sander2011. Indep of T for 245<T<300
             [[:H, :HO2], [:H2, :O2], :(0.5*6.9e-12)], # 0.5 is from Krasnopolsky suggestion to Mike
             [[:H, :HO2], [:H2O, :O1D], :(1.6e-12)], # O1D is theoretically mandated
             [[:H, :DO2], [:OH, :OD], :(7.2e-11)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HD, :O2], :(0.5*6.9e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18
             [[:H, :DO2], [:HDO, :O1D], :(1.6e-12)], # Yung88: rate same as H-ana. verified Yung89 3/28/18. Yung88 has this as yielding HDO and O, not HDO and O1D
             [[:D, :HO2], [:OH, :OD], :(0.71*7.2e-11)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base: 7.05, minor disagreement)
             [[:D, :HO2], [:HD, :O2], :(0.71*0.5*6.9e-12)], # Yung88: rate 0.71*H-ana (assumed). verified Yung89 3/28/18 (base 7.29, minor disagreement)
             [[:D, :HO2], [:HDO, :O1D], :(0.71*1.6e-12)], # Yung88: rate 0.71*H-ana (assumed). Changed to O1D to match what Mike put in 3rd line from top of this section.
             [[:H, :DO2], [:HO2, :D], :(1e-10 ./ (0.54*exp(890 ./ T)))], # Yung88 (assumed) - turn off for Case 2
             [[:D, :HO2], [:DO2, :H], :(1.0e-10)], # Yung88. verified Yung89 3/28/18 - turn off for Case 2

             ## H + H2O2. deuterated analogues added 3/29
             [[:H, :H2O2], [:HO2, :H2],:(2.81e-12*exp(-1890 ./ T))], # verified NIST 4/3/18. Only valid for T>300K. No experiment for lower.
             # [[:H, :HDO2], [:DO2, :H2], :(0)], # Cazaux2010: branching ratio = 0
             # [[:H, :HDO2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:DO2, :H2], :(0)], # Cazaux2010: BR = 0
             # [[:D, :H2O2], [:HO2, :HD], :(0)], # Cazaux2010: BR = 0
             [[:H, :H2O2], [:H2O, :OH],:(1.7e-11*exp(-1800 ./ T))], # verified NIST 4/3/18
             [[:H, :HDO2], [:HDO, :OH], :(0.5*1.16e-11*exp(-2110 ./ T))], # Cazaux2010: BR = 0.5. Rate for D + H2O2, valid 294<T<464K, NIST, 4/3/18
             [[:H, :HDO2], [:H2O, :OD], :(0.5*1.16e-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:HDO, :OH], :(0.5*1.16e-11*exp(-2110 ./ T))], # see previous line
             [[:D, :H2O2], [:H2O, :OD], :(0.5*1.16e-11*exp(-2110 ./ T))], # see previous line
             [[:D, :HDO2], [:OD, :HDO], :(0.5*1.16e-11*exp(-2110 ./ T))], # added 4/3 with assumed rate from other rxns
             [[:D, :HDO2], [:OH, :D2O], :(0.5*1.16e-11*exp(-2110 ./ T))], # sourced from Cazaux et al

             # Interconversion of odd H
             ## H + O2
             [[:H, :O2], [:HO2], threebody(:(2.0*4.4e-32*(T/300.)^-1.3), # Sander2011, 300K+. Yung89: 5.5e-32(T/300)^-1.6, 7.5e-11 valid 200-300K.
                                           :(7.5e-11*(T/300.)^0.2))],  # NIST has the temp info.
             [[:D, :O2], [:DO2], threebody(:(2.0*4.4e-32*(T/300.)^-1.3), # Yung88: rate same as H-ana.
                                           :(7.5e-11*(T/300.)^0.2))],

             ## H + O3
             [[:H, :O3], [:OH, :O2], :(1.4e-10*exp(-470 ./ T))], # verified Yung89, NIST 4/3/18
             [[:D, :O3], [:OD, :O2], :(0.71*1.4e-10*exp(-470 ./ T))], # Yung88: rate 0.71*H-ana (assumed). verified Yung89, NIST 4/3/18.
             ## O + OH
             [[:O, :OH], [:O2, :H], :(1.8e-11*exp(180 ./ T))], # Sander2011. KIDA+NIST 4/3/18 150-500K: 2.4e-11*exp(110 ./ T). Yung89: 2.2e-11*exp(120/T) for both this and D analogue.
             [[:O, :OD], [:O2, :D], :(1.8e-11*exp(180 ./ T))], # Yung88: rate same as H-ana.
             ## O + HO2
             [[:O, :HO2], [:OH, :O2], :(3.0e-11*exp(200 ./ T))], # Sander2011. KIDA (220-400K): 2.7e-11*exp(224/T)
             [[:O, :DO2], [:OD, :O2], :(3.0e-11*exp(200 ./ T))], # Yung88: rate same as H-ana. verified Yung89 4/3/18
             ## O + H2O2
             [[:O, :H2O2], [:OH, :HO2], :(1.4e-12*exp(-2000 ./ T))], # Sander2011. verified NIST 4/3/18.
             [[:O, :HDO2], [:OD, :HO2], :(0.5*1.4e-12*exp(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             [[:O, :HDO2], [:OH, :DO2], :(0.5*1.4e-12*exp(-2000 ./ T))], # Yung88: rate same as H-ana (assumed). verified Yung89 4/3/18
             ## OH + OH
             [[:OH, :OH], [:H2O, :O], :(4.2e-12*exp(-240 ./ T))], # NIST+KIDA, 200-350K: 6.2e-14*(T/300)^2.62*exp(945 ./ T) changed 4/3/18. Yung89: 4.2e-12*exp(-240/T). old rate w/mystery origin: 1.8e-12.
             [[:OD, :OH], [:HDO, :O], :(4.2e-12*exp(-240 ./ T))], # Yung88: rate same as H-ana
             [[:OH, :OH], [:H2O2], threebody(:(1.3*6.9e-31*(T/300.)^-1.0),:(2.6e-11))], # Sander2011. Why 1.3?
             [[:OD, :OH], [:HDO2], threebody(:(1.3*6.9e-31*(T/300.)^-1.0),:(2.6e-11))], # Yung88: rate same as H-ana
             ## OH + O3
             [[:OH, :O3], [:HO2, :O2], :(1.7e-12*exp(-940 ./ T))], # Sander2011, temp by NIST 220-450K. Yung89: 1.6 not 1.7 -> temp 200-300K by NIST (older info)
             [[:OD, :O3], [:DO2, :O2], :(1.7e-12*exp(-940 ./ T))], # Yung88: rate same as H-ana
             ## OH + HO2
             [[:OH, :HO2], [:H2O, :O2], :(4.8e-11*exp(250 ./ T))], # verified NIST 4/3/18. Yung89: 4.6e-11*exp(230/T) for this and next 2.
             [[:OH, :DO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ T))], # Yung88: same as H-ana.
             [[:OD, :HO2], [:HDO, :O2], :(4.8e-11*exp(250 ./ T))], # Yung88: same as H-ana.
             ## OH + H2O2
             [[:OH, :H2O2], [:H2O, :HO2], :(2.9e-12*exp(-160 ./ T))], # NIST+KIDA 4/3/18, valid 240-460K. Yung89: 3.3e-12*exp(-200/T). Sander2011 recommends an average value of 1.8e-12, but this seems too high for martian temps
             [[:OD, :H2O2], [:HDO, :HO2], :(2.9e-12*exp(-160 ./ T))], # Yung88: same as H-ana (assumed)
             [[:OD, :H2O2], [:H2O, :DO2], :(0)],  # Yung88 (assumed)
             [[:OH, :HDO2], [:HDO, :HO2], :(0.5*2.9e-12*exp(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             [[:OH, :HDO2], [:H2O, :DO2], :(0.5*2.9e-12*exp(-160 ./ T))], # Yung88: rate 0.5*H-ana.
             ## HO2 + O3
             [[:HO2, :O3], [:OH, :O2, :O2], :(1.0e-14*exp(-490 ./ T))], # Sander2011. Yung89: 1.1e-14*exp(-500/T). KIDA 250-340K: 2.03e-16*(T/300)^4.57*exp(693/T). All give comparable rate values (8.6e-16 to 1e-15 at 200K)
             [[:DO2, :O3], [:OD, :O2, :O2], :(1.0e-14*exp(-490 ./ T))], # Yung88: same as H-ana (assumed)
             ## HO2 + HO2
             [[:HO2, :HO2], [:H2O2, :O2], :(3.0e-13*exp(460 ./ T))], # Sander2011. Yung89: 2.3e-13*exp(600/T). KIDA 230-420K: 2.2e-13*exp(600/T)
             [[:DO2, :HO2], [:HDO2, :O2], :(3.0e-13*exp(460 ./ T))], # Yung88: same as H-ana (assumed)
             [[:HO2, :HO2, :M], [:H2O2, :O2, :M], :(2*2.1e-33*exp(920 ./ T))], # Sander2011.
             [[:HO2, :DO2, :M], [:HDO2, :O2, :M], :(2*2.1e-33*exp(920 ./ T))], # added 3/13 with assumed same rate as H analogue

             ## OH + D or OD + H (no non-deuterated analogues)
             [[:OD, :H], [:OH, :D], :(3.3e-9*(T^-0.63)/(0.72*exp(717 ./ T)))], # rate: Yung88. NIST (Howard82): 5.25E-11*(T/298)^-0.63  - turn off for Case 2
             [[:OH, :D], [:OD, :H], :(3.3e-9*T^-0.63)], # Yung88  - turn off for Case 2

             # CO2 recombination due to odd H (with HOCO intermediate)
             ## straight to CO2
             [[:CO, :OH], [:CO2, :H], threebodyca(:(1.5e-13*(T/300.)^0.6),:(2.1e9*(T/300.)^6.1))], # Sander2011
             [[:CO, :OD], [:CO2, :D], threebodyca(:(1.5e-13*(T/300.)^0.6),:(2.1e9*(T/300.)^6.1))], # Yung88: same as H-ana.
             ### possible deuterated analogues below
             [[:OH, :CO], [:HOCO], threebody(:(5.9e-33*(T/300.)^-1.4),:(1.1e-12*(T/300.)^1.3))], # Sander2011
             [[:OD, :CO], [:DOCO], threebody(:(5.9e-33*(T/300.)^-1.4),:(1.1e-12*(T/300.)^1.3))],

             [[:HOCO, :O2], [:HO2, :CO2], :(2.09e-12)], # verified NIST 4/3/18
             [[:DOCO, :O2], [:DO2,:CO2], :(2.09e-12)],  # assumed?

             # CO2+ attack on molecular hydrogen
             [[:CO2pl, :H2], [:CO2, :H, :H], :(8.7e-10)], # from Kras 2010 / Scott 1997
             [[:CO2pl, :HD], [:CO2pl, :H, :D], :((2/5)*8.7e-10)]
             ]
