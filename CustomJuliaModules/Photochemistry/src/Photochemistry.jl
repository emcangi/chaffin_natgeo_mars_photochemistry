__precompile__()

module Photochemistry

using PyPlot
using PyCall
using HDF5, JLD
using LaTeXStrings
using Distributed
using DelimitedFiles
using SparseArrays
using LinearAlgebra
using PlotUtils

include("../../../PARAMETERS.jl")

export deletefirst, fluxsymbol, readandskip, searchdir, search_subfolders, create_folder, input, 
       plot_bg, plotatm, plot_temp_prof,
       get_ncurrent, write_ncurrent, n_tot, getpos, 
       effusion_velocity, speciesbcs, scaleH,  
       fluxcoefs, getflux, fluxes, 
       Keddy, Dcoef,
       meanmass, lossequations, loss_rate, production_equations, production_rate, chemical_jacobian, getrate, reactionrates, ratefn, chemJmat, 
       co2xsect, h2o2xsect_l, h2o2xsect, hdo2xsect, binupO2, o2xsect, ho2xsect_l, O3O1Dquantumyield, padtosolar, quantumyield,
       Tpiecewise, Psat, Psat_HDO

# basic utility functions ==========================================================

# searches path for key, which can be a regular expression
searchdir(path, key) = filter(x->occursin(key,x), readdir(path))

function deletefirst(A, v)
    #=
    A: list
    v: any; an element that may be present in list A

    returns: list A with its first element equal to v removed.
    =#
    index = something(findfirst(isequal(v), A), 0)  # this horrible syntax introduced by Julia.
    keep = setdiff([1:length(A);],index)
    return A[keep]
end

function fluxsymbol(x)
    #= 
    Converts string x to a symbol. f for flux. 
    =#
    return Symbol(string("f",string(x)))
end

function readandskip(f, delim::Char, T::Type; skipstart=0)
    #= 
    reads in data from a file

    f: file to read
    delim: character that functions as delimiter
    T: type to cast data to. e.g., "1" could be cast to String, Char, Int, Float, etc.
    skipstart: initial number of rows in the file to ignore 

    returns: array containing file data
    =# 
    f = open(f,"r")
    if skipstart>0
        for i in [1:skipstart;]
            readline(f)
        end
    end
    f = readdlm(f, delim, T)
end

function search_subfolders(path, key)
    #=
    path: a folder containing subfolders and files.
    key: the text pattern in folder names that you wish to find.

    Searches the top level subfolders within path for all folders matching a 
    certain regex given by key. Does not search files or sub-subfolders.
    =#
    folders = []
    for (root, dirs, files) in walkdir(path)
        if root==path
            for dir in dirs
                push!(folders, joinpath(root, dir)) # path to directories
            end
        end
    end

    wfolders = filter(x->occursin(key, x), folders)
    return wfolders
end

function create_folder(foldername, parentdir)
    #=
    Checks to see if foldername exists within parentdir. If not, creates it.
    =#
    println("Checking for existence of $(foldername) folder in $(parentdir)")
    dircontents = readdir(parentdir)
    if foldername in dircontents
        println("Folder $(foldername) exists")
    else
        mkdir(parentdir*foldername)
        println("Created folder ", foldername)
    end
end

function input(prompt::String="")::String
    #=
    Prints a prompt to the terminal and accepts user input, returning it as a string.
    =#
    print(prompt)
    return chomp(readline())
end

# plot functions ===============================================================

function plot_bg(axob)
    #=
    Takes a PyPlot axis object axob and
    1) sets the background color to a light gray
    2) creates a grid for the major tick marks
    3) turns off all the axis borders

    intended to immitate seaborn
    =#
    axob.set_facecolor("#ededed")
    axob.grid(zorder=0, color="white", which="major")
    for side in ["top", "bottom", "left", "right"]
        axob.spines[side].set_visible(false)
    end
end

function plotatm(n_current; t=nothing, iter=nothing, Hsp_only=false, savepath=nothing)
    #=
    Makes a "spaghetti plot" of the species concentrations by altitude in the
    atmosphere. 

    n_current: array of species densities by altitude
    t: timestep, for plotting the atmosphere during convergence
    iter: iteration, for plotting the atmosphere during convergence
    Hsp_only: whether to plot all species or H-bearing species only (other than water)
    savepath: where to save the resulting .png file
    =#

    clf()
    
    if Hsp_only==false
        species_to_plot = fullspecieslist
    else
        species_to_plot = [:OH, :OD, :HD, :H2, :H, :D, :H2O2, :HDO2, :DO2, :HO2]
    end

    for sp in species_to_plot
        plot(n_current[sp], alt[2:end-1]/1e5, color = speciescolor[sp],
             linewidth=2, label=sp, linestyle=speciesstyle[sp], zorder=1)
    end
    tight_layout()
    ylim(0, zmax/1e5)
    xscale("log")
    xlim(1e-15, 1e18)
    if Hsp_only==true
        xlim(1e-15, 1e18)
    end
    xlabel(L"Species concentration (cm$^{-3}$)")
    ylabel("Altitude [km]")
    if t!=nothing && iter!=nothing
        title("time = $(t), iteration $(iter)")
    end
    legend(bbox_to_anchor=[1.01,1], loc=2, borderaxespad=0)# for a limited-species plot: loc="lower left", fontsize=12)#
    # convfig.canvas.draw()  # use if you want the plots in the live figure to update as the simulation converges

    if savepath != nothing
        savefig(savepath*"Hspec_only.png", bbox_inches="tight")
    end
end

function plot_temp_prof(temps_to_plot, savepath, alt)
    #=
    Creates a .png image of the tepmeratures plotted by altitude in the atmosphere

    temps_to_plot: an array of temperatures by altitude
    savepath: where to save the resulting .png image
    alt: the altitude array; passed in so that you could select other altitudes if you wanted
    =#

    fig, ax = subplots(figsize=(4,6))
    ax.set_facecolor("#ededed")
    grid(zorder=0, color="white", which="both")
    for side in ["top", "bottom", "left", "right"]
        ax.spines[side].set_visible(false)
    end

    plot(temps_to_plot, alt/1e5)

    ax.set_ylabel("Altitude (km)")
    ax.set_yticks([0, 100, alt[end]/1e5])
    ax.set_yticklabels([0,100, alt[end]/1e5])
    ax.set_xlabel("Temperature (K)")

    ax.text(temps_to_plot[end]*0.9, 185, L"T_{exo}")
    ax.text(temps_to_plot[Int64(length(temps_to_plot)/2)]+5, 75, L"T_{tropo}")
    ax.text(temps_to_plot[1], 10, L"T_{surface}")

    savefig(savepath*"/temp_profile.png", bbox_inches="tight") 
end

# Array manipulation ===========================================================
function get_ncurrent(readfile)
    #=
    Retrieves the matrix of species concentrations by altitude from an HDF5
    file, readfile, containing a converged atmosphere.
    =#
    n_current_tag_list = map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,"n_current/n_current_mat");
    n_current = Dict{Symbol, Array{Float64, 1}}()

    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]] = reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    return n_current
end

function write_ncurrent(n_current, filename)
    #=
    Writes out the current state of species density by altitude to an .h5 file
    (for converged atmosphere). 

    n_current: the array of the current atmospheric state in species density by altitude
    filename: whatever you want to call the filename when it's saved as .h5
    =# 
    n_current_mat = Array{Float64}(undef, length(alt)-2, 
                                   length(collect(keys(n_current))));
    for ispecies in [1:length(collect(keys(n_current)));]
        for ialt in [1:length(alt)-2;]
            n_current_mat[ialt, ispecies] = n_current[collect(keys(n_current))[ispecies]][ialt]
        end
    end
    h5write(filename,"n_current/n_current_mat",n_current_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string, collect(keys(n_current))))
end

function write_ncurrent(n_current, filename, alt)
    #=
    Override to allow alt (altitude array) to be passed in because of tweaks I made to other functions.
    =# 
    n_current_mat = Array{Float64}(undef, length(alt)-2, 
                                   length(collect(keys(n_current))));
    for ispecies in [1:length(collect(keys(n_current)));]
        for ialt in [1:length(alt)-2;]
            n_current_mat[ialt, ispecies] = n_current[collect(keys(n_current))[ispecies]][ialt]
        end
    end
    h5write(filename,"n_current/n_current_mat",n_current_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string, collect(keys(n_current))))
end

function n_tot(n_current, z, n_alt_index)
    #= 
    Calculates the total number density in #/cm^3 at a given altitude.

    n_current: the array of the current atmospheric state in species density by altitude
    z: altitude in question
    n_alt_index: index array, shuldn't really need to be passed in because it's global
    =#
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

function n_tot(n_current, z)
    #= 
    Override for when n_alt_index is global
    =#
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

function getpos(array, test::Function, n=Any[])
    #= 
    Get the position of keys that match "test" in an array

    array: Array; any size and dimension
    test: used to test the elements in the array
    n: Array; for storing the indices (I think)

    returns: 1D array of elements in array that match test function supplied
    =#
    if !isa(array, Array)
        test(array) ? Any[n] : Any[]
    else
        vcat([ getpos(array[i], test, Any[n...,i]) for i=1:size(array)[1] ]...)
    end
end

function getpos(array, value)
    #= 
    overloading getpos for the most common use case, finding indicies
    corresponding to elements in array that match value. 
    =#
    getpos(array, x->x==value)
end

# boundary condition functions =================================================
function effusion_velocity(Texo::Float64, m::Float64, zmax)
    #=
    Returns effusion velocity for a species in cm/s

    Texo: temperature of the exobase (upper boundary) in K
    m: mass of one molecule of species in amu
    zmax: max altitude in cm
    =#
    
    # lambda is the Jeans parameter (Gronoff 2020), basically the ratio of the 
    # escape velocity GmM/z to the thermal energy, kT.
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+zmax))
    vth = sqrt(2*boltzmannK*Texo/(m*mH))  # this one is in m/s
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)  # this is in cm/s
    return v
end

function speciesbcs(species, speciesbclist)
    #=
    Returns boundary conditions for a given species. If species not in the dict,
    flux = 0 at top and bottom will be the returned boundary condition.

    species: symbol for species passed in
    speciesbclist: It has to be passed in because it takes the H2O and HDO 
                   saturation at the surface as bcs, and these are variable.
                   Same story for effusion velocities (which depend on the
                   variable T_exo) and Of (flux of atomic O), which is not 
                   usually varied anymore but once was for testing purposes.
    
    =#
    get(speciesbclist, species, ["f" 0.; "f" 0.])
end

# scale height functions =======================================================

function scaleH(z, T::Float64, mm::Real)
    #= 
    Computes the general scale height (cm) of the atmosphere for the mean molar mass at some altitude

    z: Float or Int; unit: cm. altitude in atmosphere at which to calculate scale height
    T: temperature in Kelvin
    mm: mean molar mass, in amu. 
    =#
    return boltzmannK*T/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
    # constants are in MKS. Convert to m and back to cm.
end

function scaleH(z, T::Float64, species::Symbol)
    #=
    Overload: Scale height for a particular species
    =#
    mm = speciesmolmasslist[species]
    scaleH(z, T, mm)
end

function scaleH(z, T::Float64, n_current)
    #= 
    Overload: scale height when mean mass is not provided but general atmospheric
              state (n_current) is
    =#
    mm = meanmass(n_current, z)
    scaleH(z, T, mm)
end


# flux functions ==========================================================

function fluxcoefs(z, dz, Kv, Dv, Tv, Hsv, H0v, species)
    #= 
    base function to generate coefficients of the transport network. 

    z: Float64; altitude in cm.
    dz: Float64; altitude layer thickness ("resolution")
    Kv: Array; 3 elements (lower, this, and upper layer). eddy diffusion coefficient
    Dv: Array; 3 elements, same as K. molecular diffusion coefficient
    Tv: Array; 3 elements, same as K. temperature
    Hsv: Array; 3 elements, same as K. scale height by species
    H0v: Array; 3 elements, same as K. mean atmospheric scale height
    species: species symbol 

    v refers to "vector"
    u refers to "upper" (the upper parcel)
    l refers to "lower" (the lower parcel)

    Returns: a list of the coefficients for [downward, upward] flux.
    units are in 1/s.
    =#

    # Calculate the coefficients between this layer and the lower layer. 
    Dl = (Dv[1] + Dv[2])/2.0
    Kl = (Kv[1] + Kv[2])/2.0
    Tl = (Tv[1] + Tv[2])/2.0
    dTdzl = (Tv[2] - Tv[1])/dz
    Hsl = (Hsv[1] + Hsv[2])/2.0
    H0l = (H0v[1] + H0v[2])/2.0

    # two flux terms: eddy diffusion and gravity/thermal diffusion.
    # these are found in line 5 of Mike's transport_as_chemistry.pdf. 
    sumeddyl = (Dl+Kl)/dz/dz
    gravthermall = (Dl*(1/Hsl + (1+thermaldiff(species))/Tl*dTdzl) +
                    Kl*(1/H0l + 1/Tl*dTdzl))/(2*dz)

    # Now the coefficients between this layer and upper layer.
    Du = (Dv[2] + Dv[3])/2.0
    Ku = (Kv[2] + Kv[3])/2.0
    Tu = (Tv[2] + Tv[3])/2.0
    dTdzu = (Tv[3] - Tv[2])/dz
    Hsu = (Hsv[2] + Hsv[3])/2.0
    H0u = (H0v[2] + H0v[3])/2.0

    sumeddyu = (Du+Ku)/dz/dz  # this is the line where we divide by cm^2
    gravthermalu = (Du*(1/Hsu + (1 + thermaldiff(species))/Tu*dTdzu) +
                  Ku*(1/H0u + 1/Tu*dTdzu))/(2*dz)

    # this results in the following coupling coefficients:
    return [sumeddyl+gravthermall, # down
            sumeddyu-gravthermalu] # up
end

function getflux(n_current, dz, species)
    #=
    Returns a 1D array of fluxes in and out of a given altitude level for a 
    given species. For looking at vertical distribution of fluxes. does 
    not modify the concentrations.

    n_current: Array; species number density by altitude
    dz: Float64; layer thickness in cm
    species: Symbol

    returns: Array of flux at each altitude layer (#/cm^2/s)
    =#

    # each element in thesecoefs has the format [downward, upward]
    thesecoefs = [fluxcoefs(a, dz, species, n_current) for a in alt[2:end-1]]

    # thesebcs has the format [lower bc; upper bc], where each row contains a 
    # character showing the type of boundary condition, and a number giving its value
    thesebcs = boundaryconditions(species, dz, n_current)

    thesefluxes = fill(convert(Float64, NaN),length(intaltgrid))

    # in the following line for the lowest layer: 
    # first term is -(influx from layer above - outflux from this layer)
    # second term is (-this layer's lower bc that depends on concentration + bc that doesn't depend on concentration)
    thesefluxes[1] = (-(n_current[species][2]*thesecoefs[2][1]
                        -n_current[species][1]*thesecoefs[1][2]) 
                    +(-n_current[species][1]*thesebcs[1, 1]
                      +thesebcs[1, 2]))/2.0
    for ialt in 2:length(intaltgrid)-1
        thesefluxes[ialt] = (-(n_current[species][ialt+1]*thesecoefs[ialt+1][1]
                               -n_current[species][ialt]*thesecoefs[ialt][2])
                             +(-n_current[species][ialt]*thesecoefs[ialt][1]
                               +n_current[species][ialt-1]*thesecoefs[ialt-1][2]))/2.0
    end
    thesefluxes[end] = (-(thesebcs[2, 2]
                          - n_current[species][end]*thesebcs[2, 1])
                        + (-n_current[species][end]*thesecoefs[end][1]
                           +n_current[species][end-1]*thesecoefs[end-1][2]))/2.0
    return dz*thesefluxes
end

function fluxes(n_current, dz)
    #=
    Just runs getflux for all species
    =#
    thesefluxes = fill(convert(Float64, NaN),(length(intaltgrid),length(specieslist)))
    for isp in 1:length(specieslist)
        thesefluxes[:,isp] = getflux(n_current, dz, specieslist[isp])
    end
    return thesefluxes
end

# diffusion functions ==========================================================

function Dcoef(T, n::Real, species::Symbol)
    #=
    Calculates molecular diffusion coefficient for a particular slice of the
    atmosphere using D = AT^s/n, from Banks and Kockarts Aeronomy, part B,
    pg 41, eqn 15.30 and table 15.2 footnote

    T: temperature (K)
    n: number of molecules (all species) this altitude, usually calculated by n_tot
    species: whichever species we are calculating for
    =#
    dparms = diffparams(species)
    return dparms[1]*1e17*T^(dparms[2])/n
end

function Keddy(z::Real, nt::Real)
    #=
    eddy diffusion coefficient, stolen from Krasnopolsky (1993).
    Scales as the inverse sqrt of atmospheric number density

    z: some altitude in cm.
    nt: number total of species at this altitude (I think)
    =#
    z <= 60.e5 ? 10^6 : 2e13/sqrt(nt)
end

function Keddy(n_current, z)
    #=
    Override for when nt is not provided, and only the current atmospheric state
    (n_current) is
    =#
    z <= 60.e5 ? 10^6 : 2e13/sqrt(n_tot(n_current, z))
end

#=
molecular diffusion parameters. value[1] = A, value[2] = s in the equation
D = AT^s / n given by Banks & Kockarts Aeronomy part B eqn. 15.30 and Hunten
1973, Table 1.

molecular diffusion is different only for small molecules and atoms
(H, D, HD, and H2), otherwise all species share the same values (Krasnopolsky
1993 <- Hunten 1973; Kras cites Banks & Kockarts, but this is actually an
incorrect citation.)

D and HD params are estimated by using Banks & Kockarts eqn 15.29 (the
coefficient on T^0.5/n) to calculate A for H and H2 and D and HD, then using
(A_D/A_H)_{b+k} = (A_D/A_H)_{hunten} to determine (A_D)_{hunten} since Hunten
1973 values are HALF what the calculation provides.

s, the power of T, is not calculated because there is no instruction on how to
do so and is assumed the same as for the hydrogenated species.
=#

diffparams(species) = get(Dict(:H=>[8.4, 0.597], :H2=>[2.23, 0.75],
                               :D=>[5.98, 0.597], :HD=>[1.84, 0.75]), # 5.98 and 1.84
                               species,[1.0, 0.75])


# thermal diffusion factors (all verified by Krasnopolsky 2002)
thermaldiff(species) = get(Dict(:H=>-0.25, :H2=>-0.25, :D=>-0.25, :HD=>-0.25,
                                :He=>-0.25), species, 0)


# chemistry functions ==========================================================

function meanmass(n_current, z)
    #= 
    find the mean molecular mass at a given altitude 

    n_current: Array; species number density by altitude
    z: Float64; altitude in atmosphere in cm

    return: mean molecular mass in amu
    =#
    thisaltindex = n_alt_index[z]
    c = [n_current[sp][thisaltindex] for sp in specieslist]
    m = [speciesmolmasslist[sp] for sp in specieslist]
    return sum(c.*m)/sum(c)
end

function loss_equations(network, species)
    #=  
    given a network of equations in the form of reactionnet, this
    function returns the loss equations and rate coefficient for all
    reactions where the supplied species is consumed, in the form of an array
    where each entry is of the form [reactants, rate] 
    =#

    # get list of all chemical reactions species participates in:
    speciespos = getpos(network, species)
    # find pos where species is on LHS but not RHS:
    lhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(lhspos, rhspos)
        lhspos = deletefirst(lhspos, i)
    end

    # get the products and rate coefficient for the identified reactions.
    losseqns=map(x->vcat(Any[network[x][1]...,network[x][3]]), lhspos)
    # automatically finds a species where it occurs twice on the LHS
end

function loss_rate(network, species)
    #= 
    return a symbolic expression for the loss rate of species in the
    supplied reaction network. Format is a symbolic expression containing a sum
    of reactants * rate. 
    =#
    leqn=loss_equations(network, species) # get the equations
    lval=:(+($( # and add the products together
               map(x->:(*($(x...))) # take the product of the
                                    # concentrations and coefficients
                                    # for each reaction
                   ,leqn)...)))
end

function production_equations(network, species)
    #= 
    given a network of equations in the form of reactionnet, this
    function returns the production equations and rate coefficient for all
    reactions where the supplied species is produced, in the form of an array
    where each entry is of the form [reactants, rate] =#

    speciespos = getpos(network, species)#list of all reactions where species is produced
    # find pos where species is on RHS but not LHS
    lhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], #select the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(rhspos, lhspos)
        rhspos = deletefirst(rhspos, i)
    end

    # get the products and rate coefficient for the identified reactions.
    prodeqns = map(x->vcat(Any[network[x][1]...,network[x][3]]),
                 # automatically finds and counts duplicate
                 # production for each molecule produced
                 rhspos)

    return prodeqns
end

function production_rate(network, species)
    #= 
    return a symbolic expression for the loss rate of species in the
    supplied reaction network.
    =#

    # get the reactants and rate coefficients
    peqn = production_equations(network, species)

    # add up and take the product of each set of reactants and coeffecient
    pval = :(+ ( $(map(x -> :(*($(x...))), peqn) ...) ))
end

function chemical_jacobian(chemnetwork, transportnetwork, specieslist, dspecieslist)
    #= 
    Compute the symbolic chemical jacobian of a supplied chemnetwork and transportnetwork
    for the specified specieslist. Returns three arrays suitable for
    constructing a sparse matrix: lists of the first and second indices
    and the symbolic value to place at that index.
    =#

    # set up output vectors: indices and values
    ivec = Int64[] # list of first indices (corresponding to the species being produced and lost)
    jvec = Int64[] # list of second indices (corresponding to the derivative being taken)
    tvec = Any[] # list of the symbolic values corresponding to the jacobian

    nspecies = length(specieslist)
    ndspecies = length(dspecieslist)
    for i in 1:nspecies #for each species
        ispecies = specieslist[i]
        # get the production and loss equations
        peqn = []
        leqn = []
        if issubset([ispecies],chemspecies)
            peqn = [peqn; production_equations(chemnetwork, ispecies)] #***
            leqn = [leqn; loss_equations(chemnetwork, ispecies)]
        end
        if issubset([ispecies],transportspecies)
            peqn = [peqn; production_equations(transportnetwork, ispecies)]
            leqn = [leqn; loss_equations(transportnetwork, ispecies)]
        end
        for j in 1:ndspecies #now take the derivative with resp`ect to the other species
            jspecies = dspecieslist[j]
            #= find the places where the production rates depend on
            jspecies, and return the list rates with the first
            occurrance of jspecies deleted. (Note: this seamlessly
            deals with multiple copies of a species on either side of
            an equation, because it is found twice wherever it lives) =#
            ppos = map(x->deletefirst(peqn[x[1]],jspecies),getpos(peqn, jspecies))
            lpos = map(x->deletefirst(leqn[x[1]],jspecies),getpos(leqn, jspecies))
            if length(ppos)+length(lpos)>0 #if there is a dependence
                #make note of where this dependency exists
                append!(ivec,[i])
                append!(jvec,[j])
                #= smash the production and loss rates together,
                multiplying for each distinct equation, adding
                together the production and loss seperately, and
                subtracting loss from production. =#
                if length(ppos)==0
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($lval))
                elseif length(lpos)==0
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    tval = :(+($pval))
                else
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($pval,$lval))
                end
                # attach the symbolic expression to the return values
                append!(tvec,[tval])
            end
        end
    end
    return (ivec, jvec, Expr(:vcat, tvec...))
end

function getrate(chemnet, transportnet, species)
    #=
    Creates a symbolic expression for the rate at which a given species is
    either produced or lost. Production is from chemical reaction yields or
    entry from other atmospheric layers. Loss is due to consumption in reactions
    or migration to other layers.
    =#
    rate = :(0.0)
    if issubset([species],chemspecies)
        rate = :($rate
               + $(production_rate(chemnet, species))
               - $(      loss_rate(chemnet, species)))
    end
    if issubset([species],transportspecies)
        rate = :($rate
               + $(production_rate(transportnet, species))
               - $(      loss_rate(transportnet, species)))
    end

    return rate
end

function reactionrates(n_current)
    #=
    Creates an array of size length(intaltgrid) x (number of reactions).
    Populated with chemical reaction rates for each reaction based on species
    populations.
    =#
    
    theserates = fill(convert(Float64, NaN),(length(intaltgrid),length(reactionnet)))
    for ialt in 1:length(intaltgrid)
        theserates[ialt,:] = reactionrates_local([[n_current[sp][ialt] for sp in specieslist];
                                                [n_current[J][ialt] for J in Jratelist];
                                                Temp(alt[ialt+1]);
                                                n_tot(n_current, alt[ialt+1])]...)
    end
    return theserates
end

# photochemistry functions ========================================================
function co2xsect(co2xdata, T::Float64)
    #=
    Makes an array of CO2 cross sections at temperature T, of format 
    [wavelength in nm, xsect]. 

    Data are available at 195K and 295K, so crosssections at temperatures between 
    these are created by first calculating the fraction that T is along the 
    line from 195 to 295, and then calculating a sort of weighted average of 
    the crosssections at 195 and 295K based on that information. 
    =#
    clamp(T, 195, 295)
    Tfrac = (T-195)/(295-195)

    arr = [co2xdata[:,1]; (1-Tfrac)*co2xdata[:,2] + Tfrac*co2xdata[:,3]]
    reshape(arr, length(co2xdata[:,1]),2)
end

function h2o2xsect_l(l::Float64, T::Float64)
    #=
    from 260-350 the following analytic calculation fitting the
    temperature dependence is recommended by Sander 2011.

    Analytic calculation of H2O2 cross section using temperature dependencies
    l: wavelength in nm
    T: temperature in K
    =#
    l = clamp(l, 260, 350)
    T = clamp(T, 200, 400)

    A = [64761., -921.70972, 4.535649,
         -0.0044589016, -0.00004035101,
         1.6878206e-7, -2.652014e-10, 1.5534675e-13]
    B = [6812.3, -51.351, 0.11522, -0.000030493, -1.0924e-7]

    lpowA = map(n->l^n,[0:7;])
    lpowB = map(n->l^n,[0:4;])

    expfac = 1.0/(1+exp(-1265/T))

    return 1e-21*(expfac*reduce(+, map(*, A, lpowA))+(1-expfac)*reduce(+, map(*, B, lpowB)))
end

function h2o2xsect(h2o2xdata, T::Float64)
    #=
    stitches together H2O2 cross sections, some from Sander 2011 table and some
    from the analytical calculation recommended for 260-350nm recommended by the
    same.
    T: temperature in K
    =#
    retl = h2o2xdata[:,1]
    retx = 1e4*h2o2xdata[:,2] # factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    return reshape([retl; retx], length(retl), 2)
end

function hdo2xsect(hdo2xdata, T::Float64)
    #=
    Currently this is a direct copy of h2o2xsect because no HDO2 crosssections
    exist. In the future this could be expanded if anyone ever makes those
    measurements.
    =#
    retl = hdo2xdata[:,1]
    retx = 1e4*hdo2xdata[:,2] # factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    reshape([retl; retx], length(retl), 2)
end

function binupO2(list)
    #=
    For O2

    ...details on this missing
    =#
    ret = Float64[];
    for i in [176:203;]
        posss = getpos(list[:,1],x->i<x<i+1)
        dl = diff([map(x->list[x[1],1],posss); i])
        x0 = map(x->list[x[1],2], posss)
        x1 = map(x->list[x[1],3], posss)
        x2 = map(x->list[x[1],4], posss)
        ax0 = reduce(+,map(*,x0, dl))/reduce(+, dl)
        ax1 = reduce(+,map(*,x1, dl))/reduce(+, dl)
        ax2 = reduce(+,map(*,x2, dl))/reduce(+, dl)
        append!(ret,[i+0.5, ax0, ax1, ax2])
    end
    return transpose(reshape(ret, 4, 203-176+1))
end

function o2xsect(o2xdata, o2schr130K, o2schr190K, o2schr280K, T::Float64)
    #=
    For O2, used in quantum yield calculations,
    including temperature-dependent Schumann-Runge bands.

    ...details missing
    =# 
    
    o2x = deepcopy(o2xdata);
    # fill in the schumann-runge bands according to Minschwaner 1992
    T = clamp(T, 130, 500)
    if 130<=T<190
        o2schr = o2schr130K
    elseif 190<=T<280
        o2schr = o2schr190K
    else
        o2schr = o2schr280K
    end

    del = ((T-100)/10)^2

    for i in [176.5:203.5;]
        posO2 = something(findfirst(isequal(i), o2x[:, 1]), 0)
        posschr = something(findfirst(isequal(i), o2schr[:, 1]), 0)
        o2x[posO2, 2] += 1e-20*(o2schr[posschr, 2]*del^2
                                + o2schr[posschr, 3]*del
                                + o2schr[posschr, 4])
    end

    # add in the herzberg continuum (though tiny)
    # measured by yoshino 1992
    for l in [192.5:244.5;]
        posO2 = something(findfirst(isequal(l), o2x[:, 1]), 0)
        o2x[posO2, 2] += 1e-24*(-2.3837947e4
                            +4.1973085e2*l
                            -2.7640139e0*l^2
                            +8.0723193e-3*l^3
                            -8.8255447e-6*l^4)
    end

    return o2x
end

function ho2xsect_l(l::Float64)
    #= 
    compute HO2 cross-section as a function of wavelength l in nm, as given by 
    Sander 2011 JPL Compilation 
    =#
    a = 4.91
    b = 30612.0
    sigmamed = 1.64e-18
    vmed = 50260.0
    v = 1e7/l;
    if 190<=l<=250
        return HO2absx = sigmamed / ( 1 - b/v ) * exp( -a * log( (v-b)/(vmed-b) )^2 )
    else
        return 0.0
    end
end

function O3O1Dquantumyield(lambda, temp)
    #=
    Quantum yield for O(1D) from O3 photodissociation. 
    "The quantum yield of O1D from ozone photolysis is actually well-studied! 
    This adds some complications for processing." - Mike
    =#
    if lambda < 306. || lambda > 328.
        return 0.
    end
    temp=clamp(temp, 200, 320)#expression is only valid in this T range

    X = [304.225, 314.957, 310.737];
    w = [5.576, 6.601, 2.187];
    A = [0.8036, 8.9061, 0.1192];
    v = [0.,825.518];
    c = 0.0765;
    R = 0.695;
    q = exp.(-v/(R*temp))
    qrat = q[1]/(q[1]+q[2])

    (q[1]/sum(q)*A[1]*exp.(-((X[1]-lambda)/w[1])^4.)
     +q[2]/sum(q)*A[2]*(temp/300.)^2 .* exp.(-((X[2]-lambda)/w[2])^2.)
     +A[3]*(temp/300.)^1.5*exp.(-((X[3]-lambda)/w[3])^2.)
     +c)
end

function padtosolar(solarflux, crosssection::Array{Float64, 2})
    #=
    a function to take an Nx2 array and pad it with zeroes until it's the
    same length as the solar flux. Returns the cross sections only, as
    the wavelengths are shared by solarflux 
    =#
    positions = map(x->something(findfirst(isequal(x), solarflux[:,1]), 0), crosssection[:,1])
    retxsec = fill(0.,length(solarflux[:,1]))
    retxsec[positions] = crosssection[:,2]
    return retxsec
end

function quantumyield(xsect::Array, arr)
    #= 
    function to assemble cross-sections for a given pathway. Inputs are
    an Nx2 array xsect with wavelengths and photoabsorption cross
    sections, and arr, a tuple of tuples with a condition and a quantum
    yield multiplicative factor, either constant or a function of
    wavelength in the given regime. Return is an array with all of the
    matching wavelengths and the scaled cross-sections.
    =#
    lambdas = Float64[];
    rxs = Float64[];
    for (cond, qeff) in arr
        places = findall(cond, xsect[:,1])
        append!(lambdas, xsect[places, 1])
        #if we have a number then map to a function
        isa(qeff, Function) ? (qefffn = qeff) : (qefffn = x->qeff)
        append!(rxs, map(*,map(qefffn, xsect[places, 1]),xsect[places, 2]))
    end

    return reshape([lambdas; rxs],length(lambdas),2)
end

# TEMPERATURE ==================================================================

function Tpiecewise(z::Float64, Tsurf, Ttropo, Texo)
    #= 
    DO NOT MODIFY! If you want to change the temperature, define a
    new function or select different arguments and pass to Temp(z)

    a piecewise function for temperature as a function of altitude,
    using Krasnopolsky's 2010 "half-Gaussian" function for temperatures 
    altitudes above the tropopause, with a constant lapse rate (1.4K/km) 
    in the lower atmosphere. The tropopause width is allowed to vary
    in certain cases.

    z: altitude above surface in cm
    Tsurf: Surface temperature in K
    Tropo: tropopause tempearture
    Texo: exobase temperature
    =#
    
    lapserate = -1.4e-5 # lapse rate in K/cm
    ztropo = 120e5  # height of the tropopause top
    
    # set the width of tropopause. This code allows it to vary for surface or 
    # tropopause experiments.
    ztropo_bot = (Ttropo-Tsurf)/(lapserate)
    ztropowidth = ztropo - ztropo_bot

    if z >= ztropo  # upper atmosphere
        return Texo - (Texo - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Texo))
    elseif ztropo > z >= ztropo - ztropowidth  # tropopause
        return Ttropo
    elseif ztropo-ztropowidth > z  # lower atmosphere
        return Tsurf + lapserate*z
    end
end

# WATER ========================================================================
# 1st term is a conversion factor to convert to (#/cm^3) from Pa. Source: Marti & Mauersberger 1993
Psat(T::Float64) = (1e-6/(boltzmannK * T))*(10^(-2663.5/T + 12.537))

# It doesn't matter to get the exact SVP of HDO because we never saturate. 
# However, this function is defined on the offchance someone studies HDO.
Psat_HDO(T::Float64) = (1e-6/(boltzmannK * T))*(10^(-2663.5/T + 12.537))


end