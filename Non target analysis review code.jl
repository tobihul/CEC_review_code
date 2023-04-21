using Statistics, PubChemCrawler, CSV, DataFrames, StatsPlots, Distributions,LinearAlgebra
import HTTP
const prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
#This function converts the chemical names to PubChem CID's and is implemeted in the next two functions
function my_get_cid(; name=nothing, smiles=nothing, kwargs...)
    input = "compound/"
    name !== nothing && (input *= "name/$(HTTP.escapeuri(name))/")
    smiles !== nothing && (input *= "smiles/$((smiles))/")
    url = prolog * input * "cids/TXT"
    r = HTTP.request("GET", url; kwargs...)
    cids_string = String(r.body)
    cids = split(cids_string, "\n")
    cids = [cid for cid in cids if !isempty(cid) && !isspace(cid[1])]
    return parse(Int, cids[1])
end
#Function to convert IUPAC or other names accepted in PubChem into compound descriptors
#The inputs: The number of the paper, the compound Name/SMILES, the mass analyzer used, the column volumes, the polarity, acquisition mode and MS database size
#This was used for clustering purposes in scatter plots
function logP_MW_names(Paper, Compounds, MA, CV, Polarity, Acq, MS)

    cids::Vector{Int32} = zeros(length(Compounds))
    for i = 1:length(Compounds)
        cids[i] = my_get_cid(name=Compounds[i])
        @show i
    end
    
    cidss::Vector{Int32} = trunc.(Int32, cids)

    cidssf::Vector{Int32} = filter!(e->e≠0,cidss)
    
    Df_LogP_MW = CSV.File(get_for_cids(cidssf; properties="MolecularWeight,XLogP", output="CSV")) |> DataFrame
    Df_LogP_MW_no_missing = dropmissing(Df_LogP_MW)

    MatF = Matrix(Df_LogP_MW_no_missing)
    
    PapersF::Vector{Int32} = Vector{Int32}(ones(length(MatF[:,1])) * Paper)

    
    if MA == 0
        MAC = ["Orbitrap"]
        MAF = repeat(MAC, length(MatF[:,1]))
    elseif MA == 1
        MAC = ["Q-TOF"]
        MAF = repeat(MAC, length(MatF[:,1]))
    end

    CVs = [CV]
    CVsF = repeat(CVs, length(MatF[:,1]) )

    if Polarity == 0
        PolC = ["Both, separate"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 1
        PolC = ["Only positive"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 2
        PolC = ["Only negative"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 3
        PolC = ["Both, unclear"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 4
        PolC = ["Both, simultaneous"]
        PolF = repeat(PolC, length(MatF[:,1]))
    end

    if Acq == 0
        AcqC = ["DDA"]
        AcqF = repeat(AcqC, length(MatF[:,1]))
    elseif Acq == 1
        AcqC= ["DIA"]
        AcqF = repeat(AcqC, length(MatF[:,1]))
    elseif Acq == 2
        AcqC= ["both"]
        AcqF = repeat(AcqC, length(MatF[:,1]))
    elseif Acq == 3
        AcqC = ["Not reported"]
        AcqF = repeat(AcqC, length(MatF[:,1]))    
    end

    if MS == 0
        MSC = ["Small/In-house library"]
        MSF = repeat(MSC, length(MatF[:,1]))
    elseif MS == 1
        MSC = ["Medium"]
        MSF = repeat(MSC, length(MatF[:,1]))
    elseif MS == 2
        MSC = ["Large"]
        MSF = repeat(MSC, length(MatF[:,1]))
    elseif MS == 3
        MSC = ["Not Reported"]
        MSF = repeat(MSC, length(MatF[:,1]))    
    end
    Final_Mat = [PapersF MatF[:,1] MatF[:,2] MatF[:,3] MAF CVsF PolF AcqF MSF]
return Final_Mat
end
#Function to convert SMILES into compound descriptors, similar to the previous function
function logP_MW_smiles(Paper, Compounds, MA, CV, Polarity, Acq, MS)
    cids::Vector{Int32} = zeros(length(Compounds))
    for i = 1:length(cids)
        cids[i] = my_get_cid(smiles=Compounds[i])
        @show i
    end
    
    cidss::Vector{Int32} = trunc.(Int32, cids)

    cidssf::Vector{Int32} = filter!(e->e≠0,cidss)
    
    Df_LogP_MW = CSV.File(get_for_cids(cidssf; properties="MolecularWeight,XLogP", output="CSV")) |> DataFrame
    Df_LogP_MW_no_missing = dropmissing(Df_LogP_MW)

    MatF = Matrix(Df_LogP_MW_no_missing)
    
    PapersF::Vector{Int32} = Vector{Int32}(ones(length(MatF[:,1])) * Paper)

   
    if MA == 0
        MAC = ["Orbitrap"]
        MAF = repeat(MAC, length(MatF[:,1]))
    elseif MA == 1
        MAC = ["Q-TOF"]
        MAF = repeat(MAC, length(MatF[:,1]))
    end

    CVs = [CV]
    CVsF = repeat(CVs, length(MatF[:,1]) )

    if Polarity == 0
        PolC = ["Both, separate"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 1
        PolC = ["Only positive"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 2
        PolC = ["Only negative"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 3
        PolC = ["Both, unclear"]
        PolF = repeat(PolC, length(MatF[:,1]))
    elseif Polarity == 4
        PolC = ["Both, simultaneous"]
        PolF = repeat(PolC, length(MatF[:,1]))
    end

    if Acq == 0
        AcqC = ["DDA"]
        AcqF = repeat(AcqC, length(MatF[:,1]))
    elseif Acq == 1
        AcqC= ["DIA"]
        AcqF = repeat(AcqC, length(MatF[:,1]))
    elseif Acq == 2
        AcqC= ["both"]
        AcqF = repeat(AcqC, length(MatF[:,1]))
    elseif Acq == 3
        AcqC = ["Not reported"]
        AcqF = repeat(AcqC, length(MatF[:,1]))    
    end

    if MS == 0
        MSC = ["Small/In-house library"]
        MSF = repeat(MSC, length(MatF[:,1]))
    elseif MS == 1
        MSC = ["Medium"]
        MSF = repeat(MSC, length(MatF[:,1]))
    elseif MS == 2
        MSC = ["Large"]
        MSF = repeat(MSC, length(MatF[:,1]))
    elseif MS == 3
        MSC = ["Not Reported"]
        MSF = repeat(MSC, length(MatF[:,1]))    
    end
    Final_Mat = [PapersF MatF[:,1] MatF[:,2] MatF[:,3] MAF CVsF PolF AcqF MSF]

return Final_Mat
end



#To save any figures:
#Just change the name for each new figure and change the file destination to your PC, it will save the latest figure displayed

savefig("C:\\Users\\tehul\\OneDrive - UvA\\Final figures for paper\\Final 3D PCA NORMAN and CEC's.pdf")
##############################################################################################################################################################################
#Start by loading in the data from all papers

All_papers = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Final files for publication\\All papers.csv", DataFrame)
#In the provided data, the compound names/smiles have already been converted to CID's using the functions from above

cids = Vector{Int}((All_papers[:,2]))

#Just as an example, we can use the CID's to obtain a lot of different properties from PubChem
#If it returns an IO error run it once more and it should work

Properties = CSV.File(get_for_cids(cids; properties="MolecularWeight,XlogP,ExactMass", output="CSV")) |> DataFrame

#Using the dataset from the review we can plot in various manners:

#The MW and logP Plot 
unique_compounds = length(unique(cids))
scatter(All_papers[:,3], All_papers[:,4],  xlims = (0,1400), ylims = (-13,30), palette = :default,
        size = (1280, 720), grid = false, markerstrokewidth = 0.75, xlabel = "Molecular weight", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, label = "Detected structures n = $unique_compounds")

##If you would like to overlay the compounds from the papers onto SusDat, just scatter SusDat first and then the compounds 
#from papers but change the function to "scatter!".


#Grouping can be added to cluster based on methods (based on mass analyzer for example)
#Just change n in group = All_papers[:,n] to cluster based on the other parameters where n = 1:9
scatter(All_papers[:,3], All_papers[:,4], group = All_papers[:,5],  xlims = (0,1400), ylims = (-13,30), palette = :default,
        size = (1280, 720), grid = false, markerstrokewidth = 0.75, xlabel = "Molecular weight", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft)

#The same can be done with the NORMAN susdat database
susdat = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Final files for publication\\NORMAN susdat all chemicals.csv", DataFrame)

sus = filter!(x->x≠"NA",susdat[:,1]) #Certain chemicals had no CID and instead had an "NA" and were thus omitted
sus_2 = filter!(x->x≠"13923912, 11072636",sus) #These compounds returned errors
sus_3 = filter!(x->x≠"15538969;381122354",sus_2) #Also these
final_sus = parse.(Float64, sus_3)
susdat_cids = trunc.(Int64, final_sus)
susdat_final = CSV.File(get_for_cids(susdat_cids; properties="MolecularWeight,XLogP,ExactMass", output="CSV")) |> DataFrame
susdat_final = dropmissing(susdat_final) #There are some missing values for properties, so we remove these compounds
n_norman = length(susdat_final[:,1])
scatter(susdat_final[:,2], susdat_final[:,3],  xlims = (0,1400), ylims = (-13,30), palette = :default,
        size = (1280, 720), grid = false, markerstrokewidth = 0.75, xlabel = "Molecular weight", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, label = "NORMAN SusDat chemicals n = $n_norman", xtickfont=font(15), 
        ytickfont=font(15), 
        guidefont=font(13), 
        legendfont=font(13))



##############################################################################################################################################################################################
#The next section was used to obtain the compounds class for each chemical by inputting the InChIKeys into http://classyfire.wishartlab.com/.
#Only 1000 compounds can be run at once so 3 runs of 759 each were done. Load in the 3 parts:

Part_A = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Final files for publication\\Classes - Part A.csv", DataFrame)
Part_B = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Final files for publication\\Classes - Part B.csv", DataFrame)
Part_C = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Final files for publication\\Classes - Part C.csv", DataFrame)

#The following function is needed to obtain only the Superclass and ignore all the other information obtained from classyfire
function get_compound_class(Data)
    class_dict = Dict{String, Union{Missing, String}}()
    # Iterate over the rows of the dataframe to extract the compound classes
    for i in 1:nrow(Data)
        # Get the compound ID from the first column of the row
        cid = split(Data[i, :CompoundID], "-")[2]
        
        # Check if the "ClassifiedResults" column of the row contains the string "Class"
        if contains(Data[i, :ClassifiedResults], "Superclass:")
            # Extract the compound class from the "ClassifiedResults" column
            class = split(Data[i, :ClassifiedResults], ": ")[2]
            # Store the compound class in the dictionary
            class_dict[cid] = class
        end
    end

    # Create empty vector to store the compound classes
    classes = Vector{Union{Missing, String}}(undef, nrow(Data))

    # Iterate over the rows of the dataframe to fill in the compound classes
    for i in 1:nrow(Data)
        # Get the compound ID from the first column of the row
        cid = split(Data[i, :CompoundID], "-")[2]
        
        # Check if the compound ID is in the dictionary
        if haskey(class_dict, cid)
            # Check if the row is relevant
            if contains(Data[i, :ClassifiedResults], "Superclass:")
                # Extract the compound class from the "ClassifiedResults" column
                class = split(Data[i, :ClassifiedResults], ": ")[2]
                # Store the compound class in the vector
                classes[i] = class
            else
                classes[i] = "Row not relevant"
            end
        else
            classes[i] = "Not found"
        end
    end
    classes = filter(s -> s != "Row not relevant", classes)
    n = length(classes)
    count = 0
    i = 1
    while i <= n
        if classes[i] == "Not found"
            if count > 0
                deleteat!(classes, i)
                n -= 1
                i -= 1
            end
            count += 1
        else
            count = 0
        end
        i += 1
    end
    return classes
end


#This extracts only the superclass from the info by selecting only the relevant rows
classes_a = get_compound_class(Part_A)
classes_b = get_compound_class(Part_B)
classes_c = get_compound_class(Part_C)

#Now combine all 3 again and we get back to our original 2277 compounds (classes) 
Compound_classes = [classes_a; classes_b; classes_c]

#To make a sorted histogram of Strings in Julia use the following function:
function string_histogram(data)
    counts = Dict{String, Int}()
    for d in data
        counts[d] = get(counts, d, 0) + 1
    end
    sorted_counts = sort(collect(counts), by=x->x[2], rev=true)
    xticks = [x[1] for x in sorted_counts]
    bar(xticks, [x[2] for x in sorted_counts], legend=false, ylabel="Frequency", xrotation=40,
     size = (1280,720), grid = false, left_margin = 7Plots.mm, bottom_margin = 45Plots.mm, right_margin = 5Plots.mm, xtickfontsize = 13, ylabelfontsize = 13,)
end

#Plot the histogram of compounds classes using the function:
string_histogram(Compound_classes)

###########################################################################################################################################################
#This section was used to make the PCA 

#Start by obtaining the mass deffects for the NORMAN database and the detected CEC'S
#This is the function used:
er_masses = [27.9949,46.9689,26.003,43.972,30.9984,13.0078]
function ER_calc_ext(mz_values,m_ru)
    ER=Array{Any}(undef,length(mz_values))
    for i=1:length(mz_values)
        KM=mz_values[i].*(round.(m_ru)./m_ru) # Kendrick mass
        ER[i]=round.(round.(KM) .- KM ; digits=3)
    end
    return ER
end
#NORMAN
susdat_exact_masses = susdat_final[:,4]
K_def = ER_calc_ext(susdat_exact_masses, er_masses)
K_susdat = Matrix{Float64}(zeros(length(susdat_exact_masses),6))
for j = 1:length(K_susdat[1,:])
    for i = 1:length(susdat_exact_masses)
        K_susdat[i,j] = K_def[i][j]
    end
end
K_susdat
susdat
susdat_PCA = [susdat_final[:,2] susdat_final[:,3] K_susdat]

#Detected CEC'S
masses = Properties[:,4]

K_def = ER_calc_ext(masses, er_masses)

K_CEC = Matrix{Float64}(zeros(length(masses),6))
for j = 1:length(K_CEC[1,:])
    for i = 1:length(masses)
        K_CEC[i,j] = K_def[i][j]
    end
end

#Combine the EMD's, MW's and XlogP3 for NORMAN and detected CEC to build one PCA model
PCA_data = [All_papers[:,3] All_papers[:,4] K_CEC]

#Now we can perform PCA
using ScikitLearn
@sk_import decomposition: PCA

PCA_norman_and_papers = [susdat_PCA; PCA_data]

X = PCA_norman_and_papers
#Mean center and scale since variables all have different magnitudes
X = X .- mean(X,dims = 1)
X = X ./ std(X,dims=1)
#setup PCA model
pca = PCA(n_components = size(X,2))			#define parameters of PCA model 
pca.fit(X)				
#Variance explained
pca.explained_variance_ratio_
scatter(cumsum(pca.explained_variance_ratio_).*100, size = (1280,720))
plot!(cumsum(pca.explained_variance_ratio_).*100, size = (1280,720), legend = false, grid = false, xlabel = "Number of PC's",
        ylabel = "Explained Variance (%)", left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,
        xlims = (1,8))

#Select an appropiate number of principle components. In our case 3:
pca = PCA(n_components = 3)
pca.fit(X)


#plot Loadings
loadings = pca.components_		#loadings of the model
PCs = ["PC1", "PC2", "PC3"]
var_names = ["MolecularWeight", "logP", "CO","CN","CCl","CS","CF","CH"]

#Plot the loadings of each variable for each PC in a grouped bar plot:
groupedbar(loadings, xticks=(1:3, PCs), bar_position =:group, grid = false, label = ["MolecularWeight" "logP" "CO" "CN" "CCl" "CS" "CF" "CH"],
            size = (1280,720),
            left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,
            ylabel = "Loading")

#This will calculate the scores for each compound in each PC direction:
scores = pca.fit_transform(X)

#Then use these scores to plot:
#Plotting NORMAN data
scatter(scores[1:55793,1], scores[1:55793,2], scores[1:55793,3],
size = (1280, 720), markerstrokewidth = 0.5,
xlabel = "PC 1", ylabel = "PC 2", alpha = 0.75, label = "NORMAN SusDat n = $n_norman",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, xtickfont=font(15), 
ytickfont=font(15), 
guidefont=font(13), 
legendfont=font(13) 
)

#Plotting detected CEC's
scatter!(scores[55794:end,1], scores[55794:end,2], scores[55794:end,3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1, label = "Collected structures n = $unique_compounds",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
camera = (45,45))

#Variance explained used in our model
sum(pca.explained_variance_ratio_[1:3])


#Plotting the Organohalogens only
Organohalogens = findall(x-> x=="Organohalogen compounds", Compound_classes)
n_organohalogens = length(Organohalogens)

scatter!(scores[(Organohalogens.+55793),1], scores[(Organohalogens.+55793),2], scores[(Organohalogens.+55793),3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1, label = "Organohalogen compounds n = $n_organohalogens",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
camera = (45,45))

#Plotting the Organohalogens and Organic acids and derivatives

Organic_acids = findall(x-> x=="Organic acids and derivatives", Compound_classes)
Organohalogens_and_organic_acids = [Organohalogens;Organic_acids]
n_organohalogens_and_organic_acids =length(Organohalogens_and_organic_acids)


scatter!(scores[(Organohalogens_and_organic_acids.+55793),1], scores[(Organohalogens_and_organic_acids.+55793),2], scores[(Organohalogens_and_organic_acids.+55793),3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1, label = "Organohalogens and organic acids and derivatives n = $n_organohalogens_and_organic_acids",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
camera = (45,45))

#This plot is only the Detected structures grouped by the classes obtained in classyfire

scatter(scores[55794:end,1], scores[55794:end,2], scores[55794:end,3], group = Compound_classes,
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
camera = (45,45))


#############################################################################
#Plotting the distributions of MW and XLogP3
#MW:
histogram(susdat_final[:,"MolecularWeight"], bins = 275, label = "NORMAN SusDat n = $n_norman", xlims = (0, 1400), grid = false,
xlabel = "Molecular Weight", bottom_margin = 5Plots.mm, title = "a)", legendfont = 10)
histogram!(All_papers[:,"MW"], label = "Collected stuctures n = $unique_compounds")

#XLogP3:
histogram(susdat_final[:,"XLogP"], bins = 275, label = "NORMAN SusDat n = $n_norman", xlims = (-13, 30), grid = false,
xlabel = "XlogP3", bottom_margin = 5Plots.mm, title = "b)", legendfont = 10)
histogram!(All_papers[:,"logP"], label = "Collected stuctures n = $unique_compounds")

#################################################################################################################################################################