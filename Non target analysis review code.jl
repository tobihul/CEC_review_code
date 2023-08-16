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
function logP_MW_names(Paper, Compounds, MA, CV, Polarity, Acq, MSDIA)

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

savefig("C:\\Users\\tehul\\OneDrive - UvA\\Final figures for paper\\Transparent NORMAN vs detected CECs final.png")
savefig("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Final figures for paper/LogP Histogram final.png")


##############################################################################################################################################################################
#Start by loading in the data from all papers
All_papers = CSV.read("/Users/tobias/Downloads/All papers.csv", DataFrame)
All_papers_CV = CSV.read("/Users/tobias/Downloads/All papers_CV.csv", DataFrame)
#In the provided data, the compound names/smiles have already been converted to CID's using the functions from above

cids = Vector{Int}((All_papers[:,2]))

#Just as an example, we can use the CID's to obtain a lot of different properties from PubChem
#If it returns an IO error run it once more and it should work

Properties = CSV.File(get_for_cids(cids; properties="InChIKey,CanonicalSMILES", output="CSV")) |> DataFrame

#Using the dataset from the review we can plot in various manners:

#The MW and logP Plot 
unique_compounds = length(unique(cids))
p_together = scatter!(All_papers[:,3], All_papers[:,4],  xlims = (0,3500), ylims = (-35,75), palette = :default,
        size = (1280, 720), grid = false, xlabel = "Molecular weight", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, label = "Detected structures n = $unique_compounds",
        legendfont=font(13), markersize = 5, markerstrokewidth = 0.75, dpi = 300,xtickfont=font(13), ztickfont = font(13),
        ytickfont=font(13), 
        guidefont=font(20))
p_CEC = scatter(All_papers[:,3], All_papers[:,4],  xlims = (0,1400), ylims = (-13,30), palette = :default,
        size = (1280, 720), grid = false, xlabel = "Molecular weight", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, label = "Detected structures n = $unique_compounds",xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0.75, dpi = 300)
##If you would like to overlay the compounds from the papers onto SusDat, just scatter SusDat first and then the compounds 
#from papers but change the function to "scatter!".

#For plotting the MW and XlogP3 based on the column volumes used in each method
Column_volumes_plot = scatter(All_papers_CV[:,3], All_papers_CV[:,4], marker_z = All_papers_CV[:,6], xlims = (0,1400), 
                ylims = (-13,30), color = :plasma, size = (1280, 720), grid = false, markerstrokewidth = 0, label = false,
                xlabel = "Molecular weight (Da)", ylabel = "XLogP3",  zlabel= "Column volumes",
                left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, xtickfont=font(13), 
                ytickfont=font(13), 
                guidefont=font(20), 
                legendfont=font(13), markersize = 5, dpi = 300)



#Grouping can be added to cluster based on methods (based on mass analyzer for example)
#Mass analyzer
indices_orbitrap = findall(row -> occursin("Orbitrap", row.Mass_analyzer), eachrow(All_papers))
indices_QTOF = findall(row -> occursin("Q-TOF", row.Mass_analyzer), eachrow(All_papers))
p_orbi = scatter(All_papers[indices_orbitrap,3], All_papers[indices_orbitrap,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Orbitrap n = $(length(indices_orbitrap))", dpi = 300)
p_qtie = scatter(All_papers[indices_QTOF,3], All_papers[indices_QTOF,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), markershape = :utriangle,
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Q-TOF n = $(length(indices_QTOF))", dpi = 300, c =:darkorange1)    

scatter(All_papers[indices_orbitrap,3], All_papers[indices_orbitrap,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Orbitrap n = $(length(indices_orbitrap))", dpi = 300)
p_together = scatter!(All_papers[indices_QTOF,3], All_papers[indices_QTOF,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), c = :darkorange1, markershape = :utriangle ,
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Q-TOF n = $(length(indices_QTOF))", dpi = 300) 
        l = @layout [
            a{0.7w} [
                grid(1, 1)
                b{0.5h}
            ]
        ]
        
        plot(p_together, p_orbi, p_qtie, layout = l)
        
# Acquisition mode

indices_DDA = findall(row -> occursin("DDA", row.Acquisition), eachrow(All_papers))
indices_DIA = findall(row -> occursin("DIA", row.Acquisition), eachrow(All_papers))
indices_Not_reported = findall(row -> occursin("Not reported", row.Acquisition), eachrow(All_papers))
indices_both = findall(row -> occursin("both", row.Acquisition), eachrow(All_papers))

p_DDA = scatter(All_papers[indices_DDA,3], All_papers[indices_DDA,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "DDA n = $(length(indices_DDA))", dpi = 300, xticks = (0:500:1400))
p_DIA = scatter(All_papers[indices_DIA,3], All_papers[indices_DIA,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "DIA n = $(length(indices_DIA))", markershape = :utriangle, c = :chocolate1, dpi = 300, xticks = (0:500:1400))   

p_both = scatter(All_papers[indices_both,3], All_papers[indices_both,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "both n = $(length(indices_both))", markershape = :diamond, c = :chartreuse3, dpi = 300, xticks = (0:500:1400))
p_not_reported = scatter(All_papers[indices_Not_reported,3], All_papers[indices_Not_reported,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Not reported n = $(length(indices_Not_reported))", markershape = :rect, c = :hotpink, dpi = 300, xticks = (0:500:1400))   

scatter(All_papers[indices_DDA,3], All_papers[indices_DDA,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "DDA n = $(length(indices_DDA))", dpi = 300)
scatter!(All_papers[indices_DIA,3], All_papers[indices_DIA,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "DIA n = $(length(indices_DIA))", markershape = :utriangle, c = :chocolate1, dpi = 300)   

scatter!(All_papers[indices_both,3], All_papers[indices_both,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "both n = $(length(indices_both))", markershape = :diamond, c = :chartreuse3, dpi = 300)
p_together = scatter!(All_papers[indices_Not_reported,3], All_papers[indices_Not_reported,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Not reported n = $(length(indices_Not_reported))", markershape = :rect, c = :hotpink, dpi = 300)  
        
        l = @layout [
            a{0.5w} [
                grid(2, 2)
                

              
                
            ]
        ]
        
        plot(p_together, p_DDA, p_DIA, p_both, p_not_reported, layout = l)

#Polarity
indices_Separate = findall(row -> occursin("Both, separate", row.Polarity), eachrow(All_papers))
indices_Simultaneous = findall(row -> occursin("Both, simultaneous", row.Polarity), eachrow(All_papers))
indices_Unclear = findall(row -> occursin("Both, unclear", row.Polarity), eachrow(All_papers))
indices_negative = findall(row -> occursin("Only negative", row.Polarity), eachrow(All_papers))
indices_positive = findall(row -> occursin("Only positive", row.Polarity), eachrow(All_papers))

p_1 = scatter(All_papers[indices_Separate,3], All_papers[indices_Separate,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Both, separate n = $(length(indices_Separate))", dpi = 300)
p_3 = scatter(All_papers[indices_Unclear,3], All_papers[indices_Unclear,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Both, unclear = $(length(indices_Unclear))", markershape = :utriangle, c = :chocolate1, dpi = 300)
p_2 = scatter(All_papers[indices_positive,3], All_papers[indices_positive,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Only positive n = $(length(indices_positive))", markershape = :diamond, c = :chartreuse3, dpi = 300)   
p_4 = scatter(All_papers[indices_negative,3], All_papers[indices_negative,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Only negative n = $(length(indices_negative))", markershape = :rect, c = :hotpink, dpi = 300)   

p_5 = scatter(All_papers[indices_Simultaneous,3], All_papers[indices_Simultaneous,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Both, simultaneous n = $(length(indices_Simultaneous))", markershape = :star5, c = :red, dpi = 300)   

scatter(All_papers[indices_Separate,3], All_papers[indices_Separate,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Both, separate n = $(length(indices_Separate))", dpi = 300)
        scatter!(All_papers[indices_positive,3], All_papers[indices_positive,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Only positive n = $(length(indices_positive))", markershape = :diamond, c = :chartreuse3, dpi = 300)
scatter!(All_papers[indices_Unclear,3], All_papers[indices_Unclear,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Both, unclear = $(length(indices_Unclear))", markershape = :utriangle, c = :chocolate1, dpi = 300)   
scatter!(All_papers[indices_negative,3], All_papers[indices_negative,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Only negative n = $(length(indices_negative))", markershape = :rect, c = :hotpink, dpi = 300)   

p_together = scatter!(All_papers[indices_Simultaneous,3], All_papers[indices_Simultaneous,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Both, simultaneous n = $(length(indices_Simultaneous))", markershape = :star5, c = :red, dpi = 300)   

        l = @layout [
            a{0.6w} [
                grid(5, 1)
                

              
                
            ]
        ]
        
        plot(p_together, p_1, p_2, p_3, p_4, p_5, layout = l)
#Database
indices_Large = findall(row -> occursin("Large", row.Database), eachrow(All_papers))
indices_Not = findall(row -> occursin("Not Reported", row.Database), eachrow(All_papers))
indices_Small = findall(row -> occursin("Small/In-house library", row.Database), eachrow(All_papers))

p_large = scatter(All_papers[indices_Large,3], All_papers[indices_Large,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Both, Large n = $(length(indices_Large))", dpi = 300)
p_not = scatter(All_papers[indices_Not,3], All_papers[indices_Not,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Not Reported = $(length(indices_Not))", markershape = :utriangle, c = :chocolate1, dpi = 300)
p_small = scatter(All_papers[indices_Small,3], All_papers[indices_Small,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9), markersize = 2.5, markerstrokewidth = 0, label = "Small/In-house library n = $(length(indices_Small))", markershape = :diamond, c = :chartreuse3, dpi = 300)
scatter(All_papers[indices_Large,3], All_papers[indices_Large,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Both, Large n = $(length(indices_Large))", dpi = 300)
scatter!(All_papers[indices_Not,3], All_papers[indices_Not,4], xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Not Reported = $(length(indices_Not))", markershape = :utriangle, c = :chocolate1, dpi = 300)
p_together = scatter!(All_papers[indices_Small,3], All_papers[indices_Small,4],  xlims = (0,1400), ylims = (-13,30),
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, xtickfont=font(13), 
        ytickfont=font(13), 
        guidefont=font(20), 
        legendfont=font(13), markersize = 5, markerstrokewidth = 0, label = "Small/In-house library n = $(length(indices_Small))", markershape = :diamond, c = :chartreuse3, dpi = 300)

        l = @layout [
            a{0.7w} [
                grid(3, 1)
                

              
                
            ]
        ]
        
        plot(p_together, p_large, p_not, p_small, layout = l)
#The same can be done with the NORMAN susdat database
susdat = CSV.read("C:\\Users\\tehul\\OneDrive - UvA\\Final files for publication\\NORMAN susdat all chemicals.csv", DataFrame)
susdat = CSV.read("/Users/tobias/Library/CloudStorage/OneDrive-UvA/Final files for publication/NORMAN susdat all chemicals.csv", DataFrame)

sus = filter!(x->x≠"NA",susdat[:,1]) #Certain chemicals had no CID and instead had an "NA" and were thus omitted
sus_2 = filter!(x->x≠"13923912, 11072636",sus) #These compounds returned errors
sus_3 = filter!(x->x≠"15538969;381122354",sus_2) #Also these
final_sus = parse.(Float64, sus_3)
susdat_cids = trunc.(Int64, final_sus)
susdat_final = CSV.File(get_for_cids(susdat_cids; properties="MolecularWeight,XLogP,ExactMass", output="CSV")) |> DataFrame
susdat_final = dropmissing(susdat_final) #There are some missing values for properties, so we remove these compounds
n_norman = length(susdat_final[:,1])
scatter(susdat_final[:,2], susdat_final[:,3],  xlims = (0,3500), ylims = (-35,75), size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :topleft, label = "NORMAN SusDat chemicals n = $n_norman", xtickfont=font(13), 
        ytickfont=font(13),
    
        guidefont=font(20), 
        legendfont=font(13),markersize = 5, markerstrokewidth = 0.75, dpi = 300)
p_NORMAN = scatter(susdat_final[:,2], susdat_final[:,3],  xlims = (0,1400), ylims = (-13,30), palette = :default,
        size = (1280, 720), grid = false, xlabel = "Molecular weight (Da)", ylabel = "XLogP3", 
        left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
        legend = :false, label = "NORMAN SusDat chemicals n = $n_norman", xtickfont=font(9), 
        ytickfont=font(9), 
        guidefont=font(9), 
        legendfont=font(9),markersize = 5, markerstrokewidth = 0.75, dpi = 300)

        l = @layout [
            a{0.6w} [
                grid(1, 1)
                b{0.5h}
            ]
        ]
        
        plot(p_together, p_NORMAN, p_CEC, layout = l)

##############################################################################################################################################################################################
#The next section was used to obtain the compounds class for each chemical by inputting the InChIKeys into http://classyfire.wishartlab.com/.
#Only 1000 compounds can be run at once so 3 runs of 856 each were done. Load in the 3 parts:


All_Parts = CSV.read("/Users/tobias/Downloads/PART_ABC.csv", DataFrame)

All_classes = All_Parts[:,4]

#To make a sorted histogram of Strings in Julia use the following function:
function string_histogram(data)
    counts = Dict{String, Int}()
    for d in data
        counts[d] = get(counts, d, 0) + 1
    end
    sorted_counts = sort(collect(counts), by=x->x[2], rev=true)
    xticks = [x[1] for x in sorted_counts]
    bar(xticks, [x[2] for x in sorted_counts], legend=false, ylabel="Frequency", xrotation=40,
     size = (1280,720), grid = false, left_margin = 7Plots.mm, bottom_margin = 45Plots.mm, right_margin = 5Plots.mm, xtickfontsize = 13, dpi = 300
     ,xtickfont=font(13),
     ytickfont=font(13) ,yguidefont = font(25), 
     guidefont=font(20 )
)
end

#Plot the histogram of compounds classes using the function:
string_histogram(All_classes)
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
masses = Properties[:,2]

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
scatter(cumsum(pca.explained_variance_ratio_).*100, size = (1280,720), dpi = 300,left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
legend = :topleft, xtickfont=font(13), 
ytickfont=font(13), 
guidefont=font(20), 
legendfont=font(13), markersize = 5)
plot!(cumsum(pca.explained_variance_ratio_).*100, size = (1280,720), legend = false, grid = false, xlabel = "Number of PC's",
        ylabel = "Explained Variance (%)", left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,
        xlims = (1,8), dpi = 300)

#Select an appropiate number of principle components. In our case 3:
pca = PCA(n_components = 3)
pca.fit(X)


#plot Loadings
loadings = pca.components_		#loadings of the model
PCs = ["PC1", "PC2", "PC3"]
var_names = ["MolecularWeight", "XlogP3", "CO","CN","CCl","CS","CF","CH"]
legend_order = ["MolecularWeight", "XlogP3", "CO", "CN", "CF", "CS", "CCl", "CH"]


#Plot the loadings of each variable for each PC in a grouped bar plot:
b1 = bar(loadings[1,:], group = var_names, palette = :tab10, size = (1000,500),
left_margin = 7Plots.mm, bottom_margin = 7.5Plots.mm, right_margin = 5Plots.mm,
ylabel = "Loading", legendfont = font(7), guidefont = 15, xticks = false, title = "PC 1", legend = :topleft, dpi = 300, ylims = (-0.5,0.7))
b2 = bar(loadings[2,:], group = var_names, palette = :tab10, size = (1000,500),
left_margin = 7Plots.mm, bottom_margin = 7.5Plots.mm, right_margin = 5Plots.mm,
 legendfont = font(11), guidefont = 15, xticks = false, title = "PC 2", legend = false , dpi = 300,ylims = (-0.5,0.7))
b3 = bar(loadings[3,:], group = var_names, palette = :tab10, size = (1000,500),
left_margin = 7Plots.mm, bottom_margin = 7.5Plots.mm, right_margin = 5Plots.mm,
 legendfont = font(5), guidefont = 15, xticks = false, title = "PC 3" , dpi = 300,ylims = (-0.5,0.7),legend_order = legend_order, legend = :false)

plot(b1,b2,b3, layout = (1,3), dpi = 300)
#This will calculate the scores for each compound in each PC direction:
scores = pca.fit_transform(X)

#Then use these scores to plot:
#Plotting NORMAN data
scatter(scores[1:55793,1], scores[1:55793,2], scores[1:55793,3],
size = (1280, 720), markerstrokewidth = 0.75,
alpha = 0.75, label = "NORMAN SusDat n = $n_norman",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm, xtickfont=font(10), 
ytickfont=font(10),
ztickfont=font(10),
guidefont=font(12),
legendfont=font(15), dpi = 300)

#Plotting detected CEC's
scatter!(scores[55794:end,1], scores[55794:end,2], scores[55794:end,3], 
size = (1280, 720), markerstrokewidth = 0.75, legend = true, xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)",
 alpha = 1, label = "Collected structures n = $unique_compounds",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
camera = (45,45),  dpi = 300,xtickfont=font(13), ztickfont = font(13),
ytickfont=font(13) ,
guidefont=font(20), titlefont = font(25))

#Variance explained used in our model
sum(pca.explained_variance_ratio_[1:3])


#Plotting the Organohalogens only
Organohalogens = findall(x-> x=="Organohalogen compounds", All_classes)
n_organohalogens = length(Organohalogens)

scatter!(scores[(Organohalogens.+55793),1], scores[(Organohalogens.+55793),2], scores[(Organohalogens.+55793),3], 
size = (1280, 720), markerstrokewidth = 0.75, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1, label = "Organohalogen compounds n = $n_organohalogens",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
camera = (45,45), dpi = 300,xtickfont=font(13), 
ytickfont=font(13), 
guidefont=font(20))

#Plotting the Organohalogens and Organic acids and derivatives

Organic_acids = findall(x-> x=="Organic acids and derivatives", All_classes)
Organohalogens_and_organic_acids = [Organohalogens;Organic_acids]
n_organohalogens_and_organic_acids =length(Organohalogens_and_organic_acids)


scatter!(scores[(Organohalogens_and_organic_acids.+55793),1], scores[(Organohalogens_and_organic_acids.+55793),2], scores[(Organohalogens_and_organic_acids.+55793),3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1, label = "Organohalogens and organic acids and derivatives n = $n_organohalogens_and_organic_acids",
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
camera = (45,45), dpi = 300,xtickfont=font(13), 
ytickfont=font(13), 
guidefont=font(20))

#This plot is only the Detected structures grouped by the classes obtained in classyfire
indices_1 = findall(x -> x== "Benzenoids", All_classes)
indices_2 = findall(x -> x== "Organoheterocyclic compounds", All_classes)
indices_3 = findall(x -> x== "Organic acids and derivatives", All_classes)
indices_15 = findall(x -> x== "Organic oxygen compounds", All_classes)
indices_4 = findall(x -> x== "Lipids and lipid-like molecules", All_classes)
indices_5 = findall(x -> x== "Organohalogen compounds", All_classes)
indices_6 = findall(x -> x== "Organic nitrogen compounds", All_classes)
indices_7 = findall(x -> x== "Phenylpropanoids and polyketides", All_classes)
indices_8 = findall(x -> x== "Alkaloids and derivatives", All_classes)
indices_9 = findall(x -> x== "Nucleosides, nucleotides, and analogues", All_classes)
indices_10 = findall(x -> x== "Homogeneous non-metal compounds", All_classes)
indices_11 = findall(x -> x== "Organophosphorus compounds", All_classes)
indices_12 = findall(x -> x== "Organosulfur compounds", All_classes)
indices_13 = findall(x -> x== "Hydrocarbons", All_classes)
indices_14 = findall(x -> x== "Lignans, neolignans and related compounds", All_classes)

scatter(scores[(55793 .+indices_1),1], scores[(55793 .+indices_1), 2], scores[(55793 .+indices_1), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Benzenoids n = $(length(indices_1))")

scatter!(scores[(55793 .+indices_2),1], scores[(55793 .+indices_2), 2], scores[(55793 .+indices_2), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Organoheterocyclic compounds n = $(length(indices_2))")

scatter!(scores[(55793 .+indices_3),1], scores[(55793 .+indices_3), 2], scores[(55793 .+indices_3), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Organic acids and derivatives n = $(length(indices_3))")

scatter!(scores[(55793 .+indices_15),1], scores[(55793 .+indices_15), 2], scores[(55793 .+indices_15), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Organic oxygen compounds n = $(length(indices_15))")

scatter!(scores[(55793 .+indices_4),1], scores[(55793 .+indices_4), 2], scores[(55793 .+indices_4), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Lipids and lipid-like molecules n = $(length(indices_4))")

scatter!(scores[(55793 .+indices_5),1], scores[(55793 .+indices_5), 2], scores[(55793 .+indices_5), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Organohalogen compounds n = $(length(indices_5))")

scatter!(scores[(55793 .+indices_6),1], scores[(55793 .+indices_6), 2], scores[(55793 .+indices_6), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Organic nitrogen compounds n = $(length(indices_6))")

scatter!(scores[(55793 .+indices_7),1], scores[(55793 .+indices_7), 2], scores[(55793 .+indices_7), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Phenylpropanoids and polyketides n = $(length(indices_7))")

scatter!(scores[(55793 .+indices_8),1], scores[(55793 .+indices_8), 2], scores[(55793 .+indices_8), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Alkaloids and derivatives n = $(length(indices_8))")

scatter!(scores[(55793 .+indices_9),1], scores[(55793 .+indices_9), 2], scores[(55793 .+indices_9), 3],
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Nucleosides, nucleotides, and analogues n = $(length(indices_9))")

scatter!(scores[(55793 .+indices_10),1], scores[(55793 .+indices_10), 2], scores[(55793 .+indices_10), 3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Homogeneous non-metal compounds n = $(length(indices_10))")

scatter!(scores[(55793 .+indices_11),1], scores[(55793 .+indices_11), 2], scores[(55793 .+indices_11), 3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Organophosphorus compounds n = $(length(indices_11))")

scatter!(scores[(55793 .+indices_12),1], scores[(55793 .+indices_12), 2], scores[(55793 .+indices_12), 3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Organosulfur compounds n = $(length(indices_12))")

scatter!(scores[(55793 .+indices_13),1], scores[(55793 .+indices_13), 2], scores[(55793 .+indices_13), 3], 
size = (1280, 720), markerstrokewidth = 0.5, legend = true,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Hydrocarbons n = $(length(indices_13))")

scatter!(scores[(55793 .+indices_14),1], scores[(55793 .+indices_14), 2], scores[(55793 .+indices_14), 3], 
size = (1280, 720), markerstrokewidth = 0.5,
xlabel = "PC 1(33%)", ylabel = "PC 2(26%)", zlabel = "PC3(15%)", alpha = 1,
left_margin = 7Plots.mm, bottom_margin = 5Plots.mm, right_margin = 5Plots.mm,  
dpi = 300, markersize = 5, label = "Lignans, neolignans and related compounds n = $(length(indices_14))",
xtickfont=font(13),
ztickfont = font(13) ,
        ytickfont=font(13), 
        guidefont=font(20), xlims = (-4,6.5), ylims = (-3.5,12), legendfont = font(10), camera = (45,45)
        )


#############################################################################
#Plotting the distributions of MW and XLogP3
#MW:
histogram(susdat_final[:,"MolecularWeight"], bins = 275, label = "NORMAN SusDat n = $n_norman", xlims = (0, 1400), grid = false,
xlabel = "Molecular Weight (Da)", ylabel = "Frequency", bottom_margin = 5Plots.mm, top_margin = 7.5Plots.mm,  xtickfont=font(10), 
ytickfont=font(10), 
guidefont=font(15), dpi = 300, legendfont = font(10))
histogram!(All_papers[:,"MW"], label = "Collected stuctures n = $unique_compounds", dpi = 300,xtickfont=font(13), ztickfont = font(13),
ytickfont=font(13), 
guidefont=font(20), titlefont = font(25),
)
annotate!(50, 5500, text("a)", :left, 15, :black))

#XLogP3:
histogram(susdat_final[:,"XLogP"], bins = 275, label = "NORMAN SusDat n = $n_norman", xlims = (-13, 30), grid = false,
xlabel = "XlogP3", bottom_margin = 5Plots.mm, top_margin = 7.5Plots.mm, xtickfont=font(10), 
ytickfont=font(10), ylabel = "Frequency",
guidefont=font(15), dpi = 300, legendfont = font(10) )
histogram!(All_papers[:,"logP"], label = "Collected stuctures n = $unique_compounds", dpi = 300,xtickfont=font(13), ztickfont = font(13),
ytickfont=font(13),  
guidefont=font(20), titlefont = font(25))
annotate!(-10, 5780, text("b)", :left, 15))


################ Histogram compounds per paper

All_papers

count_paper = zeros(maximum(All_papers[:,1]))
for i = 1:maximum(All_papers[:,1])
    count_paper[i] = length(findall(x-> x==i, All_papers[:,1]))
end

count_paper

count_paper_df = DataFrame(Identified_compounds = count_paper)

CSV.write("/Users/tobias/Downloads/Compounds per paper.csv", count_paper_df)
num_papers = length(count_paper)
paper_labels = ["Paper $i" for i in 1:num_papers]

# Sort indices based on compound counts in descending order
sorted_indices = sortperm(count_paper, rev=true)

# Sort compound counts based on the sorted indices
sorted_counts = count_paper[sorted_indices]

# Create a sorted bar plot with rotated labels
bar(paper_labels[sorted_indices], sorted_counts, xlabel = "Paper", ylabel = "Identified Compounds",
    legend=false, xticks=(1:num_papers, sorted_indices),
    grid = false, size = (2000,1000), dpi = 300,xtickfont=font(13),
    ytickfont=font(13), yguidefont = font(25),
    guidefont=font(25), titlefont = font(20), left_margin = 12Plots.mm, bottom_margin = 15Plots.mm, 
    right_margin = 5Plots.mm, xtickfontsize = 13)
#################################################################################################################################################################
