#savefile <- "../data/ToolDevelopmentNovember.RData" # Save file for the program output
#load(savefile)

#(!exists("simfile")){
#  simfile <- c(strsplit(list.files(path = "../data/KneeSimulator/271021_rtkn1307f/txts", pattern = "*KP1.txt",
#                                   full.names = T), split = " "))
#}

#if(!exists("tekfile")){
#  tekfile <- c(strsplit(list.files(path = "../data/Tekscan/271021_rtkn1307f/", pattern = "Def_rtkn",
#                                   full.names = T), split = " "))
#}

#if(!exists("results")){
#  results <- vector(mode = "list", length = length(tekfile)) #An expected vector for the results to be stored
#}

#if(length(results) < length(tekfile)){
#  results <- c(results, vector(mode = "list", length = length(tekfile)-length(results))) # appending new list entries onto the existing resuls table from pts 1 and 2, only run once.
#}

#nms <- c("Intact", "OgPos", "Medial", "Posterior", "PostMed", "OgPos_IK", "Medial_IK", "Posterior_IK", "PostMed_IK",
#         "Medial_10mm", "PostMed_10mm", "Meniscectomy", "Medial_10_IK", "PostMed_10_IK", "Meniscectomy_IK")

#if(!all(names(results)==nms)){
#  names(results) <- nms #Naming the results elements according to the analysis order.
#}
