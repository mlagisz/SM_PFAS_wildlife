library(rotl)
library(ape)
library(dplyr)
library(rlist)

## Aims
#For the manually extracted species names (Latin and common names) collate phylogenetic information (tree)

#Load species list
species_all <- read.csv("./data/Species_clean.csv")
#str(species_all)
names(species_all)  

#check and clean:
table(species_all$Species_higher_taxon) #clean some records:
# species_all$Species_higher_taxon[species_all$Species_higher_taxon == "\tBranchiopoda"] <- "Branchiopoda"
# species_all$Species_higher_taxon[species_all$Species_higher_taxon == "\tMalacostraca"] <- "Malacostraca"
# species_all <- subset(species_all, Scientific_name != "") #remove row with missing species and higher taxon names
# species_all$Scientific_name <- trimws(species_all$Scientific_name) #remove leading and trailing whitespaces
# write_csv(species_all, "./data/Species_clean.csv")

mylist <- unique(species_all$Scientific_name)
length(mylist) #1047 unique Latin names (species or higher)
#table(is.na(mylist)) #no NA
#mylist <- mylist[mylist != ""] #remove missing name

taxa <- tnrs_match_names(mylist, context_name = "Animals") #find iTOL records
#View(taxa)

#check flags = potentially problematic taxon names
table(taxa$flags) #925 have no flags = should be ok
#flags such as "hidden", "infraspecific" probably ok; see https://github.com/OpenTreeOfLife/reference-taxonomy/wiki/Taxon-flags
#"sibling higher" can be ignored (93), it is a warning, but other are lethal

table(taxa$number_matches) #mostly single matches (982)
table(taxa$is_synonym) #mostly original name, only 88 used synonyms

#keep only records without any lethal flags
mylist2 <- taxa[taxa$flags == "" | taxa$flags == "sibling_higher" | taxa$flags == "infraspecific", "search_string"]
length(mylist2) #1021 out of 1047 names

#run search again using cleaned taxa names list
taxa2 <- tnrs_match_names(mylist2, context_name = "Animals") #find iTOL records

#retrieve a phylogenetic tree for this list of taxa names
mytree <- tol_induced_subtree(ott_ids = taxa2$ott_id, label_format = "name") #extract subtree of taxa2
#plot(mytree, cex=.8, label.offset =.1, no.margin = TRUE)
names(mytree)

#mytree$tip.label #clean _(species_in_domain_Eukaryota),_(genus_in_Opisthokonta), mrcaott13615ott65415, mrcaott41495ott150711, mrcaott12778ott32449

## clean tip labels
mytree$tip.label <- gsub("\\(.*", "", mytree$tip.label) #remove comments
tip <- c("mrcaott13615ott65415", "mrcaott41495ott150711", "mrcaott12778ott32449", "mrcaott235665ott306513") #labels of internal nodes at tips
mytree <- drop.tip(mytree, tip) #remove tips with the odd names above (not sure which taxa they represent - likely just internal nodes, i.e. higher taxa ranks)
mytree$tip.label <- gsub("_"," ", mytree$tip.label) #get rid of the underscores
#mytree$tip.label <- gsub("Capsicum frutescens", "Capsicum annum", mytree$tip.label) #if needed to replace labels
str(mytree)

#summarise how many names are matching and how many changed:

length(mytree$tip.label)
length(intersect(mytree$tip.label, mylist)) #812 names matching exactly to the original species list
length(setdiff(mytree$tip.label, mylist)) #92 names on the tree but not on the original species list
length(setdiff(mylist, mytree$tip.label)) #235 names on the original species list  but not on the tree - many higher taxonomic levels
length(intersect(mytree$tip.label, taxa2$unique_name)) #896 matching exactly to the adjusted unique species list from iTOL
896-812 #84 names substituted by iTOL
table(tolower(taxa2$search_string) == tolower(taxa2$unique_name)) #123 names not matching

#Make a table of how many words each tip label has:
table((lengths(gregexpr("\\W+", mytree$tip.label)) + 1)) #only 5 records with 3-word- species names (subspecies)

write.tree(mytree, file = "./data/phylogenetic_tree.tre") #save the tree of 904 taxa (mostly species)
# mytree <- read.tree(file = "./data/phylogenetic_tree.tre") #if you need to read it
