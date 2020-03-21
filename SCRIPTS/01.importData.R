# A relatively simple script to import and parse landmark data
# (geometric morphometric information) for further analysis.
## S.B. Robeson
## last updated 20/03/2020


# Make sure to call upon the appropriate libraries:
## Again, a script is provided with functions originally taken from
## Claude (2008), since I do not know of a package for these. Provide
## credit to the original...

source("../SCRIPTS/00.Claude2008fxn.R")
library(tidyverse)
library(geomorph)
library(shapes)
library(abind)


# Locate and extract names of relevant files:
## These operations find all .tps files in the folder of data,
## and list the names to pass to the next command (reading them).
morphFiles <- list.files(path = "../DATA/Landmark//",
                          pattern = ".tps",
                          full.names = FALSE,
                          recursive = FALSE)

morphFilenames <- strsplit(morphFiles,
                           ".",
                           fixed = TRUE) %>%
  sapply(., head, 1)



# Create factors that serve as lists of group membership:
## This is another manual solution for similar reasons to the
## previous combining of group data. I hope I will be able to
## return to this and revise accordingly 
VTANlist <- matrix(rep("Vermont_AN",48),48,1)
CHANlist <- matrix(rep("Chicago_AN",43),43,1)
VAANlist <- matrix(rep("Virginia_AN",44),44,1)
VTASlist <- matrix(rep("Vermont_AS",44),44,1)
CHASlist <- matrix(rep("Chicago_AS",44),44,1)
VAASlist <- matrix(rep("Virginia_AS",44),44,1)
VTgrpList <- rbind(VTANlist, VTASlist) %>% as_factor()
ChGrpList <- rbind(CHANlist, CHASlist) %>% as_factor()
VAgrpList <- rbind(VAANlist, VAASlist) %>% as_factor()
ANgrpList <- rbind(VTANlist, CHANlist, VAANlist) %>% as_factor()
ASgrpList <- rbind(VTASlist, CHASlist, VAASlist) %>% as_factor()
totalGrpList <- rbind(VTANlist, CHANlist, VAANlist,
                      VTASlist, CHASlist, VAASlist) %>% as_factor()




# Load all landmark data:
## This loops through all .tps files provided and reads them to
## arrays that are assigned the same names as their input files.

for(i in morphFilenames) {
  
  filepath <- file.path("../DATA/Landmark/", 
                        paste(i,
                              ".tps",
                              sep=""))
  
  assign(i, readland.tps(filepath,
                         readcurves = FALSE,
                         negNA = TRUE,
                         specID = "imageID"))
}




# Find and load angular data:
## The below applies to angular measurements the same logic as the
## above set of operations for landmark data.

## NOTE: These data were not used for analysis here!
## A one-way ANOVA might be appropriate, among other analyses. This
## dataset was simply one I did not have time to cover at lab meeting.
## These data also likely need cleaning due to lacking labels...

# anglesFiles <- list.files(path = "../DATA/Linear/",
#                          pattern = ".txt",
#                          full.names = FALSE)

#anglesFilenames <- strsplit(anglesFiles,
#                            ".",
#                            fixed = TRUE) %>%
#  sapply(., head, 1)


#for(i in anglesFilenames) {
#  
#  filepath <- file.path("../DATA/Linear/",
#                        paste(i,
#                              ".txt",
#                              sep=""))
#  
#  assign(i, read_csv(filepath))
#}

############################