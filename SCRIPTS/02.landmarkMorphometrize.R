# Conduct appropriate morphometric analyses of previously
# imported landmark data.
## S.B. Robeson
## last updated 20/03/2020


############################

## Here it seems necessary to choose individual datasets to combine.
## If there is a more efficient way to do this that's equally short,
## I would like to know it. However, with small numbers that need
## combining, I have stuck to doing the job manually, as clunky as
## it may look and feel.

# Create regional groups:
allChicago <- abind(Chicago_AN,Chicago_AS)
allVermont <- abind(Vermont_AN,Vermont_AS)
allVirginia <- abind(Virginia_AN,Virginia_AS)

# Create morphospecies groups:
allANgroup <- abind(Vermont_AN,Chicago_AN,Virginia_AN)
allASgroup <- abind(Vermont_AS,Chicago_AS,Virginia_AS)

# Create combined group:
allTotalGroup <- abind(allANgroup,allASgroup)


# Subject groups to PGPA:
## Partial Generalized Procrustes Analysis involves Procrustes superimposition,
## including transforming the input specimen landmark data through a combination
## of translation, rotation, and rescaling (to a unit centroid size) in order to
## superimpose all individual specimens for ready analysis.
## The output of this analysis is graphed in the next script.

allChicago_GPA <- allChicago %>%
  aligne() %>%
  pgpa()
allVermont_GPA <- allVermont %>%
  aligne() %>%
  pgpa()
allVirginia_GPA <- allVirginia %>%
  aligne() %>%
  pgpa()

allAN_GPA <- allANgroup %>%
  aligne() %>%
  pgpa()
allAS_GPA <- allASgroup %>%
  aligne() %>%
  pgpa()
allTotal_GPA <- allTotalGroup %>%
  aligne() %>%
  pgpa()


# Procrustes ANOVA -- statistical extension of Goodall's F-test,
# here an analogue of two-way ANOVA testing association between
# the Procrustes distance and the categorical variables i.e. the
# region and species of samples.
## Generates tables, yields a single P-value... interpretation-wise,
## just about similar to your garden-variety ANOVA.

# for Acer nigrum across regions:
allAN_ProcAN <- procD.lm(allAN_GPA$rotated~ANgrpList) # Procrustes ANOVA with permutation
summary(allAN_ProcAN)

# for Acer saccharum across regions:
allAS_ProcAN <- procD.lm(allAS_GPA$rotated~ASgrpList)
summary(allAS_ProcAN)

# for both species across regions:
allTotal_ProcAN <- procD.lm(allTotal_GPA$rotated~totalGrpList)
summary(allTotal_ProcAN)


# Pairwise post-hoc test applied to results of Procrustes ANOVA:
## This tests the significance of all the pairs of interactions
## from the Procrustes ANOVA. Unfortunately, it generates something
## of class "pairwise" which is very difficult to work with. You can
## view the output with the summary() command, but I haven't figured
## out how to transform it or anything. Due to this problem (maybe
## there's a quick fix that I don't know of?) I've manually reported
## the Bonferroni-corrected p-values (multiple comparisons) below, but
## I recognize that it's not an elegant solution in any sense.

# apply the pairwise post-hoc:
pairTested <- pairwise(allTotal_ProcAN,covariate = NULL,groups = totalGrpList)
summary(pairTested) # make sure to run this on its own! sourcing this script will not display the table

adjustedP <- p.adjust(c(0.001,0.108,0.001,0.012,
                        0.001,0.001,0.270,0.001,
                        0.035,0.001,0.001,0.001,
                        0.001,0.005,0.001),
                      method = "bonferroni")

adjustedP # left-to-right is the same as pairTested summary table's top-to-bottom
          # yes, this is incredibly clunky

############################