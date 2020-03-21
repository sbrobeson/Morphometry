# Conduct appropriate morphometric analyses of previously
# imported landmark data.


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


# Create factors that serve as lists of group membership:
## This is another manual solution for similar reasons to the
## previous combining of group data. I hope I will be able to
## return to this and revise accordingly 
VTANlist<-matrix(rep("Vermont_AN",48),48,1)
CHANlist<-matrix(rep("Chicago_AN",43),43,1)
VRANlist<-matrix(rep("Virginia_AN",44),44,1)
VTASlist<-matrix(rep("Vermont_AS",44),44,1)
CHASlist<-matrix(rep("Chicago_AS",44),44,1)
VRASlist<-matrix(rep("Virginia_AS",44),44,1)
ANgrpList<-as.factor(rbind(VTANlist,CHANlist,VRANlist))
ASgrpList <- as.factor(rbind(VTASlist,CHASlist,VRASlist))
totalGroupList<-as.factor(rbind(VTANlist,CHANlist,VRANlist,VTASlist,CHASlist,VRASlist))

# Procrustes ANOVA of 
allANgroup_sup <- gpagen(allANgroup)
allAN_ProcAN <- procD.lm(allANgroup_sup$coords~ANgrpList) # Procrustes ANOVA with permutation
summary(allAN_ProcAN)

allASgroup_sup <- gpagen(allASgroup)
allAS_ProcAN <- procD.lm(allASgroup_sup$coords~ASgrpList)
summary(allAS_ProcAN)

allTotalGroup2 <- gpagen(allTotalGroup)
allTotal_ProcAN <- procD.lm(allTotalGroup2$coords~totalGroupList)
summary(allTotal_ProcAN)

usage <- allTotalGroup2$Csize %>% names()
usage2 <- (strsplit(usage, split = "Acer ", fixed = TRUE))
usage2 <- usage2[1:264,2]
pairTested <- pairwise(allTotal_ProcAN,covariate = NULL,groups = totalGroupList)
realigned <- p.adjust(c(0.001,0.108,0.001,0.012,0.001,0.001,0.270,0.001,0.035,0.001,0.001,0.001,0.001,0.005,0.001),
                      method = "bonferroni")
# summary(pairTested, test.type = "dist", confidence = 0.95, stat.table = TRUE, p.adjust.methods)
# define adjustment for this if you want to do it in a single line

testmeanshapes(allANgroup,allASgroup,
               resamples = 1000,
               replace = FALSE,
               scale = TRUE)
testmeanshapes(Vermont_AN, Vermont_AS,
               resamples = 1000,
               replace = FALSE,
               scale = TRUE)
testmeanshapes(Chicago_AN, Chicago_AS,
               resamples = 1000,
               replace = FALSE,
               scale = TRUE)



#proc <- gpagen(allTotalGroup)
#coords2d <- two.d.array(proc$coords)
#consensus <- apply(proc$coords, c(1,2), mean)
#consensusvec <- apply(coords2d, 2, mean)
#resids <- t(t(coords2d)-consensusvec)
#Promet <- cov(resids)
#pca.stuff <- svd(Promet)
#eigenvalues <- pca.stuff$d
#eigenvectors <- pca.stuff$u
#scores <- resids%*%eigenvectors
#plot(scores[,1:2],asp=1, pch=20,cex=2)
#points(scores[1:132,], scores[1:132,],col="red",cex=0.5,pch=16)
#points(scores[,1,132:264], scores[,2,132:264],col="blue",cex=0.5,pch=16)
#points(scores[,1,89:132], scores[,2,89:132],col="yellow",cex=0.5,pch=16)

plot(prcomp(allTotalGroup, center = TRUE, scale. = TRUE))

#AN_umang <- cbind(Vermont_AN_ang,Chicago_AN_ang,Virginia_AN_ang)
#AS_umang <- abind(Vermont_AS_ang,Chicago_AS_ang,Virginia_AS_ang)
#aov(y~allTotalGroup, data = )


#pgpa_geomorph<-gpagen(LMdata)
PCAgeom <- plotTangentSpace(allTotalGroup2$coords, groups = totalGroupList, legend = TRUE, warpgrids = TRUE)

geomPCAscores<-PCAgeom$pc.scores

PCAgeomAN <- plotTangentSpace(allANgroup_sup$coords, groups = ANgrpList, legend = TRUE)
title("PCA of Black Maple Procrustes coordinates")

PCAgeomAS <- plotTangentSpace(allASgroup_sup$coords, groups = ASgrpList)
title("PCA of Sugar Maple Procrustes coordinates")


totalGroupListRev <- relevel(totalGroupList, "Virginia_AS")



# PCA of all warp scores (above) gives identical results to PCA of Procrustes coordinates
# in tangent space (orthogonal projection), except that the four zero PCs are also reported in the latter

allAN_PCA_t <- prcomp(allAN_GPA)


# CALCULATING WARP SCORES

# Uniform and nonuniform components (partial warps), Claude 2008 Functions
# Nonuniform scores differ from those generated by IMP because they are calculated in a different way
kp<-dim(allTotalGroup)[2]*dim(allTotalGroup)[1]
n<-dim(allTotalGroup)[3]
uniform<-uniform2D(allTotalGroup) # uniform terms
msh<-uniform$meanshape
Un<-uniform$uniform
X<-t(matrix(uniform$rotated,kp,n))
V<-X-t(t(rep(1,n)))%*%as.vector(msh)
Ben<-diag(1,kp)-Un%*%solve(t(Un)%*%Un)%*%t(Un)
LSR<-svd(V%*%Ben)
score<-LSR$u%*%diag(LSR$d) # nonuniform scores
PWscore<-score[,1:(kp-6)]
NonUnif<-LSR$v # nonuniform coefficients

# uniform$scores  # scores for two uniform components
# PWscore # partial warp scores

allwarps<-cbind(uniform$scores,PWscore)

pcaPW<-prcomp(allwarps) # WarpScores is calculated immediately above
pcscoresPW<-pcaPW$x
eigenvaluesRWA<-pcaPW$sdev^2
varexplainedRWA<-round((eigenvaluesRWA/sum(eigenvaluesRWA))*100,3)
vartableRWA<-as.matrix(varexplainedRWA); rownames(vartableRWA)<-colnames(pcscoresPW); colnames(vartableRWA)<-c("Percent Variance Explained")
vartableRWA

barplot(varexplainedRWA,ylab="Percent Variance Explained"); title(sub="PC Rank",mgp=c(0,0,0))
plot(pcscoresPW[,1],pcscoresPW[,2], pch=19, xlab="Relative Warp 1", ylab="Relative Warp 2")
plot(pcscoresPW[,1],pcscoresPW[,3], pch=19, xlab="Relative Warp 1", ylab="Relative Warp 3")






# PCA of all warp scores (above) gives identical results to PCA of Procrustes coordinates in tangent space (orthogonal projection), except that the four zero PCs are also reported in the latter


projALLAN<-orp(allAN_GPA$rotated) # project into tangent space, orthogonal projection
projANtrix<-as.matrix(projALLAN)
pca<-prcomp(projANtrix) # this uses svd rather than eigenanalysis, and uses n-1 as the divisor when calculating variances
pcscores<-pca$x # the last four PCs explain zero variance and are redundant
eigenvalues<-pca$sdev^2
varexplained<-round((eigenvalues/sum(eigenvalues))*100,3)
vartable<-as.matrix(varexplained); rownames(vartable)<-colnames(pcscores); colnames(vartable)<-c("Percent Variance Explained")
vartable

barplot(varexplained,ylab="Percent Variance Explained"); title(sub="PC Rank",mgp=c(0,0,0))
plot(pcscores[,1],pcscores[,2], pch=19, xlab="PC1", ylab="PC2")
plot(pcscores[,1], pcscores[,3], pch=19, xlab="PC1", ylab="PC3")







# For visualization of particular PCs, Claude 2008 Functions
k<-dim(allANgroup)[1]
p<-dim(allANGroup)[2]
mesh<-as.vector(mshape(projALLAN))
max1<-matrix(mesh+max(pca$x[,1])*pca$rotation[,1],k,p)
min1<-matrix(mesh+min(pca$x[,1])*pca$rotation[,1],k,p)
max2<-matrix(mesh+max(pca$x[,2])*pca$rotation[,2],k,p)
min2<-matrix(mesh+min(pca$x[,2])*pca$rotation[,2],k,p)
msh<-mshape(projALLAN)

tps(msh,min1,20)
points(min1,pch=21,bg="black")
title("PC1: left extreme")

tps(msh,max1,20)
points(max1,pch=22,bg="black")
title("PC1: right extreme")


ChREF <- readland.tps("../DATA/MapleImages/ChicANref.tps")
ChTAR <- readland.tps("../DATA/MapleImages/ChicANtarg.tps")
plotRefToTarget(msh,min1, method = "TPS")

tps(ChREF,ChTAR,20)
points(ChTAR,col="grey50")



covs <- cov()


############################