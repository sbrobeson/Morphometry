# Run this piece by piece to obtain relevant plots of
# geometric morphometric data and analyses.
## S.B. Robeson
## last updated 20/03/2020

# Clear prior graphics device settings:
if (!is.null(dev.list())) {
dev.off()
}


# Visualize results of Partial General Procrustes Analysis for
# some individual groups (species by region):

pgpa_geom_VTAN <- gpagen(Vermont_AN)
plot(pgpa_geom_VTAN$coords[,1,],pgpa_geom_VTAN$coords[,2,],
     asp=1, col="black", cex=0.2, 
     xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate")
title("Procrustes superimposition of Vermont Acer nigrum individuals")

pgpa_geom_CHAN <- gpagen(Chicago_AN)
plot(pgpa_geom_CHAN$coords[,1,],pgpa_geom_CHAN$coords[,2,],
     asp=1, col="red", cex=0.2, 
     xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate")
title("Procrustes superimposition of Chicago Acer nigrum individuals")

pgpa_geom_CHAS <- gpagen(Chicago_AS)
plot(pgpa_geom_CHAS$coords[,1,],pgpa_geom_CHAS$coords[,2,],
     asp=1, col="blue", cex=0.2, 
     xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate")
title("Procrustes superimposition of Chicago Acer saccharum individuals")



# Visualize results of Partial Generalized Procrustes Analysis
# for the individual regions and the individual species:

# Procrustes plot, Chicago-region hard maples
{plot(allChicago_GPA$rotated[,1,], allChicago_GPA$rotated[,2,],asp=1,
      xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
  points(allChicago_GPA$rotated[,1,1:43], allChicago_GPA$rotated[,2,1:43],col="black",cex=0.5,pch=16)
  points(allChicago_GPA$rotated[,1,44:87], allChicago_GPA$rotated[,2,44:87],col="red",cex=0.5,pch=16)
  legend("topright", c("Chicago Acer nigrum","Chicago Acer saccharum"),
         pch=16, col=c("black","red"),cex=0.8)
  title("Procrustes superimposition of Chicago hard maple shapes")}

# Procrustes plot, Vermont-region hard maples
{plot(allVermont_GPA$rotated[,1,], allVermont_GPA$rotated[,2,],asp=1,
      xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
  points(allVermont_GPA$rotated[,1,1:48], allVermont_GPA$rotated[,2,1:48],col="cyan",cex=0.5,pch=16)
  points(allVermont_GPA$rotated[,1,49:92], allVermont_GPA$rotated[,2,49:92],col="magenta",cex=0.5,pch=16)
  legend("topright", c("Vermont Acer nigrum","Vermont Acer saccharum"),
         pch=16, col=c("cyan","magenta"),cex=0.8)
  title("Procrustes superimposition of Vermont hard maple shapes")}

# Procrustes plot, Virginia hard maples
{plot(allVirginia_GPA$rotated[,1,], allVirginia_GPA$rotated[,2,],asp=1,
      xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
  points(allVirginia_GPA$rotated[,1,1:43], allVirginia_GPA$rotated[,2,1:43],col="green",cex=0.5,pch=16)
  points(allVirginia_GPA$rotated[,1,44:87], allVirginia_GPA$rotated[,2,44:87],col="blue",cex=0.5,pch=16)
  legend("topright", c("Virginia Acer nigrum","Virginia Acer saccharum"),
         pch=16, col=c("green","blue"),cex=0.8)
  title("Procrustes superimposition of Virginia hard maple shapes")}

# Procrustes plot, Acer nigrum everywhere
{plot(allAN_GPA$rotated[,1,], allAN_GPA$rotated[,2,],asp=1,
      xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
  points(allAN_GPA$rotated[,1,1:44], allAN_GPA$rotated[,2,1:44],col="cyan",cex=0.5,pch=16)
  points(allAN_GPA$rotated[,1,45:88], allAN_GPA$rotated[,2,45:88],col="black",cex=0.5,pch=16)
  points(allAN_GPA$rotated[,1,89:132], allAN_GPA$rotated[,2,89:132],col="green",cex=0.5,pch=16)
  legend("topright", c("Vermont Acer nigrum", "Chicago Acer nigrum","Virginia Acer nigrum"),
         pch=16, col=c("cyan","black","green"),cex=0.8)
  title("Procrustes superimposition of Acer nigrum shapes, all regions")}

# Procrustes plot, Acer saccharum everywhere
{plot(allAS_GPA$rotated[,1,], allAS_GPA$rotated[,2,],asp=1,
      xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
  points(allAS_GPA$rotated[,1,1:44], allAS_GPA$rotated[,2,1:44],col="magenta",cex=0.5,pch=16)
  points(allAS_GPA$rotated[,1,45:88], allAS_GPA$rotated[,2,45:88],col="red",cex=0.5,pch=16)
  points(allAS_GPA$rotated[,1,89:132], allAS_GPA$rotated[,2,89:132],col="blue",cex=0.5,pch=16)
  legend("topright", c("Vermont Acer saccharum", "Chicago Acer saccharum","Virginia Acer saccharum"),
         pch=16, col=c("magenta","red","blue"),cex=0.8)
  title("Procrustes superimposition of Acer saccharum shapes, all regions")}

# Procrustes plot, both maples, all regions
{plot(allTotal_GPA$rotated[,1,], allTotal_GPA$rotated[,2,],asp=1,
      xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
  points(allTotal_GPA$rotated[,1,1:48], allTotal_GPA$rotated[,2,1:48],col="cyan",cex=0.5,pch=16)
  points(allTotal_GPA$rotated[,1,48:91], allTotal_GPA$rotated[,2,48:91],col="black",cex=0.5,pch=16)
  points(allTotal_GPA$rotated[,1,91:135], allTotal_GPA$rotated[,2,91:135],col="green",cex=0.5,pch=16)
  points(allTotal_GPA$rotated[,1,135:179], allTotal_GPA$rotated[,2,135:179],col="magenta",cex=0.5,pch=16)
  points(allTotal_GPA$rotated[,1,179:223], allTotal_GPA$rotated[,2,179:223],col="red",cex=0.5,pch=16)
  points(allTotal_GPA$rotated[,1,223:267], allTotal_GPA$rotated[,2,223:267],col="blue",cex=0.5,pch=16)
  legend("topright", c("Vermont Acer nigrum", "Chicago Acer nigrum","Virginia Acer nigrum",
                       "Vermont Acer saccharum", "Chicago Acer saccharum","Virginia Acer saccharum"),
         pch=16, col=c("cyan","black","green","magenta","red","blue"),cex=0.8)
  title("Procrustes superimposition of both maples, all regions")} 



# Plot PCA of Procrustes coordinates along with thin plate
# spline visualization of leaf morphology variation implied
# by the principal component axes:
## A graphing function with extremely useful features, but a
## premade one that obeys its own ways and not those of ggplot...


# for both species across regions:
{PCAgeomOverall <- plotTangentSpace(allTotal_GPA$rotated,
                                   groups = totalGrpList,
                                   legend = TRUE,
                                   warpgrids = TRUE)
title("PCA of Procrustes coordinates of both maples, all regions")}


# for Acer nigrum across regions:
{PCAgeomAN <- plotTangentSpace(allAN_GPA$rotated,
                              groups = ANgrpList,
                              legend = TRUE,
                              warpgrids = TRUE)
title("PCA of Procrustes coordinates of Acer nigrum, all regions")}


# for Acer saccharum across regions:
{PCAgeomAS <- plotTangentSpace(allAS_GPA$rotated,
                              groups = ASgrpList,
                              legend = TRUE,
                              warpgrids = TRUE)
title("PCA of Procrustes coordinates of Acer saccharum")}


# for both species, Virginia (an example region):
{PCAgeomVA <- plotTangentSpace(allVirginia_GPA$rotated,
                              groups = VAgrpList,
                              legend = TRUE,
                              warpgrids = TRUE)
title("PCA of Virginia hard maple Procrustes coordinates")}



############################