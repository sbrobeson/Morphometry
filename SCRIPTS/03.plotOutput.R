


plot(allChicago_GPA$rotated[,1,], allChicago_GPA$rotated[,2,],asp=1,
     xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
points(allChicago_GPA$rotated[,1,1:43], allChicago_GPA$rotated[,2,1:43],col="black",cex=0.5,pch=16)
points(allChicago_GPA$rotated[,1,44:87], allChicago_GPA$rotated[,2,44:87],col="red",cex=0.5,pch=16)
legend("topright", c("Chicago Acer nigrum","Chicago Acer saccharum"),
       pch=16, col=c("black","red"),cex=0.8)
title("Procrustes coordinates of Chicago-region hard maples")

plot(allVermont_GPA$rotated[,1,], allVermont_GPA$rotated[,2,],asp=1,
     xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
points(allVermont_GPA$rotated[,1,1:48], allVermont_GPA$rotated[,2,1:48],col="cyan",cex=0.5,pch=16)
points(allVermont_GPA$rotated[,1,49:92], allVermont_GPA$rotated[,2,49:92],col="magenta",cex=0.5,pch=16)
legend("topright", c("Vermont Acer nigrum","Vermont Acer saccharum"),
       pch=16, col=c("cyan","magenta"),cex=0.8)
title("Procrustes coordinates of Vermont hard maples")

plot(allVirginia_GPA$rotated[,1,], allVirginia_GPA$rotated[,2,],asp=1,
     xlab="Procrustes x-coordinate", ylab="Procrustes y-coordinate",cex=0)
points(allVirginia_GPA$rotated[,1,1:43], allVirginia_GPA$rotated[,2,1:43],col="green",cex=0.5,pch=16)
points(allVirginia_GPA$rotated[,1,44:87], allVirginia_GPA$rotated[,2,44:87],col="blue",cex=0.5,pch=16)
legend("topright", c("Virginia Acer nigrum","Virginia Acer saccharum"),
       pch=16, col=c("green","blue"),cex=0.8)
title("Procrustes coordinates of Virginia hard maples")


ggplot(allVirginia_GPA$rotated[,1,], allVirginia_GPA$rotated[,2,]) +
  ggtitle("Procrustes coordinates of Virginia hard maples") +
  aes()
  labs()