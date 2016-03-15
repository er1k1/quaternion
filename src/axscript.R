source("axutils.R")
source("axtest.R")

# setup
if (!exists("antixdata")) {
	antixdata<-loaddata()
}

# setup a filter to sequence data by 10 when necessary
filt<-seq(1,dim(antixdata)[1],10)
# use filter for test sample
testQuaternionRotationFollowedByEarthCenteredGivesExpectedResults(antixdata[filt,])
timestep<-gettimestep()

# get Earth centred data for analysis and plots using antixutils.R functions
accs<-getEarthCentered(antixdata,2)
rots<-getEarthCentered(antixdata,5)
hzcumrots<-cumsum(as.vector(rots[,3])*timestep)
vtcumrots<-cumsum(as.vector(rots[,2])*timestep)
vels<-getVelocity(accs)
heights<-cumsum(vels[,3]-mean(vels[,3]))*timestep
fwdVel<-sqrt(vels[1]^2+vels[,2]^2)
bearing<- atan(vels[,1]/vels[,2])
i<-which(vels[,2]<0)
bearing[i]<-bearing[i]+pi
i2<-which(bearing<0)
bearing[i2]<-bearing[i2]+2*pi
pitch<-getFwdRotatn(antixdata)
pitchcumrots<-cumsum(as.vector(pitch))


	
#	first an animation of phone axes rotating within Earths axes
animateQuaternions(antixdata)
# 	and an animation of Earth centred rotations in 3d
animateAccs(cumrots,1,"Earth centred rotations")

# plots using above data
	#1. 
	par(mfrow=c(1,1),cex.main=1)
# cumulative horizontal rotation
par(new=F)
plot(hzcumrots,type="l",col="red",axes=F,xlab=NA,ylab=NA,)
box()

par(new=T)
# cumulative vertical rotation (pitch axis)
plot(pitchcumrots,type="l",col="blue",axes=F,xlab=NA,ylab=NA)
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
mtext(side = 1, "Time /100s", line = 2)
mtext(side = 2, paste("rads"), line = 2)
legend("bottomright",c("cum horiz rotn","cum vertical rotn"),col=c("red","blue"),lty = c(1, 1, 1),pch=c(1,1,NA),cex=0.65,bty="n",xjust=1)
#####
par(new=F,ask=T)
	#2.
	par(mfrow=c(1,1),cex.main=1)

# forward velocity (in bike direction)
plot(fwdVel[filt],type="l",col="green",main="forward motion",axes=F,xlab=NA,ylab=NA)
box()
par(new=T)
	
# bearing
plot(bearing[filt],type="l",col="blue",axes=F,xlab=NA,ylab=NA)
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
mtext(side = 1, "Time /100s", line = 2)
mtext(side = 2, paste("rads"), line = 2)

legend("bottomright",c("forward velocity","bearing"),col=c("green","blue"),lty = c(1, 1),pch=c(1,1),cex=0.65,bty="n",xjust=1)

# vertical
######
	#3.
	par(new=F)
par(mfrow=c(1,1),cex.main=1)

# vertical distance, velocity, acceleration
plot(heights,type="l",col="red",axes=F,xlab=NA,ylab=NA)
box()
par(new=T)

plot(vels[filt,3],type="l",col="blue",axes=F,xlab=NA,ylab=NA)
par(new=T)

#plot(accs[filt,3],type="l",col="cyan",main="z axis",axes=F,xlab=NA,ylab=NA)
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
mtext(side = 1, "Time /100s", line = 2)
mtext(side = 2, paste("heights"), line = 2)
legend("bottomright",c("rel vert velocity","rel height"),col=c("blue","red"),lty = c(1, 1, 1),pch=c(1,1,NA),cex=0.65,bty="n",xjust=1)

par(mfrow=c(2,1),cex.main=1)

	par(new=F,cex=.8)
# plot close up of rotation sequence
plot(accs[2634:2750,3],type="l", col="black",axes=F,xlab=NA,ylab=NA,main="rotation sequence at 26.34-27.5 secs")
box()
par(new=T)
plot(heights[2634:2750],type="l", col="blue",axes=F,xlab=NA,ylab=NA)
axis(side = 1, tck = -.015, labels = NA, cex=0.6)
axis(side = 2, tck = -.015, labels = NA, cex=0.6)
axis(side = 1, lwd = 0, line = -.4, cex=0.6)
axis(side = 2, lwd = 0, line = -.4, las = 1, cex=0.6)
mtext(side = 1, "Time /100s", line = 1.5, cex=0.8)
mtext(side = 2, paste("height m"), line = 2.5, cex=1)
par(new=T)
plot(vels[2634:2750],type="l", col="green",axes=F,xlab=NA,ylab=NA)
legend("bottomleft",c("rel vert accn","rel height", "rel vert velocity"),col=c("black","blue","green"),lty = c(1, 1, 1),pch=c(1,1,NA),cex=0.75,bty="n",xjust=1)


# split plot to simplify
	par(new=F,cex=.9)

plot(bearing[2634:2750],type="l", col="black",axes=F,xlab=NA,ylab=NA)
box()
par(new=T)
plot(vtcumrots[2634:2750],type="l", col="red",axes=F,xlab=NA,ylab=NA)
axis(side = 1, tck = -.015, labels = NA, cex=0.6)
axis(side = 2, tck = -.015, labels = NA, cex=0.6)
axis(side = 1, lwd = 0, line = -.4, cex=0.6)
axis(side = 2, lwd = 0, line = -.4, las = 1, cex=0.6)
mtext(side = 1, "Time /100s", line = 1.5, cex=0.8)
mtext(side = 2, paste("radians"), line = 2.5, cex=1)
par(new=T)

plot(hzcumrots[2634:2750],type="l", col="cyan",axes=F,xlab=NA,ylab=NA)
par(new=T)
plot(fwdVel[2634:2750],type="l", col="magenta",axes=F,xlab=NA,ylab=NA)
legend("bottomleft",c("bearing","vert rotn", "horiz rotn", "fwd velocity" ),col=c("black","red","cyan","magenta"),lty = c(1, 1, 1),pch=c(1,1,NA),cex=0.75,bty="n",xjust=1)
