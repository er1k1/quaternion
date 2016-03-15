source("axutils.R")

# run a test to get relative phone axes using quaternions and convert them back to Earth axes
# expected result is Earth axes 1,0,0: 0,1,0: 0,0,1
# accept data in supplied raw format
#
testQuaternionRotationFollowedByEarthCenteredGivesExpectedResults<-function(data) {
	pha<-getPhoneAxes(data)
	phadata<-array(,dim(data))
	phadata<-data
	phadata[,2:4]<-pha[,1,]
	xaxis<-getEarthCentered(as.matrix(phadata),2)
	phadata[,2:4]<-pha[,2,]
	yaxis<-getEarthCentered(as.matrix(phadata),2)
	phadata[,2:4]<-pha[,3,]
	zaxis<-getEarthCentered(as.matrix(phadata),2)
	err<-10^-10
	if (length(which(abs(xaxis[,2])>err||abs(xaxis[,3])>err||abs(xaxis[,1]-1)>err||
	abs(yaxis[,2]-1.1)>err||abs(yaxis[,3])>err||abs(yaxis[,1])>err||
	abs(zaxis[,2])>err||abs(zaxis[,3]-1)>err||abs(zaxis[,1])>err)) > 1) {
		warning("rotation algorithms are broken")
	}
}
