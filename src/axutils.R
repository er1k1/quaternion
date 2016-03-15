library(scatterplot3d)
library(rotations)
loaddata<-function(path) {
	antixdatapath="./data/jtt_data.csv"
	data=as.matrix(read.csv(antixdatapath))
	subtr=matrix(0,dim(data)[1],dim(data)[2])
	subtr[,1]<-data[1,1]
	return(data-subtr)
}

# animate the 3d data changing with time
# accept any 3d data matrix and offset to first of 3 dimensions x,y,z (2 for linear acceleration, 5 for gyro or 1 for pure 3d matrix)
# 
animateAccs<-function(data, offset, title="") {
	x=offset
	y=offset+1
	z=offset+2
	par(mfrow=c(1,1))
	amax<-max(data[,x])
	bmax<-max(data[,y])
	cmax<-max(data[,z])
	amin<-min(data[,x])
	bmin<-min(data[,y])
	cmin<-min(data[,z])
	for (i in 1:(dim(accs)[1]-3)) {
		rb<-scatterplot3d(data[c(i,i+1),x],data[c(i,i+1),y],data[c(i,i+1),z],scale.y=0.2,main=paste(round(i/100.0,digits=1),"secs","\n 3D of", title),xlab="x",ylab="z",zlab="y",xlim=c(amax,amin),ylim=c(bmax,bmin),zlim=c(cmax,cmin),color="7",type="l")
		par(new=T)
		rb<-scatterplot3d(data[c(i+1,i+2),x],data[c(i+1,i+2),y],data[c(i+1,i+2),z],scale.y=0.2,xlab="x",ylab="z",zlab="y",xlim=c(amax,amin),ylim=c(bmax,bmin),zlim=c(cmax,cmin),color="6",type="l")
		par(new=T)
		rb<-scatterplot3d(data[c(i+2,i+3),x],data[c(i+2,i+3),y],data[c(i+2,i+3),z],scale.y=0.2,xlab="x",ylab="z",zlab="y",xlim=c(amax,amin),ylim=c(bmax,bmin),zlim=c(cmax,cmin),color="5",type="l")
		par(new=F)
	}
}

# show phone axes orientation with time
# accept data in supplied raw format
#
animateQuaternions<-function(data) {
	try(Qs<-as.Q4(data[,c(11,8:10)]))
	axsS<-array(0, c(3,3))
	axsS[1,1]=axsS[2,2]=axsS[3,3]<-1
	axsE<-array(0,c(3,3))
	axsE[1,1]=axsE[2,2]=axsE[3,3]<-1
	
	for (i in seq(1,dim(Qs)[1],10)) {
		
		axsE[1,] <- rotate(Qs[i,],axsS[1,])
		axsE[2,] <- rotate(Qs[i,],axsS[2,])
		axsE[3,] <- rotate(Qs[i,],axsS[3,])
		East<-axsE[1,]
		North<-axsE[2,]
		Up<-axsE[3,]
		
		rb<-scatterplot3d(axsS[1,1],axsS[1,2],axsS[1,3],scale.y=0.2,main=paste(round(i/100.0,digits=1),"secs"),xlab="x",ylab="y",zlab="z",xlim=c(2,-2),ylim=c(2,-2),zlim=c(2,-2),color="blue",type="p")
		par(new=T)
		rb<-scatterplot3d(axsS[2,1],axsS[2,2],axsS[2,3],scale.y=0.2,main=paste(round(i/100.0,digits=1),"secs"),xlab="x",ylab="y",zlab="z",xlim=c(2,-2),ylim=c(2,-2),zlim=c(2,-2),color="green",type="p")
		par(new=T)
		rb<-scatterplot3d(axsS[3,1],axsS[3,2],axsS[3,3],scale.y=0.2,main=paste(round(i/100.0,digits=1),"secs"),xlab="x",ylab="y",zlab="z",xlim=c(2,-2),ylim=c(2,-2),zlim=c(2,-2),color="red",type="p")
		par(new=T)
		
		rb<-scatterplot3d(c(0,East[1]),c(0,East[2]),c(0,East[3]),scale.y=0.2,main=paste(round(i/100.0,digits=1),"secs"),xlab="x",ylab="y",zlab="z",xlim=c(2,-2),ylim=c(2,-2),zlim=c(2,-2),color="blue",type="l")
		par(new=T)
		rb<-scatterplot3d(c(0,North[1]),c(0,North[2]),c(0,North[3]),scale.y=0.2,main=paste(round(i/100.0,digits=1),"secs"),xlab="x",ylab="y",zlab="z",xlim=c(2,-2),ylim=c(2,-2),zlim=c(2,-2),color="green",type="l")
		par(new=T)
		rb<-scatterplot3d(c(0,Up[1]),c(0,Up[2]),c(0,Up[3]),scale.y=0.2,main=paste(round(i/100.0,digits=1),"secs"),xlab="x",ylab="y",zlab="z",xlim=c(2,-2),ylim=c(2,-2),zlim=c(2,-2),color="red",type="l")
	
		par(new=F,ask=F)
	}
	
}	

# return Hamilton product used for quaternion rotations
# accepts two quaternions in correct order
# (to do - make this accept matrices of quaternions
#
HamiltonProduct<-function(Q1, Q2) {
	#  a_1a_2 - b_1b_2 - c_1c_2 - d_1d_2
	#+ (a_1b_2 + b_1a_2 + c_1d_2 - d_1c_2)i
	#+ (a_1c_2 - b_1d_2 + c_1a_2 + d_1b_2)j
	#+ (a_1d_2 + b_1c_2 - c_1b_2 + d_1a_2)k.
	
	H <- vector("numeric",4)
	H[1] <- as.numeric(Q1[1]*Q2[1] - Q1[2]*Q2[2] - Q1[3]*Q2[3] - Q1[4]*Q2[4])
	H[2] <- as.numeric(Q1[1]*Q2[2] + Q1[2]*Q2[1] + Q1[3]*Q2[4] - Q1[4]*Q2[3])
	H[3] <- as.numeric(Q1[1]*Q2[3] - Q1[2]*Q2[4] + Q1[3]*Q2[1] + Q1[4]*Q2[2])
	H[4] <- as.numeric(Q1[1]*Q2[4] + Q1[2]*Q2[3] - Q1[3]*Q2[2] + Q1[4]*Q2[1])
	return(H[])
}

# return a vector rotated by a quaternion
# accept quaternion and vector for rotation or matrix of quaternions with same length matrix of vectors
#
rotate<-function(Q1, vec3) {
	
	if (length(Q1) > 4) {
		len<-length(is.Q4(Q1))
	} else {
		len <- 1
	}
	Qvec<-qvec(vec3)
	len2<-dim(Qvec)[1]
	if (len!= len2) {
		stop("quaternion must have same number of elements as vectors")
	}
	ax<-array(,c(len,3))
	for (i in 1:len) {
		HQ<- HamiltonProduct(as.vector(Q1[i,]), as.vector(Qvec[i,]))
		Q1conj<- -Q1[i,]
		newVec3<-HamiltonProduct(HQ, as.vector(Q1conj))
		ax[i,]<-mis.axis(as.Q4(newVec3))
	}
	return(ax)
	
}

# return relative axes of phone as timeinterval*x-axis,y-axis,z-axis*x,y,z eg 4000*3*3 for 4000 timeintervals
# accept data in supplied raw format
#
getPhoneAxes<-function(data) {
	# get Quaternions from data
	try(Qs<-as.Q4(data[,c(11,8:10)]))
	
	len <- dim(Qs)[1]
	# set up Earth axes as identity matrix
	axsS<-array(0, c(3,3))
	axsS[1,1]=axsS[2,2]=axsS[3,3]<-1
	axsE<-array(, c(len,3,3))
	j=0
	for (i in 1:len) {
		
		axsE[i,1,] <- rotate(Qs[i,],axsS[1,])
		axsE[i,2,] <- rotate(Qs[i,],axsS[2,])
		axsE[i,3,] <- rotate(Qs[i,],axsS[3,])
		
	}
	return(axsE)

}


# return rotated phone centered data as Earth centered (timeinterval*x,y,z eg 4000*3 for 4000 time intervals). 
# argument Offset is the column number for the first element (2 for linear acceleration, 5 for gyro)
#
getEarthCentered<-function(data, offset) {

	# get Quaternions from data
	len <- dim(data)[1]
	try(Qs<-as.Q4(data[,c(11,8:10)]))
	
	inverseQs<- -Qs
	
	
	# rotate phone centric data
	axsS<-array(, c(len,4))
	norm<- normalise(data,offset)
	axsS<- norm[,1:3]
	factor<- norm[,4]
	axsE<-array(, c(len,3))
	
	axsE <- rotate(inverseQs,axsS)
	
	return(axsE * factor)

}

# return normalised vector plus normalisation factor
# argument Offset is the column number for the first element (2 for linear acceleration, 5 for gyro)
#
normalise<-function(data,offset) {
	pcd<-data[,offset:(offset+2)]
	norm<-array(,c(dim(pcd)[1],dim(pcd)[2]+1))
	factor<-sqrt(rowSums(pcd^2))
	norm[,1:3]<-pcd/factor
	norm[,4]<-factor
	return(norm)
}

# return quaternion equivalent of vector
# accept vector or matrix of vectors
#
qvec<-function(vecs) {
	if (is.vector(vecs)) {
		mvecs<-t(as.matrix(vecs))
	} else {
		mvecs<-as.matrix(vecs)
	}
	qvecs<-matrix(0,dim(mvecs)[1],dim(mvecs)[2]+1)
	qvecs[,2:dim(qvecs)[2]]<-mvecs
	return(as.Q4(qvecs,1))
}



# return a timeseries index for start and finish of continuous rising or falling values in data
# accept a vector of the interesting data quantity over time 
# accept the cumulated limit of interest
#
getContinuous<-function(timeseries, limit) {

	filtBy = 1
	previ = 0
	totalRotation=0
	filt<-seq(1,length(timeseries),filtBy)
	started = F
	j = 0
	maxTotal = limit
	out<-array(,c(length(timeseries)/filtBy,3))
	for (i in filt) {
		currentRotation<-sum(timeseries[(previ+1):i])/filtBy
		if (i==1) {
			prevRotation<-currentRotation
			dirn<-T
		}
		if (((currentRotation-prevRotation)>=0)==dirn) {
			if (!started) {
				starti<-i
				dirn<-(currentRotation-prevRotation)>0
				started=T
			}
			totalRotation=currentRotation-timeseries[starti]
		} else {
			if (abs(totalRotation)>maxTotal) {
				maxTotal<-abs(totalRotation)
			}
			if (abs(totalRotation)>limit) {
				j=j+1
				out[j,]<-c(starti,i-filtBy,totalRotation)
			}
			totalRotation=0
			started=F
			dirn<-(currentRotation-prevRotation)>0
		}


		previ<- i

		prevRotation<-currentRotation
		
	}
	#out<-cbind(out,maxTotal)
	return(out[which(out[,2]!="NA"),])
}


# return integration of acceleration values corrected by mean value, as a velocity vector
# accept acceleration vector
#
getVelocity<-function(accel) {
	
	accel<-accel-mean(accel)
	vel<-array(,dim(accel))
	vel[,1]<-cumsum(accel[,1])
	vel[,2]<-cumsum(accel[,2])
	vel[,3]<-cumsum(accel[,3])
	return(vel*gettimestep())

}

# return rotation angle caused by pitch from quaternion representing rotation in x and y
# apply the rotation to the normalised x and y velocities to find vertical component
# accept data in supplied raw format
#
getFwdRotatn<-function(data) {
	accs<-getEarthCentered(antixdata,2)
	rots<-getEarthCentered(antixdata,5)
	xyvel<-getVelocity(accs)
	# normalise velocity before rotation
	normdXyvel<-normalise(xyvel,1)[,1:2]
	Qz<-matrix(0,dim(rots)[1],4)
	# create a quaternion with the x and y rotations
	Qz<-as.Q4(cbind(rots[,1:2]*gettimestep(),0))
	# apply the quaternion to the horizontal velocity vector ie in direction of forward movement
	rotatedVel<-rotate(Qz,cbind(normdXyvel,0))
	pitchRotatn <-atan(rotatedVel[,3]/sqrt((rotatedVel[,1])^2+(rotatedVel[,2])^2))

	return(pitchRotatn)

}

# return data timestep - could read this from the file

gettimestep<-function() {
	return(1/100)
}