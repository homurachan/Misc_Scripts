#!/usr/bin/env python

import math,os,sys
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(lst,boxsize,angle_step,output) =  parse_command_line()
	aa=open(lst,"r")
	instar_line=aa.readlines()
	relion30=1
	relion30=judge_relion30_or_relion31(inline=instar_line)
	print "Is relion3.0? = "+str(relion30)
	mline=-1
	if(relion30):
		mline=judge_mline0(inline=instar_line)
	else:
		#regard as relion3.1
		MLINE=judge_mline0(inline=instar_line)
		mline=judge_mline1(inline=instar_line,start=MLINE)
	if(mline<0):
		print "starfile error or this script cannot handle it. EXIT."
		quit()
	print "mline = "+str(mline)
	offset=boxsize/2.0
	R=offset
	out=open(output,"w")
	for i in range(0,mline):
		if (instar_line[i].split()):
			if (str(instar_line[i].split()[0])=="_rlnAngleRot"):
				ROT_index=int(instar_line[i].split('#')[1])-1
			if (str(instar_line[i].split()[0])=="_rlnAngleTilt"):
				TILT_index=int(instar_line[i].split('#')[1])-1
			if (str(instar_line[i].split()[0])=="_rlnAnglePsi"):
				PSI_index=int(instar_line[i].split('#')[1])-1
			
	PDF=[]
	C_PDF=0
	Min_DISTANCE=angle_step
	ROT_=[]
	TILT_=[]
	PSI_=[]
	N=0.0
	for l in range(mline,len(instar_line)):
		if(l-mline)%1000==0:
			print "Finish "+str(l-mline)+" , C_PDF="+str(C_PDF)
		if (instar_line[l].split()):
			rot=float(instar_line[l].split()[ROT_index])
			tilt=float(instar_line[l].split()[TILT_index])
			psi=float(instar_line[l].split()[PSI_index])
	#		print rot,tilt,psi
			N+=1.0
			if(C_PDF<1):
				PDF.append([])
				ROT_.append([])
				TILT_.append([])
				PSI_.append([])
				C_PDF+=1
				PDF[C_PDF-1]=0.0
				ROT_[C_PDF-1]=rot
				TILT_[C_PDF-1]=tilt
				PSI_[C_PDF-1]=psi
			if(C_PDF>=1):
				Distance=99999999.0
				for i in range(0,C_PDF):
		#			Distance_tmp=calculateAngularDistance(rot1=rot, tilt1=tilt,psi1=psi,rot2=float(ROT_[i]), tilt2=float(TILT_[i]),psi2=float(PSI_[i]))
					X1=Euler_angles2direction(alpha=rot, beta=tilt)
					X2=Euler_angles2direction(alpha=float(ROT_[i]), beta=float(TILT_[i]))
			#		print X1,X2
					dot_mult=dotProduct(v1=X1,v2=X2)
					if(math.fabs(dot_mult)>1.0):
						continue
					Distance_tmp=math.acos(dot_mult)/3.14159265359*180.0
					if(Distance_tmp<Distance):
						Distance=Distance_tmp
						REM=i
	#			print Distance
				if(Distance>Min_DISTANCE):
					PDF.append([])
					ROT_.append([])
					TILT_.append([])
					PSI_.append([])
					C_PDF+=1
					PDF[C_PDF-1]=1.0
					ROT_[C_PDF-1]=rot
					TILT_[C_PDF-1]=tilt
					PSI_[C_PDF-1]=psi
				else:
					PDF[REM]+=1.0
	for i in range(0,C_PDF):
		PDF[i]=PDF[i]/N
#		print PDF[i]
	PDF_MEAN=calc_mean(PDF)
	PDF_SIGMA=calc_sigma(PDF)
	PDF_MAX=calc_max(PDF)
	print PDF_MEAN,PDF_SIGMA,PDF_MAX
	Rmax_frac = 0.3
	width_frac = 0.5
	for i in range(0,C_PDF):
		colscale = (float(PDF[i]) - PDF_MEAN) / PDF_SIGMA
		if(colscale>5.0):
			colscale=5.0
		if(colscale<-1.0):
			colscale=-1.0
		colscale=colscale/6.0
		colscale+=1.0/6.0
		Rp = R + Rmax_frac * R * float(PDF[i]) / PDF_MAX
		AV=Euler_angles2direction(alpha=float(ROT_[i]), beta=float(TILT_[i]))
		width=Min_DISTANCE/2.0
		if(width>0.5):
			width=0.5
		tmp=""
		if (math.fabs((R - Rp) * float(AV[0])) > 0.01 or math.fabs((R - Rp) * float(AV[1])) > 0.01 or math.fabs((R - Rp) * float(AV[2])) > 0.01):
			tmp+=(".color "+str(colscale)+" 0 "+str(1. - colscale)+"\n")
			tmp+=(".cylinder "+str(R*float(AV[0])+(offset))+" "+str(R*float(AV[1])+(offset))+" "+ str(R*float(AV[2])+(offset))+ " "+ str(Rp * float(AV[0])+ offset)  + " "+ str(Rp * float(AV[1])+ offset) + " "+ str(Rp * float(AV[2])+ offset) + " "+ str(width)+"\n")
			
		out.write(tmp)
	out.close()
	aa.close()
	

	
def parse_command_line():
	usage="%prog <input list> <boxsize> <angle step> <output>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<5: 
		print "<input list> <boxsize> <angle step> <output>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	lst=str(args[0])
	boxsize=float(args[1])
	angle_step=float(args[2])
	output=str(args[3])
	return (lst,boxsize,angle_step,output)
def SQR(x):
	y=float(x)
	return(y*y)
def judge_relion30_or_relion31(inline):
	trys=3
	RELION30=1
	for i in range (0,trys):
		
		if(inline[i].split()):
			if(len(inline[i].split())>=2):
				if(inline[i].split()[0][0]=="#" and int(str(inline[i].split()[2]))>30000):
					RELION30=0
					break
	return RELION30
def judge_mline0(inline):
	trys=60
	intarget=-1
	for i in range (0,trys):
		if(inline[i].split()):
			if(inline[i].split()[0][0]!="_"):
				if(intarget==1):
					return i
					break
				else:
					continue
			if(inline[i].split()[0][0]=="_"):
				intarget=1
def judge_mline1(inline,start):
	trys=70
	intarget=-1
	for i in range (start,trys):
		if(inline[i].split()):
			if(inline[i].split()[0][0]!="_"):
				if(intarget==1):
					return i
					break
				else:
					continue
			if(inline[i].split()[0][0]=="_"):
				intarget=1				
def Euler_angles2direction(alpha, beta):

	alpha = DEG2RAD(alpha)
	beta = DEG2RAD(beta)
	v=[]
	for i in range(0,3):
		v.append([])
	ca = math.cos(alpha)
	cb = math.cos(beta)
	sa = math.sin(alpha)
	sb = math.sin(beta)
	sc = sb * ca
	ss = sb * sa

	v[0]= sc
	v[1] = ss
	v[2] = cb
	return v

def DEG2RAD(x):
	return(x/180.0*3.14159265359)
def Euler_angles2matrix(alpha, beta, gamma):
	alpha = DEG2RAD(alpha)
	beta  = DEG2RAD(beta)
	gamma = DEG2RAD(gamma)
	ca =  math.cos(alpha)
	cb =  math.cos(beta)
	cg =  math.cos(gamma)
	sa =  math.sin(alpha)
	sb =  math.sin(beta)
	sg =  math.sin(gamma)
	cc =  cb * ca
	cs =  cb * sa
	sc =  sb * ca
	ss =  sb * sa
	A=[]
	for i in range(0,3):
		A.append([])
		for j in range(0,3):
			A[i].append([])
	A[0][0] =  cg * cc - sg * sa
	A[0][1] =  cg * cs + sg * ca
	A[0][2] = -cg * sb
	A[1][0] = -sg * cc - cg * sa
	A[1][1] = -sg * cs + cg * ca
	A[1][2] = sg * sb
	A[2][0] =  sc
	A[2][1] =  ss
	A[2][2] = cb
	return A

def calculateAngularDistance(rot1, tilt1,psi1,rot2, tilt2,psi2):

#	direction1=Euler_angles2direction(alpha=rot1, beta=tilt1)
#	direction2=Euler_angles2direction(alpha=rot2, beta=tilt2)
	min_axes_dist = 3600.0

	E1=Euler_angles2matrix(alpha=rot1, beta=tilt1, gamma=psi1)
	E2=Euler_angles2matrix(alpha=rot2, beta=tilt2, gamma=psi2)
	v1=[]
	v2=[]
	axes_dist = 0;
	for i in range(0,3):
		v1=E1[i]
		v2=E2[i]
		axes_dist += math.acos(CLIP(a=dotProduct(v1, v2),b=-1., c=1.))*180.0/3.14159265359
	axes_dist=axes_dist/3.0
	if (axes_dist < min_axes_dist):
		min_axes_dist = axes_dist
	return min_axes_dist
	
def CLIP(a,b,c):
	if(float(a)<float(b)):
		return float(b)
	else:
		if(float(a)>float(c)):
			return float(c)
		return float(a)
def dotProduct(v1,v2):
	if(len(v1)!=len(v2)):
		return -9999999.0
	sum=0.0
	for i in range(0,len(v1)):
		sum+=float(v1[i])*float(v2[i])
	return sum
def calc_mean(PDF):
	X=len(PDF)
	SUM=0.0
	for i in range(0,X):	
		SUM+=float(PDF[i])
	SUM=SUM/X
	return SUM
def calc_sigma(PDF):
	X=len(PDF)
	SUM=0.0
	SUM_SQR=0.0
	for i in range(0,X):	
		SUM+=float(PDF[i])
		SUM_SQR+=SQR(float(PDF[i]))
	SUM_SQR=SUM_SQR/X
	SUM=SUM/X
	RETURN=math.sqrt(SUM_SQR-SQR(SUM))
	return RETURN
def calc_max(PDF):
	X=len(PDF)
	Y=-999999999.0
	for i in range(0,X):
		if(float(PDF[i])>Y):
			Y=float(PDF[i])
	return Y
if __name__== "__main__":
	main()
