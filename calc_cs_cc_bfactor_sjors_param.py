#!/usr/bin/env python

import math,os,sys
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(energy,delta_E,Cc,output) =  parse_command_line()
	energy=1000.0*energy
	PI=3.14159265359
	step=0.001
	
	X=4.04e-14
	Cc=Cc*1.0e-3
	CS=Cc
	print wavelength(energy/1000.0)

	out=open(output,"w")
	for i in range(0,1000):
		k=float(step*i)
		delta=Cc*math.sqrt(X+SQR(delta_E/energy))*1.0e10
		if i==0:
			print "Const_cc="+str(wavelength(energy/1000.0)*delta)
		tmpx=PI*wavelength(energy/1000.0)*delta
		Ecc=math.exp(-0.5*SQR(tmpx)*SQR(SQR(k)))
		tmp=str(k)+"\t"+str(Ecc)+"\n"
		out.write(tmp)
	out.close()
	
	
	
def parse_command_line():
	usage="%prog <input energy in kV> <input deltaE> <Cc in mm> <output txt>"
	parser = OptionParser(usage=usage, version="%1")
	
	if len(sys.argv)<5: 
		print "<input energy in kV> <input deltaE> <Cc in mm> <output txt>"
		sys.exit(-1)
	
	(options, args)=parser.parse_args()
	energy=float(args[0])
	delta_E=float(args[1])
	Cc=float(args[2])
	output=args[3]
	return (energy,delta_E,Cc,output)
def SQR(x):
	y=float(x)
	return(y*y)
def wavelength(v):
	emass=510.99906
	hc=12.3984244
	w = hc/math.sqrt( v * ( 2.0*emass + v ) )
	return(w)
if __name__== "__main__":
	main()


			
