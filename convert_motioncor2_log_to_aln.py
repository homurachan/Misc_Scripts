#!/usr/bin/env python

import math,os,sys
try:
	from optparse import OptionParser
except:
	from optik import OptionParser

def main():
	(file_index,xsize,ysize,xpatch,ypatch,first,last,frames) =  parse_command_line()
	index=open(file_index,"r")
	index_data=index.readlines()
	for N in range(0,len(index_data)):
		star=str(index_data[N].split()[0])
		aa=open(star,"r")
		instar_line=aa.readlines()
		tmp=str(star.split(".log")[0])
		outname=tmp+".aln"
		xx=open(outname,"w")
		xx.write("setting\nstackSize: "+str(xsize)+" "+str(ysize)+" "+str(frames)+" 1\npatches: "+str(xpatch)+" "+str(ypatch)+" 0\nthrow: "+str(first)+" "+str(last)+"\nglobalShift   stackID 0\n")
		for i in range(1,len(instar_line)):
			if(instar_line[i].split()):
				Serial=int(str(instar_line[i].split()[0]))
				Shift_X=float(str(instar_line[i].split()[1]))
				if(Shift_X>=0.0):
					wSX="+"+str(Shift_X)
				else:
					wSX=str(Shift_X)
				Shift_Y=float(str(instar_line[i].split()[2]))
				if(Shift_Y>=0.0):
					wSY="+"+str(Shift_Y)
				else:
					wSY=str(Shift_Y)
				tmp1=""
				tmp1+=("\t"+str(Serial)+"\t"+str(wSX)+"\t"+str(wSY)+"\n")
				xx.write(tmp1)
		aa.close()
		xx.close()
	index.close()

def parse_command_line():
	usage="%prog <input file index> <camera X,Y size> <patch X,Y> <throw First,Last> <Total frames>"
	parser = OptionParser(usage=usage, version="%1")
	if len(sys.argv)<6: 
		print "<input file index> <camera X,Y size> <patch X,Y> <throw First,Last> <Total frames>"
		sys.exit(-1)
	(options, args)=parser.parse_args()
	file_index = str(args[0])
	tmp=str(args[1])
	xsize=int(tmp.split(',')[0])
	ysize=int(tmp.split(',')[1])
	tmp=str(args[2])
	xpatch=int(tmp.split(',')[0])
	ypatch=int(tmp.split(',')[1])
	tmp=str(args[3])
	first=int(tmp.split(',')[0])
	last=int(tmp.split(',')[1])
	frames=int(args[4])
	return (file_index,xsize,ysize,xpatch,ypatch,first,last,frames)

if __name__== "__main__":
	main()
