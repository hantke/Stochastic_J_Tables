##########################################################################

# Run as python Read_Tables.py
# If you want to print the tables in an ASCII table in the Out folder try:
# python Read_Tables.py Print
import numpy as np
import sys
##########################################################################
def Table2ASCII_3D(A,Name,X1_Min = 11.,X1_Max = 14.5,X2_Min = -3.,X2_Max = 0.,X3_Min = -2.5,X3_Max = -0.5,X1_NBin = 20,X2_NBin = 20,X3_NBin = 20):
	W = file('Out/'+Name+'_3D.txt','w')
	# The 3D Table use 20 bins
	# The first dimention is Halo mass in units of log10(Mass/h^-1Msol)
	# The second dimention is the change of mass in logscale, calculated as log10((M1-M2)/M1)
	# The third dimention is the change of redshift in logscale, calculated as log10(z1-z2)
	
	print >>W,'#(1) log(M_min) log of the minimum halo mass of the selection in logscale (in units of Msol/h)'
	print >>W,'#(2) log(M_max) log of the maximum halo mass of the selection in logscale (in units of Msol/h)'
	print >>W,'#(3) log(dM_min) log of the minimum change of halo mass of the selection in logscale (calculated as (M1-M2)/M1)'
	print >>W,'#(4) log(dM_max) log of the maximum change of halo mass of the selection in logscale (calculated as (M1-M2)/M1)'
	print >>W,'#(5) log(dz_min) log of the minimum change of redshift of the selection in logscale (calculated as z1-z2)'
	print >>W,'#(6) log(dz_max) log of the maximum change of redshift of the selection in logscale (calculated as z1-z2)'
	print >>W,'#(7) '+Name
	print >>W,'#log(M_min)	log(M_max)	log(dM_min)	log(dM_max)	log(dz_min)	log(dz_max)	'+Name
	
	for i in range(X1_NBin):
		X1Min = X1_Min+(i)/float(X1_NBin)*(X1_Max-X1_Min)
		X1Max = X1_Min+(i+1)/float(X1_NBin)*(X1_Max-X1_Min)
		for j in range(X2_NBin):
			X2Min = X2_Min+(j)/float(X2_NBin)*(X2_Max-X2_Min)
			X2Max = X2_Min+(j+1)/float(X2_NBin)*(X2_Max-X2_Min)
			for k in range(X3_NBin):
				X3Min = X3_Min+(k)/float(X3_NBin)*(X3_Max-X3_Min)
				X3Max = X3_Min+(k+1)/float(X3_NBin)*(X3_Max-X3_Min)
				print >>W,X1Min,X1Max,X2Min,X2Max,X3Min,X3Max,A[i][j][k]
	W.close()
	
def Table2ASCII_4D(A,Name,X1_Min = 11.,X1_Max = 14.5,X2_Min = -3.,X2_Max = 0.,X3_Min = -2.5,X3_Max = -0.5,X4_Min = -2.5,X4_Max = 2.5,X1_NBin = 12,X2_NBin = 12,X3_NBin = 12,X4_NBin = 12):
	W = file('Out/'+Name+'_4D.txt','w')
	# The 4D Table use 20 bins
	# The first dimention is Halo mass in units of log10(Mass/h^-1Msol)
	# The second dimention is the change of mass in logscale, calculated as log10((M1-M2)/M1)
	# The third dimention is the change of redshift in logscale, calculated as log10(z1-z2)
	# The fourth dimention is the modulus of the specific angular momentum of the selection in logscale
	print >>W,'#(1) log(M_min) log of the minimum halo mass of the selection in logscale (in units of Msol/h)'
	print >>W,'#(2) log(M_max) log of the maximum halo mass of the selection in logscale (in units of Msol/h)'
	print >>W,'#(3) log(dM_min) log of the minimum change of halo mass of the selection in logscale (calculated as (M1-M2)/M1)'
	print >>W,'#(4) log(dM_max) log of the maximum change of halo mass of the selection in logscale (calculated as (M1-M2)/M1)'
	print >>W,'#(5) log(dz_min) log of the minimum change of redshift of the selection in logscale (calculated as z1-z2)'
	print >>W,'#(6) log(dz_max) log of the maximum change of redshift of the selection in logscale (calculated as z1-z2)'
	print >>W,'#(7) log(j_min) log of the minimum specific angular momentum of the selection in logscale (in physical units, Mpc/h * km/s)'
	print >>W,'#(8) log(j_max) log of the maximum specific angular momentum of the selection in logscale (in physical units, Mpc/h * km/s)'
	print >>W,'#(9) '+Name
	print >>W,'#log(M_min)	log(M_max)	log(dM_min)	log(dM_max)	log(dz_min)	log(dz_max)	log(j_min)	log(j_max)	'+Name
	
	for i in range(X1_NBin):
		X1Min = X1_Min+(i)/float(X1_NBin)*(X1_Max-X1_Min)
		X1Max = X1_Min+(i+1)/float(X1_NBin)*(X1_Max-X1_Min)
		for j in range(X2_NBin):
			X2Min = X2_Min+(j)/float(X2_NBin)*(X2_Max-X2_Min)
			X2Max = X2_Min+(j+1)/float(X2_NBin)*(X2_Max-X2_Min)
			for k in range(X3_NBin):
				X3Min = X3_Min+(k)/float(X3_NBin)*(X3_Max-X3_Min)
				X3Max = X3_Min+(k+1)/float(X3_NBin)*(X3_Max-X3_Min)
				for l in range(X4_NBin):
					X4Min = X4_Min+(l)/float(X4_NBin)*(X4_Max-X4_Min)
					X4Max = X4_Min+(l+1)/float(X4_NBin)*(X4_Max-X4_Min)
					print >>W,X1Min,X1Max,X2Min,X2Max,X3Min,X3Max,X4Min,X4Max,A[i][j][k][l]
	W.close()
	
def Get_Index(X,X_Min,X_Max,NBin):
	index = int((X-X_Min)/float(X_Max-X_Min)*NBin)
	if index < 0:
		print 'WARNING Index out of range (<0)'
		return 0
	if index >= NBin:
		print 'WARNING Index out of range (>=NBin)'
		return NBin-1
	return index
	
	
def main():
	Command = None
	if len(sys.argv) > 1:	Command = sys.argv[1]
	Angle_3D = np.load('Tables/Table_3D_Ang.npy')
	Angle_4D = np.load('Tables/Table_4D_Ang.npy')
	Dir_3D = np.load('Tables/Table_3D_Dir.npy')
	Dir_4D = np.load('Tables/Table_4D_Dir.npy')
	
	#Here we read the tables and load them in a numpy 3D and 4D array
	#The tables give the median angle and change in direction
	#in function of the Halo mass, the change of mass of the halo
	#in a time lapse, and the change of redshift of the time lapse.
	
	#To get the index, run the 'Get_Index' Function, as:
	#Get_Index(Prop,X_Min,X_Max,NBin)
	#where if Prop is the Mass, then:
		#Mass has to be in log10(Mass/h**-1Msol)
		#X_Min is 11.5
		#X_Max is 14.
	#if the Prop is change of mass:
		#dM has to be calculated as log10((M1-M2)/M1) (with M1 the final mass of the halo)
		#X_Min is -3
		#X_Max is 0
	#if the Prop is change of redshift:
		#dz has to be calculated as log10(z1-z2) (with z1 the final redshift)
		#X_Min is -2.5
		#X_Max is -0.5
		
	#for the 4D table you can also use the modulus of the specific angular momentum vector:
		#j has to be calculated as log10(j/ Mpc/h * km/s)
		#X_Min is -2.5
		#X_Max is 2.5
	#If you are using the 3D table, then NBin=20. if you are using the 4D table, then NBin=12
	#Example, for a mass M, a change of mass dm, a change of redshift dz, and a specific
	#angular momentum of j (all in logscale) the median change of angle and median change of direction
	#can be calculated as:
		#Angle = Angle_3D[Get_Index(M,11.5,14,20)][Get_Index(dm,-3.,0,20)][Get_Index(dz,-2.5,-0.5,20)]
		#Dir =     Dir_3D[Get_Index(M,11.5,14,20)][Get_Index(dm,-3.,0,20)][Get_Index(dz,-2.5,-0.5,20)]
		#Or for the 4D table:
		#Angle = Angle_4D[Get_Index(M,11.5,14,12)][Get_Index(dm,-3.,0,12)][Get_Index(dz,-2.5,-0.5,12)][Get_Index(j,-2.5,2.5,12)]
		#Dir     = Dir_4D[Get_Index(M,11.5,14,12)][Get_Index(dm,-3.,0,12)][Get_Index(dz,-2.5,-0.5,12)][Get_Index(j,-2.5,2.5,12)]
	
	if Command == 'Print':
		Table2ASCII_3D(Angle_3D,'Angle')
		Table2ASCII_3D(Dir_3D,'Dir')
		Table2ASCII_4D(Angle_4D,'Angle')
		Table2ASCII_4D(Dir_4D,'Dir')
	# Here print the Tables in an ASCII table
	
	
if __name__ == "__main__":
    main()
