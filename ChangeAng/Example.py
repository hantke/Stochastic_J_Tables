import numpy as np
from RandomAngle import RandomAngle as RAng 
from numpy.random import random as rand
from numpy.random import normal

##############################################
#How to run:
#1) Go to the Package folder and run: python setup.py build
#2) Copy Package/build/lib.linux-x86_64-2.7/RandomAngle.so into this folder (replace the old file)
#3) Run this code as: python Example.py
#This will reproduce the redshift  evolution of the direction of the angular momentum vector 
#of a dark matter halo in function of it change of mass, change  of redshift an halo mass
#############################################

def cos_ang(Vec1,Vec2):
	Mod1 = np.linalg.norm(np.array(Vec1))
	Mod2 = np.linalg.norm(np.array(Vec2))
	return (Vec1[0]*Vec2[0]+Vec1[1]*Vec2[1]+Vec1[2]*Vec2[2])/(Mod1*Mod2+1e-20)
def random_angle(Vec,Alpha):#Vec, the 3D vector, and Alpha, the angle in radians
	return RAng(Vec[0],Vec[1],Vec[2],Alpha,2.*rand())
def random_angle_Dir(Vec,Vec_Pass,Alpha,Dir,Limit = 3000):
	Vec = np.array(Vec)
	Vec_Pass = np.array(Vec_Pass)
	Phi = np.cos(Dir*2.*np.pi/180) #The prohibition angle is equal to the double of the median of the direction angle
	Diff1 = Vec-Vec_Pass
	for i in range(Limit):
		Vec_New = random_angle(Vec,Alpha/180.*np.pi)
		Diff2 = np.array(Vec_New)-np.array(Vec)
		dot = cos_ang(Diff1,Diff2)
		#if i > 0.99*Limit:	print Dir,Phi,dot,Vec_Pass,Vec,Vec_New
		if dot > Phi:	return Vec_New
	print 'WARNING ON random_angle_Dir(). Number of iterations above Limit = ',Limit
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
	Table_3D_Ang = np.load('Table_3D_Ang.npy')
	Table_3D_Dir = np.load('Table_3D_Dir.npy')
	
	# Example, lets take a halo with and specific angular momentum of 1,0,-1,
	# a halo mass of 10^12 h^-1 Msolar, z=1 and evolve it till z = 0, with an
	# increase of mass of ~1% per snapshot
	
	####################
	z = 1 #Initial Redshift
	dz = 0.05 #Change of redshift between the snapshots (can be not constant too)
	M = 1e12 #Initial halo mass
	####################
	
	j = []
	j.append(np.array([1,0,-1]))
	print '#z Mass/h^-1Msol jx jy jz alpha'
	while z > 0:
		dm = abs(normal(0.01,0.002)) #Change of mass between the snapshots
		dlogz = np.log10(dz)
		dlogM = np.log10(dm)
		logM = np.log10(M)
		ind_x,ind_y,ind_z = Get_Index(logM,11.5,14,20),Get_Index(dlogM,-3.,0,20),Get_Index(dlogz,-2.5,-0.5,20)
		Alpha = Table_3D_Ang[ind_x][ind_y][ind_z]
		print z,logM,j[-1][0],j[-1][1],j[-1][2],Alpha
		if z == 1:
			j.append(random_angle(j[-1],Alpha/180.*np.pi))
		else:
			Dir = Table_3D_Dir[ind_x][ind_y][ind_z]
			j.append(random_angle_Dir(j[-1],j[-2],Alpha,Dir))
		
		M *= (1+dm)
		z -= dz
		
	print z,logM,j[-1][0],j[-1][1],j[-1][2],Alpha
if __name__ == "__main__":
	main()