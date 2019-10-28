#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
Ce programme permet de calculer le rmsf ainsi que le B-factor.
"""

__author__ = "Emmanuel Edouard MOUTOUSSAMY && Igor Ulrich MOUTOUSSAMY"
__version__  = "1.0.0"
__copyright__ = "copyleft"
__date__ = "2014/11"

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def error_args():
	"""Cette fonction permet d'afficher une message d'erreur. Ce message permet
	à l'utilisateur de prendre connaissance des arguments à entrer.
	"""
	sys.exit("""
ERREUR! Les arguments donnés sont incorrect!
Voici les arguments à entrer:

option 				type			Description
--------------------------------------------------------------------------------
-f 				.xtc			fichier de tajectoire
-s 				.tpr			fichier de topologie
-r(optionnel) 			.pdb 			pdb(comparaison B-factor)
-p(optionnel) 			string 			chemin d'accés aux pdb   
""")

def check_args():
	"""Cette fonction permet de contrôler les arguments passer à la ligne de 
	commande.
	Retourne: une liste conteant les argument sous le format suivant : 
	[.xtc, .tpr, .pdb ,chemin d'accés au pdb]
	"""
	args = [0,0,0,0] # liste permettant de stocker les arguments
	for i in range(len(sys.argv)): # Parcourt des argument

		if sys.argv[i] == "-f": #test du fichier de trajectoire .xtc
			if sys.argv[i + 1][-4:] ==".xtc":	
				args[0] = sys.argv[i + 1]

		elif sys.argv[i] == "-s": #Test du fichier .tpr
			if sys.argv[i + 1][-4:] ==".tpr":
				args[1] = sys.argv[i + 1]

		elif sys.argv[i] == "-r": #Test du fichier .pdb
			if sys.argv[i + 1][-4:] ==".pdb":
				args[2] = sys.argv[i + 1]

		elif sys.argv[i] == "-p": #test du chemin d'accés
			args[3] = sys.argv[i + 1]

	# Affichage d'un message d'erreur si aucun fichier .xtc n'est donné
	if args[0] == 0 and args[3] == 0: 
		error_args()
	# Affichage d'un message d'erreur si aucun fichier .tpr n'est donné
	if args[1] == 0 and args[3] == 0:
		error_args()	

	return args

def mk_directory():
	"""Cette fonction permet de Créer un repertoire pour contenir les résultats.
	Retourne: le nom du repertoire.
	"""
	i = 1
	while os.path.exists("Results_%i"%i): #test de l'existence du repertoire
		i = i +1
	dir_name = "Results_%i"%i
	os.makedirs(dir_name) #Création d'un repertoire
	return dir_name

def get_pdb_trajectory(xtc,tpr,path):
	"""Cette fonction utilise une fonction de Gromacs pour obtenir les différente 
	frame d'une trajectoire.
	-Arguments : une fichier xtc, un fichier tpr et un repertoire
	-Retourne : les différente frame sour format pdb dans le repertoire
	"""
	command = 'echo 3 | trjconv -f %s -s %s -o -sep %s/out.pdb'\
	%(xtc,tpr,path)
	os.system(command) #commande gromacs pour obtenir les différentes

def extract_coord(pdb,path):
	"""Cette fonction permet d'extraire les coordonées d'un fichier pdb.
	Arguments : un fichier pdb 
	Retourne : une matrice numpy contenant les coordonées x,y,z
	"""
	with open("%s%s"%(path,pdb),"r") as input_file: #lecture du pdb
		i=0
		for line in input_file:
			if line[0:4] == "ATOM": #selection des lignes correspondant aux atomes
				if i < 1 :
					pos = line[32:55].split()
					#création d'une martice contenant les coordonées x,y,z:
					coord =np.array([pos[0],pos[1],pos[2]]) 
					i = i + 1
				else: 
					pos = line[32:55].split()
					#incrémentation de la mtrive
					coord = np.vstack((coord,pos))
		coord = coord.astype(np.float)	#convertion des données de la matrice en float
		return coord

def mean_coord(nb_frame,path):
	"""Cette fonction permet d'obtenir les coordonées moyennes.
	Arguments : le nombre de frames
	Retourne : une matrice numpy contenant les coordonées moyennes x,y,z
	pour chaque atomes.
	"""
	print"Calcul de la matrice des coordonées moyenne..."
	i = 0
	for pdb in os.listdir(path): # Parcourt des pdb
		if i == 0:
			pos = extract_coord(pdb,path) #Extraction des coordonnées du pdb
			i = i + 1
		else:
			pos2 = extract_coord(pdb,path) #Extraction des coordonées du pdb
			pos = pos + pos2 #addition des coordonnées
	pos = pos/nb_frame #divion de la matrice par le nombre de frames
	return pos

def calcul_rmsf(nb_residus,nb_frame,path,mean):
	"""Cette fonction permet de calculer le rmsf.
	Arguments : le nombre de Residus par pdb et le nombre de frame.
	Retourne : un tableau contenant les rmsf
	"""
	print"Calcul des rmsf ..."
	rmsf = [0] * nb_residus
	for pdb in os.listdir(path):
		pos = extract_coord(pdb,path)
		for i in range(nb_residus):
		#Calcul ((x-xmoy)**2) + ((y - ymoy)**2) + ((z-zmoy)**22):
 			rmsf[i] = rmsf[i] +  ((((pos[i][0] - mean[i][0])**2)+\
			((pos[i][1] - mean[i][1])**2)  + ((pos[i][2] - mean[i][2])**2)))
	for i in range(nb_residus):
		rmsf[i] = (np.sqrt(rmsf[i]/nb_frame))/10
		#la divison par 10 permet de ramener les résultats en nm. 
		rmsf[i] = round(rmsf[i],4) 
		#les résultats sont rammené à 4 chiffre significatif
	return rmsf  


def output_xvg(table,filename,title,xlab,ylab):
	"""Cette fonction permet de créer un fichier xvg contenant les rmsf.
	Arguments : un tableau contenant les rmsf.
	Retourne : un fichier xvg contenat les rmsf
	"""
	with open("%s/%s.xvg"%(dir_name,filename),"w") as output:
	#Création de l'entête:
		output.write("""#Result of the %s calculation:
@	title	"%s"
@	xaxis	label "%s"
@	yaxis	label "%s"	
@type xy\n"""%(title,title,xlab,ylab))
		j = 1
		for i in table:
			output.write("%i\t%s\n"%(j,str(i))) #Ecriture des RMSF
			j = j +1

def plot_rmsf(rmsf_table,nb_residus):
	"""Cette fonction permet tracer les rmsf en fonction des C-alpha.
	Arguments : un tableau contenant les rmsf ainsi que le nombre de résidus
	Retourne : un fichier png.
	"""
	nb_residus= range(nb_residus) 
	plt.plot(nb_residus,rmsf_table,label= 'rmsf_cal')
	plt.ylabel('RMSF(nm)')
	plt.xlabel("Residues C-alpha Atoms")
	plt.savefig("%s/RMSF2.png"%dir_name)

def calcul_bfactor(rmsf_table):
	print"Calcul des B-factor ...\n"
	"""Cette fonction permet de calculer le B-Factor.
	Arguments : un tableau conteant les rmsf
	Retourne : un tableau conteant les B-factor calculé
	"""
	bfactor_cal = [] #création d'une liste pour contenir les B-factor
	for i in rmsf_table:
		bfactor_cal.append(((i**2)*(np.pi**2)*8)/3) #calcul du B-factor
	return bfactor_cal


def extract_bfactor(pdb):
	"""Cette fonction permet d'extraire le B_factor d'un pdb.
	Arguments : un fichier pdb
	Retourne : un tableau conteant les B-factor
	"""
	with open(pdb,"r") as input_file:
		bfactor=[]
		for line in input_file:
			if line[0:4] == "ATOM":
				if line[13:15] == "CA":
					bfactor.append(float(line[61:66]))
	return bfactor

def plot_bfactor(nb_residus,bfactor_cal,bfactor,directory):
	"""Cette fonction permet de tracer le B-factor calculé et le B-factor réel
	en fonction des C-alpha.
	Arguments : un tableau conteant les B-Factor calculé et un tableau contenant
	le B-factor réel.
	Retourne : un fichier png conteant le graphique.
	"""
	plt.clf()
	plt.plot(range(nb_residus),bfactor_cal,label= 'B-Factor_cal')
	plt.plot(range(nb_residus),bfactor,label= 'B-Factor')
	plt.ylabel('B-Factor') 
	plt.xlabel("Residues C-alpha Atoms")
	plt.legend() #ajout de la legend
	plt.savefig("%s/Bfactor_cal_Compare_Bfactor.png"%dir_name)

def bfactor(nb_residus,rmsf,directory):
	"""Cette fonction permet de calculer les bfactor, d'extraire le B-factor réel
	d'un pdb et de les tracer en fonction des C-alpha.
	Arguments : un tableau conteant les B-Factor calculé et un tableau contenant
	le B-factor réel.
	Retourne : un fichier png conteant le graphique.
	"""
	bfactor_cal = calcul_bfactor(rmsf) #calcul du B_factor
	#Création d'un fichier .xvg contenant les B-factor:
	output_xvg(bfactor_cal,"bfactor_result","B-Factor","C-alpha atoms","B-factor")
	if args[2] != 0:	
		bfactor = extract_bfactor(args[2]) #obtention du B-factor
		plot_bfactor(nb_residus,bfactor_cal,bfactor,directory) #plot du B-factor et du B-Factor d'origine

def remove_pdb(directory,path):
	"""Cette fonction permet d'effacer les pdb générés los du calcul du rmsf.
	Arguments: les repertoire contenant les résultats ainsi que le chemin
	d'accés au pdb.
	"""
	if path == 0: #condition permettant de controler si des pdb ont été generés.
		choice = raw_input("Voulez vous effacer les pdb générés ? (y/n) :")
		if choice == "y":
			os.system("rm -rf %s/pdb_trajectory"%directory) #supression des pdb
		if choice != "y" and choice != "n": #limitaion du choix
			remove_pdb(directory,path) #apel recursif à la fonction
	print "\n"


def final_message(directory):
	print """Terminé!!!! Vos résultats sont dans le repertoire : %s
Merci d'avoir utilisé ce programme.
MOUTOUSSAMY Emmanuel & MOUTOUSSAMY Ulrich 2014\n\n"""%directory $!

def results(path,nb_frame):
	mean = mean_coord(nb_frame,path) #Création de la matrice de reference
	nb_residus = len(mean)	#determination du nombre de résidus
	rmsf = calcul_rmsf(nb_residus,nb_frame,path,mean) #calcul du rmsf
	#ecriture du fichier de sortie:
	output_xvg(rmsf,"rmsf_result","rms fluctuation","C-alpha atoms","nm")
	plot_rmsf(rmsf,nb_residus)	#plot du RMSF.
	bfactor(nb_residus,rmsf,dir_name) #calcul et plot du B-factor.
	remove_pdb(dir_name,args[3]) #suppression des pdb générés.
	final_message(dir_name) #permet d'afficher un message final.

if __name__ == '__main__':

	args = check_args() #contrôle des arguments
	# Création d'un répertoire pour contenir les résultats:
	dir_name = mk_directory()

	if args[3] == 0:
		#Création d'un repertoire pour contenir les pdb:
		os.makedirs("%s/pdb_trajectory"%dir_name)
		path = "%s/pdb_trajectory/"%dir_name #chemin d'accés au pdb
		get_pdb_trajectory(args[0],args[1],path) #obtention des pdb
		nb_frame = (len(os.listdir(path))-1) #determination du nombre de frames

	else:
		path = args[3] #chemin d'accés au pdb
		nb_frame = (len(os.listdir(path))-1) #determination du nombre de frames

	#calcul du rmsf,bfactor et gestion des fichier de sortie:
	results(path,nb_frame)

	

