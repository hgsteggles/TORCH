import sys
import warnings
import numpy as np

f = open("data/cornish/cornish-uchiis.txt", 'r')
content = f.readlines()

cornish_survey = np.genfromtxt("data/cornish/cornish-uchiis.txt", delimiter=',', dtype=str)
cornish_distances = np.genfromtxt("data/cornish/cornish-distances.txt", skip_header=1, dtype=str)

name_to_distance = dict()
have_distance = set()


for i in range(len(cornish_distances[:,0])):
	have_distance.add(cornish_distances[i,0])
	name_to_distance[cornish_distances[i,0]] = cornish_distances[i,1]

ofile = open("data/cornish/cornish-uchii-distances.txt", 'w')

ofile.write("#Name,l_deg,b_deg,RA_deg,Dec_deg,dRA_asec,dDec_asec,Peak_mJybm,dPeak_mJybm,Flux_mJy,dFlux_mJy,Angscale_asec,dAngscale_asec,AngscaleDecon_asec,gaussMajor_asec,dGaussMajor_asec,gaussMinor_asec,dGaussMminor_asec,gaussPosangle_deg,dGaussPosangle_deg,RMS_mJybm,Sky_mJybm,Sigma,mType,fArtefact,fCluster,fEdge,fHiNoise,fHi5Sig,fNearBright,fSmoothWeighting,fOverlap7sig,fOverlap5sig,dist\n")

for i in range(len(cornish_survey[:,0])):
	if cornish_survey[i,0] in have_distance:
		line = cornish_survey[i,0]
		for j in range(1, len(cornish_survey[i,:])):
			line += "," + cornish_survey[i,j]
		line += "," + name_to_distance[cornish_survey[i,0]]

		ofile.write(line + '\n')

ofile.close()