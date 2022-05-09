# -*- coding: utf-8 -*-

import sys
import numpy as np
import numpy.matlib
import math as ma



def userpara(Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,single_group_ele):
    
    Bandwidth=1 #MHz
    Noise_o=-104 #dBm  -174dBm/Hz
    Noise=pow(10, (Noise_o/10))/pow(10,3)
    RIS_group = int(Tx_antRIS/single_group_ele)
    
    #Area: Square
    lengA=200  # m
    
    #power limits
    
    P_max_o=50 #dBm
    P_max=pow(10, (P_max_o/10))/pow(10,3)
    P_k_o=10  #dBm
    P_k=pow(10, (P_k_o/10))/pow(10,3)
    P_R_o=10  #dBm
    P_R=pow(10, (P_R_o/10))/pow(10,3)
    P_A_o=10  #dBm
    P_A=pow(10, (P_A_o/10))/pow(10,3)
    P_B_o=3   #dBW;
    P_B=pow(10, (P_B_o/10))
    
    mu=5/4
    
    #Users
    BS_loc=np.array([lengA/2,lengA/2])
    
    RISloc=np.zeros([RIS_Lnum,2])
    
    for i in range(0,RIS_Lnum):   # i = l_RIS, i從0開始
        RISloc[i] = np.array([ma.cos(ma.pi*2*(i+1)/RIS_Lnum),ma.sin(ma.pi*2*(i+1)/RIS_Lnum)])*lengA/3+[lengA/2,lengA/2]
    
    User_loc = np.zeros([Num_User,2])
    
    User_loc[0] =np.random.rand(1,2)*(1/2)*(lengA*(i+1))
    
    for i in range(0,Num_User):   #i = userloc
            User_loc[i] = User_loc[0]
    
    # for i in range(0,Num_User):   #i = userloc
    #     User_loc[i] =np.random.rand(1,2)*((i+1)/2)*(lengA*(i+1))
        
    dismin=10 #m
    
    alpha=0.2 #-3.53
    beta=3.76
    degrade=1e6
    
    Distance_UserBS = np.maximum(np.sqrt(np.sum(np.square(np.matlib.repmat(BS_loc, Num_User, 1) - User_loc), axis=1)), dismin).reshape(2,1)
    PathLoss_UserBS1 = np.sqrt(np.divide(pow(10,alpha), np.power(Distance_UserBS, beta))).reshape(2,1)
    PathLoss_UserBS = np.matlib.repmat(PathLoss_UserBS1, 1, Tx_antBS)
    PathLoss_UserBS = np.multiply(np.multiply(PathLoss_UserBS,1/ma.sqrt(2)), np.random.rand(Num_User,Tx_antBS)+1j*np.random.rand(Num_User,Tx_antBS))/degrade
    
    
    
    
    Distance_UserRIS = np.zeros([Num_User,RIS_Lnum])
    
    for i in range(0,RIS_Lnum):
        Distance_UserRIS[:,i] = np.maximum(np.sqrt(np.sum(np.square(np.matlib.repmat(RISloc[i,:], Num_User, 1) - User_loc), axis=1)), dismin)
    
    PathLoss_UserRIS1 = np.sqrt(np.divide(pow(10,alpha), np.power(Distance_UserRIS, beta))).reshape(2,1)
    #PathLoss_UserRIS=np.zeros([Num_User,RIS_Lnum,Tx_antRIS], dtype = 'complex_')
    PathLoss_UserRIS=np.zeros([Tx_antRIS,Num_User,RIS_Lnum], dtype = 'complex_')
    
    for i in range(0,Num_User):       # i = k_user
        for j in range(0,RIS_Lnum):   # j = l_RIS
            #PathLoss_UserRIS[i,j,:] = (np.multiply(np.multiply(np.matlib.repmat(PathLoss_UserRIS1[i,j], Tx_antRIS, 1), 1/ma.sqrt(2)), np.random.randn(Tx_antRIS, 1)+1j*np.random.randn(Tx_antRIS, 1))).reshape(Tx_antRIS,)
            PathLoss_UserRIS[:,i,j] = (np.multiply(np.multiply(np.matlib.repmat(PathLoss_UserRIS1[i,j], Tx_antRIS, 1), 1/ma.sqrt(2)), np.random.randn(Tx_antRIS, 1)+1j*np.random.randn(Tx_antRIS, 1))).reshape(Tx_antRIS,)
            
    PathLoss_UserRIS_group_1 = np.zeros([RIS_group, single_group_ele], dtype = 'complex_')
    PathLoss_UserRIS_group_2 = np.zeros([RIS_group, single_group_ele], dtype = 'complex_')
    
    #flag_userris = 0
    for i in range(0, RIS_group):
        PathLoss_UserRIS_group_1[i,:] = PathLoss_UserRIS[i*single_group_ele:(i+1)*single_group_ele, 0, 0]
        PathLoss_UserRIS_group_2[i,:] = PathLoss_UserRIS[i*single_group_ele:(i+1)*single_group_ele, 1, 0]
    
    
    
    
    Distance_RISBS = np.zeros([RIS_Lnum,1])
    
    for i in range(0,RIS_Lnum):
        Distance_RISBS[i,0] = np.maximum(np.sqrt(np.sum(np.square(RISloc[i,:] - BS_loc))), dismin)
        
    PathLoss_RISBS1 = np.sqrt(np.divide(pow(10,alpha), np.power(Distance_RISBS, beta)))#.reshape(2,1)
    #PathLoss_RISBS=np.zeros([RIS_Lnum,Tx_antRIS,Tx_antBS], dtype = 'complex_')
    PathLoss_RISBS=np.zeros([Tx_antBS,RIS_Lnum,Tx_antRIS], dtype = 'complex_')
    
    
    for i in range(0,RIS_Lnum):
        #PathLoss_RISBS[i,:,:] = (np.multiply(np.multiply(np.matlib.repmat(PathLoss_RISBS1[i,0], Tx_antRIS, Tx_antBS), 1/ma.sqrt(2)), np.random.randn(Tx_antRIS, Tx_antBS)+1j*np.random.randn(Tx_antRIS, Tx_antBS)))#.reshape(Tx_antRIS,)
        PathLoss_RISBS[:,i,:] = (np.multiply(np.multiply(np.matlib.repmat(PathLoss_RISBS1[i,0], Tx_antBS, Tx_antRIS), 1/ma.sqrt(2)), np.random.randn(Tx_antBS, Tx_antRIS)+1j*np.random.randn(Tx_antBS, Tx_antRIS)))
    
    risbs = np.zeros([Tx_antRIS, Tx_antBS], dtype = 'complex_')
    
    PathLoss_RISBS_group = np.zeros([RIS_group, single_group_ele, Tx_antBS], dtype = 'complex_')
    #PathLoss_RISBS_group_2 = np.zeros([RIS_group, single_group_ele, Tx_antBS], dtype = 'complex_')
    
    for j in range(0, Tx_antBS):
        risbs[:,j] = PathLoss_RISBS[j,:,:]
    
    for i in range(0, RIS_group):
        PathLoss_RISBS_group[i,:,:] = risbs[i*single_group_ele:(i+1)*single_group_ele,:]
    
    
    
    
    return Bandwidth, Noise, P_max, P_k, P_R, P_A, P_B, mu, PathLoss_UserBS, PathLoss_UserRIS, PathLoss_RISBS, PathLoss_UserRIS_group_1,\
        PathLoss_UserRIS_group_2, PathLoss_RISBS_group