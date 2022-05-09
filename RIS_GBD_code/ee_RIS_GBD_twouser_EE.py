# -*- coding: utf-8 -*-
import sys
import numpy as np
import numpy.matlib
import math as ma


def RIS_EE_twousers(thetamar,power_P, power_S, xonoff, RIS_element_control,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,
                    PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS):
    
    #RIS_group = int(Tx_antRIS/single_group_ele)
    
    glvec_1 = np.zeros([Tx_antBS,1], dtype = 'complex_')
    glvec_2 = np.zeros([Tx_antBS,1], dtype = 'complex_')
    
    glvec_1[:,0] = PathLoss_UserBS[0,:]
    glvec_1 = glvec_1.conj().T
    glvec_1_off = glvec_1.conj().T
    glvec_2[:,0] = PathLoss_UserBS[1,:]
    glvec_2 = glvec_2.conj().T
    
    userris_1 = np.zeros([Tx_antRIS, 1], dtype = 'complex_')
    userris_2 = np.zeros([Tx_antRIS, 1], dtype = 'complex_')
    risbs = np.zeros([Tx_antRIS,Tx_antBS], dtype = 'complex_')
    ee_S=0
    ee2=0
    rate_S=0
    rate2=0
    
    for j in range(0, Tx_antBS):
        risbs[:,j] = PathLoss_RISBS[j,:,:]
    
    
    for i in range(0, RIS_Lnum):
        userris_1[:,0] = PathLoss_UserRIS[:,0,:].reshape(32,)
        
        thetal_1 = thetamar[i*Tx_antRIS:(i+1)*Tx_antRIS,0]
        
        glvec_1 = glvec_1 + np.dot(np.dot(userris_1.conj().T,(xonoff[i]*np.diag(thetal_1))), risbs)
        
    #-------SU 經過RIS反射-----------#
    
    barg1_S=np.linalg.norm(glvec_1)**2/(Noise + (np.linalg.norm(glvec_1)**2)*power_P)
    p_0 = P_k #+ P_R*np.sum(xonoff)*Tx_antRIS
    ee_S = Bandwidth*ma.log(1+barg1_S*power_S, 2)/(mu*power_S+p_0)
    rate_S = Bandwidth*ma.log(1+barg1_S*power_S, 2)
    
    
    barg1_P = np.linalg.norm(glvec_1)**2/(Noise + (np.linalg.norm(glvec_1)**2)*power_S)
    SINR_P = power_P*barg1_P
    
    #-------SU 沒有經過RIS反射-----------#
    
    barg1_S_off = np.linalg.norm(glvec_1_off)**2/(Noise + (np.linalg.norm(glvec_1)**2)*power_P)
    ee_S_off = Bandwidth*ma.log(1+barg1_S_off*power_S, 2)/(mu*power_S+p_0)
    rate_S_off = Bandwidth*ma.log(1+barg1_S_off*power_S, 2)
    
    barg1_P_Soff = np.linalg.norm(glvec_1)**2/(Noise + (np.linalg.norm(glvec_1_off)**2)*power_S)
    SINR_P_Soff = power_P*barg1_P_Soff
    
    # for index_group in range(0, RIS_group):
    #     if (RIS_element_control[index_group] == 1):
    #         # userris_1[:,0] = PathLoss_UserRIS[0,i,:]
    #         # risbs[:,:] = PathLoss_RISBS[i,:,:]
            
    #         #thetal_1 = thetamar_1[i*Tx_antRIS:(i+1)*Tx_antRIS,0]
            
    #         #glvec_1 = glvec_1 + np.dot(np.dot(userris_1.conj().T,(xonoff[i]*np.diag(thetal_1))), risbs)
    #         # print(PathLoss_UserRIS_group_1[index_group,:].shape)
    #         # print(thetamar[:,index_group].shape)
    #         # print(PathLoss_RISBS_group[index_group,:,:].shape)
    #         glvec_1 = glvec_1 + np.dot(np.dot(PathLoss_UserRIS_group_1[index_group,:],(xonoff[0]*np.diag(thetamar[:,index_group]))), PathLoss_RISBS_group[index_group,:,:])
            
    #         barg1_1=np.linalg.norm(glvec_1)**2/Noise;
    #         p_0 = 2*P_k + P_R*np.sum(xonoff)*Tx_antRIS +P_B
    #         ee1 += Bandwidth*ma.log(1+barg1_1*power, 2)/(mu*power+p_0)
    #         rate1 += Bandwidth*ma.log(1+barg1_1*power, 2)
    
            
    #     else:
    #         # userris_2[:,0] = PathLoss_UserRIS[1,i,:]
    #         # risbs[:,:] = PathLoss_RISBS[i,:,:]
            
    #         # thetal_2 = thetamar_2[i*Tx_antRIS:(i+1)*Tx_antRIS,0]
    #         # print(PathLoss_UserRIS_group_1[index_group,:].shape)
    #         # print(thetamar[:,index_group].shape)
    #         # print(PathLoss_RISBS_group[index_group,:,:].shape)
    #         glvec_2 = glvec_2 + np.dot(np.dot(PathLoss_UserRIS_group_2[index_group,:],(xonoff[0]*np.diag(thetamar[:,index_group]))), PathLoss_RISBS_group[index_group,:,:])
            
    #         barg1_2=np.linalg.norm(glvec_2)**2/Noise;
    #         p_0 = 2*P_k + P_R*np.sum(xonoff)*Tx_antRIS +P_B
    #         ee2 += Bandwidth*ma.log(1+barg1_2*power, 2)/(mu*power+p_0)
    #         rate2 += Bandwidth*ma.log(1+barg1_2*power, 2)
    
    # ee=ee1+ee2
    # rate=rate1+rate2
    
    return ee_S, rate_S, SINR_P, ee_S_off, rate_S_off, SINR_P_Soff
        