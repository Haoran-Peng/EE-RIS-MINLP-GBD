# -*- coding: utf-8 -*-
import numpy as np
import numpy.matlib
import math as ma
#from cvxopt import matrix, solvers
import cvxpy as cvx


def RIS_SCA_thetaopt(powerini,thetamarini,RIS_ele_control,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS
                     ,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,):
    
    
    Qnum = int(np.sum(xonoff)*Tx_antRIS)
    # RIS_group = int(Tx_antRIS/single_group_ele)
    
    glvec_1 = np.zeros([Tx_antBS,1], dtype = 'complex_')
    glvec_2 = np.zeros([Tx_antBS,1], dtype = 'complex_')
    glvec_1[:,0] = PathLoss_UserBS[0,:]
    glvec_2[:,0] = PathLoss_UserBS[1,:]
    
    userris_1 = np.zeros([Tx_antRIS,1], dtype = 'complex_')
    
    risbs = np.zeros([Tx_antRIS, Tx_antBS], dtype = 'complex_')
    
    Ulmar_1 = np.zeros([Qnum,Tx_antBS], dtype = 'complex_')
    Ulmar_2 = np.zeros([Qnum,Tx_antBS], dtype = 'complex_')
    
    # Ulmar_11 = np.zeros([RIS_group, single_group_ele, Tx_antBS], dtype = 'complex_')
    # Ulmar_22 = np.zeros([RIS_group, single_group_ele, Tx_antBS], dtype = 'complex_')
    
    U_diag = np.eye(Tx_antRIS)
    
    theta_initial = np.ones([Tx_antRIS, 1], dtype = 'complex_')
    thetavec_1 = np.zeros([Tx_antRIS, 1], dtype = 'complex_')
    
    
    for j in range(0, Tx_antBS):
        risbs[:,j] = PathLoss_RISBS[j,:,:]
    
    flag = 0

    for i in range(0,RIS_Lnum):
        userris_1[:,0] = PathLoss_UserRIS[:,0,:].reshape(32,)
        # userris_1[:,0] = PathLoss_UserRIS[:,i,0]
        # risbs[:,:] = PathLoss_RISBS[i,:,:]
        # #risbs[:,:] = PathLoss_RISBS[:,i,:]
        
        # print(U_diag*(PathLoss_UserRIS_group_1[i,:].conj()))
        # print((U_diag*(PathLoss_UserRIS_group_1[i,:].conj())).shape)
        # print(PathLoss_RISBS_group[i,:,:].shape)
        
        Ulmar_1[flag*Tx_antRIS : (flag+1)*Tx_antRIS,:] = np.dot(U_diag*(userris_1.conj().T), risbs)
        flag = flag+1
    
    # flag = 0
    # for i in range(0,RIS_group):
        
    #     Ulmar_11[i,:,:] = np.dot(U_diag*(PathLoss_UserRIS_group_1[i,:].conj()), PathLoss_RISBS_group[i,:,:])
    #     flag = flag+1
    
    
    # for i in range(0,RIS_group):
    #     # userris_2[:,0] = PathLoss_UserRIS[1,i,:]
    #     # risbs[:,:] = PathLoss_RISBS[i,:,:]
        
    #     # print((U_diag*(userris_1.conj().T)).shape)
        
    
    #     #Ulmar_2[flag*Tx_antRIS : (flag+1)*Tx_antRIS,:] = np.dot(U_diag*(userris_2.conj().T), risbs)
    #     flag = flag+1
    
    # flag = 0    
    # for i in range(0,RIS_group):
        
    #     Ulmar_22[i,:,:] = np.dot(U_diag*(PathLoss_UserRIS_group_2[i,:].conj()), PathLoss_RISBS_group[i,:,:])
    #     flag = flag+1
        

    flag = 0    
    for i in range(0, RIS_Lnum):
            theta_initial[:, i] = thetamarini[flag*Tx_antRIS:(flag+1)*Tx_antRIS,0]
            flag = flag+1
            
    itermax = 3e2
    itermax = int(itermax)
    epsilon = 1e-10
    itermax_1 = int(itermax+1)
    # gainvalue_1 = np.zeros([RIS_group,itermax_1, 1], dtype = 'complex_')
    gainvalue_1 = np.zeros([itermax_1, 1], dtype = 'complex_')
    
    gainvalue_1[0,0] = np.linalg.norm(glvec_1+ np.dot(Ulmar_1.conj().T, theta_initial))
    
    for index_iter in range(0, itermax):
        temp_vecupd = glvec_1+ np.dot(Ulmar_1.conj().T, theta_initial)
        # print(glvec_2.shape)
        # print(Ulmar_11[i,:,:].shape)
        # print(theta_initial[:,i].shape)
        # print(temp_vecupd.shape)
        thetavecupd = np.dot(Ulmar_1, temp_vecupd)
        thetavecupd = thetauni(thetavecupd)
        gainvalue_1[index_iter+1,0] = np.linalg.norm(glvec_1 + np.dot(Ulmar_1.conj().T, thetavecupd))
        if (abs(gainvalue_1[index_iter+1,0] - gainvalue_1[index_iter,0])/gainvalue_1[index_iter,0] < epsilon and index_iter>100):
            thetavec_1 = thetavecupd.conj()
            continue
        theta_initial = thetavecupd
    
    thetavec_1 = thetavecupd.conj()
    
    # for i in range(0, RIS_group):
    #     if (RIS_ele_control[i] == 0):
    #         gainvalue_1[i,0,0] = np.linalg.norm(glvec_2+ np.dot(Ulmar_22[i,:,:].conj().T, theta_initial[:,i].reshape(4,1)))
            
    #         for index_iter in range(0, itermax):
    #             temp_vecupd = glvec_2+ np.dot(Ulmar_22[i,:,:].conj().T, theta_initial[:,i].reshape(4,1))
    #             # print(glvec_2.shape)
    #             # print(Ulmar_11[i,:,:].shape)
    #             # print(theta_initial[:,i].shape)
    #             # print(temp_vecupd.shape)
    #             thetavecupd = np.dot(Ulmar_22[i,:,:], temp_vecupd)
    #             thetavecupd = thetauni(thetavecupd)
    #             gainvalue_1[i,index_iter+1,0] = np.linalg.norm(glvec_2+ np.dot(Ulmar_22[i,:,:].conj().T, thetavecupd))
    #             if (abs(gainvalue_1[i,index_iter+1,0] - gainvalue_1[i,index_iter,0])/gainvalue_1[i,index_iter,0] < epsilon and index_iter>100):
    #                 thetavec_1[:,i] = thetavecupd.conj().reshape(4,)
    #                 continue
    #             theta_initial[:,i] = thetavecupd.reshape(4,)
            
    #     else:
    #          gainvalue_1[i,0,0] = np.linalg.norm(glvec_1+ np.dot(Ulmar_11[i,:,:].conj().T, theta_initial[:,i].reshape(4,1)))
             
    #          for index_iter in range(0, itermax):
    #             temp_vecupd = glvec_1+ np.dot(Ulmar_11[i,:,:].conj().T, theta_initial[:,i].reshape(4,1))
    #             # print(glvec_1.shape)
    #             # print(Ulmar_11[i,:,:].shape)
    #             # print(theta_initial[:,i].shape)
    #             # print(temp_vecupd.shape)
    #             thetavecupd = np.dot(Ulmar_11[i,:,:], temp_vecupd)
    #             thetavecupd = thetauni(thetavecupd)
    #             gainvalue_1[i,index_iter+1,0] = np.linalg.norm(glvec_1+ np.dot(Ulmar_11[i,:,:].conj().T, thetavecupd))
    #             if (abs(gainvalue_1[i,index_iter+1,0] - gainvalue_1[i,index_iter,0])/gainvalue_1[i,index_iter,0] < epsilon and index_iter>100):
    #                 thetavec_1[:,i] = thetavecupd.conj().reshape(4,)
    #                 continue
    #             theta_initial[:,i] = thetavecupd.reshape(4,)
        
    #     thetavec_1[:,i] = thetavecupd.conj().reshape(4,)
        
    theta_SCAopti = np.zeros([Tx_antRIS, 1], dtype = 'complex_')
    
    
    theta_SCAopti = thetavec_1
    
    # print(type(theta_SCAopti))
    # print(theta_SCAopti.shape)
    # theta_SCAopti = theta_SCAopti
    return theta_SCAopti
            
def thetauni(thetavec_input):
    theta_uni = thetavec_input
    theta_lennum = len(thetavec_input)
    
    for i in range(0, theta_lennum):
        theta_uni[i] = thetavec_input[i]/np.linalg.norm(thetavec_input[i])
        
    return theta_uni
    
    
        