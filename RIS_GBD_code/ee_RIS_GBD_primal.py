# -*- coding: utf-8 -*-
import sys
import numpy as np
import numpy.matlib
import math as ma
import cvxpy as cvx


def RIS_GBD_primal(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoff,
                   Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS, power_P, power_S,thetamarini):
    
    
    Qnum = int(np.sum(xonoff)*Tx_antRIS)
    #RIS_Cnum = 
    
    glvec_1 = np.zeros([Tx_antBS,1], dtype = 'complex_')
    glvec_2 = np.zeros([Tx_antBS,1], dtype = 'complex_')
    glvec_1[:,0] = PathLoss_UserBS[0,:]
    glvec_2[:,0] = PathLoss_UserBS[1,:]
    
    # glvec_1[:,0] = PathLoss_UserBS[0,:]
    glvec_1_poweropt = glvec_1.conj().T
    glvec_1_poweropt_off = glvec_1.conj().T
    
    userris_1 = np.zeros([Tx_antRIS,1], dtype = 'complex_')
    userris_2 = np.zeros([Tx_antRIS,1], dtype = 'complex_')
    
    risbs = np.zeros([Tx_antRIS, Tx_antBS], dtype = 'complex_')
    
    Ulmar_1 = np.zeros([Tx_antRIS,Tx_antBS], dtype = 'complex_')
    # Ulmar_2 = np.zeros([Qnum,Tx_antBS], dtype = 'complex_')
    
    
    # Ulmar_11 = np.zeros([RIS_group, single_group_ele, Tx_antBS], dtype = 'complex_')
    # Ulmar_22 = np.zeros([RIS_group, single_group_ele, Tx_antBS], dtype = 'complex_')
    
    # U_diag = np.eye(single_group_ele)
    
    U_diag = np.eye(Tx_antRIS)
    
    theta_initial = np.ones([Tx_antRIS, 1], dtype = 'complex_')
    thetavec_1 = np.zeros([Tx_antRIS, 1], dtype = 'complex_')
    
    for j in range(0, Tx_antBS):
        risbs[:,j] = PathLoss_RISBS[j,:,:]
    

    flag = 0
    #------------------First RIS element--------------#
    
    # for i in range(0,RIS_Lnum):
    #     userris_1[:,0] = PathLoss_UserRIS[0,i,:]
    #     risbs[:,:] = PathLoss_RISBS[i,:,:]
        
    #     #print((U_diag*(userris_1.conj().T)).shape)
        
    #     if (xonoff[i]==1):
    #         Ulmar_1[flag*Tx_antRIS : (flag+1)*Tx_antRIS,:] = np.dot(U_diag*(userris_1.conj().T), risbs)
    #         flag = flag+1
    
    #------------------second RIS group--------------#
    # for i in range(0,RIS_group):
        
    #     Ulmar_11[i,:,:] = np.dot(U_diag*(PathLoss_UserRIS_group_1[i,:].conj()), PathLoss_RISBS_group[i,:,:])
    #     flag = flag+1
        
    #------------------Third PU,SU--------------#
    
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
    
    flag = 0
    # ------------------First RIS element--------------#
    
    # for i in range(0,RIS_Lnum):
    #     userris_2[:,0] = PathLoss_UserRIS[1,i,:]
    #     risbs[:,:] = PathLoss_RISBS[i,:,:]
        
    #     #print((U_diag*(userris_1.conj().T)).shape)
        
    #     if (xonoff[i]==1):
    #         Ulmar_2[flag*Tx_antRIS : (flag+1)*Tx_antRIS,:] = np.dot(U_diag*(userris_2.conj().T), risbs)
    #         flag = flag+1
    
    #------------------second RIS group--------------#
    
    # for i in range(0,RIS_group):
        
    #     Ulmar_22[i,:,:] = np.dot(U_diag*(PathLoss_UserRIS_group_2[i,:].conj()), PathLoss_RISBS_group[i,:,:])
    #     flag = flag+1
    
    
    flag = 0
    # theta_initial = np.ones([single_group_ele, RIS_group], dtype = 'complex_')
    # thetavecini_1 = np.ones([Qnum,1], dtype = 'complex_')
    # thetavecini_2 = np.ones([Qnum,1], dtype = 'complex_')
    
    #------------------First RIS element--------------#

    # for i in range(0,RIS_Lnum):
    #     if(xonoff[i] == 1):
    #         theta_initial[flag*Tx_antRIS:(flag+1)*Tx_antRIS, 0] = thetamarini[i*Tx_antRIS:(i+1)*Tx_antRIS, 0]
    #         flag = flag+1
    
    #------------------second RIS group--------------#
    
    # for i in range(0, RIS_group):
    #         theta_initial[:, i] = thetamarini[flag*single_group_ele:(flag+1)*single_group_ele,0]
    #         flag = flag+1
        
    #------------------Third PU,SU--------------#
    
    flag = 0    
    for i in range(0, RIS_Lnum):
            theta_initial[:, i] = thetamarini[flag*Tx_antRIS:(flag+1)*Tx_antRIS,0]
            flag = flag+1
    
    # RIS_element_control_inv = np.ones([RIS_group,1])
    # RIS_element_control_inv = RIS_element_control_inv - RIS_element_control
    # thetavecini_1 = thetamarini_pri1 * RIS_element_control
    # thetavecini_2 = thetamarini_pri2 * RIS_element_control_inv
    
    
    # prob_para_1 = glvec_1 + np.dot(Ulmar_1.conj().T, thetavecini_1)
    # prob_para_2 = glvec_2 + np.dot(Ulmar_2.conj().T, thetavecini_2)
    
    # prob_para_1_mid = np.array(np.linalg.norm(prob_para_1)**2).reshape(1,1)
    # prob_para_2_mid = np.array(np.linalg.norm(prob_para_2)**2).reshape(1,1)
    
    # prob_para_1_last = 2*(np.dot(np.dot(prob_para_1.conj().T, Ulmar_1.conj().T), thetavecini_1)).real
    # prob_para_2_last = 2*(np.dot(np.dot(prob_para_2.conj().T, Ulmar_2.conj().T), thetavecini_2)).real
    
    # thetaopti_real_1 = cvx.Variable(Qnum, complex = True)
    #thetaopti_img_1 = cvx.Variable(Qnum)
    
    #thetaopti_real_2 = cvx.Variable(Qnum, complex = True)
    #thetaopti_img_2 = cvx.Variable(Qnum)
    
    theta_PU_SU_RISphase = cvx.Variable(Tx_antRIS, complex = True)
    
    
    flag_obj = 0
    obj_pri_1=0
    obj_pri_2=0
    
    #-------Scond RIS group-----------#
    
    # for i in range(0, RIS_group):
    #     # print(glvec_2.shape)
    #     # print((glvec_2 + np.dot(Ulmar_22[i,:,:].conj().T, theta_initial[:,i])).shape)
    #     if (RIS_element_control[i] == 0):
    #         temp_vec_2 = (glvec_2 + (np.dot(Ulmar_22[i,:,:].conj().T, theta_initial[:,i])).reshape(8,1)).conj().T
    #         obj_pri_2 += 2*cvx.real(temp_vec_2@Ulmar_22[i,:,:].conj().T@thetaopti_real_2[i*single_group_ele:(i+1)*single_group_ele])
    #     # else:
    #     #     temp_vec = glvec_2
    #     #     obj_pri += norm(glvec_2)
    
    
    # obj_in_real = cvx.real(2*(prob_para_1.conj().T@Ulmar_1.conj().T@thetaopti_real_1) + 2*(prob_para_2.conj().T@Ulmar_2.conj().T@thetaopti_real_2))
    # obj = cvx.Maximize(obj_in_real)
    
    #-------Third PU,SU-----------#
    
    # print(glvec_2.shape)
    # print((glvec_2 + np.dot(Ulmar_22[i,:,:].conj().T, theta_initial[:,i])).shape)
    
    
    
    temp_vec_2 = (glvec_1 + (np.dot(Ulmar_1.conj().T, theta_initial))).conj().T

    
    obj_theta_phase = 2*cvx.real(temp_vec_2@Ulmar_1.conj().T@theta_PU_SU_RISphase)
    
    obj_theta = cvx.Maximize(obj_theta_phase)
    
    
    
    constraints_theta = []
    
    for i in range(0, Tx_antRIS):
        constraints_theta.append(cvx.abs(theta_PU_SU_RISphase[i]) <= 1)
    
    
    prob_theta = cvx.Problem(obj_theta, constraints_theta)
    prob_theta.solve()
    
    theta_flag = prob_theta.status
    theta_temp = theta_PU_SU_RISphase.value.conj()
    
    # else:
    #     temp_vec_2 = (glvec_1 + (np.dot(Ulmar_1.conj().T, theta_initial))).conj().T
    #     obj_SU = (cvx.norm(glvec_1)**2)*power_S
        
    #     obj_PU = 2*cvx.real(temp_vec_2@Ulmar_1.conj().T@theta_PU_SU_RISphase)*power_PU 
        
    #     obj = cvx.Maximize(obj_SU/(obj_PU + Noise))
    
        
    #--------------First RIS element---------------#
    
    # theta_num1 = 0
    # theta_num2 = 0
    # theta_close1 = 0
    # theta_close2 = 0
    
    # constraints = []
    # for i in range(0, Qnum):    
    #     if(RIS_element_control[i] == 0):
    #         constraints.append(cvx.abs(thetaopti_real_2[i]) <= 1)
    #         theta_num2 = theta_num2 + 1
    #         #constraints.append(cvx.abs(thetaopti_img_2[i]) <= 1/2)
    #         constraints.append(thetaopti_real_1[i] == 0)
    #         theta_close1 = theta_close1 + 1
    #         #constraints.append(thetaopti_img_1[i] == 0)
    #     else:
    #         constraints.append(cvx.abs(thetaopti_real_1[i]) <= 1)
    #         theta_num1 = theta_num1 + 1
    #         #constraints.append(cvx.abs(thetaopti_img_1[i]) <= 1/2)
    #         constraints.append(thetaopti_real_2[i] == 0)
    #         theta_close2 = theta_close2 + 1
    #         #constraints.append(thetaopti_img_2[i] == 0)
    
    
    
    #--------------Second RIS group---------------#
    
    
    # constraints = []
    # theta_second = 0

    # for i in range(0,RIS_group):
    #     if (RIS_element_control[i] == 0):
    #         for j in range(0, single_group_ele):
    #             constraints.append(cvx.abs(thetaopti_real_2[i*single_group_ele+j]) <= 1)
    #             theta_second +=1
            
    # for i in range(0,RIS_group):
    #     if (RIS_element_control[i] == 1):
    #         temp_vec_1 = (glvec_1 + (np.dot(Ulmar_11[i,:,:].conj().T, theta_initial[:,i])).reshape(8,1)).conj().T
    #         obj_pri_1 += 2*cvx.real(temp_vec_1@Ulmar_11[i,:,:].conj().T@thetaopti_real_2[i*single_group_ele:(i+1)*single_group_ele])
    
    # constraints.append(obj_pri_1 >= 1e-20)
    
    # prob = cvx.Problem(obj, constraints)
    # prob.solve()
    
    
    #--------------Second RIS group---------------#
    
    #-------Third PU,SU-----------#
    
    power_SU = cvx.Variable()
    power_PU = cvx.Parameter()
    power_PU = power_P
    
    for i in range(0, RIS_Lnum):
        userris_1[:,0] = PathLoss_UserRIS[:,0,:].reshape(32,)
        
        thetal_1 = theta_temp[i*Tx_antRIS:(i+1)*Tx_antRIS]
        
        glvec_1_poweropt = glvec_1_poweropt + np.dot(np.dot(userris_1.conj().T,(xonoff[i]*np.diag(thetal_1))), risbs)
        
        
    
    if (xonoff[0] == 1):
        # temp_vec_2 = (glvec_1 + (np.dot(Ulmar_1.conj().T, theta_initial))).conj().T
        # obj_SU = 2*cvx.real(temp_vec_2@Ulmar_1.conj().T@theta_PU_SU_RISphase)*power_SU
        
        # obj_PU = 2*cvx.real(temp_vec_2@Ulmar_1.conj().T@theta_PU_SU_RISphase)*power_P 
        
        # obj = cvx.Maximize(obj_SU/(obj_PU + Noise))
        
        obj_SU = power_SU
        
        obj_PU = power_PU 
        
        obj_power = cvx.Maximize(obj_SU/(obj_PU + Noise))
    
    else:
        # temp_vec_2 = (glvec_1 + (np.dot(Ulmar_1.conj().T, theta_initial))).conj().T
        obj_SU = (cvx.norm(glvec_1_poweropt_off)**2)*power_SU
        
        obj_PU = np.linalg.norm(glvec_1_poweropt)**2*power_PU 
        # obj_PU = np.linalg.norm(glvec_1_poweropt)**2
        
        obj_power = cvx.Maximize(obj_SU/(obj_PU + Noise))
    
    
    
    constraints_power = []
               
    
    constraints_power.append(power_SU >= 0)
    constraints_power.append(power_SU <= P_max)
    # print(obj_PU.is_real())
    # print(obj_SU.is_real())
    constraints_power.append(obj_PU >= 5*(obj_SU + Noise))
    
    prob_power = cvx.Problem(obj_power, constraints_power)
    prob_power.solve()
    
    #-------Third PU,SU-----------#
    
    # total_conslen = theta_num1 +theta_num2 + theta_close1 + theta_close2
    
    # thetaopt_dual_1 = np.zeros([Qnum,1], dtype = 'complex_')
    # thetaopt_dual_2 = np.zeros([Qnum,1], dtype = 'complex_')
    # theta_closedual_1 = np.zeros([Qnum,1], dtype = 'complex_')
    # theta_closedual_2 = np.zeros([Qnum,1], dtype = 'complex_')
    
    # flag_dual = 0
    # for i in range(0, int(total_conslen/2)):
    #     if(RIS_element_control[i] == 1):
    #         thetaopt_dual_1[i] = constraints[flag_dual].dual_value
    #         flag_dual = flag_dual+1
    #         theta_closedual_2[i] = constraints[flag_dual].dual_value
    #         flag_dual = flag_dual+1
    #     else:
    #         thetaopt_dual_2[i] = constraints[flag_dual].dual_value
    #         flag_dual = flag_dual+1
    #         theta_closedual_1[i] = constraints[flag_dual].dual_value
    #         flag_dual = flag_dual+1
    
    
    #--------------------Second RIS group----------------------#
    
    # theta_second_dual = np.zeros([Qnum,1], dtype = 'complex_')

    # for i in range(0, RIS_group):
    #     if (RIS_element_control[i] == 0):
    #         for j in range(0, single_group_ele):
    #             theta_second_dual[i*single_group_ele+j] = constraints[i*single_group_ele+j].dual_value
                
                
    # theta_primary_SINR = constraints[len(constraints)-1]
    
    #--------------------Third RIS PU,SU----------------------#
    if (xonoff[0] == 1):
        barg1_S = np.linalg.norm(glvec_1_poweropt)**2/(Noise + (np.linalg.norm(glvec_1_poweropt)**2)*power_P)
        p_0 = P_k #+ P_R*np.sum(xonoff)*Tx_antRIS
        ee_S = Bandwidth*ma.log(1+barg1_S*power_SU.value, 2)/(mu*power_SU.value + p_0)
        PU_SINR = (np.linalg.norm(glvec_1_poweropt)**2)*power_P
    
    else:
        barg1_S = np.linalg.norm(glvec_1_poweropt_off)**2/(Noise + (np.linalg.norm(glvec_1_poweropt)**2)*power_P)
        p_0 = P_k #+ P_R*np.sum(xonoff)*Tx_antRIS
        ee_S = Bandwidth*ma.log(1+barg1_S*power_SU.value, 2)/(mu*power_SU.value + p_0)
        # PU_SINR = 
    
    
    theta_PU_SU_phasedual = np.zeros([Tx_antRIS,1], dtype = 'complex_')
    power_SU_dual = np.zeros([2,1], dtype=(float))
    
    for i in range(0, Tx_antRIS):
       theta_PU_SU_phasedual[i] = constraints_theta[i].dual_value
    
    power_SU_dual[0] = constraints_power[0].dual_value
    power_SU_dual[1] = constraints_power[1].dual_value
    PU_SINR_dual = constraints_power[2].dual_value
    
    prob_exitflag_theta = prob_theta.status
    prob_exitflag_power = prob_power.status
    prob_LowerBound = ee_S
    
    # thetavec_primal_1 = thetaopti_real_1.value
    # thetavec_primal_2 = thetaopti_real_2.value
    
    thetavec_PU_SU_primal = theta_PU_SU_RISphase.value
    
    # theta_primal_1 = thetaopti_real_1.value.conj()
    # theta_primal_2 = thetaopti_real_2.value.conj()
    # theta_primal = theta_primal_1 + theta_primal_2
    
    theta_PU_SU_primal = theta_PU_SU_RISphase.value.conj()
    
    # theta_primaldual_1 = thetaopt_dual_1
    # theta_primaldual_2 = thetaopt_dual_2
    # theta_primaldual = thetaopt_dual_1 + thetaopt_dual_2
    
    theta_PU_SU_primal_phasedual = theta_PU_SU_phasedual
    
    # theta_primalclose_1 = theta_closedual_1
    # theta_primalclose_2 = theta_closedual_2
    # theta_primalclose = theta_primalclose_1 + theta_closedual_2
    
    power_S_opt = power_SU.value
    
    # return Ulmar_1
    return power_S_opt,thetavec_PU_SU_primal, theta_PU_SU_primal, theta_PU_SU_primal_phasedual, power_SU_dual, PU_SINR_dual, prob_LowerBound, prob_exitflag_theta, prob_exitflag_power
    # return thetavec_PU_SU_primal, theta_PU_SU_primal, theta_PU_SU_primal_phasedual, power_SU_dual, PU_SINR_dual, prob_LowerBound, prob_exitflag
    
    
    