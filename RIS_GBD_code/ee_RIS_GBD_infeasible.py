# -*- coding: utf-8 -*-
import sys
import numpy as np
import numpy.matlib
import math as ma
import cvxpy as cvx

def RIS_GBD_infeasible(PathLoss_UserBS, PathLoss_UserRIS, PathLoss_RISBS, Tx_antBS, Tx_antRIS, RIS_Lnum, xonoff, power_P, theta_initial, P_max, Noise):
    
    
    glvec_1 = np.zeros([Tx_antBS,1], dtype = 'complex_')
    glvec_1[:,0] = PathLoss_UserBS[0,:]
    
    glvec_1_poweropt = glvec_1.conj().T
    glvec_1_poweropt_off = glvec_1.conj().T
    
    userris_1 = np.zeros([Tx_antRIS,1], dtype = 'complex_')
    risbs = np.zeros([Tx_antRIS, Tx_antBS], dtype = 'complex_')
    
    Ulmar_1 = np.zeros([Tx_antRIS,Tx_antBS], dtype = 'complex_')
    U_diag = np.eye(Tx_antRIS)
    
    for j in range(0, Tx_antBS):
        risbs[:,j] = PathLoss_RISBS[j,:,:]
    
    flag = 0
    for i in range(0,RIS_Lnum):
        userris_1[:,0] = PathLoss_UserRIS[:,0,:].reshape(32,)
        
        Ulmar_1[flag*Tx_antRIS : (flag+1)*Tx_antRIS,:] = np.dot(U_diag*(userris_1.conj().T), risbs)
        flag = flag+1
    
    infeasible_var_theta = cvx.Variable()
    infeasible_var_power = cvx.Variable()
    
    
    # theta_infeasible_1 = cvx.Variable(Qnum, complex=True)
    # theta_infeasible_2 = cvx.Variable(Qnum, complex=True)
    
    infeasible_theta_PUSU_phase = cvx.Variable(Tx_antRIS, complex=True)
    infeasible_powerSU = cvx.Variable()
    
    
    obj_infeasible_theta = cvx.Minimize(infeasible_var_theta)
    
    constraints_theta_infeasible = []
    # constraints_infeasible.append(infeasible_var >= 0)
    
    # theta_infeasible_num1 = 0
    # theta_infeasible_num2 = 0
    # theta_infeasible_close1 = 0
    # theta_infeasible_close2 = 0
    
    # for i in range(0, Qnum):
    #     if (RIS_element_control[i] ==0):
    #         constraints_infeasible.append(cvx.abs(theta_infeasible_2[i]) <= 1+infeasible_var)
    #         theta_infeasible_num2 = theta_infeasible_num2 + 1
    #         constraints_infeasible.append(theta_infeasible_1[i] == 0)
    #         theta_infeasible_close1 = theta_infeasible_close1 + 1
            
    #     else:
    #         constraints_infeasible.append(cvx.abs(theta_infeasible_1[i]) <= 1+infeasible_var)
    #         theta_infeasible_num1 = theta_infeasible_num1 + 1
    #         constraints_infeasible.append(theta_infeasible_2[i] == 0)
    #         theta_infeasible_close2 = theta_infeasible_close2 + 1
    
    # theta_second_infeasible = 0
    # obj_pri_1_inf = 0
    
    # for i in range(0,RIS_group):
    #     if (RIS_element_control[i] == 0):
    #         for j in range(0, single_group_ele):
    #             constraints_infeasible.append(cvx.abs(infeasible_thetasecond_2[i*single_group_ele+j]) <= 1 + infeasible_var)
    #             theta_second_infeasible +=1
            
    # for i in range(0,RIS_group):
    #     if (RIS_element_control[i] == 1):
    #         temp_vec_inf = (glvec_1 + (np.dot(Ulmar_11[i,:,:].conj().T, theta_initial[:,i])).reshape(8,1)).conj().T
    #         obj_pri_1_inf += 2*cvx.real(temp_vec_inf@Ulmar_11[i,:,:].conj().T@infeasible_thetasecond_2[i*single_group_ele:(i+1)*single_group_ele])
    
    for i in range(0, Tx_antRIS):
        constraints_theta_infeasible.append(cvx.abs(infeasible_theta_PUSU_phase[i]) <= 1 + infeasible_var_theta)
    
    constraints_theta_infeasible.append(infeasible_var_theta >= 0)
    
    prob_infeasible_theta = cvx.Problem(obj_infeasible_theta, constraints_theta_infeasible)
    prob_infeasible_theta.solve()
    
    
    theta_initial = infeasible_theta_PUSU_phase.value.conj()
    
    
    for i in range(0, RIS_Lnum):
        userris_1[:,0] = PathLoss_UserRIS[:,0,:].reshape(32,)
        
        thetal_1 = theta_initial[i*Tx_antRIS:(i+1)*Tx_antRIS]
        
        glvec_1_poweropt = glvec_1_poweropt + np.dot(np.dot(userris_1.conj().T,(xonoff[i]*np.diag(thetal_1))), risbs)    
    
    
    obj_infeasible_power = cvx.Minimize(infeasible_var_power)
    
    constraints_power_infeasible = []
    
    
    if (xonoff[0] == 1):
        obj_inf_SU = infeasible_powerSU
        
        obj_inf_PU = power_P
        
    else:
        obj_inf_SU  = (cvx.norm(glvec_1_poweropt_off)**2)*infeasible_powerSU
        
        obj_inf_PU = np.linalg.norm(glvec_1_poweropt)**2*power_P 
    
    
    
    
    constraints_power_infeasible.append(infeasible_powerSU + infeasible_var_power>= 0)
    constraints_power_infeasible.append(infeasible_powerSU  <= P_max + infeasible_var_power)
    constraints_power_infeasible.append(obj_inf_PU  >= 5*(obj_inf_SU + Noise) + infeasible_var_power)
    constraints_power_infeasible.append(infeasible_var_power >= 0)
    
    
    
    prob_infeasible_power = cvx.Problem(obj_infeasible_power, constraints_power_infeasible)
    prob_infeasible_power.solve()
    
    # total_infeasible_conslen = theta_infeasible_num1 + theta_infeasible_num2 + theta_infeasible_close1 + theta_infeasible_close2
    
    # thetaopt_infeasible_dual_1 = np.zeros([Qnum,1], dtype = 'complex_')
    # thetaopt_infeasible_dual_2 = np.zeros([Qnum,1], dtype = 'complex_')
    # theta_infeasible_closedual_1 = np.zeros([Qnum,1], dtype = 'complex_')
    # theta_infeasible_closedual_2 = np.zeros([Qnum,1], dtype = 'complex_')
    
    # flag_infeasible_dual = 1
    # for i in range(0, Qnum):
    #     if(RIS_element_control[i] == 1):
    #         thetaopt_infeasible_dual_1[i] = constraints_infeasible[flag_infeasible_dual].dual_value
    #         flag_infeasible_dual = flag_infeasible_dual + 1
    #         theta_infeasible_closedual_2[i] = constraints_infeasible[flag_infeasible_dual].dual_value
    #         flag_infeasible_dual = flag_infeasible_dual + 1
    #     else:
    #         thetaopt_infeasible_dual_2[i] = constraints_infeasible[flag_infeasible_dual].dual_value
    #         flag_infeasible_dual = flag_infeasible_dual + 1
    #         theta_infeasible_closedual_1[i] = constraints_infeasible[flag_infeasible_dual].dual_value
    #         flag_infeasible_dual = flag_infeasible_dual + 1
    

    theta_infeasible_PUSU_dual = np.zeros([Tx_antRIS,1], dtype = 'complex_')
    power_infeasible_SU_dual = np.zeros([2,1], dtype = float)
    
       
    # for i in range(0, RIS_group):
    #     if (RIS_element_control[i] == 0):
    #         for j in range(0, single_group_ele):
    #             theta_infeasiblesecond_dual[i*single_group_ele+j] = constraints_infeasible[i*single_group_ele+j].dual_value
                
    # theta_infeasible_primaSINR = constraints_infeasible[len(constraints_infeasible)-2].dual_value
    # zzz = constraints_infeasible[len(constraints_infeasible)-1].dual_value
    
    
    for i in range(0, Tx_antRIS):
        theta_infeasible_PUSU_dual[i] = constraints_theta_infeasible[i].dual_value
        
    power_infeasible_SU_dual[0] = constraints_power_infeasible[0].dual_value
    power_infeasible_SU_dual[1] = constraints_power_infeasible[1].dual_value
    
    infeasible_PU_SINR_dual = constraints_power_infeasible[2].dual_value
    
    # thetavec_primalinfeasible_1 = theta_infeasible_1.value
    # thetavec_primalinfeasible_2 = theta_infeasible_2.value
    
    thetavec_infeasible_PUSU = infeasible_theta_PUSU_phase.value
    
    # theta_primalinfeasible_1 = theta_infeasible_1.value.conj()
    # theta_primalinfeasible_2 = theta_infeasible_2.value.conj()
    # theta_primalinfeasible = theta_primalinfeasible_1 + theta_primalinfeasible_2
    
    theta_infeasible_PUSU = infeasible_theta_PUSU_phase.value.conj()
    
    # theta_primalinfeasible_dual_1 = thetaopt_infeasible_dual_1
    # theta_primalinfeasible_dual_2 = thetaopt_infeasible_dual_2
    

    power_infeasible_SU_opt = infeasible_powerSU.value
    
    # theta_primalinfeasible_closedual_1 = theta_infeasible_closedual_1
    # theta_primalinfeasible_closedual_2 = theta_infeasible_closedual_2
    
    return thetavec_infeasible_PUSU, theta_infeasible_PUSU, theta_infeasible_PUSU_dual, power_infeasible_SU_opt, power_infeasible_SU_dual,infeasible_PU_SINR_dual
    
    