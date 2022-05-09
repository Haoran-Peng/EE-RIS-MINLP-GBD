# -*- coding: utf-8 -*-
import sys
import numpy as np
import numpy.matlib
import math as ma
import matplotlib.pyplot as plt

from ee_RIS_GBD_para import userpara
from ee_RIS_GBD_twouser_EE import RIS_EE_twousers
from ee_RIS_GBD_primal import RIS_GBD_primal
from ee_RIS_GBD_infeasible import RIS_GBD_infeasible
from ee_RIS_theta_SCAopt import RIS_SCA_thetaopt


RIS_Lnum = 1
Num_User = 2
Tx_antBS = 8
Tx_antRIS = 32
single_group_ele = 4
RIS_group = int(Tx_antRIS/single_group_ele)
power_S_max = 50


power = 30
power=pow(10, (power/10))/pow(10,3)

xonoff = np.ones([RIS_Lnum,1])
# xonoff = np.zeros([RIS_Lnum,1])
thetamarini = np.ones([Tx_antRIS,1], dtype = 'complex_')
thetamarini_pri1 = np.ones([Tx_antRIS,1], dtype = 'complex_')
thetamarini_pri2 = np.ones([Tx_antRIS,1], dtype = 'complex_')

RIS_element_control = np.zeros([RIS_group,1])
powerini = power
power_P = power
power_S = 0
Qnum = int(np.sum(xonoff)*Tx_antRIS)
theta_SCAopti_1 = np.zeros([single_group_ele, RIS_group], dtype = 'complex_')
# print(theta_SCAopti_1.shape)

[Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,PathLoss_UserRIS_group_1,PathLoss_UserRIS_group_2, PathLoss_RISBS_group]\
    =userpara(Num_User, Tx_antBS, RIS_Lnum, Tx_antRIS, single_group_ele)

#------SCA-------#

# theta_SCAopti_1 = RIS_SCA_thetaopt(powerini,thetamarini,RIS_element_control,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS,\
#                       PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)

# # print(theta_SCAopti_1)
# # print(theta_SCAopti_1.shape)

# [ee_S, rate_S, SINR_P, ee_S_off, rate_S_off, SINR_P_Soff] = RIS_EE_twousers(theta_SCAopti_1, power_P, power_S,xonoff, RIS_element_control,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,
#                     PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)
# index_increase = 0
# flag_data = np.zeros([6,51], dtype=(float))

# flag_data[0,0] = ee_S
# flag_data[1,0] = rate_S
# flag_data[2,0] = SINR_P
# flag_data[3,0] = ee_S_off
# flag_data[4,0] = rate_S_off
# flag_data[5,0] = SINR_P_Soff


# for idx_power in range(0, power_S_max):
#     # print(RIS_element_control)
#     power_S = idx_power
#     power_S=pow(10, (power_S/10))/pow(10,3)
#     # print(power_S)
    
#     theta_SCAopti_1 = RIS_SCA_thetaopt(powerini,thetamarini,RIS_element_control,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS,\
#                       PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)
        
#     [ee_S, rate_S, SINR_P, ee_S_off, rate_S_off, SINR_P_Soff] = RIS_EE_twousers(theta_SCAopti_1, power_P, power_S, xonoff, RIS_element_control,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,
#                     PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)
    
#     flag_data[0,idx_power+1] = ee_S
#     flag_data[1,idx_power+1] = rate_S
#     flag_data[2,idx_power+1] = SINR_P
#     flag_data[3,idx_power+1] = ee_S_off
#     flag_data[4,idx_power+1] = rate_S_off
#     flag_data[5,idx_power+1] = SINR_P_Soff

#------SCA-------#



# print(flag_data[0])
# print(flag_data[1])
# print(flag_data[2])
# print(flag_data[3])
# print(flag_data[4])
# print(flag_data[5])

xx = np.arange(51)


# plt.plot(xx,flag_data[0], label="ee_S")
# plt.plot(xx,flag_data[1], label="rate_S")
# plt.plot(xx,flag_data[2], label="SINR_P")
# plt.ylabel('ee')
# plt.xlabel('power_SU')
# plt.legend()
# plt.show()


# plt.plot(xx,flag_data[3], label="ee_S_off")
# plt.plot(xx,flag_data[4], label="rate_S_off")
# plt.ylabel('ee')
# plt.xlabel('power_SU')
# plt.legend()
# plt.show()



# plt.plot(xx,flag_data[2], label="SINR_P")
# plt.plot(xx,flag_data[5], label="SINR_P_Soff")
# plt.ylabel('SINR_PU')
# plt.xlabel('power_SU')
# plt.legend()
# plt.show()



# plt.plot(xx,flag_data[0], label="ee_S")
# plt.plot(xx,flag_data[3], label="ee_S_off")
# plt.plot(xx,flag_data[1], label="rate_S")
# plt.plot(xx,flag_data[4], label="rate_S_off")

# plt.ylabel('ee')
# plt.xlabel('power_SU')
# plt.legend()
# plt.show()

# plt.ylabel('ee')
# plt.xlabel('power')
# plt.legend()
# plt.show()

#------------------GBD-----------------#

power_S = power

# [power_S_opt,thetavec_PU_SU_primal, theta_PU_SU_primal, theta_PU_SU_primal_phasedual, power_SU_dual, PU_SINR_dual, prob_LowerBound, prob_exitflag_theta, prob_exitflag_power]\
# = RIS_GBD_primal(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoff,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,powerini,power_P,power_S,thetamarini,RIS_element_control)

# print(power_S_opt)
# print(thetavec_PU_SU_primal)
# print(theta_PU_SU_primal)
# print(theta_PU_SU_primal_phasedual)
# print(power_SU_dual)
# #print(PU_SINR_dual)
# print(prob_exitflag_theta)
# print(prob_exitflag_power)

[thetavec_infeasible_PUSU, theta_infeasible_PUSU, theta_infeasible_PUSU_dual, power_infeasible_SU_opt, power_infeasible_SU_dual,infeasible_PU_SINR_dual]\
  = RIS_GBD_infeasible(PathLoss_UserBS, PathLoss_UserRIS, PathLoss_RISBS, Tx_antBS, Tx_antRIS, RIS_Lnum, xonoff, power_P, thetamarini, P_max, Noise)

print(thetavec_infeasible_PUSU)
print(theta_infeasible_PUSU)
print(theta_infeasible_PUSU_dual)
print(power_infeasible_SU_opt)
print(power_infeasible_SU_dual)
print(infeasible_PU_SINR_dual)



# [thetavec_primal_1, thetavec_primal_2, theta_primal, theta_primal_1, theta_primal_2, theta_primaldual, theta_primaldual_1, theta_primaldual_2, theta_primalclose, theta_primalclose_1, theta_primalclose_2, prob_LowerBound, prob_exitflag]\
# = RIS_GBD_primal(Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,xonoff,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS,powerini,thetamarini,thetamarini_pri1,thetamarini_pri2,RIS_element_control)

# [thetavec_primalinfeasible_1, thetavec_primalinfeasible_2, theta_primalinfeasible, theta_primalinfeasible_1, theta_primalinfeasible_2, theta_primalinfeasible_dual_1, theta_primalinfeasible_dual_2, theta_primalinfeasible_closedual_1, theta_primalinfeasible_closedual_2]\
#  = RIS_GBD_infeasible(Qnum,RIS_element_control)


# # [ee,ee1,ee2,rate,rate1,rate2]= RIS_EE_twousers(thetamar,thetamar_1,thetamar_2,power,xonoff,Bandwidth,Noise,P_max,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS,PathLoss_UserRIS,PathLoss_RISBS,Num_User,Tx_antBS,RIS_Lnum,Tx_antRIS)

# print(thetavec_primalinfeasible_1)
# print(theta_primalinfeasible_1)
# print(thetavec_primalinfeasible_2)
# print(theta_primalinfeasible_2)
# print(theta_primal)
# print(theta_primaldual_1)
# print(theta_primaldual_2)


# D = np.array([3,2,6,5])
# E = np.ones([3,3])
# F = np.eye(5)
# print(D.shape)
# print(np.diag(D))