# -*- coding: utf-8 -*-

import cplex
from cplex.exceptions import CplexError
from cplex.exceptions import CplexSolverError
import sys
import numpy as np
import copy
import math
import time


class Master_class:
    my_prob = cplex.Cplex()
    init_flag = False

    def __init__(self, Bandwidth_in,Noise_in,P_max_in,P_k,P_R,P_A,P_B,mu,PathLoss_UserBS_in,PathLoss_UserRIS_in,PathLoss_RISBS_in,
                 xonoff,Num_User_in,Tx_antBS_in,RIS_Lnum_in,Tx_antRIS_in,powerini, power_P_in, power_S_in, thetamarini_in,
                 RIS_element_control, PU_SINR_min_in):
        
        self.my_prob = cplex.Cplex()
        self.init_flag == False

        self.my_prob.objective.set_sense(self.my_prob.objective.sense.maximize)
        global PathLoss_UserBS, PathLoss_UserRIS, PathLoss_RISBS, Num_User, Tx_antBS, RIS_Lnum, \
            Tx_antRIS, Bandwidth, Noise, P_max, power_P, PU_SINR_min

        PathLoss_UserBS = PathLoss_UserBS_in
        PathLoss_UserRIS = PathLoss_UserRIS_in
        PathLoss_RISBS = PathLoss_RISBS_in
        Num_User = Num_User_in
        Tx_antBS = Tx_antBS_in
        RIS_Lnum = RIS_Lnum_in
        Tx_antRIS = Tx_antRIS_in
        Bandwidth = Bandwidth_in
        Noise = Noise_in
        P_max = P_max_in
        power_P = power_P_in
        PU_SINR_min = PU_SINR_min_in
        # single_group_ele = single_group_ele_in
        # RIS_group = RIS_group_in



        

    def populatebyrow(self, my_obj, my_lb, my_ub, my_colnames, my_new_row, my_new_rhs, RIS_Lnum):
        if self.init_flag == False:
            self.init_flag = True
            # my_ctype = ""
            # for index in range(0, Tx_antRIS):
            #     my_ctype = my_ctype + "I"
            # my_ctype = my_ctype + "C"

            # self.my_prob.variables.add(
            #     obj=my_obj, lb=my_lb, ub=my_ub, names=my_colnames, types=my_ctype)
            
            # my_colnames = []

            # index_i = 1
            # index_j = 0
            # for index in range(0, K*L):
            #     index_j = index_j+1
            #     if index_j > L:
            #         index_j = 1
            #         index_i = index_i + 1
            #     str_temp = "rho"+repr(index_i)+"_"+repr(index_j)
            #     my_colnames.append(str_temp)
            
            # for index in range(0, Tx_antRIS):
            #     str_temp = "RISelement"+repr(index_i)
            #     my_colnames.append(str_temp)
            #     index_i = index_i + 1

            # my_colnames.append("Psi")
            
            # t_sense = "L"
            # t_new_rhs = [1]
            # for index_k in range(0,K):
            #     variables_coef = np.zeros(K*L+1, dtype=np.float_)
            #     for index_l in range(0,L):
            #         variables_coef[index_k*L+index_l]=1.0
                
            #     t_new_row = [my_colnames, variables_coef.tolist()]
            #     t_new_row = [t_new_row]   
                
            #     #print(t_new_row,t_sense,t_new_rhs)
            #     self.my_prob.linear_constraints.add(
            #             lin_expr=t_new_row, senses=t_sense, rhs=t_new_rhs)

        my_new_row = [my_new_row]
        my_sense = "G"

        #print("my_row:",my_new_row)        
        #print("my_sense:",my_sense)
        #print("rhs:",my_new_rhs)
        self.my_prob.linear_constraints.add(
            lin_expr=my_new_row, senses=my_sense, rhs=my_new_rhs)

    def lplex1(self, my_obj, my_lb, my_ub, my_colnames, my_new_row, my_new_rhs, Tx_antRIS, numsol, ret_flag):
        
        try:
            handle = self.populatebyrow(
                my_obj, my_lb, my_ub, my_colnames, my_new_row, my_new_rhs, RIS_Lnum)
        except CplexSolverError:
            print("Exception raised during populate")
            return

        if ret_flag == True:
            #print("test_numsol",numsol)
            self.my_prob.parameters.mip.pool.capacity.set(numsol)
            self.my_prob.parameters.mip.pool.replace.set(1)
            #time_start =time.process_time()
            self.my_prob.populate_solution_pool()
            #time_use = time.process_time() 


            num_solution = self.my_prob.solution.pool.get_num()
            #print("The solution pool contains %d solutions." % num_solution)
            # meanobjval = cpx.solution.pool.get_mean_objective_value()

            numsol_real = min(numsol, num_solution)
            sol_pool = []
        
            obj_temp = np.zeros(numsol_real)
            for i in range(numsol_real):
                obj_temp[i] = self.my_prob.solution.pool.get_objective_value(i) 
            new_index = sorted(range(len(obj_temp)), key=lambda k: obj_temp[k])
            print(obj_temp)
            print(new_index)

            for j in range(numsol_real):
                i = new_index[j]
                objval_i = self.my_prob.solution.pool.get_objective_value(i)
                x_i = self.my_prob.solution.pool.get_values(i)
                nb_vars=len(x_i)
                sol = []
                for k in range(nb_vars):
                    sol.append(x_i[k])
                sol_pool.append(sol)
                print("object:",i,objval_i)
                print("value:",i,x_i)

            print("pools =",sol_pool)

            # Print information about the incumbent
            print("Objective value of the incumbent  = ",
                self.my_prob.solution.get_objective_value())

            r_Psi = self.my_prob.solution.get_objective_value()

            #r_rho = np.ones([numsol_real,K, L], dtype=int)
            r_rho = np.ones([numsol_real, RIS_Lnum], dtype=int)
            # for i in range(numsol_real):
            #     x_i = sol_pool[i]
            #     index_i = 0
            #     index_j = -1
            #     for index in range(0, K*L):
            #         index_j = index_j+1
            #         if index_j >= L:
            #             index_j = 0
            #             index_i = index_i + 1
            #         r_rho[i,index_i,index_j] = x_i[index]
            
            for i in range(0, numsol_real):
                x_i = sol_pool[i]
                
                for index in range(0, RIS_Lnum):
                    r_rho[i, index] = x_i[index]
                    
            print("r_rho",r_rho)
            print("r_psi",r_Psi)


            return r_rho, r_Psi, numsol_real#, time_use-time_start
        else:
            return -1, -1, -1, -1

    def master_sol(self,feasible_flag_multi,thetavec_PU_SU_primal_multi,thetavec_PU_SU_primaldual_multi,
                   power_SU_primal_multi, power_SU_primal_dual_multi, PU_SINR_multi, PU_SINR_primal_dual_multi, infeasible_theta_multi,
                   infeasible_power_multi, numsol):
        #N = K*L + 1
        
        my_obj = np.zeros(RIS_Lnum+1, dtype=np.float_)
        my_obj[Tx_antRIS] = 1.0
        # print(my_obj)
        my_lb = np.zeros(RIS_Lnum+1, dtype=np.float_)
        my_ub = np.ones(RIS_Lnum+1, dtype=np.float_)

        my_lb[RIS_Lnum] = -cplex.infinity
        my_ub[RIS_Lnum] = cplex.infinity

        my_colnames = []

        # index_i = 1
        # index_j = 0
        # for index in range(0, K*L):
        #     index_j = index_j+1
        #     if index_j > L:
        #         index_j = 1
        #         index_i = index_i + 1
        #     str_temp = "rho"+repr(index_i)+"_"+repr(index_j)
        #     my_colnames.append(str_temp)
        index_i = 1
        for index in range(0, RIS_Lnum):
            str_temp = "RIS_SU_onoff"+repr(index_i)
            my_colnames.append(str_temp)
            index_i = index_i + 1

        my_colnames.append("Psi")

        this_num = feasible_flag_multi.shape[0]
        # print(feasible_flag_multi)
        # print("test")
        # print(this_num)
        for num_index in range(this_num):
            feasible_flag = feasible_flag_multi[num_index]
            # p_D = p_D_multi[num_index]
            # eta = eta_multi[num_index]
            # alpha = alpha_multi[num_index]
            # lameda = lameda_multi[num_index]
            # miu = miu_multi[num_index]
            # nu = nu_multi[num_index]
            # theta = theta_multi[num_index]
            
            thetavec_PU_SU_primal = thetavec_PU_SU_primal_multi[num_index]
            power_SU_primal = power_SU_primal_multi[num_index]
            PU_SINR = PU_SINR_multi[num_index]
            #thetavec_primal_2 = thetavec_primal_2_multi[num_index]
            thetavec_PU_SU_primaldual = thetavec_PU_SU_primaldual_multi[num_index]
            #theta_primaldual_2 = theta_primaldual_2_multi[num_index]
            power_SU_primal_dual = power_SU_primal_dual_multi[num_index]
            PU_SINR_primal_dual = PU_SINR_primal_dual_multi[num_index]
            
            infeasible_theta = infeasible_theta_multi[num_index]
            infeasible_power = infeasible_power_multi[num_index]

            # feasible_flag_multi,p_D_multi,eta_multi,alpha_multi,lameda_multi,miu_multi,nu_multi,theta_multi,numsol

            if feasible_flag:
                variables_names = copy.deepcopy(my_colnames)
                variables_coef = np.zeros(RIS_Lnum+1, dtype=np.float_)
                rho_coef = np.zeros([RIS_Lnum, 1], dtype=np.float_)
                new_rhs = 0.0

                # ψ 
                variables_coef[RIS_Lnum] = 1.0

                # constant
                # new_rhs -= eta
                
                new_rhs -= power_SU_primal_dual[0]*power_SU_primal
                
                new_rhs += power_SU_primal_dual[1]*power_SU_primal
                new_rhs -= power_SU_primal_dual[1]*P_max
                
                

                for index_l in range(0, Tx_antRIS):
                    new_rhs += thetavec_PU_SU_primal[index_l]*(abs(thetavec_PU_SU_primaldual[index_l]))
                    # new_rhs += theta_primaldual_2[index_l]*(abs(thetavec_primal_2[index_l])-1)
                    new_rhs -= thetavec_PU_SU_primaldual[index_l]*1
                    # new_rhs -= theta_primaldual_2[index_l]*1
                
                
                
                glvec_1_total = 0
                glvec_1 = np.zeros([Tx_antBS,1], dtype = 'complex_')
                glvec_1[:,0] = PathLoss_UserBS[0,:]
                
                userris_1 = np.zeros([Tx_antRIS, 1], dtype = 'complex_')
                risbs = np.zeros([Tx_antRIS,Tx_antBS], dtype = 'complex_')
                
                
                for index_bs in range(0, Tx_antBS):
                    glvec_1_total += glvec_1[index_bs]**2
                
                
                
                ris_ref = np.zeros([Tx_antBS,1], dtype = 'complex_')
                for i in range(0, RIS_Lnum):
                    userris_1[:,0] = PathLoss_UserRIS[:,0,:].reshape(32,)
                    thetal_1 = thetavec_PU_SU_primal[i*Tx_antRIS:(i+1)*Tx_antRIS,0]
                    ris_ref = np.dot(np.dot(userris_1.conj().T,(np.diag(thetal_1))), risbs)
                
                
                for index_ris in range(0,Tx_antRIS):
                    rho_coef += 2*glvec_1[index_ris]*ris_ref[index_ris] + ris_ref[index_ris]**2
                    
                
                new_rhs -= PU_SINR_primal_dual*(glvec_1_total + rho_coef)*power_P
                new_rhs += PU_SINR_primal_dual*PU_SINR_min*(glvec_1_total*power_SU_primal + Noise)
                
                rho_coef -= PU_SINR_primal_dual*(PU_SINR_min - infeasible_power)*rho_coef*power_SU_primal 
                
                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         new_rhs += lameda[index_l]*p_D[index_k][index_l]

                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         new_rhs += miu[index_k][index_l] * \
                #             (p_D[index_k][index_l])

                # for index_k in range(0, K):
                #     for index_l in range(0,L):
                #         new_rhs -= nu[index_k][index_l]*p_D[index_k][index_l]

                # for index_l in range(0, L):
                #     new_rhs += theta[index_l]*eta

                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         new_rhs -= theta[index_l]*math.log(1+p_D[index_k][index_l]/(
                #             a[index_k][index_l]+b[index_k][index_l]*p_D[index_k][index_l]), 2.0)
                
                #------element on/off cosfficient------#
                
                # for index_RIS in range(0, Tx_antRIS):
                #     rho_coef[index_RIS] = glvec_1

                # #print("rho_coef0:",rho_coef)
                # for index_l in range(0, L):
                #     for index_k in range(0, K):
                #         rho_coef[index_k][index_l] = miu[index_k][index_l] * \
                #             p_max_k_l[index_k][index_l]

                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         index_temp = index_k*L + index_l
                #         variables_coef[index_temp] = rho_coef[index_k][index_l]
                #         # print(index_k,index_l,variables_names[index_temp])
                        
                #------element on/off cosfficient------#

                new_row = [variables_names, variables_coef.tolist()]

                new_rhs=new_rhs.tolist()

                #print("feasible")
                #print(my_obj,my_lb,my_ub,my_colnames,new_row,new_rhs,K,L,numsol,num_index,num_index==this_num-1)
                r_rho, r_Psi, r_solnum, time_use = self.lplex1(
                    my_obj, my_lb, my_ub, my_colnames, new_row, new_rhs, Tx_antRIS, numsol, num_index==this_num-1)

                #res_file = open('result.txt',mode='a+')
                #res_file.writelines(['coef:\t',str(variables_coef),'\nflag:\t',str(num_index==this_num-1),'\ncut_num_sofar\t',str(self.my_prob.linear_constraints.get_num()),'\n'])
                #res_file.writelines(['\n'])
                #res_file.close()

            else:
                variables_names = copy.deepcopy(my_colnames)
                variables_coef = np.zeros(RIS_Lnum+1, dtype=np.float_)
                rho_coef = np.zeros([RIS_Lnum, 1], dtype=np.float_)
                new_rhs = 0.0
                
                
                new_rhs -= power_SU_primal_dual[0]*(power_SU_primal+infeasible_power)
                
                new_rhs += power_SU_primal_dual[1]*power_SU_primal
                new_rhs -= power_SU_primal_dual[1]*P_max
                new_rhs -= power_SU_primal_dual[1]*infeasible_power
                
                
                
                for index_l in range(0, Tx_antRIS):
                    new_rhs += thetavec_PU_SU_primal[index_l]*(abs(thetavec_PU_SU_primaldual[index_l]))
                    new_rhs -= thetavec_PU_SU_primaldual[index_l]*1
                    new_rhs -= thetavec_PU_SU_primaldual[index_l]*infeasible_theta
                
                
                
                glvec_1_total = 0
                glvec_1 = np.zeros([Tx_antBS,1], dtype = 'complex_')
                glvec_1[:,0] = PathLoss_UserBS[0,:]
                
                userris_1 = np.zeros([Tx_antRIS, 1], dtype = 'complex_')
                risbs = np.zeros([Tx_antRIS,Tx_antBS], dtype = 'complex_')
                
                
                for index_bs in range(0, Tx_antBS):
                    glvec_1_total += glvec_1[index_bs]**2
                
                
                
                ris_ref = np.zeros([Tx_antBS,1], dtype = 'complex_')
                for i in range(0, RIS_Lnum):
                    userris_1[:,0] = PathLoss_UserRIS[:,0,:].reshape(32,)
                    thetal_1 = thetavec_PU_SU_primal[i*Tx_antRIS:(i+1)*Tx_antRIS,0]
                    ris_ref = np.dot(np.dot(userris_1.conj().T,(np.diag(thetal_1))), risbs)
                
                
                for index_ris in range(0,Tx_antRIS):
                    rho_coef += 2*glvec_1[index_ris]*ris_ref[index_ris] + ris_ref[index_ris]**2
                
                
                # new_rhs -= PU_SINR_primal_dual*PU_SINR
                # new_rhs -= PU_SINR_primal_dual*infeasible_power
                new_rhs += PU_SINR_primal_dual*(PU_SINR_min-infeasible_power)*(glvec_1_total*power_SU_primal + Noise)
                new_rhs -= PU_SINR_primal_dual*(glvec_1_total + rho_coef)*power_P
                
                
                variables_coef[0] = -PU_SINR_primal_dual*(PU_SINR_min - infeasible_power)*rho_coef*power_SU_primal 

                # constant
                # for index_l in range(0, L):
                #     new_rhs -= lameda[index_l]*(P_max_D+alpha)
                
                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         new_rhs += lameda[index_l]*p_D[index_k][index_l]

                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         new_rhs += miu[index_k][index_l] * \
                #             (p_D[index_k][index_l]-alpha)

                # for index_k in range(0, K):
                #     for index_l in range(0,L):
                #         new_rhs -= nu[index_k][index_l]*(p_D[index_k][index_l]+alpha)

                # for index_l in range(0, L):
                #     new_rhs += theta[index_l]*(eta-alpha)
                
                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         new_rhs -= theta[index_l]*math.log(1+p_D[index_k][index_l]/(
                #             a[index_k][index_l]+b[index_k][index_l]*p_D[index_k][index_l]), 2.0)
                        
                # for index_l in range(0, L):
                #     for index_k in range(0, K):
                #         rho_coef[index_k][index_l] = miu[index_k][index_l] * \
                #             p_max_k_l[index_k][index_l]

                # for index_k in range(0, K):
                #     for index_l in range(0, L):
                #         index_temp = index_k*L + index_l
                #         variables_coef[index_temp] = rho_coef[index_k][index_l]
                #         # print(index_k,index_l,variables_names[index_temp])

                new_row = [variables_names, variables_coef.tolist()]

                print("infeasible")
                r_rho, r_Psi, r_solnum, time_use = self.lplex1(
                    my_obj, my_lb, my_ub, my_colnames, new_row, new_rhs, Tx_antRIS, numsol, num_index==this_num-1)
        return r_rho, r_Psi, r_solnum#, time_use