# @author: Genhui Zheng
# @institution: University of Texas at Austin
# @Email: genhuizheng@utexas.edu
from cmath import inf
import numpy as np
import pandas as pd

# multiply one numpy list by one fix calue 
def multiple_list_for_np_list(delta_output, prob):
    for i in range(len(delta_output)):
        delta_output[i] = delta_output[i]*prob
    return delta_output

# add one numpy list to the another
def update_list_for_np_list(delta_output, current_infect_POP):
    for i in range(len(current_infect_POP)):
        current_infect_POP[i] = delta_output[i] + current_infect_POP[i]
    return current_infect_POP

#whether it is longer than the duration of infection
def whether_T_more_than_recovery_period(recorder_list, duration_of_disease):
    return len(recorder_list) >= duration_of_disease

#caculate all the infected people
def all_infect_sum(infectious_population_list):
    ALL_INFECT_POPULATION = 0
    for i in range(0, len(infectious_population_list)):
        ALL_INFECT_POPULATION += np.sum(infectious_population_list[i])
    return ALL_INFECT_POPULATION

#summarize the data_record
def parameter_assemble_for_recording(S1, S2, infectious_population_list, recovery_population_list, out_list, lambda_list):
    return [[i.copy() for i in S1], S2.copy(), [i.copy() for i in infectious_population_list], [i.copy() for i in recovery_population_list], [i.copy() for i in out_list], lambda_list.copy()]

#generate a vacant record
def zero_generation_for_recording():
    return [[np.array([[0.0, 0.0, 0.0]]).copy() for i in range(0, 120)],
        np.array([[0., 0., 0.]]), 
        [[one * np.array([[0.8, 0.1, 0.1]]) * i for i in [0.5, 0.5]] for one in [0, 0, 0, 0]],
        [[one * np.array([[0.8, 0.1, 0.1]]) * i for i in [0.5, 0.5]] for one in [0, 0, 0, 0]],
        [[one * np.array([[0.8, 0.1, 0.1]]) * i for i in [0.5, 0.5]] for one in [0, 0, 0, 0]],
        [0,0,0,0,0]]

#recording the data
def recorder(recorder_list, one_record, duration_of_disease):
    if not whether_T_more_than_recovery_period(recorder_list, duration_of_disease):
        recorder_list.append(one_record.copy())
        return recorder_list, zero_generation_for_recording()
    else:
        recorder_list.append(one_record.copy())
        temp = recorder_list.pop(0)
        return recorder_list, temp.copy()

#calculate the change of the population owning maternal antibody
def delta_susceptible_for_mother_antibody(ALL_INFECT, lambda_1_m, infect_ratio_in_secrete_type, S1):
    temp = ALL_INFECT * lambda_1_m * infect_ratio_in_secrete_type.copy() * S1.copy()
    for i in range(0, len(temp)):
        for j in range(0, len(temp[i][0])):
            if temp[i][0][j] < 0 and temp[i][0][j] + S1[i][0][j] < 0:
                temp[i][0][j] = -S1[i][0][j].copy()
            elif temp[i][0][j] > 0 and temp[i][0][j] > S1[i][0][j]:
                temp[i][0][j] = S1[i][0][j].copy()
            else:
                continue
    return temp

#calculate the change of the suspectible population
def delta_susceptible_with_no_anti(ALL_INFECT, lambda_1_n, 
    infect_ratio_in_secrete_type, S2):
    temp = ALL_INFECT * lambda_1_n * infect_ratio_in_secrete_type * S2.copy()
    for i in range(0, len(temp[0])):
        if temp[0][i] < 0 and temp[0][i] + S2[0][i] < 0:
            temp[0][i] = -S2[0][i].copy()
        elif temp[0][i] > 0 and temp[0][i] > S2[0][i]:
            temp[0][i]= S2[0][i].copy()
        else:
            continue
    return temp

#calculate the sum of population owning maternal antibody
def sum_susceptible_for_mother_antibody(S1):
    SUM_S = np.array([[0.0, 0.0, 0.0]])
    for one in S1:
        SUM_S += one
    return SUM_S

#calculate the sum of the population under the 1st infection
def one_time_for_1st_infect(ALL_INFECT, lambda_1_m, lambda_1_n,infect_ratio_in_secrete_type, S1, S2, light_severe_distribution_list):
    temp = sum_susceptible_for_mother_antibody(delta_susceptible_for_mother_antibody(
            ALL_INFECT, lambda_1_m, infect_ratio_in_secrete_type, S1)) + \
        delta_susceptible_with_no_anti(ALL_INFECT, lambda_1_n, infect_ratio_in_secrete_type, S2)
    return [temp.copy() * i for i in light_severe_distribution_list]

#calculate the sum of the population under infection after recovery, not accounting for the change of different symptoms
def calculate_recovery_infected_pop(ALL_INFECT, lambda_infect, infect_ratio_in_secrete_type, recovery_pop):
    temp = ALL_INFECT * lambda_infect * infect_ratio_in_secrete_type.copy() * recovery_pop
    for j in range(0, len(temp)):
        for i in range(0, len(temp[j][0])):
            if temp[j][0][i] < 0 and temp[j][0][i] + recovery_pop[j][0][i] < 0:
                temp[j][0][i] = -recovery_pop[j][0][i].copy()
            elif temp[j][0][i] > 0 and temp[j][0][i] > recovery_pop[j][0][i]:
                temp[j][0][i]= recovery_pop[j][0][i].copy()
            else:
                continue
    return temp

#calculate the sum of the population under infection after recovery, accounting for the change of different symptoms
def calculate_infect_again_pop(ALL_INFECT, lambda_infect, infect_ratio_in_secrete_type, recovery_pop, matrix_of_light_severe):
    temp = calculate_recovery_infected_pop(ALL_INFECT, lambda_infect, infect_ratio_in_secrete_type, recovery_pop)
    final_pop = [np.array([[0., 0., 0.]]), np.array([[0., 0., 0.]])]
    for i in range(0, 2):
        for j in range(0, 2):
            final_pop[i] += temp[j] * matrix_of_light_severe[j][i]
    return final_pop

def run_one_model(input_overall_time, input_lambda, n0, n1, n2, n3, season_rate):
    ##parameters
    # mother_anti_rate: the birth rate of children with maternal antibodies evey day
    mother_anti_rate = 15
    #the duration time of maternal antibody for one child
    duration_of_mother_anti = 30*4
    #the duration of first infection
    duration_of_disease_1st = 7
    #the duration of second or above infection
    duration_of_disease_2nd_after = 4
    #the probability that every infected person infects the one with maternal antibody 
    lambda_1_m = input_lambda / n0
    #the probability that every infected person infects the suspectible person
    lambda_1_n = input_lambda
    #the probability that every infected person infects the person recovering from infection (recovery)
    """range from the 2nd to the 4th"""
    lambda_list_from_2_to_4 = [input_lambda * n1, input_lambda * n2, input_lambda * n3]
    #the probability that one person is IMMUNE to the virus after recovery
    """range from the 2nd to the 4th"""
    prob_of_non_infect_after_recovery_list = [0.52, 0.32/0.48, 0.13/0.16, 0.03/0.03]
    #the distribution proportion of severe and light in the crowd after infection
    light_severe_distribution_list = [72/85, 13/85]
    #the distribution proportion of secrete type, weak secrete and not secrete in the crowd
    secrete_type_distribution = np.array([[75/97, 14/97, 8/97]])
    #the infected proportion of secrete type, weak secrete and not secrete relative to the secrete type
    infect_ratio_in_secrete_type = np.array([[1, 0.38, 0.23]])
    #the transition matrix of different symptoms (light severe)
    transfer_p = (72 - 3/19 *85)/59
    matrix_of_light_severe= [[transfer_p, 1 - transfer_p],[1, 0]]
    #the rate between real R0 in fall or spring seaon and R0
    spring_fall_rate = season_rate
    #the rate between R0 and real R0 in summer or winter seaon
    summer_winter_rate = season_rate


    #overall lasting time of the model 
    overall_time = input_overall_time

    ##initial value of model parameters
    #the initial number of population owning the ANTIBODY from their moms 
    S1 = [mother_anti_rate * (secrete_type_distribution.copy()) for i in range(0, duration_of_mother_anti)]
    #the initial number of SUSPECTIBLE population, but not infected
    S2 = duration_of_mother_anti * mother_anti_rate * 14 *(secrete_type_distribution.copy())

    #the initial number of INFECTED population
    #ranges from 1st infected to 4th infected
    infect_list = [[one * secrete_type_distribution.copy() * i for i in light_severe_distribution_list] for one in [1, 0, 0, 0]]
    #the initial number of RECOVERY population recovering from infection
    #ranges from 1st infected to 4th infected
    recovery_list = [[one * secrete_type_distribution.copy() * i for i in light_severe_distribution_list] for one in [0, 0, 0, 0]]
    #the initial number of the population who is IMMUNE to the virus after recovery
    #ranges from 1st infected to 4th infected
    out_list = [[one * secrete_type_distribution.copy() * i for i in light_severe_distribution_list] for one in [0, 0, 0, 0]]

    disease_1st_recorder_list = []
    disease_2nd_recorder_list = []
    #the summary of data record
    data_list = []   

    for T_current in range(0, overall_time):

        current_record = parameter_assemble_for_recording(S1.copy(), S2.copy(), infect_list, recovery_list, out_list, [lambda_1_m, lambda_1_n] + lambda_list_from_2_to_4)

        ALL_INFECT = all_infect_sum(infect_list)
        disease_1st_recorder_list, N_recovery_record = recorder(disease_1st_recorder_list, current_record.copy(), duration_of_disease_1st)
        disease_2nd_recorder_list, I_2_recovery_record = recorder(disease_2nd_recorder_list, current_record.copy(), duration_of_disease_2nd_after)


        new_infect_record = [one_time_for_1st_infect(ALL_INFECT, lambda_1_m, lambda_1_n, infect_ratio_in_secrete_type.copy(), 
                current_record[0], current_record[1], light_severe_distribution_list.copy()).copy(),]
        new_recover_record = [multiple_list_for_np_list(
                one_time_for_1st_infect(all_infect_sum(N_recovery_record[2]), lambda_1_m, lambda_1_n, infect_ratio_in_secrete_type, 
                    N_recovery_record[0], N_recovery_record[1], light_severe_distribution_list),
                    1 - prob_of_non_infect_after_recovery_list[0])
                    ]
        for infect_count in range(1, 4):
            new_infect_record.append(calculate_infect_again_pop(ALL_INFECT, lambda_list_from_2_to_4[infect_count - 1], infect_ratio_in_secrete_type, 
                current_record[3][infect_count - 1],matrix_of_light_severe).copy())
            new_recover_record.append(multiple_list_for_np_list(
                    calculate_infect_again_pop(all_infect_sum(I_2_recovery_record[2]), lambda_list_from_2_to_4[infect_count - 1], infect_ratio_in_secrete_type, 
                        I_2_recovery_record[3][infect_count - 1],matrix_of_light_severe),
                        1 - prob_of_non_infect_after_recovery_list[infect_count]
                    ))
            

        data_list.append(current_record.copy() + [[ALL_INFECT]] + 
            [new_infect_record] + [new_recover_record] + [[lambda_1_m, lambda_1_n]+ lambda_list_from_2_to_4])



        S1 = update_list_for_np_list(
            delta_susceptible_for_mother_antibody(-ALL_INFECT, 
                lambda_1_m, infect_ratio_in_secrete_type.copy(), current_record[0]),
            S1)

        S1.append(mother_anti_rate * secrete_type_distribution.copy())

        
        S2 += S1.pop(0).copy()

        S2 += delta_susceptible_with_no_anti(-ALL_INFECT, lambda_1_n, infect_ratio_in_secrete_type.copy(), current_record[1])  
        
        infect_list[0] = update_list_for_np_list(
            one_time_for_1st_infect(ALL_INFECT, lambda_1_m, lambda_1_n, infect_ratio_in_secrete_type.copy(), 
                current_record[0], current_record[1], light_severe_distribution_list.copy()),
            infect_list[0])

        infect_list[0] = update_list_for_np_list(
            one_time_for_1st_infect(-all_infect_sum(N_recovery_record[2]), N_recovery_record[5][0], N_recovery_record[5][1], infect_ratio_in_secrete_type.copy(), 
                N_recovery_record[0], N_recovery_record[1], light_severe_distribution_list.copy()),
            infect_list[0])
        
        recovery_list[0] = update_list_for_np_list(
            multiple_list_for_np_list(
                one_time_for_1st_infect(all_infect_sum(N_recovery_record[2]), N_recovery_record[5][0], N_recovery_record[5][1], infect_ratio_in_secrete_type, 
                    N_recovery_record[0], N_recovery_record[1], light_severe_distribution_list),
                    1 - prob_of_non_infect_after_recovery_list[0]
                ),
            recovery_list[0])

        recovery_list[0] = update_list_for_np_list(
            calculate_recovery_infected_pop(-ALL_INFECT, 
                lambda_list_from_2_to_4[0], infect_ratio_in_secrete_type, current_record[3][0]),
            recovery_list[0])

        out_list[0] = update_list_for_np_list(
            multiple_list_for_np_list(
                one_time_for_1st_infect(all_infect_sum(N_recovery_record[2]), N_recovery_record[5][0], N_recovery_record[5][1],  infect_ratio_in_secrete_type, 
                    N_recovery_record[0], N_recovery_record[1], light_severe_distribution_list),
                    prob_of_non_infect_after_recovery_list[0]
                ),
            out_list[0])

        for infect_count in range(1, 4):
            infect_list[infect_count] = update_list_for_np_list(
                calculate_infect_again_pop(ALL_INFECT, lambda_list_from_2_to_4[infect_count - 1], infect_ratio_in_secrete_type, 
                    current_record[3][infect_count - 1],matrix_of_light_severe),
                infect_list[infect_count])

            infect_list[infect_count] = update_list_for_np_list(
                calculate_infect_again_pop(-all_infect_sum(I_2_recovery_record[2]), I_2_recovery_record[5][infect_count + 1], infect_ratio_in_secrete_type, 
                    I_2_recovery_record[3][infect_count - 1],matrix_of_light_severe),
                infect_list[infect_count])
            
            recovery_list[infect_count] = update_list_for_np_list(
                multiple_list_for_np_list(
                    calculate_infect_again_pop(all_infect_sum(I_2_recovery_record[2]), I_2_recovery_record[5][infect_count + 1], infect_ratio_in_secrete_type, 
                        I_2_recovery_record[3][infect_count - 1],matrix_of_light_severe),
                        1 - prob_of_non_infect_after_recovery_list[infect_count]
                    ),
                recovery_list[infect_count])

            recovery_list[infect_count] = update_list_for_np_list(
                calculate_recovery_infected_pop(-ALL_INFECT, 
                    lambda_list_from_2_to_4[infect_count - 1], infect_ratio_in_secrete_type, current_record[3][infect_count]),
                recovery_list[infect_count])

            out_list[infect_count] = update_list_for_np_list(
                multiple_list_for_np_list(
                    calculate_infect_again_pop(all_infect_sum(I_2_recovery_record[2]), I_2_recovery_record[5][infect_count + 1], infect_ratio_in_secrete_type, 
                        I_2_recovery_record[3][infect_count - 1],matrix_of_light_severe),
                        prob_of_non_infect_after_recovery_list[infect_count]
                    ),
                out_list[infect_count])

    return data_list

def to_text(data_list, filename):
    data_frame_dict = {}
    data_frame_dict['Time'] = [i + 1 for i in range(0, len(data_list))]
    type_list = ['Secrete', 'W-Secrete','N-Secrete']
    infect_count = ['1','2','3','4']
    level_list = ['light','severe']

    for i in range(0, len(type_list)):
        data_frame_dict['Maternal_Antibody_' + type_list[i]] = [sum_susceptible_for_mother_antibody(j[0])[0, i] for j in data_list]

    for i in range(0, len(type_list)):
        data_frame_dict['Susceptible_' + type_list[i]] = [j[1][0, i] for j in data_list]

    for j in range(0, len(infect_count)):
        for k in range(0, len(level_list)):
            for i in range(0, len(type_list)):
                data_frame_dict[
                    'No.' + infect_count[j] + ' Infection_pop' + type_list[i] + '_' + 'Symptom_' + level_list[k]] \
                    = [one_record[2][j][k][0, i] for one_record in data_list]

    for j in range(0, len(infect_count)):
        for k in range(0, len(level_list)):
            for i in range(0, len(type_list)):
                data_frame_dict[
                    'No.' + infect_count[j] + ' Recovery_pop' + type_list[i] + '_' + 'Symptom_' + level_list[k]] \
                    = [one_record[3][j][k][0, i] for one_record in data_list]   

    for j in range(0, len(infect_count)):
        for k in range(0, len(level_list)):
            for i in range(0, len(type_list)):
                data_frame_dict[
                    'No.' + infect_count[j] + ' Immune_pop' + type_list[i] + '_' + 'Symptom_' + level_list[k]] \
                    = [one_record[4][j][k][0, i] for one_record in data_list]
    
    for j in range(0, len(infect_count)):
        for k in range(0, len(level_list)):
            for i in range(0, len(type_list)):
                data_frame_dict[
                    'No.' + infect_count[j] + 'New Infected' + type_list[i] + '_' + 'Symptom_' + level_list[k]] \
                    = [one_record[7][j][k][0, i] for one_record in data_list]
    
    for j in range(0, len(infect_count)):
        for k in range(0, len(level_list)):
            for i in range(0, len(type_list)):
                data_frame_dict[
                    'No.' + infect_count[j] + 'New Recovery' + type_list[i] + '_' + 'Symptom_' + level_list[k]] \
                    = [one_record[8][j][k][0, i] for one_record in data_list]
    infect_rate_list = ['Maternal Antibody', 'Susceptible', '2nd', '3rd', '4th']
    for i in range(len(infect_rate_list)):
        data_frame_dict[
                    '{0}_lambda'.format(infect_rate_list[i])] \
                    = [one_record[9][i] for one_record in data_list]      
    data_frame_dict['All_infected'] = [one_record[6][0] for one_record in data_list]
    df = pd.DataFrame(data_frame_dict)
    df.to_excel(filename, index=False)
    return

