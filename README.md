# Function Descriptions

SIRO_model code file provides a model for predicting the outcome of epidemic pattern under specific condition.

---

### 'multiple_list_for_np_list(delta_output, prob)'

**Purpose**:  
Multiplies each element in a numpy array ('delta_output') by a constant factor ('prob').  

**Inputs**:  
- 'delta_output' (list of numpy arrays): The input list whose elements are modified in place.  
- 'prob' (float): The multiplier applied to each element in the list.  

**Returns**:  
The modified 'delta_output' list after multiplication.  

---

### 'update_list_for_np_list(delta_output, current_infect_POP)'

**Purpose**:  
Adds corresponding elements from one numpy list ('delta_output') to another ('current_infect_POP'), modifying 'current_infect_POP' in place.  

**Inputs**:  
- 'delta_output' (list of numpy arrays): The list to be added.  
- 'current_infect_POP' (list of numpy arrays): The list receiving the added values.  

**Returns**:  
The updated 'current_infect_POP'.  

---

### 'whether_T_more_than_recovery_period(recorder_list, duration_of_disease)'
**Purpose**:  
Checks if the number of entries in the 'recorder_list' exceeds or equals the 'duration_of_disease', which helps determine if the infection duration has ended.  

**Inputs**:  
- 'recorder_list' (list): The record of infected individuals over time.  
- 'duration_of_disease' (int): The threshold length of the recovery period.  

**Returns**:  
- 'True' if the list length is greater than or equal to the duration; otherwise, 'False'.  

---

### 'all_infect_sum(infectious_population_list)'
**Purpose**:  
Calculates the total number of infected individuals across all subgroups within the 'infectious_population_list'.  

**Inputs**:  
- 'infectious_population_list' (list of numpy arrays): A list containing subgroups of infected populations.  

**Returns**:  
- 'ALL_INFECT_POPULATION' (float): The total number of infected individuals.  

---

### 'parameter_assemble_for_recording(S1, S2, infectious_population_list, recovery_population_list, out_list, lambda_list)'
**Purpose**:  
Aggregates population data and parameters into a structured record for later use.  

**Inputs**:  
- 'S1', 'S2' (numpy arrays): Populations with maternal antibodies and susceptible populations without antibodies, respectively.  
- 'infectious_population_list' (list of numpy arrays): Current infected populations across subgroups.  
- 'recovery_population_list' (list of numpy arrays): Recovered populations grouped by infection rounds.  
- 'out_list' (list of numpy arrays): Immune populations across subgroups.  
- 'lambda_list' (list): Infection rate probabilities.  

**Returns**:  
- A structured list of copied parameters for recording.  

---

### 'zero_generation_for_recording()'
**Purpose**:  
Initializes an empty data record for population groups with zero values.  

**Inputs**:  
None.  

**Returns**:  
- A pre-defined list containing zeroed arrays for all population categories.  

---

### 'recorder(recorder_list, one_record, duration_of_disease)'
**Purpose**:  
Maintains a rolling record of population changes, adding a new entry ('one_record') while removing the oldest entry if the list exceeds the 'duration_of_disease'.  

**Inputs**:  
- 'recorder_list' (list): The current list of records.  
- 'one_record' (list): The new record to be added.  
- 'duration_of_disease' (int): The maximum allowable size for the 'recorder_list'.  

**Returns**:  
- 'recorder_list' (list): Updated list of records.  
- 'temp' (list): The oldest record, removed from the list if necessary.  

---

### 'delta_susceptible_for_mother_antibody(ALL_INFECT, lambda_1_m, infect_ratio_in_secrete_type, S1)'
**Purpose**:  
Calculates the change in population with maternal antibodies due to infections, constrained by the limits of existing population values.  

**Inputs**:  
- 'ALL_INFECT' (float): Total number of currently infected individuals.  
- 'lambda_1_m' (float): Infection rate for populations with maternal antibodies.  
- 'infect_ratio_in_secrete_type' (numpy array): Relative infection rates by secretion type.  
- 'S1' (numpy array): Current population with maternal antibodies.  

**Returns**:  
- 'temp' (numpy array): Adjusted changes in the population with maternal antibodies.  

---

### 'delta_susceptible_with_no_anti(ALL_INFECT, lambda_1_n, infect_ratio_in_secrete_type, S2)'
**Purpose**:  
Calculates the change in susceptible populations without maternal antibodies.  

**Inputs**:  
Same as above, but specific to 'lambda_1_n' and 'S2'.  

**Returns**:  
Similar adjustments for 'S2' populations.  

---

### 'sum_susceptible_for_mother_antibody(S1)'
**Purpose**:  
Aggregates the total population with maternal antibodies across all subgroups.  

**Inputs**:  
- 'S1' (list of numpy arrays): Populations with maternal antibodies.  

**Returns**:  
- 'SUM_S' (numpy array): Total population with maternal antibodies.  

---

### 'run_one_model(input_overall_time, input_lambda, n0, n1, n2, n3, season_rate)'
**Purpose**:  
Simulates the dynamics of infections, recoveries, and immune populations over a given timeframe.  

**Inputs**:  
- 'input_overall_time' (int): Total simulation time.  
- 'input_lambda' (float): Base infection rate.  
- 'n0', 'n1', 'n2', 'n3' (float): Scaling factors for infection rates.  
- 'season_rate' (float): Seasonal variation scaling factor.  

**Returns**:  
- 'data_list' (list): A detailed record of population dynamics over time.  

---

### 'to_text(data_list, filename)'
**Purpose**:  
Organizes simulation results into a structured Excel file, with columns representing various population categories over time.  

**Inputs**:  
- 'data_list' (list): Simulation data generated by 'run_one_model'.  
- 'filename' (str): The output filename for the Excel file.  

**Returns**:  
None (saves data to a file).  

---

This document explains the functionality of each part of the script to help understand its use in modeling population dynamics.
