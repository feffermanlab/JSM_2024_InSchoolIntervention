# JSM_2024_InSchoolIntervention
Code for observing model output under different theorized responses to intervention


The single Matlab file in this repository is used to solve the IVP of the model for spread of infectious disease through a school. 
The school is given an explicit relational structure to designate which classrooms share local resources. It is assumed that 
every classroom shares some global resources.

Once you specify an adjacency matrix for the structure of the school and classroom sizes for each classroom, and once you specify parameters for transmission rates and recovery rates, you can begin to explore how changing the control parameters, $\phi, \delta_a, \delta_b, \delta_c$ changes the model output. 

The function odefcn contains the system of ODEs itself, each additional function is used to make an observation in some region of the parameter space. 
