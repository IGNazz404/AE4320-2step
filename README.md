# AE4320-2step
The code for the aerodynamic model identification using the two-step approach.

Thank you for reading the README file for the code of my AE4320 assignment.

This README serves to guide to professor through the various files used to construct the results
for the 2-step approach assignment.

Main Files
------------
- assignmentMultiStateIteratedEKF.m
	- Performs the data-pre processing and performs the (I)EKF algorithm
	- Dependencies:
		- rk4.m: Numerical integration algorithm, Runge-Kutta
		- assignment_kf_calc_h.m: system of measurement equations
		- assignment_kf_calc_f.m: system of state derivative equations
		- assignment_kf_calc_Fx.m: Jacobian matrix of measurement system
		- assignment_kf_calc_Hx.m: Jacobian matrix of state system

- parameterEstimation.m
	- Performs the aerodynamic model identificiation and validation
	- Main Dependencies:
		- OLS_Identification.m: performs Ordinary Least Squares for the identification
			dataset
		- OLS_Validation.m: performs Ordinary Least Squares for the Validation
			dataset
		- WLS_Identification.m: performs Weighted Least Squares for the identification
			dataset
		- WLS_Validation.m: performs Weighted Least Squares for the Validation
			dataset
	- Other Dependencies:
		- dominance_Tester.m: Constructs boxplots which are used in deriving the most 
			dominant model terms in validation
		- M2structure.m: Contains the alternative model structure used in validation

Files Used for Manual Calculations
-----------------------------------
These are the files used to perform manual calculations, they are not directly used in one
of the main files.
- assignment_calc_Jac_f_syms.m, assignment_calc_Jac_h_syms.m, DHx_calc.m, der_eqs.m:
	A set of files used to find the Jacobian matrix entries for the measurement and state
	systems. These were used before I knew about the 'jacobian' function in Matlab :)

Data Files:
--------------
- Identification Dataset: "da3211_2.mat", "dadoublet_1.mat", "de3211_1.mat", "dr3211_1.mat"
- Validation Dataset: "dedoublet_1.mat", "dr3211_2.mat"

Contact Details
----------------
Name:			Ynias J.E. Prencipe 
Student Number:		4777158
Student e-mail:		y.j.e.prencipe@student.tudelft.nl
