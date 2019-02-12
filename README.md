## Chapter5 
This folder contains the scripts behind the calcualtions presented in the article:
# A risk-based approach for calibration of design codes
Michele Baravalle; Jochen Köhler.  
Structural Safety - Volume 78, May 2019, Pages 63-75. 
https://doi.org/10.1016/j.strusafe.2018.12.003  

### Info  
The scripts with names beginning with “MAIN” are the main scripts where all the inputs are defined and the calculations are performed. All the other scripts are called from the main scripts.  
•	MAIN_Optimize_lifetime_relaibility_target.m – Optimizes the lifetime target reliability as the system and at the component level, it evaluates the robustness of the design with the optimal reliability, and it optimizes each structure individually.  
•	MAIN_Optimize_yearly_relaibility_target.m – Optimizes the yearly target reliability as the system and at the component level, it evaluates the robustness of the design with the optimal reliability, and it optimizes each structure individually. 
•	The use of parallel computing is highly recommended for shortening the running time

### The scripts make use of functions in the following packages:  
•	FERUMcore and FERUMsystems_MSR packages available at http://projects.ce.berkeley.edu/ferum/.   
•	Matlab Statistics and Machine Learning Toolbox (see https://se.mathworks.com/products/statistics.html )  
