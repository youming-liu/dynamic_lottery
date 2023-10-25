
This readme.txt file was generated on 2023-10-15 by Youming Liu

-------------------
GENERAL INFORMATION
-------------------

Title: Matlab Code Supporting 'Dynamic Efficiency in Resource Allocation: Evidence from Vehicle License Lotteries in Beijing'


Author Information 

1. Name: Youming Liu
Institution: Bank of Canada, Banking and Payment
Email: youmingliu@bank-banque-canada.ca

2. Name: Shanjun Li
Institution: Dyson School of Applied Economics and Management, Cornell University and NBER
Email:SL2448@cornell.edu

3. Caixia Shen
Institution: Zhejiang University of Finance and Economics
Email: shencaixia@gmail.com; 


-----------------------------------
CODE FILES AND DESCRIPTIONS 
-----------------------------------

Main_DemandEstimation: The main demand estimation program 

Demand_ImportData : Imports data files into Matlab 

f_GMMobj: Compute GMM objective function value

f_mu: Compute consumer heterogenious preference on vehicle attributes (not including prices)

f_rho: Compute consumer heterogenious preference on vehicle prices

f_ev: value function of vehicle purchase with the absence of lottery policy

f_ev0: value function of vehicle purchase in the lottery period

f_ew: value function of lottery participation

f_CalcDeltaEvolution: estimate AR(1) process of the logit-inclusive value, delta, in vehicle purchase stage, and predicted values for the next period

f_CalcPhiEvolution: estimate AR(1) process of the inclusive value, phi, in lottery participation stage, and perdicted values for the next period

f_DynamicShareHat: compute market shares

f_DynamicMeanVal: back out mean utilities for lottery participation and vehicle models

f_moments: Compute moment conditions

f_se: Compute standard errors using numerical grandient 

