# Co-optimization-of-Distributed-Energy-Resources-Investment-and-Network-Expansion-Planning
Material Repository of the paper "Co-optimization of Distributed Energy Resources Investment and Network Expansion Planning" 
Inside you will find:
1) Code in matlab of our model to run it for veryfication reasons
2) Network setup used to our testing
3) Detailed Results Matrices and produced figures used in the paper

The proposed running process of our model is the following:
1) Run locally in you Personal Computer "single_level_sizing_v2.m"
   
   a) Choose the variable values, such RoI, Investement Yearlife, CVaR, assets size blocks etc. according to the needs of your own co-optimization problem or test.
   b) Select between the proposed Transmission Networks and Distribution Networks. Network Models are stored in folder "Network Setup". Model provides, also, the ability to users to store and run their own TN or DN Model.
   c) Choose the Number of Scenarios and Number of Timeslots you prefer according to the needs of your needs of your own co-optimization problem or test.
   d) Select between minimum_return_choice = 0 or 1 for Co-Optimization model or DSO-oriented model, respectively.
   
2) Export the results of your testing in the files investement_report and power_flow_report
