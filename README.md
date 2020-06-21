# Agent-based evolving network modeling (ABENM) with Evolving Contact Network Algorithm (ECNA)

## About
Demonstration of agent-based evolving network modeling (ABENM) - a new simulation technique that combines features of agent-based modeling with compartmental modeling
It uses a new evolving contact network algorithm (ECNA) - A new algorithm for generation of scale-free networks that follow power-law degree distribution .

It is developed for simulation of epidemic projections for infectious diseases with low prevalence, where using current agent-based network modeling can be computationaly infeasible, and where contact network structures are relevant to model that compartmental modeling is not suitable. 

ABENM simulates only infected persons and their immediate contacts (infected and suceptible) as agents. All other susceptibles are simulated as compartmental modeling, with each compartment representing the degree (log-2 bin) of the nodes. 
As new persons become newly infected, the ECNA determines who their immediate contacts are- specifcally, what is the their degree. It then draws a suceptibles persons from the compartmental model corresponding to that speicifc compartment and transitions them to the network as neighbors of the newly infected node. 

## Software
A sample computational model to demonstrate ABENM and ECNA is constructed in the Netlogo software (Netlogo software is a free softwrae- it can be downloaded here- https://ccl.northwestern.edu/netlogo/ )

## HOW TO USE IT
A basic understanding of scale-free networks would be necessary to understand the concepts of this work. 

To run the model. Open the ECNA.nlogo file in  Netlogo 
ON the interface, click on the following buttons, in ordee 
1. clear
2. setup
3. runECNA


## CREDITS AND REFERENCES

Developed by the disease modeling lab for https://blogs.umass.edu/chaitrag/projects/ 
Funding: NIH R01AI127236
