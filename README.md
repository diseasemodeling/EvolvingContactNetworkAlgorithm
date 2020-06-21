# Agent-based evolving network modeling (ABENM) with Evolving Contact Network Algorithm (ECNA)
- Demonstration in Netlogo 

## 
Demonstration of agent-based evolving network modeling (ABENM) - a new simulation technique that combines features of agent-based modeling with compartmental modeling
It uses a new evolving contact network algorithm (ECNA) - A new algorithm for generation of scale-free networks that follow power-law degree distribution .

It is developed for simulation of diseases with low prevalence, where using current agent-based network modelingcan be computationaly infeasible, and where contact network structures are relevant to model that compartmnetal modeling is not suitable. 



## HOW TO USE IT
A basic understanidng of scale-free networks would be necessary to understand the concepts of this work. 

To run the model, on the interface, click on the following button, in order 
1. clear
2. setup
3. runECNA

## THINGS TO NOTICE

ABENM simulates only infected persons and their immediate contacts (infected and suceptible) as agents. All other susceptibles are simulated as compartmental modeling, with each compartment representing the degree (log-2 bin) of the nodes. 
As new persons becomes infected, the ECNA determines who the immediate contacts are- i.e., what is the their degree, and transitions correpsonding susceptible persons from the compartment to the network as neighbors of the newly infected node. 


## CREDITS AND REFERENCES

Developed by the disease modeling lab for https://blogs.umass.edu/chaitrag/projects/ 
Funding: NIH R01AI127236
