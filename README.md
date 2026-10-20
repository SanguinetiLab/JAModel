# Computational Joint Action
This repository contains the source code for the computational model of joint action described in paper [De Vicariis C, Chackochan VT, Bandini L, Ravaschio E, Sanguineti V (2024) Computational joint action: Dynamical models to understand the development of joint coordination. PLoS Comput Biol 20(10): e1011948](https://doi.org/10.1371/journal.pcbi.1011948).

The model is summarized below.
![Graphical Model](./graphicalmodel.png)

Each agent has a sensory system, a partner model and a controller. 

-The *sensory system* perceives the partner action.
-The *partner model* combines prior experience and sensory information to predict the partner's next action.
-The *controller* generates the agent's action on the basis of the prediction of the partner model.

Directory ./Data_Paper contains the data used in the paper. 

Directory ./src contains the classes that specify the general model, which can be used for both simulating a joint action scenario or to analyze joint action time series.


