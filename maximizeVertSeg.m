%requires vani_deterministic, isres. Produces parameters, statistics, the generation when parameter set was found..

function [output, statistics, Gm, VertGoodSets] = maximizeVertSeg()

global VertGoodSet; %modify if change parameter number
global Cutoff; 
global Lb;
global Ub;
global PopSize;
global time_steps;
global eps;

VertGoodSet=ones(1,44);
Cutoff = 14;

minutes=600; %1200;
eps = 0.01; %time step to be used for Euler's method, default is 0.01
time_steps = (minutes / eps);

PopSize=100;%modify population size (100-400) - 100
generations = 300;%modify number of generations - 1000
parents = 5;%modify number of parents. Should be about 1/7 of lightPopSize. 15


Lb = [30,27,10,22,0.12, 0.11,0.11,0.15,32,31,34,31,0.2,0.13,0.28,0.12,0.25,0.1,0.16,0.11,0.26,0.12,8.8,8.6,6.1, 0.8,0.6,0.4,10,0.005,0.06,0.006,0.004,0.0006,0.03,0.001,0.05,0.007,0.03,0.002,0.07,160,200,240];
Ub = [60,57,57,59,0.37,0.39,0.4,0.38,63,62,62,65,0.38,0.39,0.4,0.39,0.4,0.36,0.34,0.34,0.4,0.4,12,11.6,12,2,1.8,1.8,18,0.03,0.3,0.029,0.18,0.009,0.28,0.016,0.29,0.03,0.3,0.024,0.3,720,920,720];


[output,statistics, Gm] = isres('vani_deterministic', 'max', [Lb; Ub], PopSize,generations,parents, 0.45, 1);
%change SRESlight4 to serial when running locally

VertGoodSets=VertGoodSet(2:end,:);

end