function f = findScoreVertSeg(x)

global time_steps;
global eps;
eps=0.01;
% global Cutoff;
% global VertGoodSet;
% global Lb;
% global Ub;
% global PopSize;

% Iterate through every parameter set
% CHUNK_SIZE = 10000;
% PARS = 1; % number of parameters to simulate, default value is 1
% x = 2; y = 1; %width and height of the tissue being simulated, default is 2x1
%toPrint = false;  % boolean marking whether or not concentrations should be printed to a text file

%mutants = ["/wt", "/delta", "/her6", "/her1", "/her7", "/her76"];

f=0;
% testPop=1; % FOR TESTING WE UNCOMMENT
%num_of_parameters = 44;

%%%%%%%%%%% NEED TO INITIALIZE GOODSET %%%%%%
%Goodset = 

%creating a 2-dimensional array to store parameter sets of each population
%param_sets = zeros(testPop,num_of_parameters);

%original parameter set from x
param_set_wt = x;
   
% For the wild type and every mutant, perform the following steps:
%  1) Adjust the appropriate protein synthesis rates to create mutants if necessary
%  2) Run the simulation
%  3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
%  4) Otherwise, test oscillation features
%  5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions
    
    
%  here we use the same parameter set param_set_wt and modify it for each
%  condition and after testing for that condition, return it to the
%  original condition. This is to save space.

    %----- Wild Type Conditions --------
    mh1_wt=deterministic_model(param_set_wt);
    if length(mh1_wt)~=time_steps 
    %deterministic_model return mh1=0 when the model is
    %inconsistent with biological possibilities
        f=0;
        return;
    end
    if satisfiesWTconditions(mh1_wt(1,:))
        f=f+5;
    end
    if f==5 %% otherwise exit this parameter set?
        [wperiod, wamplitude]= findPeriodandAmplitude(mh1_wt(1,:));
        f=f+(wperiod > 29 && wperiod < 31 ); %add a point for satisfying period conditions
    end
    
    % --------- Mutant conditions -----------
    
    % ---- Delta mutant ----
    
    wt_psd = param_set_wt(4);
    param_set_wt(4) = 0; %mutant condition : setting psd=0
    mh1_delta = deterministic_model(param_set_wt);
    if length(mh1_delta)==time_steps
        % this means the model works for the mutant conditions
        [dperiod, damplitude]= findPeriodandAmplitude(mh1_delta(1,:));
        
        %amplitude condition to be satisfied for the mutant
        delta_mutant_amplitude = (damplitude/wamplitude > 0.3) && (damplitude/wamplitude < 0.85);
        %period condition to be satisfied
        delta_mutant_period = ((dperiod / wperiod) > 1.07) && ((dperiod / wperiod) < 1.20);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f=f+delta_mutant_period + delta_mutant_amplitude;
    end
    param_set_wt(4) = wt_psd; %returning parameter set to WT conditions
    
    % ---- Her1 mutant ----
    
    wt_psh1 = param_set_wt(1);
    param_set_wt(1) = 0; %mutant condition : setting psh1=0
    mh1_h1 = deterministic_model(param_set_wt);
    if length(mh1_h1)==time_steps
        % this means the model works for the mutant conditions
        [h1period, h1amplitude]= findPeriodandAmplitude(mh1_h1(1,:));
        
        %amplitude condition to be satisfied for the mutant
        her1_mutant_amplitude = ((h1amplitude/wamplitude) > 0.85) && ((h1amplitude/wamplitude) < 1.15);
        %period condition to be satisfied
        her1_mutant_period = ((h1period/wperiod) > 0.97) && ((h1period/wperiod) < 1.03);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f=f+her1_mutant_period + her1_mutant_amplitude;
    end
    param_set_wt(1) = wt_psh1; %returning parameter set to WT conditions
    
    % ---- Her7 mutant ----
    
    wt_psh7 = param_set_wt(3);
    param_set_wt(3) = 0; %mutant condition : setting psh7=0
    mh1_h7 = deterministic_model(param_set_wt);
    if length(mh1_h7)==time_steps
        [h7period, h7amplitude]= findPeriodandAmplitude(mh1_h7(1,:));
        
        her7_mutant_amplitude = (h7amplitude/wamplitude > 0.1) && (h7amplitude/wamplitude < 0.4);
        her7_mutant_period = (h7period/wperiod > 0.97) && ((h7period/wperiod) < 1.03);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f=f+her7_mutant_period + her7_mutant_amplitude;  
    else
        disp('length(mh1_h7)~=time_steps')
    end
    param_set_wt(3) = wt_psh7; %returning parameter set to WT conditions
    
    % ---- Her6 mutant ----
    
    wt_psh6 = param_set_wt(2);
    param_set_wt(2) = 0; %mutant condition : setting psh7=0
    mh1_h6 = deterministic_model(param_set_wt);
    
    if length(mh1_h6)==time_steps
        [h6period, h6amplitude]= findPeriodandAmplitude(mh1_h6(1,:));
        
        her6_mutant_amplitude = (h6amplitude/wamplitude > 0.85) && (h6amplitude/wamplitude < 1.15);
        her6_mutant_period = ((h6period / wperiod) > 1.05) && ((h6period / wperiod) < 1.07);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f=f+her6_mutant_period + her6_mutant_amplitude;
        if ~her6_mutant_amplitude
            disp('~her6_mutant_amplitude:');
            disp(h6_mutant_amplitude);
        end
        if ~her6_mutant_period
            disp('~her6_mutant_period:')
            disp(h6period);
        end
        else
        disp('length(mh1_h6)~=time_steps')
    end
    param_set_wt(2) = wt_psh6; %returning parameter set to WT conditions
    
    % ---- Her6 and Her7 mutant ----
    
    param_set_wt(3) = 0;
    param_set_wt(2) = 0;
    mh1_h76 = deterministic_model(param_set_wt);
    
    if length(mh1_h76)==time_steps
    [h76period, h76amplitude]= findPeriodandAmplitude(mh1_h76(1,:));
    
    her76_mutant_amplitude = (h76amplitude/wamplitude > 0.85) && (h76amplitude/wamplitude < 1.15);
    her76_mutant_period = ((h76period / wperiod) > 1.05) && ((h76period / wperiod) < 1.07);
    %add a point for each mutant condition that is satisfied: period and amplitude
    f=f+her76_mutant_period + her76_mutant_amplitude;
    end

    
    %  If the paramater set created oscillatory behavior in wild type and all the mutant conditions were satisfied:
    %  1) Print the appropriate message
    %  2) Print the oscillation features into the appropriate files
    %  2) Print the parameter set into the output file containing parameter that passed the conditions.
    
end