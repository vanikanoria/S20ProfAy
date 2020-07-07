function f = find_stochastic_v3_score(x)

% global time_steps;
% global eps;
% eps=0.01;
% global Cutoff;

f=0;
%original parameter set from x
param_set_wt = repmat(x,[1,1,6]);
param_set_wt=permute(param_set_wt,[1 3 2]); 

param_set_wt(:,2,4)=0; %delta mutant
param_set_wt(:,3,1)=0; %her1 mutant
param_set_wt(:,4,3)=0; %her7 mutant
param_set_wt(:,5,2)=0; %her6 mutant
param_set_wt(:,6,3)=0;param_set_wt(:,6,2)=0;%her6 and her7 mutant

% For the wild type and every mutant, perform the following steps:
%  1) Adjust the appropriate protein synthesis rates to create mutants if necessary
%  2) Run the simulation
%  3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
%  4) Otherwise, test oscillation features
%  5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions
    
    
    %----- Wild Type Conditions --------
    mh1_wt=vani_stochastic_v3(param_set_wt(k,1,:));
    if length(mh1_wt)~=timesteps
    %vani_stochastic_v3 return mh1=0 when the model is
    %inconsistent with biological possibilities
        f(k)=0;
        return;
    end
    if satisfiesWTconditions(mh1_wt)
        f(k)=f(k)+5;
    end
    if f(k)==5 %% otherwise exit this parameter set?
        [wperiod, wamplitude]= findPeriodandAmplitude(mh1_wt);
        f(k)=f(k)+(wperiod > 29 && wperiod < 31 ); %add a point for satisfying period conditions
    end
    
    % --------- Mutant conditions -----------
    
    % ---- Delta mutant ----
    mh1_delta = vani_stochastic_v3(param_set_wt(k,2,:));
    if length(mh1_delta)==timesteps
        % this means the model works for the mutant conditions
        [dperiod, damplitude]= findPeriodandAmplitude(mh1_delta);
        
        %amplitude condition to be satisfied for the mutant
        delta_mutant_amplitude = (damplitude/wamplitude > 0.3) && (damplitude/wamplitude < 0.85);
        %period condition to be satisfied
        delta_mutant_period = ((dperiod / wperiod) > 1.07) && ((dperiod / wperiod) < 1.20);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f(k)=f(k)+delta_mutant_period + delta_mutant_amplitude;
    end
   
    % ---- Her1 mutant ----
    mh1_h1 = vani_stochastic_v3(param_set_wt(k,3,:));
    if length(mh1_h1)==timesteps
        % this means the model works for the mutant conditions
        [h1period, h1amplitude]= findPeriodandAmplitude(mh1_h1);
        
        %amplitude condition to be satisfied for the mutant
        her1_mutant_amplitude = ((h1amplitude/wamplitude) > 0.85) && ((h1amplitude/wamplitude) < 1.15);
        %period condition to be satisfied
        her1_mutant_period = ((h1period/wperiod) > 0.97) && ((h1period/wperiod) < 1.03);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f(k)=f(k)+her1_mutant_period + her1_mutant_amplitude;
    end
     
    % ---- Her7 mutant ----
    
    mh1_h7 = vani_stochastic_v3(param_set_wt(k,4,:));
    if length(mh1_h7)==timesteps
        [h7period, h7amplitude]= findPeriodandAmplitude(mh1_h7);
        
        her7_mutant_amplitude = (h7amplitude/wamplitude > 0.1) && (h7amplitude/wamplitude < 0.4);
        her7_mutant_period = (h7period/wperiod > 0.97) && ((h7period/wperiod) < 1.03);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f(k)=f(k)+her7_mutant_period + her7_mutant_amplitude;
    end
    
    % ---- Her6 mutant ----
    
    mh1_h6 = vani_stochastic_v3(param_set_wt(k,5,:));
    
    if length(mh1_h6)==timesteps
        [h6period, h6amplitude]= findPeriodandAmplitude(mh1_h6);
        
        her6_mutant_amplitude = (h6amplitude/wamplitude > 0.85) && (h6amplitude/wamplitude < 1.15);
        her6_mutant_period = ((h6period / wperiod) > 1.05) && ((h6period / wperiod) < 1.07);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f(k)=f(k)+her6_mutant_period + her6_mutant_amplitude;
    end

    % ---- Her6 and Her7 mutant ----
    
    mh1_h76 = vani_stochastic_v3(param_set_wt(k,6,:));
    
    if length(mh1_h76)==timesteps
    [h76period, h76amplitude]= findPeriodandAmplitude(mh1_h76);
    
    her76_mutant_amplitude = (h76amplitude/wamplitude > 0.85) && (h76amplitude/wamplitude < 1.15);
    her76_mutant_period = ((h76period / wperiod) > 1.05) && ((h76period / wperiod) < 1.07);
    %add a point for each mutant condition that is satisfied: period and amplitude
    f(k)=f(k)+her76_mutant_period + her76_mutant_amplitude;
    end
end
