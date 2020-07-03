function f = findScoreVertSeg(x)

global time_steps;
global eps;
eps=0.01;
% global Cutoff;
% global VertGoodSet;
% global Lb;
% global Ub;
% global PopSize;length(mh1_h76)~=time_steps

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
param_set_wt = repmat(x',[1,1,6]);
param_set_wt=permute(param_set_wt,[1 3 2]); 
param_set_wt(:,2,4)=0; % delta mutant
param_set_wt(:,3,1)=0; % her1 mutant
param_set_wt(:,4,3)=0; % her7 mutant
param_set_wt(:,5,2)=0; % her6 mutant
param_set_wt(:,6,3)=0; param_set_wt(:,6,2)=0; % her6 and her7 mutant
wperiod=0;wamplitude=0;

%parfor k=1:testPop % FOR TESTING WE COMMENT OUT

    % Booleans to see if we met mutant conditions on each iteration
    %             wt = false; her1_mutant = false; her7_mutant = false; her6_mutant = false;
    %             her76_mutant = false; delta_mutant = false;
    
    %             For the wild type and every mutant, perform the following steps:
    %              1) Adjust the appropriate protein synthesis rates to create mutants if necessary
    %              2) Run the simulation
    %              3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
    %              4) Otherwise, test oscillation features
    %              5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions
    
    
    %----- Wild Type Conditions --------
    mh1_wt=deterministic_model(param_set_wt(1,1,:));
%     If length(mh1_wt)~=time_steps 
%     % Deterministic_model return mh1=0 when the model is
%     % inconsistent with biological possibilities
%         f=0;
%         return;
%     end
    if satisfiesWTconditions(mh1_wt)
        f=f+5;
    end
    if f==5 % Otherwise exit this parameter set: Yes set the score for next one to 6.
        [wperiod, wamplitude]= findPeriodandAmplitude(mh1_wt);
        f=f+(wperiod > 29 && wperiod < 31); % Add a point for satisfying period conditions
    end
    
    % --------- Mutant conditions -----------
    
    % ---- Delta mutant ----
    
    if f>=6 % Wildtype conditions are satisfied
        mh1_delta = deterministic_model(param_set_wt(1,2,:));
        % This means the model works for the mutant conditions
        [dperiod, damplitude]= findPeriodandAmplitude(mh1_delta);
        % Amplitude condition to be satisfied for the mutant
        delta_mutant_amplitude = (damplitude/wamplitude > 0.3) && (damplitude/wamplitude < 0.85);
        % Period condition to be satisfied
        delta_mutant_period = ((dperiod / wperiod) > 1.07) && ((dperiod / wperiod) < 1.20);
        % Add a point for each mutant condition that is satisfied: period and amplitude
        f=f+delta_mutant_period + delta_mutant_amplitude;
    end
        
    % ---- Her1 mutant ----
    
    if f>=6 % Wildtype conditions are satisfied
        mh1_h1 = deterministic_model(param_set_wt(1,3,:));
        % This means the model works for the mutant conditions
        [h1period, h1amplitude]= findPeriodandAmplitude(mh1_h1);        
        %amplitude condition to be satisfied for the mutant
        her1_mutant_amplitude = ((h1amplitude/wamplitude) > 0.85) && ((h1amplitude/wamplitude) < 1.15);
        %period condition to be satisfied
        her1_mutant_period = ((h1period/wperiod) > 0.97) && ((h1period/wperiod) < 1.03);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f=f+her1_mutant_period + her1_mutant_amplitude;
    end
    
    % ---- Her7 mutant ----
    if f>=6 % Wildtype conditions are satisfied
         mh1_h7 = deterministic_model(param_set_wt(1,4,:));
        [h7period, h7amplitude]= findPeriodandAmplitude(mh1_h7);        
        her7_mutant_amplitude = (h7amplitude/wamplitude > 0.1) && (h7amplitude/wamplitude < 0.85); % Changed the amplitude condition
        her7_mutant_period = (h7period/wperiod > 0.97) && ((h7period/wperiod) < 1.03);
        % Add a point for each mutant condition that is satisfied: period and amplitude
        f=f+her7_mutant_period + her7_mutant_amplitude;    
    end
    
    % ---- Her6 mutant ----
        
    if f>=6 % Wildtype conditions are satisfied
         mh1_h6 = deterministic_model(param_set_wt(1,5,:));
        [h6period, h6amplitude]= findPeriodandAmplitude(mh1_h6);
        
        her6_mutant_amplitude = (h6amplitude/wamplitude > 0.85) && (h6amplitude/wamplitude < 1.15);
        her6_mutant_period = ((h6period / wperiod) > 1.05) && ((h6period / wperiod) < 1.07);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f=f+her6_mutant_period + her6_mutant_amplitude;
    end
    
    % ---- Her6 and Her7 mutant ----
      
    if f>=6
     mh1_h76 = deterministic_model(param_set_wt(1,6,:));
    [h76period, h76amplitude]= findPeriodandAmplitude(mh1_h76);
    
    her76_mutant_amplitude = (h76amplitude/wamplitude > 0.85) && (h76amplitude/wamplitude < 1.15);
    her76_mutant_period = ((h76period / wperiod) > 1.05) && ((h76period / wperiod) < 1.07);
    % Add a point for each mutant condition that is satisfied: period and amplitude
    f=f+her76_mutant_period + her76_mutant_amplitude;
    end
%     
%     % What failed? delta?
%     if ~delta_mutant_amplitude
%         disp('damplitude/wamplitude');
%         disp(damplitude/wamplitude);
%     end
%     if ~delta_mutant_period
%         disp('dperiod / wperiod:')
%         disp(dperiod / wperiod);
%     end
%     
%     % What failed? her1?
%     if ~her1_mutant_amplitude
%         disp('h1amplitude/wamplitude');
%         disp(h1amplitude/wamplitude);
%     end
%     if ~her1_mutant_period
%         disp('1period / wperiod:')
%         disp(h1period / wperiod);
%     end
%     
%     % What failed? her7?
%     if ~her7_mutant_amplitude
%         disp('h7amplitude/wamplitude');
%         disp(h7amplitude/wamplitude);
%     end
%     if ~her7_mutant_period
%         disp('h7period / wperiod:')
%         disp(h7period / wperiod);
%     end
%     
%     % What failed? her6?
%     if ~her6_mutant_amplitude
%         disp('h6amplitude/wamplitude');
%         disp(h6amplitude/wamplitude);
%     end
%     if ~her6_mutant_period
%         disp('h6period / wperiod:')
%         disp(h6period / wperiod);
%     end
%     
%     % What failed? her76?
%     if ~her76_mutant_amplitude
%         disp('h76amplitude/wamplitude');
%         disp(h76amplitude/wamplitude);
%     end
%     if ~her76_mutant_period
%         disp('h76period / wperiod:')
%         disp(h76period / wperiod);
%     end
%     
    %  If the paramater set created oscillatory behavior in wild type and all the mutant conditions were satisfied:
    %  1) Print the appropriate message
    %  2) Print the oscillation features into the appropriate files
    %  2) Print the parameter set into the output file containing parameter that passed the conditions.
    
end
