global constraints_satisfied;
timesteps = 60000;
 BestSets=[];
 %mh1_array = zeros(18,timesteps);
 f=zeros(length(PossibleSets),1);
 constraints_satisfied=0;
for i = 1:length(PossibleSets)
    x=PossibleSets(i,:);
    disp(i)
%     param_set_wt = repmat(x,[1,1,6]);
% param_set_wt=permute(param_set_wt,[1 3 2]); %change?
% param_set_wt(:,2,4)=0; %delta mutant
% param_set_wt(:,3,1)=0; %her1 mutant
% param_set_wt(:,4,3)=0; %her7 mutant
% param_set_wt(:,5,2)=0; %her6 mutant
% param_set_wt(:,6,3)=0;param_set_wt(:,6,2)=0;%her6 and her7 mutant
% %wperiod=0;wamplitude=0;
%parfor k=1:testPop % FOR TESTING WE COMMENT OUT
 %for k=1:testPop % FOR TESTING WE UNCOMMENT
    %wperiod=0;wamplitude=0;f=0;
    %fprintf("Simulating parameter set %i: \n",k); % Used for creating output file names specific to the paramater
    
    % Booleans to see if we met mutant conditions on each iteration
    %             wt = false; her1_mutant = false; her7_mutant = false; her6_mutant = false;
    %             her76_mutant = false; delta_mutant = false;
    
    %             For the wild type and every mutant, perform the following steps:
    %              1) Adjust job176_Outputthe appropriate protein synthesis rates to create mutants if necessary
    %              2) Run the simulation
    %              3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
    %              4) Otherwise, test oscillation features
    %              5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions
    
    
    %----- Wild Type Conditions --------
    mh1_wt = deterministic_model(x);
    mh1_array(i,:) = mh1_wt(1,:);
   
    if length(mh1_wt)~=timesteps
%     deterministic_model return mh1=0 when the model is
%     inconsistent with biological possibilities
        f(i)=0;
        continue;
    else
        f(i)=f(i)+5;
    end
    if satisfiesWTconditions(mh1_wt)
       % fprintf('satisfiesWTconditions\n');
        f(i)=f(i)+5;
    end
%     if f(i)==5 %% otherwise exit this parameter set?
%         [wperiod, wamplitude]= findPeriodandAmplitude(mh1_wt);
%         f=f+(wperiod > 29 && wperiod < 31 ); %add a point for satisfying period conditions
% %     else
% %         fprintf('~(wperiod > 29 && wperiod < 31 )')
%     end
    
    % --------- Mutant conditions -----------
    
%     % ---- Delta mutant ----
%     
% %     wt_psd = param_set_wt(k,4);
% %     param_set_wt(k,4) = 0; %mutant condition : setting psd=0
%     mh1_delta = deterministic_model(param_set_wt(k,2,:));
%     if length(mh1_delta)==timesteps
%         % this means the model works for the mutant conditions
%         [dperiod, damplitude]= findPeriodandAmplitude(mh1_delta);
%         
%         %amplitude condition to be satisfied for the mutant
%         delta_mutant_amplitude = (damplitude/wamplitude > 0.3) && (damplitude/wamplitude < 0.85);
%         %period condition to be satisfied
%         delta_mutant_period = ((dperiod / wperiod) > 1.07) && ((dperiod / wperiod) < 1.20);
%         %add a point for each mutant condition that is satisfied: period and amplitude
%         f=f+delta_mutant_period + delta_mutant_amplitude;
% %         if ~delta_mutant_period
% %             fprintf('\n~delta_mutant_period:');
% %             fprintf(delta_mutant_period);
% %         end
% %         if ~delta_mutant_amplitude
% %             fprintf('\n~delta_mutant_amplitude:');
% %             fprintf(delta_mutant_amplitude);
% %         end
%     end
% %     param_set_wt(k,4) = wt_psd; %returning parameter set to WT conditions
%     
%     % ---- Her1 mutant ----
%     
% %     wt_psh1 = param_set_wt(k,1);
% %     param_set_wt(k,1) = 0; %mutant condition : setting psh1=0
%     mh1_h1 = deterministic_model(param_set_wt(k,3,:));
%     if length(mh1_h1)==timesteps
%         % this means the model works for the mutant conditions
%         [h1period, h1amplitude]= findPeriodandAmplitude(mh1_h1);
%         
%         %amplitude condition to be satisfied for the mutant
%         her1_mutant_amplitude = ((h1amplitude/wamplitude) > 0.85) && ((h1amplitude/wamplitude) < 1.15);
%         %period condition to be satisfied
%         her1_mutant_period = ((h1period/wperiod) > 0.97) && ((h1period/wperiod) < 1.03);
%         %add a point for each mutant condition that is satisfied: period and amplitude
%         f=f+her1_mutant_period + her1_mutant_amplitude;
% %         if ~her1_mutant_amplitude
% %             fprintf('\n~her1_mutant_amplitude:');
% %             fprintf(her1_mutant_amplitude);
% %         end
% %         if ~her1_mutant_period
% %             fprintf('\n~her1_mutant_period:')
% %             fprintf(her1_mutant_period);
% %         end
%     end
% %     param_set_wt(k,1) = wt_psh1; %returning parameter set to WT conditions
% %     
%     % ---- Her7 mutant ----
%     
% %     wt_psh7 = param_set_wt(k,3);
% %     param_set_wt(k,3) = 0; %mutant condition : setting psh7=0
%     mh1_h7 = deterministic_model(param_set_wt(k,4,:));
%     if length(mh1_h7)==timesteps
%         [h7period, h7amplitude]= findPeriodandAmplitude(mh1_h7);
%         
%         her7_mutant_amplitude = (h7amplitude/wamplitude > 0.1) && (h7amplitude/wamplitude < 0.4);
%         her7_mutant_period = (h7period/wperiod > 0.97) && ((h7period/wperiod) < 1.03);
%         %add a point for each mutant condition that is satisfied: period and amplitude
%         f=f+her7_mutant_period + her7_mutant_amplitude;
% %         if ~her7_mutant_amplitude
% %             fprintf('\n~her7_mutant_amplitude:');
% %             fprintf(her7_mutant_amplitude);
% %         end
% %         if ~her7_mutant_period
% %             fprintf('\n~her7_mutant_period:')
% %             fprintf(her7_mutant_period);
% %         end
%     end
% %     param_set_wt(k,3) = wt_psh7; %returning parameter set to WT conditions
%     
%     % ---- Her6 mutant ----
%     
% %     wt_psh6 = param_set_wt(k,2);
% %     param_set_wt(k,2) = 0; %mutant condition : setting psh7=0
%     mh1_h6 = deterministic_model(param_set_wt(k,5,:));
%     
%     if length(mh1_h6)==timesteps
%         [h6period, h6amplitude]= findPeriodandAmplitude(mh1_h6);
%         
%         her6_mutant_amplitude = (h6amplitude/wamplitude > 0.85) && (h6amplitude/wamplitude < 1.15);
%         her6_mutant_period = ((h6period / wperiod) > 1.05) && ((h6period / wperiod) < 1.07);
%         %add a point for each mutant condition that is satisfied: period and amplitude
%         f=f+her6_mutant_period + her6_mutant_amplitude;
% %         if ~her6_mutant_amplitude
% %             fprintf('\n~her6_mutant_amplitude:');
% %             fprintf(her6_mutant_amplitude);
% %         end
% %         if ~her6_mutant_period
% %             fprintf('\n~her6_mutant_period:')
% %             fprintf(her6_mutant_period);
% %         end
%     end
% %     param_set_wt(k,2) = wt_psh6; %returning parameter set to WT conditions
%     
%     % ---- Her6 and Her7 mutant ----
%     
% %     wt_psh7 = param_set_wt(k,3);
% %     wt_psh6 = param_set_wt(k,2);
% %     param_set_wt(k,3) = 0;
% %     param_set_wt(k,2) = 0;
%     mh1_h76 = deterministic_model(param_set_wt(k,6,:));
%     
%     if length(mh1_h76)==timesteps
%     [h76period, h76amplitude]= findPeriodandAmplitude(mh1_h76);
%     
%     her76_mutant_amplitude = (h76amplitude/wamplitude > 0.85) && (h76amplitude/wamplitude < 1.15);
%     her76_mutant_period = ((h76period / wperiod) > 1.05) && ((h76period / wperiod) < 1.07);
%     %add a point for each mutant condition that is satisfied: period and amplitude
%     f=f+her76_mutant_period + her76_mutant_amplitude;
% 
%     end
    %returning parameter set to WT conditions - don't need to do for the
    %last iteration
%     param_set_wt(k,3) = wt_psh7;
%     param_set_wt(k,2) = wt_psh6;
%     

    
    %  If the paramater set created oscillatory behavior in wild type and all the mutant conditions were satisfied:
    %  1) Print the appropriate message
    %  2) Print the oscillation features into the appropriate files
    %  2) Print the parameter set into the output file containing parameter that passed the conditions.
    if f(i)>=5 %%adding the paramter sets that exceed the cutoff fitness score
        fprintf('The parameter set exceeded the cutoff: f = %i',f);
        BestSets = [BestSets; x]; %#ok<*AGROW>
    end
    
end

%Print the parameter set into the output file containing parameter that passed the conditions.


    

