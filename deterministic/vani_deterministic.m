% deterministic simulator for the zebrafish segmentation clock
% MATLAB version of the C code by Ahmet Ay, Jack Holland, and Adriana Sperlea
 
%  The program simulates the behavior of the zebrafish segmentation clock for two cells,
%  through a mathematical model using delay differential equations.
%  To solve the equations, the program uses Euler's method, with an adjustable stepsize.
%  The program tries to find parameter sets which replicate the behavior of the system in wild type
%  and several mutants.

function[f,g] = vani_deterministic(x)
%size of x according to isres = (PopSize,number of parameters)
global Cutoff;
global VertGoodSet;
global Lb;
global Ub;
global PopSize;
global time_steps;

% Iterate through every parameter set

testPop=PopSize;
f=zeros(PopSize,1);
timesteps = time_steps;
% testPop=1; % FOR TESTING WE UNCOMMENT
%num_of_parameters = 44;


%creating a 2-dimensional array to store parameter sets of each population
%param_sets = zeros(testPop,num_of_parameters);

%original parameter set from x
param_set_wt = repmat(x,[1,1,6]);
param_set_wt=permute(param_set_wt,[1 3 2]); %change?
param_set_wt(:,2,4)=0; %delta mutant
param_set_wt(:,3,1)=0; %her1 mutant
param_set_wt(:,4,3)=0; %her7 mutant
param_set_wt(:,5,2)=0; %her6 mutant
param_set_wt(:,6,3)=0;param_set_wt(:,6,2)=0;%her6 and her7 mutant
wperiod=0;wamplitude=0;
parfor k=1:testPop % FOR TESTING WE COMMENT OUT
% for k=1:testPop % FOR TESTING WE UNCOMMENT
    wperiod=0;wamplitude=0;
     
    %             For the wild type and every mutant, perform the following steps:
    %              1) Adjust the appropriate protein synthesis rates to create mutants if necessary
    %              2) Run the simulation
    %              3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
    %              4) Otherwise, test oscillation features
    %              5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions
    
    
    %----- Wild Type Conditions --------
    mh1_wt=deterministic_model(param_set_wt(k,1,:));
    if length(mh1_wt)~=timesteps
    %deterministic_model return mh1=0 when the model is
    %inconsistent with biological possibilities
        f(k)=0;
        continue;
    end
    if satisfiesWTconditions(mh1_wt)
       % fprintf('satisfiesWTconditions\n');
        f(k)=f(k)+5;
    end
    if f(k)==5 %% otherwise exit this parameter set?
        [wperiod, wamplitude]= findPeriodandAmplitude(mh1_wt);
        f(k)=f(k)+(wperiod > 29 && wperiod < 31 ); %add a point for satisfying period conditions
    end
    
    % --------- Mutant conditions -----------
    
    % ---- Delta mutant ----
    mh1_delta = deterministic_model(param_set_wt(k,2,:));
    if length(mh1_delta)==timesteps
        % this means the model wor1s for the mutant conditions
        [dperiod, damplitude]= findPeriodandAmplitude(mh1_delta);
        
        %amplitude condition to be satisfied for the mutant
        delta_mutant_amplitude = (damplitude/wamplitude > 0.3) && (damplitude/wamplitude < 0.85);
        %period condition to be satisfied
        delta_mutant_period = ((dperiod / wperiod) > 1.07) && ((dperiod / wperiod) < 1.20);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f(k)=f(k)+delta_mutant_period + delta_mutant_amplitude;
    end
    
    % ---- Her1 mutant ----
    
    mh1_h1 = deterministic_model(param_set_wt(k,3,:));
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
%     param_set_wt(1,1) = wt_psh1; %returning parameter set to WT conditions
%     
    % ---- Her7 mutant ----
    
    mh1_h7 = deterministic_model(param_set_wt(k,4,:));
    if length(mh1_h7)==timesteps
        [h7period, h7amplitude]= findPeriodandAmplitude(mh1_h7);
        
        her7_mutant_amplitude = (h7amplitude/wamplitude > 0.1) && (h7amplitude/wamplitude < 0.4);
        her7_mutant_period = (h7period/wperiod > 0.97) && ((h7period/wperiod) < 1.03);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f(k)=f(k)+her7_mutant_period + her7_mutant_amplitude;

    end
    
    % ---- Her6 mutant ----
    mh1_h6 = deterministic_model(param_set_wt(k,5,:));
    
    if length(mh1_h6)==timesteps
        [h6period, h6amplitude]= findPeriodandAmplitude(mh1_h6);
        
        her6_mutant_amplitude = (h6amplitude/wamplitude > 0.85) && (h6amplitude/wamplitude < 1.15);
        her6_mutant_period = ((h6period / wperiod) > 1.05) && ((h6period / wperiod) < 1.07);
        %add a point for each mutant condition that is satisfied: period and amplitude
        f(k)=f(k)+her6_mutant_period + her6_mutant_amplitude;

    end
    
    % ---- Her6 and Her7 mutant ----
    mh1_h76 = deterministic_model(param_set_wt(k,6,:));
    
    if length(mh1_h76)==timesteps
    [h76period, h76amplitude]= findPeriodandAmplitude(mh1_h76);
    
    her76_mutant_amplitude = (h76amplitude/wamplitude > 0.85) && (h76amplitude/wamplitude < 1.15);
    her76_mutant_period = ((h76period / wperiod) > 1.05) && ((h76period / wperiod) < 1.07);
    %add a point for each mutant condition that is satisfied: period and amplitude
    f(k)=f(k)+her76_mutant_period + her76_mutant_amplitude;
%     end
    end
      
    %  If the paramater set created oscillatory behavior in wild type and all the mutant conditions were satisfied:
    %  1) Print the appropriate message
    %  2) Print the oscillation features into the appropriate files
    %  2) Print the parameter set into the output file containing parameter that passed the conditions.
    
end

%Print the parameter set into the output file containing parameter that passed the conditions.

for j= 1:testPop
    if f(j)>=Cutoff %%adding the paramter sets that exceed the cutoff fitness score
        fprintf('The parameter set exceeded the cutoff: f(%i) = %f\n',j,f(j));
        VertGoodSet = [VertGoodSet; x(j,:)]; %#o1<*AGROW>
    end
%             if f(j)>= Cutoff %don't need this
%                  dlmwrite('VertSegSets0615.csv',x(j,:),'delimiter',',','-append');
%             end
end

fprintf('  ');
%Constraints. g means something is <= 0
g(:,1)=-x(:,1)+Lb(1);
g(:,2)=x(:,1)-Ub(1);
g(:,3)=-x(:,2)+Lb(2);
g(:,4)=x(:,2)-Ub(2);
g(:,5)=-x(:,3)+Lb(3);
g(:,6)=x(:,3)-Ub(3);
g(:,7)=-x(:,4)+Lb(4);
g(:,8)=x(:,4)-Ub(4);
g(:,9)=-x(:,5)+Lb(5);
g(:,10)=x(:,5)-Ub(5);
g(:,11)=-x(:,6)+Lb(6);
g(:,12)=x(:,6)-Ub(6);
g(:,13)=-x(:,7)+Lb(7);
g(:,14)=x(:,7)-Ub(7);
g(:,15)=-x(:,8)+Lb(8);
g(:,16)=x(:,8)-Ub(8);
g(:,17)=-x(:,9)+Lb(9);
g(:,18)=x(:,9)-Ub(9);
g(:,19)=-x(:,10)+Lb(10);
g(:,20)=x(:,10)-Ub(10);
g(:,21)=-x(:,11)+Lb(11);
g(:,22)=x(:,11)-Ub(11);
g(:,23)=-x(:,12)+Lb(12);
g(:,24)=x(:,12)-Ub(12);
g(:,25)=-x(:,13)+Lb(13);
g(:,26)=x(:,13)-Ub(13);
g(:,27)=-x(:,14)+Lb(14);
g(:,28)=x(:,14)-Ub(14);
g(:,29)=-x(:,15)+Lb(15);
g(:,30)=x(:,15)-Ub(15);
g(:,31)=-x(:,16)+Lb(16);
g(:,32)=x(:,16)-Ub(16);
g(:,33)=-x(:,17)+Lb(17);
g(:,34)=x(:,17)-Ub(17);
g(:,35)=-x(:,18)+Lb(18);
g(:,36)=x(:,18)-Ub(18);
g(:,37)=-x(:,19)+Lb(19);
g(:,38)=x(:,19)-Ub(19);
g(:,39)=-x(:,20)+Lb(20);
g(:,40)=x(:,20)-Ub(20);
g(:,41)=-x(:,21)+Lb(21);
g(:,42)=x(:,21)-Ub(21);
g(:,43)=-x(:,22)+Lb(22);
g(:,44)=x(:,22)-Ub(22);
g(:,45)=-x(:,23)+Lb(23);
g(:,46)=x(:,23)-Ub(23);
g(:,47)=-x(:,24)+Lb(24);
g(:,48)=x(:,24)-Ub(24);
g(:,49)=-x(:,25)+Lb(25);
g(:,50)=x(:,25)-Ub(25);
g(:,51)=-x(:,26)+Lb(26);
g(:,52)=x(:,26)-Ub(26);
g(:,53)=-x(:,27)+Lb(27);
g(:,54)=x(:,27)-Ub(27);
g(:,55)=-x(:,28)+Lb(28);
g(:,56)=x(:,28)-Ub(28);
g(:,57)=-x(:,29)+Lb(29);
g(:,58)=x(:,29)-Ub(29);
g(:,59)=-x(:,30)+Lb(30);
g(:,60)=x(:,30)-Ub(30);
g(:,61)=-x(:,31)+Lb(31);
g(:,62)=x(:,31)-Ub(31);
g(:,63)=-x(:,32)+Lb(32);
g(:,64)=x(:,32)-Ub(32);
g(:,65)=-x(:,33)+Lb(33);
g(:,66)=x(:,33)-Ub(33);
g(:,67)=-x(:,34)+Lb(34);
g(:,68)=x(:,34)-Ub(34);
g(:,69)=-x(:,35)+Lb(35);
g(:,70)=x(:,35)-Ub(35);
g(:,71)=-x(:,36)+Lb(36);
g(:,72)=x(:,36)-Ub(36);
g(:,73)=-x(:,37)+Lb(37);
g(:,74)=x(:,37)-Ub(37);
g(:,75)=-x(:,38)+Lb(38);
g(:,76)=x(:,38)-Ub(38);
g(:,77)=-x(:,39)+Lb(39);
g(:,78)=x(:,39)-Ub(39);
g(:,79)=-x(:,40)+Lb(40);
g(:,80)=x(:,40)-Ub(40);
g(:,81)=-x(:,41)+Lb(41);
g(:,82)=x(:,41)-Ub(41);
g(:,83)=-x(:,42)+Lb(42);
g(:,84)=x(:,42)-Ub(42);
g(:,85)=-x(:,43)+Lb(43);
g(:,86)=x(:,43)-Ub(43);
g(:,87)=-x(:,44)+Lb(44);
g(:,88)=x(:,44)-Ub(44);
end
