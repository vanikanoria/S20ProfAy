nrun=1;%30;
RHO1=zeros([nrun,1]); PVAL1=zeros([nrun,1]);
RHO2=zeros([nrun,1]); PVAL2=zeros([nrun,1]);
RHO3=zeros([nrun,1]); PVAL3=zeros([nrun,1]);

for run=1:nrun 
maxi=10000000; % Maximum number of iterations to run.
tend=60; %240; % Maximum time for simulation.
format compact 
% rand('state',sum(100*clock)) 
rand('state',run) % Sets random number generator to a specific state.

Time=[0]; % Simulation time  
mh1v=[0]; % Her1 mRNA levels
mh7v=[0]; % Her7 mRNA levels

% % Reaction rates in prev model
% msh1=5; msh7=5; mdh1=0.2; mdh7=0.2; % mRNA synthesis and degradation rates.
% 
% psh1=5; psh7=5; pdh1=0.3; pdh7=0.3; % Protein synthesis and degradation rates.

% ac1h1h1=0.01; dc1h1h1=400; % DNA association and dissocation rates.
% ac2h1h1=0.01; dc2h1h1=400;
% ac1h6h7=0.01; dc1h6h7=400;
% ac2h6h7=0.01; dc2h6h7=400;

% % Parameter values used in the 2013 study (from Table S1)
psh1=49.9139;
psh6=34.3117;
psh7=28.5626;
psd=37.7828;
pdh1=0.34951;
pdh6=0.14824;
pdh7=0.249715;
pdd=0.324316;
msh1=48.3084;
msh6=36.4073;
msh7=39.685;
msd=60.5577;
mdh1=0.322965;
mdh6=0.146372;
mdh7=0.381738;
mdd=0.352056;
pdh11=0.390961;
pdh16=0.29774;
pdh17=0.320157;
pdh66=0.268042;
pdh67=0.352037;
pdh77=0.251601;
nmh1=10.0213;
nmh6=10;
nmh7=10.4515;
nmd=7.74472;
nph1=1.5398;
nph6=0.886233;
nph7=0.539972;
npd=13.2661;
dah1h1=0.0179429;
ddh1h1=0.220856;
dah1h6=0.0270209;
ddh1h6=0.0917567;
dah1h7=0.00120525;
ddh1h7=0.258167;
dah6h6=0.0148271;
ddh6h6=0.251173;
dah7h6=0.0216093;
ddh7h6=0.188923;
dah7h7=0.0202756;
ddh7h7=0.161018;
critph11=587.298;
critph76=769.628;
critpd=490.254;

% Data storage for delayed reaction times
s1 = [inf]; % Time for each delayed reaction is scheduled in these vectors.
s2 = [inf];
s3 = [inf];
s4 = [inf];
s5 = [inf];
s6 = [inf];
s7 = [inf];
s8 = [inf];
Td1 = 1; % Number of scheduled delayed reactions.
Td2 = 1;
Td3 = 1;
Td4 = 1;
Td5 = 1;
Td6 = 1;
Td7 = 1;
Td8 = 1;

% Initialize the states. 
T = 0; % Time
Delta = 0; % Initialize the time step taken at each iteration.

%c1=1; %c1h1h1=0; c1h6h7=0; % Gene 1 levels
%c2=1; %c2h1h1=0; c2h6h7=0; % Gene 7 levels
mh1=0; mh7=0; mh6=0; % Her1, her7 and her 6 mRNA levels (added her6)
ph1=0; ph7=0; ph6=100; % Her1, Her7 and Hes6 protein levels

ph11=0;ph16=0;ph17=0;
ph66=0;ph76=0;ph77=0;
md=0;pd=0; 

a = zeros([34,1]); % Initialize the propensity levels for each reaction.

% Initialize Poisson times and next Poisson Jumps. Also set t s.
Tk = zeros([34,1]);
Pk = zeros([34,1]);
t = zeros([34,1]);

for i=1:34
    r=rand;
    Pk(i) = log(1/r);
end

for i=1:maxi-1 % Run the stochastic simulationm for max number of iterations.
    
  if T >= tend % If the simulation time passes tend finish the run.
    break
  end
  
%  Calculation of transcription rates of her1,her7 and deltaC mRNAs
fh1=msh1*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);
fh7=msh7*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);
fd=msd/(1+(ph11/critph11)^2+(ph76/critph76)^2);
    
% 1-8: delayed reactions
% 9-34: non-delayed reactions

a(1) = psh1*mh1; % Reaction 01: mh1 -> ph1	
a(2) = psh7*mh7;% // Reaction 09: mh7 -> ph7
a(3) = psh6*mh6; % // Reaction 15: mh6 -> ph6
a(4) = psd*md; % // Reaction 25: md -> pd	
%a(27) = fh1(msh1, ph11, ph76, pd(), critph, critpd); % // Reaction 27: -> mh1 
a(5)=fh1; % // Reaction 27: -> mh1 
%a(29) = fh7(msh7, ph11, ph76, pd(), critph, critpd); % // Reaction 29: -> mh7	
a(6)=fh7; % // Reaction 29: -> mh7
a(7) = psh6; % // Reaction 31: -> mh6	
%a(33) = fd(msd, ph11, ph76, critph, critpd); % // Reaction 33: -> md	
a(8)=fd;% // Reaction 33: -> md
a(9) = pdh1*ph1; % Reaction 02: ph1 ->
a(10) = pdh7*ph7;% // Reaction 10: ph7 -> 	
a(11) = dah7h7*ph7*(ph7-1)/2;% // Reaction 11: ph7+ph7 -> ph77
a(12) = ddh7h7*ph77; % // Reaction 12: ph77 -> ph7+ph7	
a(13) = dah7h6*ph7*ph6;% // Reaction 6: ph7+ph6 -> ph76	
a(14) = ddh7h6*ph76; % // Reaction 14: ph76 -> ph7+ph6	 
a(15) = dah1h1*ph1*(ph1-1)/2; % Reaction 03: ph1+ph1 -> ph11
a(16) = pdh6*ph6; % // Reaction 16: ph6 ->
a(17) = dah6h6*ph6*(ph6-1)/2; % // Reaction 17: ph6+ph6 -> ph66
a(18) = ddh6h6*ph66; % // Reaction 18: ph66 -> ph6+ph6
a(19) = pdh1*ph11; % // Reaction 19: ph11 -> 	 
a(20) = pdh1*ph17; % // Reaction 20: ph17 -> 
a(21) = pdh1*ph16; % // Reaction 21: ph16 ->  
a(22) = pdh1*ph77; % // Reaction 22: ph77 -> 	
a(23) = pdh1*ph76; % // Reaction 23: ph76 -> 	
a(24) = pdh1*ph66; % // Reaction 24: ph66 -> 	
a(25) = ddh1h1*ph11; %  Reaction 04: ph11 -> ph1+ph1
a(26) = pdd*pd; % // Reaction 26: pd -> 	
a(27) = dah1h7*ph1*ph7; % Reaction 05: ph1+ph7 -> ph17
a(28) = mdh1*mh1; % // Reaction 28: mh1 -> 	
a(29) = ddh1h7*ph17; % Reaction 06: ph17 -> ph1+ph7
a(30) = mdh7*mh7; % // Reaction 30: mh7 -> 
a(31) = dah1h6*ph1*ph6; % // Reaction 07: ph1+ph6 -> ph16
a(32) = mdh6*mh6; %) // Reaction 32: mh6 ->
a(33) = ddh1h6*ph16;% // Reaction 08: ph16 -> ph1+ph6
a(34) = mdd*md; %// Reaction 34: md ->	
   %need to add mdh6,mh6,pdh6,msh6,psh6 levels
  
  for j=1:34
    if a(j) ~= 0
        t(j) = (Pk(j) - Tk(j))/a(j);
    else
    t(j) = inf;
    end
  end

  % Choose the next reaction and its time.
  [Delta RN] = min([t',s1(1),s2(1),s3(1),s4(1),s5(1),s6(1),s7(1),s8(1)]); 
  % Delta will get the time before the next reaction
  % and RN will the number of the next reaction (which reaction).
  
  if RN>=1 && RN<=42 % 26 non-delayed reactions + 2*8 for delayed reactions 
     s1 = s1 - Delta; % Reaction 1
     s2 = s2 - Delta; % Reaction 2
     s3 = s3 - Delta; % Reaction 3
     s4 = s4 - Delta; % Reaction 4
     s5 = s5 - Delta; % Reaction 5
     s6 = s6 - Delta; % Reaction 6
     s7 = s7 - Delta; % Reaction 7
     s8 = s8 - Delta; % Reaction 8
     T=T+Delta;
  end
  
  % Delayed reactions being completed if necessary.
  
  if RN==35 % Reaction 1: mh1 -> ph1	
     ph1=ph1+1;
     s1=s1(1,2:Td1);
     Td1=Td1-1;
  end
  
  if RN==36 % Reaction 2: mh7 -> ph7	
     ph7=ph7+1;
     s2=s2(1,2:Td2);
     Td2=Td2-1;
  end
  
  if RN==37 %// Reaction 3: mh6 -> ph6
     ph6=ph6+1;
     s3=s3(1,2:Td3);
     Td3=Td3-1;
  end
  
  if RN==38 % // Reaction 4: md -> pd	
     pd=pd+1;
     s4=s4(1,2:Td4);
     Td4=Td4-1;
  end
  
  if RN==39 % // Reaction 5: -> mh1 
      mh1=mh1+1;
      s5=s5(1,2:Td5);
      Td5=Td5-1;
  end
  
  if RN==40 % Reaction 6: -> mh7
      mh7=mh7+1;
      s6=s6(1,2:Td6);
      Td6=Td6-1;
  end
  
  if RN==41  % Reaction 7: -> mh6
     mh6=mh6+1;
     s7=s7(1,2:Td7);
     Td7=Td7-1; 
  end
  
  if RN==42 % // Reaction 8: -> md
     md=md+1; 
     s8=s8(1,2:Td8); 
     Td8=Td8-1;
  end
  
  % For the selected reaction update Pk. 
  % (if not completing a delayed reaction)
  if RN<=34 && RN>=1
     r = rand;
     Pk(RN) = Pk(RN) + log(1/r);
  end
  
%   % Initiate delayed reactions.
%   if RN==1 % Genes her1 and her7 are paired here
%      s1 = [s1(1:Td1 - 1), nmh1, inf]; 
%      Td1 = Td1 + 1;
%      s2 = [s2(1:Td2 - 1), nmh7, inf]; 
%      Td2 = Td2 + 1;
%   end
%   
%   if RN==2 % Genes her1 and her7 are paired here
%      s2 = [s2(1:Td2 - 1), nmh7, inf]; 
%      Td2 = Td2 + 1;
%   end
  
  if RN==1 % Reaction 1: mh1 -> ph1
     s1 = [s1(1:Td1 - 1), nph1, inf];
     Td1=Td1+1;
  end
      
  if RN==2 % Reaction 2: mh7 -> ph7
     s2 = [s2(1:Td2 - 1), nph7, inf];
     Td2=Td2+1;
  end 
      
  if RN==3 %// Reaction 3: mh6 -> ph6
     s3 = [s3(1:Td3 - 1), nph6, inf]; 
     Td3 = Td3 + 1;
  end
  
  if RN==4 % // Reaction 4: md -> pd
     s4 = [s4(1:Td4 - 1), npd, inf]; 
     Td4 = Td4 + 1;
  end
  
  if RN==5 % // Reaction 5: -> mh1 
     s5 = [s5(1:Td5 - 1), nmh1, inf]; 
     Td5 = Td5 + 1;
  end
  
  if RN==6 % Reaction 6: -> mh7
     s6 = [s6(1:Td6 - 1), nmh7, inf]; 
     Td6 = Td6 + 1;
  end
  
  if RN==7  % Reaction 7: -> mh6
     s7 = [s7(1:Td7 - 1), nmh6, inf]; 
     Td7 = Td7 + 1;
  end
  
  if RN==8 % // Reaction 8: -> md 
     s8 = [s8(1:Td8 - 1), nmd, inf]; 
     Td8 = Td8 + 1;
  end
  
  if RN==9 % Reaction 9: ph1 ->
      ph1=ph1-1;
  elseif RN==10 % Reaction 10: ph7 ->  
      ph7=ph7-1;
  elseif RN==11 % // Reaction 11: ph7+ph7 -> ph77
      ph7=ph7-2;
      ph77=ph77+1;
  elseif RN==12 % // Reaction 12: ph77 -> ph7+ph7  
      ph77=ph77-1;
      ph7=ph7+2;
  elseif RN==13 % // Reaction 13: ph7+ph6 -> ph76 
      ph7=ph7-1;
      ph6=ph6-1;
      ph76=ph76+1;
  elseif RN==14 % // Reaction 14: ph76 -> ph7+ph6   
      ph76=ph76-1;
      ph7=ph7+1;
      ph6=ph6+1;
  elseif RN==15 % Reaction 03: ph1+ph1 -> ph11
      ph1=ph1-2;
      ph11=ph11+1;
  elseif RN==16 % Reaction 16: ph6 ->
      ph6=ph6-1;
  elseif RN==17 % // Reaction 17: ph6+ph6 -> ph66
      ph6=ph6-2;
      ph66=ph66+1;
  elseif RN==18 % // Reaction 18: ph66 -> ph6+ph6
      ph66=ph66-1;
      ph6=ph6+2;
  elseif RN==19 % // Reaction 19: ph11 ->
      ph11=ph11-1;
  elseif RN==20 % // Reaction 20: ph17 ->
      ph17=ph17-1;
  elseif RN==21 % // Reaction 21: ph16 ->
      ph16=ph16-1;
  elseif RN==22 % // Reaction 22: ph77 ->
      ph77=ph77-1;
  elseif RN==23 % // Reaction 23: ph76 ->
      ph76=ph76-1;
  elseif RN==24 % // Reaction 24: ph66 ->
      ph66=ph66-1;
  elseif RN==25 %  Reaction 25: ph11 -> ph1+ph1
      ph11=ph11-1;
      ph1=ph1+2;
  elseif RN==26 %  Reaction 26: pd ->     
      pd=pd-1;
  elseif RN==27 % Reaction 27: ph1+ph7 -> ph17
      ph1=ph1-1;
      ph7=ph7-1;
      ph17=ph17+1;
  elseif RN==28 % // Reaction 28: mh1 ->  
      mh1=mh1-1;
  elseif RN==29 % Reaction 29: ph17 -> ph1+ph7
      ph17=ph17-1;
      ph1=ph1+1;
      ph7=ph7+1;
  elseif RN==30 % // Reaction 30: mh7 -> 
      mh7=mh7-1;
  elseif RN==31 % // Reaction 07: ph1+ph6 -> ph16
      ph1=ph1-1;
      ph6=ph6-1;
      ph16=ph16+1;
  elseif RN==32 % // Reaction 32: mh6 ->
      mh6=mh6-1;
  elseif RN==33 % // Reaction 08: ph16 -> ph1+ph6
      ph16=ph16-1;
      ph1=ph1+1;
      ph6=ph6+1;
  elseif RN==34 %// Reaction 34: md ->  
      md=md-1;
  end
  
  % Update Tk values for each reaction.
  for k=1:34
     Tk(k)=Tk(k)+a(k)*Delta;
  end
    
    Time=[Time T]; % Store time, and mh1 and mh7 levels.
    mh1v=[mh1v mh1];
    mh7v=[mh7v mh7];
end % for i=1:maxi-1,

Data = [Time' mh1v' mh7v'];
DataTable=table(Time', mh1v', mh7v','VariableNames', {'Time', 'mh1','mh7'});
writetable(DataTable, strcat('Run', num2str(run),'.xlsx'),'WriteVariableNames', true);

figure
plot(Time, mh1v,'b')
hold on
plot(Time, mh7v,'r')
legend('Her1','Her7')
xlabel('Time')
ylabel('#mRNA')
saveas(gcf,strcat(['Run', num2str(run)]),'jpg');
close(gcf);

figure
histogram(mh1v)
saveas(gcf,strcat('Run_Histogram', num2str(run)),'jpg');
close(gcf)

[RHO1(run),PVAL1(run)] = corr(mh1v',mh7v','Type','Pearson');
[RHO2(run),PVAL2(run)] = corr(mh1v',mh7v','Type','Spearman');

end
Stats=table(RHO1, PVAL1, RHO2, PVAL2,'VariableNames', {'Pearson', 'PearSig','Spearman','SpearSig'});
writetable(Stats, strcat('Run_PearsonSpearman.xlsx'),'WriteVariableNames', true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
