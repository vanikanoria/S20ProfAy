nrun=30;
RHO1=zeros([nrun,1]); PVAL1=zeros([nrun,1]);
RHO2=zeros([nrun,1]); PVAL2=zeros([nrun,1]);
RHO3=zeros([nrun,1]); PVAL3=zeros([nrun,1]);

for run=1:nrun 
maxi=10000000; % Maximum number of iterations to run.
tend=240; % Maximum time for simulation.
format compact 
% rand('state',sum(100*clock)) 
rand('state',run) % Sets random number generator to a specific state.

Time=[0]; % Simulation time  
mh1v=[0]; % Her1 mRNA levels
mh7v=[0]; % Her7 mRNA levels

% Reaction rates.
msh1=5; msh7=5; mdh1=0.2; mdh7=0.2; % mRNA synthesis and degradation rates.

psh1=5; psh7=5; pdh1=0.3; pdh7=0.3; % Protein synthesis and degradation rates.

ac1h1h1=0.01; dc1h1h1=400; % DNA association and dissocation rates.
ac2h1h1=0.01; dc2h1h1=400;
ac1h6h7=0.01; dc1h6h7=400;
ac2h6h7=0.01; dc2h6h7=400;

nmh1=9; nmh7=9; % mRNA and protein synthesis time delays.
nph1=1; nph7=1;

% Data storage for delayed reaction times
s1 = [inf]; % Time for each delayed reaction is scheduled in these vectors.
s2 = [inf];
s3 = [inf];
s4 = [inf];
Td1 = 1; % Number of scheduled delayed reactions.
Td2 = 1;
Td3 = 1;
Td4 = 1;

% Initialize theh states. 
T = 0; % Time
Delta = 0; % Initialize the time step taken at each iteration.

c1=1; c1h1h1=0; c1h6h7=0; % Gene 1 levels
c2=1; c2h1h1=0; c2h6h7=0; % Gene 7 levels
mh1=0; mh7=0; % Her1 and her7 mRNA levels
ph1=0; ph7=0; ph6=100; % Her1, Her7 and Hes6 protein levels

a = zeros([16,1]); % Initialize the propensity levels for each reaction.

% Initialize Poisson times and next Poisson Jumps. Also set t s.
Tk = zeros([16,1]);
Pk = zeros([16,1]);
t = zeros([16,1]);

for i=1:16
    r=rand;
    Pk(i) = log(1/r);
end

for i=1:maxi-1 % Run the stochastic simulationm for max number of iterations.
    
  if T >= tend % If the simulation time passes tend finish the run.
    break
  end

  % Propensity calculations
  a(1)=msh1*c1;
  a(2)=msh7*c2;
  a(3)=psh1*mh1;
  a(4)=psh7*mh7;
  a(5)=mdh1*mh1;
  a(6)=mdh7*mh7;
  a(7)=pdh1*ph1;
  a(8)=pdh7*ph7;
  a(9)=ac1h1h1*c1*ph1*(ph1-1)/2;
  a(10)=ac1h6h7*c1*ph6*ph7;
  a(11)=ac2h1h1*c2*ph1*(ph1-1)/2;
  a(12)=ac2h6h7*c2*ph6*ph7;
  a(13)=dc1h1h1*c1h1h1;
  a(14)=dc1h6h7*c1h6h7;
  a(15)=dc2h1h1*c2h1h1;
  a(16)=dc2h6h7*c2h6h7;
  
  for j=1:16
    if a(j) ~= 0
        t(j) = (Pk(j) - Tk(j))/a(j);
    else
    t(j) = inf;
    end
  end

  [Delta RN] = min([t',s1(1),s2(1),s3(1),s4(1)]); % Choose the next reaction and its time.
  
  if RN>=1 && RN<=20 
     s1 = s1 - Delta;
     s2 = s2 - Delta;
     s3 = s3 - Delta;
     s4 = s4 - Delta;
     T=T+Delta;
  end
  
  % Delayed reactions being completed if necessary.
  if RN==17
     mh1=mh1+1;
     s1=s1(1,2:Td1);
     Td1=Td1-1;
  end

  if RN==18
     mh7=mh7+1;
     s2=s2(1,2:Td2);
     Td2=Td2-1;
  end

  if RN==19
     ph1=ph1+1;
     s3=s3(1,2:Td3);
     Td3=Td3-1;
  end
  
  if RN==20
     ph7=ph7+1;
     s4=s4(1,2:Td4);
     Td4=Td4-1;
  end
  
  % For the selected reaction update Pk. 
  if RN<=16 && RN>=1
     r = rand;
     Pk(RN) = Pk(RN) + log(1/r);
  end
  
  % Initiate delayed reactions.
  if RN==1 % Genes her1 and her7 are paired here
     s1 = [s1(1:Td1 - 1), nmh1, inf]; 
     Td1 = Td1 + 1;
     s2 = [s2(1:Td2 - 1), nmh7, inf]; 
     Td2 = Td2 + 1;
  end
  
%   if RN==2 % Genes her1 and her7 are paired here
%      s2 = [s2(1:Td2 - 1), nmh7, inf]; 
%      Td2 = Td2 + 1;
%   end
  
  if RN==3
     s3 = [s3(1:Td3 - 1), nph1, inf]; 
     Td3 = Td3 + 1;
  end
  
  if RN==4
     s4 = [s4(1:Td4 - 1), nph7, inf]; 
     Td4 = Td4 + 1;
  end
  
  % Non-delayed reactions.
  if RN==5
     mh1=mh1-1;
  elseif RN==6
      mh7=mh7-1;
  elseif RN==7
      ph1=ph1-1;
  elseif RN==8
      ph7=ph7-1;
  elseif RN==9
      c1=c1-1;
      ph1=ph1-2;
      c1h1h1=c1h1h1+1;
  elseif RN==10
      c1=c1-1;
      %ph6=ph6-1;
      ph7=ph7-1;
      c1h6h7=c1h6h7+1;
  elseif RN==11
      c2=c2-1;
      ph1=ph1-2;
      c2h1h1=c2h1h1+1;
  elseif RN==12
      c2=c2-1;
      %ph6=ph6-1;
      ph7=ph7-1;
      c2h6h7=c2h6h7+1;
  elseif RN==13
      c1=c1+1;
      ph1=ph1+2;
      c1h1h1=c1h1h1-1;
  elseif RN==14
      c1=c1+1;
      %ph6=ph6+1;
      ph7=ph7+1;
      c1h6h7=c1h6h7-1;
  elseif RN==15
      c2=c2+1;
      ph1=ph1+2;
      c2h1h1=c2h1h1-1;
  elseif RN==16
      c2=c2+1;
      %ph6=ph6+1;
      ph7=ph7+1;
      c2h6h7=c2h6h7-1; 
  end
  
  % Update Tk values for each reaction.
  for i=1:16
     Tk(i)=Tk(i)+a(i)*Delta;
  end
    
    Time=[Time T]; % Store time, and mh1 and mh7 levels.
    mh1v=[mh1v mh1];
    mh7v=[mh7v mh7];
end % for i=1:maxi-1,

Data = [Time' mh1v' mh7v'];
DataTable=table(Time', mh1v', mh7v','VariableNames', {'Time', 'mh1','mh7'});
writetable(DataTable, strcat('GenePaired/GenePaired_Run', num2str(run),'.xlsx'),'WriteVariableNames', true);

figure
plot(Time, mh1v,'b')
hold on
plot(Time, mh7v,'r')
legend('Her1','Her7')
xlabel('Time')
ylabel('#mRNA')
saveas(gcf,strcat(['GenePaired/GenePaired_Run', num2str(run)]),'jpg');
close(gcf);

figure
histogram(mh1v)
saveas(gcf,strcat('GenePaired/GenePaired_Histogram', num2str(run)),'jpg');
close(gcf)

[RHO1(run),PVAL1(run)] = corr(mh1v',mh7v','Type','Pearson');
[RHO2(run),PVAL2(run)] = corr(mh1v',mh7v','Type','Spearman');

end
Stats=table(RHO1, PVAL1, RHO2, PVAL2,'VariableNames', {'Pearson', 'PearSig','Spearman','SpearSig'});
writetable(Stats, strcat('GenePaired/GenePaired_PearsonSpearman.xlsx'),'WriteVariableNames', true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
