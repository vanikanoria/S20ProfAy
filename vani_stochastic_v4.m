% 2-cell system
% using approximations for dimer formation and degradation
% with dependency structure
% ph6 : constant, and no mh6
% Updating propensities wherever updating concentration

nrun=1;%30;
RHO1=zeros([nrun,1]); PVAL1=zeros([nrun,1]);
RHO2=zeros([nrun,1]); PVAL2=zeros([nrun,1]);
RHO3=zeros([nrun,1]); PVAL3=zeros([nrun,1]);

% % Parameter values used in the 2013 study (from Table S1)
psh1=49.9139;
% psh6=34.3117;
psh7=28.5626;
psd=37.7828;
pdh1=0.34951;
% pdh6=0.14824;
pdh7=0.249715;
pdd=0.324316;
msh1=30; % 48.3084;
% msh6=36.4073;
msh7= 20; %39.685;
msd=60.5577;
mdh1=0.322965;
% mdh6=0.146372;
mdh7=0.381738;
mdd=0.352056;
% pdh11=0.390961;
% pdh16=0.29774;
% pdh17=0.320157;
% pdh66=0.268042;
% pdh67=0.352037;
% pdh77=0.251601;
nmh1=10.0213;
nmh6=10;
nmh7=10.4515;
nmd=7.74472;
nph1=1.5398;
nph6=0.886233;
nph7=0.539972;
npd=13.2661;
% dah1h1=0.0179429;
% ddh1h1=0.220856;
% dah1h6=0.0270209;
% ddh1h6=0.0917567;
% dah1h7=0.00120525;
% ddh1h7=0.258167;
% dah6h6=0.0148271;
% ddh6h6=0.251173;
% dah7h6=0.0216093;
% ddh7h6=0.188923;
% dah7h7=0.0202756;
% ddh7h7=0.161018;
critph11=587.298;
critph76=769.628;
critpd=490.254;

for run=1:nrun 
maxi=10000000; % Maximum number of iterations to run.
tend=20; %240; % Maximum time for simulation.
format compact 
% rand('state',sum(100*clock)) 
rand('state',run) % Sets random number generator to a specific state.

Time=[0]; % Simulation time  
mh1v_c1=[0]; % Her1 mRNA levels -> cell 1
mh7v_c1=[0]; % Her7 mRNA levels -> cell 1
mh1v_c2=[0]; % Her1 mRNA levels -> cell 2
mh7v_c2=[0]; % Her7 mRNA levels -> cell 2


% Data storage for delayed reaction times
% row 1 -> cell 1
% row 2 -> cell 2

s = cell(2,6,1);
each_cell={inf};
% Time for each delayed reaction is scheduled in these vectors.
s(:,:)=each_cell;

Td1 = [1;1]; % Number of scheduled delayed reactions.
Td2 = [1;1];
Td3 = [1;1];
Td4 = [1;1];
Td5 = [1;1];
Td6 = [1;1];

% Initialize the states. 
T = 0; % Time
Delta = 0; % Initialize the time step taken at each iteration.

%States for each of the cells:
mh1=[0;0]; mh7=[0;0]; % Her1 and her7 levels 
ph1=[0;0]; ph7=[0;0]; ph6=[100;100]; % Her1, Her7 and Hes6 protein levels
md=[0;0];pd=[0;0]; 

a = zeros([2,12]); % Initialize the propensity levels for each reaction for both cells.

% Initialize Poisson times and next Poisson Jumps. Also set t s.
Tk = zeros([2,12]);
Pk = zeros([2,12]);
t = zeros([2,12]);

for ck=1:2 % for 2 cells
    for i=1:12
        r=rand;
        Pk(ck, i) = log(1/r);
    end
end

% Functions for calculation of transcription rates of her1,her7 and deltaC mRNAs

% I need: if ck=1, cn=2; if ck=2,cn=1

% Ideally need : fh1= @(pd(cn), ph11(ck), ph76(ck)) ...
fh1(1) = msh1*(1+(pd(2)/critpd))/(1+(pd(2)/critpd)+(ph1(1)*ph1(1)/critph11)^2 +(ph7(1)*ph6(1)/critph76)^2);
fh1(2) = msh1*(1+(pd(1)/critpd))/(1+(pd(1)/critpd)+(ph1(2)*ph1(2)/critph11)^2 +(ph7(2)*ph6(2)/critph76)^2);
fh7(1)= msh7*(1+(pd(2)/critpd))/(1+(pd(2)/critpd)+(ph1(1)*ph1(1)/critph11)^2 +(ph7(1)*ph6(1)/critph76)^2);
fh7(2)= msh7*(1+(pd(1)/critpd))/(1+(pd(1)/critpd)+(ph1(2)*ph1(2)/critph11)^2 +(ph7(2)*ph6(2)/critph76)^2);
fd=@(ph1,ph7,ph6) msd/(1+(ph1*ph1/critph11)^2+(ph7*ph6/critph76)^2);

for ck=1:2
    a(ck,1) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
    a(ck,2) = psh7*mh7(ck);% // Reaction 02: mh7 -> ph7
    a(ck,3)=fh7(ck); % // Reaction 3: -> mh7
    a(ck,4) = psd*md(ck); % // Reaction 4: md -> pd
    a(ck,5)=fh1(ck); % // Reaction 5: -> mh1
    a(ck,6)=fd(ph1(ck),ph7(ck),ph6(ck));% // Reaction 6: -> md
    a(ck,7) = mdh1*mh1(ck); % // Reaction 7: mh1 ->
    a(ck,8) = mdh7*mh7(ck); % // Reaction 8: mh7 ->
    a(ck,9) = pdh1*ph1(ck); % Reaction 09: ph1 ->
    a(ck,10) = pdh7*ph7(ck);% // Reaction 10: ph7 ->
    a(ck,11) = mdd*md(ck); %// Reaction 11: md ->
    a(ck,12) = pdd*pd(ck); % // Reaction 12: pd ->   
    
end

for i=1:maxi-1 % Run the stochastic simulationm for max number of iterations.
    
  if T >= tend % If the simulation time passes tend finish the run.
    break
  end
  
% 1-6: delayed reactions
% 7-12: non-delayed reactions

   for ckvalue=1:2
       for j=1:12
           if a(ckvalue,j) ~= 0
               t(ckvalue,j) = (Pk(ckvalue,j) - Tk(ckvalue,j))/a(ckvalue,j); 
    
           else
               t(ckvalue,j) = inf;
           end
       end
   end

  % Choose the next reaction and its time.
  [Delta RN] = min([t(1,:),s{1,1}(1),s{1,2}(1),s{1,3}(1), s{1,4}(1), s{1,5}(1), s{1,6}(1), t(2,:),s{2,1}(1),s{2,2}(1),s{2,3}(1),s{2,4}(1),s{2,5}(1),s{2,6}(1)]); 
  
  % = [non-delayed reactions of cell 1, delayed reactions of cell 1,
  % non-delayed reactions of cell 2, delayed reactions of cell 2)
  
  % Delta will get the time before the next reaction
  % and RN will the number of the next reaction (which reaction).
  
  % potential no. of reactions = 12+6+12+6
  
  if RN<=18 %this is if the reaction is in cell 1
      ck=1;
      cn=2; % other cell is cell 2
  else %the reaction is in cell 2 
      ck=2; 
      cn=1; %other cell is cell 1
      RN=RN-18;
  end

if RN>=1 && RN<=18
    for ckvalue=1:2
        for rn=1:6
            s{ckvalue,rn}=s{ckvalue,rn}-Delta;
        end
    end
    T=T+Delta;
end

% Delayed reactions being completed if necessary.
% cell 1
  if RN==13 % Reaction 1 for cell 1: mh1 -> ph1
     ph1=ph1+1;
     s{ck,RN-12}=s{ck,RN-12}(2:Td1(ck));
     Td1(ck)=Td1(ck)-1;
     a(ck,3)=fh7(ck); % // Reaction 3: -> mh7
     a(ck,5)=fh1(ck); % // Reaction 5: -> mh1
     a(ck,6)=fd(ph1(ck),ph7(ck),ph6(ck));% // Reaction 6: -> md
     a(ck,9) = pdh1*ph1(ck); % Reaction 09: ph1 ->
  end
  
  if RN==14 % Reaction 2: mh7 -> ph7	
     ph7=ph7+1;
     s{ck,RN-12}=s{ck,RN-12}(2:Td2(ck));
     Td2(ck)=Td2(ck)-1;
     a(ck,3)=fh7(ck); % // Reaction 3: -> mh7
     a(ck,5)=fh1(ck); % // Reaction 5: -> mh1
     a(ck,6)=fd(ph1(ck),ph7(ck),ph6(ck));% // Reaction 6: -> md
     a(ck,10) = pdh7*ph7(ck);% // Reaction 10: ph7 ->
  end
  
  if RN==15 %// Reaction 3: -> mh7
     mh7=mh7+1;
     s{ck,RN-12}=s{ck,RN-12}(2:Td3(ck));
     Td3(ck)=Td3(ck)-1;
     a(ck,2) = psh7*mh7(ck);% // Reaction 09: mh7 -> ph7
     a(ck,8) = mdh7*mh7(ck); % // Reaction 8: mh7 ->
  end
  
  if RN==16 % // Reaction 4: md -> pd
     pd=pd+1;
     s{ck,RN-12}=s{ck,RN-12}(2:Td4(ck));
     Td4(ck)=Td4(ck)-1;
     a(ck,3)=fh7(ck); % // Reaction 3: -> mh7
     a(ck,5)=fh1(ck); % // Reaction 5: -> mh1
     a(ck,12) = pdd*pd(ck); % // Reaction 12: pd ->
  end
  
  if RN==17 % // Reaction 5: -> mh1 
      mh1=mh1+1;
      s{ck,RN-12}=s{ck,RN-12}(2:Td5(ck));
      Td5(ck)=Td5(ck)-1;
      a(ck,1) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
      a(ck,7) = mdh1*mh1(ck); % // Reaction 7: mh1 ->
  end
  
  if RN==18 % Reaction 6: -> md
      md=md+1;
      s{ck,RN-12}=s{ck,RN-12}(2:Td6(ck));
      Td6(ck)=Td6(ck)-1;
      a(ck,4) = psd*md(ck);
      a(ck,11) = mdd*md(ck); %// Reaction 11: md ->
  end
    
  % For the selected reaction and cell, update Pk. 
  % (if not completing a delayed reaction)
    
  if RN<=12 && RN>=1
     r = rand;
     Pk(ck,RN) = Pk(ck,RN) + log(1/r);
  end

  if RN==1 % Reaction 1: mh1 -> ph1
     s{ck,RN} = s{ck,RN} (1:Td1(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nph1, inf];
     Td1(ck)=Td1(ck)+1;
  end
      
  if RN==2 % Reaction 2: mh7 -> ph7
     s{ck,RN} = s{ck,RN} (1:Td2(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nph7, inf];
     Td2(ck)=Td2(ck)+1;
  end 
      
  if RN==3 %// Reaction 3: -> mh7
     s{ck,RN} = s{ck,RN} (1:Td3(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nmh7, inf];
     Td3(ck) = Td3(ck) + 1;
  end
  
  if RN==4 % // Reaction 4: md -> pd
     s{ck,RN} = s{ck,RN} (1:Td4(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, npd, inf]; 
     Td4(ck) = Td4(ck) + 1;
  end
  
  if RN==5 % // Reaction 5: -> mh1 
     s{ck,RN} = s{ck,RN} (1:Td5(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nmh1, inf]; 
     Td5(ck) = Td5(ck) + 1;
  end
  
  if RN==6 % Reaction 6: -> md
     s{ck,RN} = s{ck,RN} (1:Td6(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nmd, inf];
     Td6(ck) = Td6(ck) + 1;
  end
  
  if RN==7  % Reaction 7: mh1 -> 
     mh1(ck)=mh1(ck)-1;
     a(ck,1) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
     a(ck,7) = mdh1*mh1(ck); % // Reaction 7: mh1 ->
  end
  
  if RN==8 % // Reaction 8: mh7 ->
     mh7(ck)=mh7(ck)-1;
     a(ck,2) = psh7*mh7(ck);% // Reaction 09: mh7 -> ph7
     a(ck,8) = mdh7*mh7(ck); % // Reaction 8: mh7 ->
  end
  
  if RN==9 % Reaction 9: ph1 ->
      ph1(ck)=ph1(ck)-1;
      a(ck,3)=fh7(ck); % // Reaction 3: -> mh7
      a(ck,5)=fh1(ck); % // Reaction 5: -> mh1
      a(ck,6)=fd(ph1(ck),ph7(ck),ph6(ck));% // Reaction 6: -> md
      a(ck,9) = pdh1*ph1(ck); % Reaction 09: ph1 ->
  elseif RN==10 % Reaction 10: ph7 ->  
      ph7(ck)=ph7(ck)-1;
      a(ck,3)=fh7(ck); % // Reaction 3: -> mh7
      a(ck,5)=fh1(ck); % // Reaction 5: -> mh1
      a(ck,6)=fd(ph1(ck),ph7(ck),ph6(ck));% // Reaction 6: -> md
      a(ck,10) = pdh7*ph7(ck);% // Reaction 10: ph7 ->
  elseif RN==11 % Reaction 11: md ->
      md(ck)=md(ck)-1;
      a(ck,4) = psd*md(ck);
      a(ck,11) = mdd*md(ck); %// Reaction 11: md ->
  elseif RN==12 %  Reaction 12: pd ->     
      pd(ck)=pd(ck)-1;    
      a(ck,3)=fh7(ck); % // Reaction 3: -> mh7
      a(ck,5)=fh1(ck); % // Reaction 5: -> mh1
      a(ck,12) = pdd*pd(ck); % // Reaction 12: pd ->
  end

  % Update Tk values for each reaction and each cell.
  for ckvalue=1:2
      for k=1:12
          Tk(ckvalue,k)=Tk(ckvalue,k)+a(ckvalue,k)*Delta;
      end
  end
  
    Time=[Time T]; % Store time, and mh1 and mh7 levels.
    mh1v_c1=[mh1v_c1 mh1(1)];
    mh1v_c2=[mh1v_c2 mh1(2)];
    mh7v_c1=[mh7v_c1 mh7(1)];
    mh7v_c2=[mh7v_c2 mh7(2)];
end % for i=1:maxi-1,

Data = [Time' mh1v_c1' mh7v_c1' mh1v_c2' mh7v_c2'];
% DataTable=table(Time', mh1v_c1, mh7v_c1', mh1v_c2', mh7v_c2','VariableNames', {'Time', 'mh1_cell1', 'mh7_cell1', 'mh1_cell2', 'mh7_cell2'});
% writetable(DataTable, strcat('Run_v2', num2str(run),'.xlsx'),'WriteVariableNames', true);

figure
plot(Time, mh1v_c1,'b')
hold on
plot(Time, mh7v_c1,'r')
hold on
plot(Time, mh1v_c2,'b')
hold on
plot(Time, mh7v_c2,'r')
legend('Her1 of cell 1','Her7 of cell 1','Her1 of cell 2','Her7 of cell 2')
xlabel('Time')
ylabel('#mRNA')
saveas(gcf,strcat(['Run_v2', num2str(run)]),'jpg');
close(gcf);

figure
histogram(mh1v_c1)
saveas(gcf,strcat('Run_v2_Histogram', num2str(run)),'jpg');
close(gcf)

% only covered cell 1 here
[RHO1(run),PVAL1(run)] = corr(mh1v_c1',mh7v_c1','Type','Pearson');
[RHO2(run),PVAL2(run)] = corr(mh1v_c1',mh7v_c1','Type','Spearman');

end
Stats=table(RHO1, PVAL1, RHO2, PVAL2,'VariableNames', {'Pearson', 'PearSig','Spearman','SpearSig'});
writetable(Stats, strcat('Run_v2_PearsonSpearman.xlsx'),'WriteVariableNames', true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
