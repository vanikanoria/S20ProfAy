% 2-cell system using approximations for dimer formation and degradation
function Data = vani_stochastic_v3(param_set)

psh1=param_set(1); 
psh6=param_set(2);
psh7=param_set(3);
psd=param_set(4);
pdh1=param_set(5);
pdh6=param_set(6);
pdh7=param_set(7);
pdd=param_set(8);
msh1=param_set(9);
msh6=param_set(10);
msh7=param_set(11);
msd=param_set(12);
mdh1=param_set(13);
mdh6=param_set(14);
mdh7=param_set(15);
mdd=param_set(16);
pdh11=param_set(17);
pdh16=param_set(18);
pdh17=param_set(19);
pdh66=param_set(20);
pdh76=param_set(21);
pdh77=param_set(22);
nmh1=param_set(23);
nmh7=param_set(24);
nmd=param_set(25);
nph1=param_set(26);
nph6=param_set(27);
nph7=param_set(28);
npd=param_set(29);
dah1h1=param_set(30);
ddh1h1=param_set(31);
dah1h6=param_set(32);
ddh1h6=param_set(33);
dah1h7=param_set(34);
ddh1h7=param_set(35);
dah6h6=param_set(36);
ddh6h6=param_set(37);
dah7h6=param_set(38);
ddh7h6=param_set(39);
dah7h7=param_set(40);
ddh7h7=param_set(41);
critph11=param_set(42);
critph76=param_set(43);
critpd=param_set(44);
nmh6=10;

maxi = 10000000; % Maximum number of iterations to run.
tend = 100; % 240; % Maximum time for simulation.
format compact 
% rand('state',sum(100*clock)) 
rand('state',1) % Sets random number generator to a specific state.

Time=[0]; % Simulation time  
mh1v_c1=[0]; % Her1 mRNA levels -> cell 1
mh1v_c2=[0]; % Her1 mRNA levels -> cell 2


% Data storage for delayed reaction times
% row 1 -> cell 1
% row 2 -> cell 2

s = cell(2,8,1);
each_cell={inf};
% Time for each delayed reaction is scheduled in these vectors.
s(:,:)=each_cell;

Td1 = [1;1]; % Number of scheduled delayed reactions.
Td2 = [1;1];
Td3 = [1;1];
Td4 = [1;1];
Td5 = [1;1];
Td6 = [1;1];
Td7 = [1;1];
Td8 = [1;1];

% Initialize the states. 
T = 0; % Time
Delta = 0; % Initialize the time step taken at each iteration.

% States for each of the cells:
mh1=[0;0]; mh7=[0;0]; mh6=[0;0]; md=[0;0]; % Her1, her7, hes6 and delta mRNA levels 
ph1=[0;0]; ph7=[0;0]; ph6=[100;100]; pd=[0;0];  % Her1, Her7, Hes6 and Delta protein levels

a = zeros([2,16]); % Initialize the propensity levels for each reaction for both cells.

% Initialize Poisson times and next Poisson Jumps. Also set t s.
Tk = zeros([2,16]);
Pk = zeros([2,16]);
t = zeros([2,16]);

for ckvalue=1:2 % for 2 cells
    for i=1:16
        r=rand;
        Pk(ckvalue, i) = log(1/r);
    end
end

% Functions for calculation of transcription rates of her1, her7 and deltaC mRNAs

% If ck=1, cn=2; if ck=2,cn=1
fh1=[0;0];fh7=[0;0];
% Ideally need : fh1= @(pd(cn), ph11(ck), ph76(ck)) ...
fh1(1) = msh1*(1+(pd(2)/critpd))/(1+(pd(2)/critpd)+(ph1(1)*ph1(1)/critph11)^2 +(ph7(1)*ph6(1)/critph76)^2);
fh1(2) = msh1*(1+(pd(1)/critpd))/(1+(pd(1)/critpd)+(ph1(2)*ph1(2)/critph11)^2 +(ph7(2)*ph6(2)/critph76)^2);
fh7(1)= msh7*(1+(pd(2)/critpd))/(1+(pd(2)/critpd)+(ph1(1)*ph1(1)/critph11)^2 +(ph7(1)*ph6(1)/critph76)^2);
fh7(2)= msh7*(1+(pd(1)/critpd))/(1+(pd(1)/critpd)+(ph1(2)*ph1(2)/critph11)^2 +(ph7(2)*ph6(2)/critph76)^2);
fd=@(ph1,ph7,ph6) msd/(1+(ph1*ph1/critph11)^2+(ph7*ph6/critph76)^2);

for i=1:maxi-1 % Run the stochastic simulationm for max number of iterations.
    
  if T >= tend % If the simulation time passes tend finish the run.
    break
  end
  
  if mod(i,100000)==0
      disp(T);
  end
  
% 1-8: delayed reactions
% 9-16: non-delayed reactions

for ck=1:2
    a(ck,1) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
    a(ck,2) = psh7*mh7(ck);% // Reaction 09: mh7 -> ph7
    a(ck,3) = psh6*mh6(ck); % // Reaction 15: mh6 -> ph6
    a(ck,4) = psd*md(ck); % // Reaction 25: md -> pd
    a(ck,5)= fh1(ck); % // Reaction 27: -> mh1
    a(ck,6)= fh7(ck); % // Reaction 29: -> mh7
    a(ck,7) = msh6; % // Reaction 31: -> mh6
    a(ck,8)= fd(ph1(ck),ph7(ck),ph6(ck));% // Reaction 33: -> md
    a(ck,9) = pdh1*ph1(ck); % Reaction 02: ph1 ->
    a(ck,10) = pdh7*ph7(ck);% // Reaction 10: ph7 ->
    a(ck,11) = pdh6*ph6(ck); % // Reaction 16: ph6 ->
    a(ck,12) = pdd*pd(ck); % // Reaction 26: pd ->
    a(ck,13) = mdh1*mh1(ck); % // Reaction 28: mh1 ->
    a(ck,14) = mdh7*mh7(ck); % // Reaction 30: mh7 ->
    a(ck,15) = mdh6*mh6(ck); %) // Reaction 32: mh6 ->
    a(ck,16) = mdd*md(ck); %// Reaction 34: md ->  
end

   for ckvalue=1:2
       for j=1:16
           if a(ckvalue,j) ~= 0
               t(ckvalue,j) = (Pk(ckvalue,j) - Tk(ckvalue,j))/a(ckvalue,j);
           else
               t(ckvalue,j) = Inf;
           end
       end
   end

  % Choose the next reaction and its time.
  % Non-delayed reactions of cell 1, delayed reactions of cell 1
  % Non-delayed reactions of cell 2, delayed reactions of cell 2
  [Delta RN] = min([t(1,:),s{1,1}(1),s{1,2}(1),s{1,3}(1),s{1,4}(1),s{1,5}(1),s{1,6}(1),s{1,7}(1),...
      s{1,8}(1),t(2,:),s{2,1}(1),s{2,2}(1),s{2,3}(1),s{2,4}(1),s{2,5}(1),s{2,6}(1),s{2,7}(1),s{2,8}(1)]); 
  % Delta will get the time before the next reaction
  % and RN will get the number of the next reaction (which reaction).
  % Potential no. of reactions = 16+8+16+8
  
  if RN<=24 % This is if the reaction is in cell 1
      ck=1;
      cn=2; % Other cell is cell 2
  else % The reaction is in cell 2 
      ck=2; 
      cn=1; % Other cell is cell 1
      RN=RN-24;
  end

if RN>=1 && RN<=24
    for ckvalue=1:2
        for rn=1:8
            s{ckvalue,rn}=s{ckvalue,rn}-Delta;
        end
    end
    T=T+Delta;
end

% Delayed reactions being completed if necessary.
% Cells 1 & 2
  if RN==17 % Reaction 1 for cell 1 & 2: mh1 -> ph1
     ph1(ck)=ph1(ck)+1;
     s{ck,RN-16}=s{ck,RN-16}(2:Td1(ck));
     Td1(ck)=Td1(ck)-1;
  end
  
  if RN==18 % Reaction 2: mh7 -> ph7	
     ph7(ck)=ph7(ck)+1;
     s{ck,RN-16}=s{ck,RN-16}(2:Td2(ck));
     Td2(ck)=Td2(ck)-1;
  end
  
  if RN==19 % Reaction 3: mh6 -> ph6
     ph6(ck)=ph6(ck)+1;
     s{ck,RN-16}=s{ck,RN-16}(2:Td3(ck));
     Td3(ck)=Td3(ck)-1;
  end
  
  if RN==20 % Reaction 4: md -> pd     
     pd(ck)=pd(ck)+1;
     s{ck,RN-16}=s{ck,RN-16}(2:Td4(ck));
     Td4(ck)=Td4(ck)-1;
  end
  
  if RN==21 % Reaction 5: -> mh1 
      mh1(ck)=mh1(ck)+1;
      s{ck,RN-16}=s{ck,RN-16}(2:Td5(ck));
      Td5(ck)=Td5(ck)-1;
  end
  
  if RN==22 % Reaction 6: -> mh7
      mh7(ck)=mh7(ck)+1;
      s{ck,RN-16}=s{ck,RN-16}(2:Td6(ck));
      Td6(ck)=Td6(ck)-1;
  end
  
  if RN==23  % Reaction 7: -> mh6
     mh6(ck)=mh6(ck)+1;
     s{ck,RN-16}=s{ck,RN-16}(2:Td7(ck));
     Td7(ck)=Td7(ck)-1; 
  end
  
  if RN==24 % Reaction 8: -> md
     md(ck)=md(ck)+1; 
     s{ck,RN-16}=s{ck,RN-16}(2:Td8(ck));
     Td8(ck)=Td8(ck)-1;
  end
  
  % For the selected reaction and cell, update Pk. 
  % (if not completing a delayed reaction)
  if RN<=16 && RN>=1
     r = rand;
     Pk(ck,RN) = Pk(ck,RN) + log(1/r);
  end
  
  if RN==1 % Reaction 1: mh1 -> ph1
     s{ck,RN} = s{ck,RN}(1:Td1(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nph1, inf];
     Td1(ck)=Td1(ck)+1;
  end
      
  if RN==2 % Reaction 2: mh7 -> ph7
     s{ck,RN} = s{ck,RN}(1:Td2(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nph7, inf];
     Td2(ck)=Td2(ck)+1;
  end 
      
  if RN==3 % Reaction 3: mh6 -> ph6
     s{ck,RN} = s{ck,RN}(1:Td3(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nph6, inf];
     Td3(ck) = Td3(ck) + 1;
  end
  
  if RN==4 % Reaction 4: md -> pd
     s{ck,RN} = s{ck,RN}(1:Td4(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, npd, inf]; 
     Td4(ck) = Td4(ck) + 1;
  end
  
  if RN==5 % Reaction 5: -> mh1 
     s{ck,RN} = s{ck,RN}(1:Td5(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nmh1, inf]; 
     Td5(ck) = Td5(ck) + 1;
  end
  
  if RN==6 % Reaction 6: -> mh7 
     s{ck,RN} = s{ck,RN}(1:Td6(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nmh7, inf];
     Td6(ck) = Td6(ck) + 1;
  end
  
  if RN==7  % Reaction 7: -> mh6
     s{ck,RN} = s{ck,RN} (1:Td7(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nmh6, inf];
     Td7(ck) = Td7(ck) + 1;
  end
  
  if RN==8 % Reaction 8: -> md 
     s{ck,RN} = s{ck,RN} (1:Td8(ck) - 1); 
     s{ck,RN} = [s{ck,RN}, nmd, inf];
     Td8(ck) = Td8(ck) + 1;
  end
  
  if RN==9 % Reaction 9: ph1 ->
      ph1(ck)=ph1(ck)-1;
  elseif RN==10 % Reaction 10: ph7 ->  
      ph7(ck)=ph7(ck)-1;
  elseif RN==11 % Reaction 11: ph6 ->
      ph6(ck)=ph6(ck)-1;
  elseif RN==12 %  Reaction 12: pd ->     
      pd(ck)=pd(ck)-1;
  elseif RN==13 % Reaction 13: mh1 ->  
      mh1(ck)=mh1(ck)-1;
  elseif RN==14 % Reaction 14: mh7 -> 
      mh7(ck)=mh7(ck)-1;
  elseif RN==15 % Reaction 15: mh6 ->
      mh6(ck)=mh6(ck)-1;
  elseif RN==16 % Reaction 16: md ->  
      md(ck)=md(ck)-1;
  end
  
  % Update Tk values for each reaction and each cell.
  for ckvalue=1:2
      for k=1:16
          Tk(ckvalue,k)=Tk(ckvalue,k)+a(ckvalue,k)*Delta;
      end
  end
    
   
    Time=[Time T]; % Store time, and mh1 and mh7 levels.
    mh1v_c1=[mh1v_c1 mh1(1)];
    mh1v_c2=[mh1v_c2 mh1(2)];
    
    
end % for i=1:maxi-1,

Data = [Time' mh1v_c1' mh1v_c2'];
%DataTable=table(Time', mh1v_c1', mh1v_c2','VariableNames', {'Time', 'mh1_cell1', 'mh7_cell1'});
%writetable(DataTable, strcat('Run_v4', num2str(run),'.xlsx'),'WriteVariableNames', true);

%To find the sync score, type the following into the command window:
%sync_score = corr(mh1v_c1,mh1v_c2,'Type','Pearson');

% figure
% plot(Time, mh1v_c1,'b')
% hold on
% plot(Time, mh1v_c2,'r')
% hold on
% 
% xlabel('Time')
% ylabel('#mRNA')
%saveas(gcf,strcat(['Run_v2', num2str(run)]),'jpg');
%close(gcf);

end