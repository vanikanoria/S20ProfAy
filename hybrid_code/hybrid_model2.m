function [Y] = hybrid_model2(x) % x: parameter set of size 44
% Hybrid model for a single cell system
% Set system specific parameters
num_states = 14;
num_reactions = 34;
num_cells=1; %rows*columns; %two in this case

global s;
s = cell(8,1); % Time for each delayed reaction is scheduled in these vectors.
each_cell={inf};
s(:,1)=each_cell;
global a;
global Td;
Td = ones(8,1); % Number of scheduled delayed reactions by reaction number.

% Choose time discretization of output (w.r.t. minutes)
minutes=100;%600;
step=1; % Step size at which data is stored (1 minute)
num_steps=minutes/step;
Tend=minutes; % More clear variable name;
% Initialize output variables
T=1:step:minutes; T = T(:); % Vector of time points in unit (minutes)
t=0;

% Setting partition of the reactions (0 = deterministic, 1 = stochastic)
partition = ones(num_reactions,1);
partition([11,12,13,14,15,17,18,25,27,29,31,33])=0; % Dimer association and dissociation 
partition([19,20,21,22,23,24])=0; % Dimer degradation
partition([3,7,16,32])=0; % hes6 mRNA synthesis and decay

delay = zeros(num_reactions,1);
delay(1:6) = 1;
delay(8) = 1;
delayed_set = delay; % 1 means delayed and 0 means not delayed

R = get_R();
% Net changes in species levels caused by firing of the reactions (1st dimension species, 2nd reaction)

psh1=x(1); psh6=x(2);psh7=x(3);psd=x(4);pdh1=x(5);pdh6=x(6);
pdh7=x(7);pdd=x(8);msh1=x(9);msh6=x(10);msh7=x(11);msd=x(12);
mdh1=x(13);mdh6=x(14);mdh7=x(15);mdd=x(16);pdh11=x(17);pdh16=x(18);
pdh17=x(19);pdh66=x(20);pdh76=x(21);pdh77=x(22);nmh1=x(23);
nmh7=x(24);nmd=x(25);nph1=x(26);nph6=x(27);nph7=x(28);npd=x(29);
dah1h1=x(30);ddh1h1=x(31);dah1h6=x(32);ddh1h6=x(33);dah1h7=x(34);
ddh1h7=x(35);dah6h6=x(36);ddh6h6=x(37);dah7h6=x(38);ddh7h6=x(39);
dah7h7=x(40);ddh7h7=x(41);critph11=x(42);critph76=x(43);critpd=x(44);

%functions to calculate the transcription rates of mRNA
fd=@(ph11,ph76) msd/(1+(ph11/critph11)^2+(ph76/critph76)^2);
fh1=@(ph11,ph76,pd)msh1*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);
fh7=@(ph11,ph76,pd)msh7*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);

function a = get_propensities(y)
    a = zeros([num_reactions,1]); % Initialize the propensity levels for each reaction for cell.
    
    % The different molecule types are: 
    % [ph1;ph7;ph6,pd;mh1;mh7;mh6;md;ph11;ph76;ph17;ph16;ph77;ph66]
    ph1=y(1);ph7=y(2);ph6=y(3);pd=y(4);mh1=y(5);mh7=y(6);mh6=y(7);md=y(8);
    ph11=y(9);ph76=y(10);ph17=y(11);ph16=y(12);ph77=y(13);ph66=y(14);
    
    for ck=num_cells
        a(1,ck) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
        a(2,ck) = psh7*mh7(ck);% // Reaction 2: mh7 -> ph7
        a(3,ck) = psh6*mh6(ck); % // Reaction 3: mh6 -> ph6
        a(4,ck) = psd*md(ck); % // Reaction 4: md -> pd
        a(5,ck)= fh1(ph11(ck),ph76(ck),pd(ck)); % // Reaction 5: -> mh1
        a(6,ck)= fh7(ph11(ck),ph76(ck),pd(ck)); % // Reaction 6: -> mh7
        a(7,ck) = msh6; % // Reaction 7: -> mh6
        a(8,ck)= fd(ph11(ck),ph76(ck));% // Reaction 8: -> md
        a(9,ck) = pdh1*ph1(ck); % Reaction 02: ph1 ->
        a(10,ck) = pdh7*ph7(ck);% // Reaction 10: ph7 ->
        a(11,ck) = dah7h7*ph7(ck)*(ph7(ck)-1)/2;% // Reaction 11: ph7+ph7 -> ph77
        a(12,ck) = ddh7h7*ph77(ck); % // Reaction 12: ph77 -> ph7+ph7
        a(13,ck) = dah7h6*ph7(ck)*ph6(ck);% // Reaction 6: ph7+ph6 -> ph76
        a(14,ck) = ddh7h6*ph76(ck); % // Reaction 14: ph76 -> ph7+ph6
        a(15,ck) = dah1h1*ph1(ck)*(ph1(ck)-1)/2; % Reaction 03: ph1+ph1 -> ph11
        a(16,ck) = pdh6*ph6(ck); % // Reaction 16: ph6 ->
        a(17,ck) = dah6h6*ph6(ck)*(ph6(ck)-1)/2; % // Reaction 17: ph6+ph6 -> ph66
        a(18,ck) = ddh6h6*ph66(ck); % // Reaction 18: ph66 -> ph6+ph6
        a(19,ck) = pdh1*ph11(ck); % // Reaction 19: ph11 ->
        a(20,ck) = pdh1*ph17(ck); % // Reaction 20: ph17 ->
        a(21,ck) = pdh1*ph16(ck); % // Reaction 21: ph16 ->
        a(22,ck) = pdh1*ph77(ck); % // Reaction 22: ph77 ->
        a(23,ck) = pdh1*ph76(ck); % // Reaction 23: ph76 ->
        a(24,ck) = pdh1*ph66(ck); % // Reaction 24: ph66 ->
        a(25,ck) = ddh1h1*ph11(ck); %  Reaction 04: ph11 -> ph1+ph1
        a(26,ck) = pdd*pd(ck); % // Reaction 26: pd ->
        a(27,ck) = dah1h7*ph1(ck)*ph7(ck); % Reaction 05: ph1+ph7 -> ph17
        a(28,ck) = mdh1*mh1(ck); % // Reaction 28: mh1 ->
        a(29,ck) = ddh1h7*ph17(ck); % Reaction 06: ph17 -> ph1+ph7
        a(30,ck) = mdh7*mh7(ck); % // Reaction 30: mh7 ->
        a(31,ck) = dah1h6*ph1(ck)*ph6(ck); % // Reaction 07: ph1+ph6 -> ph16
        a(32,ck) = mdh6*mh6(ck); %) // Reaction 32: mh6 ->
        a(33,ck) = ddh1h6*ph16(ck);% // Reaction 08: ph16 -> ph1+ph6
        a(34,ck) = mdd*md(ck); %// Reaction 34: md ->
    end
end

% events function for ode45:
% The output arguments value, isterminal, and direction are vectors whose 
% ith element corresponds to the ith event:
% value(i) is a mathematical expression describing the ith event. 
% An event occurs when value(i) is equal to zero. 
% isterminal(i) = 1 if the integration is to terminate when the ith event occurs.
% Otherwise, it is 0.
% direction(i) = 0 if all zeros are to be located (the default).
% A value of +1 locates only zeros where the event function is increasing, 
% and -1 locates only zeros where the event function is decreasing.
function [value,isterminal,direction] = events(t,y)
    value      = y(end); %log(rand);
    isterminal = 1;
    direction  = 1;
end
    
function s = start_delayed_reaction(r,s)
    
        if r==1 % Reaction 1: mh1 -> ph1
            s{r} = s{r} (1:Td(1) - 1);
            s{r} = [s{r}, nph1, inf];
            Td(1)=Td(1)+1;
            
        elseif r==2 % Reaction 2: mh7 -> ph7
            s{r} = s{r} (1:Td(2) - 1);
            s{r} = [s{r}, nph7, inf];
            Td(2)=Td(2)+1;
            
        elseif r==3 % // Reaction 3: mh6 -> ph6
            s{r} = s{r} (1:Td(3) - 1);
            s{r} = [s{r}, nph6, inf];
            Td(3)=Td(3) + 1;
            
        elseif r==4 % // Reaction 4: md -> pd
            s{r} = s{r} (1:Td(4) - 1);
            s{r} = [s{r}, npd, inf];
            Td(4) = Td(4) + 1;
            
        elseif r==5 % // Reaction 5: -> mh1
            s{r} = s{r} (1:Td(5) - 1);
            s{r} = [s{r}, nmh1, inf];
            Td(5) = Td(5) + 1;
        
        elseif r==6 % Reaction 6: -> mh7
            s{r} = s{r} (1:Td(6) - 1);
            s{r} = [s{r}, nmh7, inf];
            Td(6) = Td(6) + 1;
        
        elseif r==8 % // Reaction 8: -> md
            s{r} = s{r} (1:Td(8) - 1);
            s{r} = [s{r}, nmd, inf];
            Td(8) = Td(8) + 1;
        end
    end
    
    % perform stochastic reaction event
    function [y,s] = perform_stochastic_reaction(y,s)
        a = get_propensities(y);    % compute propensities of all reactions
        stochastic_a = partition.*a;  % take only stochastic reactions
    % Choose a reaction (first reaction algorithm by Gillespie)
    r = find(cumsum(stochastic_a) >= sum(stochastic_a)*rand,1); 
    % r : number of reaction that is carried out
    if delayed_set(r)==1
        s = start_delayed_reaction(r,s);
    else  %if it is not delayed...
   % update species levels accordingly and draw new random number
        y = max(y+R(:,r),0); 	
    end 
    end
    
    function [y,s] = end_delayed_reaction(r,s,y)
        s{r}=s{r}(2:Td(r));
        Td(r) = Td(r)-1;
        y = [max(y+R(:,r),0)] ;
    end
    
    function ydot = update_ode(t,y)   
        a=get_propensities(y);
        ydot=[R*((1-partition).*a); partition'*a]; %sum of stochastic propensities
    end
     
Y=zeros(num_steps, num_states);
% set ode options
options1 = odeset('Events',@events);
options2 = odeset('NonNegative',1:num_states);
options3 = odeset('NonNegative',1:num_states,'Events',@events);
% set pseudo-initial species numbers
y0      = [ones(num_states,1)];
y = y0; 
logrand=log(rand);
step=2;
while t < Tend
    % a = get_propensities(y);
    
    [tau_n, reaction] = min([s{1}(1),s{2}(1),s{3}(1),s{4}(1),s{5}(1),s{6}(1),s{8}(1)]);
    if tau_n==Inf
        tau_n=1;
    end

    [tsol,ysol,te,ye,ie] = ode45(@update_ode,[t t+tau_n],[y; logrand],options3);
    
    if ~isempty(ie)
       y=ye(1:num_states)';
       logrand=log(rand);
       tau_n=tsol(end)-t;           
       t=te;
       [y,s] = perform_stochastic_reaction(y,s);
%            if length(ie) > 1
%             error('ODEsolution:Events','\nMore than one event detected.');
%            end
    else
        logrand = ysol(end,num_states+1);
        y=ysol(end,1:num_states)';
        t=tsol(end);
        if reaction==7 
            reaction=8;
        end
        [y,s] = end_delayed_reaction(reaction,s,y);
    end
    
    %updating delays
    for rn=[1:6,8]
           s{rn}=s{rn}-tau_n;
    end      
        % write output if necesarry
        if t/T(step)<1.01&&t/T(step)>0.99
            Y(step,:)  = y(1:num_states);
            step       = step+1;
            disp(step);
        end
        
       if step > num_steps
            break;
       end
end
end