function hybrid_model(param_set)
% set system specific parameters
num_species;
num_reactants;
num_reactions;
minutes=600;
eps=0.01;
time_steps=minutes/eps;

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

%Convert the time delay values to integers, because the deterministic
% simulation uses discrete time points:
nmh1_eps=round(nmh1/eps);
nmh7_eps=round(nmh7/eps);
nmd_eps=round(nmd/eps);
nph1_eps=round(nph1/eps);
nph6_eps=round(nph6/eps);
nph7_eps=round(nph7/eps);
npd_eps=round(npd/eps);

cells=2; %rows*columns; %two in this case

%functions to calculate the transcription rates of mRNA
fd=@(ph11,ph76) msd/(1+(ph11/critph11)^2+(ph76/critph76)^2);
fh1=@(ph11,ph76,pd)msh1*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);
fh7=@(ph11,ph76,pd)msh7*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);

a = zeros([2,34]); % Initialize the propensity levels for each reaction for both cells.
function c = initial_propensities(a)
    % 1-8: delayed reactions
            % 9-34: non-delayed reactions
            for ck=1:2
                a(ck,1) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
                a(ck,2) = psh7*mh7(ck);% // Reaction 09: mh7 -> ph7
                a(ck,3) = psh6*mh6(ck); % // Reaction 15: mh6 -> ph6
                a(ck,4) = psd*md(ck); % // Reaction 25: md -> pd
                a(ck,5)=fh1(ph11(ck),ph76(ck),pd(3-ck)); % // Reaction 27: -> mh1
                a(ck,6)=fh7(ph11(ck),ph76(ck),pd(3-ck)); % // Reaction 29: -> mh7
                a(ck,7) = psh6; % // Reaction 31: -> mh6
                a(ck,8)=fd(ph11(ck),ph76(ck));% // Reaction 33: -> md
                a(ck,9) = pdh1*ph1(ck); % Reaction 02: ph1 ->
                a(ck,10) = pdh7*ph7(ck);% // Reaction 10: ph7 ->
                a(ck,11) = dah7h7*ph7(ck)*(ph7(ck)-1)/2;% // Reaction 11: ph7+ph7 -> ph77
                a(ck,12) = ddh7h7*ph77(ck); % // Reaction 12: ph77 -> ph7+ph7
                a(ck,13) = dah7h6*ph7(ck)*ph6(ck);% // Reaction 6: ph7+ph6 -> ph76
                a(ck,14) = ddh7h6*ph76(ck); % // Reaction 14: ph76 -> ph7+ph6
                a(ck,15) = dah1h1*ph1(ck)*(ph1(ck)-1)/2; % Reaction 03: ph1+ph1 -> ph11
                a(ck,16) = pdh6*ph6(ck); % // Reaction 16: ph6 ->
                a(ck,17) = dah6h6*ph6(ck)*(ph6(ck)-1)/2; % // Reaction 17: ph6+ph6 -> ph66
                a(ck,18) = ddh6h6*ph66(ck); % // Reaction 18: ph66 -> ph6+ph6
                a(ck,19) = pdh1*ph11(ck); % // Reaction 19: ph11 ->
                a(ck,20) = pdh1*ph17(ck); % // Reaction 20: ph17 ->
                a(ck,21) = pdh1*ph16(ck); % // Reaction 21: ph16 ->
                a(ck,22) = pdh1*ph77(ck); % // Reaction 22: ph77 ->
                a(ck,23) = pdh1*ph76(ck); % // Reaction 23: ph76 ->
                a(ck,24) = pdh1*ph66(ck); % // Reaction 24: ph66 ->
                a(ck,25) = ddh1h1*ph11(ck); %  Reaction 04: ph11 -> ph1+ph1
                a(ck,26) = pdd*pd(ck); % // Reaction 26: pd ->
                a(ck,27) = dah1h7*ph1(ck)*ph7(ck); % Reaction 05: ph1+ph7 -> ph17
                a(ck,28) = mdh1*mh1(ck); % // Reaction 28: mh1 ->
                a(ck,29) = ddh1h7*ph17(ck); % Reaction 06: ph17 -> ph1+ph7
                a(ck,30) = mdh7*mh7(ck); % // Reaction 30: mh7 ->
                a(ck,31) = dah1h6*ph1(ck)*ph6(ck); % // Reaction 07: ph1+ph6 -> ph16
                a(ck,32) = mdh6*mh6(ck); %) // Reaction 32: mh6 ->
                a(ck,33) = ddh1h6*ph16(ck);% // Reaction 08: ph16 -> ph1+ph6
                a(ck,34) = mdd*md(ck); %// Reaction 34: md ->
            end
end

function [] = update_propensities(y)
    
    % initialize propensities with reaction constants
    a = param_set; %reaction_constants
    
    
    for i = 1:num_reactions % first two reactions are independent from y
        a(i) = a(i)*prod(y(V(i-2).reactants)); % prod function gives the product of array elements
    end
end
    
% partition of the reactions (zero = deterministic, one = stochastic)
partition = zeros(num_reactions,1);
a_partition       = 0;  % all reactions with propensities less than 'asd_dflt' will be treated stochastically, otherwise deterministically
y_partition       = 0;	% all reactions with any reactant level less than 'ysd_dflt' will be treated stochastically, otherwise deterministically

% 1 means stochastic and 0 means deterministic
function [] = set_partition(y)

    update_propensities(y);
    for r = 1:num_reactions
        if a(r) < a_partition
            partition(r) = 1;
        else
            Reactants = V(r-2).reactants;
            for s = 1:length(Reactants)
                if R(Reactants(s),r) && y(Reactants(s)) < y_partition
                    partition(r) = 1;
                    break
                end
            end
        end
    end
end


        
function [value,isterminal,direction] = events(t,y)
    
    [v1,i1,d1] = M.reaction_event(y);
    [v2,i2,d2] = P.patient_monitoring(t);
    
    value      = [v1 v2];
    isterminal = [i1 i2];
    direction  = [d1 d2];
end

% event function of stochastic reaction
function [value,isterminal,direction] = reaction_event(y)
    value      = y(end);
    isterminal = 1;
    direction  = 1;
end


% perform stochastic reaction event
function y = perform_reaction_event(y)
    update_propensities(y);                            % compute propensities of all reactions
    a = partition.*a;                           % take only stochastic reaction
    r = find(cumsum(a) >= sum(a)*rand,1);       % and choose a reaction (first reaction algorithm by Gillespie)
    y = [max(y(1:end-1)+R(:,r)',0) log(rand)];	% update species levels accordingly and draw new random number
    
    % update partitioning of reactions
    set_partition(y);
end

options = odeset('Events',@events,'NonNegative',1:M.NrSpecies);

% choose time discretization of output (w.r.t. days)
dT;

% compute and return initial conditions of the system:
% initial state is approximated by a pre-run of length T,
% if the flag 'wt_only' is set to true the inital state will contain no mutants
function y0 = set_initial_conditions(T,wt_only)
    
    % pseudo-initial loads (obtained by pure deterministic simulations, used as initial conditions for pre-run)
    Initial_T  = 2.0e11 + 0.02*2.0e11 - 1.0e6 - 1.0e5;
    Initial_T1 = 1.0e5;
    Initial_T2 = 1.0e6;
    Initial_TL = 1.0e5;
    Initial_M  = 2.0e9 + 0.02*2.0e9 - 1.0e6 - 1.0e5;
    Initial_M1 = 1.0e3;
    Initial_M2 = 1.0e4;
    Initial_V  = 17.5e6;
    
    % set pseudo-initial species numbers
    y0      = [zeros(NrSpecies,1); 0];
    y0(1:8) = [Initial_T;Initial_M;Initial_T1;Initial_T2;Initial_TL;Initial_M1;Initial_M2;Initial_V];
    
    % compute reaction constants for the case when no drug is given
    if wt_only
        SelectionMatrixBackUp = SelectionMatrix;
        SelectionMatrix(1:AvailableTargets,2:end) = 0;
    end
    set_drugcombination(0);
    
    % initially set all reaction to deterministic
    partition = zeros(NrReactions,1);
    set_asd(0); set_ysd(0);
    
    % integrator options: all species numbers have non-negative values!
    options = odeset('NonNegative',1:NrSpecies);

    % compute initial species numbers
    [t,y] = ode15s(@RHS,[0 T],y0,options);
    
    % set initial species vector
    y0 = round(y(end,:));
    if wt_only
        y0(8:end-2) = 0;
        SelectionMatrix = SelectionMatrixBackUp;
        clear SelectionMatrixBackUp;
    end
    
    % draw first random number for stochastic reaction event
    y0(end) = log(rand);
    
    % choose first drug combination
    set_drugcombination(1);
    
    % set default partitioning thresholds and update partitioning
    set_asd(asd_dflt); set_ysd(ysd_dflt);
    set_partition(y0);
end


% functions to set partitioning thresholds
function [] = set_asd(threshold)
    asd = threshold;
end
function [] = set_ysd(threshold)
    ysd = threshold;
end


% set default partitioning thresholds
asd_dflt = 20; ysd_dflt = 20;
set_asd(asd_dflt); set_ysd(ysd_dflt);





% set event function of stochastic reactions and patient monitoring

% set ode options

while t < tend
    
    % solve system
    [t,y,te,ye,ie] = ode45(M.RHS,[t T(step)],y,options);
    
    % check for events
    if ~isempty(ie)
        
        % only one event at a particular time can be handled
        if length(ie) > 1
            error('ODEsolution:Events','\nMore than one event detected.');
        else
            
            % update time
        	t = te;
            %fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.2E days:\t%3.2E seconds',t,toc);
            
            % handle event
            switch ie
                
                % stochastic reaction event
                case 1
                    y = M.perform_reaction_event(ye);
                
                % patient check up
                case 2
                    
                    % get actual state, and
                    y = ye;
                    
                    % perform check up, i.e., measure the viral load and
                    % decide if treatment failed
                    if P.perform_checkup(y)
                        
                        % switch drug combination if not done before
                        if switch_to_nextcombo
                            M.set_drugcombination(2);       % activate next drug combination
                            switch_to_nextcombo = false;	% no further switches considered
                        
                        % otherwise: failure
                        else
                            
                            % get actual output values and stop simulation
                            T(step)    = t;
                            Y(step,:)  = y(1:M.NrSpecies);
                            VL(step)   = M.get_viral_load(y);
                            RC(step,:) = M.get_reproductive_capacities(y);
                            break
                        end
                    end
                
                otherwise
                    error('ODEsolution:Events','\nUnknown event detected.');
            end
            
            % write output if necesarry
            if t == T(step)
                Y(step,:)  = y(1:M.NrSpecies);
                VL(step)   = M.get_viral_load(y);
                RC(step,:) = M.get_reproductive_capacities(y);
                step       = step+1;
            end
        end
    
    % no event
    else
        t  = t(end);    % update time
        y  = y(end,:);	% and state
        fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.2E days:\t%3.2E seconds',t,toc);
        
        % write output
        Y(step,:)  = y(1:M.NrSpecies);
        VL(step)   = M.get_viral_load(y);
        RC(step,:) = M.get_reproductive_capacities(y);
        step       = step+1;
    end
end


end

