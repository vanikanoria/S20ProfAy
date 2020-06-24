function Y = vani_hybrid(x)
% hybrid model for a single cell system
%x: parameter set
% size of parameter set = 44 x 1

% set system specific parameters
num_states = 14;
num_reactions = 34;
cells=1; %rows*columns; %two in this case

% choose time discretization of output (w.r.t. days)
minutes=60;%600;
step=1; 
%step size at which data is stored (1 minute)
num_steps=minutes/step;
Tend=minutes; %more clear variable name;
% initialize output variables
T =1:step:minutes; T = T(:); % vector of time points in unit (minutes)
t=0;

% setting partition of the reactions (zero = deterministic, one = stochastic)
partition = zeros(num_reactions,1);
partition(1:8) = 1;
% 1 means stochastic and 0 means deterministic

R = get_R();
% net changes in species levels caused by firing of the reactions
% (1st dimension species, 2nd reaction)

psh1=x(1); psh6=x(2);psh7=x(3);psd=x(4);pdh1=x(5);pdh6=x(6);
pdh7=x(7);pdd=x(8);msh1=x(9);msh6=x(10);msh7=x(11);msd=x(12);
mdh1=x(13);mdh6=x(14);mdh7=x(15);mdd=x(16);pdh11=x(17);pdh16=x(18);
pdh17=x(19);pdh66=x(20);pdh76=x(21);pdh77=x(22);nmh1=x(23);
nmh7=x(24);nmd=x(25);nph1=x(26);nph6=x(27);nph7=x(28);npd=x(29);
dah1h1=x(30);ddh1h1=x(31);dah1h6=x(32);ddh1h6=x(33);dah1h7=x(34);
ddh1h7=x(35);dah6h6=x(36);ddh6h6=x(37);dah7h6=x(38);ddh7h6=x(39);
dah7h7=x(40);ddh7h7=x(41);critph11=x(42);critph76=x(43);critpd=x(44);

%Convert the time delay values to integers, because the deterministic
% simulation uses discrete time points:
nmh1_eps=round(nmh1/eps);nmh7_eps=round(nmh7/eps);
nmd_eps=round(nmd/eps);nph1_eps=round(nph1/eps);
nph6_eps=round(nph6/eps);nph7_eps=round(nph7/eps);
npd_eps=round(npd/eps);

%functions to calculate the transcription rates of mRNA
fd=@(ph11,ph76) msd/(1+(ph11/critph11)^2+(ph76/critph76)^2);
fh1=@(ph11,ph76,pd)msh1*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);
fh7=@(ph11,ph76,pd)msh7*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);

function a = get_propensities(y)
    a = zeros([34,1]); % Initialize the propensity levels for each reaction for cell.
    
    %The different molecule types are: [ph1;ph7;ph6,pd;mh1;mh7;mh6;md;ph11;ph76;ph17;ph16;ph77;ph66]
    ph1=y(1);ph7=y(2);ph6=y(3);pd=y(4);mh1=y(5);mh7=y(6);mh6=y(7);md=y(8);
    ph11=y(9);ph76=y(10);ph17=y(11);ph16=y(12);ph77=y(13);ph66=y(14);
    
    for ck=1
        a(1,ck) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
        a(2,ck) = psh7*mh7(ck);% // Reaction 2: mh7 -> ph7
        a(3,ck) = psh6*mh6(ck); % // Reaction 3: mh6 -> ph6
        a(4,ck) = psd*md(ck); % // Reaction 4: md -> pd
        a(5,ck)=fh1(ph11(ck),ph76(ck),pd(ck)); % // Reaction 5: -> mh1
        a(6,ck)=fh7(ph11(ck),ph76(ck),pd(ck)); % // Reaction 6: -> mh7
        a(7,ck) = psh6; % // Reaction 7: -> mh6
        a(8,ck)=fd(ph11(ck),ph76(ck));% // Reaction 8: -> md
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

% event function of stochastic reaction
function [value,isterminal,direction] = reaction_event(y)
    value      = y(end); %random number
    isterminal = 1;
    direction  = 1;
end

% perform stochastic reaction event
function y = perform_stochastic_reaction_event(y)
    a = get_propensities(y);                            % compute propensities of all reactions
    a = partition.*a;                           % take only stochastic reactions
    r = find(cumsum(a) >= sum(a)*rand,1);       % and choose a reaction (first reaction algorithm by Gillespie)
    %r : number of reaction that is carried out
    y = [max(y(1:end-1)+R(:,r)',0) log(rand)];	% update species levels accordingly and draw new random number
    %R = sparse(num_states, num_reactions);
    %so R(:,r) is the change in states due to reaction r
    %R(:,r)' is a row vector
end


    
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
    
    [v1,i1,d1] = reaction_event(y);

    value      = v1;
    isterminal = i1;
    direction  = d1;
end

    function ydot = find_next_stochastic_reaction(t,y)
    % compute right hand side
    
    % get propensities
    a = get_propensities(y);
    ydot = [R*((1-partition).*a); partition'*a];
    %partition'*a has dimension 1 x 1
    % it is the sum of propensities of all the stochastic reactions
    
    %partition = zeros(num_reactions,1);
    %R = sparse(num_states, num_reactions);
    %a = num_reactions x 1
    %ydot = [R*((1-partition).*a); partition'*a]; %the deterministic parts and then the stochastic parts
    %[R*((1-partition).*a) has dimension num_states x 1
    %partition'*a has dimension 1 x 1
    % it is the sum of propensities of all the stochastic reactions
end

Y=zeros(num_steps,num_states);
% set ode options
options = odeset('Events',@events,'NonNegative',1:num_states);
% set pseudo-initial species numbers
y0      = [zeros(num_states,1); log(rand)];
y = y0;     
step=2;
while t < Tend
    % find te until next stochastic reaction or integrate until the end of
    % the step
    [t,y,te,ye,ie] = ode45(@find_next_stochastic_reaction,[t T(step)],y,options);
    %     [Delta RN] = min([t(1,:),s{1,1}(1),s{1,2}(1),s{1,3}(1), s{1,4}(1), s{1,5}(1), s{1,6}(1), t(2,:),s{2,1}(1),s{2,2}(1),s{2,3}(1),s{2,4}(1),s{2,5}(1),s{2,6}(1)]);
    %
    %     tau_n = min(te, delayed reactions)
    % check for events
    %ie: index of triggered event function
    if ~isempty(ie)
        
        % only one event at a particular time can be handled
        if length(ie) > 1
            error('ODEsolution:Events','\nMore than one event detected.');
        else
            % update time
            t = te;
            
            % handle event
            y = perform_stochastic_reaction_event(ye);
        end
        
        % write output if necesarry
        if t == T(step)
            Y(step,:)  = y(1:num_states);
            step       = step+1;
            step
        end
        
        % no event
    else
        t  = t(end);    % update time
        y  = y(end,:);	% and state
        
        % write output
        Y(step,:)  = y(1:num_states);
        step       = step+1;
    end
end
end
