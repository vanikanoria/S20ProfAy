global a;
num_states = 14;
global num_reactions;
num_reactions = 34;
global num_cells; 
num_cells=2; 
x=create_param_set();

partition = ones(num_reactions,1);
partition([11,12,13,14,15,17,18,25,27,29,31,33])=0; % Dimer association and dissociation 
partition([19,20,21,22,23,24])=0; % Dimer degradation
partition([3,7,16,32])=0; % hes6 mRNA synthesis and decay
stochreactions=[1,2,4,5,6,8,9,10,26,28,30,34];

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



% set pseudo-initial species numbers
stocnum = size(stochreactions,2); %returns the length of dimension 2 of stochreactions
global logrand;

logrand=log(rand(stocnum,num_cells)); %here logrand has 12 rows and 2 columns
y0 = [ones(num_states,num_cells); logrand];
%y0      = [ones(num_states,1); log(rand)];
y = y0; 

t=0;
R = get_R();
options3 = odeset('NonNegative',1:num_states,'Events',@events);
tau_n=1;

[tsol,ysol,te,ye,ie] = ode45(@update_ode,[t t+tau_n],[y; logrand],options3);
ie
ye
te

 function ydot = update_ode(t,y)   
        a=get_propensities(y);
        ydot=[R*((1-partition).*a); a(stochreactions)]; %sum of stochastic propensities
 end
    
 function [value,isterminal,direction] = events(t,y)
    value      = y(num_states+1:end, num_cells); %which means the whole of the logrand function
    isterminal = ones(size(stochreactions,2),num_cells);
    direction  = ones(size(stochreactions,2),num_cells);
 end

 function a = get_propensities(y)
    num_reactions=34;
    num_cells=2;
    a = zeros([num_reactions,num_cells]); % Initialize the propensity levels for each reaction for cell.
    
    % The different molecule types are: 
    % [ph1;ph7;ph6,pd;mh1;mh7;mh6;md;ph11;ph76;ph17;ph16;ph77;ph66]
    ph1=y(1,:);ph7=y(2,:);ph6=y(3,:);pd=y(4,:);mh1=y(5,:);mh7=y(6,:);mh6=y(7,:);md=y(8,:);
    ph11=y(9,:);ph76=y(10,:);ph17=y(11,:);ph16=y(12,:);ph77=y(13,:);ph66=y(14,:);
    
     %In case of num_cells=2: we take pd levels of the other cell into
     %account in reactions number 5 and 6
    %ck_other = 3-ck; %if ck=1: ck_other=2; if ck=2: ck_other=1
    % need to be modified in case of num_cells ~= 2:
    
    for ck=num_cells
        a(1,ck) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
        a(2,ck) = psh7*mh7(ck);% // Reaction 2: mh7 -> ph7
        a(3,ck) = psh6*mh6(ck); % // Reaction 3: mh6 -> ph6
        a(4,ck) = psd*md(ck); % // Reaction 4: md -> pd
        a(5,ck)= fh1(ph11(ck),ph76(ck),pd(3-ck)); % // Reaction 5: -> mh1
        a(6,ck)= fh7(ph11(ck),ph76(ck),pd(3-ck)); % // Reaction 6: -> mh7
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