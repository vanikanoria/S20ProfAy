function mh1 = deterministic_model(param_set)

% Return mh1 levels over the timesteps for each cell
% Indicates the period and amplitude and hence if the model works with states 
% and parameters having values that are biologically possible because if the
% model attains biologically unfeasible values, it returns period and amplitude as 0
% thus indicating that the model doesn't work

global time_steps;
minutes=600; % 1200;
global eps; % Time step to be used for Euler's method, default is 0.01
eps=0.01;
time_steps = (minutes / eps); % Set the amount of time steps to be used in the simulation

% max_prop = Inf; % maximum threshold for propensity functions, default is INFINITY

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

% Runs the deterministic simulation of the model.
% For each time step:
% 1) Iterate through every cell (2 cells in this case) and update the concentrations of proteins and mRNA.
%    These concentration values are obtained by solving the differential equations for that time step, using Euler's method.
% 2) Check that the concentrations do not become negative -- a negative amount of protein is not biologically sensible
% 3) Check that the propensity functions do not go above the set threshold -- if one was specified.

% Convert the time delay values to integers, because the deterministic simulation uses discrete time points:
    
    nmh1_eps=round(nmh1/eps);
    nmh7_eps=round(nmh7/eps);
    nmd_eps=round(nmd/eps);
    nph1_eps=round(nph1/eps);
    nph6_eps=round(nph6/eps);
    nph7_eps=round(nph7/eps);
    npd_eps=round(npd/eps);
     
    cells=2; % rows*columns; %two in this case
    
    a = zeros([cells,34]); % Initialize the propensity levels for each reaction for both cells.
    % They will be checked later in the model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Number of timesteps = time_steps = time_steps
    % So we need an array of time_steps values for each state and each cell
    
    % Initializing the states: each state = array with cell no. and time at
    % which the state is calculated
    
    ph1=zeros(cells,time_steps); 
    ph7=zeros(cells,time_steps);
    ph6=zeros(cells,time_steps);
    pd=zeros(cells,time_steps);
    mh1=zeros(cells,time_steps);
    mh7=zeros(cells,time_steps);
    mh6=zeros(cells,time_steps);
    md=zeros(cells,time_steps); 
    ph11=zeros(cells,time_steps);
    ph76=zeros(cells,time_steps);
    ph17=zeros(cells,time_steps);
    ph16=zeros(cells,time_steps);
	ph77=zeros(cells,time_steps);
	ph66=zeros(cells,time_steps);

     % Functions to calculate the transcription rates of mRNA
    fd=@(ph11,ph76) msd/(1+(ph11/critph11)^2+(ph76/critph76)^2);
    fh1=@(ph11,ph76,pd) msh1*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2+(ph76/critph76)^2);
    fh7=@(ph11,ph76,pd) msh7*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2+(ph76/critph76)^2);
    
    for n = 2:1:time_steps  %% n is an integer
        for i = 1:1:cells            
            % Protein synthesis
            if n>nph1_eps
                prod_ph1 = psh1*mh1(i,n-nph1_eps);
            else
                prod_ph1 = 0;    
            end
            
            ph1(i,n) = ph1(i,n-1) + eps * (prod_ph1 - pdh1*ph1(i,n-1)-2*dah1h1*ph1(i,n-1)*ph1(i,n - 1)+2*ddh1h1*ph11(i,n - 1)-dah1h7*ph1(i,n - 1)*ph7(i,n - 1)+ddh1h7*ph17(i,n - 1)-dah1h6*ph1(i,n - 1)*ph6(i,n - 1)+ddh1h6*ph16(i,n - 1));
            
            if n > nph7_eps
                prod_ph7 = psh7*mh7(i,n - nph7_eps);
            else
                prod_ph7 = 0;
            end
            
            ph7(i,n) = ph7(i,n - 1) + eps * (prod_ph7 - pdh7*ph7(i,n - 1)-2*dah7h7*ph7(i,n - 1)*ph7(i,n - 1)+2*ddh7h7*ph77(i,n - 1)-dah1h7*ph1(i,n - 1)*ph7(i,n - 1)+ddh1h7*ph17(i,n - 1)-dah7h6*ph7(i,n - 1)*ph6(i,n - 1)+ddh7h6*ph76(i,n - 1)); 
            
            
            if n > nph6_eps
                prod_ph6 = psh6*mh6(i,n - nph6_eps);
            else
                prod_ph6 = 0;
            end
            
            ph6(i,n) = ph6(i,n-1) + eps * (prod_ph6 -pdh6*ph6(i,n - 1)-2*dah6h6*ph6(i,n - 1)*ph6(i,n - 1)+2*ddh6h6*ph66(i,n - 1)-dah1h6*ph1(i,n - 1)*ph6(i,n - 1)+ddh1h6*ph16(i,n - 1)-dah7h6*ph7(i,n - 1)*ph6(i,n - 1)+ddh7h6*ph76(i,n - 1));
               
            
            % Break out of the nested for-loops of any of them if any of the
            % checking for protein levels becoming negative
            if (ph1(i,n) < 0 || ph7(i,n) < 0 || ph6(i,n) < 0) 
                disp(ph1(i,n),ph7(i,n),ph6(i,n));
                mh1=0;
                return;
            end
            
            % Dimer proteins
            ph11(i,n) = ph11(i,n - 1) + eps * (dah1h1*ph1(i,n - 1)*ph1(i,n - 1)-ddh1h1*ph11(i,n - 1)-pdh11*ph11(i,n - 1));
            ph17(i,n) = ph17(i,n - 1) + eps * (dah1h7*ph1(i,n - 1)*ph7(i,n - 1)-ddh1h7*ph17(i,n - 1)-pdh17*ph17(i,n - 1));
            ph16(i,n) = ph16(i,n - 1) + eps * (dah1h6*ph1(i,n - 1)*ph6(i,n - 1)-ddh1h6*ph16(i,n - 1)-pdh16*ph16(i,n - 1));
            ph77(i,n) = ph77(i,n - 1) + eps * (dah7h7*ph7(i,n - 1)*ph7(i,n - 1)-ddh7h7*ph77(i,n - 1)-pdh77*ph77(i,n - 1));
            ph76(i,n) = ph76(i,n - 1) + eps * (dah7h6*ph7(i,n - 1)*ph6(i,n - 1)-ddh7h6*ph76(i,n - 1)-pdh76*ph76(i,n - 1));
            ph66(i,n) = ph66(i,n - 1) + eps * (dah6h6*ph6(i,n - 1)*ph6(i,n - 1)-ddh6h6*ph66(i,n - 1)-pdh66*ph66(i,n - 1));

            % Delta Protein
            if n > npd_eps
                prod_pd = psd*md(i,n-npd_eps);
            else
                prod_pd = 0;
            end
            
            pd(i,n) = pd(i,n-1) + eps*(prod_pd - pdd*pd(i,n - 1));
            
            
            if (ph11(i,n) < 0 || ph17(i,n) < 0 || ph16(i,n) < 0 || ph77(i,n) < 0 || ph76(i,n) < 0 || ph66(i,n) < 0 || pd(i,n) < 0)
                mh1=0;
                return;
            end
            
            if i==1
                other=2;
            else
                other=1;
            end
                
            % mRNA Synthesis
            if (n > nmh1_eps)
                prod_mh1 = fh1(ph11(i,n - nmh1_eps),ph76(i,n - nmh1_eps), pd(other,n - nmh1_eps));
            else
                prod_mh1 = 0;
            end
            mh1(i,n) = mh1(i,n-1) + eps * (prod_mh1 - mdh1*mh1(i,n - 1)); 
            
            if (n > nmh7_eps)
                prod_mh7 = fh7(ph11(i,n - nmh7_eps),ph76(i,n - nmh7_eps),pd(other,n - nmh7_eps));
            else
                prod_mh7 = 0;
            end
            mh7(i,n) = mh7(i,n-1) + eps * (prod_mh7-mdh7*mh7(i,n - 1));
            
			mh6(i,n) = mh6(i,n-1) + eps * (msh6-mdh6*mh6(i,n - 1)); 
            
            if n > nmd_eps
                prod_md = fd(ph11(i,n - nmd_eps),ph76(i,n - nmd_eps));
            else
                prod_md = 0;
            end
			md(i,n) = md(i,n-1) + eps * (prod_md - mdd*md(i,n - 1));
            
            % Checking for negative mRNA levels
            if (mh1(i,n) < 0 || mh7(i,n) < 0 || mh6(i,n) < 0 || md(i,n) < 0)
                mh1=0;
                return;
            end
            
            % Checking for negative and infinite propensities
            % C code:
%             if (max_prop != INFINITY && !checkPropensities(g, r, n, max_prop)) {
%                 return false;
%             }

            % 1-8: delayed reactions
            % 9-34: non-delayed reactions
            for ck=1:2
                a(ck,1) = psh1*mh1(ck); % Reaction 01: mh1 -> ph1
                a(ck,2) = psh7*mh7(ck);% // Reaction 09: mh7 -> ph7
                a(ck,3) = psh6*mh6(ck); % // Reaction 15: mh6 -> ph6
                a(ck,4) = psd*md(ck); % // Reaction 25: md -> pd
                a(ck,5)= fh1(ph11(ck),ph76(ck),pd(3-ck)); % // Reaction 27: -> mh1
                a(ck,6)= fh7(ph11(ck),ph76(ck),pd(3-ck)); % // Reaction 29: -> mh7
                a(ck,7) = psh6; % // Reaction 31: -> mh6
                a(ck,8)= fd(ph11(ck),ph76(ck));% // Reaction 33: -> md
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
                a(ck,19) = pdh11*ph11(ck); % // Reaction 19: ph11 ->
                a(ck,20) = pdh17*ph17(ck); % // Reaction 20: ph17 ->
                a(ck,21) = pdh16*ph16(ck); % // Reaction 21: ph16 ->
                a(ck,22) = pdh77*ph77(ck); % // Reaction 22: ph77 ->
                a(ck,23) = pdh76*ph76(ck); % // Reaction 23: ph76 ->
                a(ck,24) = pdh66*ph66(ck); % // Reaction 24: ph66 ->
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

%             if ~(a<Cutoff) && max_prop ~= inf %to satisfy this any one of the elements of a must satisfy a>=cutoff
%                 disp('x')
%                 mh1=0;
%                 return; %will return mh1 as 0
%             end
            
        end
        
    end  
end
