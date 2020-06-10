% deterministic simulator for the zebrafish segmentation clock
% MATLAB version of the C code by Ahmet Ay, Jack Holland, and Adriana Sperlea
 
%  The program simulates the behavior of the zebrafish segmentation clock for two cells,
%  through a mathematical model using delay differential equations.
%  To solve the equations, the program uses Euler's method, with an adjustable stepsize.
%  The program tries to find parameter sets which replicate the behavior of the system in wild type
%  and several mutants.


% Iterate through every parameter set

for p = 0:CHUNK_SIZE:PARS

    if PARS-p>CHUNKSIZE
        STEP = CHUNKSIZE;
    else
        STEP = PARS - p;
    end
    
    for j = 0:1:STEP 
            
            fprintf("Simulating parameter set %i: Step %i",p,j); % Used for creating output file names specific to the paramater
			
			temp_rate = rateValues(j); %// Temporary rate structure used to alter the protein synthesis rates in order to create mutants
			
            % Booleans to see if we met mutant conditions on each iteration
            wt = false; her1_mutant = false; her7_mutant = false; her6_mutant = false;
            her76_mutant = false; delta_mutant = false;
    
%             For the wild type and every mutant, perform the following steps:
%              1) Adjust the appropriate protein synthesis rates to create mutants if necessary
%              2) Run the simulation
%              3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
%              4) Otherwise, test oscillation features
%              5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions

  wt = run_mutant(&gene, t_steps, eps, temp_rate, of_wt, true, max_prop, x, y);
    
    if ~wt
         continue ;
    end
    
    
 %Inline functions for testing mutant conditions -- this is where you may change the condition ranges.

	% Wild type
    % In C: wt = fwildtype(of_wt.peaktotrough1, of_wt.peaktotrough2);
    wt= peaktotrough2 >= 1.5 && peaktotrough >= 1.5 && (peaktotrough2 / peaktotrough)<=1.5;
    if ~wt
        continue;
    end
    
    % Delta mutant
    temp_psd=0.0; %%need to incorporate this into run_mutant
    delta_mutant_model_works = run_mutant(temp_rate, of_delta, false);
    
    if ~delta_mutant_model_works
        continue;
    end
    
    delta_mutant = ((dperiod / wperiod) > 1.04) && ((dperiod / wperiod) < 1.30);
    if ~delta_mutant
        continue;
    end
    
    % >>>>>>here change psd back to its original value
    % in C: temp_rate.psd = rateValues[i].psd; 
    
	% Her1 mutant
    temp_psh1 = 0.0;
    
    
    her1_mutant = ((h1period / wperiod) > 0.97) && ((h1period / wperiod) < 1.03);
    if ~her1_mutant
        continue;
    end
    
    % Need to print?
%     if (toPrint) {
%         printForPlotting(mutants[3] + "/run0.txt", &gene, t_steps, eps);
%     }
    
    % Her7 mutant
    temp_rate.psh7 = 0.0;
    her7_mutant_model_works = run_mutant(&gene, t_steps, eps, temp_rate, of_her7, true, max_prop, x, y);
    if ~her7_mutant_model_works
        continue;
    end
    her7_mutant = ((h7period / wperiod) > 0.97) && ((h7period / wperiod) < 1.03);
    if ~her7_mutant
        continue;
    end
    %  temp_rate.psh7 = rateValues[i].psh7;
    
	% Her6 mutant
    temp_psh6 = 0.0; 
    her6_mutant_model_works = run_mutant(temp_rate, of_her13, false);
    if ~delta_mutant_model_works
        continue;
    end
    
    her6_mutant = ((h13period / wperiod) > 1.03) && ((h13period / wperiod) < 1.09);
    if ~her6_mutant
        continue;
    end
    % temp_rate.psh13 = rateValues[i].psh13; 
    
    % Her7 and Her6 mutant 
    temp_psh7 = 0.0; temp_psh6 = 0.0;
    her76_mutant_model_works = run_mutant( temp_rate, of_her713);
    if ~her76_mutant_model_works
        continue;
    end
    her76_mutant = ((h713period / wperiod) > 1.03) && ((h713period / wperiod) < 1.09);
    if ~her76_mutant
        continue;
    end
    
%     temp_rate.psh7 = rateValues[i].psh7;
%     temp_rate.psh13 = rateValues[i].psh13;
    

 
%              If the paramater set created oscillatory behavior in wild type and all the mutant conditions were satisfied:
%              1) Print the appropriate message
%              2) Print the oscillation features into the appropriate files
%              2) Print the parameter set into the output file containing parameter that passed the conditions.
%             
%             
%             cerr << terminal_blue << "Parameter set " << i << " passed." << terminal_reset << endl;
%             if (ofeat) {
%                 oft << res << "," << of_wt.period << "," << of_wt.amplitude << "," << of_wt.peaktotrough1 << ",";
%                 oft << of_delta.period << "," << of_delta.amplitude << "," << of_delta.peaktotrough1 << ",";
%                 oft << of_her1.period << "," << of_her1.amplitude << "," << of_her1.peaktotrough1 << ",";
%                 oft << of_her7.period << "," << of_her7.amplitude << "," << of_her7.peaktotrough1 << ",";
%                 oft << of_her13.period << "," << of_her13.amplitude << "," << of_her13.peaktotrough1 << ",";
%                 oft << of_her713.period << "," << of_her713.amplitude << "," << of_her713.peaktotrough1 << ",";
%             }
%             
%             allpassed<<temp_rate.psh1<<","<<temp_rate.psh7<<","<<temp_rate.psh13<<","<<temp_rate.psd<<","<<temp_rate.pdh1<<","<<temp_rate.pdh7<<",";
% 			allpassed<<temp_rate.pdh13<<","<<temp_rate.pdd<<","<<temp_rate.msh1<<","<<temp_rate.msh7<<","<<temp_rate.msh13<<","<<temp_rate.msd<<",";
% 			allpassed<<temp_rate.mdh1<<","<<temp_rate.mdh7<<","<<temp_rate.mdh13<<","<<temp_rate.mdd<<"," <<temp_rate.ddgh1h1 << "," << temp_rate.ddgh1h7 << ",";
% 			allpassed<<temp_rate.ddgh1h13<<","<<temp_rate.ddgh7h7<<","<<temp_rate.ddgh7h13<<","<<temp_rate.ddgh13h13<<",";
% 			allpassed<<temp_rate.delaymh1<<","<<temp_rate.delaymh7<<","<<temp_rate.delaymh13<<","<<temp_rate.delaymd<<","<<temp_rate.delayph1<<",";
% 			allpassed<<temp_rate.delayph7<<","<<temp_rate.delayph13<<","<<temp_rate.delaypd<<","<<temp_rate.dah1h1<<","<<temp_rate.ddh1h1<<","<<temp_rate.dah1h7<<",";
% 			allpassed<<temp_rate.ddh1h7<<","<<temp_rate.dah1h13<<","<<temp_rate.ddh1h13<<","<<temp_rate.dah7h7<<","<<temp_rate.ddh7h7<<","<<temp_rate.dah7h13<<",";
% 			allpassed<<temp_rate.ddh7h13<<","<<temp_rate.dah13h13<<","<<temp_rate.ddh13h13<<","<<temp_rate.critph1h1 << "," << temp_rate.critph7h13 <<","<<temp_rate.critpd<<endl;
% 		}

    end
 
function bool_value = model(eps, time_steps)
%returns a boolean indicting whether or not the model works

bool_value = false; %initializing the bool value as false
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Parameter values used in the 2013 study (from Table S1 of the paper)
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

CHUNK_SIZE = 10000;
PARS = 1; % number of parameters to simulate, default value is 1
x = 2; y = 1; %width and height of the tissue being simulated, default is 2x1
minutes=1200;
seed = rand;
eps = 0.01; %time step to be used for Euler's method, default is 0.01
time_steps = (minutes / eps);% Set the amount of time steps to be used in the simulation

max_prop = INFINITY; % maximum threshold for propensity functions, default is INFINITY
%toPrint = false;  % boolean marking whether or not concentrations should be printed to a text file
%ofeat = false;
%mutants = ["/wt", "/delta", "/her6", "/her1", "/her7", "/her76"];

% Runs the deterministic simulation of the model.
%      For each time step:
%      1) Iterate through every cell (2 cells in this case) and update the concentrations of proteins and mRNA.
%         These concentration values are obtained by solving the differential equations for that time step, using Euler's method.
%      2) Check that the concentrations do not become negative -- a negative amount of protein is not biologically sensible
%      3) Check that the propensity functions do not go above the set threshold -- if one was specified.

% Convert the time delay values to integers, because the deterministic
% simulation uses discrete time points:
    nmh1_eps=round(nmh1/eps);
    nmh7_eps=round(nhm7/eps);
    %nph6_eps=nph6/eps; -> not used
    nmd_eps=round(nmd/eps);
    nph1_eps=round(nph1/eps);
    nph6_eps=round(nph6/eps);
    nph7_eps=round(nph7/eps);
    npd_eps=round(npd/eps);
    
    cells=rows*columns; %two in this case
    
    a = zeros([2,34]); % Initialize the propensity levels for each reaction for both cells.
    % they will be checked later in the model
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %number of timesteps = time_steps = time_steps
    %so we need an array of time_steps values for each state and each cell
    
    %initializing the states: each state = array with cell no. and time at
    %which the state is calculated
    
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
    
    %functions to calculate the transcription rates of mRNA
    fd=@(ph11,ph76) msd/(1+(ph11/critph11)^2+(ph76/critph76)^2);
    fh1=@(ph11,ph76,pd)msh1*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);
    fh7=@(ph11,ph76,pd)msh7*(1+(pd/critpd))/(1+(pd/critpd)+(ph11/critph11)^2 +(ph76/critph76)^2);
    
    for n = 1:1:timesteps  %% n is an integer
        for i = 0:1:cells
            
            %//Protein synthesis
            if n>nph1_eps
                prod_ph1 = psh1*mh1(i,n-nph1_eps);
            else
                prod_ph1 = 0;    
            end
            ph1(i,n) = ph1(i,n-1)+ eps * (prod_ph1 - pdh1*ph1(i,n-1)-2*dah1h1*ph1(i,n-1)*ph1(i,n - 1)+2*ddh1h1*ph11(i,n - 1)-dah1h7*ph1(i,n - 1)*ph7(i,n - 1)+ddh1h7*ph17(i,n - 1)-dah1h6*ph1(i,n - 1)*ph6(i,n - 1)+ddh1h6*ph16(i,n - 1));
            
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
               
            
            %break out of the nested for-loops of any of them if any of the
            %checking for protein levels becoming negative
            if (ph1(i,n) < 0 || ph7(i,n) < 0 || ph6(i,n) < 0) 
                bool_value=false;
                break;
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
                bool_value=false;
                break;
            end
            
            % mRNA Synthesis
            if (n > nmh1_eps )
                prod_mh1 = fh1(ph11(i,n - nmh1_eps), ph76(i,n - nmh1_eps), pd(1 - i,n - nmh1_eps));
            else
                prod_mh1 = 0;
            end
            mh1(i,n) = mh1(i,n-1) + eps * (prod_mh1 - mdh1*mh1(i,n - 1)); 
            
            if (n > nmh7_eps)
                prod_mh7 = fh7(ph11(i,n - nmh7_eps), ph76(i,n - nmh7_eps), pd(1 - i,n - nmh7_eps));
            else
                prod_mh7 = 0;
            end
            mh7(i,n) = mh7(i,n-1) + eps * (prod_mh7- mdh7*mh7(i,n - 1));
            
            
			mh6(i,n) = mh6(i,n-1) + eps * (msh6-mdh6*mh6(i,n - 1)); 
            
            
            if n > nmd_eps
                prod_md = fd(ph11(i,n - nmd_eps), ph76(i,n - nmd_eps), pd(1 - i,n - nmd_eps));
            else
                prod_md = 0;
            end
			md(i,n) = md(i,n-1) + eps * (prod_md - mdd*md(i,n - 1));
            
            %checking for negative mRNA levels
            if (mh1(i,n) < 0 || mh7(i,n) < 0 || mh6(i,n) < 0 || md(i,n) < 0)
                bool_value = false;
            end
            
            %checking for negative and infinite propensities
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
                a(ck,5)=fh1(ph11(ck),ph76(ck),pd(1-ck)); % // Reaction 27: -> mh1
                a(ck,6)=fh7(ph11(ck),ph76(ck),pd(1-ck)); % // Reaction 29: -> mh7
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

            if ~(a<cutoff) && max_prop ~= inf %to satisfy this any one of the elements of a must satisfy a>=cutoff
                bool_value=false;
            end
            
        end
        
        if bool_value==false
            break;
        end
    end
end


function model_works = run_mutant(temp_rate, wild)%data &of , wild:boolean value

%  Performs the steps necessary in the simulation and analysis of the wild type or a certain mutant.
%  1) Clear the levels from the previous simulation.
%  2) Run the model for the specified duration.
%  3) Create the oscillation features of the simulation.
%  4) Return whether concentrations in the model where positive values below the propensity threshold
% 
    model_works = model(eps, t_steps, g, temp_rate, max_prop, x, y);
    ofeatures(g, eps, t_steps, wild, of); 
end

function ofeatures(mh1,wild)% data &d) 
    
%  Calculates the oscillation features -- period, amplitude and peak to trough
%  ratio for a set of concentration levels.
%  The values are calculated using the last peak and trough of the oscillations,
%  since the amplitude of the first few oscillations can be slightly unstable.
%  For the wild type, the peak and trough at the middle of the graph are also calculated
%  in order to ensure that the oscillations are sustained.
 
	tmaxlast = 0; tmaxpenult = 0; mmaxlast = 0; mmaxpenult = 0;
	tminlast = 0; tminpenult = 0; mminlast = 0; mminpenult = 0;
		
    for n = 1:1:nfinal - 1
        % check if the current point is a peak
        if mh1(0,n + 1) < mh1(0,n) && mh1(0,n) > mh1(0,n-1)
            tmaxpenult = tmaxlast;
            tmaxlast = n*eps;
            mmaxpenult = mmaxlast;
            mmaxlast = mh1(0,n);
        end
        % check if the current point is a trough
        if (mh1(0,n+1) > mh1(0,n) && mh1(0,n) < mh1(0,n-1)
            tminpenult = tminlast;
            tminlast = n * eps;
            mminpenult = mminlast;
            mminlast = mh1(0,n);
        end
        
    end
    
    if (wild)
        % calculate the peak and trough at the middle of the graph
        mmaxlast2 = 0.0; mminlast2 = 0.0;
        for m = 2:1:nfinal/2
            if(mh1(0,m + 1) < mh1(0,m) && mh1(0,m) > mh1(0,m - 1))
                mmaxlast2 = mh1(0,m);
            end
            if(mh1(0,m + 1) > mh1(0,m) && mh1(0,m) < mh1(0,m - 1))
                mminlast2 = mh1(0,m);
            end
        end
        % in order to avoid dividing by zero in case a trough is 0, set it to 1
        if(mminlast2 == 0.0 || mminlast == 0.0)
            mminlast2 = 1.0;
            mminlast = 1.0;
        end
        peaktotrough2 = mmaxlast2/mminlast2;
    end
	
	period = tmaxlast-tmaxpenult;
	amplitude = mmaxlast-mminlast;
	peaktotrough1 = mmaxlast/mminlast;
    
end