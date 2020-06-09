% deterministic simulator for the zebrafish segmentation clock
% MATLAB version of the C code by Ahmet Ay, Jack Holland, and Adriana Sperlea
 
%  The program simulates the behavior of the zebrafish segmentation clock for two cells,
%  through a mathematical model using delay differential equations.
%  To solve the equations, the program uses Euler's method, with an adjustable stepsize.
%  The program tries to find parameter sets which replicate the behavior of the system in wild type
%  and several mutants.

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

% In C code: global variables set in functions.cpp

% Also include?
% extern char* terminal_blue;
% extern char* terminal_red;
% extern char* terminal_reset;

CHUNK_SIZE = 10000;
PARS = 1; % number of parameters to simulate, default value is 1
x = 2; y = 1; %width and height of the tissue being simulated, default is 2x1
minutes=1200;
seed = rand;


eps = 0.01; %time step to be used for Euler's method, default is 0.01
max_prop = INFINITY; % maximum threshold for propensity functions, default is INFINITY
toPrint = false; ofeat = false; % boolean marking whether or not concentrations should be printed to a text file


% Iterate through every parameter set

for p = 0:CHUNK_SIZE:PARS


 
%              For the wild type and every mutant, perform the following steps:
%              1) Adjust the appropriate protein synthesis rates to create mutants if necessary
%              2) Run the simulation
%              3) Go to the next parameter set if the propensities for the current one have gone above the set threshold
%              4) Otherwise, test oscillation features
%              5) Go to the next parameter set if the current one did not produce oscillations or did not satisfy the mutant conditions


end

a = zeros([2,12]); % Initialize the propensity levels for each reaction for both cells.

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
a(5)=fh1; % // Reaction 27: -> mh1 
a(6)=fh7; % // Reaction 29: -> mh7
a(7) = psh6; % // Reaction 31: -> mh6		
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

function bool_value = model(eps, nfinal, columns, rows)

% Runs the deterministic simulation of the model.
%      For each time step:
%      1) Iterate through every cell (2 cells in this case) and update the concentrations of proteins and mRNA.
%         These concentration values are obtained by solving the differential equations for that time step, using Euler's method.
%      2) Check that the concentrations do not become negative -- a negative amount of protein is not biologically sensible
%      3) Check that the propensity functions do not go above the set threshold -- if one was specified.

    cells=rows*columns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %number of timesteps = nfinal
    %so we need an array of nfinal values for each state and each cell
    
    %initializing the states: each state = array with cell no. and time at
    %which the state is calculated
    ph1=zeros(cells,nfinal); 
    ph7=zeros(cells,nfinal);
    ph6=zeros(cells,nfinal);
    pd=zeros(cells,nfinal);
    mh1=zeros(cells,nfinal);
    mh7=zeros(cells,nfinal);
    mh6=zeros(cells,nfinal);
    md=zeros(cells,nfinal); 
    ph11=zeros(cells,nfinal);
    ph76=zeros(cells,nfinal);
    ph17=zeros(cells,nfinal);
    ph16=zeros(cells,nfinal);
	ph77=zeros(cells,nfinal);
	ph66=zeros(cells,nfinal);
   
    %for the dimers we don't need to store the values at each timestep, so
    %only store the values for each cell at any point - right?

 
    bool_value = false;
    
    for n = 1:1:nfinal  %% WHAT
        for i = 0:1:cells
            
            %//Protein synthesis
            if n>nph1
                prod_ph1 = psh1*mh1(i,n-nph1);
            else
                prod_ph1 = 0;    
            end
            ph1(i,n) = ph1(i,n-1)+ eps * (prod_ph1 - pdh1*ph1(i,n-1)-2*dah1h1*ph1(i,n-1)*ph1(i,n - 1)+2*ddh1h1*ph11(i,n - 1)-dah1h7*ph1(i,n - 1)*ph7(i,n - 1)+ddh1h7*ph17(i,n - 1)-dah1h6*ph1(i,n - 1)*ph6(i,n - 1)+ddh1h6*ph16(i,n - 1));
            
            if n > nph7
                prod_ph7 = psh7*mh7(i,n - nph7);
            else
                prod_ph7 = 0;
            end
            ph7(i,n) = ph7(i,n - 1) + eps * (prod_ph7 - pdh7*ph7(i,n - 1)-2*dah7h7*ph7(i,n - 1)*ph7(i,n - 1)+2*ddh7h7*ph77(i,n - 1)-dah1h7*ph1(i,n - 1)*ph7(i,n - 1)+ddh1h7*ph17(i,n - 1)-dah7h6*ph7(i,n - 1)*ph6(i,n - 1)+ddh7h6*ph76(i,n - 1)); 
            
            
            if n > nph6 
                prod_ph6 = psh6*mh6(i,n - nph6);
            else
                prod_ph6 = 0;
            end
            ph6(i,n) = ph6(i,n-1) + eps * (prod_ph6 -pdh6*ph6(i,n - 1)-2*dah6h6*ph6(i,n - 1)*ph6(i,n - 1)+2*ddh6h6*ph66(i,n - 1)-dah1h6*ph1(i,n - 1)*ph6(i,n - 1)+ddh1h6*ph16(i,n - 1)-dah7h6*ph7(i,n - 1)*ph6(i,n - 1)+ddh7h6*ph76(i,n - 1));
               
            
            %break out of the nested for-loops of any of them if any of the
            %checking for protein levels becoming negative
            if (ph1(i,n) < 0 || ph7(i,n) < 0 || ph6(i,n) < 0) 
                %negative_flag = true; 
                bool_value=false;
                break;
            end
            
            %//Dimer proteins
            ph11(i,n) = ph11(i,n - 1) + eps * (dah1h1*ph1(i,n - 1)*ph1(i,n - 1)-ddh1h1*ph11(i,n - 1)-ddgh1h1*ph11(i,n - 1));
            ph17(i,n) = ph17(i,n - 1) + eps * (dah1h7*ph1(i,n - 1)*ph7(i,n - 1)-ddh1h7*ph17(i,n - 1)-ddgh1h7*ph17(i,n - 1));
            ph16(i,n) = ph16(i,n - 1) + eps * (dah1h6*ph1(i,n - 1)*ph6(i,n - 1)-ddh1h6*ph16(i,n - 1)-ddgh1h6*ph16(i,n - 1));
            ph77(i,n) = ph77(i,n - 1) + eps * (dah7h7*ph7(i,n - 1)*ph7(i,n - 1)-ddh7h7*ph77(i,n - 1)-ddgh7h7*ph77(i,n - 1));
            ph76(i,n) = ph76(i,n - 1) + eps * (dah7h6*ph7(i,n - 1)*ph6(i,n - 1)-ddh7h6*ph76(i,n - 1)-ddgh7h6*ph76(i,n - 1));
            ph66(i,n) = ph66(i,n - 1) + eps * (dah6h6*ph6(i,n - 1)*ph6(i,n - 1)-ddh6h6*ph66(i,n - 1)-ddgh6h6*ph66(i,n - 1));

            %// Delta Protein
            if n > npd 
                prod_pd = psd*md(i,n-npd);
            else
                prod_pd = 0;
            end
            pd(i,n) = pd(i,n-1) + eps*(prod_pd - pdd*pd(i,n - 1));
            
            
            if (ph11(i,n) < 0 || ph17(i,n) < 0 || ph16(i,n) < 0 || ph77(i,n) < 0 || ph76(i,n) < 0 || ph66(i,n) < 0 || pd(i,n) < 0)
                %negative_flag = true; 
                bool_value=false;
                break;
            end
            
            
            avgpdh1 = 0; avgpdh7 = 0; avgpdd = 0;
            % delta from other cells - for this model, delta from the other
            % cell since it is a 2-cell model
            
            if (rows == 1 && columns == 2) % will always be true in this case
                % the delta input is coming from the other cell
                avgpdh1 = pd(1 - i,n - nmh1);
                avgpdh7 = pd(1 - i,n - nmh7);
                avgpdd = pd(1 - i,n - nmd);
            end
            
            % mRNA Synthesis
            if (n > nmh1 )
                prod_mh1 = fh1(ph11(i,n - nmh1), ph76(i,n - nmh1), avgpdh1, msh1, critph1h1, critph7h13, critpd);
            else
                prod_mh1 = fh1(0, 0, 0, msh1, critph1h1, critph7h13, critpd);
            end
            mh1(i,n) = mh1(i,n-1) + eps * (prod_mh1 - mdh1*mh1(i,n - 1)); 
            
            if (n > nmh7)
                prod_mh7 = fh7(ph11(i,n - nmh7), ph76(i,n - nmh7), avgpdh7, msh7, critph1h1, critph7h13, critpd);
            else
                prod_mh7 = fh7(0, 0, 0, msh7, critph1h1, critph7h13, critpd);
            end
            mh7(i,n) = mh7(i,n-1) + eps * (prod_mh7- mdh7*mh7(i,n - 1));
            
            
			mh6(i,n) = mh6(i,n-1) + eps * (msh6-mdh6*mh6(i,n - 1)); 
            
            
            if n > nmd 
                prod_md = fd(ph11(i,n - nmd), ph76(i,n - nmd), avgpdd, msd, critph1h1, critph7h6, critpd);
            else
                prod_md = fd(0, 0, 0, msd, critph1h1, critph7h6, critpd);
            end
			md(i,n) = md(i,n-1) + eps * (prod_md - mdd*md(i,n - 1));
            
            %checking for negative mRNA levels
            if (mh1(i,n) < 0 || mh7(i,n) < 0 || mh6(i,n) < 0 || md(i,n) < 0)
                bool_value = false;
            end
            
            %checking for negative and infinite propensities
            if (max_prop ~= inf && ~checkPropensities(g, r, n, max_prop)) 
                bool_value=false;
            end
        end
        
        if bool_value==false
            break;
        end
        
        
        
    end
end