function M = HIVmodel(Drugs,DrugCombos,Selective_Disadvantage,Resistance)
%
% This function builds the original HIV model (see von Kleist M, Menz S,
% Stocker H, Arasteh K, Huisinga W, Schuette S (2011) "HIV Quasispecies
% Dynamics during Pro-active Treatment Switching: Impact on Multi-Drug
% Resistance and Resistance Archiving in Latent Reservoirs", PlosOne, and
% von Kleist M, Menz S, Huisinga W (2010) "Drug-Class Specific Impact of 
% Antivirals on the Reproductive Capacity of HIV", PLoS Comput Biol 6(3):e1000720
%
%   Input arguments are:
%   
%       1st) Matrix of available drugs, specifying drug target (1st column)
%            and drug efficacies (2nd column)
%   
%       2nd) Matrix of drug combinations, each row specifies a combination
%            with entries referring to the index of the drug in the 1st
%            argument
%   
%       3rd) Selective disadvantage (relative to wildtype)
%   
%       4th) Resistance (relative to wildtype) of each mutant against the
%            corresponding drug
%
% 
% Copyright (C) 2011, Free University Berlin
% Contact: vkleist@zedat.fu-berlin.de
%
% FOR ACADEMIC USE this program is free software; you can redistribute 
% it and/or modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.


%%%%%%%%%% START %%%%%%%%%%


% get drug specific parameters
AvailableTargets = 4;
NrDrugs          = size(Drugs,1);         % total number of drugs
NrCombinations   = size(DrugCombos,1);    % total number of drug combinations

% check consistancy
if ~isempty(find(Drugs(:,1) < 1 | Drugs(:,1) > AvailableTargets, 1))
    error('Model:argChk',...
        ['\nAvailable indices of drug classes are:\n\n'...
            '\t1 - entry inhibitors (FI/CCR5)\n'...
            '\t2 - reverse transcriptase inhibitors (NRTI/NNRTI)\n'...
            '\t3 - integrase (InI) \n'...
            '\t4 - protease (PI)/ maturation inhibitors (MI)']);
end
if ~isempty(find(DrugCombos < 1 | DrugCombos > NrDrugs, 1))
    error('Model:argChk',...
        '\nSpecify drug combinations by reffering to valid indices of drugs.');
end


% set system specific parameters
NrMutants          = 2^NrDrugs;         % total number of different mutants
NrSpeciesPerMutant = 6;                 % number of species per mutant
NrSpecies          = 3 + NrSpeciesPerMutant*NrMutants; % total number of different species

NrReactionsPerMutant      = 3*NrMutants + 9;	% number of reactions per mutant
NrDeathReactionsPerMutant = 8;                  % number of death reactions per mutant
NrReactions               = 5 + (NrReactionsPerMutant+NrDeathReactionsPerMutant)*NrMutants; % total number of different reactions


% model parameters used for simulation (see Tabel 1 in main article)
lambda_T      = 2.0e9;      % birth rate of unifected T-cells in unit [1/day]
lambda_M      = 6.9e7;      % birth rate of unifected macrophages in unit [1/day]
death_T       = 0.02;       % death rate of unifected T-cells in unit [1/day]
death_T1      = 0.02;       % death rate of T1-cells in unit [1/day]
death_T2      = 1.0;        % death rate of T2-cells in unit [1/day]
death_TL      = 0.0001;     % death rate of latently infected T-cells (TL-cells) in unit [1/day]
death_M       = 0.0069;     % death rate of unifected macrophages in unit [1/day]
death_M1      = 0.0069;     % death rate of M1-cells in unit [1/day]
death_M2      = 0.09;       % death rate of M2-cells in unit [1/day]
deathPIC_T    = 0.35;       % intracellular degradation rate of pre-integration complex within T-cells in unit [1/day]
deathPIC_M    = 0.0035;     % intracellular degradation rate of pre-integration complex within macrophages in unit [1/day]
beta_T        = 8.0e-12;    % rate of successfull infection of T-cells in unit [1/day]
beta_M        = 1.0e-14;    % rate of successfull infection of mcarophages in unit [1/day]
k_T           = 0.35;       % integration rate constant in T-cells in unit [1/day]
k_M           = 0.07;       % integration rate constant in M-cells in unit [1/day]
n_T           = 1000;       % total number of released virus from T2-cells in unit [1/day]
n_M           = 100;        % total number of released virus from M2-cells in unit [1/day]
CL            = 23.0;       % clearance rate of virus in unit [1/day]
mu            = 2.16e-5;    % point mutation probability per base and reverse transcritpion process
rho           = 0.33;       % probability that reverse transcription is successfully completed
p_latent      = 8e-6;       % probability that a T-cell becomes latently infected
alpha         = 0.001;      % activation rate of TL-cells
frac_uninfect = 0.33;       % fraction of uninfectious virus
PlasmaFac     = 236750.0;	% blood plasma conversion factor in unit [mL plasma/HIV RNA]


% derived auxilliary variables
one_minus_select_dis   = 1.0 - Selective_Disadvantage;
frac_infect            = 1.0 - frac_uninfect;
nT_dT2                 = frac_infect*n_T/death_T2;
nM_dM2                 = frac_infect*n_M/death_M2;
rho_inv                = 1.0/rho;
p_unlatent             = 1.0 - p_latent;
one_minus_p_death_r_TL = 1.0 - p_latent*death_TL/(alpha + death_TL);
alpha_r_TL             = alpha/(alpha + death_TL);


% compute matrices of selective disadvantage and resistances of the mutants
SelectionMatrix  = ones(AvailableTargets,NrMutants);    % selective disadvantages of mutants for every target (1st dimension target, 2nd mutant)
MutationMatrix   = zeros(NrDrugs,NrMutants);            % mutation scheme of mutants against the drugs (1st dimension drug, 2nd mutant)
for i = 1:NrMutants
    
    % get possible mutation scheme
    MutationMatrix(:,i) = (int16(bitget(i-1, NrDrugs:-1:1)))';
    
    % find resistant mutations, and
    DrugIDs = find(MutationMatrix(:,i));
    % set slective disadvantage and resistance accordingly
    if ~isempty(DrugIDs)
        for j = 1:length(DrugIDs)
            actualDrug = DrugIDs(j);
            SelectionMatrix(Drugs(actualDrug,1),i)  = SelectionMatrix(Drugs(actualDrug,1),i).*one_minus_select_dis;
        end
    end
end
% resistances of mutants for every target (1st dimension drug, 2nd mutant)
ResistanceMatrix = double(MutationMatrix)*Resistance;	% compute resistances of the mutants
ResistanceMatrix(ResistanceMatrix == 0) = 1;


% compute distances of the mutations
p = zeros(NrMutants,NrMutants); % mutation probabilities between the different strains
for i = 1:NrMutants
    for j = 1:NrMutants
        p(i,j) = sum(abs(MutationMatrix(:,i) - MutationMatrix(:,j)));
    end
end
p = (mu.^p).*((1-mu).^(NrDrugs-p))';
clear MutationMatrix mu


% compute the efficacies of the considered drug combinations
Efficacies = ones(NrCombinations,AvailableTargets,NrMutants);   % drug efficacies against mutants for every drug combination and target (1st dimension drug combination, 2nd drug target, 3rd mutant)
for i = 1:NrCombinations
    for k = 1:NrMutants
        for j = DrugCombos(i,:)
            Efficacies(i,Drugs(j,1),k) = Efficacies(i,Drugs(j,1),k)*(1 - Drugs(j,2)*ResistanceMatrix(j,k));
        end
    end
end
clear ResistanceMatrix


% set system stoichiometry
R = sparse(NrSpecies,NrReactions);                  % net changes in species levels caused by firing of the reactions (1st dimension species, 2nd reaction)
V = repmat(struct('reactants',[]),1,NrReactions-2); % index structure of the reactants of every reactions

R(1,1) = 1; % birth of TU:  -> TU
R(2,2) = 1; % birth of MU:  -> MU

j = 3;      % index of reaction (row index)
for i = 1:NrMutants     % iterate through mutants

    % infection: TU + V(i) -> T1(i)
    R(1,j)                          = -1;	% TU
    R(8+NrSpeciesPerMutant*(i-1),j) = -1;	% V(i)
    R(3+NrSpeciesPerMutant*(i-1),j) = +1;	% T1(i)
    V(j-2).reactants = [1 8+NrSpeciesPerMutant*(i-1)];	% reactants are TU and V(j)
    j = j+1;

    % clearance of PIC: T1(i) -> TU
    R(3+NrSpeciesPerMutant*(i-1),j) = -1;	% T1(i)
    R(1,j)                          = +1;	% TU
    V(j-2).reactants = 3+NrSpeciesPerMutant*(i-1);      % reactant is T1(i)
    j = j+1;

    % I. transition/mutation: T1(q) -> T2(i)
    k = 0;
    for q = j:j+NrMutants-1
        R(3+NrSpeciesPerMutant*k,q) = -1;	% T1(q)
        V(q-2).reactants = 3+NrSpeciesPerMutant*k;      % reactant is T1(q)
        k = k+1;
    end
    R(4+NrSpeciesPerMutant*(i-1),j:j+NrMutants-1) = +1;	% T2(i)
    j = j+NrMutants;

    % II. transition/mutation: T1(q) -> TL(i)
    k = 0;
    for q = j:j+NrMutants-1
        R(3+NrSpeciesPerMutant*k,q) = -1;	% T1(q)
        V(q-2).reactants = 3+NrSpeciesPerMutant*k;      % reactant is T1(q)
        k = k+1;
    end
    R(5+NrSpeciesPerMutant*(i-1),j:j+NrMutants-1) = +1;	% TL(i)
    j = j+NrMutants;

    % activation of latent cells: TL(i) -> T2(i)
    R(5+NrSpeciesPerMutant*(i-1),j) = -1;	% TL(i)
    R(4+NrSpeciesPerMutant*(i-1),j) = +1;	% T2(i)
    V(j-2).reactants = 5+NrSpeciesPerMutant*(i-1);      % reactant is TL(i)
    j = j+1;

    % infection: MU + V(i) -> M1(i)
    R(2,j)                          = -1;	% MU
    R(8+NrSpeciesPerMutant*(i-1),j) = -1;	% V(i)
    R(6+NrSpeciesPerMutant*(i-1),j) = +1;	% M1(i)
    V(j-2).reactants = [2 8+NrSpeciesPerMutant*(i-1)];	% reactants are MU and V(i)
    j = j+1;

    % clearance of PIC: M1(i) -> MU
    R(6+NrSpeciesPerMutant*(i-1),j) = -1;	% M1(i)
    R(2,j)                          = +1;	% MU
    V(j-2).reactants = 6+NrSpeciesPerMutant*(i-1);      % reactant is M1(i)
    j = j+1;

    % III. transition/mutation: M1(q) -> M2(i)
    k = 0;
    for q = j:j+NrMutants-1
        R(6+NrSpeciesPerMutant*k,q) = -1;	% M1(q)
        V(q-2).reactants = 6+NrSpeciesPerMutant*k;      % reactant is M1(q)
        k = k+1;
    end
    R(7+NrSpeciesPerMutant*(i-1),j:j+NrMutants-1) = +1;	% M2(i)
    j = j+NrMutants;

    % Ia. production: T2(i) -> V(i) + T2(i)
    R(8+NrSpeciesPerMutant*(i-1),j) = +1;           % V(i)
    V(j-2).reactants = 4+NrSpeciesPerMutant*(i-1);	% reactant is T2(i)
    j = j+1;

    % Ib. production: T2(i) -> V_NI + T2(i)
    R(end,j) = +1;                                  % V_NI
    V(j-2).reactants = 4+NrSpeciesPerMutant*(i-1);	% reactant is T2(i)
    j = j+1;

    % IIa. production: M2(i) -> V(i) + M2(i)
    R(8+NrSpeciesPerMutant*(i-1),j) = +1;           % V(i)
    V(j-2).reactants = 7+NrSpeciesPerMutant*(i-1);	% M2(i)
    j = j+1;

    % IIb. production: M2(i) -> V_NI + M2(i)
    R(end,j) = +1;                                  % V_NI
    V(j-2).reactants = 7+NrSpeciesPerMutant*(i-1);	% M2(i)
    j = j+1;  
end

% death reactions
r = 2 + NrReactionsPerMutant*NrMutants;

% death of TU: TU ->  
R(1,r+1) = -1;          % TU
V(r-1).reactants = 1;   % reactant is TU

% death of MU: MU ->  
R(2,r+2) = -1;          % MU
V(r).reactants = 2;     % reactant is MU

j = 3; l = 3;
for i = 1:NrMutants
    n = l;	% all natural death events
    for s = j:j+NrSpeciesPerMutant-1
        R(s,r+n)            = -1;
        V(r+n-2).reactants	= s;
        n = n+1;
    end
    j = j+NrSpeciesPerMutant;
    l = l+NrSpeciesPerMutant;

    % Ia. clearance of virus by infection: TU + V(i) -> TU
    R(j-1,r+l)          = -1;       % V(i)
    V(r+l-2).reactants	= [1 j-1];  % reactants are TU and V(i)

    % clearance of virus by infection: MU + V(i) -> MU
    R(j-1,r+l+1)        = -1;       % V(i)
    V(r+l-1).reactants	= [2 j-1];  % reactants are MU and V(i)
    l = l+2;
end

% death of V_NI: V_NI ->  
R(NrSpecies,NrReactions)   = -1;        % V_NI
V(NrReactions-2).reactants = NrSpecies; % reactant is V_NI


% reactions constants (only birth/death reactions)
c = zeros(NrReactions,1);

% set fixed values of c
c(1) = lambda_T;    %  -> TU
c(2) = lambda_M;	%  -> MU
for i = 0:NrMutants-1
    c(4 + i*NrReactionsPerMutant)										 = deathPIC_T;	% T1(i) -> TU(i)
    c(5 + 2*NrMutants + i*NrReactionsPerMutant)                          = alpha;       % TL(i) -> T2(i)
    c(7 + 2*NrMutants + i*NrReactionsPerMutant)							 = deathPIC_M;  % M1(i) -> MU(i)
    c(5 + NrReactionsPerMutant*NrMutants + i*NrDeathReactionsPerMutant)  = death_T1;    % T1(i) ->  
    c(6 + NrReactionsPerMutant*NrMutants + i*NrDeathReactionsPerMutant)  = death_T2;    % T2(i) ->  
    c(7 + NrReactionsPerMutant*NrMutants + i*NrDeathReactionsPerMutant)  = death_TL;    % TL(i) ->  
    c(8 + NrReactionsPerMutant*NrMutants + i*NrDeathReactionsPerMutant)  = death_M1;    % M1(i) ->  
    c(9 + NrReactionsPerMutant*NrMutants + i*NrDeathReactionsPerMutant)  = death_M2;    % M2(i) ->  
    c(10 + NrReactionsPerMutant*NrMutants + i*NrDeathReactionsPerMutant) = CL;          %  V(i) ->  
end
c(3 + NrReactionsPerMutant*NrMutants) = death_T;    % TU ->  
c(4 + NrReactionsPerMutant*NrMutants) = death_M;    % MU ->  
c(NrReactions)                        = CL;         % V_NI ->  


% reaction propensities
a = zeros(NrReactions,1);
function [] = propensities(y)
    
    % initialize propensities with reaction constants
    a = c;
    
    % first two reactions are independent from statse y
    for i = 3:NrReactions % first two reactions are independent from y
        a(i) = a(i)*prod(y(V(i-2).reactants));
    end
end


% partition of the reactions (zero = deterministic, one = stochastic)
partition = zeros(NrReactions,1);
asd       = 0;  % all reactions with propensities less than 'asd_dflt' will be treated stochastically, otherwise deterministically
ysd       = 0;	% all reactions with any reactant level less than 'ysd_dflt' will be treated stochastically, otherwise deterministically
function [] = set_partition(y)

    propensities(y);
    for r = 3:NrReactions
        if a(r) < asd
            partition(r) = 1;
        else
            Reactants = V(r-2).reactants;
            for s = 1:length(Reactants)
                if R(Reactants(s),r) && y(Reactants(s)) < ysd
                    partition(r) = 1;
                    break
                end
            end
        end
    end
end


% return value of the right hand side (ydot) at time t for given state y
function ydot = RHS(t,y)
    
    % get propensities
    propensities(y);
    
    % compute right hand side
    ydot = [R*((1-partition).*a); partition'*a];
end


% set all reaction rate constants for specified drug combination dc
function [] = set_drugcombination(dc)
    
    % get efficacies that correspond to the drug combination dc
    if dc
        efficacies = reshape(Efficacies(dc,:,:),AvailableTargets,NrMutants);
    else
        efficacies = ones(AvailableTargets,NrMutants);
    end

    r = 3;  % rate constants of first two reactions never change
    
    % iterate for mutants
    for i = 1:NrMutants

        c(r) = beta_T*SelectionMatrix(1,i)*efficacies(1,i)*SelectionMatrix(2,i)*efficacies(2,i); % TU + VI(i) -> T1(i)  [successful infection of T-cells]
        r    = r+2;

        for j = 1:NrMutants
            c(r)               = p_unlatent*p(j,i)*k_T*SelectionMatrix(3,j)*efficacies(3,j);	% T1(i) -> T2(j) [proviral integration + phenotypic mutation]
            c(r+NrMutants)     = p_latent*p(j,i)*k_T*SelectionMatrix(3,j)*efficacies(3,j);      % T1(i) -> TL(j) [proviral integration + phenotypic mutation]
            c(r+2*NrMutants+3) = p(j,i)*k_M*SelectionMatrix(3,j)*efficacies(3,j);               % M1(i) -> M2(j) [proviral integration + phenotypic mutation]
            r = r+1;
        end
        r = r+ NrMutants + 1;

        c(r) = beta_M*SelectionMatrix(1,i)*efficacies(1,i)*SelectionMatrix(2,i)*efficacies(2,i); % MU + VI(i) -> M1(i) [successful infection of M-cells]
        r = r+NrMutants+2;

        c(r) = n_T*frac_infect*SelectionMatrix(4,i)*efficacies(4,i); % T2(i) -> VI(i) + T2(i) [release of infectious Virus from productively infect. T-cells]
        r = r+1;

        c(r) = n_T*(1.0 - frac_infect*SelectionMatrix(4,i)*efficacies(4,i)); % T2(i) -> VNI + T2(i) [release of non-infectious Virus from productively infect. T-cells]
        r = r+1;

        c(r) = n_M*frac_infect*SelectionMatrix(4,i)*efficacies(4,i);  % M2(i) -> VI(i) + M2(i) [release of infectious Virus from productively infect. M-cells] 
        r = r+1;

        c(r) = n_M*(1.0 - frac_infect*SelectionMatrix(4,i)*efficacies(4,i)); % M2(i) -> VNI + M2(i) [release of non-infectious Virus from productively infect. M-cells]
        r = r+1;

        c(11 + NrReactionsPerMutant*NrMutants + (i-1)*NrDeathReactionsPerMutant) = beta_T*SelectionMatrix(1,i)*efficacies(1,i)*(rho_inv - SelectionMatrix(2,i)*efficacies(2,i)); % TU + VI(i) -> TU  [clearance of virus by unsuccessful infection of T-cells]
        c(12 + NrReactionsPerMutant*NrMutants + (i-1)*NrDeathReactionsPerMutant) = beta_M*SelectionMatrix(1,i)*efficacies(1,i)*(rho_inv - SelectionMatrix(2,i)*efficacies(2,i)); % MU + VI(i) -> MU  [clearance of virus by unsuccessful infection of M-cells]
    end
end


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


% event function of stochastic reaction
function [value,isterminal,direction] = reaction_event(y)
    value      = y(end);
    isterminal = 1;
    direction  = 1;
end


% perform stochastic reaction event
function y = perform_reaction_event(y)
    propensities(y);                            % compute propensities of all reactions
    a = partition.*a;                           % take only stochastic reaction
    r = find(cumsum(a) >= sum(a)*rand,1);       % and choose a reaction (first reaction algorithm by Gillespie)
    y = [max(y(1:end-1)+R(:,r)',0) log(rand)];	% update species levels accordingly and draw new random number
    
    % update partitioning of reactions
    set_partition(y);
end


% compute and return viral load at state y
function vl = get_viral_load(y)
    
    vl = y(end-1);  % V_NI

    for m = 1:NrMutants
		vl = vl+y((m-1)*NrSpeciesPerMutant + 8);
    end
    
	vl = vl/PlasmaFac;
end


% compute and return reproductive capacities at state y
function rc = get_reproductive_capacities(y)
    
    rc = zeros(1,NrCombinations);
    for i = 1:NrCombinations

        n = 3; r = 6;
        for m = 1:NrMutants
            
            CL_by_infT = SelectionMatrix(1,m)*Efficacies(i,1,m);
            TUbeta     = CL_by_infT*SelectionMatrix(2,m)*Efficacies(i,2,m);
            CL_by_infT = CL_by_infT/rho - TUbeta;
            CL_by_infM = y(2)*beta_M*CL_by_infT;
            CL_by_infT = CL_by_infT*y(1)*beta_T;
            MUbeta     = y(2)*beta_M*TUbeta;
            TUbeta     = TUbeta*y(1)*beta_T;
            
            trans_rateT = SelectionMatrix(3,m)*Efficacies(i,3,m);
            trans_rateM = k_M*trans_rateT;
            trans_rateT = trans_rateT*k_T;
            
            trans_deathT = trans_rateT/(death_T1 + deathPIC_T + trans_rateT);
            trans_deathM = trans_rateM/(death_M1 + deathPIC_M + trans_rateM);
            
            rc(i) = rc(i) + (TUbeta + MUbeta) * (y(n)*nT_dT2*trans_deathT*one_minus_p_death_r_TL + y(r)*nM_dM2*trans_deathM);
            n = n+1; r = r+1;
            
            rc(i) = rc(i) + (TUbeta*trans_deathT + MUbeta*trans_deathM) * (y(n)*nT_dT2 + y(n+1)*nT_dT2*alpha_r_TL + y(r)*nM_dM2);
            r = r+1;
            
            rc(i) = rc(i) + y(r)*(TUbeta*nT_dT2*trans_deathT*one_minus_p_death_r_TL + MUbeta*nM_dM2*trans_deathM);

            rc(i) = rc(i) / (CL + CL_by_infT + TUbeta + CL_by_infM + MUbeta)*SelectionMatrix(4,m)*Efficacies(i,4,m);
            n = n+5; r = r+4;
        end
    end
end


% return functions for interaction with the model
M.NrSpecies                   = NrSpecies;
M.RHS                         = @RHS;
M.set_drugcombination         = @set_drugcombination;
M.set_initial_conditions      = @set_initial_conditions;
M.set_asd                     = @set_asd;
M.set_ysd                     = @set_ysd;
M.set_partition               = @set_partition;
M.reaction_event              = @reaction_event;
M.perform_reaction_event      = @perform_reaction_event;
M.get_viral_load              = @get_viral_load;
M.get_reproductive_capacities = @get_reproductive_capacities;

end