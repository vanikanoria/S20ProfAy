function [T,Y,VL,RC] = HAART(WriteOutput)
%
% This function can be used to simulate the kinetics of HIV after application 
% of conventional HAART treatment (see von Kleist M, Menz S, Stocker H,
% Arasteh K, Huisinga W, Schuette S (2011) "HIV Quasispecies Dynamics during
% Pro-active Treatment Switching: Impact on Multi-Drug Resistance and Resistance
% Archiving in Latent Reservoirs", PlosOne)
%
%   Input: 
%
%       Optional boolean to write output. If set to true all output will be
%       written to seperate files.
%
%   Returns:
%   
%       1st) Vector of simulated time points in days
%   
%       2nd) Matrix of species levels (columns) at every time point (rows)
%   
%       3rd) Vector of viral load at every time point
%   
%       4th) Matrix of reproductive capacities,
%            w.r.t. the drug combinations (columns) and time points (rows)
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


% specify available drugs and drug combinations

% drugs and drug combination considered in the main article
% Drugs      = [2 0.9; 2 0.9; 2 0.9; 1 0.9; 1 0.9; 4 0.9];
% DrugCombos = [1 2 3; 4 5 6];
%
% each row in 'Drugs' defines a drug with parameters (columns):
%   1st - target:
%       set to 1 for entry inhibitors
%              2 for reverse transcriptase inhibitors
%              3 for integrase inhibitors
%              4 for protease / maturation inhibitors
%   2nd - efficacy: value from 0 (not effective) to 1 (full effective)
%
% each row in 'DrugCombos' specifies a drug combination by
% referring to the corresponding indices in 'Drugs'

Drugs = [   2 0.8;...
            2 0.8;...
            1 0.8;...
            4 0.8];
DrugCombos = [1 2; 3 4];

% set selective disadvantage and resistance of the mutants (relative to wildtype)
Selective_Disadvantage = 0.2;
Resistance             = 0.01;  % 0.01 = 99% resistance; 0.2 = 80% resistance, etc
if Resistance <= 0              % value has to be greater than 0
    Resistance = 1e-6;
end

% set duration of simulation (in days) 
% Tend = 735;  % value used in the main article
Tend = 10;

% construct HIV-model
fprintf(1,'\nConstructing HIV-model:\t\t'); tic;
M = HIVmodel(Drugs,DrugCombos,Selective_Disadvantage,Resistance);
fprintf(1,'%3.2E seconds',toc);

% thresholds for stochastic-deterministic partitioning 
% ysd = number of reactant below which reaction is handled stochastically 
% asd = value of the propentity function below which reaction is handled stochastically  
% smaller values-> faster simulation, less accurate
% M.set_asd(20);	% value used in the main article
% M.set_ysd(20);    % value used in the main article
M.set_asd(2);
M.set_ysd(5);

% approximate initial state by a deterministic pre-run:
%   1st argument: length of pre-run in days
%   2nd argument: boolean, specifying if the initial state should contain wild type only
fprintf(1,'\nComputing initial state:\t'); tic;
y = M.set_initial_conditions(500,false);
fprintf(1,'%3.2E seconds',toc);

% initialize patient monitoring
dT_monitoring = Tend;   % time interval between monitorings in days
% dT_monitoring = 30;	% value used in the main article
P             = PatientMonitoring(Tend,y,M.get_viral_load,dT_monitoring);

% set event function of stochastic reactions and patient monitoring
function [value,isterminal,direction] = events(t,y)
    
    [v1,i1,d1] = M.reaction_event(y);
    [v2,i2,d2] = P.patient_monitoring(t);
    
    value      = [v1 v2];
    isterminal = [i1 i2];
    direction  = [d1 d2];
end

% set ode options
options = odeset('Events',@events,'NonNegative',1:M.NrSpecies);

% choose time discretization of output (w.r.t. days)
dT = 1;

% initialize output variables
T = 0:dT:Tend; T = T(:); % vector of time points in unit [days]
if T(end) ~= Tend
    T = [T; Tend];
end
NrSteps = length(T);

Y  = zeros(NrSteps,M.NrSpecies); Y(1,:) = y(1:M.NrSpecies);                         % species levels at every time point
VL = zeros(NrSteps,1); VL(1) = M.get_viral_load(y);                                 % viral load at every time point
RC = zeros(NrSteps,size(DrugCombos,1)); RC(1,:) = M.get_reproductive_capacities(y);	% reproductive capacities at every time point

% activate switch to next drug combination and start simulation
switch_to_nextcombo = true;
fprintf(1,'\nSimulated '); tic;
t = 0; step = 2; fprintf(1,'t = %3.2E days:\t%3.2E seconds',t,toc);
while t < Tend
    
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
            fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3.2E days:\t%3.2E seconds',t,toc);
            
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

% check if treatment finished succssful
if t == Tend
	fprintf(1,'\n\n\tFinished successful!\n');
else
    fprintf(1,'\n\n\tFailure detected!\n');
    
    % delete unused entries in output variables
    T  = T(1:step);
    Y  = Y(1:step,:);
    VL = VL(1:step);
    RC = RC(1:step,:);
end

% plot results
figure(); clf;
semilogy(T,Y,'-o');
xlabel('time in days');
ylabel('level of species');

figure(); clf;
subplot(1,2,1);
semilogy(T,VL,'-o');
xlabel('time in days');
ylabel('viral load [HIV RNA/mL plasma]');
subplot(1,2,2);
semilogy(T,RC,'-o');
xlabel('time in days');
ylabel('reproductive capacities');

% write output files if required
if nargin
    if WriteOutput
        fprintf(1,'\nWrite output:\t\t\t'); tic;
        dlmwrite('Time.csv',T,'delimiter',',','precision','%.2f');
        dlmwrite('SpeciesLevels.csv',Y,'delimiter',',','precision','%.2E');
        dlmwrite('ViralLoad.csv',VL,'delimiter',',','precision','%.2E');
        dlmwrite('RepCaps.csv',RC,'delimiter',',','precision','%.2E');
        fprintf(1,'%3.2E seconds\n\n',toc);
    end
end
         
end