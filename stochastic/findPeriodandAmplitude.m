function [period, amplitude]= findPeriodandAmplitude(Data)

%the input consists of time and mh1 concentrations
Time = Data(:,1);
time_steps = size(Data,1); %this is the number of time_steps in the stochastic simulation

% if Data has two columns it has the conc of 1 cell
if size(Data,2) == 2
    num_cells = 1;
else % if Data has three columns it has the conc of 2 cells
    num_cells = 2;
end
mh1 = Data(:,2:end);

period_array = zeros(1, num_cells);
amplitude_array = zeros(1, num_cells);

%use a sliding average to smoothen the data
for each_cell = 1:num_cells
    mh1(:,each_cell) = movmean(mh1(:,each_cell),40);
end

for cell_num = 1:num_cells 
    %  Calculates the oscillation features -- period, amplitude and peak to trough
%  ratio for a set of concentration levels.
%  The values are calculated using the average of the oscillations after 
% the first three oscillations
%  since the amplitude of the first few oscillations can be slightly unstable.

    tmaxlast = 0;  mminlast = 0;

    periods =[];
    amplitudes = [];
		
    for n = 1:1:time_steps - 1
        %check if the current point is a peak
        %checking that it's greater than 5 before and 5 after
        if (mh1(n + 1,cell_num) < mh1(n,cell_num) && mh1(n,cell_num) > mh1(n-1,cell_num))
            
            if (mh1(n+2,cell_num) < mh1(n,cell_num) && mh1(n+3,cell_num) < mh1(n,cell_num) &&...
                    mh1(n+4,cell_num) < mh1(n,cell_num) && mh1(n+5,cell_num) < mh1(n,cell_num) &&...
                    mh1(n,cell_num) > mh1(n-2,cell_num)  && mh1(n,cell_num) > mh1(n-3,cell_num) &&...
                    mh1(n,cell_num) > mh1(n-4,cell_num)  && mh1(n,cell_num) > mh1(n-5,cell_num))
                
                tmaxpenult = tmaxlast;
                tmaxlast = Time(n);
                mmaxlast = mh1(n,cell_num);
                current_period = tmaxlast-tmaxpenult;
                periods = [periods current_period];
                current_amplitude = mmaxlast-mminlast;
                amplitudes = [amplitudes current_amplitude];
            end
            
        end
        % check if the current point is a trough
        %checking that it's smaller than 5 before and 5 after
        if (mh1(n+1,cell_num) > mh1(n,cell_num) && mh1(n,cell_num) < mh1(n-1,cell_num))
            
            if (mh1(n+2,cell_num) > mh1(n,cell_num) && mh1(n+3,cell_num) > mh1(n,cell_num) &&...
                    mh1(n+4,cell_num) > mh1(n,cell_num) && mh1(n+5,cell_num) > mh1(n,cell_num) &&...
                    mh1(n,cell_num) < mh1(n-2,cell_num)  && mh1(n,cell_num) < mh1(n-3,cell_num) &&...
                    mh1(n,cell_num) < mh1(n-4,cell_num)  && mh1(n,cell_num) < mh1(n-5,cell_num))
                mminlast = mh1(n,cell_num);
            end
        end
    end
    
    num_periods = length(periods);
    num_amplitudes = length(amplitudes);
    
    period_array(cell_num) = mean(periods(4:num_periods));
    amplitude_array(cell_num) = mean(amplitudes(4:num_amplitudes));
end

period = mean(period_array);
amplitude = mean(amplitude_array);
end