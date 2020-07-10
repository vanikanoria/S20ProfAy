function [period, amplitude]= findPeriodandAmplitude(mh1)

global time_steps;
global eps;
eps=0.01;

%use a sliding average to smoothen the data
mh1(1,:) = movmean(mh1(1,:),40);

    
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
        if (mh1(1,n + 1) < mh1(1,n) && mh1(1,n) > mh1(1,n-1))
            
            if (mh1(1,n + 2) < mh1(1,n) && mh1(1,n + 3) < mh1(1,n) &&...
                    mh1(1,n + 4) < mh1(1,n) && mh1(1,n + 5) < mh1(1,n) &&...
                    mh1(1,n) > mh1(1,n-2)  && mh1(1,n) > mh1(1,n-3) &&...
                    mh1(1,n) > mh1(1,n-4)  && mh1(1,n) > mh1(1,n-5))
                
                tmaxpenult = tmaxlast;
                tmaxlast = n*eps;
                mmaxlast = mh1(1,n);
                period = tmaxlast-tmaxpenult;
                periods = [periods period];
                amplitude = mmaxlast-mminlast;
                amplitudes = [amplitudes amplitude];
            end
            
        end
        % check if the current point is a trough
        %checking that it's smaller than 5 before and 5 after
        if (mh1(1,n + 1) > mh1(1,n) && mh1(1,n) < mh1(1,n-1))
            
            if (mh1(1,n + 2) > mh1(1,n) && mh1(1,n + 3) > mh1(1,n) &&...
                    mh1(1,n + 4) > mh1(1,n) && mh1(1,n + 5) > mh1(1,n) &&...
                    mh1(1,n) < mh1(1,n-2)  && mh1(1,n) < mh1(1,n-3) &&...
                    mh1(1,n) < mh1(1,n-4)  && mh1(1,n) < mh1(1,n-5))
                mminlast = mh1(1,n);
            end
        end
    end
    num_periods = length(periods);
    num_amplitudes = length(amplitudes);
    
    period = mean(periods(4:num_periods));
    amplitude = mean(amplitudes(4:num_amplitudes));
end