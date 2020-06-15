function [period, amplitude]= findPeriodandAmplitude(mh1)

global time_steps;
    
    %  Calculates the oscillation features -- period, amplitude and peak to trough
%  ratio for a set of concentration levels.
%  The values are calculated using the last peak and trough of the oscillations,
%  since the amplitude of the first few oscillations can be slightly unstable.
%  For the wild type, the peak and trough at the middle of the graph are also calculated
%  in order to ensure that the oscillations are sustained.
    
    tmaxlast = 0; tmaxpenult = 0; mmaxlast = 0; mmaxpenult = 0;
	tminlast = 0; tminpenult = 0; mminlast = 0; mminpenult = 0;
		
    for n = 1:1:time_steps - 1
        % check if the current point is a peak
        if mh1(1,n + 1) < mh1(1,n) && mh1(1,n) > mh1(1,n-1)
            tmaxpenult = tmaxlast;
            tmaxlast = n*eps;
            mmaxpenult = mmaxlast;
            mmaxlast = mh1(1,n);
        end
        % check if the current point is a trough
        if (mh1(1,n+1) > mh1(1,n) && mh1(1,n) < mh1(1,n-1))
            tminpenult = tminlast;
            tminlast = n * eps;
            mminpenult = mminlast;
            mminlast = mh1(1,n);
        end
    end
    
	period = tmaxlast-tmaxpenult;
	amplitude = mmaxlast-mminlast;
end