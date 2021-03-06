function [period, amplitude, sustainedOscillations]= checkPeriodAndAmplitude(mh1,is_wildtype)
    
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
        if (mh1(0,n+1) > mh1(0,n) && mh1(0,n) < mh1(0,n-1))
            tminpenult = tminlast;
            tminlast = n * eps;
            mminpenult = mminlast;
            mminlast = mh1(0,n);
        end
        
    end
    
    if (is_wildtype)
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