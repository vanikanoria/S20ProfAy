function isTrue = satisfiesWTconditions(mh1)
global time_steps;
 % Calculates the oscillation feature -- peak to trough
%  ratio for a set of concentration levels.
%  The values are calculated using the last peak and trough of the oscillations,
%  since the amplitude of the first few oscillations can be slightly unstable.
%  For the wild type, the peak and trough at the middle of the graph are also calculated
%  in order to ensure that the oscillations are sustained.
    
    mmaxlast = 0; % tmaxlast = 0; tmaxpenult = 0; mmaxpenult = 0;
    mminlast = 0; % tminlast = 0; tminpenult = 0; mminpenult = 0;
		
    for n = 1:1:time_steps - 1
        % Check if the current point is a peak
        if mh1(1,n + 1) < mh1(1,n) && mh1(1,n) > mh1(1,n-1)
            % tmaxpenult = tmaxlast;
            % tmaxlast = n*eps;
            % mmaxpenult = mmaxlast;
            mmaxlast = mh1(1,n);
        end
        % Check if the current point is a trough
        if (mh1(1,n+1) > mh1(1,n) && mh1(1,n) < mh1(1,n-1))
            % tminpenult = tminlast;
            % tminlast = n * eps;
            % mminpenult = mminlast;
            mminlast = mh1(1,n);
        end
    end
    
% Calculate the peak and trough at the middle of the graph

        mmaxlast2 = 0.0; mminlast2 = 0.0;
        for m = 2:1:time_steps/2
            if(mh1(1,m + 1) < mh1(1,m) && mh1(1,m) > mh1(1,m - 1))
                mmaxlast2 = mh1(1,m);
            end
            if(mh1(1,m + 1) > mh1(1,m) && mh1(1,m) < mh1(1,m - 1))
                mminlast2 = mh1(1,m);
            end
        end
        % In order to avoid dividing by zero in case a trough is 0, set it to 1
        if(mminlast2 == 0.0 || mminlast == 0.0)
            mminlast2 = 1.0;
            mminlast = 1.0;
        end
        
        peaktotrough1 = mmaxlast/mminlast; 
        peaktotrough2 = mmaxlast2/mminlast2;
        
        isTrue = (peaktotrough2 >= 1.5 && peaktotrough1 >= 1.5 && (peaktotrough2 / peaktotrough1)<=1.5);
end