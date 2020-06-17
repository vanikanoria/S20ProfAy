function [amplitude] = checkAmplitude(sol)%provide vector of solution @ which to eval amplitude
    peaks=0;
    troughs=0;
    value_peaks=0;
    value_troughs=0;
    for z=2:size(sol,1)-1
        if sol(z-1)<sol(z)&&sol(z)>sol(z+1)%Find peaks.
            peaks=peaks+1;
            if peaks>2
                value_peaks=value_peaks+sol(z);
            end
        end
        if sol(z-1)>sol(z)&&sol(z)<sol(z+1)%Find troughs.
            troughs=troughs+1;
            if troughs>2
                value_troughs=value_troughs+sol(z);
            end
        end
    end
    if peaks>=3 && troughs>=3 %Must have at least 3 peaks and at least 3 troughs
        average_peaks=value_peaks/(peaks-2);
        average_troughs=value_troughs/(troughs-2);
        amplitude=(average_peaks-average_troughs)/2;
    else
        amplitude = 0;
    end
end