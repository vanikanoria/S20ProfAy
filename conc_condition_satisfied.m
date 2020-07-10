function isTrue = satisfies_this_conc_condition(conc, upper_limit)
%returns a boolean that indicates whether the given concentration is below
%the given upper limit
isTrue=true;
global time_steps;
%eps=0.01
    
    for n = 10000:1:time_steps - 1
        % check if the current point is a peak
        if conc(1,n + 1) < conc(1,n) && conc(1,n) > conc(1,n-1)
            if conc(1,n) > upper_limit
                isTrue = false;
                return;
            end     
        end

    end
end