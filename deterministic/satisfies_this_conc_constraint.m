function isTrue = satisfies_this_conc_constraint(conc, upper_limit)
%returns a boolean that indicates whether the given concentration is below
%the given upper limit
global time_steps;
if length(conc) ~= time_steps
    isTrue = false;
    return;
end
%each conc has the dimensions: cells x time_steps

isTrue=max(conc(1, 10000:time_steps)) <= upper_limit;

if ~isTrue
    max(conc(1, 10000:time_steps))
end
%     for n = 10000:1:time_steps - 1
        % check if the current point is a peak
%         if conc(1,n + 1) < conc(1,n) && conc(1,n) > conc(1,n-1)
%             if conc(1,n) > upper_limit
%                 isTrue = false;
%                 return;
%             end     
%         end
        

%     end
end