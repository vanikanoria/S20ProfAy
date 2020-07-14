Set1 = [];
BetterSets=[];
global time_steps;
time_steps = 60000;
min_total = Inf;
total=0;

for i = 1:1001
   
    [mh1,mh7,md,ph1,ph7,pd,ph11,ph17,ph16,ph77,ph66,ph76] = deterministic_model_extended(SetEnd(i,:)');
conc_array = [mh1,mh7,md,ph1,ph7,pd,ph11,ph17,ph16,ph77,ph66,ph76];
num_mRNA = 3;
num_protein = 2;
%we count pd separately
num_dimers = 6;
total_num = num_mRNA + num_protein + num_dimers + 1;

for j = 1:total_num
    conc = conc_array(1,(j-1)*time_steps+1:(j-1)*time_steps+time_steps);
    
    total = total +  max(conc(1, 10000:time_steps));

    %find max conc of each of the states and add them for all the sets
    % the set with the least sum wins
end

if total < min_total
   
        min_total = total;
        BetterSets = [BetterSets ; SetEnd(i,:)];
        Set1 = SetEnd(i,:);
end
     disp(i);
end
