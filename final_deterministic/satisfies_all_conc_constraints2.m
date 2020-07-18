function isTrue = satisfies_all_conc_constraints2(param_set, mRNA_limit,protein_limit,dimer_limit, pd_limit)
%calls satisfies_this_conc_constraint(conc, upper_limit) for mRNA, proteins and dimers
%returns a boolean that indicates whether or not the concentrations of the
%above are below the upper limit
global time_steps;
isTrue=true;
mRNA_upper_limit = mRNA_limit;
protein_upper_limit = protein_limit;
dimer_upper_limit = dimer_limit;
pd_upper_limit = pd_limit;
upper_limit_array = [mRNA_upper_limit, mRNA_upper_limit,mRNA_upper_limit, ...
    protein_upper_limit, protein_upper_limit, pd_upper_limit, ...
    dimer_upper_limit, dimer_upper_limit, dimer_upper_limit, ...
    dimer_upper_limit, dimer_upper_limit, dimer_upper_limit];

[mh1,mh7,md,ph1,ph7,pd,ph11,ph17,ph16,ph77,ph66,ph76] = deterministic_model_extended(param_set);
conc_array = [mh1,mh7,md,ph1,ph7,pd,ph11,ph17,ph16,ph77,ph66,ph76];
num_mRNA = 3;
num_protein = 2;
%we count pd separately
num_dimers = 6;
total_num = num_mRNA + num_protein + num_dimers + 1;

for i = 1:total_num
    conc = conc_array(1,(i-1)*time_steps+1:(i-1)*time_steps+time_steps);
    upper_limit = upper_limit_array(i);
    isTrue = satisfies_this_conc_constraint(conc, upper_limit);
    if ~isTrue
        
        i
        return;
    end
end

