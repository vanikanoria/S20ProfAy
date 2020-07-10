function isTrue = satisfies_all_conc_constraints(param_set)
isTrue=true;
mRNA_upper_limit = 100;
protein_upper_limit = 1000;
dimer_upper_limit = 1000;
upper_limit_array = [mRNA_upper_limit, mRNA_upper_limit,mRNA_upper_limit, ...
    protein_upper_limit, protein_upper_limit, protein_upper_limit, ...
    dimer_upper_limit, dimer_upper_limit, dimer_upper_limit, ...
    dimer_upper_limit, dimer_upper_limit, dimer_upper_limit];

conc_array = deterministic_model_extended(param_set);
% = [mh1,mh7,md,ph1,ph7,ph6,ph11,ph17,ph16,ph77,ph66,ph76]

num_mRNA = 3;
num_protein = 3;
num_dimers = 6;
total_num = num_mRNA + num_protein + num_dimers;

for i = 1:total_num
    conc = conc_array(i);
    upper_limit = upper_limit_array(i);
    isTrue = conc_condition_satisfied(conc, upper_limit);
    if ~isTrue
        return;
    end
end

