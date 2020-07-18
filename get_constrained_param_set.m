for i = 1:length(BestSets)

if (n>10000) && (mh1(i,n) > 100 || mh7(i,n) > 100 ||md(i,n) > 100)
    mh1=0;
    continue;
elseif (n>10000) && (ph1(i,n) > 1000 ||ph7(i,n) > 1000 || pd(i,n) > 1000)
    %             disp('protein too low:');
    mh1=0;
    continue;
elseif (n>10000) && (ph11(i,n) > 1000 ||ph17(i,n) > 1000 ||ph76(i,n) > 1000 ...
        || ph77(i,n) > 1000 || ph66(i,n) > 1000|| ph16(i,n) > 1000)
    %             disp('dimer level too low:');
    mh1=0;
    continue;
end

end