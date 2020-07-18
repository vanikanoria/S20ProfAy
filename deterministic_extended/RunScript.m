%%
ppn = 16;
totalProcs = 15;
%ClusterInfo.setQueueName('matlab');
%ClusterInfo.setProcsPerNode(ppn);
c = parcluster;

%%
job183 = c.batch(@maximizeVertSeg,4,{},'Pool',totalProcs,'AttachedFiles',{'/Users/vani/Documents/MATLAB/VertSeg/S20ProfAy'});


%%
job181_State = job181.State;
%job172_Diary = job171.diary;
job181_Output = job181.fetchOutputs;

%%


dlmwrite('BestSets2.csv',job71_output{4},'delimiter',',','-append');