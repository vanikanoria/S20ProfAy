%%
ppn = 16;
totalProcs = 15;
%ClusterInfo.setQueueName('matlab');
%ClusterInfo.setProcsPerNode(ppn);
c = parcluster;

%%
job171 = c.batch(@maximizeVertSeg,4,{},'Pool',totalProcs,'AttachedFiles',{'/Users/vani/Documents/MATLAB/VertSeg/S20ProfAy'});

%%
job171_State = job171.State;
%job171_Diary = job171.diary;
job171_Output = job171.fetchOutputs;

%%


dlmwrite('VertSegSets.csv',job71_output{4},'delimiter',',','-append');