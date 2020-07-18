%%
ppn = 16;
totalProcs = 15;
%ClusterInfo.setQueueName('matlab');
%ClusterInfo.setProcsPerNode(ppn);
c = parcluster;

%%
job189 = c.batch(@maximizeVertSeg,4,{},'Pool',totalProcs,'AttachedFiles',{'/Users/vani/Documents/MATLAB/VertSeg/S20ProfAy/hybrid_code'});

%%
job189_State = job189.State;
%job171_Diary = job171.diary;
job189_Output = job189.fetchOutputs;

%%


dlmwrite('VertSegSets.csv',job71_output{4},'delimiter',',','-append');