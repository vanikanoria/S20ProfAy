%%
ppn = 16;
totalProcs = 15;
%ClusterInfo.setQueueName('matlab');
%ClusterInfo.setProcsPerNode(ppn);
c = parcluster;

%%
job174 = c.batch(@maximizeVertSeg,4,{},'Pool',totalProcs,'AttachedFiles',{'/Users/vani/Documents/MATLAB/VertSeg/S20ProfAy'});


%%
job174_State = job174.State;
%job172_Diary = job171.diary;
job174_Output = job174.fetchOutputs;

%%


dlmwrite('VertSegSets.csv',job71_output{4},'delimiter',',','-append');