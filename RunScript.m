%%
ppn = 16;
totalProcs = 15;
%ClusterInfo.setQueueName('matlab');
%ClusterInfo.setProcsPerNode(ppn);
c = parcluster;

%%
job176 = c.batch(@maximizeVertSeg,4,{},'Pool',totalProcs,'AttachedFiles',{'/Users/vani/Documents/MATLAB/VertSeg/S20ProfAy'});


%%
job176_State = job176.State;
%job172_Diary = job171.diary;
job176_Output = job176.fetchOutputs;

%%


dlmwrite('VertSegSets.csv',job71_output{4},'delimiter',',','-append');