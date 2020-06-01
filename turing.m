oldpath=path 
path('~/Documents/MATLAB/SupportPackages/R2018a/parallel/pbs/nonshared', oldpath);

email='vkanoria@colgate.edu'; 
clusterHost='turing.colgate.edu'; 
jobStorageLocation='/Users/vani/Desktop/VertSeg/S20ProfAy'; 
remoteJobStorageLocation='/home/vkanoria/Somewhere on the Turing Cluster';

cluster=parallel.cluster.Generic( 'JobStorageLocation', jobStorageLocation);

set(cluster, 'HasSharedFilesystem', false); 
set(cluster, 'ClusterMatlabRoot', '/shared/MATLAB/R2018a'); 
set(cluster, 'OperatingSystem', 'unix'); 
set(cluster, 'IntegrationScriptsLocation', '~/Documents/MATLAB/SupportPackages/R2018a/parallel/pbs/nonshared'); 
set(cluster, 'NumWorkers', 8);

cluster.AdditionalProperties.ClusterHost=clusterHost; 
cluster.AdditionalProperties.RemoteJobStorageLocation=remoteJobStorageLocation; 
cluster.AdditionalProperties.QueueName = 'matlab'; 
cluster.AdditionalProperties.EmailAddress = email; 
cluster.AdditionalProperties.UseUniqueSubfolders = true;

saveAsProfile(cluster, 'Turing'); 
parallel.defaultClusterProfile('Turing');