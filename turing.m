oldpath=path 
path('~/Documents/MATLAB/SupportPackages/R2018a/parallel/pbs/nonshared', oldpath);

email='vkanoria@colgate.edu'; 
clusterHost='turing.colgate.edu'; 
jobStorageLocation='/Users/vani/Documents/MATLAB/VertSeg/S20ProfAy'; 
remoteJobStorageLocation='/home/vkanoria/Documents';

cluster=parallel.cluster.Generic( 'JobStorageLocation', jobStorageLocation);

set(cluster, 'HasSharedFilesystem', false); 
set(cluster, 'ClusterMatlabRoot', '/shared/MATLAB/R2018a'); 
set(cluster, 'OperatingSystem', 'unix'); 
set(cluster, 'IntegrationScriptsLocation', '~/Documents/MATLAB/SupportPackages/R2018a/parallel/pbs/nonshared'); 
set(cluster, 'NumWorkers', 32);

cluster.AdditionalProperties.ClusterHost=clusterHost; 
cluster.AdditionalProperties.RemoteJobStorageLocation=remoteJobStorageLocation; 
cluster.AdditionalProperties.QueueName = 'matlab'; 
cluster.AdditionalProperties.EmailAddress = email; 
cluster.AdditionalProperties.UseUniqueSubfolders = true;

saveAsProfile(cluster, 'Turing'); 
parallel.defaultClusterProfile('Turing');