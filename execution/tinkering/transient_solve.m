% % ############ BEFORE SOLVING TO-DOs###############
% % Make sure to manually change the path for the following files:
% % marshall.m
% % toolkits.m
% % generic.m
md.inversion.iscontrol=0;
md.miscellaneous.name='KNS'
fid=fopen('C:\Users\domjon\Desktop\ISSM_KNS\Runs\1.bin','wb');
% % Set cluster #md.cluster
% 	% generic parameters #help generic
% 	% set only the name and number of process
md.cluster=generic('name',oshostname(),'np',4);
% 	% Set which control message you want to see #help verbose
md.verbose=verbose('convergence',true);

% 	% ########### Solve #help solve ##########
% 	% we are solving a StressBalance
% 	% md=solve(md,'Stressbalance');
% % plotmodel(md,'data',md.results.StressbalanceSolution.Vel);
% 
% % Change the working directory:
% % pwd
% cd('C:\Users\domjon\Desktop\ISSM_KNS\Runs')


% % ########## Transient solve ###########
% % set the transient model to ignore the thermal model
% % #md.transient
md.transient.isthermal=0;
% define the timestepping scheme
% everything here should be provided in years #md.timestepping
% give the length of the time_step
md.timestepping.time_step=1;
% give final_time 
md.timestepping.final_time=100;

md=solve(md,'Transient');

%% 

plotmodel(md, ...
    'data', md.inversion.vel_obs,'title','Observed Velocity',...
    'data',md.results.TransientSolution(1).Vel,'title','Modeled Velocity (first)',...
    'data',md.results.TransientSolution(15).Vel,'title','Modeled Velocity (middle)',...
    'data',md.results.TransientSolution(33).Vel,'title','Modeled Velocity (last)',...
    'colorbar#all','on','colorbartitle#1-4','(m/yr)');
    % 'caxis#1-4',([1.5,]),...
    % 'log#1-4',10


