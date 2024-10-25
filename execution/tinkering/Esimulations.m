
    md=loadmodel('pddforce.mat')

    %%
x = [0+1/12:1/12:1];

    pr=[];
    for i=1:12
        test= [];
        test= [precip(:,end-72+i), precip(:,end-60+i), precip(:,end-48+i), precip(:,end-36+i), precip(:,end-24+i), precip(:,end-12+i)];
        test= mean(test,2);
        pr= [pr, test];
    end

    md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
	md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,1*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
    md.smb.precipitations_reconstructed(1:end-1,:) = pr;
    md.smb.precipitations_reconstructed(end,:) = x;

    md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

t=[];
    for i=1:12
        test= [];
        test= [tmp(:,end-72+i), tmp(:,end-60+i), tmp(:,end-48+i), tmp(:,end-36+i), tmp(:,end-24+i), tmp(:,end-12+i)];
        test= mean(test,2);
        t= [t, test];
    end

    md.smb.istemperaturescaled=0;
	md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,1*12);
	md.smb.temperatures_reconstructed(1:end-1,:) = t;
    md.smb.temperatures_reconstructed(end,:) = x;

    md.transient.issmb=1;
    md.transient.ismasstransport=1;
    md.transient.isstressbalance=1;
    md.transient.isthermal=0;
    md.transient.isgroundingline=1;
    % md.transient.isgia=0;
    md.transient.isesa=0;
    md.transient.isdamageevolution=0;
    md.transient.ismovingfront=1;
    md.transient.ishydrology=0;
    % md.transient.isslr=0;
    md.transient.isoceancoupling=0;
    % md.transient.iscoupler=0;

    md.groundingline.migration='SubelementMigration';
    

    % Set timestepping options, run for 50 years, saving every year
md.timestepping = timesteppingadaptive(md.timestepping);
md.timestepping.final_time=250;
md.timestepping.time_step_min=0.0001;
md.settings.output_frequency=5; % output every Nth timestep
    % md.settings.recording_frequency=0; % frequency of recording for a restart

    %playing around with mass transport
md.masstransport.isfreesurface=0;
md.masstransport.stabilization=5;

    % Additional options
    md.inversion.iscontrol=0;
    md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

	% Solve
	md.toolkits=toolkits;
    md.verbose=verbose('all');
    %md.cluster=generic('name',oshostname,'np',2);
    md.cluster=barkla();

md.calving.calvingrate=zeros(md.mesh.numberofvertices,1);
md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);

loadonly = 0;

   %Make sure jobs are submitted without MATLAB waiting for job completion 
   md.settings.waitonlock = 0;
   md.cluster.interactive = 0; %only needed if you are using the generic cluster
   md.miscellaneous.name = 'pddtestts_250'; 

  %Submit job or download results, make sure that there is no runtime name (that includes the date)
   md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);


%%
   %Save model if necessary
loadonly = 1;

   %Make sure jobs are submitted without MATLAB waiting for job completion 
   md.settings.waitonlock = 0;
   md.cluster.interactive = 0; %only needed if you are using the generic cluster
   md.miscellaneous.name = 'pddsetup'; 

  %Submit job or download results, make sure that there is no runtime name (that includes the date)
   md=solve(md,'Stressbalance','runtimename',false,'loadonly',loadonly);
   md=loadresultsfromcluster(md);
   save('pddsmallts250.mat', 'md', '-v7.3');
%%
   md=loadmodel('pddsmallts250.mat')
   %%
   
    plotmodel(md,'data',md.results.TransientSolution(4456).Thickness,'title','modelled thickness', ...
        'data',md.geometry.thickness,'title','observed thickness',...
        'colorbar#all','on','colorbartitle#all','(m)',...
	    'caxis#all',([1.5,40000]),...
	    'log#all',10)

    rmse_val = calculate_rmse(md.results.TransientSolution(26).Vel, md.initialization.vel);
%%
    plotmodel(md,'data',md.mask.ice_levelset,'title','Ice', ...
        'colorbar#all','on','caxis#all',([-1,1]));

    %%
    custom_cmap = [linspace(1, 0, 256)', linspace(1, 1, 256)', linspace(1, 1, 256)'];
colormap(custom_cmap);
axis equal
hold on

gt = shaperead('./trunk/KNS/paleo-positions/paleo-positions.shp');
mapshow(gt);
S = shaperead('./trunk/KNS/fjord.shp');
mapshow(S, 'LineWidth', 1.5, 'EdgeColor', 'none', 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);

colorbar;

% Add and customize the color bar
c = colorbar;
%%
% Set the background color of the figure to 'none'
set(gcf, 'Color', 'none');

% Optional: You can also set the axis background to 'none' if needed
set(gca, 'Color', 'none');

export_fig('plot.png', '-png', '-transparent');


    %rmse_val = calculate_rmse(md.results.TransientSolution(26).Thickness, md.geometry.thickness);
   
