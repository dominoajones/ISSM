md=model;
disp('Making mesh');
md.miscellaneous.name='KNS';
md.miscellaneous.notes='Setting up initialization';
md.mesh.epsg=3413;

% Creating a mesh
md=triangle(model,'./domain/domain.exp',1000000);

disp(' --Refining Mesh');
ncdatavx = './trunk/datasets/Velocity/Clipped/vx_vel-CL.nc';
ncdatavy = './trunk/datasets/Velocity/Clipped/vy_vel-CL.nc' ;

velx= ncread(ncdatavx,'Band1');
xx= ncread(ncdatavx,'x');
yx= ncread(ncdatavx,'y');

vely= ncread(ncdatavy,'Band1');
xy= ncread(ncdatavx,'x');
yy= ncread(ncdatavx,'y');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./domain/no-ice-mask.exp',1);
vx(find(in))=0;
vy(find(in))=0;
vel(find(in))=0;


% Refine mesh using surface velocities as metric
md=bamg(md,'hmin',500,'hmax',1000000,'field',vel,'err',5); 
[md.mesh.lat,md.mesh.long] = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

disp('  --Refine mesh at TW margin');
md.private.bamg=struct();

h=NaN*ones(md.mesh.numberofvertices,1);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./domain/refinement.exp',1);
h(find(in))=500; % may be too fine
% in=ContourToNodes(md.mesh.x,md.mesh.y,'./domain/no-ice-mask.exp',1);
% h(find(in))=2500;
% % plotmodel(md,'data',in,'edgecolor','w');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./domain/no-ice-mask.exp',1);
vx(find(in))=0;
vy(find(in))=0;
vel(find(in))=0;

md=bamg(md,'hmin',500,'hmax',50000,'field',vel,'err',5,'hVertices',h);
plotmodel(md,'data','mesh');

clear h in ncdatavx ncdatavy vel velx vely vx vy xx xy yx yy %tidy

save ./runs/mesh.mat md

%%
disp('Parameterization');
md = loadmodel('./runs/mesh.mat');
    md = setmask(md,'','');
	md = parameterize(md,'./KNS.par');
	md = setflowequation(md,'SSA','all');

    save ./runs/parameterized.mat md

%%
disp('Setting up inversion');
md = loadmodel('./runs/parameterized.mat');
md.inversion=taoinversion();

md.groundingline.migration='SubelementMigration';

md.inversion.vx_obs=md.initialization.vx;
md.inversion.vy_obs=md.initialization.vy;
md.inversion.vel_obs=md.initialization.vel;
md.inversion.algorithm='lmvm'; % minimization algorithm: 'blmvm', 'cg', or 'lmvm'

md.stressbalance.loadingforce = NaN*ones(md.mesh.numberofvertices, 3);
md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);

md.stressbalance.isnewton=0;
md.initialization.pressure=md.constants.g*md.materials.rho_ice*md.geometry.thickness;

disp('Inversion parameters');
	md.inversion.iscontrol=1;
	md.inversion.maxsteps=300;
	md.inversion.maxiter=600;
	md.inversion.fatol=0; %convergence criterion: f(X)-f(X*) (X: current iteration, X*: "true" solution, f: cost function)
    md.inversion.frtol=0; %convergence criterion: |f(X)-f(X*)|/|f(X*)|
    md.inversion.gatol=0; %convergence criterion: ||g(X)|| (g: gradient of the cost function)
    md.inversion.grtol=0; %convergence criterion: ||g(X)||/|f(X)|
    md.inversion.gttol=0.0000001; %convergence criterion: ||g(X)||/||g(X0)|| (g(X0): gradient at initial guess X0)
	md.verbose=verbose('control',true);
    md.inversion.thickness_obs=md.geometry.thickness;
    md.inversion.surface_obs=md.geometry.surface;

	% Cost functions
   % 101: SurfaceAbsVelMisfit
   % 102: SurfaceRelVelMisfit
   % 103: SurfaceLogVelMisfit
   % 104: SurfaceLogVxVyMisfit
   % 105: SurfaceAverageVelMisfit
	md.inversion.cost_functions=[101 103 104 501]; 
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,4);
    md.inversion.cost_functions_coefficients(:,1)=200;
	md.inversion.cost_functions_coefficients(:,2)=1;
    md.inversion.cost_functions_coefficients(:,3)=1;
    md.inversion.cost_functions_coefficients(:,4)=1e-7;

    %changing cost function coefficients where mesh is refined
h=NaN*ones(md.mesh.numberofvertices,1);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./domain/refinement.exp',1);
h(find(in))=1;
pos=find(h==1);
md.inversion.cost_functions_coefficients(pos,1)=400;
md.inversion.cost_functions_coefficients(pos,2)=15;
md.inversion.cost_functions_coefficients(pos,3)=15;

    % exclude regions where velocity is 0.0 m/yr from cost functions
    pos = find(md.inversion.vel_obs==0.0);
    md.inversion.cost_functions_coefficients(pos,1:2) = 0;

	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.min_parameters=0.1*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=150*ones(md.mesh.numberofvertices,1); 

    md.stressbalance.restol=0.01; 
    md.stressbalance.reltol=0.1; 
	md.stressbalance.abstol=NaN;

    disp('Running inversion');
	% Solve
	md.toolkits=toolkits;
    md.verbose=verbose('all');
    md.cluster=generic('name',oshostname(),'np',2);
    md.miscellaneous.name='invert';
    md=solve(md,'Stressbalance');

	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

    save ./runs/inverted.mat md

    %%
    disp('Relaxation period initialization');
    md = loadmodel('./runs/inverted.mat');
    
disp('  --Regress basal friction coefficient');
in = ContourToNodes(md.mesh.x, md.mesh.y, './domain/refinement.exp', 1) & md.friction.coefficient>1.5;
indepbase = md.geometry.base(find(in));
depfriction = md.friction.coefficient(find(in));

% Perform linear regression  on the filtered data
p = polyfit(indepbase, depfriction, 1);
fittedValues = polyval(p, indepbase);

% Plot the data and the regression line
% figure;
% scatter(indepbase, depfriction, 'b'); % Scatter plot of the filtered data
% hold on;
% plot(indepbase, fittedValues, 'r', 'LineWidth', 2); % Regression line
% hold off;

% coefficient  =0.0634*elevation + 63.45
in = md.mask.ice_levelset>0;
md.friction.coefficient(find(in))= p(1)*md.geometry.base(find(in)) + p(2);

in=ContourToNodes(md.mesh.x, md.mesh.y, './domain/smoothing.exp', 1) ;
% md.friction.coefficient(find(in))=smoothdata(md.friction.coefficient(find(in)),'SmoothingFactor',0.2);
% md.friction.coefficient(find(in))=smoothdata(md.friction.coefficient(find(in)),'rlowess',4);
md.friction.coefficient(find(in))=smoothdata(md.friction.coefficient(find(in)),"movmedian",10);

disp('  --Tranisent run setup');
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

    md.groundingline.migration='None';

    md.timestepping = timesteppingadaptive(md.timestepping);
    md.timestepping.final_time=50;
    md.timestepping.time_step_min=1/12;
    md.timestepping.time_step_max=1/12;
    md.settings.output_frequency=5;

    % Additional options
    md.inversion.iscontrol=0;
    md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

	% Solve
	md.toolkits=toolkits;
    md.verbose=verbose('all');
    md.cluster=barkla();

disp('  --Other initialization parameters');
    md.calving=calvinghab();
    md.calving.flotation_fraction=0.8;
    md.masstransport.isfreesurface=0;

    md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
    md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
    md.masstransport.stabilization = 1;

disp ('   --SOLVE');
name='steadystate'
loadonly = 0;

   md.settings.waitonlock = 0;
   md.cluster.interactive = 0; %only needed if you are using the generic cluster
   md.miscellaneous.name = name; 

   md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

   clear depfriction fittedValues in indepbase loadonly name p

%%
disp('Save result');
md=loadresultsfromdisk(md,'steadystate.outbin');

disp('  --Restart after relaxation period');
    % reinitialise model from end of the relaxation
    md = transientrestart(md);
    % remove inversion results...
    md.results = rmfield(md.results,'StressbalanceSolution');
    % remove relaxation results...
    md.results = rmfield(md.results,'TransientSolution3');

    save ./runs/steadystate.mat md

    %%
    disp('Extrude model test');
md = loadmodel('./runs/steadystate.mat');
   md=extrude(md,5,1);
   md=setflowequation(md,'MOLHO','all');
   md.stressbalance.spcvx_shear = nan(md.mesh.numberofvertices,1);
   md.stressbalance.spcvy_shear = nan(md.mesh.numberofvertices,1);
   md.stressbalance.spcvx_base = nan(md.mesh.numberofvertices,1);
   md.stressbalance.spcvy_base = nan(md.mesh.numberofvertices,1);

   md.smb.pddfac_snow = 4.3*ones(md.mesh.numberofvertices,1); % Notice these are spatially constant, but if you wanted could be spatially varying since they are on the vertices.
   md.smb.pddfac_ice = 8.3*ones(md.mesh.numberofvertices,1);
%%
    md.calving.flotation_fraction=0;
   disp ('   --SOLVE');
   name='test'
   loadonly = 0;
   %md.timestepping.final_time=md.timestepping.start_time+50;
   md.settings.waitonlock = 0;
   md.cluster.interactive = 0; %only needed if you are using the generic cluster
   md.miscellaneous.name = name; 
   md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

   %%
   md=loadresultsfromdisk(md,'test.outbin');

   figure('WindowState', 'maximized')
    plotmodel(md,'data','transient_movie','transient_movie_field','Vel','log#all',10,'caxis#all',([1.5,6000]),'transient_movie_output','smb_NOsmooth_STILL_thick.gif')
