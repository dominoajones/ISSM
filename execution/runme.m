md=model;
disp('Making mesh');
md.miscellaneous.name='KNS';
md.miscellaneous.notes='Setting up initialization';
md.mesh.epsg=3413;

% Creating a mesh
md=triangle(model,'/projects/domino/ISSM/execution/domain/domain/domain.exp',1000000);

disp(' --Refining Mesh');
ncdatavx = '/projects/domino/ISSM/execution/trunk/datasets/Velocity/Clipped/vx_vel-CL.nc';
ncdatavy = '/projects/domino/ISSM/execution/trunk/datasets/Velocity/Clipped/vy_vel-CL.nc' ;

velx= ncread(ncdatavx,'Band1');
xx= ncread(ncdatavx,'x');
yx= ncread(ncdatavx,'y');

vely= ncread(ncdatavy,'Band1');
xy= ncread(ncdatavx,'x');
yy= ncread(ncdatavx,'y');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);
in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/ISSM/execution/domain/domain/no-ice-mask.exp',1);
vx(find(in))=0;
vy(find(in))=0;
vel(find(in))=0;


% Refine mesh using surface velocities as metric
md=bamg(md,'hmin',500,'hmax',1000000,'field',vel,'err',5); 
[md.mesh.lat,md.mesh.long] = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

disp('  --Refine mesh at TW margin');
md.private.bamg=struct();

h=NaN*ones(md.mesh.numberofvertices,1);
in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/ISSM/execution/domain/domain/refinement.exp',1);
h(find(in))=500; % may be too fine
% in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/ISSM/execution/domain/domain/no-ice-mask.exp',1);
% h(find(in))=2500;
% % plotmodel(md,'data',in,'edgecolor','w');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);
in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/ISSM/execution/domain/domain/no-ice-mask.exp',1);
vx(find(in))=0;
vy(find(in))=0;
vel(find(in))=0;

md=bamg(md,'hmin',500,'hmax',5000,'field',vel,'err',5,'hVertices',h);
plotmodel(md,'data','mesh');
           hold on
gt = shaperead('/projects/domino/ISSM/execution/trunk/KNS/paleo-positions/paleo-positions.shp');
mapshow(gt);      
S = shaperead('/projects/domino/ISSM/execution/trunk/KNS/fjord.shp');
mapshow(S, 'LineWidth', 1.5, 'EdgeColor', 'none', 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);

clear h in ncdatavx ncdatavy vel velx vely vx vy xx xy yx yy %tidy

save /projects/domino/ISSM/execution/runs/mesh.mat md

%%
disp('Parameterization');
md = loadmodel('/projects/domino/ISSM/execution/runs/mesh.mat');
    md = setmask(md,'','');
	md = parameterize(md,'/projects/domino/ISSM/execution/KNS.par');
	md = setflowequation(md,'SSA','all');

    save /projects/domino/ISSM/execution/runs/parameterized.mat md

%%
disp('Setting up inversion');
md = loadmodel('/projects/domino/ISSM/execution/runs/parameterized.mat');
md.inversion=taoinversion();

md.groundingline.migration='SubelementMigration';

md.inversion.vx_obs=md.initialization.vx;
md.inversion.vy_obs=md.initialization.vy;
md.inversion.vel_obs=md.initialization.vel;
md.inversion.algorithm='blmvm'; % minimization algorithm: 'blmvm', 'cg', or 'lmvm'

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
in=ContourToNodes(md.mesh.x,md.mesh.y,'/projects/domino/ISSM/execution/domain/domain/refinement.exp',1);
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

    save /projects/domino/ISSM/execution/runs/inverted.mat md

    clear h in pos

    %%
    md.timestepping.start_time=1972;
    md.timestepping.final_time=2017;
    md.mask.ice_levelset=interpMonthlyIceMaskGreene(md.mesh.x, md.mesh.y, [md.timestepping.start_time, md.timestepping.final_time]);
    %%
    md.smb=SMBforcing();
    ncdatasmb = '/projects/domino/ISSM/execution/trunk/datasets/RACMO23p2/Clipped/SMB_clipped_3413.nc';
    % finfo = ncinfo(ncdatasmb)
    xs = ncread(ncdatasmb, 'x');
    ys = ncread(ncdatasmb, 'y');

SMB= [];
 
    for i = 168:708 % 1972 to 2017, note data is given in bands each month since 1958-01-15
        varName = sprintf('Band%d', i);                         % Construct the variable name (e.g., 'Band25', 'Band26', ...)
        SMBBand = ncread(ncdatasmb, varName); % Read the variable data
        md.smb.mass_balance = InterpFromGridToMesh(xs, ys, SMBBand', md.mesh.x, md.mesh.y, 0); % Perform interpolation
        md.smb.mass_balance=md.smb.mass_balance / 917 * 12; %Convert from mm/month WE to m/yr IE
        SMB = cat(2, SMB, md.smb.mass_balance); % Concatenate along the third dimension
    end

    x = [1972+1/12:1/12:2017+1/12];
    md.smb.mass_balance = NaN(md.mesh.numberofvertices+1,length(x));
    md.smb.mass_balance(1:end-1,:)=SMB;
    md.smb.mass_balance(end,:)=x;
