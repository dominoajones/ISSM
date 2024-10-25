% Tell MATLAB where to find ISSM and external libraries by selecting
% "ISSM-earlydays" and adding subfolders to path.

% Verify that ISSM works by executing:
issmversion

% If ISSM is runnintg correctly, you will see this message: 
% % % Ice-sheet and Sea-level System Model (ISSM) Version 4.24
% % % (website: http://issm.jpl.nasa.gov contact: https://issm.ess.uci.edu/forum/)
% % % 
% % % Build date: Mon Nov 20 01:27:52 PST 2023
% % % Compiled on pc mingw64 x86_64 by jenkins
% % % 
% % % Copyright (c) 2009-2023 California Institute of Technology
% % % 
% % %     to get started type: issmdoc

cd('C:\Users\domjon\Desktop\ISSM-earlydays')

%%
% S = shaperead('trunk\KNS\KNS-AS\KNS-AS.shp');
% geoshow(S);
% exptool('KNS-setup\final-domain\domain.exp');
%   Delineate KNS catchment as contour, quit exptool and file will automatically be saved.
md=model;

% Name and Coordinate system
md.miscellaneous.name='KNS';
md.miscellaneous.notes='Initialization with final domain - basal friction inversion and model spin up';
md.mesh.epsg=3413;

%%
% Creating a mesh

md=triangle(model,'KNS-setup\final-domain\domain.exp',100000);
% plotmodel(md,'data','mesh');

disp(' Refining Mesh');
ncdatavx = 'trunk\datasets\Velocity\vel_vx.nc';
ncdatavy = 'Trunk\Data\Velocity\vel_vy.nc';

velx= ncread(ncdatavx,'Band1');
xx= ncread(ncdatavx,'x');
yx= ncread(ncdatavx,'y');

vely= ncread(ncdatavy,'Band1');
xy= ncread(ncdatavx,'x');
yy= ncread(ncdatavx,'y');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);

% Refine mesh using surface velocities as metric
md=bamg(md,'hmin',500,'hmax',25000,'field',vel,'err',5); % all values here are arbitrary
[md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,+1,51.72,64.18);

% plotmodel(md,'data','mesh');

% ### Notes on Mesh refinement from ISSM user guide ### 
% 1. 'field' followed by vel : the field we want to adapt the mesh to
% 2. err : the allowed interpolation error
% 3. hmin : minimum edge length
% 4. hmax : maximum edge length
% 5. gradation : ratio between two consecutive edges
% 6. Force the triangles to be equilateral by using the Fanisomax  option, which specifies the maximum level of anisotropy (between 0 and 1, 1 being fully isotropic).
% 7. You can further refine mesh in a specific region, see user guide pages 45-48 for instructions.

 md.private.bamg=struct();

%%
% plotmodel(md,'data','mesh','edgecolor','w');
% exptool('refinement.exp')

h=NaN*ones(md.mesh.numberofvertices,1);
in=ContourToNodes(md.mesh.x,md.mesh.y,'refinement.exp',1);
h(find(in))=350;
% plotmodel(md,'data',in,'edgecolor','w');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);

md=bamg(md,'hmin',500,'hmax',25000,'field',vel,'err',5,'hVertices',h);
% vel=shock(md.mesh.x,md.mesh.y);
plotmodel(md,'data','mesh');

%%
nctopg='trunk\datasets\Bedmachine\Clipped\KNS-clipped_BedMachine-topography.nc';
x2= ncread(nctopg,'x');
y2= ncread(nctopg,'y');
topg  = ncread(nctopg,'Band1')';
disp('   Interpolating bedrock topography');
md.geometry.base = InterpFromGridToMesh(x2,y2,topg,md.mesh.x,md.mesh.y,0);
md.geometry.bed=md.geometry.base;
%   Sanity check:
%plotmodel(md,'data',md.geometry.base)

ncsurf='Trunk\Data\Clipped\Bedmachine\KNS-clipped_BedMachine-surface.nc';
x3= ncread(ncsurf,'x');
y3= ncread(ncsurf,'y');
surf  = ncread(ncsurf, 'Band1')';
disp('   Interpolating surface elevation');
md.geometry.surface=InterpFromGridToMesh(x3,y3,surf,md.mesh.x,md.mesh.y,0);

md.geometry.thickness=md.geometry.surface-md.geometry.base;

%Set min thickness to 1 meter
pos0=find(md.geometry.thickness<=0);
md.geometry.thickness(pos0)=1;
md.geometry.surface=md.geometry.thickness+md.geometry.base;

md.inversion.thickness_obs=md.geometry.thickness;
md.inversion.surface_obs=md.geometry.surface;
%   Sanity check
%plotmodel(md,'data',md.geometry.surface)
% plotmodel(md,'data',md.geometry.thickness);

% Alternatively,
% ncthk='C:\Users\domjon\Desktop\ISSM_KNS\Trunk\Data\BedMachine\KNS-clipped_BedMachine-thickness.nc'
% x1= ncread(ncthk,'x');
% y1= ncread(ncthk,'y');
% thk= ncread(ncthk,'Band1');
% disp('   Interpolating thicknesses');
% md.geometry.thickness=InterpFromGridToMesh(x1,y1,thk',md.mesh.x,md.mesh.y,0);

% Extrude layers
%md=extrude(md,5,1); % number of layers here is arbitrary

% Define all elements as SSA:
md=setflowequation(md,'SSA','all');
%md=setflowequation(md,'HO','all');
% Consider more informative model ice flow physics later:
% SSA, SIA, Hybrid, Higher-order, Stokes
% https://www.antarcticglaciers.org/glaciers-and-climate/numerical-ice-sheet-models/hierarchy-ice-sheet-models-introduction/

% Set mask
% Assume all nodes are grounded
md=setmask(md,'','');

%save C:\Users\domjon\Desktop\ISSM_KNS\KNSGeometry md;

%%
% Parameters
disp('   Interpolating velocities');
% ncdatavv = 'Trunk\Data\Clipped\Velocity\rast.calc.vel.nc';
ncdatavv = 'Trunk\Data\Velocity\calculated-vel.nc';
% finfo = ncinfo(ncdatavv)
velv= ncread(ncdatavv,'Band1');
xv= ncread(ncdatavv,'x');
yv= ncread(ncdatavv,'y');
md.inversion.vel_obs = InterpFromGridToMesh(xv,yv,velv',md.mesh.x,md.mesh.y,0);
md.initialization.vel= md.inversion.vel_obs;

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
md.inversion.vx_obs  = vx;
md.inversion.vy_obs  = vy;
md.initialization.vx = md.inversion.vx_obs;
md.initialization.vy = md.inversion.vy_obs;
md.initialization.vz=zeros(md.mesh.numberofvertices,1);

disp('   Interpolating temperatures');
ncdatatemp='Trunk\Data\Clipped\racmo\t2m.kns-1958-2017.reprojected.3413.nc';
% finfo = ncinfo(ncdatatemp)
temp= ncread(ncdatatemp,'Band560'); %Aug. 2004
xt= ncread(ncdatatemp,'x');
yt= ncread(ncdatatemp,'y');
md.initialization.temperature=InterpFromGridToMesh(xt,yt,temp',md.mesh.x,md.mesh.y,0); %already in Kelvin
% Sanity check:
% plotmodel(md,'data',md.initialization.temperature)

disp('   Construct ice rheological properties'); %#md.materials
%n has one value per element ????????
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
%B has one value per vertex ????????
md.materials.rheology_B=paterson(md.initialization.temperature); %!!!!!!!!!!!!!!!!!!!!!!!!!!!
% md.materials.rheology_B=paterson((273-20)*ones(md.mesh.numberofvertices,1));
md.damage.D=zeros(md.mesh.numberofvertices,1);

% % %   From Nias et al. (2023), BedMachine v3 geometry and velocity is used to invert for basal fricion coefficient (citing Morlighem et al. 2010).

disp('   Inverse for basal friction parameters'); %#md.friction
friction_coefficient = 40;
md.friction.coefficient=friction_coefficient*ones(md.mesh.numberofvertices,1);
md.friction.p=ones(md.mesh.numberofelements,1);
md.friction.q=ones(md.mesh.numberofelements,1);

%no friction applied on floating ice
pos=find(md.mask.ocean_levelset<0);
md.friction.coefficient(pos)=0;
md.groundingline.migration='SubelementMigration';

md.inversion=m1qn3inversion();
md.inversion.vx_obs=md.initialization.vx;
md.inversion.vy_obs=md.initialization.vy;
md.inversion.vel_obs=md.initialization.vel;

% dirichlet boundary condition are known as SPCs: ice frozen to the base, no velocity	#md.stressbalance

% SPCs are initialized at NaN one value per vertex
	md.stressbalance.spcvx=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvy=NaN*ones(md.mesh.numberofvertices,1);
	md.stressbalance.spcvz=NaN*ones(md.mesh.numberofvertices,1);

md.stressbalance.referential = NaN*ones(md.mesh.numberofvertices, 6);
md.stressbalance.loadingforce = NaN*ones(md.mesh.numberofvertices, 3);

% Basal melting
% md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
% md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);

	% Control general
	md.inversion.iscontrol=1;
	md.inversion.maxsteps=40;
	md.inversion.maxiter=80;
	md.inversion.dxmin=1;
	md.inversion.gttol=1.0e-4;
	md.verbose=verbose('control',true);

	% Cost functions
	md.inversion.cost_functions=[101 103 501]; %md.inversion for the different cost functions
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=1;
	md.inversion.cost_functions_coefficients(:,2)=1;
    md.inversion.cost_functions_coefficients(:,3)=2e-6;

	% Controls
	md.inversion.control_parameters={'FrictionCoefficient'};
	md.inversion.min_parameters=0.001*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1); %!!!!!!!!!!!!!!!!!!!!!!

	% Additional parameters
	md.stressbalance.restol=0.01;
	md.stressbalance.reltol=0.1;
	md.stressbalance.abstol=NaN;

	% Solve
    md.miscellaneous.name='KNS';
	md.toolkits=toolkits;
	md.cluster=generic('name',oshostname,'np',2);
	md=solve(md,'Stressbalance');

	% Update model friction fields accordingly
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

	% plotmodel(md,'data',md.friction.coefficient)

    %Turn off inversion
   md.inversion.iscontrol=0;

cd('C:\Users\domjon\Desktop\KNS_refined_domain\')
% fid=fopen('C:\Users\domjon\Desktop\ISSM_KNS\Runs\KNS.bin','wb');
md.miscellaneous.name='KNS';
md.verbose=verbose(0);
md = solve(md,'Stressbalance');

%%
plotmodel(md,'axis#all','tight','data',md.friction.coefficient,'caxis',[0 100],'title','inverted friction coeff',...
	'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities', ...
    'data',md.initialization.vel,'title','observed velocities',...
    'data',(md.results.StressbalanceSolution.Vel-md.initialization.vel),'title','(modeled-observed) velocities',...
    'colorbar#all','on','colorbartitle#2-4','(m/yr)',...
	'caxis#2-3',([1.5,6000]),...
	'log#2-3',10);

%%

disp('      creating boundary conditions');
% %   % ### General ###
% % % help SetIceSheetBC
% % % help SetIceShelfBC
% % % help SetMarineIceSheetBC
% %   % SetIceSheetBC creates the boundary conditions for stressbalance and thermal models for an IceSheet with no Ice Front.
% %   % SetIceShelfBC: Neumann BC are used on the ice front (an ANRGUS contour around the ice front must be given in input). Dirichlet BC are used elsewhere for stressbalance
% %   % SetMarineIceSheetBC: Neumann BC are used on the ice front (an ARGUS contour around the ice front can be given in input, or it will be deduced as onfloatingice & onboundary). Dirichlet BC are used elsewhere for stressbalance

% %   % ### For now ###
% %   % At the beginning, simplify model with Dirichlet BC everywhere:
md=SetIceSheetBC(md);

% %   % ### For later ###
% % % Look at the research question and objectives. If appropriate, add Neumann BC at the front to track rate of change.
% % % Making the ice front:
% % % S = shaperead("C:\Users\domjon\Desktop\ISSM_KNS\Trunk\KNS-Catchment\kns_basin\kns_basin.shp");
% % % geoshow(S);
% % % exptool('C:\Users\domjon\Desktop\ISSM_KNS\Front.exp');
% % % Example of adjusted tidewater margin BC:
% % md=SetMarineIceSheetBC(md,'./Front.exp');
% % md.basalforcings.floatingice_melting_rate=zeros(md.mesh.numberofvertices,1);
% % md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);

% % %   For SMB, Nias et al. (2023) used a 30-year relaxation period with a 1960-1989 mean from RACMO2.3p2 

disp('   Interpolating surface mass balance');
ncdatasmb = 'C:\Users\domjon\Desktop\ISSM_KNS\Trunk\Data\RACMO23p2\smb-reprojected.nc';
% finfo = ncinfo(ncdatasmb)
smb= ncread(ncdatasmb, 'Band1');
% smb=smb(:,:,560); %Aug. 2004
xs= ncread(ncdatasmb,'x');
ys= ncread(ncdatasmb,'y');
md.smb.mass_balance=InterpFromGridToMesh(xs,ys,smb',md.mesh.x,md.mesh.y,0);
% This was in the tutorials, don't know if it's needed for RACMO data:
% md.smb.mass_balance=md.smb.mass_balance*md.materials.rho_water/md.materials.rho_ice;

md.stressbalance.maxiter=100;

% % %   Sanity check:
% % %plotmodel(md,'data',md.smb.mass_balance)

disp('Mass Transport parameters');
md.masstransport.stabilization= 4
%0: no stabilization, 1: artificial diffusion, 2: streamline upwinding 
%3: discontinuous Galerkin, 4: flux corrected transport (FCT), 5: streamline upwind Petrov-Galerkin (SUPG)
md.masstransport.spcthickness=4000*ones(md.mesh.numberofvertices,1);
%thickness constraints (NaN means no constraint)

disp('Calving Parameterization'); %following Cuzzone et al. 2022
md.calving=calvingvonmises();
md.calving.stress_threshold_groundedice=600000;
md.calving.stress_threshold_floatingice=200000;
md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);

md.levelset.reinit_frequency=1; %this bit is all me, baby
md.levelset.spclevelset=NaN*ones(md.mesh.numberofvertices,1);
md.levelset.migration_max=10000;
md.transient.ismovingfront=1;

%%
% % ############ BEFORE SOLVING TO-DOs###############
% % Make sure to manually change the path for the following files:
% % marshall.m
% % toolkits.m
% % generic.m
md.inversion.iscontrol=0;
md.miscellaneous.name='KNS'
fid=fopen('C:\Users\domjon\Desktop\KNS_refined_domain','wb');
% % Set cluster #md.cluster
% 	% generic parameters #help generic
% 	% set only the name and number of process
md.cluster=generic('name',oshostname(),'np',4);
% 	% Set which control message you want to see #help verbose
% md.verbose=verbose('convergence',true);
md.verbose=verbose('all'); 

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
md.timestepping.final_time=10;

md=solve(md,'Transient');

%% 

plotmodel(md, ...
    'data', md.inversion.vel_obs,'title','Observed Velocity',...
    'data',md.results.TransientSolution(1).Vel,'title','Modeled Velocity (Year 1)',...
    'data',md.results.TransientSolution(10).Vel,'title','Modeled Velocity (Year 10)',...
    'colorbar#all','on','colorbartitle#1-3','(m/yr)',...
	'caxis#1-2',([1.5,7000]),...
    'caxis#3',([1.5,100000]),...
	'log#1-3',10);

