disp('Initial guess');

disp('   Initialize basal friction using driving stress');
    disp('      -- Compute surface slopes and use 10 L2 projections');
    [sx,sy,s]=slope(md); sslope=averaging(md,s,10);
    disp('      -- Process surface velocity data');
    vel = md.inversion.vel_obs;
    flags=(vel==0); pos1=find(flags); pos2=find(~flags);
    vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
    vel=max(vel,0.1);
    disp('      -- Calculate effective pressure');
    Neff = md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base;
    Neff(find(Neff<=0))=1;
    disp('      -- Deduce friction coefficient');
    md.friction.coefficient=sqrt(md.materials.rho_ice*md.geometry.thickness.*(sslope)./(Neff.*vel/md.constants.yts));
    md.friction.coefficient=min(md.friction.coefficient,150);
    md.friction.p = 1.0 * ones(md.mesh.numberofelements,1);
    md.friction.q = 1.0 * ones(md.mesh.numberofelements,1);
    disp('      -- Extrapolate on ice free and floating ice regions');
    flags=(md.mask.ice_levelset>0) | (md.mask.ocean_levelset<0); pos1=find(flags); pos2=find(~flags);
    %md.friction.coefficient(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2),...
    %md.mesh.x(pos1),md.mesh.y(pos1),'natural');
    md.friction.coefficient(pos1) = 1;
    pos=find(isnan(md.friction.coefficient));
    md.friction.coefficient(pos)  = 1;

% friction_coefficient = 25;
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

disp('Inversion parameters);')
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
in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/refinement.exp',1);
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

	% Solve
	md.toolkits=toolkits;
    md.verbose=verbose('all');
    md.cluster=generic('name',oshostname(),'np',2);
    md.miscellaneous.name='invert';
   md=solve(md,'Stressbalance');

	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
    save invert md
%%

    % update initial velocity
    md.initialization.vel=md.results.StressbalanceSolution.Vel;
    md.initialization.vx=md.results.StressbalanceSolution.Vx;
    md.initialization.vy=md.results.StressbalanceSolution.Vy;

