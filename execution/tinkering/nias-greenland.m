% Script for running the SSA initialisation:
%       - Step 1: stress balance inversion with initial guess derived 
%               from driving stress 
%       - Step 2: Relaxation: with initialisation date is 2007 (Dec?), 
%       even though the climate is from RACMO 1960-1989, other data used 
%       in the initialisation is from ~2007. Run for 50 years.
%       - Step 3: Estimate floating ice basal melt rate, based on the



if any(steps==1)
        disp('   Step 1: M1QN3 control method friction -- initial guess from driving stress');

   md.inversion=m1qn3inversion();

    % Control general
    md.inversion.iscontrol=1;
    md.inversion.maxsteps=20;   %%%% increase? in Helenes code nsteps=300
    md.inversion.maxiter=40;
    md.inversion.dxmin=0.1;
    md.inversion.gttol=0.0001;
	
    % Cost functions
    md.inversion.cost_functions=[101 103 501];
    md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
    md.inversion.cost_functions_coefficients(:,1)= 200;     % 100; 
    md.inversion.cost_functions_coefficients(:,2)= 1;       % 1;
    md.inversion.cost_functions_coefficients(:,3)= 1e-7;    % 1e-8;
    
    % exclude regions where velocity is 0.0 m/yr from cost functions
    pos = find(md.inversion.vel_obs==0.0);
    md.inversion.cost_functions_coefficients(pos,1:2) = 0;
    
    disp('   Initialize basal friction using driving stress');
    disp('      -- Compute surface slopes and use 10 L2 projections');
    [sx,sy,s]=slope(md); sslope=averaging(md,s,10);
    disp('      -- Process surface velocity data');
    vel = md.inversion.vel_obs;
    flags=(vel==0); pos1=find(flags); pos2=find(~flags);
    vel(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),vel(pos2),md.mesh.x(pos1),md.mesh.y(pos1));
    %velmax = max(vel);
    %vel(vel==0 & md.mask.ice_levelset<0) = velmax;
    %disp('      -- Filling in missing ice velocity with MEaSUREs mosaic');
    %[velx, vely] = interpJoughinCompositeGreenland(md.mesh.x,md.mesh.y);
    %vel = sqrt( velx.^2 + vely.^2 );
    %idx = md.mask.ice_levelset < 0 & isnan(vel);
    %vel(idx) = sqrt( velx(idx).^2 + vely(idx).^2 );
    vel=max(vel,0.1);
    disp('      -- Calculate effective pressure');
    Neff = md.materials.rho_ice*md.geometry.thickness+md.materials.rho_water*md.geometry.base;
    Neff(find(Neff<=0))=1;
    disp('      -- Deduce friction coefficient');
    md.friction.coefficient=sqrt(md.materials.rho_ice*md.geometry.thickness.*(sslope)./(Neff.*vel/md.constants.yts));
    md.friction.coefficient=min(md.friction.coefficient,200);
    md.friction.p = 1.0 * ones(md.mesh.numberofelements,1);
    md.friction.q = 1.0 * ones(md.mesh.numberofelements,1);
    disp('      -- Extrapolate on ice free and floating ice regions');
    flags=(md.mask.ice_levelset>0) | (md.mask.ocean_levelset<0); pos1=find(flags); pos2=find(~flags);
    %md.friction.coefficient(pos1) = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2),...
    %md.mesh.x(pos1),md.mesh.y(pos1),'natural');
    md.friction.coefficient(pos1) = 1;
    pos=find(isnan(md.friction.coefficient));
    md.friction.coefficient(pos)  = 1;
    
    %%%% test regularization term (L-curve test with uob-umod)

    % Controls on inverted values
    md.inversion.control_parameters={'FrictionCoefficient'};
    md.inversion.min_parameters=1*ones(md.mesh.numberofvertices,1);
    md.inversion.max_parameters=200*ones(md.mesh.numberofvertices,1);

    % Additional parameters
    % % For stress balance
    md.stressbalance.restol=0.01; 
    md.stressbalance.reltol=0.1; 
    md.stressbalance.abstol=NaN;

    % Go solve
    md.miscellaneous.name='KNS';
	md.toolkits=toolkits;
	md.cluster=generic('name',oshostname,'np',2);
    md.verbose=verbose(1);
	md=solve(md,'Stressbalance');

	% Update model friction fields accordingly:
          % Budd
	md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;

	% plotmodel(md,'data',md.friction.coefficient)

    %Turn off inversion
   md.inversion.iscontrol=0;

cd('C:\Users\domjon\Desktop\KNS-ISSM\')
% fid=fopen('C:\Users\domjon\Desktop\ISSM_KNS\Runs\KNS.bin','wb');
md.miscellaneous.name='KNS';
md.verbose=verbose(0);
md = solve(md,'Stressbalance');

plotmodel(md,'axis#all','tight','data',md.friction.coefficient,'caxis',[0 100],'title','inverted friction coeff',...
	'data',md.results.StressbalanceSolution.Vel,'title','modeled velocities', ...
    'data',md.initialization.vel,'title','observed velocities',...
    'data',(md.results.StressbalanceSolution.Vel-md.initialization.vel),'title','(modeled-observed) velocities',...
    'colorbar#all','on','colorbartitle#2-4','(m/yr)',...
	'caxis#2-3',([1.5,6000]),...
	'log#2-3',10);

%% 

    % update initial velocity
    md.initialization.vel=md.results.StressbalanceSolution.Vel;
    md.initialization.vx=md.results.StressbalanceSolution.Vx;
    md.initialization.vy=md.results.StressbalanceSolution.Vy;
    
    % Save
    save ./models/gris.cmmtt.control_drag.ssa.sb md;

end



if any(steps==2)
        disp('   Step 2: Relaxation SSA (i.e. 50 year transient run)');

    % load friction inversion model
    md = loadmodel('./models/gris.cmmtt.control_drag.ssa.sb');

    % set which components of the transient solution to run
    md.transient.issmb=0;
    md.transient.ismasstransport=1;
    md.transient.isstressbalance=1;
    md.transient.isthermal=0;
    md.transient.isgroundingline=1;
    md.transient.isgia=0;
    md.transient.isesa=0;
    md.transient.isdamageevolution=0;
    md.transient.ismovingfront=0;
    md.transient.ishydrology=0;
    md.transient.isslr=0;
    md.transient.isoceancoupling=0;
    md.transient.iscoupler=0;

    md.groundingline.migration='SubelementMigration';
    
    % set flow equations
    md = setflowequation(md,'SSA','all');

    % Set timestepping options, run for 50 years, saving every year
    md.timestepping.time_step=0.05;%0.01; % need to adjust for CFL
    md.timestepping.start_time=0; %years
    md.timestepping.final_time=50; %years
    md.settings.output_frequency=20; % output every Nth timestep
    md.settings.recording_frequency=0; % frequency of recording for a restart

    % Additional options
    md.inversion.iscontrol=0;
    md.transient.requested_outputs={'default','IceVolumeAboveFloatation','IceVolume','TotalSmb','SmbMassBalance'}; % other options?
    
    % Solve
    md.cluster=clusters_gsfc('oibserve');
    md.cluster.np=8;
    md.miscellaneous.name='gris_ssa_relax';
    md.cluster.interactive=0; %runs in background on cluster (adds & to end of *.queue)
    md.toolkits=toolkits; 
    md.verbose=verbose('solution',true,'control',true,'convergence',true);
    md.settings.waitonlock=0; % Model results must be loaded manually with md=loadresultsfromcluster(md);
    md=solve(md,'Transient');%,'restart','gris_ssa-03-14-2019-15-51-23-60987');

    % Paterson
    %md=loadresultsfromcluster(md,'runtimename','gris_ssa_relax-10-13-2019-14-50-30-754');
    % Cuffey
    md=loadresultsfromcluster(md,'runtimename','gris_ssa_relax-10-18-2019-13-18-07-65142');
    
    plotmodel(md,'data',(md.results.TransientSolution(10).Thickness - md.results.TransientSolution(10-5).Thickness)/5,...
        'data',(md.results.TransientSolution(20).Thickness - md.results.TransientSolution(20-5).Thickness)/5,...
        'data',(md.results.TransientSolution(30).Thickness - md.results.TransientSolution(30-5).Thickness)/5,...
        'data',(md.results.TransientSolution(40).Thickness - md.results.TransientSolution(40-5).Thickness)/5,...
        'data',(md.results.TransientSolution(50).Thickness - md.results.TransientSolution(50-5).Thickness)/5,...
        'caxis#all',[-10 10])
    
    test = zeros(50,1);
    for i=2:51
        test(i-1) = sum((md.results.TransientSolution(i).Thickness - md.results.TransientSolution(i-1).Thickness).^2);
    end
    plot(test)
    
    % From this analysis, I think a 50 year relaxation longer than
    % necessary. Therefore I will save a shorter version (30 years)
    md.results.TransientSolution(32:51) = [];
    save ./models/gris.cmmtt.relax30yr.ssa.smb6089 md;
    
    % Save
    %save ./models/gris.cmmtt.relax50yr.ssa.smb6089 md;

end

 

if any(steps==3)
        disp('   Step 3: Final step of initialisation: finding floating ice basal melt rate)');

    % load friction inversion model
    %md = loadmodel('./models/gris.cmmtt.relax50yr.ssa.smb6089');
    md = loadmodel('./models/gris.cmmtt.relax30yr.ssa.smb6089');

    % find modeled dh/dt for final 5 years of the relaxation
    %dhdt_mod = (md.results.TransientSolution(51).Thickness - md.results.TransientSolution(46).Thickness)/5;
    dhdt_mod = (md.results.TransientSolution(31).Thickness - md.results.TransientSolution(26).Thickness)/5;
    
    % subtract smb from dhdt_mod
    dhdt_mod_dyn = dhdt_mod - md.smb.mass_balance;
    
    % Load in 03-09 mean dhdt_dyn data (observed)
    ncdata='/Users/inias/Documents/science/data/Greenland_dhdt/fromBea/Icedyndhdtave0309nc.nc';
    x = ncread(ncdata,'X');
    y = ncread(ncdata,'Y');
    dhdt = ncread(ncdata,'Uniform Lattice #0');
    dhdt(dhdt<-1000) = 0; % to remove the fill value of -1e32
    % % get md vertices lat lon locations onto the netcdf's projection
    [xi,yi]= ll2xy(md.mesh.lat,md.mesh.long,+1,39,71);
    % % interp data from grid to mesh
    dhdt_dyn_obs = InterpFromGridToMesh(x,y,dhdt,xi,yi,0);
    clear x y xi yi dhdt ncdata;

    % find dhdt_obs - dhdt_mod (for dynamic component only)
    % % for multi-year average dhdt_obs (2003-2009)
    diff = dhdt_dyn_obs - dhdt_mod_dyn;
    % % negative diff in floating areas mean that the smb isn't enough to
    % % account for the thinning, and so account for this difference with
    % % basal melt.
    % % set all positive values to zero (probably not required, but don't want to force ice shelf accretion) 
    diff(diff>0)=0;
    % for floating ice, we can assume this difference can be accounted for
    % in the basal melt rate (positive if melting) [m/yr] so change in sign
    % needed
    md.basalforcings.floatingice_melting_rate = -diff;
    
    % reinitialise model from end of the relaxation
    md = transientrestart(md);
    % remove inversion results...
    md.results = rmfield(md.results,'StressbalanceSolution');
    % remove relaxation results...
    md.results = rmfield(md.results,'TransientSolution3');
    % these results fields are removed to save space and speed up loading
    % of the initial state when configuring the ensemble
    
    % Save
    save ./models/gris.cmmtt.initialstate.ssa.smb6089 md;
    
    
    

end

