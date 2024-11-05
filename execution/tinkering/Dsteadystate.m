
%%
        % Basal friction
in = ContourToNodes(md.mesh.x, md.mesh.y, './execution/domain/refinement.exp', 1) & md.friction.coefficient>1.5;
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

in=ContourToNodes(md.mesh.x, md.mesh.y, './execution/domain/smoothing.exp', 1) ;
% md.friction.coefficient(find(in))=smoothdata(md.friction.coefficient(find(in)),'SmoothingFactor',0.2);
% md.friction.coefficient(find(in))=smoothdata(md.friction.coefficient(find(in)),'rlowess',4);
md.friction.coefficient(find(in))=smoothdata(md.friction.coefficient(find(in)),"movmedian",10);
plotmodel(md,'data',md.friction.coefficient,'title','Friction Coefficient')
    %%
    md.smb = SMBd18opdd();
    x = [0+1/12:1/12:100+1/12];
	a = x(1:1:end-1);

    pr=[];
    for i=1:12
        test= [];
        test= [precip(:,0+i), precip(:,12+i), precip(:,24+i), precip(:,36+i), precip(:,48+i), precip(:,60+i)];
        test= mean(test,2);
        pr= [pr, test];
    end

    pr1=[];
    for i=1:100
        pr1=[pr1 pr];
    end

    md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
	md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,100*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
    md.smb.precipitations_reconstructed(1:end-1,:) = pr1;
    md.smb.precipitations_reconstructed(end,:) = a;

    md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

t=[];
    for i=1:12
        test= [];
        test= [tmp(:,0+i), tmp(:,12+i), tmp(:,24+i), tmp(:,36+i), tmp(:,48+i), tmp(:,60+i)];
        test= mean(test,2);
        t= [t, test];
    end

    t1=[];
    for i=1:100
        t1=[t1 t];
    end

    md.smb.istemperaturescaled=0;
	md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,100*12);
	md.smb.temperatures_reconstructed(1:end-1,:) = t1;
    md.smb.temperatures_reconstructed(end,:) = a;
    md.smb.rlaps=0;

%     %md.smb = SMBd18opdd(); 
%     n=168;
%     m=180;
%     x = [0+1/12:1/12:100+1/12];
% 	a = x(1:1:end-1);
% 
%     pr=[];
%     for i=1:12
%         test= [];
%         test= [precip(:,n+i), precip(:,m+i)];
%         test= mean(test,2);
%         pr= [pr, test];
%     end
% 
%     pr1=[];
%     for i=1:100
%         pr1=[pr1 pr];
%     end
% 
%     md.smb.isprecipscaled=0; %  This allows us to set the precip ourselves to what we have done above.
% 	md.smb.precipitations_reconstructed = NaN(md.mesh.numberofvertices+1,100*12); % Needs to be size numberofvertices+1, so the last row can hold the timestep
%     md.smb.precipitations_reconstructed(1:end-1,:) = pr1;
%     md.smb.precipitations_reconstructed(end,:) = a;
%     md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.
% 
% t=[];
%     for i=1:12
%         test= [];
%         test= [tmp(:,n+i), tmp(:,m+i)];
%         test= mean(test,2);
%         t= [t, test];
%     end
% 
%     t1=[];
%     for i=1:100
%         t1=[t1 t];
%     end
% 
%     md.smb.istemperaturescaled=0;
% 	md.smb.temperatures_reconstructed= NaN(md.mesh.numberofvertices+1,100*12);
% 	md.smb.temperatures_reconstructed(1:end-1,:) = t1;
%     md.smb.temperatures_reconstructed(end,:) = a;
%     md.smb.rlaps=0;

% md.smb=SMBforcing();
% ncdatasmb = './execution/trunk/datasets/RACMO23p2/Clipped/SMB_clipped_3413.nc';
% % finfo = ncinfo(ncdatasmb)
% xs = ncread(ncdatasmb, 'x');
% ys = ncread(ncdatasmb, 'y');
% 
% SMB_month = [];
% 
% for i = 1:268 % 30 year SMB average from 1960 to 1990, note data is given in bands each month since 1958-01-15
%     varName = sprintf('Band%d', i);                         % Construct the variable name (e.g., 'Band25', 'Band26', ...)
%     SMBBand = ncread(ncdatasmb, varName); % Read the variable data
%     md.smb.mass_balance = InterpFromGridToMesh(xs, ys, SMBBand', md.mesh.x, md.mesh.y, 0); % Perform interpolation
%     SMB_month = cat(3, SMB_month, md.smb.mass_balance); % Concatenate along the third dimension
% end
% 
% % Calculate the average value for each cell across all time steps
% md.smb.mass_balance = mean(SMB_month, 3);
% 
% md.smb.mass_balance=md.smb.mass_balance / 917 * 12 %Convert from mm/month WE to m/yr IE

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
    md.timestepping.final_time=100;
    %md.timestepping.final_time=md.timestepping.start_time+100;
    md.timestepping.time_step_min=0.01;
    md.timestepping.time_step_max=1/12;
    md.settings.output_frequency=5;

    % 
    % % Set timestepping options, run for 50 years, saving every year
    % md.timestepping.time_step=0.05;%0.01; % need to adjust for CFL
    % md.timestepping.start_time=0; %years
    % md.timestepping.final_time=100; %years
    % md.settings.output_frequency=20; % output every Nth timestep
    % md.settings.recording_frequency=0; % frequency of recording for a restart

    %playing around with mass transport
md.masstransport.isfreesurface=0;

md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);
md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);

md.masstransport.stabilization = 1;

    % Additional options
    md.inversion.iscontrol=0;
    md.transient.requested_outputs={'default','IceVolume','TotalSmb','SmbMassBalance'}; % other options?

	% Solve
	md.toolkits=toolkits;
    md.verbose=verbose('all');
    md.cluster=barkla();

    %%

    md.calving=calvinghab();
    md.calving.flotation_fraction=4;
%%
    % reinitialise model from end of the relaxation
    md = transientrestart(md,600);
    % remove inversion results...
    %md.results = rmfield(md.results,'StressbalanceSolution');
    % remove relaxation results...
    md.results = rmfield(md.results,'TransientSolution3');

    md.results=rmfield(md.results,'TransientSolution');
%%

name='calv_9'
loadonly = 0;

   %Make sure jobs are submitted without MATLAB waiting for job completion 
   md.settings.waitonlock = 0;
   md.cluster.interactive = 0; %only needed if you are using the generic cluster
   md.miscellaneous.name = name; 

  %Submit job or download results, make sure that there is no runtime name (that includes the date)
   md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);

%% PASTE IN GIT BASH, then:
% scp domjon@barkla6.liv.ac.uk:/mnt/data1/users/domjon/ISSM/execution/calv_9/calv_9.outbin C:/Users/domjon/Desktop/KNS-ISSM
   md=loadresultsfromdisk(md,'calv_9.outbin')

%%
   surfmb=[];
   pos = find(md.mask.ice_levelset<0);
for i=1:24;
surfmb=[surfmb md.results.TransientSolution(i).SmbMassBalance(pos)];
end
%Plot surface mass balance time series in first subplot
plot([1:24],mean(surfmb));
%Title this plot Mean surface mass balance
title('Transient mass balance from PDD forcing');
ylabel('Mass Balance (m/yr IE)')
xlabel('Time step')
   %%
   tt=[];
   pos = find(test.mask.ice_levelset<0);
for i=1:24;
tt=[tt test.results.TransientSolution(i).SmbMassBalance(pos)];
end
%Plot surface mass balance time series in first subplot
plot([1:24],mean(tt));
%Title this plot Mean surface mass balance
title('Transient mass balance from PDD forcing');
ylabel('Mass Balance (m/yr IE)')
xlabel('Time step')

   %%
plotmodel(md,'data',md.geometry.base);
z = [-1000:10:1000];demcmap(z)
   %%
       plotmodel(md,'data',md.results.TransientSolution(end).Vel,'title','modelled thickness', ...
        'colorbar#all','on','colorbartitle#all','(m)',...
	    'caxis#all',([1.5,6000]),'log#all',10)

       hold on
gt = shaperead('./execution/trunk/KNS/paleo-positions/paleo-positions.shp');
mapshow(gt);
S = shaperead('./execution/trunk/KNS/fjord.shp');
mapshow(S, 'LineWidth', 1.5, 'EdgeColor', 'none', 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);

%%
plotmodel(md,'data','transient_movie','caxis#all',([1.5,4000]))
   %%
    % reinitialise model from end of the relaxation
    md = transientrestart(md);
    % remove inversion results...
    md.results = rmfield(md.results,'StressbalanceSolution');
    % remove relaxation results...
    md.results = rmfield(md.results,'TransientSolution3');
    save ./execution/tinkering/Shared/steadystate.mat md
%%
    plotmodel(md,'data',md.results.TransientSolution(18).Vel,'title','modelled velocities', ...
        'data',md.initialization.vel,'title','observed velocities',...
        'colorbar#all','on','colorbartitle#all','(m/yr)',...
	    'caxis#all',([1.5,6000]),'log#all',10)

    rmse_val = calculate_rmse(md.results.TransientSolution(26).Vel, md.initialization.vel);
%%
    plotmodel(md,'data',md.results.TransientSolution(700).Vel,'title','modelled thickness', ...
	    'caxis#all',([1.5,6000]),'log#all',10)

    % rmse_val = calculate_rmse(md.results.TransientSolution(26).Thickness, md.geometry.thickness);

    %%
    figure('WindowState', 'maximized')
    plotmodel(md,'data','transient_movie','transient_movie_field','Vel','log#all',10,'caxis#all',([1.5,6000]),'transient_movie_output','smb_NOsmooth_STILL_thick.gif')
        %plotmodel(md,'data','transient_movie','transient_movie_field','Vx','transient_movie_output','600-900_part2_vx.gif')

%            hold on
% gt = shaperead('./execution/trunk/KNS/paleo-positions/paleo-positions.shp');
% mapshow(gt);      
% S = shaperead('./execution/trunk/KNS/fjord.shp');
% mapshow(S, 'LineWidth', 1.5, 'EdgeColor', 'none', 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.5);

%%
plotmodel(md,'data',md.smb.mass_balance,'title','Mass Balance')
