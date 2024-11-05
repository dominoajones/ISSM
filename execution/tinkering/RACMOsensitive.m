md=loadmodel('./tinkering/Shared/steadystate.mat')

in = ContourToNodes(md.mesh.x, md.mesh.y, '.\KNS-setup\domain\refinement.exp', 1) & md.friction.coefficient>1.5;
indepbase = md.geometry.base(find(in));
depfriction = md.friction.coefficient(find(in));

% Perform linear regression on the filtered data
p = polyfit(indepbase, depfriction, 1);
fittedValues = polyval(p, indepbase);

% Step 2: Calculate R-squared
SSresid = sum((depfriction - fittedValues).^2);
SStotal = sum((depfriction - mean(depfriction)).^2);
R_squared = 1 - SSresid/SStotal;

% Step 3: Calculate p-value and statistical significance
n = length(indepbase); % Number of data points
x = indepbase;
y = depfriction;

% Standard error of the regression
s_err = sqrt(SSresid / (n-2));

% Standard error of the slope
x_mean = mean(x);
SE_slope = s_err / sqrt(sum((x - x_mean).^2));

% t-statistic for the slope
t_stat = p(1) / SE_slope;

% p-value from t-statistic
p_value = 2 * (1 - tcdf(abs(t_stat), n-2)); % two-tailed p-value

% Display results
disp(['R-squared: ', num2str(R_squared)]);
disp(['p-value: ', num2str(p_value)]);
if p_value < 0.05
    disp('The result is statistically significant at the 5% significance level.');
else
    disp('The result is not statistically significant at the 5% significance level.');
end

%%

X = [ones(length(indepbase), 1), indepbase]; % Add a column of ones for the intercept
[b,~,~,~,stats] = regress(depfriction, X);

% Step 2: Extract R-squared and p-value
R_squared = stats(1); % R-squared value
F_stat = stats(2);    % F-statistic
p_value = stats(3);   % p-value
% Note: stats(4) contains the estimate of the error variance

% Step 3: Display results
disp(['R-squared: ', num2str(R_squared)]);
disp(['p-value: ', num2str(p_value)]);
if p_value < 0.05
    disp('The result is statistically significant at the 5% significance level.');
else
    disp('The result is not statistically significant at the 5% significance level.');
end
%%

% Plot the data and the regression line
% figure;
% scatter(indepbase, depfriction, 'b(1)'); % Scatter plot of the filtered data
% hold on;
% plot(indepbase, fittedValues, 'b(2)', 'LineWidth', 2); % Regression line
% hold off;
% Step 2: Calculate fitted values
fittedValues = X * b; % Compute the fitted values using the regression coefficients

% Step 3: Create the scatter plot of the original data
figure; % Create a new figure
scatter(indepbase, depfriction,10, 'b', 'filled'); % Scatter plot in blue with filled circles
hold on; % Keep the scatter plot so we can add the regression line

% Step 4: Plot the regression line
plot(indepbase, fittedValues, 'r-', 'LineWidth', 2); % Regression line in red with a thicker line width

%%

% coefficient  =0.0634*elevation + 63.45
in = md.mask.ice_levelset>0;
md.friction.coefficient(find(in))= p(1)*md.geometry.base(find(in)) + p(2);

in=ContourToNodes(md.mesh.x, md.mesh.y, '.\KNS-setup\domain\smoothing.exp', 1) ;
md.friction.coefficient(find(in))=smoothdata(md.friction.coefficient(find(in)),'SmoothingFactor',0.9);
plotmodel(md,'data',md.friction.coefficient)

in = md.mask.ice_levelset>0;
md.smb.mass_balance=0.4+md.smb.mass_balance;
md.smb.mass_balance(find(in))=0;
% plotmodel(md,'data',md.smb.mass_balance)

md.levelset.spclevelset=NaN(md.mesh.numberofvertices,1);
md.frontalforcings.meltingrate=zeros(md.mesh.numberofvertices,1);

%Better to make levelset a signed distance instead of ±1
md.mask.ice_levelset = reinitializelevelset(md, md.mask.ice_levelset);


%Set constraints to levelset equations (levelset cannot change along boundary)
pos = find(md.mesh.vertexonboundary);
md.levelset.spclevelset(pos) = md.mask.ice_levelset(pos);


md.calving.calvingrate=zeros(md.mesh.numberofvertices,1);


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
   

    % Additional options
    md.inversion.iscontrol=0;
	md.transient.requested_outputs={'IceVolume','TotalSmb','IceVolumeAboveFloatation', ...
		'SmbMassBalance'};


%Use adaptive time stepping to make sure the model is stable
md.timestepping = timesteppingadaptive(md.timestepping);
md.timestepping.final_time=100;



% 	md.amr.groundingline_resolution=500;
% 	md.amr.groundingline_distance=10000;
% 	md.amr.hmin=500; % the same resolution used around the grounding line
% 	md.amr.hmax=10000; % the same coase resolution used to generate the coarse mesh
% 	md.amr.gradation=1.7; % this controls the ratio between two consecutive edges
% 	md.amr.fieldname='None'; % no field used here
% 	md.amr.keepmetric=0; % no field, no metric
% md.transient.amr_frequency=1;
% md.levelset.fe='P1';

% md=extrude(md,5,1)
md=setflowequation(md,'SSA','all')

loadonly = 0;

   %Make sure jobs are submitted without MATLAB waiting for job completion 
   md.settings.waitonlock = 0;
   md.cluster.interactive = 0; %only needed if you are using the generic cluster
   md.miscellaneous.name='sensitive_SSA'; 

  %Submit job or download results, make sure that there is no runtime name (that includes the date)
   md=solve(md,'Transient','runtimename',false,'loadonly',loadonly);                         

%%
   %Save model if necessary
loadonly = 1;

   %Make sure jobs are submitted without MATLAB waiting for job completion 
   md.settings.waitonlock = 0;
   md.cluster.interactive = 0; %only needed if you are using the generic cluster
   md.miscellaneous.name = 'sensitive_SSA'; 

  %Submit job or download results, make sure that there is no runtime name (that includes the date)
   md=solve(md,'Stressbalance','runtimename',false,'loadonly',loadonly);
%save sensitive_SSA md
save sensitive_SSA md

%%
md=loadmodel('sensitive_HO.mat')
plotmodel(md,'data',md.results.TransientSolution(15).Vel,'title','Transient Model Velocity', ...
        'data',md.results.TransientSolution(15).Thickness,'title','Transient Model Thickness',...
        'data', md.smb.mass_balance,'title','SMB (RACMO plus 0.4 m/yr)',...
        'data',md.friction.coefficient,...
        'caxis#1',([0.5,6500]),'caxis#2',([0,2500]),'log#1',10,...
        'colorbar#all','on','colorbartitle#1-2','m','colorbartitle#3','m/yr ice equivelant');

%%
plotmodel(md,'data','amr',1,'title','t=1 yr','fontsize',12)
%%
    plotmodel(md,'data',md.results.TransientSolution(11).Vel,'title','Transient Model Velocity', ...
        'data',md.results.TransientSolution(11).Thickness,'title','Transient Model Thickness',...
        'data', md.smb.mass_balance,'title','SMB (RACMO plus 0.4 m/yr)',...
        'data',md.friction.coefficient,...
        'caxis#1',([0.5,6500]),'caxis#2',([0,2500]),'log#1',10,...
        'colorbar#all','on','colorbartitle#1-2','m','colorbartitle#3','m/yr ice equivelant');

    %%
    plotmodel(md,'data','transient_movie','caxis#all',([0,2500]))