% IMPORTANT IMPORTANT IMPORTANT
% Please make "KNS-ISSM" you're working directory

%%
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
%%
md=model;
% Name and Coordinate system
md.miscellaneous.name='KNS';
md.miscellaneous.notes='Setting up initialization';
md.mesh.epsg=3413;

% Creating a mesh
md=triangle(model,'./execution/domain/domain.exp',1000000);

disp(' Refining Mesh');
ncdatavx = './execution/trunk/datasets/Velocity/Clipped/vx_vel-CL.nc';
ncdatavy = './execution/trunk/datasets/Velocity/Clipped/vy_vel-CL.nc' ;

velx= ncread(ncdatavx,'Band1');
xx= ncread(ncdatavx,'x');
yx= ncread(ncdatavx,'y');

vely= ncread(ncdatavy,'Band1');
xy= ncread(ncdatavx,'x');
yy= ncread(ncdatavx,'y');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/KNS-setup/domain/no-ice-mask.exp',1);
vx(find(in))=0;
vy(find(in))=0;
vel(find(in))=0;


% Refine mesh using surface velocities as metric
md=bamg(md,'hmin',500,'hmax',1000000,'field',vel,'err',5); 
[md.mesh.lat,md.mesh.long] = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

disp('Refine mesh at TW margin');
md.private.bamg=struct();

h=NaN*ones(md.mesh.numberofvertices,1);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/refinement.exp',1);
h(find(in))=500; % may be too fine
% in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/no-ice-mask.exp',1);
% h(find(in))=2500;
% % plotmodel(md,'data',in,'edgecolor','w');

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
vel	= sqrt(vx.^2+vy.^2);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/no-ice-mask.exp',1);
vx(find(in))=0;
vy(find(in))=0;
vel(find(in))=0;

md=bamg(md,'hmin',500,'hmax',50000,'field',vel,'err',5,'hVertices',h);
% % vel=shock(md.mesh.x,md.mesh.y);
plotmodel(md,'data','mesh');

disp('Interpolate geometry');
nctopg='./execution/trunk/datasets/Bedmachine/Clipped/bed-bedmachine-CL.nc';
x2= ncread(nctopg,'x');
y2= ncread(nctopg,'y');
topg  = ncread(nctopg,'bed')';
disp('   Bedrock topography');
md.geometry.base = InterpFromGridToMesh(x2,y2,topg,md.mesh.x,md.mesh.y,0);
md.geometry.bed=md.geometry.base;
%   Sanity check:
%plotmodel(md,'data',md.geometry.base)

ncsurf='./execution/trunk/datasets/Bedmachine/Clipped/surface-bedmachine-CL.nc';
x3= ncread(ncsurf,'x');
y3= ncread(ncsurf,'y');
surf  = ncread(ncsurf, 'surface')';
disp('   Surface elevation');
md.geometry.surface=InterpFromGridToMesh(x3,y3,surf,md.mesh.x,md.mesh.y,0);
in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/no-ice-mask.exp',1);
md.geometry.surface(find(in))=md.geometry.base(find(in));

md.geometry.thickness=md.geometry.surface-md.geometry.base;
md.inversion.thickness_obs=md.geometry.thickness;

%Set min thickness to 10 meters
% pos0=find(md.geometry.thickness<=10);
% md.geometry.thickness(pos0)=10;
% md.geometry.surface(find(in))=0;
% md.geometry.surface=md.geometry.thickness+md.geometry.base;

% % CHANGE FOR NO ICE ZONES TO HAVE 0 THICKNESS (IMPORTANT FOR TRANSIENT
% RUN)

md.inversion.surface_obs=md.geometry.surface;

disp('Interpolate other parameters');
disp('    Velocity');
ncdatavv = './execution/trunk/datasets/Velocity/Clipped/cal_vel-CL.nc';
% finfo = ncinfo(ncdatavv)
% novel=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/no-velocity-values.exp',1);
velv= ncread(ncdatavv,'Band1');
xv= ncread(ncdatavv,'x');
yv= ncread(ncdatavv,'y');
md.inversion.vel_obs = InterpFromGridToMesh(xv,yv,velv',md.mesh.x,md.mesh.y,0);
md.inversion.vel_obs(find(in))=0;
% md.inversion.vel_obs(find(novel))=0;
md.initialization.vel= md.inversion.vel_obs;

vx	= InterpFromGridToMesh(xx,yx,velx',md.mesh.x,md.mesh.y,0);
vy	= InterpFromGridToMesh(xy,yy,vely',md.mesh.x,md.mesh.y,0);
md.inversion.vx_obs  = vx;
md.inversion.vx_obs(find(in))=0;
% md.inversion.vx_obs(find(novel))=0;
md.inversion.vy_obs  = vy;
md.inversion.vy_obs(find(in))=0;
% md.inversion.vy_obs(find(novel))=0;
md.initialization.vx = md.inversion.vx_obs;
md.initialization.vy = md.inversion.vy_obs;
md.initialization.vz=zeros(md.mesh.numberofvertices,1);

disp('   Temperature');
ncdatatemp=./execution/trunk/datasets/RACMO23p2/Clipped/t2m_clipped.nc;
% finfo = ncinfo(ncdatatemp);
xtemp= ncread(ncdatatemp,'x');
ytemp= ncread(ncdatatemp,'y');
temp_month = [];
for i = 25:385 % 30 year SMB average from 1960 to 1990, note data is given in bands each month since 1958-01-15
    varName = sprintf('Band%d', i); 
    tempBand = ncread(ncdatatemp, varName); 
    md.initialization.temperature = InterpFromGridToMesh(xtemp, ytemp, tempBand', md.mesh.x, md.mesh.y, 0);
    temp_month = cat(3, temp_month, md.initialization.temperature); 
end
% Calculate the average value for each cell across all time steps
md.initialization.temperature = mean(temp_month, 3);

disp('   Construct ice rheological properties'); %md.materials
md.materials.rheology_n=3*ones(md.mesh.numberofelements,1);
md.materials.rheology_B=cuffey(md.initialization.temperature); 
md.damage.D=zeros(md.mesh.numberofvertices,1);
md.materials.rheology_law='None';
% law for the temperature dependance of the rheology: 'None', 'BuddJacka', 'Cuffey', 'CuffeyTemperate', 'Paterson', 'Arrhenius', 'LliboutryDuval', 'NyeH2O', or 'NyeCO2'

md=setflowequation(md,'SSA','all');

disp('Set mask');
in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/no-ice-mask.exp',1);
md=setmask(md,'','');
md.mask.ice_levelset(find(in))=1;

disp('Set boundary conditions');
% md=SetMarineIceSheetBC(md,'./execution/domain/front.exp'); 
md=SetIceSheetBC(md);

pos = find(md.mesh.vertexonboundary);
md.masstransport.spcthickness(pos) = md.geometry.thickness(pos);

%Make mass transport more stable
md.masstransport.stabilization = 1; %5: SUPG, 2: SU, 1: art diff

in=ContourToNodes(md.mesh.x,md.mesh.y,'./execution/domain/adjustBC.exp',1);
pos=find(in);
md.stressbalance.spcvx(pos)=NaN;
md.stressbalance.spcvy(pos)=NaN;
md.stressbalance.spcvz(pos)=NaN;

%%
plotmodel(md,'data',md.stressbalance.spcvx,'shading','flat','nan',0);

%%
plotmodel(md,'data',md.inversion.vx_obs,'caxis',[-160,0]);
