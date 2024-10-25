md=loadmodel('invert.mat')

[md.mesh.lat,md.mesh.long] = xy2ll(md.mesh.x,md.mesh.y,+1,39,71);

% Set up SMB D180 class and set PDD factors
	
	md.smb = SMBd18opdd();  % Turn on Class
	md.smb.isd18opd=1;	% 1 means we want to load our monthly climatologies of temperature and precipitation 
	
	md.smb.delta18o = [-40.0110;0];  % This is meaningless, but needs to be set.  Basically, if you set md.smb.isd18opd=0, that would mean you want to take a modern climatology and scale it back through time based upon an ice core record.  For now, we just need to set this, but it will not be used since we set md.smb.isd18opd=1.
	
	md.smb.rlaps=6.0; % This is the Lapse rate you want to use
	md.smb.desfac = 1; % This is the elevation desertification effect.  0 is no change, 1 is maximum change.  Which is a reduction by a factor of 2 for every 1km change in surface elevation
	md.smb.rlapslgm = 5.5; % This just needs to be set, but is not used for your purposes
	
	% Set PDD Factors
	md.smb.issetpddfac = 1;
	md.smb.pddfac_snow = 4.3*ones(md.mesh.numberofvertices,1); % Notice these are spatially constant, but if you wanted could be spatially varying since they are on the vertices.
	md.smb.pddfac_ice = 8.3*ones(md.mesh.numberofvertices,1);

    disp('   ACC');
ncdataacc='./trunk/datasets/Box1840thr2012/accbox.nc';
finfo = ncinfo(ncdataacc);
xacc= ncread(ncdataacc,'x');
yacc= ncread(ncdataacc,'y');

md.smb.precipitations_presentday = NaN(md.mesh.numberofvertices,12);
q=1
for i = 133:144 % 12 months
    index = i:12:((i+11)+12*(149)); %index across 1850 to 2000 for each month
    acc_month=[];
    for w= 1:150
        varName = sprintf('Band%d', index(w)); 
        accBand = ncread(ncdataacc, varName); 
        acc = InterpFromGridToMesh(xacc, yacc, accBand', md.mesh.x, md.mesh.y, 0);
        acc_month = cat(3, acc_month, acc); 
    end
    md.smb.precipitations_presentday(:,q)=mean(acc_month,3);
    q=q+1;
end

    disp('   DEM');
ncdatadem='./trunk/datasets/Box1840thr2012/dembox.nc';
% finfo = ncinfo(ncdatadem);
xdem= ncread(ncdatadem,'x');
ydem= ncread(ncdatadem,'y');
dem=ncread(ncdatadem,'Band1');
box_dem=InterpFromGridToMesh(xdem, ydem, dem', md.mesh.x, md.mesh.y, 0);
md.smb.s0p = (max(box_dem,0));
md.smb.s0t = (max(box_dem,0));

    disp('   Temp');
ncdatatemp='./trunk/datasets/Box1840thr2012/tempbox.nc';
%finfo = ncinfo(ncdatatemp);
xtemp= ncread(ncdatatemp,'x');
ytemp= ncread(ncdatatemp,'y');

md.smb.temperatures_presentday = NaN(md.mesh.numberofvertices,12);
q=1
for i = 133:144 % 12 months
    index = i:12:((i+11)+12*(149)); %index across 1850 to 2000 for each month
    temp_month=[];
    for w= 1:150
        varName = sprintf('Band%d', index(w)); 
        tempBand = ncread(ncdatatemp, varName); 
        temp = InterpFromGridToMesh(xtemp, ytemp, tempBand', md.mesh.x, md.mesh.y, 0);
        temp_month = cat(3, temp_month, temp); 
    end
    md.smb.temperatures_presentday(:,q)=mean(temp_month,3);
    q=q+1;
end


ncbadg='./trunk/datasets/paleoran/briner2020recons.nc';
finfo = ncinfo(ncbadg);

	pre1 = ncread(ncbadg,'P_moderate');
    pre1 = permute(pre1,[4 3 2 1]);
    
	% Time for Jessica's recon go back to -21880 ka.  Index 373 is 650 CE, 379 is 950CE and index 399 is 1850CE
	pre1 = pre1(373:399,:,:,:);

	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;

	pre_jessica = NaN(md.mesh.numberofvertices+1,26,12);
	disp('Interpolating precipitation output onto grid');
	for i = 1:26
		for k = 1:12
			pre_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGridToMesh(latB,lonB',squeeze(pre1(i,k,:,:))',md.mesh.lat,md.mesh.long,0));
			%temp_jessica(end,:,k)=((k)/12);
		end
	end

	% For precipitation, we multiply the Box climatology by the fraction from the Badgeley Reconstruction
	precip = [];
	for i = 1:26
		for j = 1:12
			precip = [precip;md.smb.precipitations_presentday(:,j)'.*pre_jessica(1:end-1,i,j)'];
		end
	end

	
	% For the md.smb.precipitations_reconstructed, we need to set the last row to represent time.  I do this a weird way (I'm not sure why), but basically you need time to be increasing.  For example, I start my simulations at 12.5 ka years ago, which is represented by 112580 and go to the year 1850, which is represented by 124981.  Either way, time has to be increasing.  Also, for every 50 year chunk, we have monthly outputs.  Therefore, the time needs to have monthly steps.  That is the reason for the stride every 1/12.  	
	x = [0+1/12:1/12:1300+1/12];
	% Stride every 50 years - Time Resolution of UW Product
	a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

	xx = [];
	for i =1:26;
		xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
	end

	precip = permute(precip,[2 1]);
	
	md.smb.precipitations_reconstructed(md.smb.precipitations_reconstructed<0) = 0.01; % There is a possibility for negative precipitation. values.  Set anything less than 0 to low value.

	temp1 = ncread(ncbadg,'T_moderate');
    temp1 = permute(temp1,[4 3 2 1]);

	% Time for Jessica's recon go back to -21880 ka.  You will have to find your respective time steps, for me, 153 is ~12.5 ka and 401 is 1850CE
	temp1 = temp1(373:399,:,:,:);

	latB  = ncread(ncbadg,'lat');
	lon  = ncread(ncbadg,'lon');
	x = find(lon>180);
	lonB(x) = lon(x)-360;

	temp_jessica = NaN(md.mesh.numberofvertices+1,26,12);
	disp('Interpolating temperature output onto grid');
	for i = 1:26
		for k = 1:12
			temp_jessica(1:md.mesh.numberofvertices,i,k)=(InterpFromGridToMesh(latB,lonB',squeeze(temp1(i,k,:,:))',md.mesh.lat,md.mesh.long,0));
			%temp_jessica(end,:,k)=((k)/12);
		end
	end

	% For temperature, we add the Box climatology to the anomalies from the Badgeley Reconstruction
	tmp = [];
	for i = 1:26
		for j = 1:12
			tmp = [tmp;md.smb.temperatures_presentday(:,j)'+temp_jessica(1:end-1,i,j)'];
		end
	end

	x = [0+1/12:1/12:1300+1/12];

	% Stride every 50 years - Time Resolution of UW Product
	a = x(1:600:end);b = x(2:600:end);c = x(3:600:end);d = x(4:600:end);e = x(5:600:end); f = x(6:600:end); g = x(7:600:end); h = x(8:600:end); ii = x(9:600:end);j = x(10:600:end); k = x(11:600:end); l = x(12:600:end);

	xx = [];
	for i =1:26;
		xx = [xx;a(i);b(i);c(i);d(i);e(i);f(i);g(i);h(i);ii(i);j(i);k(i);l(i)];
	end

	tmp = permute(tmp,[2 1]);

    %Convert from Celsius to Kelvin!
    md.smb.temperatures_presentday = md.smb.temperatures_presentday+274.15;
    tmp = tmp + 274.15;

    %Convert from mm/month WE to m/yr WE!
    md.smb.precipitations_presentday = md.smb.precipitations_presentday * (1/1000) * 12;
    precip = precip * (1/1000) * 12;

    md.smb.rlaps=0;

    save pddforce md
    %%

filename = 'PDDforcing_precip.xlsx';  
writematrix(test, filename);
