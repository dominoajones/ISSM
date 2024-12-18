function ncexport = export2nc(dataset,elements,x,y,filename);
%export2nc - Saves data in model as netcdf file, compatible with QGIS
%
%   Usage:
%      ncexport = export2nc(data,filename)
%
%   dataset:    values of data to be saved
%   elements:   md.mesh.elements
%   x:          md.mesh.x
%   y:          md.mesh.y
%   filename:   name of resulting file
%   
%   You can change grid resolution or CRS as needed.

% Define your mesh and grid parameters
index = elements; % Delaunay triangulation defining the mesh
x = x; % x-coordinates of the mesh vertices
y = y; % y-coordinates of the mesh vertices
data = dataset; % vertex values of data to be interpolated

% Determine the bounds of the mesh
x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

% Define the resolution of the grid (e.g., 1000 meters)
grid_resolution = 1000; % in meters

% Generate the grid points
xgrid = x_min:grid_resolution:x_max;
ygrid = y_min:grid_resolution:y_max;
default_value = NaN; % value of points located out of the mesh

% Interpolate the data onto the grid
grid = InterpFromMeshToGrid(index, x, y, data, xgrid, ygrid, default_value);

% Define the NetCDF file name
nc_filename = filename;

% Create dimensions
nccreate(nc_filename, 'x', 'Dimensions', {'x', length(xgrid)});
nccreate(nc_filename, 'y', 'Dimensions', {'y', length(ygrid)});

% Create the variable and write the data
nccreate(nc_filename, 'friction_coefficient', 'Dimensions', {'x', length(xgrid), 'y', length(ygrid)});
ncwrite(nc_filename, 'x', xgrid);
ncwriteatt(nc_filename, 'x', 'units', 'meters');
ncwriteatt(nc_filename, 'x', 'standard_name', 'projection_x_coordinate');
ncwriteatt(nc_filename, 'x', 'long_name', 'x coordinate of projection');
ncwrite(nc_filename, 'friction_coefficient', grid');

ncwrite(nc_filename, 'y', ygrid);
ncwriteatt(nc_filename, 'y', 'units', 'meters');
ncwriteatt(nc_filename, 'y', 'standard_name', 'projection_y_coordinate');
ncwriteatt(nc_filename, 'y', 'long_name', 'y coordinate of projection');

% Add global attributes for CRS
ncwriteatt(nc_filename, '/', 'Conventions', 'CF-1.6');
ncwriteatt(nc_filename, '/', 'title', 'Friction Coefficient Data');
ncwriteatt(nc_filename, '/', 'institution', 'Your Institution');
ncwriteatt(nc_filename, '/', 'source', 'Model Data');
ncwriteatt(nc_filename, '/', 'references', 'Your References');

% Add variable attributes for friction_coefficient
ncwriteatt(nc_filename, 'friction_coefficient', 'long_name', 'Friction Coefficient');
ncwriteatt(nc_filename, 'friction_coefficient', 'units', 'unitless');
ncwriteatt(nc_filename, 'friction_coefficient', 'coordinates', 'x y');
ncwriteatt(nc_filename, 'friction_coefficient', 'grid_mapping', 'crs');

nccreate(nc_filename, 'crs', 'Datatype', 'char', 'Dimensions', {'string', 1});
ncwriteatt(nc_filename, 'crs', 'grid_mapping_name', 'polar_stereographic');
ncwriteatt(nc_filename, 'crs', 'longitude_of_projection_origin', 0.0);
ncwriteatt(nc_filename, 'crs', 'latitude_of_projection_origin', 90.0);
ncwriteatt(nc_filename, 'crs', 'standard_parallel', 70.0);
ncwriteatt(nc_filename, 'crs', 'straight_vertical_longitude_from_pole', -45.0);
ncwriteatt(nc_filename, 'crs', 'false_easting', 0.0);
ncwriteatt(nc_filename, 'crs', 'false_northing', 0.0);
ncwriteatt(nc_filename, 'crs', 'EPSG_code', 'EPSG:3413');

disp(['NetCDF file ' nc_filename ' created successfully.']);