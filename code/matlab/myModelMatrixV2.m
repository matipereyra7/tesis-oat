% Creating Model-based Matrix (MM)

%% GRID
Ns = 64; % number of sensors
Nt = 512; % temporal samples
Dx = 125e-6; % pixel size (Dy = Dz) [m]
nx = 64; % number of pixels in the 2-D image region
N = nx*nx; % total number of pixels in the image region
vs = 1485; % speed of sound [m/s]
cr = 12.512e-3; % radius of the circumference where the sensors are arranged
las = 45e-3; % length of the integrating line detector
rho = 1000; % soft tissue or water density [kg/m^3]
to = 0;     % initial time [s]
tf = 15e-6; % final time [s]
sensor_arc = 348; % [�]

% Create the computational grid
Nx = round(cr*2/Dx); % number of grid points in the x direction (row)
Ny = round(cr*2/Dx); % number of grid points in the y direction (column)
Nz = round(ls/Dx);   % number of grid points in the z direction
% Origin in the center of the volume
if mod(Nx,2)==0
    Nx = Nx + 1;
end
if mod(Ny,2)==0
    Ny = Ny + 1;
end
if mod(Nz,2)==0
    Nz = Nz + 1;
end
dx = Dx;             % grid point spacing in the x direction [m]
dy = Dx;             % grid point spacing in the y direction [m]
dz = Dx;             % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);          % for 2-D simulations
%kgrid = KWaveGrid(Nx, dx, Ny, dy, Nz, dz); % for 3-D simulations
gridorigin2D = round(Nx/2); 

%% MEDIUM
% Define the properties of the propagation medium
medium.sound_speed = vs * ones(Nx, Ny);      % for 2-D simulations [m/s]
medium.density = rho * ones(Nx, Ny);         % for 2-D simulations [kg/m^3]

%% USING RADIAL SENSOR DISTRIBUTION
% Define a Cartesian sensor mask of a centered circle with Ns sensor elements
sensor_radius = cr;         % [m]
num_sensor_points = Ns;
th = linspace(0, sensor_arc*pi/180, Ns + 1); th = th(1:Ns);
sensor.mask = sensor_radius*[cos(th);sin(th)];
%% USING LINEAR ARRAY
% The mask is 64*

%%
% Set the time stepping
dt = Dx/vs;             % [s]
%kgrid.t_array = (0:(Nt-1))*dt;
kgrid.t_array = linspace(to,tf,Nt);

% Initialise system matrix
A = zeros(Ns*Nt,N);

cont = 0; % MM column number
the_minx = 100;
the_miny = 100;
the_maxx = 0;
the_maxy = 0;
for kky=1:nx
    for kkx=1:nx
        cont = cont + 1;
        % Create initial pressure distribution
        p0 = zeros(Nx,Ny);
        corx=gridorigin2D-round(nx/2)+kkx;
        cory=gridorigin2D-round(nx/2)+kky;
        p0(corx,cory)=1; % [Pa]
        
        the_minx = min(the_minx, corx);
        the_miny = min(the_miny, cory);
        the_maxx = max(the_maxx, corx);
        the_maxy = max(the_maxy, cory);
        continue
        % Set the initial condition (PA source term)
        source.p0 = p0;
        
        % Run the k=Wave simulation
        input_args = {'PMLInside',false, 'PMLSize', 8, 'Smooth', false, 'PlotPML', false, 'PlotSim', false};
        % PLMInside: if false perfectly matched layer is outside the grid (the grid are enlarged)
        % PMLSize: size of the perfectly matched layer in grid points
        % Smooth: controlling whether source.p0, medium.sound_speed, medium.density are smoothed before computation
        % PlotPML: show perfectly matched layer in the simulation plots
        % PlotSim: plot progressively the simulation iteration
        
        % Run the simulation
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
        
        % Record column cont
        % if the matrix is going to be used in a Matlab environment do this:
        %A(:,cont) = sensor_data(:);
        % if the matrix is going to be used in a Python environment do this:
        A(:,cont) = reshape(transpose(reshape(sensor_data,[Ns,Nt])),[Ns*Nt,1]);
        
        display(['Iteration: ' num2str(cont) ' of ' num2str(N)])
    end
end

disp(the_minx)
disp(the_maxx)
disp(the_miny)
disp(the_maxy)
% Threshold the matrix to remove small entries and make it more sparse
thresh = 0; %6
if thresh > 0
    A(abs(A)<10^(-thresh)) = 0;
end

% save for python
% A=A'; %in python, when loaded using h5py, the matrix is transposed. In case of the module mat73, it is not necessary this line.
save('../data/forwMat1.mat','-v7.3','A')
save('forwMat_delete_me.mat','-v7.3','A')

