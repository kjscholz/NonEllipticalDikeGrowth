% Computes the velocity of observations points along the tipline of a dike.
% Fluid is assumed to travel strictly radially away from injection point,
% which is fixed at some depth and radial location with respect to the base
% of the volcanic edifice. Constant excess pressure at dike inlet.
% Written by M. Townsend July 2024

clearvars -except F param HC RC x0
keep_stress_model = 0;
sparse_save = 0;
% close all
if keep_stress_model==0
    tic
    % load & rescale stress model
    stress_fpath = 'StressModels/CrustalStressGridded_SummerCoon2600gg_lithostatic_a.mat';
    load(stress_fpath)
    toc
    param.x0=0;
    HC = 2600; % Elevation of edifice

    func = @(R)topo_profile(R)-100; % creating function handle for when topo is close enough to zero
    RC = round(fzero(func,13e3),-3)
    % get rid of stress zeros at the edifice base
    % SXX(2,:)=SXX(3,:).*1.0156;
    % SXX(1,:)=SXX(2,:);
    % SXX=SXX.*1; % Scale stresses??
    F = SxxInterp;
end
% Set up parameter grid
depth_list = 0e3; % m
density_list = 2200:200:2600;
overpressure_list = 1e6:3e6:10e6;

[depth_grid, density_grid, overpressure_grid] = meshgrid(depth_list,density_list,overpressure_list);
grid_elements = numel(depth_grid);
batchnum=26;
% Preallocate Outputs Based on Grid Size
output_time = cell(1,grid_elements);
output_yo = cell(1,grid_elements);
output_zo = cell(1,grid_elements);
output_param = cell(1,grid_elements);

vent_height = zeros(size(depth_grid));
vent_radius = zeros(size(depth_grid));

% Assign injection point
% **Injection point is where fluid comes into the dike. Excess pressure here
% is assumed to remain constant in this version of the model**
param.yi = 0e3; % y-coordinate of injection point (m)

% create structure for constants
%param.youngs       = 1700;      % Young's modulus (Pa)
%param.pr           = 0.5;       % Poisson's ratio crust
param.mu           = param.sm;       % Shear modulus crust (Pa)
param.rho_c        = 2700;      % density crust (kg/m3)
param.g            = 9.81;      % gravitational acceleration (m/s2)
param.gamma_litho  = param.rho_c*param.g; % lithostatic gradient (Pa/m)

% constants related to dike propagation
param.eta          = 1e5;     % fluid viscosity (Pa s)
param.f_d          = 1;       % aperture of dike tip as fraction of max opening (Rubin, 1995)
param.f_p          = 1;       % fraction of pressure gradient at dike tip (Rubin, 1995)
param.f_v          = 1;       % fraction of average flow velocity at entrance (Rubin, 1995)

n   = 1001; % number of observation points
a_0 = 50; % initial dike radius (m)

% boundary conditions -- inflow conditions
%q_in_mL_min = 1.05; % mL/min
%param.q_in  = q_in_mL_min*1e-6/60; % convert to m^3/s
tic
for ii=1:grid_elements
    % variables
    Pe  = overpressure_grid(ii); % constant excess pressure (Pa)
    param.zi = depth_grid(ii); % z-coordinate of injection point (m)
    param.rho_m        = density_grid(ii);     % density fluid (kg/m3)
    param.gamma_magma  = param.rho_m*param.g; % magmastatic gradient (Pa/m)

    % time steps
    begin_time     = 0;       % initialize time
    end_time       = (6e4*3*param.eta^2)/((Pe/a_0)+(param.gamma_litho-param.gamma_magma))^2;      % maximum simulation time in seconds
    time_vector    = linspace(begin_time,end_time,8000); % time steps at which to compute (s)
    output_time{ii} = time_vector;
    [output_param{ii}, output_yo{ii}, output_zo{ii}] = mainDikePropagate(F, a_0,Pe,time_vector,n,param); % run mainChamber.m
    output_param{ii}.stop_reason
    if strcmp(output_param{ii}.stop_reason,"topo")
        z_topo_final = -topo_profile(output_yo{ii}(output_param{ii}.time_last,:))+1;
        inds = output_zo{ii}(output_param{ii}.time_last,:)<=z_topo_final;
        vent_h_opts = output_zo{ii}(output_param{ii}.time_last,inds);
        vent_r_opts = output_yo{ii}(output_param{ii}.time_last,inds);
        [vent_h,v_ind] =min(vent_h_opts);
        vent_r = vent_r_opts(v_ind);
        vent_height(ii)=vent_h;
        vent_radius(ii)=vent_r;
    else
        vent_height(ii)=NaN;
        vent_radius(ii)=NaN;
    end
    
    if sparse_save ==1
    output_yo{ii}=output_yo{ii}(output_param{ii}.time_last,:);
    output_zo{ii}=output_zo{ii}(output_param{ii}.time_last,:);
    end
end
compute_time = toc
%%
name = sprintf('batch_out/runs_%03d.mat',batchnum);
save(name,"density_grid","overpressure_grid","depth_grid","vent_height","vent_radius",...
    "output_time", "output_yo", "output_zo",...
    "output_param","compute_time", "HC","RC","stress_fpath","sparse_save",'-v7.3');
%% plots
for jj=1:1:grid_elements %floor(grid_elements/3)
    time_vector=output_time{jj};
    store_yo = output_yo{jj};
    store_zo = output_zo{jj};
    figure
    hold on
    for i = 1:floor(output_param{ii}.time_last/20):output_param{ii}.time_last %floor(length(time_vector)/20)

        plot(store_yo(i,:),store_zo(i,:),'o')
   

    end
    plot(store_yo(output_param{ii}.time_last,:),store_zo(output_param{ii}.time_last,:),'o')

    surface_line_x  = -RC:RC/100:RC;
    surface_line_y  = -topo_profile(surface_line_x);

    hold on
    plot(surface_line_x,surface_line_y,'k','LineWidth',2)

    xline(0,'k--')
    yline(0,'LineWidth',2)

    set(gca, 'Ydir','reverse')

    axis equal ;
end
