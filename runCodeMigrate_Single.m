% Computes the velocity of observations points along the tipline of a dike.
% Fluid is assumed to travel strictly radially away from injection point,
% which is fixed at some depth and radial location with respect to the base
% of the volcanic edifice. Constant excess pressure at dike inlet.
% Written by M. Townsend July 2024

clearvars -except F param HC RC x0
keep_stress_model = 0;
% close all
if keep_stress_model==0
    tic
    % load & rescale stress model
    load('StressModels/CrustalStressGridded_SummerCoon2600gg_lithostatic_c.mat')
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

% Assign injection point
% **Injection point is where fluid comes into the dike. Excess pressure here
% is assumed to remain constant in this version of the model**
param.yi = 1e3; % y-coordinate of injection point (m)

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
% variables
Pe  = 2e6; % constant excess pressure (Pa)
param.zi = 6e3; % z-coordinate of injection point (m)
param.rho_m        = 2400;     % density fluid (kg/m3)
param.gamma_magma  = param.rho_m*param.g; % magmastatic gradient (Pa/m)

% time steps
begin_time     = 0;       % initialize time
end_time       = 3e3/(Pe*((Pe/param.mu)^2*a_0)/(param.eta));      % maximum simulation time in seconds
time_vector    = linspace(begin_time,end_time,4000); % time steps at which to compute (s)


[param_out, store_yo, store_zo] = mainDikePropagate(F, a_0,Pe,time_vector,n,param); % run mainChamber.m
compute_time = toc
param_out.stop_reason
param_out.time_last
%% plots
figure
hold on
for i = 1:floor(param_out.time_last/20):param_out.time_last

plot(store_yo(i,:),store_zo(i,:),'o')

end


surface_line_x  = -RC:RC/100:RC;
surface_line_y  = -topo_profile(surface_line_x);

hold on
plot(surface_line_x,surface_line_y,'k','LineWidth',2)

xline(0,'k--')
yline(0,'LineWidth',2)

set(gca, 'Ydir','reverse')

axis equal ;

