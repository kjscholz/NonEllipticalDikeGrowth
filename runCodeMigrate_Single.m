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
% variables
Pe  = 2e6; % constant excess pressure (Pa)
param.zi = -1e3; % z-coordinate of injection point (m)
param.rho_m        = 2600;     % density fluid (kg/m3)
param.gamma_magma  = param.rho_m*param.g; % magmastatic gradient (Pa/m)

% time steps
begin_time     = 0;       % initialize time
end_time       = (6e4*3*param.eta^2)/((Pe/a_0)+(param.gamma_litho-param.gamma_magma))^2;%3e3/(Pe*((Pe/param.mu)^2*a_0)/(param.eta));      % maximum simulation time in seconds
time_vector    = linspace(begin_time,end_time,4000); % time steps at which to compute (s)


[param_out, store_yo, store_zo] =... %store_dh_dt,store_del,store_p_g_total,store_p_g_ex,store_p_g_grav,store_p_g_volcano,store_lo
    mainDikePropagate(F, a_0,Pe,time_vector,n,param); % run mainChamber.m
compute_time = toc
param_out.stop_reason
param_out.time_last
%% plots
figure
hold on
for i = 1:floor(param_out.time_last/20):param_out.time_last%

plot(store_yo(i,:),store_zo(i,:),'o')

end
plot(store_yo(param_out.time_last,:),store_zo(param_out.time_last,:),'o')

surface_line_x  = -RC:RC/100:RC;
surface_line_y  = -topo_profile(surface_line_x);

hold on
plot(surface_line_x,surface_line_y,'k','LineWidth',2)

xline(0,'k--')
yline(0,'LineWidth',2)

set(gca, 'Ydir','reverse')

axis equal ;
%%
% %% plots
% figure
% hold on
% for i = 1:2:param_out.time_last
% 
% plot(store_yo(i,200:300),store_zo(i,200:300),'o')
% 
% end
% plot(store_yo(param_out.time_last,200:300),store_zo(param_out.time_last,200:300),'o')
% 
% surface_line_x  = -RC:RC/100:RC;
% surface_line_y  = -topo_profile(surface_line_x);
% 
% hold on
% plot(surface_line_x,surface_line_y,'k','LineWidth',2)
% 
% xline(0,'k--')
% yline(0,'LineWidth',2)
% 
% set(gca, 'Ydir','reverse')
% 
% axis equal ;
% %% plots
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% plot(store_yo(i,:),store_zo(i,:),'o')
% 
% end
% plot(store_yo(param_out.time_last,:),store_zo(param_out.time_last,:),'o')
% 
% surface_line_x  = -RC:RC/100:RC;
% surface_line_y  = -topo_profile(surface_line_x);
% 
% hold on
% plot(surface_line_x,surface_line_y,'k','LineWidth',2)
% 
% xline(0,'k--')
% yline(0,'LineWidth',2)
% 
% set(gca, 'Ydir','reverse')
% 
% axis equal ;
% %%
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% scatter(i*ones(size(store_zo(i,200:300))),store_lo(i,200:300),40,log10(abs(store_p_g_total(i,200:300))),'filled')
% 
% end
% colorbar
% title("absolute value of total pressure  driving gradient log scale for observations points at 70^o -100^o ")
% ylabel("distance from source")
% xlabel("iteration")
% %%
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% scatter(i*ones(size(store_zo(i,200:7:300))),store_lo(i,200:7:300),40,log10(abs(store_p_g_total(i,200:7:300))),'filled','MarkerEdgeColor','k')
% 
% end
% %clim([-2e3,2e3])
% ylim([40,60])
% xlim([200,param_out.time_last+1])
% colorbar
% title("total pressure  driving gradient for observations points at 70^o -100^o ")
% ylabel("distance from source")
% xlabel("iteration")
% %%
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% scatter(i*ones(size(store_zo(i,200:7:300))),store_lo(i,200:7:300),40,(store_p_g_ex(i,200:7:300)),'filled','MarkerEdgeColor','k')
% 
% end
% clim([1.5e4,2.5e4])
% ylim([40,60])
% xlim([200,param_out.time_last+1])
% colorbar
% title("excess pressure  driving gradient for observations points at 70^o -100^o ")
% ylabel("distance from source")
% xlabel("iteration")
% %%
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% scatter(i*ones(size(store_zo(i,200:7:300))),store_lo(i,200:7:300),40,(store_p_g_volcano(i,200:7:300)),'filled','MarkerEdgeColor','k')
% 
% end
% %clim([1.5e4,2.5e4])
% ylim([40,60])
% xlim([200,param_out.time_last+1])
% colorbar
% title("volcano stress model pressure  driving gradient for observations points at 70^o -100^o ")
% ylabel("distance from source")
% xlabel("iteration")
% 
% %%
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% scatter(i*ones(size(store_zo(i,200:7:300))),store_lo(i,200:7:300),40,(store_p_g_grav(i,200:7:300)),'filled','MarkerEdgeColor','k')
% 
% end
% %clim([1.5e4,2.5e4])
% ylim([40,60])
% xlim([200,param_out.time_last+1])
% colorbar
% title("fluid gravity/buoyancy pressure  driving gradient for observations points at 70^o -100^o ")
% ylabel("distance from source")
% xlabel("iteration")
% %%
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% scatter(i*ones(size(store_zo(i,200:7:300))),store_lo(i,200:7:300),40,log10(abs(store_dh_dt(i,200:7:300))),'filled','MarkerEdgeColor','k')
% 
% end
% %clim([-2e-2,2e-2])
% ylim([40,60])
% xlim([200,param_out.time_last+1])
% colorbar
% title("log10(abs(dh/dt)) for observations points at 70^o -100^o ")
% ylabel("distance from source")
% xlabel("iteration")
% %%
% figure
% hold on
% 
% 
% scatter(1:length(store_del),store_del,'filled','MarkerEdgeColor','k')
% 
% %clim([-2e3,2e3])
% %ylim([40,60])
% xlim([200,param_out.time_last+1])
% %colorbar
% title("Dike opening (del)")
% ylabel("Dike opening (del)")
% xlabel("iteration")
% %%
% figure
% hold on
% for i = 1:1:param_out.time_last
% 
% scatter(i*ones(size(store_zo(i,220:270))),store_lo(i,220:270),'filled')
% 
% end
% % %%
% % for i = param_out.time_last-5:1:param_out.time_last
% % figure
% % hold on
% % 
% % 
% % scatter(store_yo(i,:),store_zo(i,:),50,store_p_g_total(i,:),'filled')
% % 
% % 
% % colorbar
% % 
% % set(gca, 'Ydir','reverse')
% % 
% % axis equal ;
% % end
% % %%
% % for i = param_out.time_last-5:1:param_out.time_last
% % figure
% % hold on
% % 
% % 
% % scatter(store_yo(i,:),store_zo(i,:),50,store_p_g_volcano(i,:),'filled')
% % 
% % 
% % colorbar
% % 
% % set(gca, 'Ydir','reverse')
% % 
% % axis equal ;
% % end
% % 
% % %%
% % for i = param_out.time_last-5:1:param_out.time_last
% % figure
% % hold on
% % 
% % 
% % scatter(store_yo(i,:),store_zo(i,:),50,store_p_g_grav(i,:),'filled')
% % 
% % 
% % colorbar
% % 
% % set(gca, 'Ydir','reverse')
% % 
% % axis equal ;
% % end
% % %%
% % for i = param_out.time_last-10:1:param_out.time_last
% % figure
% % hold on
% % 
% % 
% % scatter(store_yo(i,:),store_zo(i,:),50,store_p_g_ex(i,:),'filled')
% % 
% % 
% % colorbar
% % 
% % set(gca, 'Ydir','reverse')
% % 
% % axis equal ;
% % end
% % %%
% % for i = param_out.time_last-4:1:param_out.time_last
% % figure
% % hold on
% % 
% % 
% % scatter(store_yo(i,:),store_zo(i,:),50,store_p_g_ex(i,:),'filled')
% % 
% % 
% % colorbar
% % 
% % set(gca, 'Ydir','reverse')
% % 
% % axis equal ;
% % end
% % %% plots
% % figure
% % hold on
% % for i = param_out.time_last-5:1:param_out.time_last
% % 
% % plot(store_yo(i,:),store_zo(i,:),'o')
% % 
% % end
% % scatter(store_yo(param_out.time_last-1,:),store_zo(param_out.time_last-1,:),...
% %         500,1:1:1001,'filled','s')
% % colormap(parula(8*2))
% % colorbar
% % 
% % surface_line_x  = -RC:RC/100:RC;
% % surface_line_y  = -topo_profile(surface_line_x);
% % 
% % hold on
% % plot(surface_line_x,surface_line_y,'k','LineWidth',2)
% % 
% % xline(0,'k--')
% % yline(0,'LineWidth',2)
% % 
% % set(gca, 'Ydir','reverse')
% % 
% % axis equal ;
% % xlim([-3000,3000])
% % ylim([-3000,3000])
% 
% % 
% % %% plots
% % figure
% % hold on
% % plot(1:400, store_dh_dt(1:400,255));%))
% % 
% % %% plots
% % figure
% % hold on
% % plot(1:4000, store_del(1:4000));%))