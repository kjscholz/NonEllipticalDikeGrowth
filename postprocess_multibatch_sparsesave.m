clear all
close all
batchnum = 16;
loadBackgroundStress = false; % T/F For loadign stress in
name = sprintf('batch_out/runs_%03d.mat',batchnum);
load(name)

grid_elements = numel(density_grid);
%%

if loadBackgroundStress
    y = linspace(-RC,RC,1000);
    z = linspace(max(depth_grid(:))*1.1,-HC,1000);
    [Y,Z] = meshgrid(y,z);
    X = output_param{1}.x0*ones(size(Y));
    Z_topo = topo_profile(Y);
    topo_mask = find(-Z>Z_topo);
    stress_fpath = "StressModels/CrustalStressGridded_SummerCoon2600gg_lithostatic_a.mat";
    load(stress_fpath);
    SXX = SxxInterp(X,Y,Z);
    SXX(topo_mask)=nan;
end
%%
figure
%scatter(overpressure_grid(:),density_grid(:),500,-vent_height(:),'filled','s')

contourf(squeeze(overpressure_grid)/1e6,squeeze(density_grid),...
    -squeeze(vent_height))
colormap(parula(15))
cbar=colorbar;
xlim([0,20])
ylim([2050,2650])
clim([0,1500])
%%
figure
%scatter(overpressure_grid(:),density_grid(:),500,-vent_height(:),'filled','s')

contour(squeeze(overpressure_grid)/1e6,squeeze(density_grid),...
    -squeeze(vent_height), 'LineColor', 'k')


xlim([0,20])
ylim([2050,2650])
%clim([0,1500])
%%
set(0, 'DefaultAxesFontName', 'Arial')
set(0, 'DefaultTextFontName', 'Arial')
figure
scatter(overpressure_grid(:)/1e6,density_grid(:),200,-vent_height(:),'filled','o')
colormap(parula(15))

cbar=colorbar;
xlim([0.75,20.25])
ylim([2075,2625])
clim([0,1500])
hold on
[C, h]=contour(squeeze(overpressure_grid)/1e6,squeeze(density_grid),...
    -squeeze(vent_height),0:250:1500, 'LineColor', 'k','LineWidth',2);
clabel(C, h);

pbaspect([1 1 1])
ylabel("Magma Density (kg/m^3)")
xlabel("Excess Magma Pressure (MPa)")
ylabel(cbar,"Vent Elevtion (m)")
title(['Initial Depth ' num2str(depth_grid(1)) ' km']);
fontsize(gcf,scale=1.5)
ax=gca;
ax.YTick = 2100:100:2600;
ax.XTick = 0:5:20;

ax.XMinorTick = 'on';
%% plots
surface_line_x  = -12e3:12e3/100:12e3;
surface_line_y  = -topo_profile(surface_line_x);
domain_limy = [-8e3,8e3];
domain_limz = [-2600, 12e3];
fig_size=800
figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size])
count =1;
tiledlayout(2,3,"TileSpacing","compact")
ax={}
density_title={};
for jj=1:1:grid_elements
    fignum=mod(jj,6)+1;
    time_vector=output_time{jj};
    store_yo = output_yo{jj};
    store_zo = output_zo{jj};
    param = output_param{jj};
    
    if jj<7
        ax{fignum}=nexttile(fignum);
        density_title{fignum}=density_grid(jj);
    end
            plot(ax{fignum},store_yo,-store_zo, LineWidth=1)



            count=count+1;
            hold on

axis equal ;
            %scatter(vent_radius(jj),-vent_height(jj),40,'magenta','filled')
            if strcmp(param.stop_reason,"time")
                disp("time")
            end

   end
for i=1:6

    plot(ax{i},surface_line_x,-surface_line_y,'k','LineWidth',2)

    title(ax{i},['\rho_m = ' num2str(density_title{i}) ' kg/m^3']);
    xline(ax{i},0,'k--')
    %yline(0,'LineWidth',2)
    xlabel(ax{i},'Radial Position (m)')
    ylabel(ax{i},'Vertical Position (m)')
    %set(gca, 'Ydir','reverse')

    

    xlim(ax{i},[-12e3, 12e3])
    ylim(ax{i},[-10e3, 2600])
    
end
%fontsize(gcf,scale=2.5)

    ylim(ax{1},[-10e3, 2600])