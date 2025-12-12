clear all
close all
batchnum = 12;
loadBackgroundStress = false; % T/F For loadign stress in
name = sprintf('batch_out/runs_%03d.mat',batchnum);
load(name)

grid_elements = numel(density_grid);


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
scatter3(overpressure_grid(:),density_grid(:),depth_grid(:),...
    500,-vent_height(:),'filled','s')
colormap(parula(8*2))
colorbar
xlim([0,11e6])
ylim([2150,2650])
clim([0,1500])
%% plots
surface_line_x  = -12e3:12e3/100:12e3;
surface_line_y  = -topo_profile(surface_line_x);
domain_limy = [-8e3,8e3];
domain_limz = [-2600, 12e3];
fig_size=800
figure('Renderer', 'painters', 'Position', [10 10 fig_size*1.5 fig_size])
count =1;
let=['a','b','c','d'];
tiledlayout(2,2,"TileSpacing","compact")
for jj=1:1:grid_elements
    time_vector=output_time{jj};
    store_yo = output_yo{jj};
    store_zo = output_zo{jj};
    param = output_param{jj};

    ylim_min = max(min(store_yo(:)), domain_limy(1));
    ylim_max = max(max(store_yo(:)), domain_limy(2));

    zlim_min = domain_limz(1);
    zlim_max = min(max(store_zo(:)), domain_limz(2));
    % if abs(depth_grid(jj)-6e3)<0.1
    %     if (density_grid(jj)==2200)||(density_grid(jj)==2600)
    %         if (overpressure_grid(jj)==1e6)||(overpressure_grid(jj)==7e6)
                nexttile
                if loadBackgroundStress
                    contourf(Y,Z,SXX)
                end

                hold on
                step = 100;%floor(param.time_last/34);
                if count ==1
                    step = 10;
                end
                base = floor(param.time_last/step);
                start = rem(param.time_last,step);
                param.time_last;
                if start ==0
                    start =1;
                end


                for i = start:step:param.time_last+1
                    plot(store_yo(i,:),-store_zo(i,:), LineWidth=2)

                end

                if count==4
                    plot(store_yo(param.time_last,:),-store_zo(param.time_last,:),LineWidth=2)
                end
                
                count=count+1;
                hold on
                plot(surface_line_x,-surface_line_y,'k','LineWidth',2)
                plot([0e4, 5e4],[-4e3, 5e4*tand(45)],'g')
                title([let(jj) ') \rho_m = ' num2str(density_grid(jj)) 'kg/m^3 P_e = ' num2str(overpressure_grid(jj)/1e6) ' MPa'])
                xline(0,'k--')
                %yline(0,'LineWidth',2)
                xlabel('Radial Position (m)')
                ylabel('Vertical Position (m)')
                %set(gca, 'Ydir','reverse')

                axis equal ;
                xlim([-12e3, 12e3])
                ylim([-9e3, 2600])

                %scatter(vent_radius(jj),-vent_height(jj),40,'magenta','filled')


            end
%         end
%     end
% end
fontsize(gcf,scale=2.5)