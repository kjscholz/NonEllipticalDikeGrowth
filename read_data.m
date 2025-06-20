
% Read in experimental data
large_cone_data = xlsread('large_cone_constant_flux.xlsx');
t_large   = large_cone_data(1:end,1); % time in minutes
b_large   = large_cone_data(1:end,2); % height of fluid-filled portion (m)
h_large   = large_cone_data(1:end,3); % position of upper dike tip from original injection point (m)
a2_large  = large_cone_data(1:end,4); % breadth of fluid-filled portion (m)


% Read in experimental data
medium_cone_data = xlsread('medium_cone_constant_flux.xlsx');
t_medium   = medium_cone_data(1:end,1); % time in minutes
b_medium   = medium_cone_data(1:end,2); % height of fluid-filled portion (m)
h_medium   = medium_cone_data(1:end,3); % position of upper dike tip from original injection point (m)
a2_medium = medium_cone_data(1:end,4); % breadth of fluid-filled portion (m)


% Read in experimental data
small_cone_data = xlsread('small_cone_constant_flux.xlsx');
t_small   = small_cone_data(1:end,1); % time in minutes
b_small   = small_cone_data(1:end,2); % height of fluid-filled portion (m)
h_small   = small_cone_data(1:end,3); % position of upper dike tip from original injection point (m)
a2_small= small_cone_data(1:end,4); % breadth of fluid-filled portion (m)




