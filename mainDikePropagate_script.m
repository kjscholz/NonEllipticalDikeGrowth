
% Prescribe observation points along dike tip-line
% **points are initially linearly spaced in their angle wrt injection point
b_0 = a_0;
param.thetas = linspace(pi,-pi,n); % angle of each point in radians

param.yo_0 = param.yi+a_0.*cos(param.thetas); % y-coordinates of observations points (m)
param.zo_0 = param.zi-a_0.*sin(param.thetas); % z-coordinates of observation points (m)

[~,param.ind_upper] = min(abs(param.thetas-pi/2));
[~,param.ind_lower] = min(abs(param.thetas+pi/2));
[~,param.ind_left]  = min(abs(param.thetas-pi));
[~,param.ind_right] = min(abs(param.thetas-0));

% Initial dike volume, assumed to start as penny shape (m^3)
major = a_0; minor = b_0;
k_0     = real(sqrt(1-(minor^2)/(major^2))); % modulus
[~,E_0] = ellipke(k_0^2);
del_0 = (2*a_0*(1-param.pr)*Pe)/(param.mu*E_0); % m
Vd_0       = (2*pi/3)*a_0*b_0*del_0; %initial dike volume (m^3)


% initial observations points (m)
yo = param.yo_0;
zo = param.zo_0;
x0_o = zeros(size(zo));

% initialize the array to store observations points
store_yo = NaN(numel(time_vector),numel(yo));
store_zo =  NaN(numel(time_vector),numel(yo));

for i = 1:length(time_vector)-1

    dt = time_vector(i+1)-time_vector(i); % time step in seconds
    lo = sqrt((yo-param.yi).^2+(zo-param.zi).^2); % distances between injection point and observations points (m)

    % % % COMPUTE FLUID PRESSURE GRADIENTS AT EACH POINT % % %

    % First term related to runout of excess pressure
    p_g_ex = param.f_p*Pe./lo; % (Pa/m)

    % Second term related to gravity on Fluid
    p_g_grav    = sin(param.thetas).*(param.gamma_litho-param.gamma_magma); % gravitational gradient

    % Third term from load of volcanic edifice + lithostatic
    Sv_at_zo_yo = F(x0_o,yo,zo); % value of volcano stress at observation points
    Sv_at_zi_yi = F(x0,param.yi,param.zi); % value of volcano stress at injection point
    dSv = Sv_at_zo_yo-Sv_at_zi_yi;
    p_g_volcano = -dSv./lo; 

    % new points at which to evaluate volcano stress  
    %yn = yo+param.dl.*cos(param.thetas); % y-coordinates of points for evaluating stress gradient (m)
    %zn = zo-param.dl.*sin(param.thetas); % z-coordinates of points for evaluating stress gradient (m)
    %Sv_at_zn_yn = F(zn,yn); % value of stress at points zn, yz
    %dSv = Sv_at_zn_yn-Sv_at_zo_yo;
    %p_g_volcano = -dSv./param.dl;

    p_g_total = p_g_ex+p_g_grav+p_g_volcano; % Total fluid pressure gradients (Pa/m)
    % disp(mean(p_g_total))

    % % % COMPUTE VELOCITIES OF OBSERVATION POINTS % % %

    % 1. Compute max dike opening (pretend it's an ellipsoid for this)

    % compute dike height near injection point

    % use min and max lengths from injection point as semi-minor and semi-major axis of
    % ellipse
    %minor = min(lo); % semi-minor axis of dike (m)
    %major = max(lo); % semi-major axis of dike (m)
    b = zo(param.ind_lower) - zo(param.ind_upper); % height of dike near injection point (m)
    a = yo(param.ind_right) - yo(param.ind_left);  % horizontal breadth of dike near injection point (m)

    % Elliptic integrals
    if b > a
        major = b/2; minor = a/2;
    else
        major = a/2; minor = b/2;
    end

    k     = real(sqrt(1-(minor^2)/(major^2))); % modulus
    [F2,E2] = ellipke(k^2);
    A2 = ((E2-F2)/E2)*(1/(((major^2)/(minor^2))-1));
    % max dike opening (ish)
    del = (2*(b/2)*(1-param.pr)*Pe)/(param.mu*E2); % m

    % 2. Compute velocities of observation points (m/s)
    dh_dt = param.f_v*(1/(3*param.eta))*p_g_total*((param.f_d*del)^2);
    db_dt = (dh_dt(param.ind_upper)+dh_dt(param.ind_lower))/2;
    da_dt = (dh_dt(param.ind_left)+dh_dt(param.ind_right))/2;

    % % % ADVANCE POINTS BASED ON TIME STEP % % %
    dh = dh_dt*dt; % increment of growth at each point (m)
    lo = lo+dh; % new distance of observation points from injection site (m)
    yo = param.yi+lo.*cos(param.thetas); % y-coordinates of observations points (m)
    zo = param.zi-lo.*sin(param.thetas); % z-coordinates of observation points (m)

    store_yo(i,:) = yo;
    store_zo(i,:) = zo;

end



