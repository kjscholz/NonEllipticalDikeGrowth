clear all
%fname = 'stress_models/step01a_materials_cauchy_stress.csv';
fname = 'StressModels/large_cone_gelatin_highRes.mat';
pylith_model = 0; % T/F for pylith model which has different stress convention
x0 = 0; % coordinate position of X cross section
g = 9.807;


%%% Crustal parameters
param.rho_c = 2700; % Crust density (kg/m^3)
param.sm    = 1e3;  % Shear modulus (MPa)
param.pr    = 0.25; % Poisson's ratio

if pylith_model
    HC = 2600; % Elevation of edifice

    func = @(R)topo_profile(R)-100; % creating function handle for when topo is close enough to zero
    RC = round(fzero(func,13e3),-3)

    data = readmatrix(fname);

    % Assign Apprsopriate Field Names
    X = data(:,1);
    Y = data(:,2);
    % CHANGE SIGN CONVENTION FROM Down neg to Down Positive
    Z = -1*data(:,3);

    % CHANGE SIGN CONVENTION FROM Compression Negative to Compression Positive
    Sxx = -1*data(:,4); Syy = -1*data(:,5); Szz = -1*data(:,6);
    Sxy = -1*data(:,7); Syz = -1*data(:,8); Sxz = -1*data(:,9);

    [SxxInterp, ~, ~,~, ~,~] = StressInterp(X,Y,Z,Sxx,Syy,Szz, Sxy, Sxz, Syz);
else
    % load & rescale stress model
    load(fname)

    % get rid of stress zeros at the edifice base
    SXX(2,:)=SXX(3,:).*1.0156;
    SXX(1,:)=SXX(2,:);
    tic
    [SxxInterp, ~, ~,~, ~,~] = StressInterp(Y',Z',SXX');
    toc

end

fname_save= "StressModels/CrustalState_LargeGelatinCone.mat";

save(fname_save,"SxxInterp",...
    "param", "HC","RC","x0","g","fname",'-mat')
%% Functions
function [SxxInterp, SyyInterp, SzzInterp,SxyInterp, SxzInterp, SyzInterp] = StressInterp(varargin)
% inputs:
%   (X, Y, Z, Sxx, Syy, Szz, Sxy, Sxz, Syz)
%   or
%   (R, Z, Sxx, Syy, Szz, Sxy, Sxz, Syz) if axis symmetric
% outputs:
%   [SxxInterp, SyyInterp, SzzInterp,SxyInterp, SxzInterp, SyzInterp)
extrap_method = "nearest";
if isequal(nargin,9)
    % varargin is Stensor (X, Y, Z, Sxx, Syy, Szz, Sxy, Sxz, Syz)
    X = varargin{1}; Y = varargin{2}; Z = varargin{3};
    Sxx = varargin{4}; Syy = varargin{5}; Szz = varargin{6};
    Sxy = varargin{7}; Sxz = varargin{8}; Syz = varargin{9};

    SxxInterp = scatteredInterpolant(X(:), Y(:), Z(:), Sxx(:), "linear", extrap_method);
    SyyInterp = scatteredInterpolant(X(:), Y(:), Z(:), Syy(:), "linear", extrap_method);
    SzzInterp = scatteredInterpolant(X(:), Y(:), Z(:), Szz(:), "linear", extrap_method);
    SxyInterp = scatteredInterpolant(X(:), Y(:), Z(:), Sxy(:), "linear", extrap_method);
    SxzInterp = scatteredInterpolant(X(:), Y(:), Z(:), Sxz(:), "linear", extrap_method);
    SyzInterp = scatteredInterpolant(X(:), Y(:), Z(:), Syz(:), "linear", extrap_method);

elseif isequal(nargin,8)
    % varargin is Stensor (R, Z, Sxx, Syy, Szz, Sxy, Sxz, Syz)
    R = varargin{1}; Z = varargin{2};
    Sxx = varargin{3}; Syy = varargin{4}; Szz = varargin{5};
    Sxy = varargin{6}; Sxz = varargin{7}; Syz = varargin{8};

    SxxInterp = scatteredInterpolant(R(:), Z(:), Sxx(:), "linear", extrap_method);
    SyyInterp = scatteredInterpolant(R(:), Z(:), Syy(:), "linear", extrap_method);
    SzzInterp = scatteredInterpolant(R(:), Z(:), Szz(:), "linear", extrap_method);
    SxyInterp = scatteredInterpolant(R(:), Z(:), Sxy(:), "linear", extrap_method);
    SxzInterp = scatteredInterpolant(R(:), Z(:), Sxz(:), "linear", extrap_method);
    SyzInterp = scatteredInterpolant(R(:), Z(:), Syz(:), "linear", extrap_method);
elseif isequal(nargin,3)
    % varargin is Stensor (R, Z, Sxx)
    R = varargin{1}; Z = varargin{2};
    Sxx = varargin{3};
    SxxInterp = griddedInterpolant(R, Z, Sxx, "linear", extrap_method);
    SyyInterp = [];
    SzzInterp = [];
    SxyInterp = [];
    SxzInterp = [];
    SyzInterp = [];

end

end