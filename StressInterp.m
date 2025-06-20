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

end

end