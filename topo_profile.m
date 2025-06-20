function [Z] = topo_profile(varargin)
% topo_profile(R): 
%   assumes R is a radial position or array of radial positions
% topo_profile(X, Y):
%   assumes X = x_coordinates and Y = y_coordinates
%   X and Y must have the same shape
% RETURNS    z which has size(R) and contains elevation above free surface
%
% Hard coded edifice shape constants and heigh using values from
% edifice_shape.csv

if nargin == 2
    X = varargin{1};
    Y = varargin{2};
    if ~isequal(size(X),size(Y))
        error("X and Y must have same dimensions")
    end
    R = sqrt(X.^2+Y.^2);
elseif nargin == 1
    R = varargin{1};
else
    error("Incorrect Number of args passed in: pass 1 argument is radial position, 2 arguments is x and y position")
end
% constants for fitting Mt Fuji
a1 = 0.3609019150579505; b1 = 2.524518954184922; 
a2 = 0.620889711991542; b2 = 0.2647231053801891;

H = 2600; % height of edifice in m 

R_nd = R/H; % nondimensional heigh
Z_nd = a1.*exp(-b1.*R_nd.^2)+a2.*exp(-b2.*R_nd.^2);
Z = Z_nd*H;
end