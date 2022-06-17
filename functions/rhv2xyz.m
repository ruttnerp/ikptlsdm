function xyz = rhv2xyz(rhv, Type)
% Function for calculating polar coordinates into cartesian coordinates.
%
% input:  xyz...  [nx3] cartesian coordinates in [m]
%         Type... [string] defines orientation of system
% output: rhv...  [nx3] polar coordinates (range, horizontal, vertical) in [m,rad,rad]
%
% author:   Manuel Weber
% date:     2020-09-02
%
%% Function begins
%initialize
xyz = zeros(size(rhv));

switch Type
    case 'left'     % left-hand system
        xyz(:,1) = rhv(:,1).*sin(rhv(:,3)).*sin(rhv(:,2));
        xyz(:,2) = rhv(:,1).*sin(rhv(:,3)).*cos(rhv(:,2));
        xyz(:,3) = rhv(:,1).*cos(rhv(:,3));
        
    case 'right'    % right-hand system
        xyz(:,1) = rhv(:,1).*sin(rhv(:,3)).*cos(rhv(:,2));
        xyz(:,2) = rhv(:,1).*sin(rhv(:,3)).*sin(rhv(:,2));
        xyz(:,3) = rhv(:,1).*cos(rhv(:,3));
end