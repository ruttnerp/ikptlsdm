function [fitresult, gof] = createFit(z_all, time_all, z_v2dp)
%CREATEFIT(Z_ALL,TIME_ALL,Z_V2DP)
%  Create a fit.
%
%  Data for 'z_time_v2dp_2D' fit:
%      X Input: z_all
%      Y Input: time_all
%      Z Output: z_v2dp
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 16-Jun-2022 16:36:35


%% Fit: 'z_time_v2dp_2D'.
[xData, yData, zData] = prepareSurfaceData( z_all, time_all, z_v2dp );

% Set up fittype and options.
ft = fittype( 'poly33' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'LAR';

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

% Plot fit with data.
figure( 'Name', 'z_time_v2dp_2D' );
h = plot( fitresult, [xData, yData], zData );
legend( h, 'z_time_v2dp_2D', 'z_v2dp vs. z_all, time_all', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'z_all', 'Interpreter', 'none' );
ylabel( 'time_all', 'Interpreter', 'none' );
zlabel( 'z_v2dp', 'Interpreter', 'none' );
grid on
view( -176.3, 17.7 );

