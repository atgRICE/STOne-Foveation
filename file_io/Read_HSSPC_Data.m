function [wavelength_binned,M_binned] = Read_HSSPC_Data(filename,background_filename,numbins,lambda_lo,lambda_hi)
%Read_HSSPC_Data Reads in data from HyperSpectral Imager.

% Get Data
M = importdata(filename);
% Get Background
bg = dlmread(background_filename);
% Pull out wavelengths
wavelength = M(2:end,1);
% Remove wavelengths and timestamps
M = M(2:end,2:end);
% Background-subtraction
M = M - bg(:,2);
% Bin Data
[wavelength_binned,M_binned] = Bin_Data(wavelength,M,numbins,lambda_lo,lambda_hi);

end
