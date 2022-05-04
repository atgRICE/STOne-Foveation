function [wavelength_binned,M_binned] = Bin_Data(wavelength,M,numbins,lambda_lo,lambda_hi)
%BIN_DATA Downsamples spectrometer data by binning.

% Interpolation Method
    % Smooth Data
M_smooth = M; %smoothdata(M,'Gaussian',40);
% M_smooth = smoothdata(M,'movmedian',100);
    % Interpolate Smoothed Measurements
n = 0;
while n*numbins < length(wavelength(wavelength<lambda_hi & wavelength>lambda_lo))
    n = n + 1;
end
% Put onto grid that aligns with target binning.
wavelength_interp = linspace(lambda_lo,lambda_hi,n*numbins);
M_interp = interp1(wavelength,M_smooth,wavelength_interp);
if size(M_interp,1)==1
    M_interp = M_interp';
end

for ii = 1:numbins
    M_binned(:,ii) = sum(M_interp((ii-1)*n+1:ii*n,:),1,'omitnan');
    wavelength_binned(ii) = mean(wavelength_interp((ii-1)*n+1:ii*n),'omitnan');
end

end
