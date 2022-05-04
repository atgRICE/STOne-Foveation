function RGB = HSI2RGB_convert_vis(Im_recon, band_w, scale)
% convert HS data cube to RGB image
% input: Im_recon n*n*numBands

load('C:\Users\agilj\Desktop\Lab_Stuff\Extra_Matlab_Code\HSI2RGB\xyzbar.mat'); % 33 channels for 400:10:720 nm
visband = band_w(band_w>=400 & band_w<=720);
xyzbar_interp = interp1([400:10:720]',xyzbar,visband);

[r,c,w] = size(Im_recon);
Im_tmp = reshape(Im_recon, r*c, w);
Im_tmp = Im_tmp(:,1:length(visband));

XYZ = Im_tmp*xyzbar_interp;
XYZ = reshape(XYZ, r, c, 3);
XYZ = max(XYZ, 0);
%XYZ = XYZ/max(XYZ(:));
XYZ = XYZ/scale;

RGB = XYZ2sRGB_exgamma(XYZ);
% RGB = max(RGB, 0);
% RGB = min(RGB, 1);

end