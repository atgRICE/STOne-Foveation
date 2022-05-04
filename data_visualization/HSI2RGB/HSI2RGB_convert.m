function RGB = HSI2RGB_convert(Im_recon, band_w, scale)
% convert HS data cube to RGB image
% input: Im_recon n*n*numBands

load('C:\Users\Track\Desktop\Lab_Stuff\KentOptronics_Stuff\Matlab\HSI2RGB\xyzbar.mat'); % 33 channels for 400:10:720 nm
xyzbar_interp = interp1([400:10:720]',xyzbar,band_w);
% xyzbar_interp = interp1([400:10:720]',xyzbar,band_w(18:65)); % 100 band

[r,c,w] = size(Im_recon);
Im_tmp = reshape(Im_recon, r*c, w);
Im_tmp = Im_tmp(:,:);    % 33 Band


XYZ = Im_tmp*xyzbar_interp;
XYZ = reshape(XYZ, r, c, 3);
XYZ = max(XYZ, 0);
%XYZ = XYZ/max(XYZ(:));
XYZ = XYZ/scale;

RGB = XYZ2sRGB_exgamma(XYZ);
% RGB = max(RGB, 0);
% RGB = min(RGB, 1);

end