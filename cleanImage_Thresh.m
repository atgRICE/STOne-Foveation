function [pic_clean] = cleanImage_Thresh(pic,perc)
%CLEANIMAGE_Thresh Cleans L1 reconstructed images. Support for hyperspectral.
% Finds the top and bottom (perc)% intensities in the full datacube, then hard
% thresholds them.

% Sort data
pic_sorted = sort(pic(:),'ascend');
% Find number of points to threshold on each end.
num_to_thresh = round(perc*length(pic(:)));
% Find low and high threshold values.
low_val = pic_sorted(num_to_thresh);
high_val = pic_sorted(end-num_to_thresh);
% Threshold Image
pic_clean = pic;
pic_clean(pic_clean>high_val) = high_val;
pic_clean(pic_clean<low_val) = low_val;

end

