% Create a preview from undersampled data.  The method automatically
% calculates the downsampling factor by using the number of measurements.
% modified by Liyang on 07/15/2015
% modified for spectrometer-based hyperspectral by ATG on 2020 March 16.

function [preview vec]= makePreview(R, b, order )

N = size(R,1);  %  Number of pixels in high resolution
n = sum(R(:,1));     % pixels in low resolution


delta = N/n;

if delta==1
    for i = 1:size(R,2)
        vec = STO(b(:,i));
        preview(:,:,i) = nestedVectorToImage(vec, order);
    end
    return;
end

% added by Liyang
k = log2(delta); k = k/2; % delta = 4^k
k = ceil(k);
delta = 4^k;
%

for i = 1:size(R,2)
    blockData(:,:,i) = reshape(b(:,i),[delta, N/delta]);
    blockR(:,:,i) = reshape(R(:,i),[delta, N/delta]);
end
means = squeeze(sum(blockR.*blockData,1)./sum(blockR,1));
if size(means,1)==1
    means = means';
end

blocks = kron(means,ones(delta,1));

% vec = blocks(:);

vec = STO(blocks);

for i = 1:size(vec,2)
    preview(:,:,i) = nestedVectorToImage(vec(:,i), order);
end

rootDelta = round(sqrt(delta));

preview = preview(1:rootDelta:end,1:rootDelta:end,:);

end

