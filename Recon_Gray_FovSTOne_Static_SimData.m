clear
close all

%% Parameters
% Remember, the larger the difference between c and c_prev, the greater the
% error introduced by the downsampling approximation. This can be partially
% alleviated by reducing the value of mu for large differences in
% resolution.
res_full = 512; % Cannot be changed. This is the DMD pattern resolution.
res_prev = 256; % Can be changed. Resolution outside ROI. Must be 2^n for int. n.
c = log2(res_full);
c_prev = log2(res_prev);
mu = 10; % Lower mu -> Don't trust meas, trust regularization.
        % Higher mu -> Don't trust regularization, trust meas.
csRat = 0.1; % Percentage of the total number of meas. to use.
numMeas = floor(csRat*res_full^2); % Number of measurements to reconstruct.
foveate = 1;

%% Read in Data
pic = imread('monarch_test_img.png');
pic = rgb2gray(pic);
pic = pic(:,1:512);
pic = double(pic)/255;

%% Make CS Measurements
% Note that the procedure used here to acquire compressive measurements in
% simulation is equivalent to sampling the scene using a permuted sensing
% matrix.
% y = (SamplePermMat)*Phi*(PixelPermMat)*x
order = createOrderingData(res_full,'full');
pic_lin = imageToNestedVector(pic,order); % Pixelwise permutation
data = STO(pic_lin); % Apply full-res sensing matrix
data = data(order.samplingOrder); % Sampling order permutation
data = data(1:numMeas); % Subsampling.

%% Data Preprocessing
rowsToSample = order.samplingOrder(1:numMeas);
% Define b and R
b = zeros(res_full^2,1);
R = b;
b(rowsToSample) = data(1:numMeas);
R(rowsToSample) = 1;

% Account for >50% light acceptance of STO pats and 1/0 -> +/- conversion.
count_1 = 1; count_0 = 0;
for ii = 1:round(log2(res_full))
    count_1_new = count_1*3 + count_0;
    count_0_new = count_1 + count_0*3;
    count_1 = count_1_new; count_0 = count_0_new;
end
tmp = b;
tmp = tmp(R==1);
tmp_mean = mean(tmp);
tmp_mean = tmp_mean/count_1*(count_1+count_0)/2;
tmp = tmp - tmp_mean;
b(R==1) = tmp;
b = b*2^(c-1);

clear tmp tmp_mean count*

%% Select ROIs
if foveate
    %%% Uncomment to pick your own set of ROIs! %%%
%     % Make preview
%     tmp = makePreview(R,b,order);
%     prev = kron(tmp,ones(res_full/size(tmp,1)));
%     
%     % Set up GUI
%     clear TL BR
%     global done;
%     global ready;
%     done = 0;
%     ready = 0;
%     roiSelectionPic = figure; hold on
%     imshow(mat2gray(prev));
%     H1 = uicontrol('Style','pushbutton','String','Done',...
%                   'position',[25 5 50 30],'Callback',@(src,evnt)isDone);
%     H2 = uicontrol('Style','pushbutton','String','Select Region',...
%                   'position',[85 5 80 30],'Callback',@(src,evnt)isReady);
%     % User selects multiple ROIs. Get TL and BR from that.
%     count = 1;
%     while ~done
%         pause(0.01);
%         if ready
%             hold on
%             % Get region.
%             rect = getrect(roiSelectionPic);
%             TL(count,:) = round([max(rect(1),1),max(rect(2),1)]);
%             BR(count,:) = round([min(rect(1)+rect(3),res_full),...
%                                  min(rect(2)+rect(4),res_full)]);
%             % Ensure you don't split pixels.
%             for i = 1:2
%                 modval = mod(TL(count,i),res_full/res_prev);
%                 if modval ~= 1
%                     TL(count,i) = TL(count,i)-modval+1;
%                 end
%                 modval = mod(BR(count,i),res_full/res_prev);
%                 if modval ~=0
%                     BR(count,i) = BR(count,i)-modval+(res_full/res_prev);                
%                 end
%             end
%             % Mark selected ROI on the image.
%             plot(TL(count,1),TL(count,2),'ro');
%             plot(TL(count,1),BR(count,2),'ro');
%             plot(BR(count,1),TL(count,2),'ro');
%             plot(BR(count,1),BR(count,2),'ro');
%             plot([TL(count,1),TL(count,1)],[TL(count,2),BR(count,2)],'r-');
%             plot([TL(count,1),BR(count,1)],[TL(count,2),TL(count,2)],'r-');
%             plot([BR(count,1),BR(count,1)],[BR(count,2),TL(count,2)],'r-');
%             plot([BR(count,1),TL(count,1)],[BR(count,2),BR(count,2)],'r-');
%             
%             count = count + 1;
%         end
%         ready = 0;
%     end
%     numROI = size(TL,1);
%     TL(:,[1,2]) = TL(:,[2,1]);
%     BR(:,[1,2]) = BR(:,[2,1]);
    %%% End "Uncomment to pick your own set of ROIs!" %%%

    % For Reproducing Paper Result.
    % Comment this out if own ROI selected.
    % This requires res_full = 512, res_prev = 256 for reproducing the
    % paper code.
    % Other values can be used, but remember that you cannot have half of
    % a low-resolution pixel in the reconstruction. This is handled for
    % you in the above code -- that is, the selected foveated regions are
    % shifted to remove any split low-resolution pixels.
    numROI = 1;
    TL = [41,21];
    BR = [84,64];

    %% Get Foveated Region.
    tic
    nestvals = [];
    for kk = 1:numROI
        [xmesh,ymesh] = meshgrid(TL(kk,1):BR(kk,1),TL(kk,2):BR(kk,2));
        subvals = [xmesh(:),ymesh(:)];
        linvals = sub2ind([res_full,res_full],subvals(:,1),subvals(:,2));
        nestim = nestedVectorToImage(1:res_full^2,order);
        nestim = nestim(:);
        nestvals = [nestvals;nestim(linvals)]; % nestvals is our foveated region here.
        clear linvals nestim subvals xmesh ymesh
    end
    nestvals = unique(nestvals);

    jump = 4^(c-c_prev); % We will downsample all pixels in [i,i+jump].
    count = 1;
    clear record_ds
    for i = 1:jump:res_full^2
        % Downsample if we are outside of the foveated region. Otherwise, keep
        % all values.
        if anyEq(i,nestvals) % In Foveated Region
            record_ds(count:count+jump-1) = 0;
            count = count + jump;
        else % Outside Foveated Region
            % We already have the High-Res Downsample, so just take first
            % value.
            record_ds(count) = 1;
            count = count + 1;
        end
    end
    order.n1 = res_full^2;
    order.n2 = count-1;
    LRreg = setdiff(1:res_full^2,nestvals);
    fovParams = uint32([res_full,res_prev,sort(LRreg)]); % MUST have this form.
    times.PatternDownsample = toc;
    fprintf(['Time to Get Foveation Parameters: ',num2str(times.PatternDownsample),' seconds.\n']);

    %% Reconstruction
    tic
    w = [1,1];
    [vec_recon,~] = pdhg_image_l1_fov(R,b,mu,w,order,res_prev,TL,BR,...
                                    record_ds,fovParams);
%     % This uses an L2 measurement consistency term in reconstruction
%     % instead of an L1 measurement consistency term. As such, it does not
%     % handle the noise as well.
%     [vec_recon,~] = pdhg_image_fov(R,b,mu,order,res_prev,TL,BR,...
%                                     record_ds,fovParams);
    times.FoveaRecon = toc;
    fprintf(['Time for DS Foveated Reconstruction: ',num2str(times.FoveaRecon),' seconds.\n']);
    % Need to upsample the LR region before going back to an image.
    vec_recon_us = UndoRLE(vec_recon,record_ds,fovParams);
else
    tic
%     [vec_recon_us,~] = pdhg_image_l1(R,b,mu,order);
    [vec_recon_us,~] = pdhg_image(R,b,mu,order);
    times.FullResRecon = toc;
    fprintf(['Time for Full-Res Reconstruction: ',num2str(times.FullResRecon),' seconds.\n']);
end
% Now convert upsampled vec_recon into an image.
tic
rec = nestedVectorToImage(vec_recon_us,order);
times.NE2I = toc;
fprintf(['Time to Unpermute Image: ',num2str(times.NE2I),' seconds.\n']);

% Post-Processing (optional)
rec_clean = cleanImage_Thresh(rec,0.001);

% Display
figure;imshow(mat2gray(rec_clean));
title([num2str(res_full),'/',num2str(res_prev),' DS Foveated Reconstruction'])
set(gca,'fontsize',14);

%% --- FUNCTION DECLARATIONS --- %%
function isDone(~,~)
global done
done = 1;
end

function isReady(~,~)
global ready
ready = 1;
end