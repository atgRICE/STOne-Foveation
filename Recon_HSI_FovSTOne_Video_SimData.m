clear
%close all

%% Parameters
res_full = 128; % Cannot be changed. Resolution of DMD patterns.
res_prev = 64; % Can be changed. Resolution outside of ROI. Must be 2^n.
c = log2(res_full);
c_prev = log2(res_prev);
dataPerFrame = 1024; % Number of measurements to use per frame.
shiftPerFrame = dataPerFrame; % How many meas. to shift by between frames.
% In practice, shiftPerFrame can be different from dataPerFrame. In
% simulation, it introduces artifacts due to the discreteness of the
% object's motion. Thus, we set it to dataPerFrame here.
foveate = 1; % Set to 0 for full-res reconstruction.
mu = 10;

%% Load Simulation Video
tic
vidLoc = './Utilities/ForHSVid_Simulation/ball/HSI/';
sz = [207,471,16];
% Get all files in folder.
files = dir(vidLoc);
files = files(3:end);
for i = 1:length(files)
    if strcmp(files(i).name(end-2:end),'png')
        files_img(i) = files(i);
    else
        file_boxLocs = files(i);
    end
end
% Read in images, convert to datacube.
hsi = zeros(sz(1),sz(2),sz(3),length(files_img));
for i = 1:length(files_img)
    img = imread([vidLoc,'/',files_img(i).name]);
    hsi(:,:,:,i) =  X2Cube(img);
end
times.loadVid = toc;
fprintf(['Time to Load Simulation Video: ',num2str(times.loadVid),' seconds.\n']);

%% Get Part of Video
vid = hsi(30:157,235-64:235+63,:,151:190);
numFrames = size(vid,4);
numBands = size(vid,3);
lambda = linspace(470,620,16);

% Display
for i = 1:numFrames
    for j = 1:numBands
        vid_rs(:,:,j,i) = kron(vid(:,:,j,i),ones(4));
    end
end
for i = 1:numFrames
    vid_rgb(:,:,:,i) = HSI2RGB_convert_vis(vid_rs(:,:,:,i),lambda,1);
end

%% Simulate CS-Measurement Acquisition
tic
% Encode 4D-Video as matrix.
for i = 1:numFrames
    tmp(:,:,(i-1)*numBands+1:i*numBands) = vid(:,:,:,i);
end
order = createOrderingData(res_full,'full');
for i = 1:size(tmp,3)
    vid_enc(:,i) = imageToNestedVector(tmp(:,:,i),order);
end
% STOne Transform
meas = STO(vid_enc);
% Sampling Order Permutation
meas = meas(order.samplingOrder,:);
% Sort out the measurements we actually acquired.
data = zeros(numFrames*dataPerFrame,numBands);
for i = 1:numFrames
    FrameX = meas(:,(i-1)*numBands+1:i*numBands);
    FrameX = circshift(FrameX,-(i-1)*dataPerFrame,1);
    data((i-1)*dataPerFrame+1:i*dataPerFrame,:) = FrameX(1:dataPerFrame,:);
end
times.simData = toc;
fprintf(['Time to Simulate Measurement: ',num2str(times.simData),' seconds.\n']);

%% Data Preprocessing
tic
% Define b and R
b = zeros(res_full^2,numBands*numFrames);
R = b;
rowsToSample = order.samplingOrder;
start = 1;
stop = dataPerFrame;
for f = 1:numFrames
    rowsInThisFrame = rowsToSample(1:dataPerFrame);
    b(rowsInThisFrame,(f-1)*numBands+1:f*numBands) = data(start:stop,:);
    R(rowsInThisFrame,(f-1)*numBands+1:f*numBands) = 1;
    start = start+shiftPerFrame;
    stop = stop+shiftPerFrame;
    rowsToSample = circshift(rowsToSample,-shiftPerFrame,1);
end
rowsToSample = order.samplingOrder;

clear tmp tmp_mean count*
times.preprocess = toc;
fprintf(['Time to Preprocess Data: ',num2str(times.preprocess),' seconds.\n']);

%% Generate Preview
tic
prev = makePreview(R,b,order);
% Decode Temporal-Spectral Info
prev = reshape(prev,size(prev,1),size(prev,2),numBands,numFrames);
for i = 1:numFrames
    for j = 1:numBands
        prev_rs(:,:,j,i) = kron(prev(:,:,j,i),ones(res_full/size(prev,1)));
    end
end
for i = 1:numFrames
    prev_rgb(:,:,:,i) = HSI2RGB_convert_vis(prev_rs(:,:,:,i),lambda,1);
    for j = 1:3
        prev_rgb_rs(:,:,j,i) = kron(prev_rgb(:,:,j,i),res_full/size(prev,1));
    end
end
times.getPreview = toc;
fprintf(['Time to Get Preview: ',num2str(times.getPreview),' seconds.\n']);

%% Select ROIs
if foveate
%     %%% Uncomment to pick your own set of ROIs! %%%
%     % Set up GUI
%     clear TL BR
%     global done;
%     global ready;
%     done = 0;
%     ready = 0;
%     roiSelectionPic = figure; hold on
%     imshow(mat2gray(prev_rgb_rs(:,:,:,1)));
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
%     %%% End "Uncomment to pick your own set of ROIs!" %%%
    
    % For Reproducing Paper Result.
    % Comment this out if own ROI selected.
    numROI = 5;
    TL = [1,39; 23,71; 43,53; 47,75; 67,47];
    BR = [30,98; 50,102; 70,82; 82,102; 88,80];
    
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
    times.fovParams = toc;
    fprintf(['Time to Get Foveation Parameters: ',num2str(times.fovParams),' seconds.\n']);
    
    %% Reconstruction
    tic
    % Foveated Reconstruction
    [vec_recon,~] = pdhg_HSVid_l1_fov(R,b,mu,numBands,numFrames,order,...
                                      res_prev,TL,BR,record_ds,fovParams);
    % Need to upsample the LR region before going back to an image.
    for j = 1:size(vec_recon,2)
        count = 1;
        for i = 1:length(record_ds)
            if record_ds(i)
                vec_recon_us(count:count+jump-1,j) = vec_recon(i,j);
                count = count + jump;
            else
                vec_recon_us(count,j) = vec_recon(i,j);
                count = count + 1;
            end
        end
    end
    times.FovRecon = toc;
    fprintf(['Time for DS Foveated Reconstruction: ',num2str(times.FovRecon),' seconds.\n']);
else
    % Full-Res Reconstruction
    [vec_recon_us,~] = pdhg_hypSpecVideo_l1(R,b,mu,numBands,numFrames,order);
    times.FullResRecon = toc;
    fprintf(['Time for DS Foveated Reconstruction: ',num2str(times.FullResRecon),' seconds.\n']);
end

% Now convert upsampled vec_recon into an image.
for i = 1:size(vec_recon_us,2)
    rec_vid_mixed(:,:,i) = nestedVectorToImage(vec_recon_us(:,i),order);
end
for f = 1:numFrames
    rec(:,:,:,f) = rec_vid_mixed(:,:,(f-1)*numBands+1:f*numBands);
end

% Post-Processing (optional)
rec_clean = cleanImage_Thresh(rec,0.02);
% Display
for i = 1:numFrames
    for j = 1:numBands
        rec_clean_rs(:,:,j,i) = kron(rec_clean(:,:,j,i),ones(4));
    end
end
for i = 1:numFrames
    rec_rgb(:,:,:,i) = HSI2RGB_convert_vis(rec_clean_rs(:,:,:,i),lambda,1);
end
figure;
for i = 1:numFrames
    imshow(mat2gray(rec_rgb(:,:,:,i)));
    pause(0.1)
end
title([num2str(res_prev),'/',num2str(res_full),' DS Foveated Reconstruction'])
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