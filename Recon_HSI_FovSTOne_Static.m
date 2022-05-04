clear
close all

%% Parameters
res_full = 128; % Cannot be changed. This is the DMD pattern resolution.
res_prev = 32; % Can be changed. Resolution outside ROI. Must be 2^n for int. n.
c = log2(res_full);
c_prev = log2(res_prev);
numBins = 64; % Can be changed. Gives number of spectral bins.
mu = 5; % Lower mu -> Don't trust meas, trust regularization.
        % Higher mu -> Don't trust regularization, trust meas.
specWeight = 1; % How much to regularize each pixel's spectrum.
lambda_lo = 500; % Lowest wavelength to reconstruct.
lambda_hi = 700; % Highest wavelength to reconstruct.
numMeas = 3277; % Number of measurements to reconstruct. Must be <=3277.
foveate = 1;
filepath = './';
filedesc = 'T1_ChartAndFan_128sto_50ms_Flame.txt';
filename(1).name = [filepath,filedesc];
background_filename = [filepath,'T1_bg.txt'];
norm_filename = [filepath,'T1_allBlack.txt'];

%% Read in Data
data = [];
for i = 1:length(filename)
    % This I/O code reads data from the Flame spectrometer.
    % It verifies that there are no missing measurements from the
    % timestamps in the first row, or else replaces missing measurements
    % with 0's and returns the indices of the missing measurements. It is
    % also formatted specifically for the Synchronous Trigger mode of the
    % Flame spectrometer from OceanOptics.
    [lambda,meas] = Read_HSSPC_Data(filename(i).name,...
                                                background_filename,...
                                                numBins,lambda_lo,...
                                                lambda_hi);
    data = cat(3,data,meas);
end

%% Data Normalization
normData = load(norm_filename);
bg = load(background_filename);
normData(:,2) = normData(:,2) - bg(:,2);
[~,normData] = Bin_Data(normData(:,1)',normData(:,2)',numBins,lambda_lo,lambda_hi);
data = data./normData;
clear bg

%% Data Preprocessing
order = createOrderingData(res_full,'full');
rowsToSample = order.samplingOrder(1:numMeas);
% Define b and R
b = zeros(res_full^2,numBins);
R = b;
b(rowsToSample,:) = data(1:numMeas,:);
R(rowsToSample,:) = 1;

% Account for >50% light acceptance of STO pats and 1/0 -> +/- conversion.
count_1 = 1; count_0 = 0;
for ii = 1:round(log2(res_full))
    count_1_new = count_1*3 + count_0;
    count_0_new = count_1 + count_0*3;
    count_1 = count_1_new; count_0 = count_0_new;
end
for ii = 1:size(b,2)
    tmp = b(:,ii);
    tmp = tmp(R(:,ii)==1);
    tmp_mean = mean(tmp);
    tmp_mean = tmp_mean/count_1*(count_1+count_0)/2;
    tmp = tmp - tmp_mean;
    b(R(:,ii)==1,ii) = tmp;
    b(:,ii) = b(:,ii)*2^(c-1);
end

clear tmp tmp_mean count*

%% Select ROIs
if foveate
    %%% Uncomment to pick your own set of ROIs! %%%
    % % Make preview
    % tmp = makePreview(R,b,order);
    % tmp = HSI2RGB_convert_vis(tmp,lambda,1);
    % for i = 1:3
    %     prev(:,:,i) = kron(tmp(:,:,i),ones(res_full/size(tmp,1)));
    % end
    % 
    % % Set up GUI
    % clear TL BR
    % global done;
    % global ready;
    % done = 0;
    % ready = 0;
    % roiSelectionPic = figure; hold on
    % imshow(mat2gray(prev));
    % H1 = uicontrol('Style','pushbutton','String','Done',...
    %               'position',[25 5 50 30],'Callback',@(src,evnt)isDone);
    % H2 = uicontrol('Style','pushbutton','String','Select Region',...
    %               'position',[85 5 80 30],'Callback',@(src,evnt)isReady);
    % % User selects multiple ROIs. Get TL and BR from that.
    % count = 1;
    % while ~done
    %     pause(0.01);
    %     if ready
    %         hold on
    %         % Get region.
    %         rect = getrect(roiSelectionPic);
    %         TL(count,:) = round([max(rect(1),1),max(rect(2),1)]);
    %         BR(count,:) = round([min(rect(1)+rect(3),res_full),...
    %                              min(rect(2)+rect(4),res_full)]);
    %         % Ensure you don't split pixels.
    %         for i = 1:2
    %             modval = mod(TL(count,i),res_full/res_prev);
    %             if modval ~= 1
    %                 TL(count,i) = TL(count,i)-modval+1;
    %             end
    %             modval = mod(BR(count,i),res_full/res_prev);
    %             if modval ~=0
    %                 BR(count,i) = BR(count,i)-modval+(res_full/res_prev);                
    %             end
    %         end
    %         % Mark selected ROI on the image.
    %         plot(TL(count,1),TL(count,2),'ro');
    %         plot(TL(count,1),BR(count,2),'ro');
    %         plot(BR(count,1),TL(count,2),'ro');
    %         plot(BR(count,1),BR(count,2),'ro');
    %         plot([TL(count,1),TL(count,1)],[TL(count,2),BR(count,2)],'r-');
    %         plot([TL(count,1),BR(count,1)],[TL(count,2),TL(count,2)],'r-');
    %         plot([BR(count,1),BR(count,1)],[BR(count,2),TL(count,2)],'r-');
    %         plot([BR(count,1),TL(count,1)],[BR(count,2),BR(count,2)],'r-');
    %         
    %         count = count + 1;
    %     end
    %     ready = 0;
    % end
    % numROI = size(TL,1);
    % TL(:,[1,2]) = TL(:,[2,1]);
    % BR(:,[1,2]) = BR(:,[2,1]);
    %%% End "Uncomment to pick your own set of ROIs!" %%%

    % For Reproducing Paper Result.
    % Comment this out if own ROI selected.
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
    [vec_recon,~] = pdhg_HSI_l1_fov(R,b,mu,specWeight,order,res_prev,TL,BR,...
                                    record_ds,fovParams);
    times.FoveaRecon = toc;
    fprintf(['Time for DS Foveated Reconstruction: ',num2str(times.FoveaRecon),' seconds.\n']);
    % Need to upsample the LR region before going back to an image.
    vec_recon_us = UndoRLE(vec_recon,record_ds,fovParams);
else
    tic
    [vec_recon_us,~] = pdhg_HSI_l1(R,b,mu,specWeight,order);
    times.FullResRecon = toc;
    fprintf(['Time for Full-Res Reconstruction: ',num2str(times.FullResRecon),' seconds.\n']);
end
% Now convert upsampled vec_recon into an image.
tic
for i = 1:length(lambda)
    rec(:,:,i) = nestedVectorToImage(vec_recon_us(:,i),order);
end
times.NE2I = toc;
fprintf(['Time to Unpermute Image: ',num2str(times.NE2I),' seconds.\n']);

% Post-Processing (optional)
rec_clean = cleanImage_Thresh(rec,0.001);

% Display
figure;imshow(mat2gray(HSI2RGB_convert_vis(rec_clean,lambda,1)));
title([num2str(res_full),'/',num2str(res_prev),' DS Foveated Reconstruction'])
set(gca,'fontsize',14);

%% Get Spectra
spectra_fov = Extract_TestPoint_Spectra(rec_clean,lambda);

% To select your own spectra, use HSI_SpectralExtractor.
% When using this script, click "Ready" to select pixels. Once the cursor
% changes to a crosshair, click all pixels whose spectra you would like to
% average. Shift click to select the final pixel for a single curve. The
% average of all selected spectra will plotted and a new set can be chosen.
% Click "Done" when you are finished to output all selected spectra to the
% variable "spectra_fov". This output structure contains all selected
% points before averaging, as well as their respective spectra before
% averaging.
% See the comments at the beginning of the HSI_SpectralExtractor script for
% more details.

% spectra_fov = HSI_SpectralExtractor(rec_clean,lambda,'pseudo');

%% --- FUNCTION DECLARATIONS --- %%
function isDone(~,~)
global done
done = 1;
end

function isReady(~,~)
global ready
ready = 1;
end