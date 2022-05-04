function [infostruct] = HSI_SpectralExtractor(hsi,band_w,pictype)
%HSI_SPECTRALEXTRACTOR Easy hyperspectral image spectrum extraction.
%   Accepts the hyperspectral datacube, a vector of wavelengths associated
%   with the bandwidth of the datacube, and a string describing whether the
%   image displayed is 'pseudocolor' or 'intensity'.
%   ---
%   First, a displayable image is generated from the datacube by either:
%   1) Computing a pseudocolor RGB version of the image, or
%   2) Summing along the spectral dimension to generate an intensity image.
%   This image is displayed and 3 buttons are shown at the bottom (2 if a
%   pseudocolor image is generated). The "Ready" button is clicked
%   immediately before a pixel (or group of pixels) is selected for 
%   spectral extraction. Any zooming should be done before pressing 
%   "Ready". The "Done" button is pressed when you are finished extracting
%   spectra. Once this is pressed, the code returns the archive of selected
%   points and a structure of their locations in the image.
%   ---
%   To extract spectra, click "Ready", then click on a group of pixels, the
%   spectra of which you would like to average. Shift-click (or
%   double-click) on the final pixel. If only a single pixel is selected,
%   only that spectrum is plotted and saved. If multiple are selected, then
%   they are averaged, plotted, and saved.
%   --- 
%   pictype = 'intensity' makes an intensity image.
%   pictype = 'pseudo' makes a pseudocolor image from visible data.
%
%   By Anthony Giljum, Kelly Lab, Rice University, February 2019.

objnum = 0;
count = 0;
global done
global notready
clear x y

% Make Plot
    % What type of plot do we want?
if strcmp(pictype,'intensity') == 1
        % Sum bands
    Im_sumup = sum(hsi,3);
        % Plot Sum
    im = figure;
    imagesc(Im_sumup), colorbar, colormap hot, axis image
    title('Intensity Image')
elseif strcmp(pictype,'pseudo') == 1
    % Generate pseudocolor image from visible spectral data
    colorpic = HSI2RGB_convert_vis(hsi,band_w,1);
    im = figure;
    imshow(mat2gray(colorpic)); axis image; title('PseudoColor Image')
end

% Make buttons.
H1 = uicontrol('Style','pushbutton','String','Ready',...
               'position',[25 5 50 30],'Callback',@(src,evnt)MakeReady);
           
if strcmp(pictype,'intensity') == 1
    H2 = uicontrol('Style','pushbutton','String','Change Caxis',...
                   'position',[450,5,100,30],...
                   'Callback',@(src,evnt)Change_caxis);
end

H3 = uicontrol('Style','pushbutton','String','Done',...
               'position',[200,5,75,30],'Callback',@(src,evnt)AmDone);

clims  = caxis;

% Primary Code
    % While we want to extract spectrum...
while(1)
    % Get correct axis handle.
    figure(im);
    % Do any zooming
    notready = 1;
    done = 0;
    while notready == 1 && done == 0
        pause(0.001);
    end
    if done == 1
        break
    end
    % Get point coordinates
    [y,x,~] = impixel;
    hold on
    plot(y,x,'go')
    count = count + 1;
    
    % Label the point coordinates so we know their associated plots
    for i = 1:length(x)
        text(y(i),x(i),num2str(objnum+1),'Color','c');
    end
    
    % Define particles by grouping
    for i = 1:length(x)
        p(i,:) = hsi(x(i),y(i),:);
        infostruct(objnum+1).location(i,1) = x(i);
        infostruct(objnum+1).location(i,2) = y(i);
        infostruct(objnum+1).spectrum(i,:) = p(i,:);
    end
    
    % Average up grouped particles
    spec = sum(p,1)/length(x);
    
    % Get correct axis handle.
    if mod(objnum,6) == 0
        specplot = figure;
    else
        figure(specplot)
    end
    
    % Plot Spectrum
    subplot(2,3,mod(objnum,6)+1)
    axishandles(count) = subplot(2,3,mod(objnum,6)+1);
    % Set axes for all plots
    if objnum == 0
        currentmax = max(spec);
        currentmin = min(spec);
    end
    if max(max(spec)) > currentmax
        currentmax = max(spec);
    end
    if min(min(spec)) < currentmin
        currentmin = min(spec);
    end
    plot(band_w,spec);
    title(['Point ',num2str(objnum+1)])
    axis(axishandles,[min(band_w),max(band_w),currentmin,currentmax])

    % Refresh Variables, Update object number.
    clear x y p
    spec_archive(:,objnum+1) = spec;
    objnum = objnum + 1;
end
end

%% --- Function Dependency Declarations --- %%
% When button clicked, calls this function, which tells the program you are
% ready to extract a nanoparticle's spectrum.
function MakeReady(~, ~)
    global notready
    notready = 0;
end

function Change_caxis(~,~)
    cur_lims = caxis;
    prompt   = {'Enter Lower ColorAxis Limit','Enter Upper ColorAxis Limit'};
    title    = 'Change ColorAxis';
    dims      = [1 35];
    definput = {num2str(cur_lims(1)),num2str(cur_lims(2))};
    
    answer = inputdlg(prompt,title,dims,definput);
    answer = [str2double(answer{1}),str2double(answer{2})];

    caxis([answer(1),answer(2)])
end

function AmDone(~,~)
    global done
    done = 1;
end