function MakeSpectralMontage(hsvid,lambda,bandsToUse,doSave,filename,delayTime)
%MAKECOLORMONTAGE Makes a (video) montage where each image is the
%appropriate color.

% Name important values
numFrames = size(hsvid,4);
rgbvid = zeros(size(hsvid,1),size(hsvid,2),3,length(bandsToUse),numFrames);

% Convert video to RGB
% Formatting on the rgbvid variable: 
% [x, y, rgb channels, HS lambda associated with rgb, Frame]
for i = 1:numFrames
    count = 1;
    for j = bandsToUse
        OneSpec = zeros(size(hsvid,1),size(hsvid,2),size(hsvid,3));
        OneSpec(:,:,j) = hsvid(:,:,j,i);
        rgbvid(:,:,:,count,i) = HSI2RGB_convert_vis(OneSpec,lambda,1);
        count = count + 1;
    end
end

% Make montage of gifs, where each one represents one of the bandsToUse
% wavelengths. In this way, each gif is an RGB representation of that
% particular wavelength's intensity across the image.

% Display Montage
global ready
ready = 0;
figure;
H1 = uicontrol('Style','pushbutton','String','Exit',...
               'position',[25 5 50 30],'Callback',@(src,evnt)makeReady);
while(1)
    for i = 1:numFrames
        montage(mat2gray(rgbvid(:,:,:,:,i)))
        drawnow
        pause(0.001);
        if ready == 1
            break
        end
    end
    if ready == 1
        break
    end
end
close;

%% Save as gif
if doSave
    SaveColorMontageVid_gif(rgbvid,filename,delayTime)
end

end

%% --- DEPENDENCY FUNCTION DECLARATIONS --- %%
% Tells the code that you are finished viewing the video.
function makeReady(~,~)
global ready;
ready = 1;
end