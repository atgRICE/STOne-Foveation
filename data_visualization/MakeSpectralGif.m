function MakeSpectralGif(hspic,lambda,doSave,filename)
%MAKESPECTRALGIF Makes a gif that scans over each spectrum in the correct
%color.

% Name important values
numSpectra = size(hspic,3);
rgbpic = zeros(size(hspic,1),size(hspic,2),3,numSpectra);

% Convert video to RGB
% Formatting on the rgbvid variable: 
% [x, y, rgb channels, HS lambda associated with rgb]
for i = 1:numSpectra
    OneSpec = zeros(size(hspic));
    OneSpec(:,:,i) = hspic(:,:,i);
    rgbpic(:,:,:,i) = HSI2RGB_convert_vis(OneSpec,lambda,1);
end

% Display Montage
global ready
ready = 0;
figure;
H1 = uicontrol('Style','pushbutton','String','Exit',...
               'position',[25 5 50 30],'Callback',@(src,evnt)makeReady);
while(1)
    for i = 1:numSpectra
        imshow(mat2gray(rgbpic(:,:,:,i)))
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
    SaveSpectralSliceGif(rgbpic,filename)
end

end

%% --- DEPENDENCY FUNCTION DECLARATIONS --- %%
% Tells the code that you are finished viewing the video.
function makeReady(~,~)
global ready;
ready = 1;
end