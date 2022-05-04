function SaveGrayVid_gif(vid,filename)
%SAVEGRAYVID_GIF Saves a grayscale video as a .gif file.

% Get number of frames in video
numFramesToMake = size(vid,3);
delaytime = 0.5;
% Save as gif
for i = 1:numFramesToMake
    % Mat2gray for display
    vid_clean(:,:,i) = uint8(255*mat2gray(vid(:,:,i)));
    if i == 1
        imwrite(vid_clean(:,:,i),filename,'gif','LoopCount',Inf,'DelayTime',delaytime);
    else
        imwrite(vid_clean(:,:,i),filename,'gif','WriteMode','append','DelayTime',delaytime);
    end
end

end

