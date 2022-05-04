function SaveColorMontageVid_gif(rgbvid,filename,delayTime)
%SaveHSvid_gif Saves rgb video as a .gif file.

numFrames = size(rgbvid,5);
h = figure;
axis tight manual
for i = 1:numFrames
    montage(mat2gray(rgbvid(:,:,:,:,i)));
%     montage(rgbvid(:,:,:,:,i));
    set(gca,'color','none')
    frame = getframe(h);
    pic = frame2im(frame);
    [pic_ind,cmap] = rgb2ind(pic,256);
    
    if i == 1
        imwrite(pic_ind,cmap,filename,'gif','LoopCount',inf,'DelayTime',delayTime);
    else
        imwrite(pic_ind,cmap,filename,'gif','DelayTime',delayTime,'WriteMode','append');
    end
end
close(h)
end

