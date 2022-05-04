function SaveSpectralSliceGif(rgbpic,filename)
%SaveHSvid_gif Saves rgb video as a .gif file.

numSpectra = size(rgbpic,4);
h = figure;
% axis tight manual
for i = 1:numSpectra
    imshow(mat2gray(rgbpic(:,:,:,i)));
%     imshow(rgbpic(:,:,:,i));
    set(gca,'color','none')
    frame = getframe(h);
    pic = frame2im(frame);
    [pic_ind,cmap] = rgb2ind(pic,256);
    
    if i == 1
        imwrite(pic_ind,cmap,filename,'gif','LoopCount',inf,'DelayTime',0.3);
    else
        imwrite(pic_ind,cmap,filename,'gif','DelayTime',0.3,'WriteMode','append');
    end
end
close(h)
end

