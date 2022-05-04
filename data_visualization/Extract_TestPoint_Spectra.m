function spectra = Extract_TestPoint_Spectra(hsi,lambda)
%% Extract Spectra of Test Points
% close all
points(1,:) = [63,33];
points(2,:) = [77,36];
points(3,:) = [62,88];

implot = figure;
% imagesc(sum(hsi,3));axis image;colormap hot;colorbar
rgbpic = HSI2RGB_convert_vis(hsi,lambda,1);
imshow(mat2gray(rgbpic));


for i = 1:size(points,1)
    % Mark where Point was Taken
    figure(implot)
    hold on
    plot(points(i,2),points(i,1),'go')
    text(points(i,2)+5,points(i,1)+1,num2str(i),'Color','c','fontsize',16)
    
    % Spectrum Extraction
    spec = squeeze(hsi(points(i,1),points(i,2),:));
    spectra(:,i) = spec;
    % Plot Spectrum
    if mod(i-1,3) == 0
        specplot = figure;
    else
        figure(specplot)
    end
    subplot(3,1,mod(i,3)+1);
    axishandles(i) = subplot(3,1,mod(i-1,3)+1);
    if i == 1
        currentmax = max(spec);
        currentmin = min(spec);
    end
    if max(max(spec)) > currentmax
        currentmax = max(spec);
    end
    if min(min(spec)) < currentmin
        currentmin = min(spec);
    end
    plot(lambda,spec);
    title(['Point ',num2str(i),' - (',num2str(points(i,1)),',',num2str(points(i,2)),')'])
    axis(axishandles,[min(lambda),max(lambda),currentmin,currentmax]);
%     axis(axishandles,[min(lambda),max(lambda),0,1]);
    set(gca,'fontsize',14)

end

end