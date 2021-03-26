clear
clc
load('SWF')
sizeia=size(SWF);
pixelsizem=(0.000188288458*10);%10 times smaller images

%% removing irregularities with SW threshold
for k=1:2
    for i=1:sizeia(1)
        for j=1:sizeia(2)
            if SWF(i,j,k)<=15
                SWF(i,j,k)=0;
            elseif SWF(i,j,k)>=80
                SWF(i,j,k)=100;
            end
        end
    end
end
%% removing irregularities with relation to wedge position 
for k=1:2
    for j=1:sizeia(2)
        i=50;
        edge=0;
        while edge==0 && i<=sizeia(1)
            if SWF(i,j,k)>=80
                edge=1;
                SWF(i+1:end,j,k)=100;
            end
            i=i+1;
        end
     end
end

%% Plot SW
for k=1:2
subplot(2,1,k)
imagesc([0 sizeia(1,2)]*pixelsizem,[0 sizeia(1,1)]*pixelsizem,flipud(SW(:,:,k)))
   set(gca,'YDir','Normal')

   axis equal
   axis tight
   
   colormap(jet(256))
   caxis([0 100])
   c = colorbar;
   xlabel('X(m)')
   ylabel('Z(m)')
   text('Units','points','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'String','SW concentration %',...
    'Position',[360 50 0]);
end