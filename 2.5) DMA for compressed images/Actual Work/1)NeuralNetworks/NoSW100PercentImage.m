SWtest=PredictionC0;
for k=1:sizeia(3)
for i=1:sizeia(1)
for j=1:sizeia(2)
if SWtest(i,j,k)>=75
if SWtest(i,j,k)>=85
SWtest(i,j,k)=100;
else
SWtest(i,j,k)=(SWtest(i,j,k)/85)*100;
end
end
end
end
end
for i=1:sizeia(3)
figure (i)
imagesc([0 sizeia(1,2)]*pixelsizem,[0 sizeia(1,1)]*pixelsizem,flipud(SWtest(:,:,i)))
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