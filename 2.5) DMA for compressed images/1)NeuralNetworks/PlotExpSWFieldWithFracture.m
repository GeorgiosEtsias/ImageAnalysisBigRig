load('SWRightTest2')
SW=SWRightTest2;

%% Rectangel characteristics
x1=2134;
y1=729;
y1=1185-y1;
length=1373;
width=80;
rectangleChar=[x1 y1 length width];

sizeia=size(SW);
pixelsizem=0.0001054384;
%% Plot SW field
figure (1)
%imagesc([0 sizeia(1,2)]*pixelsizem,[0 sizeia(1,1)]*pixelsizem,flipud(SW(:,:,end)))
imagesc(flipud(SW(:,:,end)))
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
hold on
%% Plotting the fracture posssition 
 r=rectangle('Position',rectangleChar,'FaceColor',[0 0 0],'EdgeColor','w',...
    'LineWidth',2);
%r.LineStyle=':'; %dotted lines do not look good
