%-------------------G.Etsias February 28-2019-----------------------------%
%----------------Data Modification Approach (DMA)-------------------------%
% -------Deriving Prediction Data using perfect C=0% & C=100% ----------- %
% ---------------------- New Camera, Green LI --------------------------- %
% ======================================================================= %
%----Based on: G.Etsias November 20-2018 version for the old camera-------%
% ======================================================================= %
%---Derive data for training Neural Network to do SW regression analysis--%
%------------------ using 3*8=24 calibration images-----------------------%
%-----------Input variables X,Y,BS,LI(G) output SW C%---------------------%
% ======================================================================= %
clear
clc
%% --- Variables ---

% Load the 3 datasets
load('SWFreshOnlyG')
load('SPTest20G')
load('SPTest15G')
load('SWOnlyG')
% Load mean values of LI
load('M') %from the training data DMA
% Load Homogenisation Factor
load ('MeanHomoFactor') %from the training data DMA
test=1;
%% Chose Test case chareacterisitcs
BD1=780; %Bead size of the investigated homogeneous aquifer
TestCase(:,:,1)=SPTest20G;
TestCase(:,:,2)=SPTest15G;
TestNumber=size(TestCase);

% === Variables to manually set ===========================================
npts = 8; % No. of calibration images in which ANN was trained
% Number of pixels that are not full of beads and are to be
% ignored. (If sandbox full of glass pixlim=1, also pixlim=1 for intial investigations!!)
pixlim=1; % foul aquifers
% If homogenization procedur occures homog=1 if not homog=0
homog=1;
% Set plot perf =1 to plot the perfect C=0% and C=100% aquifers
plotperf=1 ;

% ======================================================================= %
sizeia = size(TestCase);
npixels = sizeia(1,1)*sizeia(1,2);

%Derive the LI data for C= 0% - 100% and bead sizes 780um-1090um-1325um
avgimgarea= zeros(sizeia(1),sizeia(2),npts); 
avgimgarea(:,:,1)= SWFreshOnlyG;
avgimgarea(:,:,2:1+TestNumber(3))= TestCase(:,:,1:end);
avgimgarea(:,:,npts)= SWOnlyG;

%% Image Homogenization
for i=1:npts
        avgimgarea(:,:,i)=avgimgarea(:,:,i)./MeanHomoFactor;
end

%% --- Cange Factor Caluclation ---- %

change1C0=zeros(size(1,1),size(1,2));
change2C0=zeros(size(1,1),size(1,2));
change3C0=zeros(size(1,1),size(1,2));
change1C100=zeros(size(1,1),size(1,2));
change2C100=zeros(size(1,1),size(1,2));
change3C100=zeros(size(1,1),size(1,2));



%Perfect C=0Applies only for C=0% images, calculating the change factor
for i=1:sizeia(1,1)
for j=1:sizeia(1,2)
    changeC0(i,j)=(avgimgarea(i,j,1)./M(1));%LI is divided by the factor of change..
    changeC100(i,j)=(avgimgarea(i,j,npts)./M(8));
end
end   

%% --- Data Modification --- % 
% Perect C=0%
for i=1:npts
   avgimgareaC0(:,:,i)=avgimgarea(:,:,i)./changeC0;
end
% Perect C=100%
for i=1:npts
   avgimgareaC100(:,:,i)=avgimgarea(:,:,i)./changeC100;
end

%% --- Plot Modified Data ---- %
if plotperf==1
 for i = 1:npts
   figure(i)
   imagesc(avgimgareaC0(pixlim:end,:,i));  
   
    axis equal
    axis tight
    xlabel('pixels')
    ylabel('pixels')
    
    caxis([0 255])
    c = colorbar;
    %colormap gray
    colormap(jet(256))
    % Creating Colorbar tittle with proper possition and orienation.
    text('Units','points','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'String','Light Intensity',...
    'Position',[350 30 0]);
    title ('Perfect C=0%')
 
   figure(i+npts)
   imagesc(avgimgareaC100(pixlim:end,:,i));  
   
    axis equal
    axis tight
    xlabel('pixels')
    ylabel('pixels')
    
    caxis([0 255])
    c = colorbar;
    %colormap gray
    colormap(jet(256))
    % Creating Colorbar tittle with proper possition and orienation.
    text('Units','points','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'String','Light Intensity',...
    'Position',[350 30 0]);
    title ('Perfect C=100%')  
end
end


%% --- DMA prediction Data ---
for i=1:npts
  LI1=avgimgareaC0(:,:,i);
  LI1=reshape(LI1,[],1); % with reshape we make LI a single column.
  LI2=avgimgareaC100(:,:,i);
  LI2=reshape(LI2,[],1); % with reshape we make LI a single column.
  DATAC0(:,i)=LI1;
  DATAC100(:,i)=LI2;
  clear LI1 LI2
end
DATAAC0=[DATAC0(:,1);DATAC0(:,2);DATAC0(:,3);DATAC0(:,4);...
    DATAC0(:,5);DATAC0(:,6);DATAC0(:,7);DATAC0(:,8)];
DATAAC100=[DATAC100(:,1);DATAC100(:,2);DATAC100(:,3);DATAC100(:,4);...
    DATAC100(:,5);DATAC100(:,6);DATAC100(:,7);DATAC100(:,8)];
