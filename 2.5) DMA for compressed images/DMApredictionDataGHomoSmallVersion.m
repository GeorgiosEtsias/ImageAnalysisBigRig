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
load('SPCalibrationsGsmall')


% Load mean values of LI
load('M') %from the training data DMA

% Load Homogenisation Factor
load ('MeanHomoFactor') %from the training data DMA

test=0;

% === Variables to manually set ===========================================
npts = 8; % No. of calibration images
% Number of pixels that are not full of beads and are to be
% ignored. (If sandbox full of glass pixlim=1, also pixlim=1 for intial investigations!!)
pixlim=1; % foul aquifers
% If homogenization procedur occures homog=1 if not homog=0
homog=1;
% Set plot perf =1 to plot the perfect C=0% and C=100% aquifers
plotperf=1 ;

% ======================================================================= %
sizeia = size(SPCalibrationsGsmall);
npixels = sizeia(1,1)*sizeia(1,2);

%Derive the LI data for C= 0% - 100% and bead sizes 780um-1090um-1325um
avgimgarea= zeros(sizeia(1),sizeia(2),npts); 
avgimgarea(:,:,1:8)= SPCalibrationsGsmall(:,:,1:end);


%% Image Homogenization
for i=1:npts
        avgimgarea(:,:,i)=avgimgarea(:,:,i)./MeanHomoFactor;
end

%% --- Cange Factor Caluclation ---- %

change1C0=zeros(size(1,1),size(1,2));
change1C100=zeros(size(1,1),size(1,2));




%Perfect C=0Applies only for C=0% images, calculating the change factor
for i=1:sizeia(1,1)
        for j=1:sizeia(1,2)
                change1C0(i,j)=(avgimgarea(i,j,1)./M(1));%LI is divided by the factor of change..
                change1C100(i,j)=(avgimgarea(i,j,8)./M(8));
        end
end
   

%% --- Data Modification --- % 
% Perect C=0%
for i=1:8
   avgimgareaC0(:,:,i)=avgimgarea(:,:,i)./change1C0;
end
% Perect C=100%
for i=1:8
   avgimgareaC100(:,:,i)=avgimgarea(:,:,i)./change1C100;
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
 
   figure(i+8)
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
%===============================780um=====================================%
% X's correct column
x=[1:sizeia(1,2)];
X=repmat(x,sizeia(1,1)-pixlim+1,1);
X=reshape(X,[],1);
% Y's correct column
y=[1:sizeia(1,1)-pixlim+1];
y=y';
Y=repmat(y,sizeia(1,2),1);

%Entering the number of rows in ONE column dataset
nmodpixels= (sizeia(1,1)-pixlim+1)*sizeia(1,2);
%If pixlim=1, nmodpixels==npixels!!
% Bead Size = 780um
BD1=[780];
%BD1=[estbsize(:,1:1194)];
BD1=repmat(BD1,nmodpixels,1);
%BD1=reshape(BD1,[],1);

%-----------C=0%---------------------------
%%Light Intensities column
A=avgimgareaC0(:,:,1);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,1);
% with reshape we make A a single column.
LI100=reshape(A,[],1);

% SW concentration column 
SW=[0];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA0C0=[X,Y,BD1,LI0,SW];
DATA0C100=[X,Y,BD1,LI100,SW];

%------------C=5%---------------------
%%Light Intensities column
A=avgimgareaC0(:,:,2);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,2);
% with reshape we make A a single column.
LI100=reshape(A,[],1);
%SW concentration column 
SW=[5];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA5C0=[X,Y,BD1,LI0,SW];
DATA5C100=[X,Y,BD1,LI100,SW];
%
%------------C=10%---------------------
%%Light Intensities column
A=avgimgareaC0(:,:,3);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,3);
% with reshape we make A a single column.
LI100=reshape(A,[],1);
% SW concentration column 
SW=[10];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA10C0=[X,Y,BD1,LI0,SW];
DATA10C100=[X,Y,BD1,LI100,SW];

%------------C=20%---------------------
%%Light Intensities column
A=avgimgareaC0(:,:,4);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,4);
% with reshape we make A a single column.
LI100=reshape(A,[],1);
% SW concentration column 
SW=[20];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA20C0=[X,Y,BD1,LI0,SW];
DATA20C100=[X,Y,BD1,LI100,SW];

%------------C=30%---------------------
%%Light Intensities column
A=avgimgareaC0(:,:,5);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,5);
% with reshape we make A a single column.
LI100=reshape(A,[],1);
% SW concentration column 
SW=[30];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA30C0=[X,Y,BD1,LI0,SW];
DATA30C100=[X,Y,BD1,LI100,SW];

%------------C=50%---------------------
%%Light Intensities column
A=avgimgareaC0(:,:,6);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,6);
% with reshape we make A a single column.
LI100=reshape(A,[],1);
% SW concentration column 
SW=[50];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA50C0=[X,Y,BD1,LI0,SW];
DATA50C100=[X,Y,BD1,LI100,SW];

%------------C=70%---------------------
%%Light Intensities column
A=avgimgareaC0(:,:,7);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,7);
% with reshape we make A a single column.
LI100=reshape(A,[],1);
% SW concentration column 
SW=[70];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA70C0=[X,Y,BD1,LI0,SW];
DATA70C100=[X,Y,BD1,LI100,SW];

%------------C=100%---------------------
%%Light Intensities column
A=avgimgareaC0(:,:,8);
% with reshape we make A a single column.
LI0=reshape(A,[],1);

A=avgimgareaC100(:,:,8);
% with reshape we make A a single column.
LI100=reshape(A,[],1);
% SW concentration column 
SW=[100];
SW=repmat(SW,nmodpixels,1);
% Combining all 4 variable together 
DATA100C0=[X,Y,BD1,LI0,SW];
DATA100C100=[X,Y,BD1,LI100,SW];

%------Creating the final (sizeia(1)*sizeia(2)*sizeia(3)) x 1 array for regresion analysis 780
DATAAC0=[DATA0C0;DATA5C0;DATA10C0;DATA20C0;DATA30C0;DATA50C0;DATA70C0;DATA100C0];
DATAAC100=[DATA0C100;DATA5C100;DATA10C100;DATA20C100;DATA30C100;DATA50C100;DATA70C100;DATA100C100];

DATAAC0 = DATAAC0(:,4);
DATAAC100 = DATAAC100(:,4);
