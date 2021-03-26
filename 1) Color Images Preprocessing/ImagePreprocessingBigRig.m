% ---------------G.Etsias November 10th 2020----------------------------- %
% =======================Purpose========================================= %
% Image Bounding and Spatial Synchronisation of colored images----------- %
% =======================Notes=========================================== %
% Designed to work with big rig images ---------------------------------- %
% All images belong to a single folder, 'ColorImages'-------------------- %
% No need for file name alterations-------------------------------------- %
% Functions used: spatoriginBigRig,spatsync2, spatsyncColor--------------- %
% Pixel size derived by the previously calulated value------------------- %
% ======================================================================= %
% ---------- Based on ColorImagePreprocessing2 (9-JAN-2019) ------------- %
% ======================================================================= %
% Final results include RGB, R,G,B and Greyscale image/matrices---------- %
tic
%% --- Variables ---
clc
clear
% === Variables to set manually  ==========================================
npts =3; % No. of printed  timesteps( equals to the number of photos )
fps = '/'; % File path separator, '\' for windows, '/' for Unix 
%
% --- Pixel sizing --------------------------------------------------------
scaledistmm = 10; % Distance on rule used to calculate pixel size
ifpixel = 0; % Determine pixel size? If yes, set to 1
psts = 1; % Pixel sizing time step (image of ruler)

% --- Domain bounding -----------------------------------------------------
ifdomain = 1; % Determine IA domain bounds? If yes, set to 1
dbts = 1; % Domain bounding time step (usually 1)
spatialsync = 1; % Spatial Synchronisation required? If yes, set to 1
dl = 1; % Domain length(m)-big rig
dh = 0.46; % Domain height(m)minus the unsaturated zone(0.54-0.08)
fullwedge = 0; % Set to 1 to analyse full wedge inc. tank sides            
offset = 60; % pixels (similar to the old rig)

%% ---User Defined Options---- %
medfilt = 15; % Median filtering level of 3,5,7,9 etc one greyscale images
medfiltRGB=15;  %Median filtering level of 3,5,7,9 etc one greyscale images
%If no set medfilt=0 or medfiltRGB=0
medfiltColor=0; %Set to 1 to have colored images with median filtering(not avdised)
plotim=0;%Set to 1 to imshow the processed colored images 
% ================(case specific saving)================================= %
imagedistribution=0; % Set to 1 to distribute all the images according to:
%commands in the final section of the script

%% --- Read Data ---
path = pwd; %(definition) pwd = Identify current folder
cd(path); %(definition) cd = Change current folder
    
    % Defining folder name, 'fname'    
    fname = 'ColorImages';  
    % cd to directory
    cd(strcat(path, fps, fname));
    % Read in image names
    fimgnames = dir('*.jpg');%dir name lists files and folders that match name.
   % cd back to original folder
    cd(path)

for i = 1:npts % number of aquired images 
    % cd to directory----------------------------------------------------
    cd(strcat(path, fps, fname));
   
    % Read the Image
      avgimgcolor(:,:,:,i)= imread(fimgnames(i,1).name);% Read colored image
      A = rgb2gray(avgimgcolor(:,:,:,i)); % Turn RGB image into Grayscale 
      avgimg(:,:,i) = round(mean(A,3));
      %Convert greyscale to double, essential for 'spatoriginColor'
        
    % cd back to original folder----------------------------------------
    cd(path)
    
  
   %% --- Determine Pixel Size ---
    % Function 'pixelsize.m'
    if ifpixel == 1 && i == psts
        pixelsizem = pixelsize(scaledistmm,(avgimg(:,:,psts)));
    elseif i == psts
       pixelsizem = 0.000188288458;%Pixelsize calculated manually through image analysis
    end 
    
    %% --- Spatial synchronisation co-ordinate finder ---
    % Select same point in each time step image and shift so that the
    % images are exactly overlapping each other
    % Function 'pixpoint.m'
    % Pick co-ordinate at the bottom of the SW side, therefore can use this
    % later to determine how far the bounded domain is away from the side
    if spatialsync == 1
            % only need the br co-ordinate to specify origin and syncpoints
           
                [tl,tr,bl,br] = spatoriginBigRig(avgimg(:,:,i),pixelsizem,dh);
            xss(:,1,i) = br(1,1);
            yss(:,1,i) = br(2,1);
            
    end
    
    %% --- Bound IA domain ---
    % Bottom left corner coords
    % Requires function 'pixpoint.m'
    if ifdomain == 1 && i == dbts
        
        [tl,tr,bl,br] = spatoriginBigRig(avgimg(:,:,i),pixelsizem,dh);    
         
        
        % Collect all the x and y terms for the bounding box corners
        if fullwedge == 1 % includes sides of tank
            % can drastically affect wmz calculation as overcomplicates
            % contour map
            xbounds = [bl(1,1) tl(1,1) tr(1,1) br(1,1)];
            ybounds = [bl(2,1) tl(2,1) tr(2,1) br(2,1)];
        else % only bottom and x co-ords need offset
            xbounds = [bl(1,1)+offset tl(1,1)+offset tr(1,1)-offset ...
                br(1,1)-offset ];
            ybounds = [bl(2,1)-offset/2 tl(2,1) tr(2,1)...
                br(2,1)+offset/2 ];                 
        end
    end
    
    %% ---Median Filtering----
    % Find pixel values bounded by the determined IA domain. Apply median
    % filtering to reduce LI noise
    % Fixing pixel positions specifically for the calibration relationship
    % finalresults folder
    
    if medfilt > 0
%         avgimgareafull(:,:,i) = medfilt2(avgimg(ybounds(1,2):ybounds(1,1),...
%             xbounds(1,2):xbounds(1,3),i),[medfilt medfilt]);
        avgimgareafull(:,:,i) = medfilt2(avgimg(:,:,i),[medfilt medfilt]); %Median filt: Greyscale
        
        %Notes on RGB images:-------------------------------------------- %
        %(1)applying median filtering in an RGB image significantly 
        %increases execution time.
        %(2)medfilt3(image,[medfilt medfilt medfilt) does filetering in 3
        %dimensions and returns a greyscale image
        %(3)medfilt3(image) without specification about the m-by-n-by-p 
        %neighborhood returns a modified RGB image.
        if medfiltColor==1
        avgimgcolorfilt(:,:,:,i)= medfilt3(avgimgcolor(:,:,:,i));%Not used
        end
        %---------------------------------------------------------------- %
                 
    elseif medfilt == 0
%         avgimgareafull(:,:,i) = avgimg(ybounds(1,2):ybounds(1,1),...
%             xbounds(1,2):xbounds(1,3),i);
        avgimgareafull(:,:,i) = avgimg(:,:,i);
    end 
end

%% --- Spatial Synchronisation ---
% The results of this section super-impose images, so they have larger
% dimensions than the original images
% Requires function 'spatsync.m' 
if spatialsync == 1
    [ssavgimgarea] = spatsync2(avgimgareafull,xss,yss);
    [avgimgcolor] = spatsyncColor(avgimgcolor,xss,yss);
    if medfiltColor==1
    [avgimgcolorfilt] = spatsyncColor(avgimgcolorfilt,xss,yss);
    end
end

%% --- Trim Image Edges ---
% -- Bound domain --
for i = 1:npts
avgimgarea(:,:,i) = ssavgimgarea(ybounds(1,2):ybounds(1,1),...
            xbounds(1,1):xbounds(1,4),i);
end
%% ----- Creating a final Array Structure ------
% Including the processed RGB, R,G,B median-fil.RGB  and Greyscale images 
for i=1:npts
avgimgNew(i).color=avgimgcolor(ybounds(1,2):ybounds(1,1),...
           xbounds(1,1):xbounds(1,4),:,i); % RGB image (3-D matrix)
avgimgNew(i).R=avgimgcolor(ybounds(1,2):ybounds(1,1),...
           xbounds(1,1):xbounds(1,4),1,i); % R image (2-D matrix)
avgimgNew(i).G=avgimgcolor(ybounds(1,2):ybounds(1,1),...
           xbounds(1,1):xbounds(1,4),2,i); % G image (2-D matrix)
avgimgNew(i).B=avgimgcolor(ybounds(1,2):ybounds(1,1),...
           xbounds(1,1):xbounds(1,4),3,i); % B image (2-D matrix)
       if medfiltColor==1
avgimgNew(i).filt=avgimgcolorfilt(ybounds(1,2):ybounds(1,1),...
           xbounds(1,1):xbounds(1,4),:,i); %(filtered RGB image 3-D matrix)
       end
%avgimgNew(i).greyscale=avgimgarea(:,:,i); % Greyscale image (2-D matrix)
avgimgNew(i).greyscale=avgimg(ybounds(1,2):ybounds(1,1),...
           xbounds(1,1):xbounds(1,4),i); % Greyscale image (2-D matrix)
end

%% ----------- Apply median filtering to 2-D R, G, B images -------------- %
if medfiltRGB >=1
%Median Filtering options:
% i)medfilt2(A,[m n]): Pixels within one-half the width of the neighborhood 
%([m n]/2) of the edges can appear distorted(3 pixels at every corner = 0)
% ii)B = J = medfilt2(___,padopt)
%avgimgNew(i).R = medfilt2(avgimgNew(i).R,[medfiltRGB medfiltRGB], 'symmetric');
%Symmetrically extend the image at the boundaries, thus giving different
%values

for i=1:npts
avgimgNew(i).R = medfilt2(avgimgNew(i).R,[medfiltRGB medfiltRGB], 'symmetric');
avgimgNew(i).G = medfilt2(avgimgNew(i).G,[medfiltRGB medfiltRGB], 'symmetric');
avgimgNew(i).B = medfilt2(avgimgNew(i).B,[medfiltRGB medfiltRGB], 'symmetric');
avgimgNew(i).greyscale = medfilt2(avgimgNew(i).greyscale,[medfiltRGB medfiltRGB], 'symmetric');
end
 
end
%% ----- Turning R, G, B from unit8 into diuble ------- %
for i=1:npts
avgimgNew(i).R =  im2double(avgimgNew(i).R)*255;
avgimgNew(i).G =  im2double(avgimgNew(i).G)*255;
avgimgNew(i).B = im2double(avgimgNew(i).B)*255;
end

%% Plotting processed images - not necessary
if plotim==1
%Plotting the RGB images
for i=1:npts
    figure(i)
    imshow(avgimgNew(i).color)
end
%Plotting the RGB images with median filtering
if medfiltColor==1
 for i=1:npts
    figure(i+4)
    imshow(avgimgNew(i).filt)
 end   
end 
% Plotting the R, G, B images with pseudo-colouring
for i=1:npts
    figure(i+8)
    imagesc(avgimgNew(i).R); 
    axis equal
    axis tight
    xlabel('pixels')
    ylabel('pixels')
    caxis([0 255])
    c = colorbar;
    colormap autumn
    %title('R')
end
for i=1:npts
    figure(i+12)
    imagesc(avgimgNew(i).G); 
    axis equal
    axis tight
    xlabel('pixels')
    ylabel('pixels')
    caxis([0 255])
    c = colorbar;
    colormap summer
    %title('G')
end
for i=1:npts
    figure(i+16)
    imagesc(avgimgNew(i).B); 
    axis equal
    axis tight
    xlabel('pixels')
    ylabel('pixels')
    caxis([0 255])
    c = colorbar;
    colormap winter
    %title('B')
end
for i=1:npts
    figure(i+20)
    imagesc(avgimgNew(i).greyscale); 
    axis equal 
    axis tight
    xlabel('pixels')
    ylabel('pixels')
    caxis([0 255])
    c = colorbar;
    %colormap gray
    colormap jet(256)
    %title('Greyscale')
end
end

%% Case-specific image distribution - not necessary %
% ---Connected in specific to the january 2019 experimental procedure---- %
if imagedistribution==1
%1)3 Homogeneous Aquifers 1090-780-1325um freshwater only
onlyfresh1=avgimgNew(1:3);
%2)3 Homogeneous Aquifers 1090-780-1325um freshwater only
onlyfresh2=avgimgNew(4:6);
%3)3 Homogeneous Aquifers 1090-780-1325um freshwater only
onlyfresh3=avgimgNew(7:9);
%4)3 Homogeneous Aquifers 1090-780-1325um freshwater only
onlyfresh4=avgimgNew(10:12);
% ----------------------------------------------------------------------- %
%5)1090um 8 Calibration Images
CalHomo1090=avgimgNew(13:21);
%6)1090um f-6mm Test Images
Homo1090Test1=avgimgNew(13:21);
%7)1090um 6-4mm Test Images
Homo1090Test2=avgimgNew(13:21);
%8)1090um 4-5mm Test Images
Homo1090Test3=avgimgNew(13:21);
% ----------------------------------------------------------------------- %
%9)780um 8 Calibration Images
CalHomo780=avgimgNew(13:21);
%10)1090um f-6mm Test Images
%11)1090um 6-4mm Test Images
%12)1090um 4-5mm Test Images
% ----------------------------------------------------------------------- %
%13)1325um 8 Calibration Images
CalHomo1325=avgimgNew(13:21);
%14)1090um f-6mm Test Images
%15)1090um 6-4mm Test Images
%16)1090um 4-5mm Test Images
% ----------------------------------------------------------------------- %
%17)1325-1090-780um 8 Calibration Images
CalLayerOne=avgimgNew(13:21);
%18)325-1090-780um f-6mm Test Images
%19)325-1090-780um 6-4mm Test Images
%20)325-1090-780um 4-5mm Test Images
% ----------------------------------------------------------------------- %
%21)780-1090-1325um 8 Calibration Images
CalLayerTwo=avgimgNew(13:21);
%22)780-1090-1325um f-6mm Test Images
%23)780-1090-1325um 6-4mm Test Images
%24)780-1090-1325um 4-5mm Test Images
% ----------------------------------------------------------------------- %
%25)1325-780-1325um 8 Calibration Images
CalLayerThree=avgimgNew(13:21);
%26)1325-780-1325um f-6mm Test Images
%27)1325-780-1325um 6-4mm Test Images
%28)1325-780-1325um 4-5mm Test Images
% ----------------------------------------------------------------------- %
  
end   
toc    
beep
