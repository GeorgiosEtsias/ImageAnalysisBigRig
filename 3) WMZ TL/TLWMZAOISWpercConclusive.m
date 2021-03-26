%-------------------------- G Etsias 30/08/2019 --------------------------%
% --- Based on ic2p31 by G.Robison --- %
% Utilized functions : tlwmz7mine & linear1
% =================== Algorithm products ================================ %
% Calculates ALL the SWI variables necessary for data evaluation -------- %
% 1)Calculates TL of the SW fields generated via ANN -------------------- %
% 2)Determines AOI in heterogeneus and homogeneous aquifers alike ------- %
% 3)Calculates WMZ for each layer of a heterogeneous aquifer (NOVEL) ---- %
% 4)WMZ is calculated for the aquifer as a whole both with and withpou the
% 20%-80% limit imposed in G.Robinson's calulations --------------------- %
% 5)Calculates the saltwater% freshwater% mixing% of pixels in each image %
% ============== Necessary pre-liminary actions ========================= %
% 1)Upload investigated aquifer structure and SW concentration field ---- %
% 2)Determine aquifer characteristics a)heterogeneity b)layers c)bead sizes
% ======================================================================= %
clear
clc
%% ---------------------- Loading SW Fields ------------------------------%
load('SWFmod')%loading te specific concentration image
SW=SWFmod;
set=char('SWFmod');
SingleField = 0; %Set SingleField =1 if the loaded SW field file includes only ONE (SUTRA Simulations)
if SingleField == 1
    sizeia = size(SW);
    sizeia(1,3) = 1;
else
sizeia=size(SW);
end
npixels = sizeia(1,1)*sizeia(1,2);

%% Pixelsize - multily the value by *10 for the smaller (1/10) images) 
pixelsizem=(0.000188288458*10);%10 times smaller images

%% ------------------- Aquifer Characteristics -------------------------- %
Heterogenous=0;%If heterogeneous aquifer is investigated set Heterogenous=1
Layers=1;%Number of aquifer layers
Beads=1;%Number of utilized bead sizes
position=800; % Posiiton (in pixels) of the middle layer 
%% ----------------------- Load Bead Structure--------------------------- %
if Heterogenous==1
    load ('Str13257801325')
    Structure=Str13257801325;
    StructureOriginal=Structure;
end
%% ----------------- Checking Imported Variables ------------------------ %
    if Layers<Beads 
    beep
    input !!!incorrectBeads!!!
    elseif Beads<Layers
    repeat=1;%Set repat=1 for repeated layers 
    else 
    repeat=0;
    end 
%% Turning the BDSR into 1,2,3s is essential for the AOI algorithm to work
% NO MANUAL CHANGE REQUIRED
% -----------------LowerLayer=1, MiddleLayer=2, UpperLayer=3--------------%
% 3 Aquifer Categories Investigated
% 1) Homogeneous or Random Heterogeneous
% 2) Heterogeneous with varying bead sizes 
% 3) Heterogeneous with repeating bead sizes 
if Heterogenous==1
% Calculating dominant sizes in upper midle and lower aquifer areas
   upperl=Structure(1,:);
   upperl=mode(upperl);
   lowerl=Structure(sizeia(1),:);
   lowerl=mode(lowerl);
   midlel=Structure(position,:);
   midlel=mode(midlel);
% Assigning the essential values   
 for i=1:sizeia(1)
   for j=1:sizeia(2)
% -------------------- 2 Layered Aqifers -------------------------------- %
       if Layers==2 
           if Structure(i,j)==lowerl
               Structure(i,j)=1;
           else 
               Structure(i,j)=2;
           end
% ------------------- 3 Layered Aquifers -------------------------------- %
       elseif Layers==3 
           if repeat==0 % No Sandwich
           if Structure(i,j)==lowerl
               Structure(i,j)=1;
           elseif Structure(i,j)==midlel
               Structure(i,j)=2;
           else
               Structure(i,j)=3;
           end
           else %Sandwich
             if Structure(i,j)==midlel
               Structure(i,j)=2;
             elseif Structure(i,j)==lowerl && i>position
               Structure(i,j)=1;
             else 
               Structure(i,j)=3;
             end
           end       
       end
   end
 end
else
    Structure=ones(sizeia(1),sizeia(2));%In the absence of any heterogenity
    StructureOriginal=Structure;
end

%% Plotting the Original and Modified Heterogeneous Structure
% -------------------------- Original Bead Structure ----------------------    
figure (1)
    %imagesc(prediction(:,:,i))
    imagesc( StructureOriginal)
    axis equal
    axis tight
    xlabel('X(pixels)')
    ylabel('Y(pixels)')
    colormap(jet(256))
   caxis([780 1325])
   %c = colorbar;
  % legend('780','1090','1325')
  title('Original Aquifer Structure(780um~1325um)')
   text('Units','points','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'String','Bead Size (\mum)',...
    'Position',[352 50 0]);
% ------------------- Modified Bead Structure ----------------------------- 
    figure (2)
    %imagesc(prediction(:,:,i))
    imagesc(Structure)
    axis equal
    axis tight
    xlabel('X(pixels)')
    ylabel('Y(pixels)')
    colormap(jet(256))
   caxis([1 3])
   %c = colorbar;
  % legend('780','1090','1325')
  title('Modified Aquifer Structure(1~3 increasing from bottom to top)')
   text('Units','points','VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'Rotation',90,...
    'String','Bead Size (1~3)',...
    'Position',[352 50 0]);
figure (3)
imagesc (Structure)

%% Dimensions in meters
xmatpre=zeros(sizeia(1,1),sizeia(1,2));
ymatpre=zeros(sizeia(1,1),sizeia(1,2));
 for i = 1:sizeia(1,1)
 for j = 1:sizeia(1,2)
  xmatpre(i,j) = j*pixelsizem;
  ymatpre(i,j) = i*pixelsizem;
 end
end

xmat = flipdim(xmatpre,2); %SW shoulds always be at the left of the image.
ymat = flipud(flipdim(ymatpre,2)); %That is why the contour plots are reversed!!

%% Toe Length and WMZ(for the whole aquifer) calculation  
for i=1:sizeia(1,3) 
cdata1=[25 47.5 75]; %specific concentrations I want to contour plot
img = SW(:,:,i);
edge=1;%remove edge effects when calculating toe length
%value of 1 - yes (experimental), value of 0 = no (numerical) [uses the algorithm for his lab and numerical investigations?]
offset=0;
tl_predicted=0; %extrapolate to predict TL at tank bottom boundary = 1, additional parameters: offset

[tlfull1,wmzfull,npointsfull,c25full,c50full,c75full] = tlwmzFract(xmat,ymat,...
img,cdata1,edge,offset,tl_predicted);

cdata2=[25 52.5 75]; %specific concentrations I want to contour plot
img = SW(:,:,i);
edge=1;%remove edge effects when calculating toe length
%value of 1 - yes (experimental), value of 0 = no (numerical) [uses the algorithm for his lab and numerical investigations?]
offset=0;
tl_predicted=0; %extrapolate to predict TL at tank bottom boundary = 1, additional parameters: offset

[tlfull2,wmzfull,npointsfull,c25full,c50full,c75full] = tlwmzFract(xmat,ymat,...
img,cdata2,edge,offset,tl_predicted);

tlfull=(tlfull1(2)+tlfull2(2))/2;

TL(i)=tlfull;
WMZ(i)=wmzfull;
end
TL=TL';
WMZ=WMZ';


%% Turning dimensions back in pixels
xmatpre=zeros(sizeia(1,1),sizeia(1,2));
ymatpre=zeros(sizeia(1,1),sizeia(1,2));
 for i = 1:sizeia(1,1)
 for j = 1:sizeia(1,2)
  xmatpre(i,j) = j;
  ymatpre(i,j) = i;
 end
 end
  
xmat = flipdim(xmatpre,2); %SW shoulds alway be at the left of the image.
ymat = flipud(flipdim(ymatpre,2)); %That is why the contour plots are reversed!!
mflip=1;

%% AOI and layer specific WMZ Calculation
%Execute tlwmzmine equation
for i=1:sizeia(1,3) %Executed for ALL images
cdata=[25 50 75]; %specific concentrations I want to contour plot
img = SW(:,:,i);
edge=1;%remove edge effects when calculating toe length
%value of 1 - yes (experimental), value of 0 = no (numerical) [uses the algorithm for his lab and numerical investigations?]
offset=0;
tl_predicted=0; %extrapolate to predict TL at tank bottom boundary = 1, additional parameters: offset
[tlfull,wmzfull,npointsfull,c25full,c50full,c75full] = tlwmzFract(xmat,ymat,...
img,cdata,edge,offset,tl_predicted);

%Only interested in the contour plots, not TL & WMZ
                    roundedc25 = round(c25full);
                    roundedc50 = round(c50full);
                    roundedc75 = round(c75full);
 % convert contour coordinates to image space coords
                    contour_imgspace(sizeia(1,1),sizeia(1,2)) = zeros;
                    
                    if sum(roundedc25(1,:)) == 0 % if no 25% contour exists
                        contour_imgspace = contour_imgspace;
                    else
                        for k = 1:length(roundedc25)
                            % rows = y co-ordinates, columns = x co-ordintaes
                            contour_imgspace(roundedc25(2,k),roundedc25(1,k)) = 1;
                        end
                    end
                    
                    if sum(roundedc50(1,:)) == 0 % if no 50% contour exists
                        contour_imgspace = contour_imgspace;
                    else
                        for k = 1:length(roundedc50)
                            % rows = y co-ordinates, columns = x co-ordintaes
                            contour_imgspace(roundedc50(2,k),roundedc50(1,k)) = 2;
                        end
                    end
                    if sum(roundedc75(1,:)) == 0 % if no 75% contour exists
                        contour_imgspace = contour_imgspace;
                    else
                        for k = 1:length(roundedc75)
                            % rows = y co-ordinates, columns = x co-ordintaes
                            contour_imgspace(roundedc75(2,k),roundedc75(1,k)) = 3;
                        end
                    end

 if mflip == 1;
                        contour_imgspace = flipud(flipdim(contour_imgspace,2));
                    else
                        contour_imgspace = flipud(contour_imgspace,2);
 end
                    % figure(i)
                    % imagesc(contour_imgspace(:,:,i))
                    % axis equal
                    % axis tight

  % Repick contour lines in image space and relate to bead
                    % boundary identifier
                    % 25% = 1, 50% = 2, 75% = 3
%At this point the contour lines are separated according to each bead layer
                    for j = 1:length(cdata)% Is 3, more concentrations can be investigated if needed
                        scratch1 = contour_imgspace == j;
                        identif_imgspace(:,:,j) = Structure.*(scratch1);
                        clear scratch1
                        %                     figure
                        %                     imagesc(identif_imgspace(:,:,j))
                        %                     axis equal
                        %                     axis tight
                    end
% identif_imgspace contains the contours with appropriate
                    % bead identifiers for each concentration
                    for kk = 1:max(max(Structure))
                        if sum(sum(identif_imgspace(:,:,1) == kk)) == 0 || ...
                                sum(sum(identif_imgspace(:,:,2) == kk)) == 0 || ...
                                sum(sum(identif_imgspace(:,:,3) == kk)) == 0
                            % if all 3 contours do not exist in zone
                            angle_intrusion(kk,1:4) = nan;
                            presence(kk,1) = nan;
                            wmzp(kk,1) = nan;
                            npointsp(kk,1) = nan;
                        else % if all 3 contours do exist, calculate intrusion angle for each
                            % determine linear gradient of data
                            % c50
                            for jj = 1:3
                                [yloc,xloc] = find(identif_imgspace(:,:,jj) == kk);
                                % gradient = ydiff/xdiff
                                x0 = [1;0]; % Starting guess
                                options = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000);
                                [x,resnorm,residual,exitflag] = lsqcurvefit(@linear1,x0,...
                                    xloc,yloc,[],[],options);
                                % figures for de-bugging
                                %                        figure
                                %                        imagesc(img)
                                %                        hold on
                                %                        plot(xloc,yloc,'kk','linewidth',1)
                                %                        f = x(1)*(min(xloc):max(xloc))+x(2)
                                %                        plot((min(xloc):max(xloc)),f,'kk','linewidth',2)
                                %                        axis equal
                                %                        axis tight
                                %                        hold off
                                % angle is gradient*45 degrees
                                angle_intrusion(kk,jj) = -(x(1)*45);
                                
                                % determine 'presence' of contour in each bead
                                % zone. Presence is a measure of how much of the
                                % contour is in each zone. Low presence could
                                % provide angles and widths that are not
                                % representative of the bead class. Presence is
                                % the percentage ratio of the length of the
                                % contour in the x direction to the total length
                                % of the bead zone and is an arbitrary
                                % indicator.
                                % 50% contour presence representative for all
                                %if jj == 2
                                %    presence(kk,1) = round(((max(xloc)-min(xloc))/...
                                %        sizeia(2))*100);
                                %end
                            end
                        end
                            % average angle of intrusion for 3 contour lines
                            angle_intrusion(kk,j+1) = mean(angle_intrusion(kk,1:j)); 
                            % average angle of all bead layers (of the
                            % average of the 3 concentrations 25,50,75)
                            if kk==max(max(Structure))
                                %kk
                                %angle_intrusion(1,j+1)
                                %angle_intrusion(kk,j+1)
                            interlayerAOI=mean(angle_intrusion(1:kk,j+1));
                            end
          %% -------------------- WMZ calculation for each bead layer ------------- %
                      if sum(sum(identif_imgspace(:,:,1) == kk)) == 0 || ...
                                sum(sum(identif_imgspace(:,:,2) == kk)) == 0 || ...
                                sum(sum(identif_imgspace(:,:,3) == kk)) == 0
                            % if all 3 contours do not exist in zone
                             wmz(kk)=0;
                            
                       else % if all 3 contours do exist, calculate intrusion angle for each
                            % determine linear gradient of data
                                
                            %[yloc25(kk),xloc25(kk)] = find(identif_imgspace(:,:,1) == kk);
                            %[yloc75(kk),xloc75(kk)] = find(identif_imgspace(:,:,3) == kk); 
                            
                            [xx25,yy25] = find(identif_imgspace(:,:,1) == kk);
                            [xx75,yy75] = find(identif_imgspace(:,:,3) == kk);
                            c25lim=[xx25,yy25]; %Part of 25% contour in the current layer
                            c75lim=[xx75,yy75]; %Part of 75% contour in the current layer
                            sizec25=size(c25lim); %Number of pixels in contour
                            countwmz(kk)=0;
                            for hh=1:(sizec25(1))
                            xf(hh).xx = find(c25lim(hh,2)==c75lim(:,2));
                            % finding the pixels with same x in both contours
                            if xf(hh).xx~=0
                               countwmz(kk)=countwmz(kk)+1;
                               samepixels=size(xf(hh).xx);
                               diff=0;
                               for zz=1:samepixels(1)%If more than one pixels have smae x's 
                               diff=(c75lim(xf(hh).xx(zz))-c25lim(hh,1))+diff;
                               end
                               pixelwmz(countwmz(kk))=diff/samepixels(1);
                            end
                            end
                            wmz1=mean(pixelwmz(:));
                            wmz(kk)=wmz1*pixelsizem;
                            clear wmz1
                            clear xf
                            clear pixelwmz
                      end
                    end
                    %% WMZ for all layers without the 20%-80% limitation
                    %5 different cases are possible
                    if Layers==2
                        if wmz(2)~=0
                            wmz(3)=(wmz(1)*countwmz(1)+wmz(2)*countwmz(2))/(countwmz(1)+countwmz(2));
                        else
                            if wmz(1)~=0 
                            wmz(3)=wmz(1);
                            else
                            wmz(3)=0;
                            end
                        end
                    elseif Layers==3
                        if wmz(3)~=0
                            wmz(4)=(wmz(1)*countwmz(1)+wmz(2)*countwmz(2)...
                                +wmz(3)*countwmz(3))/(countwmz(1)+countwmz...
                                (2)+countwmz(3));
                        else
                            if wmz(2)~=0
                            wmz(4)=(wmz(1)*countwmz(1)+wmz(2)*countwmz(2))/(countwmz(1)+countwmz(2));
                            elseif wmz(1)~=0
                                wmz(4)=wmz(1);
                            else
                                wmz(4)=0;
                            end
                        end
                    end
                                                
%% Create Structure that Incorporates all AOIs and WMZs for all investigated images
if Layers==1%-------------------------------------------------------------%
    AOI(i).Figure=i;
    AOI(i).Mean=angle_intrusion(1,4);
    AOI(i).C25=angle_intrusion(1,1);
    AOI(i).C50=angle_intrusion(1,2);
    AOI(i).C75=angle_intrusion(1,3);
    AOI(i).ALL= angle_intrusion(1,:);
    WMZperLayer(i).ALLNolimit=wmz(1);
elseif Layers==2%---------------------------------------------------------%
    AOI(i).Figure=i;
    % Upper Layer
    AOI(i).UpperLayerMean=angle_intrusion(2,4);
    AOI(i).UpperLayer25=angle_intrusion(2,1);
    AOI(i).UpperLayer50=angle_intrusion(2,2);
    AOI(i).UpperLayer75=angle_intrusion(2,3);
    % Lower Layer
    AOI(i).LowerLayerMean=angle_intrusion(1,4);
    AOI(i).LowerLayer25=angle_intrusion(1,1);
    AOI(i).LowerLayer50=angle_intrusion(1,2);
    AOI(i).LowerLayer75=angle_intrusion(1,3);
    % Interlayer
    AOI(i).Interlayer=interlayerAOI;
    % Summary
    AOI(i).UpperLayerALL=angle_intrusion(2,:);
    AOI(i).LowerLayerALL=angle_intrusion(1,:);
    % WMZ
    WMZperLayer(i).UpperLayerALL=wmz(2);
    WMZperLayer(i).LowerLayerALL=wmz(1);
    WMZperLayer(i).ALLNolimit=wmz(3);
elseif Layers==3%---------------------------------------------------------%
    AOI(i).Figure=i;
    % Upper Layer
    AOI(i).UpperLayerMean=angle_intrusion(3,4);
    AOI(i).UpperLayer25=angle_intrusion(3,1);
    AOI(i).UpperLayer50=angle_intrusion(3,2);
    AOI(i).UpperLayer75=angle_intrusion(3,3);
    % Middle Layer
    AOI(i).MiddleLayerMean=angle_intrusion(2,4);
    AOI(i).MiddleLayer25=angle_intrusion(2,1);
    AOI(i).MiddleLayer50=angle_intrusion(2,2);
    AOI(i).MiddleLayer75=angle_intrusion(2,3);
    % Lower Layer
    AOI(i).LowerLayerMean=angle_intrusion(1,4);
    AOI(i).LowerLayer25=angle_intrusion(1,1);
    AOI(i).LowerLayer50=angle_intrusion(1,2);
    AOI(i).LowerLayer75=angle_intrusion(1,3);
    % Interlayer
    AOI(i).Interlayer=interlayerAOI;
    % Summary
    AOI(i).UpperLayerALL=angle_intrusion(3,:);
    AOI(i).MiddleLayerALL=angle_intrusion(2,:);
    AOI(i).LowerLayerALL=angle_intrusion(1,:);
    % WMZ
    WMZperLayer(i).UpperLayer=wmz(3);
    WMZperLayer(i).MiddleLayer=wmz(2);
    WMZperLayer(i).LowerLayer=wmz(1);
    WMZperLayer(i).ALLNolimit=wmz(4);
end
%clearvars -except xmat ymat Structure sizeia SW set npixels mflip Layers AOI
clear contour_imgspace % If not cleared causes problems to subsequent executions
end

%% Saltwater - Mixng - Freshwater Percentage
% Percentages = Salt (>75%), Mixi(5%-75%), Fresh(<5%)
Perc=zeros(sizeia(1,3),3);

for k=1:sizeia(1,3)
    salt=0;
    mix=0;
    fresh=0;
    for i=1:sizeia(1)
        for j=1:sizeia(2)
            if SW(i,j,k)<5
                fresh=fresh+1;
            elseif SW(i,j,k)>=75
                salt=salt+1;
            else
                mix=mix+1;
            end
        end
    end
    Perc(k,1)=(salt/(sizeia(1)*sizeia(2)))*100;
    Perc(k,2)=(mix/(sizeia(1)*sizeia(2)))*100;
    Perc(k,3)=(fresh/(sizeia(1)*sizeia(2)))*100;
end

clearvars -except  TL WMZ WMZperLayer AOI Perc set

name1 = strcat('TL',set);
name2 = strcat('WMZ',set);
name3 = strcat('WMZperLayer',set)
name4 = strcat('AOI',set);
name5 = strcat('Perc',set);

beep