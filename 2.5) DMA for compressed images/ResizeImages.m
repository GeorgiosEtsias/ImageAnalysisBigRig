%Resizes laboratory images to face limitations imposed by the lack of RAM
clear
clc
%% Load orignal images
load('SPCalibrationsG')
load('SPTest15G')
load('SPTest20G')
load('SWFreshOnlyG')
load('SWOnlyG')
%% Imresize and saving
for i=1:8
   SPCalibrationsGsmall(:,:,i)= imresize(SPCalibrationsG(:,:,i),0.1);
end
save('SPCalibrationsGsmall','SPCalibrationsGsmall')

SPTest15Gsmall= imresize(SPTest15G,0.1);
save('SPTest15Gsmall','SPTest15Gsmall')
SPTest20Gsmall= imresize(SPTest20G,0.1);
save('SPTest20Gsmall','SPTest20Gsmall')
SWFreshOnlyGsmall= imresize(SWFreshOnlyG,0.1);
save('SWFreshOnlyGsmall','SWFreshOnlyGsmall')
SWOnlyGsmall= imresize(SWOnlyG,0.1);
save('SWOnlyGsmall','SWOnlyGsmall')
