tic
%----------------- July 27th 2018-----------------------%
%++Generate Prediction from Swallow Neural network trained for ALL
%the calibration images TOGETHER+++
% The variable names are random, only one that HAS to remain the same is
% ShalowNN and DATTA generated from diferent codes respectively.
%1) In case you want to get a prediction for an intrusion image change one or
%multiple image files, better to just cjange file 8.
% 2) If all the data you input are the training data, "perform" can give
% you the performance of the NN.
clc
load('ANNbig')
par=1; %set par=1 to calculate the pareto solution
%% Perfect C=0%
inputs=DATAAC0;

%targets=DATAA(:,4);
inputs=inputs';
%targets=targets';
outputs = ANNbig(inputs);
%errors = gsubtract(targets, outputs);
%performance = perform(ShalowNN, targets, outputs)
OMG=outputs;

garea=reshape(OMG,sizeia(1,1)-pixlim+1,sizeia(1,2),npts);

for k=1:npts
for i=1:sizeia(1)
    for j=1:sizeia(2)
        if garea(i,j,k)>=100
            garea(i,j,k)=100;
        end
        if garea(i,j,k)<=0
            garea(i,j,k)=0;
        end
    end
end
end

PredictionC0=garea;
%% Perfect C=100%
inputs=DATAAC100;

%targets=DATAA(:,4);
inputs=inputs';
%targets=targets';
outputs = ANNbig(inputs);
%errors = gsubtract(targets, outputs);
%performance = perform(ShalowNN, targets, outputs)
OMG=outputs;
garea=reshape(OMG,sizeia(1,1)-pixlim+1,sizeia(1,2),npts);

for k=1:npts
for i=1:sizeia(1)
    for j=1:sizeia(2)
        if garea(i,j,k)>=100
            garea(i,j,k)=100;
        end
        if garea(i,j,k)<=0
            garea(i,j,k)=0;
        end
    end
end
end

PredictionC100=garea;

%% Give the correct dimensions of Predicted SW fields
if test == 1
 PredictionC0=PredictionC0(:,:,2:sizeia(3)+1);
 PredictionC100=PredictionC100(:,:,2:sizeia(3)+1);
end
%% Pareto

sizeia=size(PredictionC0);
pixelsizem=(0.000188288458);%10 times smaller images
ALL=cat(4,PredictionC0,PredictionC100);
PredictionParetoLinAv2055=mean(ALL,4);
for k=1:sizeia(3)
    for i=1:sizeia(1)
       for j=1:sizeia(2)
            if PredictionC100(i,j,k)>=95
                PredictionParetoLinAv2055(i,j,k)=PredictionC100(i,j,k);
            elseif PredictionC0(i,j,k)<=5
                PredictionParetoLinAv2055(i,j,k)=PredictionC0(i,j,k);
            elseif PredictionC100(i,j,k)>=50
                   PredictionParetoLinAv2055(i,j,k)=((PredictionC100(i,j,k)/100)*PredictionC100(i,j,k))+...
                       ((1-(PredictionC100(i,j,k)/100))*PredictionC0(i,j,k));
            elseif PredictionC0(i,j,k)<=20
                PredictionParetoLinAv2055(i,j,k)= (((1-(40-PredictionC0(i,j,k))/40))*PredictionC100(i,j,k))+...
                       ((40-PredictionC0(i,j,k))/40*PredictionC0(i,j,k));
            end
        end
    end
end
for i=1:sizeia(3)
figure (i)
%subplot(2,1,i)
imagesc([0 sizeia(1,2)]*pixelsizem,[0 sizeia(1,1)]*pixelsizem,flipud(PredictionParetoLinAv2055(:,:,i)))
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
toc
SWF=PredictionParetoLinAv2055;