%--------------------------July 31st 2018---------------------------------%
%==============Georgios Etsias==============%
% Using Neural Net Fitting App one cannot use parallel programing toolbox.
% Using Big Data demands Parallel Computing to optimise efficiency. Current
% Script Creates a feedforward Neural Network with 10 Neurons and on hidden
% layer.
%Written using autogenerated script adding the parallel programing command
clc
trainn=DATAA(:,1);
goall=DATAA(:,2);
trainn=trainn';
goall=goall';
%trainn=tall(trainn);
%goall=tall(goall);
%If network has only 1 hidden layer with X number of neurons
%feedforwardnet(X)
%If network has 2 hidden layers  with X&Y number of neurons
%feedforwardnet([X Y])
net1 = feedforwardnet([10 10 10]);
net2 = train(net1,trainn,goall,'useParallel','no','showResources','yes');

