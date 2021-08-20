
%Ben and Syafiq ML project
%for tryty = 1:10
clear all;
close all;
filename = sprintf('TenRX_rec_x1_y0p5_Nx32_Ny16_k0p4_1_2_%d_T',2);
imp = load(filename);
X0 = (imp.Z0 - imp.Z0')/(2i); %imaginary part of the impedance matrix
Rr = (imp.Z0 + imp.Z0')/(2); %radiation resistance
Rs = 1e-3; %surface resistance
[V,D] = eig(Rr);
D(D<0) = 0;
Rr = V*D/V; %make Rr PSD
Rd = Rr + imp.Zid*Rs; %real part of the impedance matrix
Z = Rd + 1i*X0; %impedance matrix
[RE,index] = max(eig(Rr,Rd));  %max radiation efficiency (I^HRrI)/(I^HRdI)
N = length(Z); %number of elements some are in x and others in y direction
V = zeros(N,1);
V(floor(N/4)-17) = 1; %set the input voltage
V(floor(N/4)+14) = 1;

%start plot currents
[I_temp,eig_val] = eig(Rr,Rd);
I = I_temp(:,index);
opts.im = 1;
opts.re = 1;
opts.geoc = 1;
figure;
%[Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I,1,1,1,opts);
%end plot currents

I_MoM = Z\V; %method of moments currents
RE_MoM = real((I_MoM'*Rr*I_MoM)/(I_MoM'*Rd*I_MoM));

%start define geometry
geo_var = zeros(32,16);



[ind1,index_geo,GeoArea] = Rec_region2bas(imp.bas,imp.meshp,geo_var); %This function returns the index values left over after cutting it is
opts.geoind = geo_var;
figure;
%[Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I_MoM,1,1,1,opts);


%start plot currents
I_old = zeros(length(Rr),1);
[I_temp,eig_val] = eig(Rr(index_geo,index_geo),Rd(index_geo,index_geo));
[RE,index] = max(eig(Rr(index_geo,index_geo),Rd(index_geo,index_geo)));  %max radiation efficiency (I^HRrI)/(I^HRdI)
I = I_temp(:,index);
opts.im = 1;
opts.re = 1;
opts.geoc = 1;
figure;
I_old(index_geo) = I;
%[Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I_old,1,1,1,opts);
%end plot currents

I_MoM = Z\V; %method of moments currents
RE_MoM = real((I_MoM'*Rr*I_MoM)/(I_MoM'*Rd*I_MoM));


%end define geometry

%here the new variables are defined
V_new = V(index_geo);
Z_new = Z(index_geo,index_geo);
I_MoM_new = Z_new\V_new;
Rr_new = Rr(index_geo,index_geo);
Rd_new = Rd(index_geo,index_geo);
RE_MoM_new = real((I_MoM_new'*Rr_new*I_MoM_new)/(I_MoM_new'*Rd_new*I_MoM_new));

%generate random geometry data set
N_x = length(geo_var(:,1));
N_y = length(geo_var(1,:));
ind_feed = find(V == 1);
y_feed = floor(ind_feed/(N_x-1)) + 1;
x_feed = ind_feed - (N_x-1)*(y_feed-1) + 1;




NN = 10000;

reinforced_var = ones(N_x,N_y);



N_classes = 10;

class_geo = zeros(N_classes,32,16);
look_up_class = zeros(N_classes,1);

N_data = zeros(N_classes,1);
N_geo = zeros(N_classes, N_x, N_y);
%N_geo_data = zeros(NN,N_x, N_y);

N_geo_data_vec = zeros(NN,N_x*N_y);
label_sample = zeros(NN,1);

label_count = zeros(N_classes,1);
class_boundary = zeros(N_classes,1);
class_boundary(1) = 0.45;
class_boundary(2) = 0.4;
class_boundary(3) = 0.35;
class_boundary(4) = 0.3;
class_boundary(5) = 0.25;
class_boundary(6) = 0.2;
class_boundary(7) = 0.15;
class_boundary(8) = 0.1;
class_boundary(9) = 0.05;
class_boundary(10) = 0;
% % class_boundary(1) = 4000;
% % class_boundary(2) = 3000;
% % class_boundary(3) = 2500;
% % class_boundary(4) = 2000;
% % class_boundary(5) = 1300;
% % class_boundary(6) = 1000;
% % class_boundary(7) =500;
% % class_boundary(8) = 250;
% % class_boundary(9) = 100;
% % class_boundary(10) = 0;
% class_boundary(5) = 0.15;
% class_boundary(6) = 0.07;
% class_boundary(7) = 0.04;
% class_boundary(8) = 0;

sto_randvar = zeros(NN,32,16);
%for n = 1:10
i = 0;


ii = 0;
while min(label_count) < 1000 || sum(label_count) < 10000
    i = i + 1;
    % rng shuffle
    Rand_var = randi([0 1],N_x,N_y);
    Rand_var(x_feed,y_feed) = 0;
    Rand_var(x_feed-1,y_feed) = 0;
    Rand_var(x_feed,y_feed+1) = 0;
    Rand_var(x_feed-1,y_feed+1) = 0;
    
    Rand_var(:,9:16) = flip(Rand_var(:,1:8),2); %enforce symmetry
    
    %     Rand_var(x_feed,y_feed+1) = randi([0 1],1,1); %break symmetry
    %     Rand_var(x_feed-1,y_feed+1) = randi([0 1],1,1); %break symmetry
    
    [ind1,index_geo,GeoArea] = Rec_region2bas(imp.bas,imp.meshp,Rand_var);
    V_new = V(index_geo);
    Z_new = Z(index_geo,index_geo);
    I_MoM_new = Z_new\V_new;
    Rr_new = Rr(index_geo,index_geo);
    Rd_new = Rd(index_geo,index_geo);
    %Xe_new = Xe(index_geo,index_geo);
    %Xm_new = Xm(index_geo,index_geo);
    RE_MoM_new(i) = real((I_MoM_new'*Rr_new*I_MoM_new)/(I_MoM_new'*Rd_new*I_MoM_new));
   % RE_MoM_new(i) = real(max(I_MoM_new'*Xe_new*I_MoM_new,I_MoM_new'*Xm_new*I_MoM_new)/(I_MoM_new'*Rd_new*I_MoM_new));
    opts.geoind = Rand_var;
    
    
    
    
    for j = 1:N_classes
        if RE_MoM_new(i) >= class_boundary(j)
            if label_count(j) < 1000
                ii = ii + 1;
                N_data(j) = N_data(j) + 1;
                N_geo(j,:,:) = squeeze(N_geo(j,:,:)) + Rand_var;
                N_geo_data_vec(ii,:) = Rand_var(:);
                label_sample(ii) = j-1;
                label_count(j) = label_count(j) + 1;
                sto_randvar(ii,:,:) = Rand_var;
%                if look_up_class(j) == 0
%                    look_up_class(j) = 1;
%                    class_geo(j,:,:) = Rand_var;
%                end
%                 
%                 sto_randvar(ii,:,:) = class_geo(j,:,:);
                
%                 sto_randvar(ii,1:4,1:4) = ones(4,4)*Rand_var(7,9);
%                 sto_randvar(ii,5:8,1:4) = ones(4,4)*Rand_var(8,9);
%                 sto_randvar(ii,9:12,1:4) = ones(4,4)*Rand_var(9,9);
%                 sto_randvar(ii,13:16,1:4) = ones(4,4)*Rand_var(10,9);
%                 sto_randvar(ii,17:20,1:4) = ones(4,4)*Rand_var(11,9);
%                 sto_randvar(ii,21:24,1:4) = ones(4,4)*Rand_var(12,9);
%                 sto_randvar(ii,25:28,1:4) = ones(4,4)*Rand_var(13,9);
%                 sto_randvar(ii,29:32,1:4) = ones(4,4)*Rand_var(14,9);
%                 
%                 
%                 sto_randvar(ii,1:4,5:8) = ones(4,4)*Rand_var(7,10);
%                 sto_randvar(ii,5:8,5:8) = ones(4,4)*Rand_var(8,10);
%                 sto_randvar(ii,9:12,5:8) = ones(4,4)*Rand_var(9,10);
%                 sto_randvar(ii,13:16,5:8) = ones(4,4)*Rand_var(10,10);
%                 sto_randvar(ii,17:20,5:8) = ones(4,4)*Rand_var(11,10);
%                 sto_randvar(ii,21:24,5:8) = ones(4,4)*Rand_var(12,10);
%                 sto_randvar(ii,25:28,5:8) = ones(4,4)*Rand_var(13,10);
%                 sto_randvar(ii,29:32,5:8) = ones(4,4)*Rand_var(14,10);
%                 
%                 sto_randvar(ii,1:4,9:12) = ones(4,4)*Rand_var(7,11);
%                 sto_randvar(ii,5:8,9:12) = ones(4,4)*Rand_var(8,11);
%                 sto_randvar(ii,9:12,9:12) = ones(4,4)*Rand_var(9,11);
%                 sto_randvar(ii,13:16,9:12) = ones(4,4)*Rand_var(10,11);
%                 sto_randvar(ii,17:20,9:12) = ones(4,4)*Rand_var(11,11);
%                 sto_randvar(ii,21:24,9:12) = ones(4,4)*Rand_var(12,11);
%                 sto_randvar(ii,25:28,9:12) = ones(4,4)*Rand_var(13,11);
%                 sto_randvar(ii,29:32,9:12) = ones(4,4)*Rand_var(14,11);
%                 
%                 sto_randvar(ii,1:4,13:16) = ones(4,4)*Rand_var(7,12);
%                 sto_randvar(ii,5:8,13:16) = ones(4,4)*Rand_var(8,12);
%                 sto_randvar(ii,9:12,13:16) = ones(4,4)*Rand_var(9,12);
%                 sto_randvar(ii,13:16,13:16) = ones(4,4)*Rand_var(10,12);
%                 sto_randvar(ii,17:20,13:16) = ones(4,4)*Rand_var(11,12);
%                 sto_randvar(ii,21:24,13:16) = ones(4,4)*Rand_var(12,12);
%                 sto_randvar(ii,25:28,13:16) = ones(4,4)*Rand_var(13,12);
%                 sto_randvar(ii,29:32,13:16) = ones(4,4)*Rand_var(14,12);
                
                
                
                
            end
            break;
        end
    end
    
    
    
    %end
    %reinforced_var = ...;
end

NN = ii;

figure;
plot(sort(RE_MoM_new));

for i = 1:N_classes
    figure;
    contour(squeeze(N_geo(i,:,:)).');
    axis equal
end




Sigma_c = zeros(N_classes, N_x*N_y,N_x*N_y);
mean_data = zeros(N_classes, N_x,N_y);
mean_data_vec = zeros(N_classes, N_x*N_y);
prior_c = zeros(N_classes,1);
N_total = sum(N_data);

for i = 1:N_classes
    mean_data(i,:,:) = N_geo(i,:,:)/N_data(i);
    mean_data_vec(i,:) = mean_data(i,:);
    prior_c(i) = N_data(i)/N_total;
end


% % % % % for i = 1:NN
% % % % %     Sigma_c(label_sample(i),:,:) = squeeze(Sigma_c(label_sample(i),:,:)) + (squeeze(N_geo_data_vec(i,:)) - mean_data_vec(label_sample(i),:))*(squeeze(N_geo_data_vec(i,:)) - mean_data_vec(label_sample(i),:)).';
% % % % % end

% digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
%     'nndatasets','DigitDataset');
% imds = imageDatastore(digitDatasetPath, ...
%     'IncludeSubfolders',true,'LabelSource','foldernames');


%imds.Folders = 'Images';
for i = 1:NN
    %str = compose('Images/fig%d.png', i);
    imwrite(squeeze(sto_randvar(i,:,:)),sprintf('/home/syafiq/Documents/course/ML/ML_LTH/matlab/Images/fig%d.png',i))
    %     imds.Files(i) = sprintf('Images/fig%d.png',i);
    %     imds.Labels(i) = label_sample(i);
end

f = fullfile('/home/syafiq/Documents/course/ML/ML_LTH/matlab/Images');
imds = imageDatastore(f, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');


cats = label_sample(1:NN);

try_cat = categorical(cats);
for i = 1:NN
    imds.Labels(i) = try_cat(i);
end


figure;
perm = randperm(10000,20);
for i = 1:20
    subplot(5,4,i);
    imshow(imds.Files{perm(i)});
end

labelCount = countEachLabel(imds);

img = readimage(imds,1);
size(img)

numTrainFiles = 999;
[imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');

CC = length(unique(label_sample));

layers = [
    imageInputLayer([32 16 1])
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,64,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
%     
    convolution2dLayer(3,128,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    
%     convolution2dLayer(4,128,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     
    
    
    fullyConnectedLayer(CC+1)
    softmaxLayer
    classificationLayer];




options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'GradientThresholdMethod', 'absolute-value',...
    'GradientThreshold',6,...
    'MaxEpochs',4, ...
    'Momentum',0.65,...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');


% options = trainingOptions('sgdm', ...
%     'LearnRateSchedule','piecewise', ...
%     'LearnRateDropFactor',0.2, ...
%     'LearnRateDropPeriod',5, ...
%     'MaxEpochs',20, ...
%     'ValidationData',imdsValidation, ...
%     'ValidationFrequency',30, ...
%     'MiniBatchSize',64, ...
%     'Plots','training-progress')


net = trainNetwork(imdsTrain,layers,options);


recycle('off');
delete('Images\*');

% %classify
% class_c = 1;
% P_c = 0;
% for i = 1:10
%     class_c = 1;
%     for j = 1:N_classes
%         P_c = log(prior_c(j)) - 0.5*log(abs(2*pi*Sigma_c(,:,:)));
%     end
% end
%
%
%end
