%Ben and Syafiq ML project

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

%start plot currents
[I_temp,eig_val] = eig(Rr,Rd);
I = I_temp(:,index);
opts.im = 1;
opts.re = 1;
opts.geoc = 1;
figure;
[Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I,1,1,1,opts);
%end plot currents

I_MoM = Z\V; %method of moments currents
RE_MoM = real((I_MoM'*Rr*I_MoM)/(I_MoM'*Rd*I_MoM));

%start define geometry
geo_var = zeros(32,16);
geo_var(11,8) = 1; %by setting entry equal to one geometry is cut should try to not cut the voltage source
[ind1,index_geo,GeoArea] = Rec_region2bas(imp.bas,imp.meshp,geo_var); %This function returns the index values left over after cutting it is
opts.geoind = geo_var;
figure;
[Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I_MoM,1,1,1,opts);
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



N_classes = 7;
N_data = zeros(N_classes,1);
N_geo = zeros(N_classes, N_x, N_y);
%N_geo_data = zeros(NN,N_x, N_y);

N_geo_data_vec = zeros(NN,N_x*N_y);
label_sample = zeros(NN,1);
class_boundary = zeros(N_classes,1);
class_boundary(1) = 0.6;
class_boundary(2) = 0.5;
class_boundary(3) = 0.4;
class_boundary(4) = 0.3;
class_boundary(5) = 0.2;
class_boundary(6) = 0.1;
class_boundary(7) = 0;


%for n = 1:10
for i = 1:NN
    % rng shuffle
    Rand_var = randi([0 1],N_x,N_y);
    Rand_var(x_feed,y_feed) = 0;
    Rand_var(x_feed-1,y_feed) = 0;
    [ind1,index_geo,GeoArea] = Rec_region2bas(imp.bas,imp.meshp,Rand_var);
    V_new = V(index_geo);
    Z_new = Z(index_geo,index_geo);
    I_MoM_new = Z_new\V_new;
    Rr_new = Rr(index_geo,index_geo);
    Rd_new = Rd(index_geo,index_geo);
    RE_MoM_new(i) = real((I_MoM_new'*Rr_new*I_MoM_new)/(I_MoM_new'*Rd_new*I_MoM_new));
    opts.geoind = Rand_var;
    
    for j = 1:N_classes
        if RE_MoM_new(i) >= class_boundary(j)
            N_data(j) = N_data(j) + 1;
            N_geo(j,:,:) = squeeze(N_geo(j,:,:)) + Rand_var;
            N_geo_data_vec(i,:) = Rand_var(:);
            label_sample(i) = j;
            break;
        end
    end
    
    
    
    %end
    %reinforced_var = ...;
end

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


for i = 1:NN
    Sigma_c(label_sample(i),:,:) = squeeze(Sigma_c(label_sample(i),:,:)) + (squeeze(N_geo_data_vec(i,:)) - mean_data_vec(label_sample(i),:))*(squeeze(N_geo_data_vec(i,:)) - mean_data_vec(label_sample(i),:)).';
end



%classify 
class_c = 1;
P_c = 0;
for i = 1:10
    class_c = 1;
    for j = 1:N_classes
        P_c = log(prior_c(j)) - 0.5*log(abs(2*pi*Sigma_c(,:,:)));
    end
end 



