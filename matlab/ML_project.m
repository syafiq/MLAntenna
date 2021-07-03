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
%[Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I,1,1,1,opts);
%end plot currents

I_MoM = Z\V; %method of moments currents
RE_MoM = real((I_MoM'*Rr*I_MoM)/(I_MoM'*Rd*I_MoM));

%start define geometry
geo_var = zeros(32,16);
geo_var(11,8) = 1; %by setting entry equal to one geometry is cut should try to not cut the voltage source
[ind1,index_geo,GeoArea] = Rec_region2bas(imp.bas,imp.meshp,geo_var); %This function returns the index values left over after cutting it is
opts.geoind = geo_var;
figure;
%[Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I_MoM,1,1,1,opts);
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
N1 = 0; %number of samples in class
N2 = 0;
N3 = 0;
N4 = 0;
N5 = 0;
N6 = 0;
N7 = 0;

N1_geo = zeros(N_x,N_y);
N2_geo = zeros(N_x,N_y);
N3_geo = zeros(N_x,N_y);
N4_geo = zeros(N_x,N_y);
N5_geo = zeros(N_x,N_y);
N6_geo = zeros(N_x,N_y);
N7_geo = zeros(N_x,N_y);

reinforced_var = ones(N_x,N_y);
%for n = 1:10
for i = 1:10000
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
    if RE_MoM_new(i) > 0.6 %class 1
        N1_geo = N1_geo + Rand_var;
        N1 = N1 + 1;
        figure;
        [Jx,Jy,xx,yy,rho]=Jplot4(imp.bas,imp.meshp,I_MoM,1,1,1,opts);
    elseif RE_MoM_new(i) > 0.5 %class 2
        N2_geo = N2_geo + Rand_var;
        N2 = N2 + 1;
    elseif RE_MoM_new(i) > 0.4 %class 3
        N3_geo = N3_geo + Rand_var;
        N3 = N3 + 1;
    elseif RE_MoM_new(i) > 0.3 %class 4
        N4_geo = N4_geo + Rand_var;
        N4 = N4 + 1;
    elseif RE_MoM_new(i) > 0.2 %class 5
        N5_geo = N5_geo + Rand_var;
        N5 = N5 + 1;
    elseif RE_MoM_new(i) > 0.1 %class 6
        N6_geo = N6_geo + Rand_var;
        N6 = N6 + 1;
    else   %class 7
        N7_geo = N7_geo + Rand_var;
        N7 = N7 + 1;
    end
%end
%reinforced_var = ...;
end

%figure;
%plot(sort(RE_MoM_new));
%figure;
%contour(N1_geo.');
%axis equal
%figure;
%contour(N2_geo.');
%axis equal
%figure;
%contour(N3_geo.');
%axis equal
%figure;
%contour(N4_geo.');
%axis equal
%figure;
%contour(N5_geo.');
%axis equal
%figure;
%contour(N6_geo.');
%axis equal
%figure;
%contour(N7_geo.');
%axis equal
%figure;
%contour(N1_geo.'+N2_geo.'+N3_geo.'+N4_geo.'+N5_geo.'+N6_geo.'+N7_geo.');
%axis equal