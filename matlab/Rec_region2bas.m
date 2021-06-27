function [ind1,ind2,A,ind3] = Rec_region2bas(bas,mesh,geoind)
%
% uses
% bas: basis parameters
% mesh: mesh parameters
% x1,x2,y1,y2: coordinates for the active region

B = bas.B;      % bases function values in Gauss points times Gauss weights
Bp = bas.Bp;    % bases function derivative values in Gauss points times Gauss weights
B0 = bas.GW;    % orthogonal part

BxNN = bas.BxNN; % total number of bases functions in x-direction
BxN = bas.BxN;  % total number of basis function types in the x-direction per area
Bxn = bas.Bxn;  % x-number of a basis function Bxn(b#,t#)
Bxt = bas.Bxt;  % type of bases function

ByNN = bas.ByNN;
ByN = bas.ByN;
Byn = bas.Byn;
Byt = bas.Byt;

txN = mesh.txN;     % number of rectangles in x-direction
ixN = bas.ixN;     % number of integration points
xi = bas.xi;       % x-coordinates xi(int-#,rec-#)
dx = mesh.dx;
tyN = mesh.tyN;
iyN = bas.iyN;
yi = bas.yi;
dy = mesh.dy;


Ax1 = zeros(BxNN*tyN,1);
Ay1 = zeros(ByNN*txN,1);
Ax2 = zeros(BxNN*tyN,1);
Ay2 = zeros(ByNN*txN,1);
%dt = singdist(center);


A = 0; % area
for nx = 1:txN  % över (basfunktions) områden x
    for nxi = 1:1%ixN  % över integrationspunkter x
        x = xi(nxi,nx);  % x-coord
        dx1 = dx(nx); % dx/??
        for ny = 1:tyN  % över (basfunktions) områden y
            if geoind(nx,ny) == 1               
                dy1 = dy(ny); % dy/??
                A = A + dx1*dy1;
                for nyi = 1:1%iyN  % över integrationspunkter y
                    y = yi(nyi,ny);  % y-coord
                    
                    for bb1 = 1:BxN % över basfunktioner i x-led
                        jj = Bxn(bb1,nx);
                        if jj>0
                            jj = jj + (ny-1)*BxNN;
                            B1 = B(nxi,Bxt(bb1,nx))*B0(nyi)*dx1*dy1;
                            Bp1 = Bp(nxi,Bxt(bb1,nx))*B0(nyi)*dy1;
                            Ax1(jj) = Ax1(jj) + abs(B1);
                            Ax2(jj) = Ax2(jj) + Bp1;
%                             disp(strcat('x:',num2str([nx ny jj])));
                        end
                    end
                    for bb1 = 1:ByN % över basfunktioner i y-led
                        jj = Byn(bb1,ny);
                        if jj>0
                            jj = jj + (nx-1)*ByNN;
                            B1 = B(nyi,Byt(bb1,ny))*B0(nxi)*dy1*dx1;
                            Bp1 = Bp(nyi,Byt(bb1,ny))*B0(nxi)*dx1;
                            Ay1(jj) = Ay1(jj) + abs(B1);
                            Ay2(jj) = Ay2(jj) + Bp1;
%                             disp(strcat('y:',num2str([nx ny jj])));
                        end
                    end
                end
            end
        end
    end
end % (basfunktioner) omraden

A2 = [Ax2; Ay2];
Axp = [Ax1; Ay1];
Ayp = [Ay1; Ax1];
ind2 = find(Axp==0);
ind1 = find(Axp>0);
temp = max(Axp);
ind3 = find(Axp>0.9*temp);
% figure(2); clf
% imagesc(Z21)
% colorbar
