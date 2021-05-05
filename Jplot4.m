function [Jx,Jy,xx,yy,rho]=Jplot4(bas,meshp,Jxopt,px,ixN,iyN,opts)
% plotts current density from MoM evaluted current
% ixN, iyN number of evaluation points in each rectangle
% px quiver points, 0 non, 1 all, n, each nth-rec. vector for quiver
% 0 < px < 1 quiver in points with

% colormap(ho)
hott=hot(256);
hottt=hott(size(hott,1):-1:1,:);
ho=[1-hottt(2:end,:); hottt];
if nargin < 7
    opts.arrowwidth = 0.75;
    opts.re = 1;
    opts.im = 1;
    opts.arrowcolor = 'k';
    if exist('myseqcmap')
        opts.colormap = 'myseqcmap';
    else
        opts.colormap = hottt;
    end
    opts.geo = 0;
    opts.geoc = 0;  % restrict to geo if 1
    if nargin < 6
        ixN = 1;  % center point
        if nargin < 5
            iyN = 1;
            if nargin < 4
                px = 1;
            end
        end
    end
else
    if ~isfield(opts,'re') && ~isfield(opts,'im')
        opts.re = 1;
        opts.im = 0;
    end
    if ~isfield(opts,'re')
        opts.re = 0;
    end
    if ~isfield(opts,'im')
        opts.im = 0;
    end
    if ~isfield(opts,'geoc')
        opts.geoc = 0;
    end
    if isfield(opts,'ind')
        Jxopt1 = zeros(1,meshp.txN*(meshp.tyN-1)+meshp.tyN*(meshp.txN-1));
        Jxopt1(opts.ind) = Jxopt;
        Jxopt = Jxopt1;
    end
    
    if ~isfield(opts,'arrowwidth')
        opts.arrowwidth = 0.75;
    end
    if ~isfield(opts,'arrowcolor')
        opts.arrowcolor = 'k';
    end
    if ~isfield(opts,'colormap')
        if exist('myseqcmap')
            opts.colormap = 'myseqcmap';
        else
            opts.colormap = 'hottt';
        end
    else
        if opts.colormap == 'hottt'
            opts.colormap = hottt;
        end
    end
end

BxNN = bas.BxNN; % total number of bases functions in x-direction
BxN = bas.BxN;  % total number of basis function types in the x-direction per area
Bxn = bas.Bxn;  % x-number of a basis function Bxn(b#,t#)
Bxt = bas.Bxt;  % type of bases function

ByNN = bas.ByNN;
ByN = bas.ByN;
Byn = bas.Byn;
Byt = bas.Byt;

txN = meshp.txN;     % number of rectangles in x-direction
dx = meshp.dx;
tyN = meshp.tyN;
dy = meshp.dy;


x = meshp.x;
xm = meshp.xm;
y = meshp.y;
ym = meshp.ym;

% local point
ix = linspace(0,1,ixN+1);
ix = ix(2:end)-ix(2)/2;
iy = linspace(0,1,iyN+1);
iy = iy(2:end)-iy(2)/2;

Bx = zeros(ixN,8);
Bxp = zeros(ixN,8);
By = zeros(iyN,8);
Byp = zeros(iyN,8);
for n=1:8
    [Bx(:,n),Bxp(:,n),~,~] = basfunktion1(ix,n);
    [By(:,n),Byp(:,n),~,~] = basfunktion1(iy,n);
end

% size(Bx)
% global evaluation points
xi = zeros(ixN,txN); % iN = # of integrationspunkter
yi = zeros(iyN,tyN);
for n=1:txN
    xi(:,n) = xm(n)+dx(n)*(ix-0.5);
end
for n=1:tyN
    yi(:,n) = ym(n)+dy(n)*(iy-0.5);
end


Jx = zeros(txN*ixN,tyN*iyN); % current evaluated in the evaluation points
Jy = zeros(txN*ixN,tyN*iyN); % current evaluated in the evaluation points
Jxp = zeros(txN*ixN,tyN*iyN); % current evaluated in the evaluation points
Jyp = zeros(txN*ixN,tyN*iyN); % current evaluated in the evaluation points
Jxc = zeros(txN,tyN); % current evaluated in the center points
Jyc = zeros(txN,tyN); % current evaluated in the center points
rho = zeros(txN*ixN,tyN*iyN);
rhoc = zeros(txN,tyN);
xx = zeros(txN*ixN,1);
yy = zeros(tyN*iyN,1);
m = 0;
rhotot = 0;
for nx = 1:txN  % regions in x
    for nxi = 1:ixN  % over evaluation points in x
        m = m+1; % x-pos
        xx(m) = xi(nxi,nx);
        n = 0;
        for ny = 1:tyN
            if isfield(opts,'geoind') && opts.geoc
                ingeo = ~opts.geoind(nx,ny);
            else
                ingeo = 1;
            end
            for nyi = 1:iyN  % over evaluation points in y
                n = n+1;
                yy(n) = yi(nyi,ny);
                if ingeo
                    for bb2 = 1:BxN % över basfunktioner i x-led
                        ii = Bxn(bb2,nx);
                        if ii>0
                            ii = ii + (ny-1)*BxNN;
                            %                         [m n]
                            Jx(m,n) = Jx(m,n) + Jxopt(ii)*Bx(nxi,Bxt(bb2,nx))/dy(ny);
                            rho(m,n) = rho(m,n) + Jxopt(ii)*Bxp(nxi,Bxt(bb2,nx))/dy(ny)/dx(nx);
                            Jxc(nx,ny) = Jxc(nx,ny) + Jxopt(ii)/2/dy(ny);
                            rhotot = rhotot + Jxopt(ii)*Bxp(nxi,Bxt(bb2,nx))/dy(ny)/dx(nx);
                        end
                    end
                    for bb2 = 1:ByN % över basfunktioner i y-led
                        ii = Byn(bb2,ny);
                        if ii>0
                            ii = ii + (nx-1)*ByNN + tyN*BxNN;
                            Jy(m,n) = Jy(m,n) + Jxopt(ii)*By(nyi,Byt(bb2,ny))/dx(nx);
                            rho(m,n) = rho(m,n) + Jxopt(ii)*Byp(nyi,Byt(bb2,ny))/dy(ny)/dx(nx);
                            Jyc(nx,ny) = Jyc(nx,ny) + Jxopt(ii)/2/dx(nx);
                            rhotot = rhotot + Jxopt(ii)*Byp(nyi,Byt(bb2,ny))/dx(nx)/dy(ny);
                        end
                    end
                else
                    Jyp(m,n) = nan(1);
                    Jxp(m,n) = nan(1);
                end
            end
        end
    end
end

Jamp = sqrt(abs(Jx).^2+abs(Jy).^2);
Jnan = sqrt(abs(Jxp).^2+abs(Jyp).^2);
% Jamp = sqrt(abs(Jx).^2+abs(Jy).^2);
Jam = max(max(Jamp));

if Jam>0
    imAlpha=ones(size(Jnan));
    imAlpha(isnan(Jnan))=0;
    imagesc(xx,yy,Jamp','AlphaData',imAlpha',[0 Jam]);
    set(gca,'color',[1 1 1]);
else
%     Jamp = sqrt(abs(Jxc).^2+abs(Jyc).^2);
%     Jam = max(max(Jamp));
%     pcolor(meshp.xm,meshp.ym,Jamp')
    W = [1 2 1; 2 4 2; 1 2 1]/8;
    Jamp = conv2(Jamp,W,'same');
    pcolor(xx,yy,Jamp')
    shading interp
end
% quiver plots
Jamp = sqrt(abs(Jxc).^2+abs(Jyc).^2);
Jam = max(max(Jamp));

if px~=0 % add quiver plot
    hold on
    if length(px)>1
        % to be done        
    elseif px<1 % takes only large values
%         ind = find(Jamp>px*Jam);
%         Jx1 = Jxc(ind);
%         Jy1 = Jyc(ind);
%         [Xp,Yp] = ndgrid(meshp.xm,meshp.ym);
%         Xp1 = Xp(ind);
%         Yp1 = Yp(ind);
%         dd = min(meshp.dx);
        Ja1 = zeros(size(Jamp)+2);
        Jx1 = [];%zeros(size(Jamp));
        Jy1 = [];%zeros(size(Jamp));
        Xp1 = [];
        Yp1 = [];
        Ja1(2:end-1,2:end-1) = Jamp/Jam;
        [Ja1m,i1,i2] = maxm(Ja1);
        while Ja1m>px
            Jx1 = [Jx1 Jxc(i1-1,i2-1)];
            Jy1 = [Jy1 Jyc(i1-1,i2-1)];
            Xp1 = [Xp1 meshp.xm(i1-1)];
            Yp1 = [Yp1 meshp.ym(i2-1)];
            for m1=-1:1
                for m2=-1:1
                    Ja1(i1+m1,i2+m2) = 0;
                end
            end
            [Ja1m,i1,i2] = maxm(Ja1);            
        end
        dd = 2*min(meshp.dx);        
    else
        [xN,yN] = size(Jxc);
        xqN = floor(xN/(px)); % number of arrows
        i0x = round((xN-xqN*px)/2);
        if i0x < px/2
            xqN = xqN-1;
            i0x = round((xN-xqN*px)/2);
        end        
        yqN = floor(yN/(px)); % number of arrows
        i0y = round((yN-yqN*px)/2);
        if i0y < px/2
            yqN = yqN-1;
            i0y = round((yN-yqN*px)/2);
        end        
        Jx1 = Jxc(i0x:px:end-i0x+1,i0y:px:end-i0y+1);
        Jy1 = Jyc(i0x:px:end-i0x+1,i0y:px:end-i0y+1);
        [Xp1,Yp1] = ndgrid(meshp.xm(i0x:px:end-i0x+1),meshp.ym(i0y:px:end-i0y+1));
        Jamp = sqrt(abs(Jx1).^2+abs(Jy1).^2);
        Jam = max(max(Jamp));
        ind = find(Jamp>Jam/20);  % remove small values
        Jx1 = Jx1(ind);
        Jy1 = Jy1(ind);
        Xp1 = Xp1(ind);
        Yp1 = Yp1(ind);
        dd = min(meshp.dx)*px;
    end
    ma = sqrt(max(max(abs(Jx1).^2+abs(Jy1).^2)));  % max length
    S  = 0.9*dd/ma;
    Jx1 = S*Jx1;
    Jy1 = S*Jy1;
    %     if isfield(opts,'arrowcolor')
    if opts.re
        [Jm,ind] = max(abs(Jxopt));
        Jv = Jxopt(ind);
        Jv = Jv/abs(Jv);
        if opts.im
            Jv=1;
        end
        Jx1 = Jx1/Jv;
        Jy1 = Jy1/Jv;
        quiver(real(Xp1-Jx1/2).',real(Yp1-Jy1/2).',real(Jx1).',real(Jy1).',0,opts.arrowcolor,'filled','linewidth',opts.arrowwidth)
    end
    if opts.im
        quiver(Xp1'-imag(Jx1/2).',Yp1'-imag(Jy1/2).',imag(Jx1).',imag(Jy1).',0,'b','filled','linewidth',opts.arrowwidth)
    end
    %     else
    %         quiver(Xp1'-Jx1.'/2,Yp1'-Jy1.'/2,Jx1.',Jy1.',0,'-k','filled','linewidth',0.73)
    %     end
end
box off
axis tight
axis equal

%colormap(opts.colormap);
if isfield(opts,'geoind')
    xN = txN; yN = tyN;
    geoind = opts.geoind;
    ig = []; % in geo
    eg = []; % out geo
    for nx = 1:xN  % regions in x
        for ny = 1:yN
            ag = [nx nx+1 ny ny; nx nx ny ny+1; nx+1 nx+1 ny ny+1; nx nx+1 ny+1 ny+1];
            if ~geoind(nx,ny)
                ig = [ig; ag];
            else
                eg = [eg; ag];
            end
        end
    end
    % add boundary
    for ny = 1:yN
        ag = [1 1 ny ny+1; xN+1 xN+1 ny ny+1];
        eg = [eg; ag];
    end
    for nx = 1:xN
        ag = [nx nx+1 1 1; nx nx+1 yN+1 yN+1];
        eg = [eg; ag];
    end
    
    bg = intersect(ig,eg,'rows');
    hold on
    for m=1:size(bg,1)
        plot(meshp.x(bg(m,1:2)),meshp.y(bg(m,3:4)),'k')
    end
end

if isfield(opts,'geobas')
    gb = sort(opts.geobas);
    gbxi = find(gb<=meshp.tyN*bas.BxNN);
    gbx = gb(gbxi);
    gbyi = find(gb>meshp.tyN*bas.BxNN);
    gby = gb(gbyi);
%     gbx
%     gby
    for nx = 1:txN  % regions in x
        for ny = 1:tyN
            for bb2 = 1:2 % över basfunktioner i x-led
                ii = Bxn(bb2,nx);
                if ii>0
                    ii = ii + (ny-1)*BxNN;
                    if min(abs(ii-gbx))<1e-5
                        if bas.Bxt(bb2,nx) == 6 % increasing
                            plot([meshp.x(nx+1) meshp.x(nx+1)],[meshp.y(ny) meshp.y(ny+1)],'g','linewidth',1)
                        else
                            plot([meshp.x(nx) meshp.x(nx)],[meshp.y(ny) meshp.y(ny+1)],'g','linewidth',1)
                        end
                    end
                end
            end
            for bb2 = 1:2 % över basfunktioner i y-led
                ii = Byn(bb2,ny);
                if ii>0
                    ii = ii + (nx-1)*ByNN + tyN*BxNN;
                    if min(abs(ii-gby))<1e-5
                        if bas.Byt(bb2,ny) == 6 % increasing
                            plot([meshp.x(nx) meshp.x(nx+1)],[meshp.y(ny+1) meshp.y(ny+1)],'g','linewidth',1)
                        else
                            plot([meshp.x(nx) meshp.x(nx+1)],[meshp.y(ny) meshp.y(ny)],'g','linewidth',1)
                        end
                    end
                end
            end
        end
    end    
end

if isfield(opts,'feed')
    gb = sort(opts.feed);
    gbxi = find(gb<=meshp.tyN*bas.BxNN);
    gbx = gb(gbxi);
    gbyi = find(gb>meshp.tyN*bas.BxNN);
    gby = gb(gbyi);
%     gbx
%     gby
    for nx = 1:txN  % regions in x
        for ny = 1:tyN
            for bb2 = 1:2 % över basfunktioner i x-led
                ii = Bxn(bb2,nx);
                if ii>0
                    ii = ii + (ny-1)*BxNN;
                    if min(abs(ii-gbx))<1e-5
                        if bas.Bxt(bb2,nx) == 6 % increasing
                            plot([meshp.x(nx+1) meshp.x(nx+1)],[meshp.y(ny) meshp.y(ny+1)],'c','linewidth',4)
                        else
                            plot([meshp.x(nx) meshp.x(nx)],[meshp.y(ny) meshp.y(ny+1)],'c','linewidth',4)
                        end
                    end
                end
            end
            for bb2 = 1:2 % över basfunktioner i y-led
                ii = Byn(bb2,ny);
                if ii>0
                    ii = ii + (nx-1)*ByNN + tyN*BxNN;
                    if min(abs(ii-gby))<1e-5
                        if bas.Byt(bb2,ny) == 6 % increasing
                            plot([meshp.x(nx) meshp.x(nx+1)],[meshp.y(ny+1) meshp.y(ny+1)],'c','linewidth',4)
                        else
                            plot([meshp.x(nx) meshp.x(nx+1)],[meshp.y(ny) meshp.y(ny)],'c','linewidth',4)
                        end
                    end
                end
            end
        end
    end        
end

if isfield(opts,'rembas')
    gb = sort(opts.rembas);
    gbxi = find(gb<=meshp.tyN*bas.BxNN);
    gbx = gb(gbxi);
    gbyi = find(gb>meshp.tyN*bas.BxNN);
    gby = gb(gbyi);
%     gbx
%     gby
    for nx = 1:txN  % regions in x
        for ny = 1:tyN
            for bb2 = 1:2 % över basfunktioner i x-led
                ii = Bxn(bb2,nx);
                if ii>0
                    ii = ii + (ny-1)*BxNN;
                    if min(abs(ii-gbx))<1e-5
                        if bas.Bxt(bb2,nx) == 6 % increasing
                            plot([meshp.x(nx+1) meshp.x(nx+1)],[meshp.y(ny) meshp.y(ny+1)],'w','linewidth',3)
                        else
                            plot([meshp.x(nx) meshp.x(nx)],[meshp.y(ny) meshp.y(ny+1)],'w','linewidth',3)
                        end
                    end
                end
            end
            for bb2 = 1:2 % över basfunktioner i y-led
                ii = Byn(bb2,ny);
                if ii>0
                    ii = ii + (nx-1)*ByNN + tyN*BxNN;
                    if min(abs(ii-gby))<1e-5
                        if bas.Byt(bb2,ny) == 6 % increasing
                            plot([meshp.x(nx) meshp.x(nx+1)],[meshp.y(ny+1) meshp.y(ny+1)],'w','linewidth',3)
                        else
                            plot([meshp.x(nx) meshp.x(nx+1)],[meshp.y(ny) meshp.y(ny)],'w','linewidth',3)
                        end
                    end
                end
            end
        end
    end    
end


axis ij