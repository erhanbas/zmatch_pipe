function features=detector_3D(I,varargin)
% The rows of features represent:
% 1 x
% 2 y
% 3 z
% 4 interpolated index of detection scale
% 5 detected scale
% 6 (half) first semi-axis of ellipsoid
% 7 (half) second semi-axis of ellipsoid
% 8 (half) third semi-axis of ellipsoid
% 9 contrast
% 10 offset
% 11 -19 t: hessian
sigma0 = .6;
thresh =  0.025;
S      =  3 ;
presmooth = sigma0 * 2^(1/S);
min_scale= 0;
max_scale      =  min(size(I))/6;
operator = 'doh';
features=[];
checkboundary=1;
asp = 10;
like = 14;
min_max = 0.5;
for k=1:2:length(varargin)
    switch lower(varargin{k})
        case 'like'
            like = varargin{k+1} ;
        case 'min_max'
            min_max= varargin{k+1} ;

        case 'asp'
            asp = varargin{k+1} ;
        case 'min_scale'
            min_scale = varargin{k+1} ;
        case 'max_scale'
            max_scale = varargin{k+1} ;
        case 'numlevels'
            S = varargin{k+1} ;
        case 'sigma0'
            sigma0 = varargin{k+1} ;
        case 'presmooth'
            presmooth = varargin{k+1} ;
        case 'threshold'
            thresh = varargin{k+1} ;
        case 'operator'
            operator = lower(varargin{k+1}) ;
        case 'checkboundary'
            checkboundary= varargin{k+1} ;
        otherwise
            error(['Unknown parameter ''' varargin{k} '''.']) ;
    end
end
if presmooth~=0
    I = imsmooth3(I,presmooth);
end
gss = gaussianss_3D(I,min_scale,max_scale,S,sigma0,operator);
%%
if strcmp(operator,'doh')
    TT = 1e-10;
elseif strcmp(operator,'log') || strcmp(operator,'dog')
    TT = 1e-4;
end
[M,N,T,S]=size(gss.op);
jdx = localmax_3D(gss.op, TT) ;
jdx = [jdx localmax_3D(-gss.op,TT)];
idx =[];
for jj=1:length(jdx)
    ii = jdx(jj);
    hess = [gss.Lxx(ii) gss.Lxy(ii) gss.Lxz(ii);gss.Lxy(ii) gss.Lyy(ii) gss.Lyz(ii);gss.Lxz(ii) gss.Lyz(ii) gss.Lzz(ii)];
    [~,D]=eig(hess);
    D=diag(D);
    if sum(D<0)~=3 && sum(D>0)~=3
        continue;
    end
    idx(end+1) = ii;
end

[i,j,k,s] = ind2sub( size( gss.op), idx ) ;
oframes = [j(:)';i(:)';k(:)';s(:)'] ;
oframes(1:4,:) = oframes(1:4,:) - 1;
oframes = refinemx_3D(oframes,gss.op) ;
oframes(1:4,:) = oframes(1:4,:) + 1;
oframes = oframes(:,oframes(1,:)>1 & oframes(2,:)>1 & oframes(3,:)>1 & oframes(4,:)>1 &oframes(1,:)<N & oframes(2,:)<M & oframes(3,:)<T & oframes(4,:)<S);
cord0=floor(oframes(1:4,:));
cord1=ceil(oframes(1:4,:));
cordr=round(oframes(1:4,:));
x=oframes(1,:);
y=oframes(2,:);
z=oframes(3,:);
s=oframes(4,:);
x1=cord0(1,:);
y1=cord0(2,:);
z1=cord0(3,:);
s1=cord0(4,:);
x2=cord1(1,:);
y2=cord1(2,:);
z2=cord1(3,:);
s2=cord1(4,:);
xr=cordr(1,:);
yr=cordr(2,:);
zr=cordr(3,:);
sr=cordr(4,:);
ind=[
    sub2ind(size(gss.op),y1,x1,z1,s1);
    sub2ind(size(gss.op),y1,x1,z2,s1);
    sub2ind(size(gss.op),y2,x1,z1,s1);
    sub2ind(size(gss.op),y2,x1,z2,s1);
    sub2ind(size(gss.op),y1,x2,z1,s1);
    sub2ind(size(gss.op),y1,x2,z2,s1);
    sub2ind(size(gss.op),y2,x2,z1,s1);
    sub2ind(size(gss.op),y2,x2,z2,s1);
    sub2ind(size(gss.op),y1,x1,z1,s2);
    sub2ind(size(gss.op),y1,x1,z2,s2);
    sub2ind(size(gss.op),y2,x1,z1,s2);
    sub2ind(size(gss.op),y2,x1,z2,s2);
    sub2ind(size(gss.op),y1,x2,z1,s2);
    sub2ind(size(gss.op),y1,x2,z2,s2);
    sub2ind(size(gss.op),y2,x2,z1,s2);
    sub2ind(size(gss.op),y2,x2,z2,s2)
    ];
indr=sub2ind(size(gss.op),yr,xr,zr,sr); %20130606

x20=x2-x;
y20=y2-y;
z20=z2-z;
s20=s2-s;
x01=x-x1;
y01=y-y1;
z01=z-z1;
s01=s-s1;
x21=x2-x1;
y21=y2-y1;
z21=z2-z1;
s21=s2-s1;

m=[
    x20 .* y20 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
    x20 .* y20 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21 ;
    x20 .* y01 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
    x20 .* y01 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y20 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y20 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y01 .* z20.*s20 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y01 .* z01.*s20 ./ x21 ./ y21 ./ z21./s21;
    x20 .* y20 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
    x20 .* y20 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21 ;
    x20 .* y01 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
    x20 .* y01 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y20 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y20 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y01 .* z20.*s01 ./ x21 ./ y21 ./ z21./s21 ;
    x01 .* y01 .* z01.*s01 ./ x21 ./ y21 ./ z21./s21;
    ];
LXX = zeros(16,numel(indr));
LYY = zeros(16,numel(indr));
LZZ = zeros(16,numel(indr));
LXY = zeros(16,numel(indr));
LXZ = zeros(16,numel(indr));
LYZ = zeros(16,numel(indr));
OP = zeros(16,numel(indr));
LL = zeros(16,numel(indr));
for kk=1:16
    LXX(kk,:)=gss.Lxx(ind(kk,:));
    LYY(kk,:)=gss.Lyy(ind(kk,:));
    LZZ(kk,:)=gss.Lzz(ind(kk,:));
    LXY(kk,:)=gss.Lxy(ind(kk,:));
    LXZ(kk,:)=gss.Lxz(ind(kk,:));
    LYZ(kk,:)=gss.Lyz(ind(kk,:));
    OP(kk,:)=gss.op(ind(kk,:));
    LL(kk,:)=gss.L(ind(kk,:));
end
% %noise suppression
OPR=repmat(gss.op(indr),16,1);
f1 = (OP./OPR)>0;
f2 = sum(f1)>=like;
OP(~f1)=NaN;
mt=min(abs(OP));
Mt=max(abs(OP));
f3= min_max<(mt./Mt);
f4=f2&f3;
Lxx = zeros(size(indr));
Lyy = zeros(size(indr));
Lzz = zeros(size(indr));
Lxy = zeros(size(indr));
Lxz = zeros(size(indr));
Lyz = zeros(size(indr));
op = zeros(size(indr));
L = zeros(size(indr));
for kk=1:16
    Lxx = Lxx + LXX(kk,:).*m(kk,:);
    Lyy = Lyy + LYY(kk,:).*m(kk,:);
    Lzz = Lzz + LZZ(kk,:).*m(kk,:);
    Lxy = Lxy + LXY(kk,:).*m(kk,:);
    Lxz = Lxz + LXZ(kk,:).*m(kk,:);
    Lyz = Lyz + LYZ(kk,:).*m(kk,:);
    op = op + OP(kk,:).*m(kk,:);
    L = L + LL(kk,:).*m(kk,:);
end
nanidx =isnan(Lxx);
Lxx(nanidx)=gss.Lxx(indr(nanidx));
Lyy(nanidx)=gss.Lyy(indr(nanidx));
Lzz(nanidx)=gss.Lzz(indr(nanidx));
Lxy(nanidx)=gss.Lxy(indr(nanidx));
Lxz(nanidx)=gss.Lxz(indr(nanidx));
Lyz(nanidx)=gss.Lyz(indr(nanidx));
op(nanidx)=gss.op(indr(nanidx));
L(nanidx)=gss.L(indr(nanidx));

oframes1=[];
for ii=1:length(Lxx)
    if ~f4(ii)
        continue;
    end
    hess = [Lxx(ii) Lxy(ii) Lxz(ii);Lxy(ii) Lyy(ii) Lyz(ii);Lxz(ii) Lyz(ii) Lzz(ii)];
    [V,D]=eig(hess);
    D=diag(D);
    if sum(D<0)~=3 && sum(D>0)~=3
        continue;
    end
    [B,IX]=sort(abs(D),'descend');
    D=D(IX);
    V= V(:,IX);
    R1=D(1)/D(2);
    R2=D(1)/D(3);
    if strcmp(operator,'doh')
        K=R1 * (R2 - 5 * R1 * (1 + R2)) / (R1 * (-5 + R2) - 5 * R2);
        G=R2 * (5 * R2 + R1 * (-1 + 5 * R2)) / R1 / (-R2 + 5 * R1 * (1 + R2));
        H=-0.1e1 / 0.6e1 + 0.5e1 / 0.6e1 / R1 + 0.5e1 / 0.6e1 / R2;
        c = (-1)^(D(1)>0)*abs(op(ii) ^ (0.1e1 / 0.3e1) * ((1 + H) * (1 + H * K) * (1 + G * H * K)) ^ (0.5e1 / 0.6e1) * G ^ (-0.1e1 / 0.2e1) * H ^ (-0.3e1 / 0.2e1) / K);
    elseif strcmp(operator,'log') || strcmp(operator,'dog')
        K=R1 * (R2 ^ 2 + R1 ^ 2 * (3 + 2 * R2 + 3 * R2 ^ 2)) / (2 * R1 * R2 + 3 * R2 ^ 2 + R1 ^ 2 * (3 + R2 ^ 2));
        G=R2 * (3 * R2 ^ 2 + 2 * R1 * R2 ^ 2 + R1 ^ 2 * (1 + 3 * R2 ^ 2)) / R1 / (R2 ^ 2 + R1 ^ 2 * (3 + 2 * R2 + 3 * R2 ^ 2));
        H=((2 * R1 * R2 + 3 * R2 ^ 2 + R1 ^ 2 * (3 + R2 ^ 2)) / R1 / R2 / (R1 + R2 + R1 * R2)) / 0.2e1;
        c =  -op(ii) * G ^ (-0.1e1 / 0.2e1) / K * (H / (1 + H) / (1 + H * K) / (1 + G * H * K)) ^ (-0.3e1 / 0.2e1) / (3 + H ^ 2 * K * (1 + G + G * K) + 2 * H * (1 + K + G * K));
    end
    if strcmp(operator,'doh')
        d = -c *abs( G * H ^ (0.3e1 / 0.2e1) * K * (G * (1 + H) * (1 + H * K) * (1 + G * H * K)) ^ (-0.1e1 / 0.2e1)) + L(ii);
    elseif strcmp(operator,'log') || strcmp(operator,'dog')
        d =-c * G * H ^ (0.3e1 / 0.2e1) * K * (G * (1 + H) * (1 + H * K) * (1 + G * H * K)) ^ (-0.1e1 / 0.2e1) + L(ii);
    end
    if K<0 ||G<0||H<0 || sqrt(K)>asp || sqrt(G)>asp
        continue;
    end
    if abs(c)<thresh || abs(c)> 1.2
        continue;
    end
    oframes1(:,end+1)=[oframes(:,ii);zeros(size(K));H;K;G;c;d;reshape(V,9,1)];
end
if isempty(oframes1)
    return;
end
oframes1(5,:)=interp1(1:length(gss.scales),gss.scales,oframes1(4,:),'spline');
oframes1(6,:) = oframes1(5,:).*sqrt(oframes1(6,:));
oframes1(7,:) = oframes1(6,:).*sqrt(oframes1(7,:));
oframes1(8,:) = oframes1(7,:).*sqrt(oframes1(8,:));
for ii=1:size(oframes1,2)
    V=reshape(oframes1(11:19,ii),3,3);
    R=oframes1(6:8,ii);
    D=diag(1./R.^2);
    t=V*D*V';
    u=t(1,1);
    v=t(2,2);
    w=t(3,3);
    e=t(1,2);
    f=t(1,3);
    g=t(2,3);
    oframes1(11:19,ii) = t(:);
    if checkboundary
        zboundary=((-2 * e * f * g + g ^ 2 * u + e ^ 2 * w + v * (f ^ 2 - u * w)) / (e ^ 2 - u * v)) ^ (-0.1e1 / 0.2e1);
        yboundary=((-2 * e * f * g + g ^ 2 * u + e ^ 2 * w + v * (f ^ 2 - u * w)) / (f ^ 2 - u * w)) ^ (-0.1e1 / 0.2e1);
        xboundary=((-2 * e * f * g + g ^ 2 * u + e ^ 2 * w + v * (f ^ 2 - u * w)) / (g ^ 2 - v * w)) ^ (-0.1e1 / 0.2e1);
        if oframes1(3,ii) + zboundary >T ||oframes1(3,ii) -zboundary <1 || oframes1(1,ii) + xboundary >T ||oframes1(1,ii) -xboundary <1 || oframes1(2,ii) + yboundary >T ||oframes1(2,ii) -yboundary <1
            continue;
        end
    end
    features(:,end+1)=oframes1(:,ii);
end
end