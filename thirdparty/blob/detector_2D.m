function features=detector_2D(I,varargin)
% The lines of oframes represents:
% 1 x
% 2 y
% 3 interpolated index of detection scale
% 4 detected scale
% 5 (minor semi-axis) width
% 6 (major semi-axis) length
% 7 contrast
% 8 offset
% 9 angle, orientation of the minor semi axis
% 10 u
% 11 e
% 12 e
% 12 v
% ellipse: u x^2 + 2 v x y + w y ^2 = 1

sigma0 = .6;
thresh =  0.025;
S      =  3 ;
presmooth = sigma0 * 2^(1/S);
min_scale= 0;
max_scale      =  min(size(I))/6;
operator = 'doh';
checkboundary=1;
asp = 10;
iso=0;
features = [];
like = 6;
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
        case 'iso'
            iso = varargin{k+1} ;
        otherwise
            error(['Unknown parameter ''' varargin{k} '''.']) ;
    end
end
if presmooth~=0
    I = imsmooth2(I,presmooth);
end
gss = gaussianss_2D(I,min_scale,max_scale,S,sigma0,operator);
if strcmp(operator,'doh')
    TT = 1e-10;
elseif strcmp(operator,'log') || strcmp(operator,'dog')
    TT =1e-4;
end
[M,N,S]=size(gss.op);
jdx = localmax_2D(gss.op, TT) ;
if ~strcmp(operator,'doh')
    jdx = [jdx localmax_2D(-gss.op,TT)];
end
idx =[];
for jj=1:length(jdx)
    ii = jdx(jj);
    hess = [gss.Lxx(ii) gss.Lxy(ii); gss.Lxy(ii) gss.Lyy(ii)];
    [~,D]=eig(hess);
    D=diag(D);
    if sum(D<0)~=2 && sum(D>0)~=2
        continue;
    end
    idx(end+1) = ii;
end

[i,j,s] = ind2sub( size( gss.op), idx ) ;
oframes = [j(:)';i(:)';s(:)'] ;
oframes(1:3,:) = oframes(1:3,:) - 1;
oframes = refinemx_2D(oframes,gss.op) ;
oframes(1:3,:) = oframes(1:3,:) + 1;
oframes = oframes(:,oframes(1,:)>1 & oframes(2,:)>1 & oframes(3,:)>1 &oframes(1,:)<M & oframes(2,:)<N & oframes(3,:)<S);
cord0=floor(oframes(1:3,:));
cord1=ceil(oframes(1:3,:));
cordr=round(oframes(1:3,:));
x=oframes(1,:);
y=oframes(2,:);
s=oframes(3,:);
x1=cord0(1,:);
y1=cord0(2,:);
s1=cord0(3,:);
x2=cord1(1,:);
y2=cord1(2,:);
s2=cord1(3,:);
xr=cordr(1,:);
yr=cordr(2,:);
sr=cordr(3,:);
ind=[
    sub2ind(size(gss.op),y1,x1,s1);
    sub2ind(size(gss.op),y1,x1,s2);
    sub2ind(size(gss.op),y2,x1,s1);
    sub2ind(size(gss.op),y2,x1,s2);
    sub2ind(size(gss.op),y1,x2,s1);
    sub2ind(size(gss.op),y1,x2,s2);
    sub2ind(size(gss.op),y2,x2,s1);
    sub2ind(size(gss.op),y2,x2,s2);
    ];
indr=sub2ind(size(gss.op),yr,xr,sr); %20130606

x20=x2-x;
y20=y2-y;
s20=s2-s;
x01=x-x1;
y01=y-y1;
s01=s-s1;
x21=x2-x1;
y21=y2-y1;
s21=s2-s1;

m=[
    x20 .* y20 .* s20 ./ x21 ./ y21 ./ s21 ;
    x20 .* y20 .* s01 ./ x21 ./ y21 ./ s21 ;
    x20 .* y01 .* s20 ./ x21 ./ y21 ./ s21 ;
    x20 .* y01 .* s01 ./ x21 ./ y21 ./ s21 ;
    x01 .* y20 .* s20 ./ x21 ./ y21 ./ s21 ;
    x01 .* y20 .* s01 ./ x21 ./ y21 ./ s21 ;
    x01 .* y01 .* s20 ./ x21 ./ y21 ./ s21 ;
    x01 .* y01 .* s01 ./ x21 ./ y21 ./ s21;
    ];
LXX = zeros(8,numel(indr));
LYY = zeros(8,numel(indr));
LXY = zeros(8,numel(indr));
OP = zeros(8,numel(indr));
LL = zeros(8,numel(indr));
for kk=1:8
    LXX(kk,:)=gss.Lxx(ind(kk,:));
    LYY(kk,:)=gss.Lyy(ind(kk,:));
    LXY(kk,:)=gss.Lxy(ind(kk,:));
    OP(kk,:)=gss.op(ind(kk,:));
    LL(kk,:)=gss.L(ind(kk,:));
end
% %noise suppression
OPR=repmat(gss.op(indr),8,1);
f1 = (OP./OPR)>0;
f2 = sum(f1)>=like;
OP(~f1)=NaN;
mt=min(abs(OP));
Mt=max(abs(OP));
f3= min_max<(mt./Mt);
f4=f2&f3;
Lxx = zeros(size(indr));
Lyy = zeros(size(indr));
Lxy = zeros(size(indr));
op = zeros(size(indr));
L = zeros(size(indr));
for kk=1:8
    Lxx = Lxx + LXX(kk,:).*m(kk,:);
    Lyy = Lyy + LYY(kk,:).*m(kk,:);
    Lxy = Lxy + LXY(kk,:).*m(kk,:);
    op = op + OP(kk,:).*m(kk,:);
    L = L + LL(kk,:).*m(kk,:);
end
nanidx =isnan(Lxx);
Lxx(nanidx)=gss.Lxx(indr(nanidx));
Lyy(nanidx)=gss.Lyy(indr(nanidx));
Lxy(nanidx)=gss.Lxy(indr(nanidx));
op(nanidx)=gss.op(indr(nanidx));
L(nanidx)=gss.L(indr(nanidx));

oframes1=[];
for ii=1:length(Lxx)
    if ~f4(ii)
        continue;
    end
    hess = [Lxx(ii) Lxy(ii);Lxy(ii) Lyy(ii)];
    [V,D]=eig(hess);
    D=diag(D);
    if sum(D<0)~=2 && sum(D>0)~=2
        continue;
    end
    [B,IX]=sort(abs(D),'descend');
    D=D(IX);
    V= V(:,IX);
    t1=atan(V(2,1)/V(1,1));
    r=D(1)/D(2);
    if strcmp(operator,'doh')
        K = r.^2;
        H = 1./r;
        c =   sqrt(abs(op(ii))) .* ((1 + r) .^ 2) ./ r;
        c(D(1)>0)=c(D(1)>0)*-1;
    elseif strcmp(operator,'log') || strcmp(operator,'dog')
        K = (r + 3 * r .^ 3) ./ (3 + r .^ 2);
        H = ((3 + r .^ 2) ./ r ./ (1 + r)) / 0.2e1;
        c =   -(op(ii).* (3 + 2 * r + 3 * r .^ 2) .^ 2 ./ (1 + r) .^ 2 .* (3 + 10 * r .^ 2 + 3 * r .^ 4) .^ (-0.1e1 / 0.2e1)) / 0.2e1;
    end
    if strcmp(operator,'doh')
        d = -c .* sqrt(r) ./ (0.1e1 + r) + L(ii);
    elseif strcmp(operator,'log') || strcmp(operator,'dog')
        d =  (-c .* sqrt((3 + 10 * r .^ 2 + 3 * r .^ 4)) + ((3 + 2 * r + 3 * r .^ 2) .* L(ii))) ./ (3 + 2 * r + 3 * r .^ 2);
    end
    if K<0 ||H<0 || sqrt(K)>asp
        continue;
    end
    if abs(c)<thresh || abs(c)> 1.2
        continue;
    end
    oframes1(:,end+1)=[oframes(:,ii);zeros(size(K));H;K;c;d;t1;reshape(V,4,1)];
end
if isempty(oframes1)
    return;
end
oframes1(4,:)=interp1(1:length(gss.scales),gss.scales,oframes1(3,:),'spline');
oframes1(5,:) = oframes1(4,:).*sqrt(oframes1(5,:));
oframes1(6,:) = oframes1(5,:).*sqrt(oframes1(6,:));
for ii=1:size(oframes1,2)
    V=reshape(oframes1(10:13,ii),2,2);
    R=oframes1(5:6,ii);
    D=diag(1./R.^2);
    t=V*D*V';
    u=t(1,1);
    v=t(2,2);
    e=t(1,2);
    oframes1(10:13,ii) = t(:);
    if checkboundary
        ybound=abs(u .* (u .* (-e .^ 2 + u .* v)) .^ (-0.1e1 / 0.2e1));
        xbound = abs(sqrt(v) .* (-e .^ 2 + u .* v) .^ (-0.1e1 / 0.2e1));
        if (oframes1(1,ii)+xbound)>N || (oframes1(1,ii) - xbound)<1 || (oframes1(2,ii)+ ybound)>M || (oframes1(2,ii) - ybound)<1
            continue;
        end
    end
    features(:,end+1)=oframes1(:,ii);
end
end