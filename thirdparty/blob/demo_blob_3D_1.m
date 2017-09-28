clear;
NOISY = 0;
DISPLAY = 1;
RELERR = 1;
V = -1+2*rand(9,1000);
%2-10
r = 2  + 8.*rand(3,1000);
r = sort(r,1);
d=rand(1000,1);
c=-1+2*rand(1000,1);
W=[61 61 61];
ori =   W(1)*rand(3,1000)/2;
k = 1;
RE=[];
for j=1:1000
    v=orth(reshape(V(:,j),3,3));
    D=diag(1./r(:,j).^2);
    t=v*D*v';
    
    I1=gauss3(t,ori(:,j),W);
    if DISPLAY
        figure(1);clf;
        subplot(1,2,1);
        D=I1;
        D=uint8(D*255);
        Ds = smooth3(D);
        hiso = patch(isosurface(Ds,5),...
            'FaceColor',[1,.75,.65],...
            'EdgeColor','none','FaceAlpha',0.1);
        isonormals(Ds,hiso)
        hcap = patch(isocaps(D,5),...
            'FaceColor','red',...
            'EdgeColor','none');
        view(35,30);
        axis([1 W(1) 1 W(1) 1 W(1)]);
        daspect([1,1,1])
        lightangle(45,30);
        camlight; lighting gouraud;
        alpha(0.5);
        title('expected');
    end
    I1=I1-min(I1(:)) ;
    I1=I1/max(I1(:)) ;
    I1=I1*c(j) + d(j) ;
    if NOISY
        I1=imnoise(I1,'gaussian');
        sigma0 = 1.6;
    else
        sigma0 = .6;
    end
    
    [frames1] = detector_3D( I1,'sigma0',sigma0,'checkboundary',0) ;
    if isempty(frames1)
        continue;
    end
    I=zeros(size(I1));
    for ii=1:size(frames1,2)
        m2=frames1([11 12 13],ii);
        t=reshape(frames1(11:19,ii),3,3);
        I1=gauss3(t,[frames1(1,ii)-1; frames1(2,ii)-1;frames1(3,ii)-1],W);
        flg=I1>1e-6;
        I(flg)=I1(flg);
    end
    if DISPLAY
        figure(1);
        subplot(1,2,2);
        hold on;
        D=I;
        D=uint8(D*255);
        Ds = smooth3(D);
        hiso = patch(isosurface(Ds,5),...
            'FaceColor',[1,.75,.65],...
            'EdgeColor','none','FaceAlpha',0.1);
        isonormals(Ds,hiso)
        hcap = patch(isocaps(D,5),...
            'FaceColor','red',...
            'EdgeColor','none');
        view(35,30);
        axis([1 W(1) 1 W(1) 1 W(1)]);
        daspect([1,1,1])
        lightangle(45,30);
        camlight; lighting gouraud;
        alpha(0.5);
        title('detected');
    end
end