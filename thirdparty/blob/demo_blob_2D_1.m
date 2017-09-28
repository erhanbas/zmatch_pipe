clear;
NOISY = 1;
ISO=0;
RELERR = 1;
baseline=rand(1000,1);


WINSIZ=127;
O = WINSIZ*2/3*rand(2,1000);
r = 2  + 8.*rand(2,1000);
r = sort(r,1);
V = -1+2*rand(4,1000);

k = 1;
for j=1:100
    d = baseline(j);
    c = -d + rand();
    v=orth(reshape(V(:,j),2,2));
    D=diag(1./r(:,j).^2);
    t=v*D*v';
    I1=gauss21(t,O(:,j),[WINSIZ WINSIZ]);
    I1=I1-min(I1(:)) ;
    I1=I1/max(I1(:)) ;
    I1=I1*c + d ;
    if NOISY
        I1= imnoise(I1,'gaussian');
        sigma0=1.6;
    else
        sigma0=.6;
    end
    figure(1);
    clf;
    subplot(1,2,1);
    imshow(I1);
    title('expected');
    [frames1] = detector_2D( I1,'min_max',.9,'sigma0',sigma0) ;
    if isempty(frames1)
        continue;
    end
    I=zeros(size(I1));
    subplot(1,2,2);
    imshow(I1);
    hold on;
    for ii=1:size(frames1,2)
        x0=frames1(1,ii);
        y0=frames1(2,ii);
        Mi=reshape(frames1(10:13,ii),2,2);
        I1=gauss21(Mi,[x0,y0],[WINSIZ WINSIZ]);
        drawellipse(Mi,x0,y0);
        flg=I1>1e-6;
        I(flg)=I1(flg);
    end
%     imshow(I);
    drawnow;
    title('detected');
    hold off;
    pause(.5);
end