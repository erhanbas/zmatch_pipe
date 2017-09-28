fname = 'mri.tif';
info = imfinfo(fname);
num_images = numel(info);
I1=[];
for k = 1:num_images
    A = imread(fname, k);
    I1(:,:,end+1)=im2double(A);
end
% I1=imnoise(I1,'gaussian');
[M,N,K]=size(I1);
W=[N M K];

figure(1);
subplot(1,2,1);
hold on;
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
axis([1 W(1) 1 W(2) 1 W(3)]);
daspect([1,1,1])
lightangle(45,30);
camlight; lighting gouraud;
alpha(0.75);
[frames1] = detector_3D(I1,'checkboundary',0) ;
I=zeros(size(I1));
figure(1);
subplot(1,2,2);
hold on;
for ii=1:size(frames1,2)
    t=reshape(frames1(11:19,ii),3,3);
    I1=gauss3(t,[frames1(1,ii)-1; frames1(2,ii)-1;frames1(3,ii)-1],W);
    flg=I1>1e-6;
    I(flg)=I1(flg);
end
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
axis([1 W(1) 1 W(2) 1 W(3)]);
daspect([1,1,1])
lightangle(45,30);
camlight; lighting gouraud;
alpha(0.75);