I1=imreadbw('cameraman.tif');
figure(1);clf;
imshow(I1);
hold on;
[frames1] = detector_2D( I1,'threshold',0.3) ;
for ii=1:size(frames1,2)
    x0=frames1(1,ii);
    y0=frames1(2,ii);
    Mi=reshape(frames1(10:13,ii),2,2);
    drawellipse(Mi,x0,y0);
end

hold off;