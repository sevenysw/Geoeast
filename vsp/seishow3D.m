function seishow3D(D,c,dmin,dmax)

if (~exist('c' ,'var')) 
    c=100;
end

cm = max(D(:));
[n1,n2,n3] = size(D);

di = 0.05; dc = 0.05; dt = 0.004;
% xt = [1, 50, 100, 129, 178, 228];
% xtl= [0, 1, 2, 0, 1 ,2];
% yt = [1, 50, 100, 129, 178, 228];
% ytl= [ 0, 1, 2, 0, 0.2, 0.4];
xt = [1, 40, 65, 104]*2;
xtl= [[0, 40]*di, [0, 40]*dc]*2;
yt = [1, 40, 65, 104]*2;
ytl= [[0 40]*di, [100 140]*dt]*2;

I                   = zeros(n1+n3+1,n2+n3+1);
I(1:n3,1:n2)        = squeeze( D(round(n1/2),:,:) )';
I(:,n2+1)  = cm;
I(1:n3,n2+2:n2+n3+1)  = mean(D(:))*ones(n3,n3);
I(n3+2:n3+n1+1,1:n2)  = squeeze(D(:,:, round(n3/2)));
I(n3+1,:)  = cm;
I(n3+2:n3+n1+1,n2+2:n2+n3+1) = squeeze( D(:,round(n2/2),:));
%figure,
I = clip(I,c,c);

if (~exist('dmin' ,'var')) 
    dmin = min(I(:)) ;
end
if (~exist('dmax' ,'var')) 
    dmax = max(I(:)) ;
end

figure,imagesc(I,[dmin,dmax]);%axis image;
colorbar; colormap(seismic(2)); colormap gray;
set(gca,'XTick',xt,'YTick',yt);
set(gca,'XTickLabel',xtl,'YTickLabel',ytl);
set(gca,'FontSize',12);
xlabel('X (km)                       Y (km)');
ylabel('Time (s)                            Y (km)');
%imagesc(I);axis image; colormap gray; axis off;colorbar;
% figure, plot(I(:,64));
end