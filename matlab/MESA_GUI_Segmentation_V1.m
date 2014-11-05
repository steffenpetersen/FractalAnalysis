function [xd,yd] = MESA_GUI_Segmentation_V1(Image)


imagesc(Image)
[mask, xi, yi]=roipoly;
% load('tmp');
Im_mask=mask.*Image;
pix_zeros=find(Im_mask==0);
edge=sub2ind(size(Im_mask),round(yi),round(xi));
av_edge=mean(Image(edge));
Im_mask(pix_zeros)=av_edge;
Img=Im_mask/(max(max(Im_mask)));% normalization between 0 and 1.

A=255;
Imgmin  = min(Img(:));
Imgmax  = max(Img(:));
f = (Img-Imgmin)/(Imgmax-Imgmin);
Img=A*f; % normalize to the range of [0,1].
nu=0.001*A^2; % calculates the coefficient of arc length term.
 
%%   III   Level set function.
sigma = 4;% scale parameter that specifies the size of the neighbourhood.
iter_outer=50; 
iter_inner=10;% inner iteration for level set evolution.
timestep=.1;
mu=1;% coefficient for distance regularization term (regularizes the level set function).
c0=1;


initialLSF = c0*ones(size(Img));% initializes level set function.
initialLSF(30:90,50:90) = -c0;
% initialLSF(5:end-5,5:end-5) = -c0;
% initialLSF(size(Img,1)/2-5:size(Img,1)/2+5,size(Img,2)/2-5:size(Img,2)/2+5) = -c0;
u=initialLSF;


%%   IV    Bias field estimation for CMR.
epsilon=1;
b=ones(size(Img));  
K=fspecial('gaussian',round(2*sigma)*2+1,sigma);% gaussian kernel.
KI=conv2(Img,K,'same');
KONE=conv2(ones(size(Img)),K,'same');
[row,col]=size(Img);
N=row*col;

for n=1:iter_outer
    [u, b, C]= lse_bfe(u,Img, b, K,KONE, nu,timestep,mu,epsilon, iter_inner);% calls the lse_bfeg function. 
end

% contour extraction
imagesc(Image)
hold on
[c,ch]=contour(u,[0 0],'r');
     ch=get(ch,'children');
     nc=numel(ch);
     xd=[];
     yd=[];
     
for i=1:nc
     xdraw=get(ch(i),'xdata');
     ydraw=get(ch(i),'ydata');
     xd=cat(1,xd,xdraw);
     yd=cat(1,yd,ydraw);
end
hold off






