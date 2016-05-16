

ARA_X    ARA_Z  ARA_Y     Plate#   x_pxl  z_pxl  yz_pxl
-40      520    20       54        238    62     214
0     2090    20       54        231   118     214
0     4260    20       54        230   195     214

500      860  -5995     113        241    76     449
4360     4840  -5995     113        382   212     449
10     6770  -5995     113        229   282     449

a25 = [ 238    62   -214
    231   118   -214
    230   195   -214
    241    76   -449
    382   212   -449
    229   282   -449]';


ARA = [  -40         520          20
    0        2090          20
    0        4260          20
    500         860       -5995
    4360        4840       -5995
    10        6770       -5995]';

for j0=5%1:6  % j0=5  & ii=[1:6] optimize pxl error
    for k=1:3
        switch k
            case 1
                ii=setdiff([1:6],j0);
            case 2
                ii=[2 4 6];
            case 3
                ii=[1:6];
        end
        a=(ARA-repmat(ARA(:,j0),[1 size(ARA,2)]));
        b=(a25-repmat(a25(:,j0),[1 size(a25,2)]));
        
        M=b(:,ii)/a(:,ii); B=round(mean((M*ARA-a25)')');
        a25_mdl = round(M*ARA-repmat(B,[1,size(ARA,2)]));
        err(:,:,k) = round(a25_mdl-a25)./a25;
    end;
    squeeze(sum(round(sqrt(sum(err.^2,1))*1000),2))
end;

%% The transformation matrix that maps ARA coordinates into the pixel coords of Atlas25 image cube
M = [  0.0358   -0.0000    0.0012;
    -0.0007    0.0350   -0.0001;
    0         0       0.0391];

B = [ -233  -45  215]';

% pxl error round(a25_mdl-a25)
err_pxl=[  -6     2     3     3     0    -3;
    1     0    -1     0     0     1;
    0     0     0     0     0     0]

% pmds_Oct=ubu(:,1); x_Oct=ubu(:,2); y_Oct=ubu(:,3); z_Oct=ubu(:,4);
% xzy_a25_Oct = round(M*[x_Oct z_Oct y_Oct]'-repmat(B,[1,numel(x_Oct)]));
% xzy_a25_Oct(3,:)=-xzy_a25_Oct(3,:);
%
% pmds_new=ubu(:,1); x_new=ubu(:,2); y_new=ubu(:,3); z_new=ubu(:,4);
% xzy_a25_new = round(M*[x_new z_new y_new]'-repmat(B,[1,numel(x_new)]));
% xzy_a25_new(3,:)=-xzy_a25_new(3,:);
%
%%
load Atlas25.mat
load AtlasMask25.mat
load Annotation.mat
load colormap_Annotation
% load the matlab struct containing the most recent Hits analysis result,e.g.
load 'InjHits_Analysis_18-Jan-2016.mat'

%pmds=ubu(:,1); x=ubu(:,2); y=ubu(:,3); z=ubu(:,4);  injNo=ubu(:,5);
InjPar_AAV.DoubleInj=logical(InjPar_AAV.DoubleInj);
pmds={InjPar_AAV.brnID{:},InjPar_AAV.brnID{InjPar_AAV.DoubleInj}};
injNo=[InjPar_AAV.Ainj, InjPar_AAV.Ainj2(InjPar_AAV.DoubleInj)];
x=[InjPar_AAV.x, InjPar_AAV.x2(InjPar_AAV.DoubleInj)]';
x=abs(x); % map RHemi injections into LHemi
y=[InjPar_AAV.y, InjPar_AAV.y2(InjPar_AAV.DoubleInj)]';
z=[InjPar_AAV.z, InjPar_AAV.z2(InjPar_AAV.DoubleInj)]';

%*******************************************************
% May-01 special:    
% The t-date Ntot = 847 injection coverage still has holes.
% But the final coverage is desired to be shown for the pre-SFN abstract presentation
% To this end, the coverage density is to show the final coverage AFTER we 
% add in n=25 Augm injection already in pipeline (n=12) or 
% planned to inject to fill holes in coverage(n=13)
% 
x847=x; y847=y; z847=z;
ii=[6 7 9 10 11 14:20 22:34];
xAugm=Augm_x(ii); yAugm=Augm_y(ii); zAugm=Augm_z(ii);
x=[x; xAugm]; y=[y; yAugm]; z=[z; zAugm]; 
%*******************************************************

xzy_a25 = round(M*[x z y]'-repmat(B,[1,numel(x)]));
% The pixel coordinate in the rostrocausal Y dim has to be positive so flip it
xzy_a25(3,:)=-xzy_a25(3,:);
% Note that the mediolateral coordinates xzy_a25(1,:) are positive to reflect LHemi injections
% Note that the Annotation is made of the RHemi (and reflected back to the Left Hemi (while keeping the symmetric original RHemi un-annotated) for the color plates in the printed atlas)
% Thus we flip the injection x pxl-coordinate to map into the RHemi (the annotated half) of the mask
[ny nz nx]=size(AtlasMask25);
xzy_a25(1,:)=nx+1-xzy_a25(1,:);

for i=1:ny
    ii=find(xzy_a25(3,:)==i);
    clf;
    % sagittal
    % imagesc(squeeze(Atlas25(:,:,i))'); colormap gray; axis image;
    % coronal
    %imagesc(squeeze(Atlas25(i,:,:))); colormap gray; axis image;
    imagesc(squeeze(AtlasMask25(i,:,:))); colormap gray; axis image;
    hold on;
    aa=xzy_a25([1 2],ii);
    plot(aa(1,:),aa(2,:),'r.','markersize',[24]);
    title(sprintf('Y = %d  um',(i-213)*25));
    pause(0.1);
    %waitforbuttonpress;
end;

%% Calculate the Inj Density in the RHemi
iiBrn=find(AtlasMask25)';
[yMsk,zMsk,xMsk] = ind2sub(size(AtlasMask25),iiBrn);
x_pxl=1; %28.5;
z_pxl=1; %28.5;
y_pxl=1; %26.0;
D=Inf(size(iiBrn));
% density calculation ~ 2sec/inj (~30 mins for the AAV full data set N~825)
for i=1:size(xzy_a25,2)
    dx=x_pxl*(xMsk-xzy_a25(1,i));
    dy=y_pxl*(yMsk-xzy_a25(3,i));
    dz=z_pxl*(zMsk-xzy_a25(2,i));
    dd=sqrt(dx.^2+dy.^2+dz.^2);
    D=min([D;dd]);
    if mod(i,25)==0 fprintf('...done injection %d\n',i); end;
end
D(find(D==0))=x_pxl/2;
save(sprintf('InjDensity_AAV_%s.mat',date),'D');

%%
% load the matlab struct containing the most recent injection density analysis,e.g.
load('InjDensity_AAV_19-Jan-2016');
iiBrn=find(AtlasMask25)';

% colormap_Annotation
clr_bkg=0;         % neutral gray or black
clr_tissue=5;       % bright cyan color   (for colormap_Annotation)
clr_tissue=250;       % bright cyan color (for colormap gray)
clr_tissue_RH=22;       % a distinct color for Right Hemi
clr_grd=200;         % high contrast with both background and tissue
clr_inj_cntr=80;      % center patch dark color
clr_inj_rim=60;       % rim of 1mm outer radius -- lighter color of same hue as center

InjDensity = AtlasMask25;
ii_tissue=AtlasMask25==255;
ii_bkg=AtlasMask25==0;
InjDensity(ii_tissue)=clr_tissue;
InjDensity(ii_bkg)=clr_bkg;

%InjDensity(InjDensity==255) = 128;
%InjDensity(iiBrn')=uint8(D);
nexp=1;
%InjDensity(iiBrn)=uint8(255*(1-(D/max(D(:))).^(nexp)));
InjDensity(iiBrn)=uint8(255*(min(D(:))./D).^(nexp));
InjDensity=InjDensity(:,:,end:-1:1);

x_pxl=28.5;
z_pxl=28.5;
y_pxl=26.0;

iy_Breg=213;
ix_Breg=230;
iz_Breg=20;


InjDensity0=InjDensity;
II=InjDensity0;

%for i=1:4 ii=find(II>=Icrit(i) & II<Icrit(i-1)); II(ii)=Icrit(i); end;
%II(II>8)=192; II(II<=8&II>5)=32; II(II<5)=4;

II(II>7)=clr_inj_cntr;      % corresponds to r=521 um radius around inj center
II(II<=7&II>4)=clr_inj_rim; % corresponds to r=912 um radius around inj center
ii_farpxl = II<=4;
ii_fartissue=ii_farpxl & ~ii_bkg;
II(ii_fartissue)=clr_tissue;
%II(:,:,1:230)=clr_tissue_RH;
II(ii_bkg)=clr_bkg;

II(:,:,20:35:450)=clr_grd;
II(:,5:35:320,:)=clr_grd;

InjDensity=II;
clear II;

figure(1); clf;
set(gcf,'PaperOrientation','portrait')
set(gcf,'PaperPosition',[ 0.2500    0.2500   10.5000    8.0000]);
colormap(colormap_Annotation);
%colormap(spring(256));
for j=1:1:501 %25:1:525
    i=j+24;
    A=squeeze(InjDensity(i,:,:));
    aL=squeeze(Annotation(i,:,end:-1:1));
    aR=squeeze(Annotation(i,:,:));
    aaR=unique(aR(:));
    aaL=unique(aL(:));
    image(A);
    axis image; hold on;
    contour(aR,aaR,'w-','linewidth',2);
    contour(aL,aaL,'w-','linewidth',2);
    hold off;
    %colormap(colormap_Annotation);
    axis off;
    Ycoronal=-(i-iy_Breg)*y_pxl;
    titlestr=sprintf('rostrocaudal  Y= %4d um',Ycoronal);
    save_image_flnm=sprintf('InjDensAAV_Coronal%03d.png',j);
    title(titlestr);
    set(gca,'xtick',[20:35:450],'xticklabel',[-6000:1000:6000]);
    %set(gca,'xdir', 'reverse');
    set(gca,'ytick',[5:35:320],'yticklabel',[-1000:1000:8000]);
    xlabel('mediolateral  X [um]');
    ylabel('dorsoventral  Z [um]');
    axis on;
    saveas(gcf,save_image_flnm);
    pause(0.3);
    %waitforbuttonpress;
end;

command=['ffmpeg -y -r 10 -f image2 -i ' sprintf('OST%d',k) filesep '%4d_' sprintf('OST%d-F.png',k) ' -r 10 -vb 20M -aspect ' num2str(1.331,3) sprintf(' OST%d_C.mp4',k)]
command=['ffmpeg -y -r 10 -f image2 -i ' sprintf('InjDensAAV_Coronal') filesep 'InjDensAAV_Coronal%3d.png' ' -r 10 -vb 20M -aspect ' num2str(1.3125,3) sprintf(' InjDensAAV_Coronal.mp4')]
system(command)
for i=1:501
    imfl0=sprintf('InjDensAAV_Coronal%03d.png',i-1+25);
    imfl1=sprintf('InjDensAAV_Coronal%03d.png',i);
    A=imread(imfl0);
    imwrite(A,imfl1,'png');
    if mod(i,25)==0 fprintf(' ...done img %4d\n',i); end;
end;


%%

figure(2); clf;
colormap(colormap_Annotation);
% coronal-stack movie
clear A
i_frm=0;
for iy=30:1:510 % 1:ny
    i_frm=i_frm+1;
    counter_str = sprintf('y= %d um',round((iy_Breg-iy)*y_pxl));
    png_flnm=sprintf('%04d_Coronal_AAV-InjDensity.png',i_frm);
    a=squeeze(uint8(InjDensity(iy,:,:)));
    %a(300,350-[1:round(1000/x_pxl)])=220; % scale bar;
    for i=1:3 A(:,:,i)=a; end;
    % A = A(:,end:-1:1,:);   % Horiz-Flip so that LH is shown on the left
    A1 = myTextInsert(A,[340 20],counter_str,10);
    image(A1);
    %colormap gray;
    colormap(colormap_Annotation);
    axis image; axis off;
    %set(gca,'xdir','reverse');  % LH to be shown on the left
    title(counter_str);
    imwrite(A1,png_flnm,'png');
    pause(0.1);
end;

clear A
% sagittal-stack movie
n0=40;
for ix=nx-n0:-1:nx/2-20
    counter_str = sprintf('x= %d um',round((ix-ix_Breg)*x_pxl));
    png_flnm=sprintf('%04d_Sagittal_AAV-InjDensity.png',nx-n0-ix+1);
    a=squeeze(uint8(InjDensity(:,:,ix)))';
    a(:,3:35:528)=clr_grd;
    %a(300,350-[1:round(1000/y_pxl)])=220; % scale bar;
    for i=1:3 A(:,:,i)=a; end;
    A1 = myTextInsert(A,[40 20],counter_str,10);
    %image(A);
    image(A1);
    %colormap gray;
    axis image; axis off;
    title(counter_str);
    imwrite(A1,png_flnm,'png');
    pause(0.1);
end;

clear A
% transverse-stack movie
n0=20;
for iz=n0:nz
    counter_str = sprintf('z= %d um',round((iz-iz_Breg)*z_pxl));
    png_flnm=sprintf('%04d_Transverse_AAV-InjDensity.png',iz-n0);
    a=squeeze(uint8(InjDensity(:,iz,:)));
    a(3:35:528,:)=clr_grd;
    %a(510,40+[1:round(1000/x_pxl)])=220; % scale bar;
    for i=1:3 A(:,:,i)=a; end;
    A1 = A(:,end:-1:1,:);   % Horiz-Flip so that LH is shown on the left
    A1 = myTextInsert(A1,[40 20],counter_str,10);
    image(A1);
    %colormap gray;
    axis image; axis off;
    %set(gca,'xdir','reverse');  % LH to be shown on the left
    title(counter_str);
    imwrite(A1,png_flnm,'png');
    pause(0.1);
end;

% CORONAL AspectRatio:  (nx*x_pxl)/(nz*z_pxl) = (456*28.5)/(320*28.5) = 456/320 = 1.425
% SAGGITAL AspectRatio: (ny*y_pxl)/(nz*z_pxl) = (528*26)/(320*28.5) = 1.5053
% TRANSVERSE AspectRatio: (nx*x_pxl)/(ny*y_pxl) = (528*26)/(456*28.5) = 1.056
ffmpeg -y -r 10 -f image2 -i AAV-InjDensity-Coronal/%04d_Coronal_AAV-InjDensity.png -r 10 -vb 20M -aspect 1.425  AAV-InjDensity_C.mp4
ffmpeg -y -r 10 -f image2 -i AAV-InjDensity-Sagittal/%04d_Sagittal_AAV-InjDensity.png -r 10 -vb 20M -aspect 1.505  AAV-InjDensity_S.mp4
ffmpeg -y -r 10 -f image2 -i AAV-InjDensity-Transverse/%04d_Transverse_AAV-InjDensity.png -r 10 -vb 20M -aspect 1.056  AAV-InjDensity_T.mp4

ffmpeg -y -r 10 -f image2 -i AAV-InjDensity-Coronal_thresholded/%04d_Coronal_AAV-InjDensity.png -r 10 -vb 20M -aspect 1.425  AAV-InjDensity_thresholded_C.mp4
ffmpeg -y -r 10 -f image2 -i AAV-InjDensity-Sagittal_thresholded/%04d_Sagittal_AAV-InjDensity.png -r 10 -vb 20M -aspect 1.505  AAV-InjDensity_thresholded_S.mp4
ffmpeg -y -r 10 -f image2 -i AAV-InjDensity-Transverse_thresholded/%04d_Transverse_AAV-InjDensity.png -r 10 -vb 20M -aspect 1.056  AAV-InjDensity_thresholded_T.mp4

