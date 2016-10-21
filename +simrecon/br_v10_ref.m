clc
clear
tic

%% inpout parameters
cd('C:\Users\Avizo\Downloads\sim_test_data\');
dir_rowdata='RAWcell1'; %Input the folder of the raw dataset
PSF='ygbead_zp1um_NAp55nap44_ExpPsf.tif'; %Input the file of PSF
dir_SI=['SI_',[dir_rowdata],'PbyP_v10initalDeckKpExpPSF-100_Abs_centerPSFttt']; %Specify the reconstructed folder
orientation=[0]; %Corresponding to the file names
phase=[0,-2/5,-4/5,-6/5,-8/5]; %Check the phase shift
ss=250; %step size in nm%
firstz='100.000'; %name of first folder in z
kp=[337 400,337 463]; %Specify the Kp value if known
%kp=[]
apor=336; %apodization radius
px=102 %Input the pixel size for creating SIMpsf
lambda=488 %Excitation wavelength if known
n=1.33 %Check the medium refractive index
mr=39; %specify the radius of mask

mkdir(dir_SI) 
mkdir([dir_SI,'\summedWF']);
sz=size(dir(dir_rowdata),1)-2;
fileID=fopen([dir_SI,'\info.txt'],'w');

%% load PSF
psf=single(imread([PSF]));
[num idx]=max(psf(:));
[x y] = ind2sub(size(psf),idx);
psf=psf(1:2*(x-1),1:2*(y-1));
psf=abs(psf-100);
psf=psf./sum(sum(psf));

%% SIM reconstruction
kpcal=[];
E0001=single(imread([dir_rowdata,'\',firstz,'\',num2str(orientation(1),'%05.3d'),'1','.tif']));
dim=size(E0001,1);

%% Apodization function
%[vx, vy] = meshgrid(1:dim,1:dim);
%A1=cos(pi/(2*apor)*sqrt((vx-(floor(dim/2)+1))^2+(vy-(floor(dim/2)+1))^2));
%A1(A1<0) = 0;
for x=1:dim;
    for y=1:dim
A1(x,y)=cos(pi/(2*apor)*sqrt((x-(floor(dim/2)+1))^2+(y-(floor(dim/2)+1))^2));

if A1(x,y)<0
    A1(x,y)=0;
end
    end
end

%% mask for crosscorrelation
for x=1:2*dim-1
    for y=1:2*dim-1
        
        if (x-dim)^2+(y-dim)^2<(0.9*mr)^2|(x-dim)^2+(y-dim)^2>(1.1*mr)^2;
            m(x,y)=0;
        else 
            m(x,y)=1;
        end

    end
end

 A000= [3 exp(-i*phase(1,1)*2*pi) exp(+i*phase(1,1)*2*pi) exp(-i*phase(1,1)*pi) exp(+i*phase(1,1)*pi);
        3 exp(-i*phase(1,2)*2*pi) exp(+i*phase(1,2)*2*pi) exp(-i*phase(1,2)*pi) exp(+i*phase(1,2)*pi);
        3 exp(-i*phase(1,3)*2*pi) exp(+i*phase(1,3)*2*pi) exp(-i*phase(1,3)*pi) exp(+i*phase(1,3)*pi);
        3 exp(-i*phase(1,4)*2*pi) exp(+i*phase(1,4)*2*pi) exp(-i*phase(1,4)*pi) exp(+i*phase(1,4)*pi);
        3 exp(-i*phase(1,5)*2*pi) exp(+i*phase(1,5)*2*pi) exp(-i*phase(1,5)*pi) exp(+i*phase(1,5)*pi)];
 A000=A000/9; %A matrix can be improved with the precise intensity of each beam!%

for z=1:1
 z  
    subdir=num2str(str2num(firstz)+(ss/1000)*(z-1),'%5.3f');
    summedWF=zeros(dim,dim);


for r=1:size(orientation,2)
r
    %% load row images
    E0001=single(imread([dir_rowdata,'\',subdir,'\',num2str(orientation(r),'%05.3d'),'1','.tif']));
    E0002=single(imread([dir_rowdata,'\',subdir,'\',num2str(orientation(r),'%05.3d'),'2','.tif']));
    E0003=single(imread([dir_rowdata,'\',subdir,'\',num2str(orientation(r),'%05.3d'),'3','.tif']));
    E0004=single(imread([dir_rowdata,'\',subdir,'\',num2str(orientation(r),'%05.3d'),'4','.tif']));
    E0005=single(imread([dir_rowdata,'\',subdir,'\',num2str(orientation(r),'%05.3d'),'5','.tif']));
    E0001=E0001./4;
    E0002=E0002./4;
    E0003=E0003./4;
    E0004=E0004./4;
    E0005=E0005./4;    
    
    summedWF=summedWF+(E0001+E0002+E0003+E0004+E0005)./5;

    %% deconvolution
    E0001=deconvlucy(E0001,psf,10);
    E0002=deconvlucy(E0002,psf,10);
    E0003=deconvlucy(E0003,psf,10);
    E0004=deconvlucy(E0004,psf,10);
    E0005=deconvlucy(E0005,psf,10);
    %dampar or ringing reduction can be added further%
    
    %% Fourier transform

EE0001=fftshift(fft2(E0001));
EE0002=fftshift(fft2(E0002));
EE0003=fftshift(fft2(E0003));
EE0004=fftshift(fft2(E0004));
EE0005=fftshift(fft2(E0005));

EE0001=EE0001.*A1;
EE0002=EE0002.*A1;
EE0003=EE0003.*A1;
EE0004=EE0004.*A1;
EE0005=EE0005.*A1;

    %% retieve domains
    %D0001=zeros(dim,dim);
    %D0002=zeros(dim,dim);
    %D0003=zeros(dim,dim);
    %D0004=zeros(dim,dim);
    %D0005=zeros(dim,dim);
    
   %EE1=[EE0001(1,1),EE0001(1,2)];
   %EE2=[EE0002(1,1),EE0002(1,2)];
   %EE3=[EE0003(1,1),EE0003(1,2)];
   %EE4=[EE0004(1,1),EE0004(1,2)];
   %EE5=[EE0005(1,1),EE0005(1,2)];
s=[dim^2,1];

EE0001=reshape(EE0001,s);
EE0002=reshape(EE0002,s);
EE0003=reshape(EE0003,s);
EE0004=reshape(EE0004,s);
EE0005=reshape(EE0005,s);
B000=[EE0001,EE0002,EE0003,EE0004,EE0005];
%D000=linsolve(A000,B000.');
A000
sum(B000(:))
D000=mldivide(A000,B000.');
   % parfor x=1:dim;
   %     for y= 1:dim;           
   %         B000=[EE0001(x,y) EE0002(x,y) EE0003(x,y) EE0004(x,y) EE0005(x,y)];
              %B000=[EE1(x),EE2(x),EE3(x),EE4(x),EE5(x)];
            %A000= [3 exp(-i*phase(1,1)*2*pi) exp(+i*phase(1,1)*2*pi) exp(-i*phase(1,1)*pi) exp(+i*phase(1,1)*pi);
            %       3 exp(-i*phase(1,2)*2*pi) exp(+i*phase(1,2)*2*pi) exp(-i*phase(1,2)*pi) exp(+i*phase(1,2)*pi);
            %       3 exp(-i*phase(1,3)*2*pi) exp(+i*phase(1,3)*2*pi) exp(-i*phase(1,3)*pi) exp(+i*phase(1,3)*pi);
            %       3 exp(-i*phase(1,4)*2*pi) exp(+i*phase(1,4)*2*pi) exp(-i*phase(1,4)*pi) exp(+i*phase(1,4)*pi);
            %       3 exp(-i*phase(1,5)*2*pi) exp(+i*phase(1,5)*2*pi) exp(-i*phase(1,5)*pi) exp(+i*phase(1,5)*pi)];
            %A000=A000/9; %A matrix can be improved with the precise intensity of each beam!%
            %D000=linsolve(A000,B000.');   
            %D000=mldivide(A000,B000.');
            %D0001(x,y)=D000(1);
            %D0002(x,y)=D000(2);
            %D0003(x,y)=D000(3);
            %D0004(x,y)=D000(4);
            %D0005(x,y)=D000(5);
       % end
    %end
   
   D0001=D000(1,:);
   D0002=D000(2,:);
   D0003=D000(3,:);
   D0004=D000(4,:);
   D0005=D000(5,:);
D0001=reshape(D0001,[dim,dim]);
D0002=reshape(D0002,[dim,dim]);
D0003=reshape(D0003,[dim,dim]);
D0004=reshape(D0004,[dim,dim]);
D0005=reshape(D0005,[dim,dim]);
     
    %% cross correlation to find kp (use m1 term)
 S0001=zeros(2*dim,2*dim);
 S0002=zeros(2*dim,2*dim);
 S0003=zeros(2*dim,2*dim);
 S0004=zeros(2*dim,2*dim);
 S0005=zeros(2*dim,2*dim);
 
 S0001(floor(dim/2)+1:floor(dim/2)+dim,floor(dim/2)+1:floor(dim/2)+dim)=D0001;

 if size(kp,1)==0
 x1=xcorr2(D0001,D0004);
 x2=xcorr2(D0001,D0005);
 
 x1=x1.*m;
 x2=x2.*m;
 %figure(2)
 %imagesc(abs(x1))
  
 [M,I]=max(abs(x1(:)));
 [y_1,x_1]=ind2sub(size(abs(x1)),I);
 x_1=x_1+1;
 y_1=y_1+1;

 [M,I]=max(abs(x2(:)));
 [y_2,x_2]=ind2sub(size(abs(x2)),I);
 x_2=x_2+1;
 y_2=y_2+1;
 
 shift1=[y_1-(floor(dim/2)+1),x_1-(floor(dim/2)+1)];
 shift2=[y_2-(floor(dim/2)+1),x_2-(floor(dim/2)+1)];
 
 displace1=[y_1-(floor(dim)+1),x_1-(floor(dim)+1)];
 y_3=2*displace1(1)+floor(dim)+1;
 x_3=2*displace1(2)+floor(dim)+1;
 displace2=[y_2-(floor(dim)+1),x_2-(floor(dim)+1)];
 y_4=2*displace2(1)+floor(dim)+1;
 x_4=2*displace2(2)+floor(dim)+1;

 shift3=[y_3-(floor(dim/2)+1),x_3-(floor(dim/2)+1)];
 shift4=[y_4-(floor(dim/2)+1),x_4-(floor(dim/2)+1)];
 kpcal=vertcat(kpcal,[x_1-(floor(dim/2)),y_1-(floor(dim/2)),x_3-(floor(dim/2)),y_3-(floor(dim/2))])
 
 else
 shift1=[kp(r,2)-1,kp(r,1)-1]
 shift2=[dim-shift1(1),dim-shift1(2)]
 shift3=[kp(r,4)-1,kp(r,3)-1]
 shift4=[dim-shift3(1),dim-shift3(2)]
 kpcal=kp;
 end
 eval(['shift',num2str(orientation(r),'%05.3d'),'1','=shift1',';']);
 eval(['shift',num2str(orientation(r),'%05.3d'),'2','=shift2',';']);
 eval(['shift',num2str(orientation(r),'%05.3d'),'3','=shift3',';']);
 eval(['shift',num2str(orientation(r),'%05.3d'),'4','=shift4',';']);

    %% find initial phase for m1 term (somehow it is a little strange, can be further improved)
 WF2x=abs(ifft2(ifftshift(S0001))); %wide-field image for comparison%
 
 figure; imagesc(WF2x); axis equal tight;
 

 
 for t=1:1:20
 ip=(t-1)*0.1*pi;
 S0004(shift1(1)+1:shift1(1)+dim,shift1(2)+1:shift1(2)+dim)=D0004*exp(i*ip);
 S0005(shift2(1)+1:shift2(1)+dim,shift2(2)+1:shift2(2)+dim)=D0005*exp(-i*ip);
 %S0002(shift2(1)+1:shift2(1)+dim,shift2(2)+1:shift2(2)+dim)=D0002*exp(i*ip);
 %S0003(shift1(1)+1:shift1(1)+dim,shift1(2)+1:shift1(2)+dim)=D0003*exp(-i*ip);
 S1=S0001+S0004+S0005;
 I1=abs(ifft2(ifftshift(S1))); 
 intphas1(t)=sum(sum(I1.*WF2x));
 figure(1)
 imagesc(abs(S1))
 end
 [M,I]=max(intphas1);
 optintphas1=(I-1)*0.1*pi;
 %S0002(shift2(1)+1:shift2(1)+dim,shift2(2)+1:shift2(2)+dim)=D0002*exp(i*optintphas);
 %S0003(shift1(1)+1:shift1(1)+dim,shift1(2)+1:shift1(2)+dim)=D0003*exp(-i*optintphas);
 
 S0004(shift1(1)+1:shift1(1)+dim,shift1(2)+1:shift1(2)+dim)=D0004*exp(i*optintphas1);
 S0005(shift2(1)+1:shift2(1)+dim,shift2(2)+1:shift2(2)+dim)=D0005*exp(-i*optintphas1);
 
 S1=S0001+S0004+S0005;
 I1=abs(ifft2(ifftshift(S1)));
 
 figure('Name', ['BJ m', num2str(2)], 'NumberTitle', 'off');
            imagesc(I1);
            axis equal tight;
            
 name=['SIm1_',num2str(orientation(r),'%05.3d')];
 mkdir([dir_SI,'\',name]);
 imwrite(uint16(I1),[dir_SI,'\',name,'\',num2str(z),'.tif']);

    %% find initial phase for m2 term 
 for t=1:1:20
 ip=(t-1)*0.1*pi;
 S0002(shift3(1)+1:shift3(1)+dim,shift3(2)+1:shift3(2)+dim)=D0002*exp(i*ip);
 S0003(shift4(1)+1:shift4(1)+dim,shift4(2)+1:shift4(2)+dim)=D0003*exp(-i*ip);
 %S0002(shift2(1)+1:shift2(1)+dim,shift2(2)+1:shift2(2)+dim)=D0002*exp(i*ip);
 %S0003(shift1(1)+1:shift1(1)+dim,shift1(2)+1:shift1(2)+dim)=D0003*exp(-i*ip);
 S1=S0001+S0002+S0003;
 I1=abs(ifft2(ifftshift(S1))); 
 intphas2(t)=sum(sum(I1.*WF2x));
 figure(1)
 imagesc(abs(S1))
 end
 [M,I]=max(intphas2);
 optintphas2=(I-1)*0.1*pi;
 %S0002(shift2(1)+1:shift2(1)+dim,shift2(2)+1:shift2(2)+dim)=D0002*exp(i*optintphas);
 %S0003(shift1(1)+1:shift1(1)+dim,shift1(2)+1:shift1(2)+dim)=D0003*exp(-i*optintphas);
 
 S0002(shift3(1)+1:shift3(1)+dim,shift3(2)+1:shift3(2)+dim)=D0002*exp(i*optintphas2);
 S0003(shift4(1)+1:shift4(1)+dim,shift4(2)+1:shift4(2)+dim)=D0003*exp(-i*optintphas2);%% reconstruct SIM image
 S1=S0001+S0002+S0003;
 
 I1=abs(ifft2(ifftshift(S1)));
 name=['SIm2_',num2str(orientation(r),'%05.3d')];
 mkdir([dir_SI,'/',name]);
 imwrite(uint16(I1),[dir_SI,'/',name,'/',num2str(z),'.tif']);
 %imwrite(uint16(I1),[dir_SI,'\','I',num2str(orientation(r),'%05.3d'),'.tif'])

  toc;
  
    %% save frequency domains 
 eval(['F',num2str(orientation(r),'%05.3d'),'1','=S0001',';']); 
 eval(['F',num2str(orientation(r),'%05.3d'),'2','=S0002',';']); 
 eval(['F',num2str(orientation(r),'%05.3d'),'3','=S0003',';']);
 eval(['F',num2str(orientation(r),'%05.3d'),'4','=S0004',';']); 
 eval(['F',num2str(orientation(r),'%05.3d'),'5','=S0005',';']);

    %% test if we need to create new SIMpsf or Apodization
 %test=[eval(['shift1','~=','shift',num2str(orientation(r),'%05.3d'),'1'])
 %      eval(['shift2','~=','shift',num2str(orientation(r),'%05.3d'),'2'])
 %      eval(['shift3','~=','shift',num2str(orientation(r),'%05.3d'),'3'])
 %      eval(['shift4','~=','shift',num2str(orientation(r),'%05.3d'),'4'])];

 %test=vertcat(test,test)
end

%% test if we need to create new SIMpsf or Apodization
    if z==1
        test=1
    elseif size(kp,1)~=0
        test=0
    elseif kpcal((z-1)*size(orientation,2)+1:(z-1)*size(orientation,2)+size(orientation,2),1:4)~=kpcal((z-2)*size(orientation,2)+1:(z-2)*size(orientation,2)+size(orientation,2),1:4)
    %elseif kpcal((z-1)*size(orientation,2)+2,1:4)~=kpcal((z-2)*size(orientation,2)+2,1:4)
    %elseif kpcal((z-1)*size(orientation,2)+3,1:4)~=kpcal((z-2)*size(orientation,2)+3,1:4)
        test=1
    else
        test=0
    end

    %% Create SIM psf

if test==1

    k=2*pi/(lambda/n);
    py=px;
    
for r=1:size(orientation,2)
 
    %% create pattern and rawdata
kpx=floor(dim/2)+1-kpcal(r,3);
kpy=floor(dim/2)+1-kpcal(r,4);
alpha=atan(kpy/kpx)
kx=k*cos(alpha)
ky=k*sin(alpha)
period=px*dim/sqrt(kpx^2+kpy^2)
theta=asin(488/(2*1.33*period))
[vx, vy] = meshgrid(1:dim,1:dim);
P1=(3+2*cos(2*(kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,1)*2*pi)+4*cos((kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,1)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
P2=(3+2*cos(2*(kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,2)*2*pi)+4*cos((kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,2)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
P3=(3+2*cos(2*(kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,3)*2*pi)+4*cos((kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,3)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
P4=(3+2*cos(2*(kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,4)*2*pi)+4*cos((kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,4)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
P5=(3+2*cos(2*(kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,5)*2*pi)+4*cos((kx*(vx-(floor(dim/2)+1))*px+ky*(vy-(floor(dim/2)+1))*py)*sin(theta)+phase(1,5)*pi)*cos(k*0*0*(cos(theta)-1)))/9;

%for x=1:dim;
%    for y=1:dim;
%P1(y,x)=(3+2*cos(2*(kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,1)*2*pi)+4*cos((kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,1)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
%P2(y,x)=(3+2*cos(2*(kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,2)*2*pi)+4*cos((kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,2)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
%P3(y,x)=(3+2*cos(2*(kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,3)*2*pi)+4*cos((kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,3)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
%P4(y,x)=(3+2*cos(2*(kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,4)*2*pi)+4*cos((kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,4)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
%P5(y,x)=(3+2*cos(2*(kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,5)*2*pi)+4*cos((kx*(x-(floor(dim/2)+1))*px+ky*(y-(floor(dim/2)+1))*py)*sin(theta)+phase(1,5)*pi)*cos(k*0*0*(cos(theta)-1)))/9;
%    end
%end

SP=zeros(dim,dim);
SP(floor(dim/2)+1,floor(dim/2)+1)=10000;

SP1=conv2(SP.*P1,psf,'same');
SP2=conv2(SP.*P2,psf,'same');
SP3=conv2(SP.*P3,psf,'same');
SP4=conv2(SP.*P4,psf,'same');
SP5=conv2(SP.*P5,psf,'same');

%% reconstruct

SP1=fftshift(fft2(SP1));
SP2=fftshift(fft2(SP2));
SP3=fftshift(fft2(SP3));
SP4=fftshift(fft2(SP4));
SP5=fftshift(fft2(SP5));

%DSP0001=zeros(dim,dim);
%DSP0002=zeros(dim,dim);
%DSP0003=zeros(dim,dim);
%DSP0004=zeros(dim,dim);
%DSP0005=zeros(dim,dim);

SP1=reshape(SP1,s);
SP2=reshape(SP2,s);
SP3=reshape(SP3,s);
SP4=reshape(SP4,s);
SP5=reshape(SP5,s);
B000=[SP1,SP2,SP3,SP4,SP5];    
%DSP000=linsolve(A000,B000.');
DSP000=mldivide(A000,B000.');
%    for x=1:dim;
%        for y= 1:dim;           
%            B000=[SP1(x,y) SP2(x,y) SP3(x,y) SP4(x,y) SP5(x,y)];
%            DSP000=linsolve(A000,B000.');          
%            DSP0001(x,y)=DSP000(1);
%            DSP0002(x,y)=DSP000(2);
%            DSP0003(x,y)=DSP000(3);
%            DSP0004(x,y)=DSP000(4);
%            DSP0005(x,y)=DSP000(5);
%        end
%    end

 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'1','=zeros(2*dim,2*dim)',';'])
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'2','=zeros(2*dim,2*dim)',';'])
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'3','=zeros(2*dim,2*dim)',';'])
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'4','=zeros(2*dim,2*dim)',';'])
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'5','=zeros(2*dim,2*dim)',';'])

   DSP0001=DSP000(1,:);
   DSP0002=DSP000(2,:);
   DSP0003=DSP000(3,:);
   DSP0004=DSP000(4,:);
   DSP0005=DSP000(5,:);
DSP0001=reshape(DSP0001,[dim,dim]);
DSP0002=reshape(DSP0002,[dim,dim]);
DSP0003=reshape(DSP0003,[dim,dim]);
DSP0004=reshape(DSP0004,[dim,dim]);
DSP0005=reshape(DSP0005,[dim,dim]);
SP0001=zeros(2*dim,2*dim);
SP0002=zeros(2*dim,2*dim);
SP0003=zeros(2*dim,2*dim);
SP0004=zeros(2*dim,2*dim);
SP0005=zeros(2*dim,2*dim);
SP0001(floor(dim/2)+1:floor(dim/2)+dim,floor(dim/2)+1:floor(dim/2)+dim)=DSP0001;
SPWF=abs(ifft2(ifftshift(SP0001))); %wide-field image for comparison%
     %% find initial phase for m1 term

 for t=1:1:20
 ip=(t-1)*0.1*pi;
 eval(['SP0004(shift',num2str(orientation(r),'%05.3d'),'1(1)+1:shift',num2str(orientation(r),'%05.3d'),'1(1)+dim,shift',num2str(orientation(r),'%05.3d'),'1(2)+1:shift',num2str(orientation(r),'%05.3d'),'1(2)+dim)','=DSP0004*exp(i*ip)',';'])
 eval(['SP0005(shift',num2str(orientation(r),'%05.3d'),'2(1)+1:shift',num2str(orientation(r),'%05.3d'),'2(1)+dim,shift',num2str(orientation(r),'%05.3d'),'2(2)+1:shift',num2str(orientation(r),'%05.3d'),'2(2)+dim)','=DSP0005*exp(i*-ip)',';'])
 
 SPm1=SP0001+SP0004+SP0005;
 ISPm1=abs(ifft2(ifftshift(SPm1))); 
 intphasSPm1(t)=sum(sum(ISPm1.*SPWF));
 figure(1)
 imagesc(abs(SPm1))
 end
 [M,I]=max(intphasSPm1);
 optintphasSPm1=(I-1)*0.1*pi;
 eval(['SP0004(shift',num2str(orientation(r),'%05.3d'),'1(1)+1:shift',num2str(orientation(r),'%05.3d'),'1(1)+dim,shift',num2str(orientation(r),'%05.3d'),'1(2)+1:shift',num2str(orientation(r),'%05.3d'),'1(2)+dim)','=DSP0004*exp(i*optintphasSPm1)',';'])
 eval(['SP0005(shift',num2str(orientation(r),'%05.3d'),'2(1)+1:shift',num2str(orientation(r),'%05.3d'),'2(1)+dim,shift',num2str(orientation(r),'%05.3d'),'2(2)+1:shift',num2str(orientation(r),'%05.3d'),'2(2)+dim)','=DSP0005*exp(-i*optintphasSPm1)',';'])
 %SP0004(shift1(1)+1:shift1(1)+dim,shift1(2)+1:shift1(2)+dim)=DSP0004*exp(i*optintphasSPm1);
 %SP0005(shift2(1)+1:shift2(1)+dim,shift2(2)+1:shift2(2)+dim)=DSP0005*exp(-i*optintphasSPm1);
 
    %% find initial phase for m2 term 
 for t=1:1:20
 ip=(t-1)*0.1*pi;
 eval(['SP0002(shift',num2str(orientation(r),'%05.3d'),'3(1)+1:shift',num2str(orientation(r),'%05.3d'),'3(1)+dim,shift',num2str(orientation(r),'%05.3d'),'3(2)+1:shift',num2str(orientation(r),'%05.3d'),'3(2)+dim)','=DSP0002*exp(i*ip)',';'])
 eval(['SP0003(shift',num2str(orientation(r),'%05.3d'),'4(1)+1:shift',num2str(orientation(r),'%05.3d'),'4(1)+dim,shift',num2str(orientation(r),'%05.3d'),'4(2)+1:shift',num2str(orientation(r),'%05.3d'),'4(2)+dim)','=DSP0003*exp(i*-ip)',';'])
 %SP0002(shift3(1)+1:shift3(1)+dim,shift3(2)+1:shift3(2)+dim)=DSP0002*exp(i*ip);
 %SP0003(shift4(1)+1:shift4(1)+dim,shift4(2)+1:shift4(2)+dim)=DSP0003*exp(-i*ip);
 
 SPm2=SP0001+SP0002+SP0003;
 ISPm2=abs(ifft2(ifftshift(SPm2))); 
 intphasSPm2(t)=sum(sum(ISPm2.*SPWF));
 figure(1)
 imagesc(abs(SPm2))
 end
 [M,I]=max(intphasSPm2);
 optintphasSPm2=(I-1)*0.1*pi;
 eval(['SP0002(shift',num2str(orientation(r),'%05.3d'),'3(1)+1:shift',num2str(orientation(r),'%05.3d'),'3(1)+dim,shift',num2str(orientation(r),'%05.3d'),'3(2)+1:shift',num2str(orientation(r),'%05.3d'),'3(2)+dim)','=DSP0002*exp(i*optintphasSPm2)',';'])
 eval(['SP0003(shift',num2str(orientation(r),'%05.3d'),'4(1)+1:shift',num2str(orientation(r),'%05.3d'),'4(1)+dim,shift',num2str(orientation(r),'%05.3d'),'4(2)+1:shift',num2str(orientation(r),'%05.3d'),'4(2)+dim)','=DSP0003*exp(i*-optintphasSPm2)',';'])
 %SP0002(shift3(1)+1:shift3(1)+dim,shift3(2)+1:shift3(2)+dim)=DSP0002*exp(i*optintphasSPm2);
 %SP0003(shift4(1)+1:shift4(1)+dim,shift4(2)+1:shift4(2)+dim)=DSP0003*exp(-i*optintphasSPm2);
    %% save frequency domains 
 eval(['FSP',num2str(orientation(r),'%05.3d'),'1','=SP0001',';']); 
 eval(['FSP',num2str(orientation(r),'%05.3d'),'2','=SP0002',';']); 
 eval(['FSP',num2str(orientation(r),'%05.3d'),'3','=SP0003',';']);
 eval(['FSP',num2str(orientation(r),'%05.3d'),'4','=SP0004',';']); 
 eval(['FSP',num2str(orientation(r),'%05.3d'),'5','=SP0005',';']);
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'1(floor(dim/2)+1:floor(dim/2)+dim,floor(dim/2)+1:floor(dim/2)+dim)','=SP0001',';']); 
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'4(shift',num2str(orientation(r),'%05.3d'),'1(1)+1:shift',num2str(orientation(r),'%05.3d'),'1(1)+dim,shift',num2str(orientation(r),'%05.3d'),'1(2)+1:shift',num2str(orientation(r),'%05.3d'),'1(2)+dim)','=SP0004',';']); 
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'5(shift',num2str(orientation(r),'%05.3d'),'2(1)+1:shift',num2str(orientation(r),'%05.3d'),'2(1)+dim,shift',num2str(orientation(r),'%05.3d'),'2(2)+1:shift',num2str(orientation(r),'%05.3d'),'2(2)+dim)','=SP0005',';']); 
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'2(shift',num2str(orientation(r),'%05.3d'),'3(1)+1:shift',num2str(orientation(r),'%05.3d'),'3(1)+dim,shift',num2str(orientation(r),'%05.3d'),'3(2)+1:shift',num2str(orientation(r),'%05.3d'),'3(2)+dim)','=SP0002',';']); 
 %eval(['FSIMP',num2str(orientation(r),'%05.3d'),'3(shift',num2str(orientation(r),'%05.3d'),'4(1)+1:shift',num2str(orientation(r),'%05.3d'),'4(1)+dim,shift',num2str(orientation(r),'%05.3d'),'4(2)+1:shift',num2str(orientation(r),'%05.3d'),'4(2)+dim)','=SP0003',';']); 

end

    %% reconstruct SIM point spread function
WFSP=zeros(2*dim,2*dim);
for r=1:size(orientation,2)
WFSP=WFSP+eval(['FSP',num2str(orientation(r),'%05.3d'),'1']);
end
SIMSP=WFSP./size(orientation,2);
for r=1:size(orientation,2)
 SIMSP=SIMSP+eval(['FSP',num2str(orientation(r),'%05.3d'),'2'])+eval(['FSP',num2str(orientation(r),'%05.3d'),'3']) ...
            +eval(['FSP',num2str(orientation(r),'%05.3d'),'4'])+eval(['FSP',num2str(orientation(r),'%05.3d'),'5']);
end

SIMpsf=abs(ifft2(ifftshift(SIMSP)));
mkdir([dir_SI,'\SIMpsf']);
imwrite(uint16(SIMpsf),[dir_SI,'\SIMpsf\',num2str(z),'.tif'])

end
%% Export final image
    summedWF=summedWF./size(orientation,2);
    imwrite(uint16(summedWF),[dir_SI,'\summedWF\',num2str(z),'.tif'])

    %% SIM image with all orientations
WF2x=zeros(2*dim,2*dim);
for r=1:size(orientation,2)
WF2x=WF2x+eval(['F',num2str(orientation(r),'%05.3d'),'1']);
end

SIm1=WF2x./size(orientation,2);
SIm2=WF2x./size(orientation,2);
SIall=WF2x./size(orientation,2);
for r=1:size(orientation,2)
 SIm1=SIm1+(eval(['F',num2str(orientation(r),'%05.3d'),'4'])+eval(['F',num2str(orientation(r),'%05.3d'),'5']))./3;
 SIm2=SIm2+(eval(['F',num2str(orientation(r),'%05.3d'),'2'])+eval(['F',num2str(orientation(r),'%05.3d'),'3']))./3;      
 SIall=SIall+(eval(['F',num2str(orientation(r),'%05.3d'),'2'])+eval(['F',num2str(orientation(r),'%05.3d'),'3']))./3 ...
            +(eval(['F',num2str(orientation(r),'%05.3d'),'4'])+eval(['F',num2str(orientation(r),'%05.3d'),'5']))./3;
end
imagesc(abs(SIall))

%WF2xA=WF2x.*eval(['FA',num2str(orientation(1),'%05.3d'),'1']);
%SIallA=SIall.*FAall;

WF2x=abs(ifft2(ifftshift(WF2x)));
mkdir([dir_SI,'\recWF']);
imwrite(uint16(WF2x),[dir_SI,'\recWF\',num2str(z),'.tif'])
SIm1=abs(ifft2(ifftshift(SIm1)));
mkdir([dir_SI,'\SIm1']);
imwrite(uint16(SIm1),[dir_SI,'\SIm1\',num2str(z),'.tif'])
SIm2=abs(ifft2(ifftshift(SIm2)));
mkdir([dir_SI,'\SIm2'])
imwrite(uint16(SIm2),[dir_SI,'/SIm2/',num2str(z),'.tif'])
SIall=abs(ifft2(ifftshift(SIall)));
mkdir([dir_SI,'\SIall'])
imwrite(uint16(SIall),[dir_SI,'\SIall\',num2str(z),'.tif'])

%WF2xA=abs(ifft2(ifftshift(WF2xA)));
%mkdir([dir_SI,'\recWFA']);
%imwrite(uint16(WF2xA),[dir_SI,'\recWFA\',num2str(z),'.tif'])
%SIallA=abs(ifft2(ifftshift(SIallA)));
%mkdir([dir_SI,'\SIallA'])
%imwrite(uint16(SIallA),[dir_SI,'\SIallA\',num2str(z),'.tif'])

%% second deconvolution
SIMpsf=SIMpsf./sum(sum(SIMpsf));
SIallSecondDec=deconvlucy(SIall,SIMpsf,10);
SIallSecondDec=SIallSecondDec./4;
%SIallSecondDecA=fftshift(fft2(SIallSecondDec));
%SIallSecondDecA=SIallSecondDecA.*FAall;
%SIallSecondDecA=abs(ifft2(ifftshift(SIallSecondDecA)));
mkdir([dir_SI,'\SIallSecondDec']);
imwrite(uint16(SIallSecondDec),[dir_SI,'\SIallSecondDec\',num2str(z),'.tif'])
%mkdir([dir_SI,'\SIallSecondDecA'])
%imwrite(uint16(SIallSecondDecA),[dir_SI,'\SIallSecondDecA\',num2str(z),'.tif'])

%% Save Kp information
fprintf(fileID,'%s\r\n','------------------');
fprintf(fileID,'%s\r\n','------------------');
fprintf(fileID,'%s\r\n',num2str(z));
fprintf(fileID,'%5.2f\n',kpcal);
%fprintf(fileID,'%s\r\n %5.2f\n %5.2f\r\n','kp1=',kpcal(1,1),kpcal(1,2));
%fprintf(fileID,'%s\r\n %5.2f\n %5.2f\r\n','kp2=',kpcal(2,1),kpcal(2,2));
%fprintf(fileID,'%s\r\n %5.2f\n %5.2f\r\n','kp3=',kpcal(3,1),kpcal(3,2));

end

fclose(fileID);
toc