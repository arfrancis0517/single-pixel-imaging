%% single pixel camera _to learn the sharpness of the mask on the OS and OS effects
clc;clear all;close all;
frequency = 0.350e12;% in THz
c = physconst('LightSpeed'); % speed of light in meter per sec
lambda = (c/frequency)*1;% in meter or ~1 mm
wo = 0.00125*1; % beam waist 1.5mm
k = (2*pi)/lambda;% wavenumber
mask = "NNLS"; % "NNLS" or  "hadamard" for mask generation
Measurements =  15000; % Number of measurements
Blk_Size =8; % Block size
h=zeros;t=zeros; % 
%% Sampling
%provide grid space, sampling period and z for gaussian beam propagation
Lx = 0.050; % grid space in meter
Ly = 0.050;
Nx = 1024;Ny =1024; % no of samples
deltax = Lx/Nx; deltay = Ly/Ny; %sampling period
[x,y,X,Y,FX,FY,r,theta,fx,fy,rf] = sampling(Nx,Ny,deltax,deltay);
%% CREATING OBJECT
n= 200;% number of pixels
w = 20; % number of 1spokes
p_n =412; % number of pixels on both sides for zeropadding
[U] = siemensstar(n,w,p_n);
uf=fftshift(fft2(U));
imagesc(FY,FX,((abs(((uf))))),'CDataMapping','scaled'),axis square,colormap jet
imagesc(X,Y,((abs(((U))))),'CDataMapping','scaled'),axis square,colormap parula
%%  Gaussian beam Parameters
z= 0.012; %  propagation distance in meters0.001
paraxial =1; % 1 for Paraxial & 0 for non paraxial
wox= 0.01;% gaussian beam for tilt beam
w1=0.005;% beam waist radius for illumination beam
Zr = (pi*wox^2)/lambda;wx=wox*sqrt(1+(z/Zr)^2); % rayleigh length & beam waist radius
R= z+(Zr^2/z);% radius of curvature
phi = atan(z/Zr);% Gouy phase
azi = 38*pi/180; ele = 50*pi/180; % Azimuthal & Elevation angle in radians
%% Gaussian beam illumination
gau= exp(-(((x-0.002).^2+(y-0.0002).^2))/wx^2);%.*exp(-1i*phi+1i*k*sin(ele).*sin(azi).*Y-1i*k*sin(ele)*cos(azi).*X);
gau1 = gaussian(w1,lambda,r,k,0.08,paraxial);
imagesc(X,Y,((abs(((gau1))))),'CDataMapping','scaled'),axis square,colormap jet
%% tilted plane wave
Z_hologram = 0.00040;
prop = Propagator(Nx, lambda, Lx, Z_hologram);
pz=IFT2Dc(FT2Dc(gau.*exp(1i*k*sin(ele).*sin(azi).*y+1i*k*sin(ele)*cos(azi).*x)).*prop);
imagesc(abs(pz))
%% SIMULATING HOLOGRAM
U1 = IFT2Dc(FT2Dc(U.*gau1).*prop);
U2 = IFT2Dc(FT2Dc(gau1).*prop);
hologram = abs((U1+pz).*conj(U1+pz));%abs(U+pz).^2;fftshift(pz).*prop
imagesc(abs(hologram./max(hologram(:))))
hologram1 = abs((U2+pz).*conj(U2+pz));
imagesc(abs(hologram1./max(hologram1(:))))
h1_norm = hologram1./(max(hologram1(:)));
h_norm = hologram./max(hologram(:));
figure, imshow(hologram, []);
colormap(gray)
t1 =abs(h1_norm-h_norm);
imagesc(abs(t1./max(t1(:))))
%% Mask generation
[Image,A] = Maskgen(Measurements,Blk_Size,Nx,Ny,mask);
%% illumination from halogen lamp on OS
F_num = 2.53; % Fnumber for Projection lensof masks
%% transfer functions related to propagation
Lens= @(r,k,F)(F^2./(r.^2 + F^2)).*exp(-1i.*(k*(sqrt(r.^2+F^2)-F)));% Lens function
%% Propagation from OS to detector, Angular spectrum method
F=0.200;F1 = 0.150;F2 = 0.025;%% Focal length
Dia = 0.050; Dia1 = 0.081; Dia2 = 0.030; % Diameter of lens in mm
z_prop = 0.0001; % propagation distance to first lens
%% from OS to Detector
if(mask=="hadamard")
    t=zeros;h=zeros;
    for i = 1:1024
        U_objp =  Image{i}.*t1; % mask and the object
        U1p=ifft2(fftshift(fft2(U_objp)).*Propagator(Nx,lambda,Lx,z_prop));
                 U_Fp = ifft2((fft2(Lens(r,k,F).*pupil(r,Dia,Nx,Ny).*U1p)).*Propagator(Nx,Ny,deltax,deltay,lambda,F));
        %          U_F1p = ifft2((fft2(Lens(r,k,F1).*pupil(r,Dia1,Nx,Ny).*U_Fp)).*prop(k,fx,fy,F1));
        %          U_F2p = ifft2((fft2(Lens(r,k,F).*pupil(r,Dia2,Nx,Ny).*U_F1p)).*prop(k,fx,fy,F2));
        t(i)=sum(sum(abs(U_objp(:))));
        %% negatuve
        U_objn =  A{i}.*t1;
        U1n=ifft2(fftshift(fft2(U_objn)).*Propagator(Nx,lambda,Lx,z_prop));
                 U_Fn =ifft2((fft2(Lens(r,k,F).*pupil(r,Dia,Nx,Ny).*U1n)).*Propagator(Nx,Ny,deltax,deltay,lambda,F));
        %          U_F1n = ifft2((fft2(Lens(r,k,F1).*pupil(r,Dia1,Nx,Ny).*U_Fn)).*prop(k,fx,fy,F1));
        %          U_F2n = ifft2((fft2(Lens(r,k,F2).*pupil(r,Dia2,Nx,Ny).*U_F1p)).*prop(k,fx,fy,F2));
        h(i)=sum(sum(abs(U_objn(:))));
        %           imagesc(abs(U1p))
    end
elseif(mask=="NNLS")
    h = zeros; t=A;
    for i = 1:Measurements
        U_obj = Image{i}.*t1;
                 U1=ifft2(fftshift(fft2(U_obj).*Propagator(Nx,lambda,Lx,z_prop)));
%                  U_F = ifft2(fftshift(fft_U(Lens(r,k,F).*pupil(r,Dia,Nx,Ny).*U1)).*Propagator(Nx,Ny,deltax,deltay,lambda,F));
        %          U_F1 =  ifft2((fft2(Lens(r,k,F1).*pupil(r,Dia1,Nx,Ny).*U_F)).*prop(k,fx,fy,F1));
        %          U_F2 = ifft2((fft2(Lens(r,k,F2).*pupil(r,Dia2,Nx,Ny).*U_F1)).*prop(k,fx,fy,F2));
        h(i)=sum(sum(abs(U_obj)));
        imagesc(abs(U_obj))
    end
end
%% Reconstruction
[img_rec,y] = reconstruction(t,h',Nx/Blk_Size,mask);% Blk_Size for hadarmd & (Nx/Blk_Size)for NNLS
plot(y)
imagesc(abs(img_rec./max(img_rec(:)))),axis square, colormap jet
%% Hologram reconstruction 


