clear all;close all;clc;
tic
% Start a simulation of the reconstruction

%%
% General options


% Testing Options
results = 'C:\Computational Imaging Simulations\Simulation Results\20 - comparison convol vs nnls - c\';
proj =20730; % Projections
ng = 0; %Gaussian noise in dB SNR
np = 0; %poisson noise enabled
off = 0;

use_poisson_noise = 0;
use_gaussian_noise = 0;
use_offset = 0;

use_quantization = 0;
quantbits = 16;

use_countscaling = 0;
maxcount = 1000000;
% store masks and where
save_masks = 1;
use_existing_masks = 0;
masks_dir = 'C:\Computational Imaging Simulations\masks\cs-comp 16II\';%%'C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7';
use_hadamard = 0;

% Noise given as normally distributed with fixed snr
use_randn_noise = 0;
noise_snr = 0;


%gaussian_snr = 10; % in dB

% Specify the measurement dimensions and target ans kernel and indices
%Set Dimensions
height = 1080;%48;%for Bs 16   768;
width = 1920;%64;% for Bs 16  1024;
%nmasks = 400000;

%Percentage of total measurement energy as constant offset

offset_percentage = 0;

%Kernel properties
ke_probability = 0.0004;
ke_blockSize = 5;
N = height*width;
%A= zeros(proj,N);
% use a callback method -> write the corresponding things withing callback
% itself
use_callback = 1;
callback_func = @callback;

%% 

%Get Target ready
%target = mat2gray(imread('USAF-3.tif'));
%target = logical(target);
%target = phantom(1024);
%target = imresize(target, [height, width], 'Method', 'nearest');

%Create Kernel and indices
kernel = generate_kernel(0.5, 10, height, width, 'stream');
imagesc(kernel)
if use_hadamard
    
    %kernel = generate_kernel_hadamard (256, 4, height, width); 
else 
    
    %kernel=imread('Kernel16.tif');
    %kernel=im2double(kernel);
    %kernel = generate_kernel(ke_probability, ke_blockSize, height, width, 'circles');
end    




%%
% Reconstruction options

opt = struct;
opt.tol = 1e-18; %tolerance
opt.max_time = 10 * 60; %maximum calculation time in seconds
opt.max_iter = 5000; %maximimum iterations

% For flexbox
opt.method ='flexbox';%'flexbox';
opt.flexbox_internal_iter = 500000;  %How many iterations before checkpoint -> callback
opt.flexbox_nnls = 1; %use nnls or nor
opt.flexbox_l1grad = 0;%; %1e5; %factor for ||grad(f)||_1 regularization (this is TV-1-Norm)
opt.flexbox_l2grad = 0; %factor for ||grad(f)||_2^2 regularization (TV-2-Norm)
opt.flexbox_wname = 'none'; %use some wavelet or 'none'

%opt.flexbox_wname = 'db45';
%opt.flexbox_wfactor = 1; %weighting in minimization

%Alternative Method: Least Squarest
% opt.method = 'lsqr';

%%
% Rest of the code -> so the simulation
for o = off% vari offset
    offset_percentage = o;
    for p = np
        use_poisson_noise = p;
        for g = ng%700000
            gaussian_snr = g;
            for i = proj%700000
                nmasks = i;
                sampling = i/height * width; 
                indices = randsample(height*width, nmasks);
                if (use_existing_masks)  % load maskset, kernel and variables
                    nmaskstemp = nmasks;
                    load (fullfile(masks_dir,'workspace.mat'),'nmasks','height','width','indices','kernel','y') 
                    X = ['max. number of existing masks ',num2str(nmasks)];
                    disp(X);
                    if (nmaskstemp < nmasks)
                        nmasks = nmaskstemp;

                    end
                    y = y(1:nmasks);
                    indices = indices(1:nmasks);
                    %load (fullfile(masks_dir,'workspace.mat')) 
                else

                    if save_masks == 0
                        %Generate y1, which is yield using convolution directly
                        %y = ifft2(fft2(kernel) .* fft2(target));
                        %y = y(indices);
                        
                        y = zeros(nmasks, 1);
                        if exist(masks_dir, 'dir') == 0
                            [~,~,~] = mkdir(masks_dir);
                        end
                        for i = 1:nmasks
                            mask = generate_mask(kernel, [height, width], indices(i));
                            %imwrite(mask, fullfile(masks_dir, [int2str(i),'.png']));
                            y(i) = abs(sum(sum(mask .* target)));

                        end
                    else
                        %Generate y from the masks themselves
                        y = zeros(nmasks, 1);
                        if exist(masks_dir, 'dir') == 0
                            [~,~,~] = mkdir(masks_dir);
                        end
                        for i = 1:nmasks
                            mask = generate_mask(kernel, [height, width], indices(i));
                            imwrite(mask, fullfile(masks_dir, [int2str(i),'.png']));
                            temp=imread([masks_dir num2str(i) '.png']);
                            temp=reshape(temp,height*width,1);
                            
                            A(i,:) = temp;
                            y(i) = abs(sum(sum(mask .* target)));

                        end
                    end

                    save(fullfile(masks_dir,'workspace.mat')); 
                end
                
                if use_offset == 1
                    offset = (offset_percentage/100) * max(y);
                    y = y + offset;
                else
                    offset = 0;
                end

                if use_countscaling == 1 % normalize to max number of counts
                    y_max = max (y);
                    y = y * (maxcount/y_max);
                end
                
              
                if use_randn_noise == 1
                    noise_r = randn(nmasks,1);
                    noise_r = noise_snr .* norm(y) .* (noise_r ./ norm(noise_r));
                else
                    noise_r = zeros(nmasks, 1);
                end

                if use_poisson_noise == 1
                   noise_p = poissrnd(y)-y;
                else
                   noise_p = zeros(nmasks, 1);
                end

                if use_gaussian_noise == 1
                    %noise_g = y + gauss_std*randn -y;
                    noise_g = awgn(y,gaussian_snr,'measured')-y;
                    %n_ref = wgn(y,gaussian_snr,'measured');
                else
                    noise_g = zeros(nmasks, 1);
                end    


                %else
                %    noise = zeros(nmasks, 1);
                %end
                noise_sum = sqrt (noise_r.^2 + noise_p.^2 + noise_g.^2);
                y = y + noise_sum;
                
                if use_quantization == 1  % use quantization
                    y  = xquant(y, quantbits, maxcount,'mid-rise');
                end
                y = max(y,0);
                
                % Apply the noise and offset and clip negative values
                %y = abs(y + offset + noise);

                %Reconstruct
                if use_callback == 1
                    opt.callback = callback_func;
                end
                
                if use_countscaling == 1 % normalize to max number of counts  %scale y for recon
                    %y_max = max (y);
                    y = y * (y_max/maxcount);
                end
                
                if use_offset == 1 && offset_percentage > 0
                    offset = max(y)/(offset_percentage/100);                 %scale y for recon
                end    
                convolve_flag = 1;
                %[x1, mu] = reconstruct(kernel,A,transpose(A) ,indices, y, convolve_flag, opt);% reconstruction 
                %imagesc(x)
                [x, mu] = reconstruct(kernel,A,transpose(A) ,indices, y, convolve_flag, opt);
%                 convolve_flag = 0;
%                 [x2, mu] = reconstruct(kernel,A,transpose(A) ,indices, y, convolve_flag, opt);% reconstruction 
                imagesc(x)
                %x3 = lsqr(A,y,opt.tol,opt.max_iter);
               % x3 = reshape(x3,height,width);
%                 nnopt  = struct;
%                 nnopt.Tol = opt.tol;
%                 nnopt.Iter = opt.max_iter;
%                 x4=nnls(A,y,nnopt);
%                 x4 = reshape(x4,height,width);
                %imagesc(x)

                %conv_recon(kernel,indices,y);

                %Lets have a look
%                 figure;
%                 subplot(4,1,1);
%                 imagesc(x1);
%                 title('conv');
%                 subplot(4,1,2);
%                 imagesc(x2);
%                 title('lsqr + afun');
%                 subplot(4,1,3);
%                 imagesc(x3);
%                 title('lsqr');

%                 subplot(4,1,4);
%                 imagesc(x4);
%                 title('nnls');
%                 
                filenamesum = [results int2str(nmasks)];
                if use_poisson_noise == 1
                    filenamesum = [filenamesum '_P' int2str(maxcount)];
                end
                if use_gaussian_noise == 1
                    filenamesum = [filenamesum '_G' num2str(gaussian_snr) 'dB'];
                end
                if use_quantization == 1
                    filenamesum = [filenamesum '_' num2str(quantbits) 'bit'];
                end
                if use_offset == 1
                    filenamesum = [filenamesum '_O' num2str(offset_percentage)];
                end
                filenamesum = [filenamesum  '2.tif'];
                %x = imadjust(x);
                imwrite(x, filenamesum);
            end
        end
    end
end
imwrite(kernel,[results 'kernel_' num2str(ke_probability) '_' int2str(ke_blockSize) '.png']);

numMasks=45000;
%
temp=imread(['C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7w\3.png']);
wmask=temp.*255;
imwrite(wmask,['C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7w\3.png']);
%}
state=1;
st=0;

 for i=1:numMasks;  
     temp=imread(['C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7\' int2str(i),'.png']);
     mask=temp.*255;
        %mask=scaleimg(temp,[1080,1920]).*255; %For rescale
        %imwrite(mask,['C:\Users\Cristos\Desktop\Christos\Masks\Cmasks w\' num2str(2*i+1) '.png']);
   
    
    if (state==4)
        
        copyfile(['C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7w\3.png'],['C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7w\' num2str(i+2+st) '.png'])
        st=st+1;
        imwrite(mask,['C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7w\' num2str(i+2+st) '.png']);
        state=2;
    else
                
        imwrite(mask,['C:\Users\Cristos\Desktop\Sven\Masks\Cmasks 7w\' num2str(i+2+st) '.png']);
        state=state+1;
    end
    
  
    %imwrite(wmask,['C:\Users\Cristos\Desktop\Cristos\Masks\Cmasks w\' num2str(2*i+2) '.png']);
    %copyfile(['C:\Users\Cristos\Desktop\Cristos\Masks\Cmasks w\2.png'],['C:\Users\Cristos\Desktop\Christos\Masks\Cmasks w\kernel' num2str(ke_) '.png'])
         
 end

toc
