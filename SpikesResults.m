
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This code constitutes an important tool for obtaining the results in
%    the manuscript entitled "Fourier-Legendre cryo-EM to generate 
%    SARS-CoV-2 S protein volumetric object from 2D Legendre bases 
%    Interpolator", submitted to The JOSA (Journal of the Optical Society
%     of America) A .
%
%    Recalling that the ROI are built from 3D tomographic by cryo electron 
%    microscopy (cryo-EM) image data from the research work "In situ 
%    structural analysis of SARS-CoV-2 spike reveals flexibility mediated
%    by three hinges" with DOI: 10.1126/science.abd5223 [28].
%
%    From such work, we selected several frames (video image captures from
%    the tomographic image) and for each frame we obtained the phase of the
%    observed spikes or Protein S from the SARS-CoV-2 cells, segmented by 
%    multiples ciclos de filtrado por Fourier transform Method. Subsequently,
%    we integrate the phase of the different Spikes frames de una misma 
%    SARS-CoV-2 VP (virak particle) in a 3D matrix, where we obtain the 3D 
%    volumetric object (3D solid) of the SARS-CoV-2 virion particle or VP.
%
%    Correspondings Authors:
%    Dr. Jesus Alonso Arriaga Hernandez
%    jesus.arriagahdz@correo.buap.mx;    dr.j.a.arriaga.hernandez@gmail.com
%    Dra. Bolivia Teresa Cuevas Otahola
%    bolivia.cuevasotahola@viep.com.mx;          b.cuevas.otahola@gmail.com
%
%    This algorithm contains a simple routine to call the main files 
%    required to obtain the results as well as to apply our main function
%    "FourierPN.m".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    In a first step, we start by reading and analysing the frames in the 
%    video "abd5223s1.mov" from Cryo-EM obtained by Turonová et al. [28], 
%    which we included as supplementary material labelled as Video-1.  
%    Subsequently, we apply the Fourier transform method by multiple cycles
%    and perform a segmentation with a simple rectangular region.

%    Reading the video in the storage directory in the PC
Video = VideoReader('C:\Your PC\Spikes\Videos\Visualization 1.mov');

%    Previous instruction interpreting the video as an structure where 
%    the frames sizes and number of frames are obtained.
N = Video.NumFrames;
A = Video.Height;
B = Video.Width; 

%    Auxiliary variable "frames" independently sotres the frames for their
%    subsequent analysis
frames = zeros(A,B,N);

for n = 1:N
    hasFrame(Video);   
    aux = double(rgb2gray(readFrame(Video)));
    frames(:,:,n) = aux;
end

%    T is a number increasing the resolution in Matlab, however, it causes a
%    sharpness loss for T > 5. Hence, we choose T=3 considering the computing
%    power of the user.
T = 3;
Frames = imresize(frames, T);
Frames0 = uint8(Frames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    First segmentation section, where we determine the ROI in the frame
%    in such manner that it has a suitable size for the rest of the frames
%    and in order to have data from a single spike.

%    For Spike 1
Spike1 = Frames(:,:,50:97);
[A, B, C1] = size(Spike1);

%    ROI segmentation mask
MaxSpike1 = Frames0(:,:,72);

%    Basic segmentation of the mask with real size in comparison with the 
%    video, for Spike 1
%Mask1 = imcrop(MaxSpike1, [782.5 804.5 84 122]);
%figure; imshow(Mask1)

%    ROI application in all frames of Spike 1
for k = 1:C1
    Aux = Spike1(:,:,k);
    roiSpike1 = imcrop(Aux, [782.5 804.5 84 122]);
    ROISpike1(:,:,k) = roiSpike1;
end
ROISpike10 = uint8(ROISpike1);

%    For Spike 2
Spike2 = Frames(:,:,81:125);
[A, B, C2] = size(Spike2);

%    ROI segmentation mask 
MaxSpike2 = Frames0(:,:,102);

%    Basic segmentation of the mask with real size in comparison with the 
%    video, for Spike 2
%Mask2 = imcrop(MaxSpike2, [1301.5 849.5 103 158]);
%figure; imshow(Mask2)

%    ROI application in all frames of Spike 2
for k = 1:C2    
    Aux = Spike2(:,:,k);
    roiSpike2 = imcrop(Aux, [1301.5 849.5 103 158]);
    ROISpike2(:,:,k) = roiSpike2;
end
ROISpike20 = uint8(ROISpike2);

%    t is a number increasing the resolution in Matlab. We choose T=4
t = 6;
ROISpike1 = imresize(ROISpike10, t);
ROISpike2 = imresize(ROISpike20, t);

%    Tones inversion. SARS-CoV-2 virion data in bright tones

Spike1 = imcomplement(ROISpike1);
Spike2 = imcomplement(ROISpike2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    In this part we focus on removing the Spikes noise. To this aim, we 
%    apply and construct the mask to perfectly delimit the Spike information
%    in the largest area frame (BigSpike). We use the Fourier 2D transform 
%    in the "FourierPN" function. Such filter also reduces the periodic 
%    noise and induces a certain continuity in the image data, removing the
%    pixelating errors.

bspike1 = Spike1(:,:,23);
bspike2 = Spike2(:,:,25);

%    Reading the BigSpike masks 
load('C:\Your PC\Spikes\Programs\BMS1.mat');
bigSpike1 = uint8(BMS1);
load('C:\Your PC\Spikes\Programs\BMS2.mat');
bigSpike2 = uint8(BMS2);

%    Application of the BigSpike mask to all frames for each Spike

for k = 1:C1    
    aux = Spike1(:,:,k);
    Aux = aux.*bigSpike1;
    BigSpike1(:,:,k) = Aux;
end

for k = 1:C2    
    aux = Spike2(:,:,k);
    Aux = aux.*bigSpike2;
    BigSpike2(:,:,k) = Aux;
end

%    Construction of the squared matrices.

[A1, B1] = size(bigSpike1);
Aux0 = zeros(A1);
aux0 = B1/2;
for k = 1:C1    
    aux = BigSpike1(:,:,k);
    AA = aux;
    for i = 1 : B1    
        Aux0(:, ((A1/2)-1-aux0+i)) = AA(:,i);
    end
    sqBigSpike1(:,:,k) = Aux0;
end
SQsqBigSpike1 = uint8(sqBigSpike1);

[A2, B2] = size(bigSpike2);
Aux0 = zeros(A2);
aux0 = B2/2;
for k = 1:C2    
    aux = BigSpike2(:,:,k);
    AA = aux;
    for i = 1 : B2    
        Aux0(:, ((A2/2)-1-aux0+i)) = AA(:,i);
    end
    sqBigSpike2(:,:,k) = Aux0;
end
sqsqBigSpike2 = uint8(sqBigSpike2);
for i = 1 : C2
    SQsqBigSpike2(:,:,i) = imrotate(sqsqBigSpike2(:,:,i),180);
end

%    Periodic noise removal and characterization of the base noise, where
%    the base noise is removed from all the frames in order to subtract it 
%    to obtain only the Spike information by means of segmentation masks 
%    in a semi-automatic process. 
%    We use the FourierPN function, based on the analysis by means of the 
%    2D Fourier Transform (FFT2) and simple filtering using pass-bands.
%    This function has two input and two output parameters 
%    {data, maxiter con [MaskSeg, WiNoise] = FourierPN(data, maxiter)}.
%    Where MaskSeg is the binary segmentation mask and WiNoise the data 
%    analyzed using Fourier by means of our proposal for mean periodic noise
%    removal, "maxiter" is the number of iterations, where we recommend 
%    setting maxiter=100 to start and to obtain optimal results for the  
%    segmentation mask we use maxiter=100. Considering maxiter=800 or more, 
%    we optimize the results of the segmentation masks, with larger execution
%    maps and requiring more computing power. Hence, such value should be 
%    set to a value between 100 and 500 in the following lines.

load('C:\Your PC\Spikes\GitHub Files\PerNoiSpike1.mat');
load('C:\Your PC\Spikes\GitHub Files\PerNoiSpike2.mat');
%     Reading automatic masks for 500 iterations
load('C:\Your PC\Spikes\GitHub Files\MaskSegSpike1.mat')
load('C:\Your PC\Spikes\GitHub Files\MaskSegSpike2.mat')

%    Segmentation masks and mean noise for Spike 1, where the 
%    variable MaskSegSpike1 are the masks of the semiautomatic process
%    for Spike 1
%[MaskSegSpike1, PerNoiSpike1] = FourierPN(SQsqBigSpike1, 100);

%    Segmentation masks and mean noise for Spike 2, where the
%    variable MaskSegSpike2 are the masks of the semiautomatic process
%    for Spike 2
%[MaskSegSpike2, PerNoiSpike2] = FourierPN(SQsqBigSpike2, 100);

%    Mask data for Spike 1 as Double variables
Mask1 = double(MaskSegSpike1);

%    Mask data for Spike 2 as Double variables
Mask2 = double(MaskSegSpike2);

%    In the case of using the binary masks manually built, it is necessary
%    to use for each frame the following lines. The required files can be
%    found in the GitHub repository, which are found in files for manual
%    and automatic masks (MaskSegSpike)
%
%    Mask1
read_mask = dir('C:\Your PC\Spikes\Masks\Masks Spike 1\*.mat'); 

%    Mask2
%read_mask = dir('C:\Your PC\Spikes\Masks\Masks Spike 2\*.mat'); 

J = length(read_mask);

for k = 1:J
    %    Instruction for the Spike 1 frames
    mask0 = strcat('C:\Your PC\Spikes\Masks\Masks Spike 1\BiMaskSpike1_', num2str(k), '.mat');
    
    %    Instruction for the Spike 2 frames
    %mask0 = strcat('C:\Your PC\Spikes\Masks\Masks Spike 2\BiMaskSpike2_', num2str(k), '.mat');
    
    mask = load(mask0);
    masks(:,:,k) = cell2mat(struct2cell(mask));
end
masks0 = double(masks);
Masks = uint8(masks0);
[a, b, J] = size(Masks);

%    Frames filtered for the masks, data only from
for k = 1 : J
    %    Spike1
    Spike1_data(:,:,k) = Masks(:,:,k).*PerNoiSpike1(:,:,k);
    %Spike1_data(:,:,k) = MaskSegSpike1(:,:,k).*PerNoiSpike1(:,:,k);
    %MetaData1(:,:,k) = Masks(:,:,k).*SQsqBigSpike1(:,:,k);
    %MetaData1(:,:,k) = MaskSegSpike1(:,:,k).*SQsqBigSpike1(:,:,k);
    
    %    Spike2   
    %Spike2_data(:,:,k) = Masks(:,:,k).*PerNoiSpike2(:,:,k);
    %Spike2_data(:,:,k) = MaskSegSpike2(:,:,k).*PerNoiSpike2(:,:,k);
    %MetaData2(:,:,k) = Masks(:,:,k).*SQsqBigSpike2(:,:,k);
    %MetaData2(:,:,k) = MaskSegSpike2(:,:,k).*SQsqBigSpike2(:,:,k);
    
end

%    Selection of the files according to the analyzed Spike (1 or 2)

% We consider manual or automatic masks for each Spike, according to the
% the choice in the previous lines, as mentioned in the GitHub repository

Data =  Spike1_data;
%Data = Spike2_data;

% Metadata data without filtering is considered with the addition of manual
% or automatic masks for each Spike according to the choice in the previous 
% lines, , as mentioned in the GitHub repository

%Data = MetaData1;
%Data = MetaData2;

% Unfiltered metadata considering only the ROI by the mean BigMask in 
% all frames, without Fourier analysis.  BiigSpike Mask is included in the
% GitHub repository with the variable BMS1 or BMS2  

%Data = SQsqBigSpike1;
%Data = SQsqBigSpike2;

%    The following lines allow to save the data sequence, which is the
%    multiplication of the mask by the data without periodic noise
%    (PerNoiSpike1 or PerNoiSpike2) and the metadata (untreated data
%    segmented at most by a binary mask of simple region, as ROISpike
%    and bigspike; of variable SQsqBigSpike1 or SQsqBigSpike2).
%    We recall that these data subsequently saved are treated with the
%    Python algorithm 3DLegendreS.
%    We bear in mind that the results are compared for both Spikes and
%    for the filtered data and metadata, to finally obtain a total of
%    four 3D volumetric Objects or 3D Solids.  

%    The line to create the variables where the frames are stored and
%    that will be executed in Python with 3DLegendreS

for k= 1 : J
  %    Frames Spike
  
  aux = Data(:,:,k);
  for i = 1 : a
      for j = 1 : b
          if aux(i,j) == 0
              aux1(i,j) = NaN;
          else
              aux1(i,j) = aux(i,j);
        end
      end
  end
  A = aux1;

%   It is important to bear in mind that "FrameSpike1" or "FrameSpike2" 
%   need to be modified according to the Spike under analysis since when
%   saving the files, only the memory space will be stored instead of the
%   data
  eval(sprintf('FrameSpike1_%d = A;', k));

%    Spike 1, these lines only store the data, which should be previously
%    created, which is the reason of the previous line

%    In this line, we indicate the destination directory where the file
%    sequence will be stored, as shown in the example    
%   filename = ['C:\Your PC\Spikes\Frames\FrameSpike1_', int2str(k), '.mat'];
    filename = ['FrameSpike1_', int2str(k), '.mat'];

%    This line saves the files, only in the workspace directory 
    save(filename, 'A','-mat');
end


%    The variables saved in previous lines are imported and opened in Python.
%    We consider the size of FrameSpike1_T (for each T with T = 48
%    for Spike1 and T = 45 for Spike2). To this aim, we generate intermediate
%    frames from the sizes of each frame (MxMxT), intending to obtain a 
%    a matrix MxMxM by means of the 2D Legendre base interpolator algorithm
%    integrated into the function 3DLegendreS in Python. The procedure to 
%    generate the frames from number T to M (M = 738 for Spike1 and T = 954 
%    for Spike2) requires generating intermediate frames until the desired 
%    number M is reached. We recall that the file size is significantly large
%    ~6.5Gb, which reflects also the reading and manipulation procedures.

%    The following lines read the M frames generated for each Spike and 
%    build a 3D matrix with MxMxM-shape. Considering the sizes of the Spikes
%    and the SARS-CoV-2 viral particles, finally, a 3D solid is built,
%    observed with different filters to recognize and identify the 
%    RNA and other elements of the S Protein.

%    Images for Cell 2 fit
M_Files = dir('C:\Your PC\Spikes\GitHub Files\3DModelsPython\Spike_1\Filter Data Manual Masks\*.mat'); 
%M_Files = dir('C:\Your PC\Spikes\GitHub Files\3DModelsPython\Spike_2\Filter Data Manual Masks\*.mat'); 
M = length(M_Files);

for m = 1:M
    
    %    Instruction for opening the 3DLegendreS fit files of Spike
    file = strcat('C:\Your PC\Spikes\GitHub Files\3DModelsPython\Spike_1\Filter Data Manual Masks\FrameSpike1_', num2str(m), '.mat');        
    %file = strcat('C:\Your PC\Spikes\GitHub Files\3DModelsPython\Spike_2\Filter Data Manual Masks\FrameSpike1_', num2str(m), '.mat');        
    
    EstrAux = load(file);
    fitimage = EstrAux.mat;
    %    Identification of the images generated by 3DLegendreS
    aux_0 = fitimage;
    FitImage(:,:,m) = aux_0;

    %    Images with edges and zeros in NaN variables
    Aux = FitImage(:,:,m);
    for i=1:a
        for j=1:b
            if Aux(i,j) == 0
                NanImages(i,j,m) = NaN;
            else
                NanImages(i,j,m) = Aux(i,j);
            end
        end
    end
end

%    Instruction to see the 3D model of the SARS-CoV-2 cell.
%    Variables contained in the 3D Solids folder of the GitHub repository
volumeViewer(NanImages);
%volumeViewer(Volum3DObject_S1);


%    Measurements of the Spike obtained from the video metadata "In situ 
%    structural analysis of SARS-CoV-2 spike reveals flexibility mediated
%    by three hinges". To this aim we consider the Figs. 1 and 2  

% considering as nucleus diameter of 76.7nm for SARS-CoV-2 VP 1 (spike 1)
% considering as nucleus diameter of 81.1nm for SARS-CoV-2 VP 2 (spike 2)
% Considering a length of Spike 1 of 1 25.5nm for spike 1
% Considering a length of Spike 2 of 29.5nm for spike 2

%   Rescaling the measurements, 59nm for Spike 1 or 2 for a cubic matrix, 
%   which is plotted in several of its isophotes without integrating it into
%   the 3D solid, without the unnecesary data.
%   Spike 1
XX = (A1/(A1-(0.01*A1)))*25;

%    Spike 2
%XX = (A2/(A2-(0.1*A2)  ))*30;

[X, Y]=meshgrid(0:(XX/a):XX);
X=imresize(X,[M M]);    Y=imresize(Y,[M M]);

figure;
hold on
 for ii = 1 : M
     zz = (J) - (ii);
     [~,h] = contour(X,Y,NanImages(:,:,ii),30);% plot contour at the bottom
     h.ContourZLevel = zz;
     
 end
hold off
view(3);
xlabel('X [nm]'); ylabel('Y [nm]'); zlabel('Z [nm]', 'Rotation',0)





