function [MaskSeg, WiNoise] = FourierPN(data, maxiter)

%    This code constitutes an important tool for obtaining the results in
%    the manuscript entitled "Fourier and bi-dimensional Legendre base
%    interpolator applied to cryo-EM to generate a SARS-CoV-2 S protein 3D
%    volumetric object", submitted to The Journal of Medical Systems (JMS):
%    the home of clinical informatics Research by SPRINGER.
%
%    In this code we search for base periodic noise in all frames under the
%    cryo-EM study, Turonová et al. [28], in order to remove it from the 
%    frames under study in our proposal and hence obtain metadata free from 
%    such noise.
%
%    This code is based on the Fourier transform (multiple cycles) and the 
%    Distance Regularized Level Set Evolution Method (DRLSEM). 
%
%    Correspondings Authors:
%    Dr. Jesus Alonso Arriaga Hernandez
%    jesus.arriagahdz@correo.buap.mx;    dr.j.a.arriaga.hernandez@gmail.com
%
%    Dra. Bolivia Teresa Cuevas Otahola
%    bolivia.cuevasotahola@viep.com.mx;          b.cuevas.otahola@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    Frames of the data under study

[A, B, C] = size(data);

%    Identification of elements of interest of the n frames with base noise
%    or persistently present noise, from the first to the last frames (multiple
%    cycles) in order to be compared with the frame where the spike covers 
%    the largest area (filtering by multiple cycles).

CC = floor( (C / 2));

aux0 = data(:,:,1);
aux01 = data(:,:,2);
auxtop = data(:,:,C);
auxtop1 = data(:,:,C-1);
auxmax = data(:,:,CC);

%    Corresponding Fourier transform

ftaux0 = fftshift(fft2(aux0));
ftaux01 = fftshift(fft2(aux01));
ftauxtop = fftshift(fft2(auxtop));
ftauxtop1 = fftshift(fft2(auxtop1));
ftauxmax = fftshift(fft2(auxmax));

%    Construction of the ideal filter to determine and identify the center 
%    of the frequencies, to subsequently remove the low average frequencies
%    in every frame.

[x,y]=meshgrid(-((A/2)-1):(A/2),-((B/2)-1):(B/2));
x0=(A/2);
y0=(B/2);

%    Filter width
sigma=20;          

%    Filter use for periodic noise

filter = exp(-1*((x).^2+(y).^2)/(2*(1.5*sigma)^2));

%    Filtered data

fil_ftaux0 = ftaux0.*filter;
fil_ftaux01 = ftaux01.*filter;
fil_ftauxtop = ftauxtop.*filter;
fil_ftauxtop1 = ftauxtop1.*filter;
fil_ftauxmax = ftauxmax.*filter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Mean noise in all frames

freprom = (fil_ftaux0 + fil_ftaux01 + fil_ftauxtop + fil_ftauxtop1 + fil_ftauxmax)./5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    Frames without mean noise

BaseNoise = abs(ifft2(freprom));

for i = 1 : C
    aux = data(:,:,i);
    ftaux = fftshift(fft2(aux)).*filter;
    invftaux = (ifft2(ftaux));
    InvFTaux = abs(invftaux) - BaseNoise;
    ft_data(:,:,i) = abs(InvFTaux);
end

%    Data without the background noise
FT_data = uint8(ft_data);

%    First output variable of the present function, data without noise
WiNoise = FT_data;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    At this point the function obtains results, where the Spike information
%    can be easily identified.  Subsequently, the code performs the automation
%    of the segmentation process.
%    Distance Regularized Level Set Evolution Method (DRLSEM).

for j = 1 : C
    Img=FT_data(:,:,j);
    frm=0;
    
    % Regularization parameters, length, Phi area.  Elements of the Delta 
    % function in the regularization. Maximum number of iterations "maxiter"
    % and sigma of the Gaussian (different from the filter used in the 
    % background noise.

    timestep=1;  
    mu=0.2;  
    lambda=5; 
    alfa= -3;  
    epsilon=1.5; 
    c0=2;        % this values cannot be larger than 10

    maxiter = maxiter;  %iterations number for the construction of the phi mask 
    sigma=2.0;    % scale parameter in Gaussian kernel

    G=fspecial('gaussian',30,sigma); 
    Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussian convolution
    
    %   Computation of the border indices in the function "g" based on the 
    %   image gradient
    [Ix,Iy]=gradient(Img_smooth);
    f=Ix.^2+Iy.^2;
    %   g represent edge indicator function.
    g=exp(-f);

    %   Construction of the mask with spike data only. This is the base for
    %   the construction of the 3D solid of the Spike using our proposal 
    %   "2D Legendre base interpolator"

    phi = c0*ones(size(Img));

    %   Position of the starting point to search for the segmentation 
    %   regions, to remove all the background noise and delimit only the 
    %   Spike information in a quasi-autonomous process.
    
    if j < 5    
        phi(298:301,373:376)=-c0; 
        
    elseif j < 6
        phi(290:300,360:370)=-c0;
        phi(335:345,325:335)=-c0;
    
    elseif j < 10
        phi(290:300,360:370)=-c0;
        phi(335:345,325:335)=-c0;
        phi(194:200,371:377)=-c0;   
        phi(440:450,305:311)=-c0;
    
    elseif j < 22
        phi(290:300,360:370)=-c0;
        phi(335:345,325:335)=-c0;
        phi(194:200,371:377)=-c0;   
        phi(440:450,305:311)=-c0;    
        phi(114:120,500:506)=-c0;   
        phi(194:200,250:256)=-c0;   
        phi(430:435,420:425)=-c0;   
        phi(205:210,460:470)=-c0;
    
    elseif j < 28
        phi(360:540,400:415)=-c0;
        phi(290:300,325:335)=-c0;
        phi(210:220,455:465)=-c0;
        phi(134:139,350:355)=-c0;

    elseif j < 32
       phi(420:470,380:390)=-c0;
       phi(160:215,370:380)=-c0;

    elseif j < 35
       phi(440:450,390:400)=-c0;
       phi(200:225,360:390)=-c0;

    elseif j < C
       phi(208:217,364:373)=-c0;
    end

    [vx, vy]=gradient(g);
    
    for k=1:maxiter
        %   Use of phi as a Newmann border condition, described
        %   in each point different to zero in phi

        [nrow,ncol] = size(phi);
        phi([1 nrow],[1 ncol]) = phi([3 nrow-2],[3 ncol-2]);
        phi([1 nrow],2:end-1) = phi([3 nrow-2],2:end-1);
        phi(2:end-1,[1 ncol]) = phi(2:end-1,[3 ncol-2]);
        
        %   Computation of the regularization elements considering phi
        %   as the level set function that minimizes the regularization 
        %   functional.

        [phi_x,phi_y]=gradient(phi);
        ss=sqrt(phi_x.^2 + phi_y.^2);
        aa=(ss>=0) & (ss<=1);
        bb=(ss>1);
        ps=aa.*sin(2*pi*ss)/(2*pi)+bb.*(ss-1);  
        dps=((ps~=0).*ps+(ps==0))./((ss~=0).*ss+(ss==0));  

        [nxx,junk]=gradient(dps.*phi_x - phi_x);
        [junk,nyy]=gradient(dps.*phi_y - phi_y);
        ff1 = nxx+nyy;
        distRegTerm = ff1 + 4*del2(phi);

        %   Computation in differential form of the are where the 
        %   regularization functional is minimized in terms of phi

        diracphi=(1/2/epsilon)*(1+cos(pi*phi/epsilon));
        b = (phi<=epsilon) & (phi>=-epsilon);
        diracPhi = diracphi.*b;

        areaTerm=diracPhi.*g;
    
        %   Computation of the regularization elements that properly 
        %   approximate the distance regularized level set evolution in
        %   terms of phi
        [phi_x,phi_y]=gradient(phi);
        s=sqrt(phi_x.^2 + phi_y.^2);

        %   add a small positive number to avoid division by zero
        Nx=phi_x./(s+1e-10); 
        Ny=phi_y./(s+1e-10);

        [nxx,junk]=gradient(Nx);
        [junk,nyy]=gradient(Ny);
        ff2 = nxx+nyy;
        edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*ff2;
    
        %   Definitive computation of the distance regularized level set 
        %   evolution according to phi
        phi=phi + timestep*(mu/timestep*distRegTerm + lambda*edgeTerm + alfa*areaTerm);
    
        %   show result in every 50 iteration regarding "maxiter"
        if mod(k,50)==1
            frm=frm+1;
            h=figure(5);
            set(gcf,'color','w');
            subplot(1,2,1);
            II=Img;
            II(:,:,2)=Img;II(:,:,3)=Img;

            imshow(II); axis off; axis equal; hold on;  
            q=contour(phi, [0,0], 'r');
            msg=['Contour frame number=' num2str(j)];  title(msg);
           
            subplot(1,2,2);
            mesh(-phi); 
            hold on;  contour(phi, [0,0], 'r','LineWidth',2);
            view([180-30 65]);
            msg=['Iteration number=' num2str(k)];
            title(msg);
            %pause(0.1)
        
            %   Instructions to perform the the animation of the mask 
            %   construction
        
            frame = getframe(h);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            %Write to the GIF File
            if frm == 1
                imwrite(imind,cm,'C:\Users\jesus\Music/femur.gif','gif', 'Loopcount',inf);
            else        
                imwrite(imind,cm,'C:\Users\jesus\Music/femur.gif','gif','WriteMode','append');
            end
        end
    end
    close(h)
    
    Phi = uint8(phi);
    PHi = (Phi/max(max(Phi)));
    PHI = imcomplement(PHi);
    mask = PHI-min(min(PHI)); 
    spikeSeg(:,:,j) = mask.*Img;
end

%    Second output variable of the present function, semiautomatic masks

%    Final result of the function, segmentation masks.
%    Information of the Spike of interest only.
MaskSeg=spikeSeg;

end


