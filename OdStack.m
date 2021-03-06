
function S=OdStack(structin,method,plotflag,manualiocheck)
%function S=OdStack(structin,method)
%
%Automated STXM raw data to optical density (OD) conversion.
%Particle regions are detected by a constant threshold condition or by using Otsu's thresholding algorithm.
%Particle-free image zones are used to calculate the normalization intensity Izero
%by averaging the intensity data over the particle-free zones.
%R.C. Moffet, T.R. Henn February 2009
%
%Manual Io selection option added by Matthew Fraund (2016)
%
%Inputs
%------
%structin       aligned raw data stack structure array (typically the output of the
%               FTAlignStack.m script).
%method         string defining the thresholding method:    -'O' thresholding using Otsu's method
%                                                           -'C' thresholding using a constant value
%                                                                (98.5% of maximal intensity)
%															
%plotflag		1 or 0 which controls the plotting of figures
%manualiocheck	string deciding whether to manually select Io region
%
%Outputs
%-------
%S              structure array containing the OD converted STXM data

%Making third input (plotflag) optional, defaulting to no figure
if nargin < 3
	plotflag = 0;
	manualiocheck = 'no';
elseif nargin < 4
	manualiocheck = 'no';
end


% create temporary variables
stack=structin.spectr;
eVlength=length(structin.eVenergy);

% construct the output struct
S=structin;
clear S.spectr;
S.spectr=zeros(size(structin.spectr,1),size(structin.spectr,2),size(structin.spectr,3));

% calculate imagesc limits
xAxislabel=[0,S.Xvalue];
yAxislabel=[0,S.Yvalue];

% particle masking & thresholding with constant threshold condition
if strcmp(method,'C')==1
    
    imagebuffer=mean(stack,3);
    imagebuffer=medfilt2(imagebuffer); % median filtering of the stack mean
    GrayImage=mat2gray(imagebuffer); % Turn into a greyscale with vals [0 1]
    Mask=zeros(size(GrayImage));
    Mask(GrayImage>=0.85)=1; % Detect particle free image regions
    
    % particle masking & thresholding using Otsu's method
elseif strcmp(method,'O')==1
    
    imagebuffer=mean(stack,3);  %% Use average of all images in stack
    GrayImage=mat2gray(imagebuffer); %% Turn into a greyscale with vals [0 1]
    GrayImage=imadjust(GrayImage,[0 1],[0 1],15); %% increase contrast
    %     GrayImage=imadjust(GrayImage,[0 1],[0 1],2); %%Contrast of 2 for maps
    
    Thresh=graythresh(GrayImage); %% Otsu thresholding
    Mask=im2bw(GrayImage,Thresh); %% Give binary image
    
    % Thresholding for maps -- doesnt need as much contrast adjustment
elseif strcmp(method,'map')==1
    imagebuffer=mean(stack,3);  %% Use average of all images in stack
    GrayImage=mat2gray(imagebuffer); %% Turn into a greyscale with vals [0 1]
    %     GrayImage=imadjust(GrayImage,[0 1],[0 1],15); %% increase contrast
    GrayImage=imadjust(GrayImage,[0 1],[0 1],2); %%Contrast of 2 for maps
    
    Thresh=graythresh(GrayImage); %% Otsu thresholding
    Mask=im2bw(GrayImage,Thresh); %% Give binary image
    
    
    
    % Thresholding method not defined
else
    
    display('Error! No thresholding method defined! Input structure not converted!')
    return
    
end
 S.mask=Mask;
% Izero extraction

%%%%% this section replaced with an optional input
%If you get annoyed with this prompt, uncomment the manualiocheck line and comment
%the line with "inputdlg" line. 
% manualiocheck = 'no';
% manualiocheck = inputdlg('Do you want to manually define an Io region? (Cancel continues with automatic method)','Manual Io Check',1,{'yes'});

if strcmp(manualiocheck,'yes') == 1
    avgstackfig = figure('Name','Define an Io region','NumberTitle','off');
    imagesc(imagebuffer);
    colormap gray;
    
    manualrect = getrect(avgstackfig);
    close(avgstackfig);
    
    cliprowstart = round(manualrect(2));
    cliprowend = round(cliprowstart + manualrect(4));
    clipcolstart = round(manualrect(1));
    clipcolend = round(clipcolstart + manualrect(3));
    manualIomatrix = imagebuffer((cliprowstart):(cliprowend),(clipcolstart):(clipcolend));
    
    Izero = zeros(eVlength,2);
    stdIzero = zeros(eVlength,2);
    errIzero = zeros(eVlength,2);
    
    Izero(:,1)=S.eVenergy;
    stdIzero(:,1) = S.eVenergy;
    errIzero(:,1) = S.eVenergy;
    
    for cnt=1:eVlength
        
        buffer=stack((cliprowstart):(cliprowend),(clipcolstart):(clipcolend),cnt); %This selects the region selected previously for energy "cnt"
        Izero(cnt,2)=mean(mean(buffer));
        stdIzero(cnt,2) = std(std(buffer));
        numIzero = numel(buffer);
        errIzero(cnt,2) = 1.96 .* stdIzero(cnt,2) ./ sqrt(numIzero); %1.96 comes from t-distribution with an alpha level of 0.05
        clear buffer
        
    end

    
else
    
    Izero=zeros(eVlength,2);
    stdIzero = zeros(eVlength,2);
    errIzero = zeros(eVlength,2);
    
    Izero(:,1)=S.eVenergy;
    stdIzero(:,1) = S.eVenergy;
    errIzero(:,1) = S.eVenergy;
    
    % loop over energy range of stack, calculate average vor each energy -> return_spec
    for cnt=1:eVlength
        
        buffer=stack(:,:,cnt);
        Izero(cnt,2)=mean(buffer(Mask==1));
        stdIzero(cnt,2) = std(buffer(Mask==1));
        numIzero = sum(sum(Mask));
        errIzero(cnt,2) = 1.96 .* stdIzero(cnt,2) ./ sqrt(numIzero); %1.96 comes from t-distribution with an alpha level of 0.05
        clear buffer
        
    end
end

S.Izero=Izero;
S.stdIzero = stdIzero;
S.errIzero = errIzero;
% stack conversion Intensity --> Optical Density

S.errOD = zeros(size(S.spectr));

for k=1:eVlength
    
    S.spectr(:,:,k)= -log(stack(:,:,k)/Izero(k,2)); 
    S.errOD(:,:,k) = (errIzero(k,2)) .* sqrt((1./stack(:,:,k).^2) + (1./Izero(k,2).^2));
end

% S.position = structin.position;

%Plot results
if plotflag==1
    figure
    subplot(2,2,1)
    imagesc(xAxislabel,yAxislabel,imagebuffer)
    axis image
    colorbar
    title('Raw Intensity Stack Mean')
    colormap gray
    xlabel('X-Position (µm)','FontSize',11,'FontWeight','normal')
    ylabel('Y-Position (µm)','FontSize',11,'FontWeight','normal')
    
    subplot(2,2,2)
    imagesc(xAxislabel,yAxislabel,mean(S.spectr,3));
    axis image
    colorbar
    colormap gray
    title('Optical Density Stack Mean')
    xlabel('X-Position (µm)','FontSize',11,'FontWeight','normal')
    ylabel('Y-Position (µm)','FontSize',11,'FontWeight','normal')
    
    
    subplot(2,2,3)
    imagesc(xAxislabel,yAxislabel,Mask)
    colorbar
    axis image
    title('Izero Region Mask')
    xlabel('X-Position (µm)','FontSize',11,'FontWeight','normal')
    ylabel('Y-Position (µm)','FontSize',11,'FontWeight','normal')
    
    subplot(2,2,4)
    plot(Izero(:,1),Izero(:,2))
    title('Izero')
    xlabel('Photon energy (eV)','FontSize',11,'FontWeight','normal')
    ylabel('Raw Counts','FontSize',11,'FontWeight','normal')
    
    if length(S.eVenergy)>1
        xlim([min(S.eVenergy),max(S.eVenergy)])
        ylim([0.9*min(S.Izero(:,2)),(max(S.Izero(:,2))+0.1*min(Izero(:,2)))])
    end
end
return
