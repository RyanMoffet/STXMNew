function [Sout]=DirLabelOrgVolFrac(Snew)

%% Calculates organic volume fraction
%% Need to run CarbonMaps.m first

%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sout.VolFrac is a vector of individual particle volume fractions

%% Sout.ThickMap is a map of inorganic and inorganic thicknesses as well as a volume fraction map
%%                    ThickMap(:,:,1)=torg;
%%                    ThickMap(:,:,2)=tinorg;
%%                    ThickMap(:,:,3)=tsoot;
%%                    ThickMap(:,:,4)=volFraction;

%% To do list:
%% 1. Incorporate to DirLabelMapsStruct.m (which feeds the ParticleAnalysis2 app)

%% Define fundamental assumptions
orgComp='adipic';
Nc=6; % number of carbons 
inorgComp='NaCl';

inorgDens=2.16;
orgDens=1.36;
sootDens=1.8;

sootMW=12.01;
orgMW=146.1412;

%% Define varibles from CarbonMaps
Mask=Snew.binmap;
LabelMat=Snew.LabelMat;
%% Find energies
[~,preidx] = min(abs(Snew.eVenergy - 278)); %find(Snew.eVenergy<283);
ccidx=find(Snew.eVenergy>284 & Snew.eVenergy<286);
[~,postidx] = min(abs(Snew.eVenergy - 320));%find(Snew.eVenergy>310 & Snew.eVenergy<330); %length(Snew.eVenergy);
cidx=find(Snew.eVenergy>288 & Snew.eVenergy<289);
sp2idx=find(Snew.eVenergy>284.5 & Snew.eVenergy<285.6);
% if isempty(cidx)
%     cidx=postidx;
% end

postim=Snew.spectr(:,:,postidx);
preim=Snew.spectr(:,:,preidx);

%%
[uorgpre,uorgpost,uinorgpre,uinorgpost]=PreToPostRatioVolFrac(inorgComp,orgComp); %% get calculated cross sections
xin=uinorgpost./uinorgpre;
xorg=uorgpre./uorgpost; % added 1/9/17
% torg=(postim-xin.*preim)./(uorgpost.*orgDens+xin.*uorgpre.*orgDens); %
% old def before 1/9/2017 - note sign difference with the following:
torg=(postim-xin.*preim)./(uorgpost.*orgDens-xin.*uorgpre.*orgDens);
torg(torg<0)=0;
torg=torg./100; % convert to meters
% figure,imagesc(torg),colorbar;
% tinorg=(preim-uorgpre.*((postim-xin.*preim)./(uorgpost+xin.*uorgpre)))./...
%     (uinorgpre.*inorgDens);
tinorg=(preim-xorg.*postim)./(uinorgpre.*inorgDens-xorg.*uorgpost.*inorgDens);% updated due to error 1/9/17

tinorg(tinorg<0)=0;
tinorg=tinorg./100; % convert to meters
% figure,imagesc(tinorg),colorbar;
% figure,imagesc(volFraction),colorbar
MatSiz=size(LabelMat);
XSiz=Snew.Xvalue/MatSiz(1);
YSiz=Snew.Yvalue/MatSiz(2);
xdat=[0:XSiz:Snew.Xvalue];
ydat=[0:YSiz:Snew.Yvalue];
pixelsize=mean([XSiz,YSiz]);
%% Correct for thick inorganic regions assuming thickness is that of a cube
if sum(preim(preim>1.5))>0 % if the particle contains thick regions with OD>1.5
    for i=1:max(max(LabelMat)) % loop over each particle
        npix=length(preim(LabelMat==i & preim>1.5)); % count pixels in ith particle that have OD>1.5
        thickness=sqrt(npix).*pixelsize.*1e-6; % calculate inorganic the thickness based on that of a cube (works for big NaCl inclusions)
        tinorg(LabelMat==i & preim>1.5)=thickness; % for the ith particle, replace OD>1.5 thicknesses with geometric correction
    end
end
%% Correct for soot containing particles. 
tsoot=zeros(size(torg));
if ~isempty(sp2idx)
    tsoot=(Snew.sp2.*Snew.BinCompMap{3}.*torg.*orgDens.*Nc.*sootMW)./(sootDens.*orgMW);
    torg(Snew.BinCompMap{3}==1)=(1-Snew.sp2(Snew.BinCompMap{3}==1)).*torg(Snew.BinCompMap{3}==1);
end

%% calculate organic volume fraction
volFraction=torg./(torg+tinorg+tsoot);
volFraction(Mask==0)=0;

%% Integrate volume fractions for individual particles
VolFrac=zeros(max(max(LabelMat)),1);
for i=1:max(max(LabelMat)) % loop over particles
    sumOrgThick=nansum(torg(LabelMat==i));
    sumInorgThick=nansum(tinorg(LabelMat==i));
    sumSootThick=nansum(tsoot(LabelMat==i));
    VolFrac(i)=sumOrgThick./(sumOrgThick+sumInorgThick+sumSootThick);
end

%% Do figures
figure('Name',Snew.particle,'NumberTitle','off','Position',[1,1,715,869]);
subplot(2,2,1),imagesc(xdat,ydat,torg),colorbar, 
axis image, 
title('organic thickness (m)'), 
xlabel('X (\mum)');
ylabel('Y (\mum)'); 
subplot(2,2,2),imagesc(xdat,ydat,tinorg),colorbar,axis image,
title('inorganic thickness (m)'), 
xlabel('X (\mum)');
ylabel('Y (\mum)');
subplot(2,2,3),imagesc(xdat,ydat,tsoot),colorbar,
axis image,
title('soot thickness'),
xlabel('X (\mum)');
ylabel('Y (\mum)');
subplot(2,2,4),imagesc(xdat,ydat,volFraction),colorbar,
axis image,
title('organic volume fraction'),
xlabel('X (\mum)');
ylabel('Y (\mum)');

figure, hist(VolFrac), 
title('Organic Volume Fraction'),
xlabel('volume fraction'),
ylabel('#');
export_fig([Snew.particle,'OVFa'],'-png');

%% prepare outputs
Snew.VolFrac=VolFrac;
ThickMap(:,:,1)=torg;
ThickMap(:,:,2)=tinorg;
ThickMap(:,:,3)=tsoot;
ThickMap(:,:,4)=volFraction;
Snew.ThickMap=ThickMap;
Snew.VolFrac=VolFrac;
Sout=Snew;