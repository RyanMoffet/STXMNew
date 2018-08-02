function multpartimSootCarbox(InDir)


%%  [LabelCnt,OutSpec,size,label,CmpSiz,OutOCSpec,OutRad,RadStd,OutCmp]...
%%     =DirLabel(InDir,sample,saveflg,savepath,sp2)
%% INPUT: InDir - string containing directory where stacks are saved as
%%                .mat files
%%        sample - cell array of strings indicating sample names
%%
%% OUTPUT: LabelCnt - matrix containing counts for particle classes
%%                    columns: [OC,ECOC,INECOC,INOC] rows: [sample1,
%%                    sample2, sample3...]
%%         OutSpec - cell array containing grand average spectra for the
%%                   different samples. [energy,OD]
%%         size - cell array containing sizes for individual particles
%%                within the different samples
%%         label - cell array containing labels for individual particles
%%         CmpSiz - cell array of sample matrices containing component
%%                  areas (in um^2). Rows correspond to particles and columns are as
%%                  follows: ['OCArea','InArea','Sp2Area','KArea','CO3Area','Total
%%                  Area']
%%         OutOCSpec - Cell array of spectra for OC particles only
%%         OutRad - Cell array of average radial particle scans for each
%%                  class. Row is relative dist from center and column is
%%                  particle type.
%%         OutCmp - cell array of average spectra for component specified
%%                  in call to ComponentSpec.m
%% OCT 2009, RCM
cd(InDir)
startingdir=dir;
LabelCnt=zeros(4,length(InDir));
figsiz=1;
scrsz=get(0,'ScreenSize');
figure('units','normalized','position',[0,0,figsiz.*scrsz(4)/scrsz(3),figsiz]);

%% begin loop over sample directories
directory=dir;
subnum=ceil(sqrt(length(directory)-2));
nstack=length(directory)-2;
msize=ceil(sqrt(nstack));

x=zeros(msize,msize);y=x;
ctr=1;
for j=1:length(directory)
    ind=strfind(directory(j).name,'.mat');
    if ~isempty(ind)
        load(directory(j).name)
        x(ctr)=Snew.Xvalue;
        y(ctr)=Snew.Yvalue;
        ctr=ctr+1;
    end
end
maxdim=max([sum(x')+x(:,end)',sum(y)+y(end,:)]);
x=x./maxdim;y=y./maxdim;
xmax=max(x');
ymax=max(y);
rcnt=1;
ccnt=1;
ybot=1-ymax(1);
xpos=0;
xmove=0;
ctr=1;
for j=1:length(directory)  %%
    ind=strfind(directory(j).name,'.mat');
    if ~isempty(ind)
        load(directory(j).name);
        sprintf(directory(j).name)
        %% Run DiffMaps and label particles
        Sinp=SootCarboxMap(Snew,0);
        if ~isfield(Sinp,'SootCarbox')
            beep
            disp('Warning: SootCarbox doesnt exist, image scaling may be screwed');
            continue
        else
            clear varitest
            %set up axes
            MatSiz=size(Sinp.SootCarbox);
            XSiz=x(ctr)/MatSiz(1);
            YSiz=y(ctr)/MatSiz(2);
            xdat=[0:XSiz:x(ctr)];
            ydat=[0:YSiz:y(ctr)];
            if ccnt<=msize
                ccnt=ccnt+1;
                xpos=xpos+xmove;
                ypos=ybot+(ymax(rcnt)-y(ctr))./2;
            else
                ccnt=1;
                xpos=0;
                rcnt=rcnt+1;
                ybot=ybot-ymax(rcnt);
                ypos=ybot+(ymax(rcnt)-y(ctr))./2;
                %             xpos=xpos+x(j)./2;
                ccnt=ccnt+1;
            end
            axes('position',[xpos,ypos,x(ctr),y(ctr)])
            imagesc([0, Snew.Xvalue],[0,Snew.Yvalue],Sinp.SootCarbox),...
                caxis([0,5]),hold on,
            [B,L] = bwboundaries(Sinp.binmap);
            [ymaxT,xmaxT,emaxT]=size(S.spectr);
            for k = 1:length(B)
                boundary = B{k};
                plot(boundary(:,2).*(S.Xvalue./xmaxT), ...
                    boundary(:,1).*(S.Yvalue./ymaxT), 'w', 'LineWidth', 2)
            end
            
            set(gca, 'XTick', [], 'YTick', [])
            %         axis('image')
            xmove=x(ctr);
            ctr=ctr+1;
        end
    end
end

return


