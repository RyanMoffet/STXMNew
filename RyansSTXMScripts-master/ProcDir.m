function ProcDir(RawDir,FinDir,plotflag,varargin)

%% This script processes stxm data found in "RawDir" and places it in P1Dir

datafolder=RawDir;
cd(datafolder) %% move to raw data folder
foldstruct=dir;

if ~isempty(varargin)
    badE=varargin{1};
else
    badE=[];
end

cnt=1;
numobj=length(dir);
if numobj==2
    beep
    disp('empty directory');
end
for i = 3:numobj %% loops through stack folders in raw data folder
    try cd(fullfile(datafolder,foldstruct(i).name)); %% move to data folder
            S=LoadStackRaw(pwd); %% load stack data
            S.particle=sprintf('%s',foldstruct(i).name); %% print particle name
            S=DeglitchStack(S,badE);
            %             filename=sprintf('%s%s%s',P1Dir,'\S',foldstruct(i).name); %% define directory to save file in
            %             cd(P1Dir)
            %             save(sprintf('%s%s','S',foldstruct(i).name)) %% save stack data in .mat file
%             figure,plot(S.eVenergy,squeeze(mean(mean(S.spectr))))
            S=AlignStack(S);
            if length(S.eVenergy)<5
                Snew=OdStack(S,'map',1);
            else
                Snew=OdStack(S,'O',1);
            end
            cd(FinDir)
            save(sprintf('%s%s','F',foldstruct(i).name))
            cd(FinDir);
    catch
        %         cd(fullfile(datafolder,foldstruct(i).name)); %% move back to raw data folder
        cd(datafolder); %% move back to raw data folder
        %         cd ..
        beep
        cnt=cnt+1;
    end
    
end
close all
