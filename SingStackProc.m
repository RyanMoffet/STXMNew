function SingStackProc(path,plotflag,Name,varargin)

%% This script processes stxm data found in "path" and places it in
%% "path". In path, the stack data will be saved in the array Snew.
%% INPUT: path - string containing path to raw stxm data, processed datat will be saved here.
%%        Name - name of raw stxm data folder to be processed
%%        varargin{1} - vector of energies to be removed from stack
%% OCT 2009 RCM

datafolder=path;
cd(datafolder) %% move to raw data folder
foldstruct=dir;

if ~isempty(varargin)
    badE=varargin{1};
else
    badE=[];
end

ImpTest=0;
cnt=1;
numobj=length(dir);
for i = 3:numobj %% loops through stack folders in raw data folder
    if strcmp(foldstruct(i).name,Name)
        try cd(fullfile(datafolder,foldstruct(i).name)); %% move to data folder
            S=LoadStackRaw(pwd); %% load stack data
            S.particle=sprintf('%s',foldstruct(i).name); %% print particle name
            S=DeglitchStack(S,badE);
            %             filename=sprintf('%s%s%s',P1Dir,'\S',foldstruct(i).name); %% define directory to save file in
            %             cd(P1Dir)
            %             save(sprintf('%s%s','S',foldstruct(i).name)) %% save stack data in .mat file
            figure,plot(S.eVenergy,squeeze(mean(mean(S.spectr))))
            S=AlignStack(S);
            if length(S.eVenergy)==1
                Snew=OdStack(S,'O',plotflag);
            else
                Snew=OdStack(S,'C',plotflag);
            end
            cd(path)
            save(sprintf('%s%s','F',foldstruct(i).name))
            ImpTest=1;
        catch 
            
            cd(datafolder); %% move back to raw data folder
            %         cd ..
            disp('wrong path?')
            cnt=cnt+1;
        end
    else
        continue
    end
end


if ImpTest==0
    error('No import performed: Wrong filename or path?');
end
% cd(P1Dir)
% numobj=length(dir);
% for i = 3:numobj %% loops through stack matfiles
%     if strcmp(foldstruct(i).name,Name)
%         foldstruct=dir;
%         load(sprintf('%s',foldstruct(i).name));
%         Snew=AlignStack(S);
%         Snew=OdStack(Snew,'C');
%         S
%         cd(path)
%         save(sprintf('%s%s','F',foldstruct(i).name))
%         cd(P1Dir);
%     else
%         continue
%     end
% end
cd(path)