%% This code was edited/adapted from StrainOutput in OpenXY, from the brilliant team over at BYU. https://github.com/BYU-MicrostructureOfMaterials/OpenXY
% Credits include @ZachClayburn, @bjack205, @eretnek, @dtfullwood, Kathryn Small, and others.

function meanshiftEff=StrainOutput_1(Settings,Components,DoShowGB,smin,smax,MaxMisorientation,IQcutoff)
%STRAINOUTPUT
%StrainOutput(Settings,Components,DoShowGB,smin,smax,MaxMisorientation)
%Plot strain components
%code bits for this function taken from Step2_DisloDens_Lgrid_useF_2.m
%authors include: Collin Landon, Josh Kacher, Sadegh Ahmadi, and Travis Rampton
%modified for use with HROIM GUI code, Jay Basinger 4/20/2011

data = Settings.data;
r = data.rows;
c = data.cols;

%%
if iscell(data.g)
    Angles=cell2mat(data.g);
    n=length(Angles)/3;
    angles=zeros(n,3);
    n=[1:n];
    for i=1:3
        angles(:,i)=Angles((n-1)*3+i)';
    end
else
    angles = Settings.data.g;
end

% Strain plots ************************

FArray = data.F;

%Ignore obviously bad points

%assignin('base','BadIndex',BadIndex)
%cutoff = (mean(SSE(SSE ~= inf)) + 15*std(SSE(SSE ~= inf)));
%BadPoints = SSE > cutoff;
%assignin('base','BadPoints',BadPoints)
%BadPoints = squeeze(BadPoints(1,:,:))';
%BadIndex = BadIndex(BadPoints);
%for j = 1:length(BadIndex)
%    data.F(:,:,BadIndex(j)) = eye(3);
%end
%%Ignore obviously bad points
%SSE = [Settings.fitMetrics.SSE];
%SSE = cell2mat(struct2cell(SSE));
%BadIndex = 1:length(SSE);
%assignin('base','BadIndex',BadIndex)
%cutoff = (mean(SSE(SSE ~= inf)) + 15*std(SSE(SSE ~= inf)));
%BadPoints = SSE > cutoff;
%assignin('base','BadPoints',BadPoints)
%BadPoints = squeeze(BadPoints(1,:,:))';
%BadIndex = BadIndex(BadPoints);
%for j = 1:length(BadIndex)
%    data.F(:,:,BadIndex(j)) = eye(3);
%end
%Ignore obviously bad points
SSE = Settings.SSE;
for i=1:length(SSE); Fn(i)=norm(squeeze(FArray(:,:,i)));end
BadIndex = 1:length(SSE);
cutoff = (mean(SSE(SSE ~= inf)) + 2*std(SSE(SSE ~= inf)));
BadPoints = SSE > cutoff;
BadIndex = BadIndex(BadPoints);
for j = 1:length(BadIndex)
    data.F(:,:,BadIndex(j)) = eye(3);
end
%FList = [data.F{:}];
%FArray=reshape(FList,[3,3,length(FList(1,:))/3]);


thisF=zeros(3);
FSample=FArray;
U = zeros(size(FArray));
for i=1:length(FArray(1,1,:))
    g=euler2gmat(angles(i,1),angles(i,2),angles(i,3));% this give g(sample to crystal)
    thisF(:,:)=FArray(:,:,i);
    %[R,U(1:3,1:3,i)]=poldec(thisF);
    FSample(:,:,i)=g'*thisF*g; %check this is crystal to sample
end
grainmap=Settings.GrainVals.grainID; % grain numbers for averaging F within each grain
ng=unique(grainmap);
zeroF=zeros(1,length(Fn));
zeroF(Fn==0)=1;
meanshiftF=FSample;
meanshiftE=zeros(3,3,length(Fn));
meanshiftEff=zeros(1,length(Fn));
Findex=[1:length(Fn)];
for kk=1:length(ng)
    onemap=ones(1,length(Fn));
    sumF=sum(FArray(:,:,grainmap==ng(kk)),3);
    onemap(grainmap~=ng(kk))=0;
    sumg=sum(grainmap==ng(kk))-sum(sum((onemap.*zeroF)));
    meanF=sumF/sumg;
    Thisindex=Findex(grainmap==ng(kk));
    for kkk=1:length(Thisindex)
        %meanshiftF(:,:,Thisindex(kkk))=inv(squeeze(meanF))*squeeze(FArray(:,:,Thisindex(kkk)));
        meanshiftF(:,:,Thisindex(kkk))=squeeze(meanF)/squeeze(FArray(:,:,Thisindex(kkk)));
        meanshiftE(:,:,Thisindex(kkk))=((squeeze(meanshiftF(:,:,Thisindex(kkk))))'*squeeze(meanshiftF(:,:,Thisindex(kkk)))-eye(3))/2;
        %meanshiftEff(Thisindex(kkk))=sqrt(2/9*((meanshiftE(1,1,Thisindex(kkk))-meanshiftE(2,2,Thisindex(kkk))).^2+(meanshiftE(2,2,Thisindex(kkk))-meanshiftE(3,3,Thisindex(kkk))).^2+(meanshiftE(1,1,Thisindex(kkk))-meanshiftE(3,3,Thisindex(kkk))).^2));
        exx=(2*meanshiftE(1,1,Thisindex(kkk))-meanshiftE(2,2,Thisindex(kkk))-meanshiftE(3,3,Thisindex(kkk)))/3;
        eyy=(-meanshiftE(1,1,Thisindex(kkk))+2*meanshiftE(2,2,Thisindex(kkk))-meanshiftE(3,3,Thisindex(kkk)))/3;
        ezz=(-meanshiftE(1,1,Thisindex(kkk))-meanshiftE(2,2,Thisindex(kkk))+2*meanshiftE(3,3,Thisindex(kkk)))/3;
        meanshiftEff(Thisindex(kkk))=2/3*sqrt(3*(exx^2+eyy^2+ezz^2)/2+3*(meanshiftE(1,2,Thisindex(kkk))^2+meanshiftE(1,3,Thisindex(kkk))^2+meanshiftE(2,3,Thisindex(kkk))^2));
    end
end
if strcmp(Settings.ScanType,'Hexagonal')
    NColsOdd = c;
    NColsEven = c-1;
    NRows = r;
    count = 1;
    
    NewstrainEff = zeros(NRows,NColsEven);
    if i==1 && j==1
        gmap = zeros(NRows,NColsEven);
    end
    for rr = 1:NRows
        if bitget(abs(rr),1)~=0 %odd
            
            for cc = 1:NColsOdd
                NewstrainEff(rr,cc) = meanshiftEff(count);
                if i==1 && j==1
                    gmap(rr,cc)=Settings.grainID(count);
                end
                count = count + 1;
            end
        else
            for cc = 1:NColsEven
                NewstrainEff(rr,cc) = meanshiftEff(count);
                if i==1 && j==1
                    gmap(rr,cc)=Settings.grainID(count);
                end
                count = count + 1;
            end
        end
    end
    meanshiftEff = NewstrainEff;
    
else
    meanshiftEff=reshape(meanshiftEff, [c r])';
    
end
cMap = jet(128);
cMap(1,:) = cMap(1,:)./3;
cMap(end,:) = cMap(end,:)./3;
% %     region=[113 15 120 150]           %if you want to isolate an area
figure;imagesc(meanshiftEff)
%     rectangle('position', region, 'edgecolor', 'r', 'linewidth', 2)
%     %draw a rectangle around specified region
AverageStrain=sum(meanshiftEff)/sum(meanshiftEff~=0);
title(['\epsilon_E_f_f, Average Strain: ' num2str(AverageStrain)],'fontsize',14)
%shading flat
%axis equal tight
% view(2)
%colorbar
%colormap(cMap)
%caxis([smin smax])
%    lines = [];
%    oldLines = [];
%     hold off
%     
%     condition = inpolygon(meanshiftEff, region);
%     effstrain_reg = meanshiftEff(condition);
%     figure; plot(effstrain_reg)
%     AverageStrain_reg=sum(effstrain_reg)/sum(effstrain_reg~=0);
%     title(['\epsilon_E_f_f, Average Strain: ' num2str(AverageStrain_reg)],'fontsize',14)
% %     plot(ebsd_reg, ebsd_reg('Nickel').orientations)
%     shading flat
%     axis equal tight
%     % view(2)
%     colorbar
%     colormap(cMap)
%     caxis([smin smax])
         lines = [];
         oldLines = [];
    
if DoShowGB && ~strcmp(Settings.ScanType,'Hexagonal')
    if isempty(lines)% Cobbled together way to speed things up
        lines = PlotGBs(Settings.grainID,[Settings.Nx Settings.Ny],Settings.ScanType);
        if isfield(Settings, 'grainsHaveBeenSplit') && ...
                Settings.grainsHaveBeenSplit
            redlines = PlotGBs(Settings.oldGrains.grainID, [Settings.Nx Settings.Ny], Settings.ScanType,gca,1,'k'); % was 1.5 and 'r'**
        end
    else
        hold on
        for ii = 1:size(lines,1)
            plot(lines{ii,1},lines{ii,2},'LineWidth',1,'Color','k')
        end%ii = 1:size(lines,1)
        if isfield(Settings, 'grainsHaveBeenSplit') && ...
                Settings.grainsHaveBeenSplit
            for ii = 1:size(redlines,1)
                plot(redlines{ii,1},redlines{ii,2},'LineWidth',1,'Color','k')% was 1.5 and 'r'**
            end%ii = 1:size(redlines,1)
        end
        hold off
    end
end 
