%%Loading in datafiles and setting parameters
clear
clc
%OpenXY Strain Output .MAT%
OpenXYout = load('AnalysisParams_Ti64_HT650C_T=0_scan2b_p15step_1200x_20kV_SS5_WD18 cleaned-Dil_remapped.mat');
Settings = OpenXYout.Settings;
%EBSD ang/csv scan/grain file (copy and paste MTEX import wizard data here)
%Not needed, pulling grainID from OpenXYout.Settings.GrainVals
%%MTEX loading(copy paste from import_wizard('EBSD'))
% crystal symmetry, repleace and use correct values for your material
CS = {... 
  'notIndexed',...
  crystalSymmetry('432', [3.3 3.3 3.3], 'mineral', 'Ti (beta)', 'color', [0.53 0.81 0.98])};
% plotting convention
setMTEXpref('xAxisDirection','north');
setMTEXpref('zAxisDirection','outOfPlane');
pname = 'Z:\S11 Drop\Ti64 Remapping';
fname = [pname '\Corr_Ti64_HT650C_T=0_scan2b_p15step_1200x_20kV_SS5_WD18 cleaned-Dil_remapped.ang'];
ebsd = EBSD.load(fname,CS,'interface','ang',...
  'convertEuler2SpatialReferenceFrame');

[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',15*degree);

%% Calculate strain from OpenXY MAT, plot to XY grid, map unique grain ID from MTEX. 
%ONLY WORKS WITH SQUARE SCAN DATA RIGHT NOW BECAUSE HEX SUX
equivstrain_avg = StrainOutput_1(Settings,1,0,0.08,5,0)';
%defining grains, OpenXY style
%grainID = Settings.GrainVals.grainID;
%defining grains, MTEX style
grainID = ebsd.grainId;
grainID_xy_map = reshape(grainID,[Settings.Ny, Settings.Nx])';

%%Find XY coordinates of unique grain IDs, sum all strains at those
%%equivalent points, reassign an average value to those XY points.
unique_grains = unique(grainID_xy_map);
%grain_ID_point_store = zeros(max(unique_counts),max(unique_counts),size(unique_grains,1));
unique_counts = histc(grainID_xy_map(:),unique_grains);
grain_avg_strain = zeros(Settings.Ny, Settings.Nx)';
for i = 1:size(unique_grains,1)
    point_count = 1;
    overall_percent = 100*(i/size(unique_grains,1));
    temp_grain = zeros(unique_counts(i),2);
    %Scan thru grainID array, finding points that equal a given unique grainID
    for x = 1:1:Settings.Nx
        for y = 1:1:Settings.Ny
            
            if grainID_xy_map(x,y) == unique_grains(i)
                %grain_temp_store = grain_ID_point_store(:,:,i);
                temp_grain(point_count,:)=[x y];
                point_count = point_count + 1;
                grain_percent =  100*(point_count/unique_counts(i));
                percent_out = ['Overall Done: ', num2str(overall_percent), ', Grain Done: ', num2str(grain_percent)];
                disp(percent_out)
            end
        end
    end
    temp_grain_strain = zeros(size(temp_grain,1),1);
    for j = 1:1:size(temp_grain,1)      
        temp_grain_strain(j) = equivstrain_avg(temp_grain(j,1),temp_grain(j,2));  
    end
    temp_grain_mean_strain = mean(temp_grain_strain);
    for a = 1:1:size(temp_grain,1)      
        grain_avg_strain(temp_grain(a,1),temp_grain(a,2)) = temp_grain_mean_strain;  
    end
end
figure;
imagesc(grain_avg_strain);
lines = PlotGBs(grainID,[Settings.Ny Settings.Nx],Settings.ScanType);
AverageStrain=sum(grain_avg_strain)/sum(grain_avg_strain~=0);
title(['\epsilon_E_f_f, Average Strain: ' num2str(AverageStrain)],'fontsize',14)
colormap(jet);
colorbar
caxis([0 .02]);
