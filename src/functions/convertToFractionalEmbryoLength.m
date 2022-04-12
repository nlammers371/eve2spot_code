function [spot_struct] = convertToFractionalEmbryoLength(resultsFolder,spot_struct)
        
% NL: this is adapted from the Garcia lab pipeline
% EllipsePosDV = [];


% DVExperiment = strcmpi(liveExperiment.experimentAxis, 'DV');

load([resultsFolder,filesep,'APDetection.mat'])
% Ellipses = getEllipses(liveExperiment);
% correctDV = exist([resultsFolder, filesep,'DV',filesep,'DV_correction.mat'], 'file');
% if correctDV
%     load([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'DV_correction.mat'], 'DV_correction');
% end

%Angle between the x-axis and the AP-axis
if exist('coordPZoom', 'var')
    APAngle=atan2((coordPZoom(2)-coordAZoom(2)),...
        (coordPZoom(1)-coordAZoom(1)));
else error('coordPZoom not defined. Was AddParticlePosition.m run?'); end

APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
% DVLength = APLength/2;

%The information in Ellipses is
%(x, y, a, b, theta, maxcontourvalue, time, particle_id)
for i=1:length(spot_struct)    
    spot_struct(i).APPosNucleus = NaN(size(spot_struct(i).yPosNucleus));
    for j=1:length(spot_struct(i).xPosNucleus)
        
        %Angle between the x-axis and the nucleus using the A position as a
        %zero
        try
            nucleusAngles_AP=atan2((double(spot_struct(i).yPosNucleus(j))-coordAZoom(2)),...
            (double(spot_struct(i).xPosNucleus(j))-coordAZoom(1)));
        catch
            error('goddamit')
        end
            
        %Distance between the points and the A point
        Distances=sqrt((coordAZoom(2)-spot_struct(i).yPosNucleus(j)).^2 ...
            + (coordAZoom(1)-spot_struct(i).xPosNucleus(j)).^2);
        
        APPositions=Distances.*cos(nucleusAngles_AP-APAngle);
        spot_struct(i).APPosNucleus(j)=APPositions/APLength*100;
        
%         if DVExperiment && correctDV
%             DVPositions=Distances.*sin(nucleusAngles_AP-APAngle);
%             EllipsePosDV{i}(j)=abs(DVPositions-DV_correction)/DVLength;
%         else
%             DVPositions=Distances.*sin(nucleusAngles_AP-APAngle);
%             EllipsePosDV{i}(j)=DVPositions/DVLength;
%         end
        
    end
end



end