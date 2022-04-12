function [ ECOpar ] = ECOpuckCals( sn )

% function [ ECOpar ] = ECOpuckCals( sn )
%
% Log and retirieve ECOpuck calibr ations.  Listed are the scale factor and the dark
% count so chlorophyll concentration, particle concentration and CDOM
% concentration can be calculated.
%
% KIM 08.12

switch sn
    case 'FLBBCD2K-2710'
        
        %chlorophyll fluorescence to chlorophyll concentration in microg/l
        ECOpar.chl.dark = 50; % [counts]
        ECOpar.chl.sf = 0.0073; % microg/l/count
        
        % backscattering to particle concentration in (m^-1 sr^-1)
        ECOpar.backscat.dark = 51; % [counts]
        ECOpar.backscat.sf = 1.888e-6; % (m^-1 sr^-1)/count
        
        % change CDOM fluorescence to CDOM concentration ppb
        ECOpar.CDOM.dark = 49; % [counts]
        ECOpar.CDOM.sf = 0.0905; % ppb/count
        
    case 'FLBBCD2K-3383'
        
        %chlorophyll fluorescence to chlorophyll concentration in microg/l
        ECOpar.chl.dark = 46; % [counts]
        ECOpar.chl.sf = 0.0073; % microg/l/count
        
        % backscattering to particle concentration in (m^-1 sr^-1)
        ECOpar.backscat.dark = 47; % [counts]
        ECOpar.backscat.sf = 1.659e-6; % (m^-1 sr^-1)/count
        
        % change CDOM fluorescence to CDOM concentration ppb
        ECOpar.CDOM.dark = 41; % [counts]
        ECOpar.CDOM.sf = 0.0906; % ppb/count
end