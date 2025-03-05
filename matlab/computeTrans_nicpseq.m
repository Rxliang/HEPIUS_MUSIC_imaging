function Trans = computeTrans_nicpseq(varargin)
% computeTrans - Assign array specifications for known transducer types
% Arguments:
%    One argument:
%        1st arg - Trans
%        if Trans is a structure, return the full Trans structure. If ID string, return transducer name only.
%    Two arguments:
%        1st arg should be transducer name as string.
%        2nd arg should be parameter desired (eg 'maxHighVoltage').
%        Returns value of parameter.
%
% Use built-in defaults, but allow prior specification of some attributes, such as Trans.frequency.
%
% The Vantage version is modified to shift all transducer frequencies to the nearest value supported with 4X sampling using the 
% 250 MHz Vantage system clock. Note that for SW Versions 3.0+, the demodulation frequency can be specified independently
% of the transducer's nominal center frequency. Nevertheless, because the default is to set Receive.demodF = Trans.frequency,
% specifying a hardware supported frequency here is helpful. When several options are available near the nominal 
% center frequency, Trans.frequency is chosen to provide the highest ADC sampling rate. 
%
% The elevation aperture and elevation focal depth are not currently used in Verasonics software, 
% and when known are provided here only for convenience. A value of zero indicates that the value is unknown.
% Note that these 2 quanities will not be converted to wavelengths, regardless of Trans.units.
% 
% Trans.id is a double that holds the logical 32-bit value "0x00FFIIDD", where "FF" is the format 
% version and "IIDD" is the scanhead ID. The format versions known are:
%   0 - ATL/Philips transducers
%   1 - Verasonics vendor 1
%   2 - Verasonics vendor 2
% For example, "Trans.id = hex2dec('00000250');" specifies format version 0
% and the L7-4's ID.  This is specified in the code below as"Trans.id = hex2dec('0250');". 
%
% Trans.impedance is used for estimating transmit output current.
% The output current value is used in checking power dissipation limits for
% transmit devices within the system, but not for the transducer itself.
% It is the user's responsibility to keep transmit power levels below the
% safe operating limits of their transducer.
% Note that a non-realistic default impedance of 20 ohms has been set for MUX probes
%   (which should NEVER be used for push, and perhaps not for any extended transmit bursts), 
%   as a means of severely limiting transmit output levels when long duration profile 5 transmit bursts are used. 
%   In fact, the impedance of an HVMux switch is high enough that the switch will overheat during an
%   extended burst, and the artificially low impedance value is used by software to prevent most operation in profile 5.

% Copyright 2001-2016 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.


% Known transducers and their corresponding IDs and high voltage limits. *** Do not modify these values. ***
KnownTransducers = {'L7-4',      '000250',  96,   0;...
                    'L10-5',     '00074C',  70,   1;...
                    'L11-4v',    '01ABB4',  75,   0;...
                    'L11-5',     '000351',  96,   0;...
                    'L12-3v',    '01BBC3',  75,   1;...
                    'L12-5 38mm','000755',  75,   1;...
                    'L12-5 50mm','000B5B',  75,   1;...
                    'L22-8v',    '018A8A',  35,   1;...  
                    'L22-14v',   '02AB18',  30,   0;...
                    'L22-14vLF', '02AB17',  30,   0;...
                    'L38-22v',   '03BB30',  25,   1;...
                    'GEC1-5D',   '270225',  50,   0;...                    
                    'GEC1-5DPlus',   '270225',  50,   0; 
                    'GEC1-5DSingle',   '270225',  50,   0;                      
                    'GEC1-6D',   '27023E',  50,   0;...
                    'GE4CD',     '270200',  50,   0;...
                    'GEIC5-9D',  '270210',  50,   0;... 
                    'GEL3-12D',  '270243',  50,   1;...
                    'GEL3-12D_64',  '000000',  50,   1;...                    
                    'GEM5ScD',   '27022E',  50,   0;... 
                    'GEM5ScD_64','000000',  50,   0;...                                                             
                    'GEM5SD',    '270231',  50,   0;...                     
                    'GE6SD',     '270223',  50,   0;...                     
                    'GE6SD_64',  '000000',  50,   0;...                                         
                    'CL10-5',    '00034D',  96,   0;...
                    'CL15-7',    '00035C',  96,   0;...
                    'C4-2',      '0020D1',  96,   0;...
                    'C5-2',      '0020D9',  96,   0;...
                    'C5-2v',     '01AC52',  96,   0;...
                    'C7-4',      '00224E',  96,   0;... % no data
                    'C8-4V',     '00228C',  96,   0;... % no data
                    'C8-5',      '0022DE',  96,   0;... % no data
                    'C9-5ICT',   '00228B',  96,   0;...
                    'P3-2',      '004428',  96,   0;... % no data
                    'P4-1',      '00483E',  96,   0;...
                    'P4-2',      '004439',  96,   0;...
                    'P4-2v',     '01AA42',  96,   0;...
                    'P5-3',      '004529',  96,   0;... % no data
                    %'P6-3',      '004D3B',  96,   0;...
                    'P6-3_64',      '004D3B',  96,   0;...                    
                    'P7-4',      '00462A',  96,   0;... 
                    'H235',      '04235A',  96,   0;... % Sonic Concepts H-235
                    'Adapter Embedded S.T.E', '01FFFA',  96,   0;... % no data
                    '260 ZIF Backshell Kit',  '0BAD00',  96,   0;... % no data
                    '408 ZIF Backshell Kit',  '0BAD02',  96,   0;... % no data
                    '260 ZIF S.T.E. Fixture', '01FFFC',  96,   0;... % no data
                    '408 ZIF S.T.E. Fixture', '01FFFB',  96,   0; ...
                    'MUSIC_Biobox64',         '000000',  50,   0}; % no data
                

% The following "known" probes include ID info only (no specs): C7-4,  C8-4V,  C8-5,  P3-2,  P5-3,  P7-4  
% (Specifications for these probes are planned for a future release).
%                    'GE6SD',     '270223',  50,   0;...                     

switch nargin
    
    case 1
        Trans = varargin{1};
        if ~isstruct(Trans)  % if a structure is not provided as input, assume input is ID to translate into string.
            % input argument may be ID as a hex string, or the actual ID
            % numeric value
            if ischar(Trans)
                % convert hex string to number so it won't matter how many
                % leading zeros were provided
                probeID = hex2dec(Trans);
            else
                % not a string so presumably a numeric value
                probeID = Trans;
            end
            probeIDhex = num2str(probeID, '%06X');
            probenum = find(strcmpi(probeIDhex, KnownTransducers(:, 2)), 1);
            if isempty(probenum), Trans = 'Unknown';
            else Trans = KnownTransducers{probenum, 1}; % return the probe name (string value)
            end
            return
        end
        if ~isfield(Trans,'name'), error('computeTrans: Trans.name must be provided in input structure.'); end
        probenum = find(strcmpi(Trans.name, KnownTransducers(:, 1)), 1);
        if isempty(probenum), error('computeTrans: Trans.name not recognized as known transducer.'); end
        speedOfSound = 1.540;  % default speed of sound in mm/usec
        verbose = 2;
        if evalin('base','exist(''Resource'',''var'')&&isfield(Resource,''Parameters'')')
            if evalin('base','isfield(Resource.Parameters,''speedOfSound'')')
                speedOfSound = evalin('base','Resource.Parameters.speedOfSound')/1000; % speed of sound in mm/usec
            end
            if evalin('base','isfield(Resource.Parameters,''verbose'')')
                verbose = evalin('base','Resource.Parameters.verbose');
            end
        end

        % check for user-specified units, and print warning message if not found
        if ~isfield(Trans,'units') || isempty(Trans.units)
            fprintf(2, 'Warning: Trans.units not specified; selecting default units of mm.\n');
            fprintf(2, 'If script requires wavelength units, add an explicit definition of\n');
            fprintf(2, '"Trans.units = ''wavelengths'';" before calling computeTrans.\n');
            Trans.units = 'mm';
        end
        if ~strcmp(Trans.units, 'mm') && ~strcmp(Trans.units, 'wavelengths')
            error('computeTrans: Unrecognized value for Trans.units.  Must be ''mm'' or ''wavelengths''.');
        end

        % if Trans.frequency value has already been specified, we will use
        % it as is.  VSX and update() will confirm the value matches the
        % A/D sample rate constraints, and will exit with an error message
        % to the user if not.  Therefore we do not need to validate the
        % Trans.frequency value here (and could not, since we don't know
        % intended use of 4/3 sampling or interleave, etc.).
        if isfield(Trans,'frequency')
            if isempty(Trans.frequency)
                % if empty, remove it so cases below will assign default frequency
                Trans = rmfield(Trans, 'frequency');
            end
        end
        % also allow user-specified Bandwidth to override the default:
        if isfield(Trans,'Bandwidth')
            if isempty(Trans.Bandwidth)
                % if empty, remove it so cases below will assign default
                % Bandwidth
                Trans = rmfield(Trans, 'Bandwidth');
            end
        end

        Trans.lensCorrection = 0; % specify default value, in case it is not set for a particular transducer;
        Trans.id = hex2dec(KnownTransducers(probenum, 2));
        
        if strcmp(KnownTransducers{probenum, 1}, 'GEM5SD')
          KnownTransducers{probenum, 1} = 'GEM5ScD';
        end
        
        switch KnownTransducers{probenum, 1}
            case 'L7-4'
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                % Vantage:  5.208 is closest supported frequency to 5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 7]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector=1
                Trans.numelements = 128;
                Trans.elementWidth = .250; % width in mm
                Trans.spacingMm = .298;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 25; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0.887; % in mm units; was 3 wavelengths;
                Trans.impedance = [3 23.9-125i;  3.25 25.4-116i;  3.5 26-106i;  3.75 25.4-98.9i;  4 25.9-89.4i;  4.25 27-79.7i;...
                    4.5 32.8-72.6i;  4.75 39.2-66.2i;  5 46.1-69.6i;  5.25 46.5-72.4i;  5.5 41.9-71.6i;  5.75 43.2-69.8i;...
                    6 42.3-69.8i;  6.25 38.2-71i;  6.5 33.5-66.2i;  6.75 32-59.8i;  7 34.4-54.2i;  7.25 37.4-50.3i;...
                    7.5 42.3-48.2i;  7.75 47.8-47.9i;  8 53-51.3i];

            case 'L10-5'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .1729; % width in mm
                Trans.spacingMm = .1979;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 0; % nominal elevation focus depth from lens on face of transducer (unknown)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.183; % in mm units; was 6 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 10.5, ...
                                     'clock', 5, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,65), ...
                                     'VDASAperture', zeros(17,65,'uint8'));
                for i = 0:64, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                ApertureDataFilename = 'L10-5ApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:65
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L11-4v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % nominal frequency in MHz
                % Vantage:  6.25 is closest supported frequency to 6.4286 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.5, 10.5]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = 0.270; % width in mm
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 5; % active elevation aperture in mm (estimae)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.4785; % in mm units; was 5 wavelengths
                Trans.impedance = [ 3 17.4-123i;  3.25 18.6-110i;  3.5 20.5-100i;  3.75 20.4-91i;  4 23.3-83.3i;  4.25 21.1-78.3i;...
                    4.5 21.2-69.8i;  4.75 20.8-65.3i;  5 19.6-57.8i;  5.25 21.5-52.2i;  5.5 20.3-47.6i;  5.75 20.4-40.8i;...
                    6 21.3-37.2i;  6.25 20.4-31.5i;  6.5 22.7-25.8i;  6.75 23.1-22.8i;  7 23.5-17.3i;  7.25 26-13.8i;...
                    7.5 25.7-9.62i;  7.75 28.3-4.22i;  8 31.2-2.05i;  8.25 32.2+1.5i;  8.5 36.8+4.35i;  8.75 37.5+3.61i;...
                    9 38.7+7.27i;  9.25 39.7+7.01i;  9.5 38.5+11.5i;  9.75 42+15.4i;  10 42.9+18.2i;  10.25 47+24i;...
                    10.5 53.9+25.6i;  10.75 60.1+27.2i;  11 70.7+25.1i;  11.25 77.2+16.7i;  11.5 82.2+5.99i;  11.75 77-9.58i;...
                    12 65.4-15.8i;  12.25 54.3-17.7i;  12.5 44-15.6i;  12.75 35.2-10.6i;  13 29.1-3.17i;  13.25 26.2+4.47i;...
                    13.5 25.3+11i;  13.75 25.2+16.5i;  14 25.3+20.8i;  14.25 25+24.8i;  14.5 24.6+29.1i;  14.75 24.6+33.6i;...
                    15 25+38.1i;  15.25 25.6+42.4i;  15.5 26.5+46.5i;  15.75 27.5+50.7i;  16 28.9+54.9i];   

            case 'L11-5'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = 0.235; % width in mm
                Trans.spacingMm = 0.260;   % Spacing between elements in mm.
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.impedance = 50; % using default value 


            case 'L12-3v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 12]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .170; % width in mm
                Trans.spacingMm = .200;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.183; % in mm units; was 6 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 5.0, ...
                                     'clock', 5, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,65), ...
                                     'VDASAperture', zeros(33,65,'uint8'));
                for i = 0:64, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                ApertureDataFilename = 'L12-3v_ApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:65
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L12-5 38mm'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [5, 11]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 192;
                Trans.elementWidth = .1729; % width in mm
                Trans.spacingMm = .1979;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.365; % in mm units; was 12 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 100, ...
                                     'logicRail', 10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,65), ...
                                     'VDASAperture', zeros(17,65,'uint8'));
                for i = 0:64, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,64-i)]'; end
                ApertureDataFilename = 'L12-5_38mmApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:65
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L12-5 50mm'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 is closest supported frequency to 8.18 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [5, 11]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 256;
                Trans.elementWidth = .1703; % width in mm
                Trans.spacingMm = .1953;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 7.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 20; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.365; % in mm units; was 12 wavelengths;
                Trans.impedance = 51; % using default value for MUX probe
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,129), ...
                                     'VDASAperture', zeros(17,129,'uint8'));
                for i = 0:128, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                ApertureDataFilename = 'L12-5_50mmApertureData.txt';
                [fid, errmsg] = fopen(ApertureDataFilename);
                if ~isempty(errmsg), error ([errmsg ':  Filename = ' ApertureDataFilename ' is not in the path. See the Example Script directory for ' Trans.name]), end
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:129
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);

            case 'L22-8v'  % Kolo CMUT probe
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [8, 21.5]; end % approx. 90% relative bandwidth
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 256;
                Trans.elementWidth = 0.0703; % width in mm
                Trans.spacingMm = 0.108;   % Spacing between elements in mm.
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0; % in mm units; was 12 wavelengths;
                Trans.impedance = 5; % using default value for MUX probe
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 35; end
                Trans.HVMux = struct('highVoltageRails', 90, ...
                                     'logicRail', 6.5, ... %10.5, ...
                                     'clock', 8, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,129), ...
                                     'VDASAperture', zeros(17,129,'uint8'));
                for i = 0:128, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                fid = fopen('L22-8vApertureData.txt');
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:129
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);
                                        
            case 'L22-14v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate), 
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                    Trans.elevationApertureMm = 1.5; % active elevation aperture in mm
                    Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [10.00, 11.16-54.28i; 11.00, 12.02-44.73i; 12.00, 14.38-39.52i; 13.00, 14.19-33.50i;...
                    14.00, 14.43-27.21i; 15.00, 16.01-21.87i; 16.00, 16.82-17.73i; 17.00, 17.81-13.12i;...
                    18.00, 18.65-8.77i; 19.00, 20.90-4.20i; 20.00, 23.60-1.48i; 21.00, 25.54+0.61i; 22.00, 27.02+1.89i;...
                    23.00, 26.97+2.95i; 24.00, 25.99+5.39i; 25.00, 25.24+9.19i];
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit
                                        
            case 'L22-14vLF'
                if ~isfield(Trans,'frequency'), Trans.frequency = 15.625; end % nominal frequency in MHz
                % Note: use 15.625 MHz for 4X sampling (62.5 MHz sample rate), 
                % or 18.75 MHz for 4/3 sampling (25.0 MHz sample rate, 50 MHz A/D rate)
                % Manufacturer specified center frequency is 18.0 MHz +/- 10%
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [14, 22]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                    Trans.elevationApertureMm = 3; % Long Focus active elevation aperture in mm
                    Trans.elevationFocusMm = 20; % Long Focus nominal elevation focus depth from lens on face of transducer
                Trans.elementWidth = 0.08; % element width in mm; assumes 20 micron kerf
                Trans.spacingMm = 0.100;   % Spacing between elements in mm.
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                % Lens Correction: From mfr data sheet, matching layers are 0.05 mm thick at 2145 m/sec
                % average velocity, and lens is 0.48 mm thick at 1147 m/sec.
                % Thus the net effective lens thickness in mm is given by the following
                % expression, which evaluates to 0.6804 mm for 1540 m/sec velocity
                Trans.lensCorrection = 1000 * speedOfSound * (0.05/2145 + 0.48/1147); % velocities in m/sec; result in mm
                Trans.impedance = [10.00, 11.16-54.28i; 11.00, 12.02-44.73i; 12.00, 14.38-39.52i; 13.00, 14.19-33.50i;...
                    14.00, 14.43-27.21i; 15.00, 16.01-21.87i; 16.00, 16.82-17.73i; 17.00, 17.81-13.12i;...
                    18.00, 18.65-8.77i; 19.00, 20.90-4.20i; 20.00, 23.60-1.48i; 21.00, 25.54+0.61i; 22.00, 27.02+1.89i;...
                    23.00, 26.97+2.95i; 24.00, 25.99+5.39i; 25.00, 25.24+9.19i];
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % data sheet lists 30 Volt limit
                
            case 'L38-22v' % Kolo CMUT probe 
                if ~isfield(Trans,'frequency'), Trans.frequency = 30; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [22 38]; end 
                Trans.type = 0;     % linear=0   Array geometry is linear (x values only).
                Trans.connType = 1; % HDI=1 
                Trans.numelements = 256;
                Trans.elementWidth = 0.065; % width in mm
                Trans.spacingMm = 0.069;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 8; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:256,1) = Trans.spacingMm*(-((256-1)/2):((256-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0.15; % in mm units
                Trans.impedance = 5;  % artificially low to prevent overdriving the array
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 25; end % do not exceed 35 V!
                Trans.HVMux = struct('highVoltageRails', 90, ...
                    'logicRail', 6.5, ... %10.5, ...
                    'clock', 8, ...
                    'clockInvert', 0, ...
                    'polarity', 0, ...
                    'latchInvert', 0, ...
                    'Aperture', zeros(Trans.numelements,129), ...
                    'VDASAperture', zeros(17,129,'uint8'));
                for i = 0:128, Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),mod((i:i+127),128)+1,zeros(1,128-i)]'; end
                fid = fopen('L38-22vApertureData.txt');
                C = textscan(fid, '%*s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'CollectOutput',1);
                for i = 1:129
                    Trans.HVMux.VDASAperture(:,i) = uint8(hex2dec(C{1}(i,:)));
                end
                fclose(fid);
                
            case 'GEL3-12D'
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % 6.5 nominal frequency in MHz from GE data sheet (results in 50 MHz ADCrate)
                if ~isfield(Trans,'Bandwidth'), BW = 0.85; Trans.Bandwidth = 6.5*[1-BW/2, 1+BW/2]; end;  % 85% relative bandwidth from GE data sheet
                Trans.type = 0;     % =0 Array geometry is linear (x values only).
                Trans.connType = 7; % GE 408 pin connector on UTA (GE "D" family of probes)
                Trans.numelements = 256;
                Trans.spacingMm = 0.200;   % element pitch in mm from GE data sheet
                Trans.elementWidth = 0.9 * Trans.spacingMm; % width in mm (wild guess of 10% kerf; not from GE data sheet)
                Trans.elevationApertureMm = 5; % active elevation aperture in mm (spec)
                Trans.elevationFocusMm = 22; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1; % in mm units  (wild guess not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Add HVMux control structure; Voltage rail settings are based on notes
                % from Marc for the GE L3-12-D.  The clockInvert, polarity, and
                % latchInvert signals are not used by the GE L3-12-D so the values
                % assigned here are meaningless placeholders that will be ignored by
                % the UTA baseboard
                Trans.HVMux = struct('highVoltageRails', 100, ...
                                     'logicRail', 3.8, ...
                                     'clock', 10, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements,129), ...
                                     'VDASAperture', zeros(3,129,'uint8'));

                % Define Trans.HVMux.Aperture for each of the 129 HVMux apertures,
                % to select the desired set of elements assuming a 1:1 mapping to
                % channels for the first 128 elements.  Then redefine the
                % HVMux.Aperture array to include the mapping from elements at the
                % connector to HW system channels done by the GE UTA adapter
                % module, to map element positions at the connector to the
                % connector channel of the HW system that is wired to that
                % connector pin.  This mapping is equivalent to the Trans.Connector
                % element to channel mapping used for non-HVMux probes, as given
                % here, for the GE UTA module while using 128 active connector
                % channels
                EL2CH = [ 33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                         113   117   114   118   115   119   116   120   100   101    99   102    98   103    97   104, ...
                          49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                         108   112   107   111   106   110   105   109   124   128   123   127   122   126   121   125, ...
                           8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                          77    69    78    70    65    71    66    72    67    81    68    82    79    83    80    84, ...
                          17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                          88    96    87    95    86    94    85    93    76    89    75    90    74    91    73    92 ]';
                for i = 0:128
                    % first define the Aperture column as if there was 1:1 mapping, and
                    % then remap the Element entries using the EL2CH array
                    Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),EL2CH(mod((i:i+127),128)+1)',zeros(1,128-i)]';
                end

                % Defining the VDASAperture array:  The GE probe requires a 16 bit
                % control word to be sent to a CPLD in the probe, to tell it which
                % aperture to select through the HVMux chips.  Thus each
                % VDASAperture column will contain that 16 bit control word, with
                % the LS 8 bits in byte 3 and MS 8 bits in byte 2.  byte 1 is
                % set to two, the length of the table in bytes.  From the GE
                % documentation, the bit assignments are as follows in the control
                % word:
                % bits 15:12 Command, set to 0000 for legacy non-pipeline mode (9
                % usec load time)
                % bit ll Echo, set to 0 to disable echo back of the command received
                % bits 10:8 Row select (elevation aperture control) zero disables
                % all rows, 001 enables only the center row so we will assume that
                % is the value to use for L3-12 since it has only one elevation
                % row.
                % bits 7:0 Aperture select:  must be one less than the aperture
                % value we use in the script, since GE counts elements from zero
                % while we start at one.
                pldCmd = 0; % command value for bits 15:12
                pldEcho = 0; % value for bit 11
                pldRow = 1; % command valiue for bits 10:8
                % Set first byte for all apertures to two, the length of the table with two data bytes per aperture
                Trans.HVMux.VDASAperture(1, :) = uint8(2*ones(1, 129));
                % Set second byte to the 8 MSbits of the control word, which are
                % the same for all apertures since the only thing that changes is
                % the aperture select in bits 0:7
                MSbyte = 16*pldCmd + 8*pldEcho + pldRow;
                Trans.HVMux.VDASAperture(2, :) = uint8(MSbyte*ones(1, 129));
                % Set third byte to the aperture select value, one less than the
                % aperture index value we use which is in the range 1:129
                Trans.HVMux.VDASAperture(3, :) = uint8(0:128);
         case 'GEL3-12D_64'
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % 6.5 nominal frequency in MHz from GE data sheet (results in 50 MHz ADCrate)
                if ~isfield(Trans,'Bandwidth'), BW = 0.85; Trans.Bandwidth = 6.5*[1-BW/2, 1+BW/2]; end;  % 85% relative bandwidth from GE data sheet
                Trans.type = 0;     % =0 Array geometry is linear
                                    % (x values only).
                Trans.connType = 1; % UTA260
                Trans.numelements = 64;
                Trans.spacingMm = 0.200;   % element pitch in mm from GE data sheet
                Trans.elementWidth = 0.9 * Trans.spacingMm; % width in mm (wild guess of 10% kerf; not from GE data sheet)
                Trans.elevationApertureMm = 5; % active elevation aperture in mm (spec)
                Trans.elevationFocusMm = 22; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1; % in mm units  (wild guess not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Add HVMux control structure; Voltage rail settings are based on notes
                % from Marc for the GE L3-12-D.  The clockInvert, polarity, and
                % latchInvert signals are not used by the GE L3-12-D so the values
                % assigned here are meaningless placeholders that will be ignored by
                % the UTA baseboard
                Trans.HVMux = struct('highVoltageRails', 100, ...
                                     'logicRail', 3.8, ...
                                     'clock', 10, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(256,129), ...
                                     'VDASAperture', zeros(3,129,'uint8'));

                % Define Trans.HVMux.Aperture for each of the 129 HVMux apertures,
                % to select the desired set of elements assuming a 1:1 mapping to
                % channels for the first 128 elements.  Then redefine the
                % HVMux.Aperture array to include the mapping from elements at the
                % connector to HW system channels done by the GE UTA adapter
                % module, to map element positions at the connector to the
                % connector channel of the HW system that is wired to that
                % connector pin.  This mapping is equivalent to the Trans.Connector
                % element to channel mapping used for non-HVMux probes, as given
                % here, for the GE UTA module while using 128 active connector
                % channels
                EL2CH = [ 33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                         113   117   114   118   115   119   116   120   100   101    99   102    98   103    97   104, ...
                          49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                         108   112   107   111   106   110   105   109   124   128   123   127   122   126   121   125, ...
                           8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                          77    69    78    70    65    71    66    72    67    81    68    82    79    83    80    84, ...
                          17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                          88    96    87    95    86    94    85    93    76    89    75    90    74    91    73    92 ]';
                for i = 0:128
                    % first define the Aperture column as if there was 1:1 mapping, and
                    % then remap the Element entries using the EL2CH array
                    Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),EL2CH(mod((i:i+127),128)+1)',zeros(1,128-i)]';
                end

                % Defining the VDASAperture array:  The GE probe requires a 16 bit
                % control word to be sent to a CPLD in the probe, to tell it which
                % aperture to select through the HVMux chips.  Thus each
                % VDASAperture column will contain that 16 bit control word, with
                % the LS 8 bits in byte 3 and MS 8 bits in byte 2.  byte 1 is
                % set to two, the length of the table in bytes.  From the GE
                % documentation, the bit assignments are as follows in the control
                % word:
                % bits 15:12 Command, set to 0000 for legacy non-pipeline mode (9
                % usec load time)
                % bit ll Echo, set to 0 to disable echo back of the command received
                % bits 10:8 Row select (elevation aperture control) zero disables
                % all rows, 001 enables only the center row so we will assume that
                % is the value to use for L3-12 since it has only one elevation
                % row.
                % bits 7:0 Aperture select:  must be one less than the aperture
                % value we use in the script, since GE counts elements from zero
                % while we start at one.
                pldCmd = 0; % command value for bits 15:12
                pldEcho = 0; % value for bit 11
                pldRow = 1; % command valiue for bits 10:8
                % Set first byte for all apertures to two, the length of the table with two data bytes per aperture
                Trans.HVMux.VDASAperture(1, :) = uint8(2*ones(1, 129));
                % Set second byte to the 8 MSbits of the control word, which are
                % the same for all apertures since the only thing that changes is
                % the aperture select in bits 0:7
                MSbyte = 16*pldCmd + 8*pldEcho + pldRow;
                Trans.HVMux.VDASAperture(2, :) = uint8(MSbyte*ones(1, 129));
                % Set third byte to the aperture select value, one less than the
                % aperture index value we use which is in the range 1:129
                Trans.HVMux.VDASAperture(3, :) = uint8(0:128);
                
            case 'CL10-5'
                if ~isfield(Trans,'frequency'), Trans.frequency = 6.25; end % nominal frequency in MHz (this is the closest allowed frequency)
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4.5, 9]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.spacingMm = 0.200;   % Spacing between elements in mm.
                Trans.elementWidth = Trans.spacingMm - 0.05; % width in mm
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.1; % in mm units; was 9 wavelengths;
                Trans.impedance = 50; % using default value


            case 'CL15-7'
                if ~isfield(Trans,'frequency'), Trans.frequency = 8.929; end % nominal frequency in MHz
                % Vantage:  8.929 is closest supported frequency to 9.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 9*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = 0.16; % width in mm
                Trans.spacingMm = 0.178;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 8; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 15; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0.6899; % in mm units; was 4 wavelengths;
                Trans.impedance = 50; % using default value


            case 'C4-2'
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.976; end % nominal frequency in MHz
                % Vantage:  2.976 is closest supported frequency to 3.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                scanangle = 74.95 * (pi/180);    % degrees converted to radians
                radiusMm = 41.219;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .050;   % guess (in mm)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end

                Trans.impedance = 50; % using default value

            case 'C5-2'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.125; end % nominal frequency in MHz
                % Vantage:  3.125 is closest supported frequency to 3.2143 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4.5]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                scanangle = 74.95 * (pi/180);    % degrees converted to radians
                radiusMm = 41.219;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .050;   % guess (in mm)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (estimate)
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.impedance = [ 1 19.4-255i;  1.25 20.1-186i;  1.5 22.4-136i;  1.75 29-94.8i;  2 43.6-62i;  2.25 64.8-48.9i;...
                    2.5 70.2-54.9i;  2.75 54.3-49.1i;  3 45.4-30.9i;  3.25 42.9-13.3i;  3.5 42.5+1.9i;  3.75 42.1+16.4i;...
                    4 41.5+31.1i;  4.25 41.2+46.3i;  4.5 42.3+62.7i;  4.75 44.8+78.7i;  5 45.2+93i;  5.25 40.9+114i;...
                    5.5 42.6+144i;  5.75 51.2+177i;  6 66.5+214i];

            case 'C5-2v'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.7; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2.4, 5]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                radiusMm = 49.57;  % radius in mm.
                spacingMm = .508; % spacing in mm.
                kerf = .048;   % in mm.
                scanangle = 128*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (spec)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.035; % in mm units; was 3 wavelengths;
                Trans.impedance = [ 1 39.2-216i;  1.25 41-162i;  1.5 41.1-118i;  1.75 47.5-87i;  2 50.8-69.7i;  2.25 44.7-49.5i;...
                    2.5 43.3-27i;  2.75 48.6-8.53i;  3 50.4+4.4i;  3.25 50.7+20.7i;  3.5 56.5+38.2i;  3.75 69.1+51.8i;...
                    4 83.2+53.8i;  4.25 91.9+45.7i;  4.5 81.7+34.6i;  4.75 66+40.1i;  5 56.6+47i;  5.25 42.4+57.2i;...
                    5.5 30.7+75i;  5.75 25.3+95.5i;  6 24.3+115i;  6.25 24.8+132i;  6.5 26+148i;  6.75 27.8+164i;  7 30+179i];
         
         case 'GEC1-5DPlus'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.4; end % nominal frequency in MHz is 3.4
                if ~isfield(Trans,'Bandwidth'), BW = 0.85; 
                % WARNING: does not apply to single el txd
                Trans.Bandwidth = 3.4*[1-BW/2, 1+BW/2]; end 
                
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 7; % =7 GE connector
                Trans.numelements = 256;
                Trans.numelementsActual = 256; %192                
                radiusMm = 56.8;  % radius in mm.
                spacingMm = 5; % 0.350; % spacing in mm.
                kerf = .9*spacingMm;   % in mm. (guess)
                scanangle = Trans.numelementsActual*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                Trans.elevationApertureMm = 11.5; % active elevation aperture in mm (estimate)
                Trans.elevationFocusMm = 63; % nominal elevation focus depth from lens on face of transducer (estimate)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
%                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
 %               Trans.ElementPos(:,2) = 0;
  %              Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
   %             Trans.ElementPos(:,4) = Angle; % Orientation of
                                               % element with
                                               % respect to z axis.
%keyboard
                Trans.ElementPos(:,1) = zeros(size(Angle)).*(Trans.radiusMm*sin(Angle));
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = zeros(size(Angle)).*(Trans.radiusMm*cos(Angle)-Trans.radiusMm);
                Trans.ElementPos(:,4) = 0*Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                %Trans.lensCorrection = 1.; % in mm units; (guess)
                Trans.lensCorrection = 0; % in mm units; (guess)
                Trans.impedance = 50;    % Needs to be measured before using in profile 5
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 5; end
                
                % Now add the Trans.Connector map that is required, to reflect the mapping
                % from elements to system channels thrugh the UTA module and assuming
                % element numbering is as shown in the GE connector pinout (from Marc's GE
                % UTA module specification)
                Trans.Connector = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                
         case 'GEC1-5DSingle'
          if ~isfield(Trans,'frequency'), Trans.frequency = 1.9578430; end % nominal frequency
          %if ~isfield(Trans,'Bandwidth'), BW = 0.20; % from Cowe
                                                      % 2007 TCD paper 
            Trans.Bandwidth = [1.9186480 1.997038]; % measured -3dB from Evaldas %end 
            Trans.type = 1;     % Array geometry is curved linear (x and z values only).
            Trans.connType = 7; % =7 GE connector
            Trans.numelements = 2;
            radiusMm = 0;  % radius in mm.
            Trans.elementWidth = 6;  % width in mm (12 is total width)
            Trans.spacingMm = Trans.elementWidth; 
            Trans.radiusMm = 0;
            Trans.elevationApertureMm = 11.5; % active elevation aperture in mm (estimate)
            Trans.elevationFocusMm = 20;                 
            %   Set default element positions (units in mm).
            Trans.ElementPos = zeros(Trans.numelements,4);
            Trans.ElementSens = ones(1, 101);
            Trans.lensCorrection = 0; % in mm units; (guess)
            Trans.impedance = 20.7260-j/(2*pi*Trans.frequency*1e6*6.87128e-9); % measured
            if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 16.5; end
            
            % Now add the Trans.Connector map that is required, to reflect the mapping
            % from elements to system channels thrugh the UTA module and assuming
            % element numbering is as shown in the GE connector pinout (from Marc's GE
            % UTA module specification)
            Trans.Connector = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                                
            Trans.Connector = Trans.Connector([256 32]);
            % end GEC1-5DSingle
         case 'GEC1-5D'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.4; end % nominal frequency in MHz is 3.4
                if ~isfield(Trans,'Bandwidth'), BW = 0.85; Trans.Bandwidth = 3.4*[1-BW/2, 1+BW/2]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 7; % =7 GE connector
                Trans.numelements = 192;
                radiusMm = 56.8;  % radius in mm.
                spacingMm = .350; % spacing in mm.
                kerf = .9*spacingMm;   % in mm. (guess)
                scanangle = Trans.numelements*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                Trans.elevationApertureMm = 11.5; % active elevation aperture in mm (estimate)
                Trans.elevationFocusMm = 63; % nominal elevation focus depth from lens on face of transducer (estimate)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.; % in mm units; (guess)
                Trans.impedance = 50;    % Needs to be measured before using in profile 5
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Now add the Trans.Connector map that is required, to reflect the mapping
                % from elements to system channels thrugh the UTA module and assuming
                % element numbering is as shown in the GE connector pinout (from Marc's GE
                % UTA module specification)
                Trans.Connector = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                
                Trans.Connector = Trans.Connector(33:224);
                
           case 'GEC1-6D'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.9; end % nominal frequency in MHz is 3.4
                if ~isfield(Trans,'Bandwidth'), BW = 0.95; Trans.Bandwidth = 3.4*[1-BW/2, 1+BW/2]; end 
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 7; % =7 GE connector
                Trans.numelements = 192;
                radiusMm = 56.8;  % radius in mm.
                spacingMm = .350; % spacing in mm.
                kerf = .9*spacingMm;   % in mm. (guess)
                scanangle = Trans.numelements*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 11.5; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 66; % nominal elevation focus depth from lens on face of transducer (estimate)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.; % in mm units; (guess)
                Trans.impedance = 50;    % Needs to be measured before using in profile 5
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Now add the Trans.Connector map that is required, to reflect the mapping
                % from elements to system channels thrugh the UTA module and assuming
                % element numbering is as shown in the GE connector pinout (from Marc's GE
                % UTA module specification)
                Trans.Connector = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                
                Trans.Connector = Trans.Connector(33:224);
                
           case 'GE4CD'
                if ~isfield(Trans,'frequency'), Trans.frequency = 3.125; end % nominal frequency = 3.2 MHz
                if ~isfield(Trans,'Bandwidth'), BW = 0.7; Trans.Bandwidth = 3.2*[1-BW/2, 1+BW/2]; end 
                Trans.type = 1;     % 1= Array geometry is curved linear (x and z values only).
                Trans.connType = 7; % =7 GE connector
                Trans.numelements = 128;
                radiusMm = 60.;  % radius in mm.
                spacingMm = .478; % spacing in mm.
                kerf = 0.9*spacingMm;   % in mm. (guess)
                scanangle = 128*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                    Trans.elevationApertureMm = 13; % active elevation aperture in mm (spec)
                    Trans.elevationFocusMm = 67; % nominal elevation focus depth from lens on face of transducer (spec)
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.; % in mm units; (guess)
                Trans.impedance = 50;  % Needs to be measured before using in profile 5
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Now add the Trans.Connector map that is required, to reflect the mapping
                % from elements to system channels thrugh the UTA module and assuming
                % element numbering is as shown in the GE connector pinout (from Marc's GE
                % UTA module specification)
                Trans.Connector = [ 33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                         113   117   114   118   115   119   116   120   100   101    99   102    98   103    97   104, ...
                          49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                         108   112   107   111   106   110   105   109   124   128   123   127   122   126   121   125, ...
                           8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                          77    69    78    70    65    71    66    72    67    81    68    82    79    83    80    84, ...
                          17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                          88    96    87    95    86    94    85    93    76    89    75    90    74    91    73    92 ]';                

            case 'C9-5ICT'
                if ~isfield(Trans,'frequency'), Trans.frequency = 7.813; end % nominal frequency in MHz
                % Vantage:  7.813 and 6.944 are closest supported frequencies to 7.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = 7.5*[0.7, 1.3]; end % default assumed value of 60% of center frequency
                Trans.type = 1;     % Array geometry is curved linear (x and z values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                scanangle = 146.677 * (pi/180);    % degrees converted to radians
                radiusMm = 8.511;  % radius in mm.
                spacingMm = radiusMm * scanangle/(Trans.numelements-1); % spacing in mm.
                kerf = .025;   % guess (in mm)
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 0; % nominal elevation focus depth from lens on face of transducer (unknown)
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                firstangle = -(scanangle/2); %   first element angle = -0.65405 radians
                deltatheta = scanangle/(Trans.numelements-1);
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.impedance = 50; % using default value
            
            case 'GEIC5-9D'
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.682; end % nominal frequency = 5.8 MHz
                if ~isfield(Trans,'Bandwidth'), BW = 0.75; Trans.Bandwidth = 5.8*[1-BW/2, 1+BW/2]; end 
                Trans.type = 1;     % =1 Array geometry is curved (x,z values only).
                Trans.connType = 7; % =7 GE 408 pin connector on UTA (GE "D" family of probes)
                Trans.numelements = 192;
                radiusMm = 10.1;  % radius in mm.
                spacingMm = .138; % spacing in mm.
                kerf = spacingMm*0.2;   % in mm. (guess)
                    Trans.elevationApertureMm = 6; % active elevation aperture in mm (spec)
                    Trans.elevationFocusMm = 35; % nominal elevation focus depth from lens on face of transducer (spec)
                scanangle = Trans.numelements*spacingMm/radiusMm;    % arc length/radius
                Trans.radiusMm = radiusMm; % radius in mm.
                Trans.spacingMm = spacingMm;  % Spacing in mm.
                Trans.elementWidth = (spacingMm - kerf);  % width in mm
                deltatheta = spacingMm/radiusMm;
                firstangle = -(scanangle/2) + 0.5*deltatheta; % first element angle
                %   Set default element positions (units in mm).
                Trans.ElementPos = zeros(Trans.numelements,4);
                Angle = firstangle:deltatheta:-firstangle;
                Trans.ElementPos(:,1) = Trans.radiusMm*sin(Angle);
                Trans.ElementPos(:,2) = 0;
                Trans.ElementPos(:,3) = Trans.radiusMm*cos(Angle)-Trans.radiusMm;
                Trans.ElementPos(:,4) = Angle; % Orientation of element with respect to z axis.
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1; % in mm units  (wild guess not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Now add the Trans.Connector map that is required, to reflect the mapping
                % from elements to system channels thrugh the UTA module and assuming
                % element numbering is as shown in the GE connector pinout (from Marc's GE
                % UTA module specification)
                Trans.Connector = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                
                Trans.Connector = Trans.Connector(33:224);
            
 case 'GEM5ScD_64' 
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.841; end % nominal frequency is 2.8 MHz from GE data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.7 4.2]; end % from frequency response plot on GE data sheet
                Trans.type = 0;  % specify 1.5D array as linear array with single row of elements.
                % Note this 1.5 D array cannot be steered in the Y direction, since the outer rows
                % of elements are connected together, allowing only the elevation focus to be adjusted.  
                % For the Verasonics implementation, there is no dynamic focusing in the elevation
                % dimension - instead a dynamic aperture approach is used, with the near field using
                % acquisitions with the central row only, and the far field using all three rows.
                Trans.connType = 1; % HDI
                Trans.numelements = 64;
                Trans.ConnectorES = (1:Trans.numelements)'; % 1:1 mapping from transducer elements to element signals at connector
                % Note the actual transducer element array is 80 X 3 for a total of 240 active
                % elements, but the pairs of side elements are connected together within the probe
                % so only 160 element signals are present at the probe connector.
                Trans.spacingMm = 0.270;   % element pitch in mm; total active aperture is 21.6 mm (80 * 0.270)  
                Trans.elementWidth = 0.230; % width in mm from GE data sheet; kerf is 0.04 mm
                Trans.elevationApertureMm = 13; % active elevation aperture in mm from GE data sheet
                elevOffsetMm = 4.875; % distance in mm from center of the array to center of the side elements
                % total elevation aperture is 13 mm with a 6.5 mm center element and 3.25 mm side
                % element on each side, connected in parallel (and thus same total area; approx.
                % same impedance).  Center of a side element is therefore 4.875 mm from center of
                % the elevation aperture
                Trans.elevationFocusMm = 77; % nominal elevation focus depth from lens on face of 
                % transducer, from GE spec sheet. Note that by adjusting the delay and weighting of
                % the side element rows relative to the center row, the elevation focus can be
                % adjusted over some range.
                Trans.ElementPos = zeros(Trans.numelements,5);
                % If treated as a full 2D array, Trans.type would be set to 2 and the ElementPos array
                % would require five entries for each element (X, Y, Z position coordinates plus angular
                % orientation in two dimensions).  When treated as a type 0 linear, all elements are
                % at z = 0 and facing straight ahead, so both angular directions are zero. In addition,
                % elements 81:160, which represent the outer row pairs, are specified to be located
                % at y = 0 (similar to the inner row), since in the far field the 3 element rows are
                % treated as if they were a single element.
                %xpos =
                %Trans.spacingMm*(-((Trans.numelements-2)/4):((Trans.numelements-2)/4));
                xpos = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                %Trans.ElementPos(:,1) = [xpos, xpos]; % same xpos
                %for both rows 
                Trans.ElementPos(:,1) = [xpos]; 
                % If desiring to treat the GEMScD as a full 2D array (type 2), uncomment the lines below
                % to offset the outer row of elements in the y direction.  To provide zero delay focus
                % at elevationFocusMM,the outer rows can be offset in the z direction as well. This
                % allows computeTXDelays to place the elevation focus at the range focal point.
                % Trans.ElementPos(81:160,2) = elevOffsetMm; % add the y-offset for outer row
                % Trans.ElementPos(81:160,3) = Trans.elevationFocusMm - sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2);
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1; % in mm units  (wild guess; not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Notes on element numbering and element to channel mapping:  GE has defined element
                % signal numbers at the probe connector for 256 elements, with the convention that
                % probes using less than 256 elements will use the elements in the center of the
                % signal assignments for a 256 element array.  Thus for this 160 element probe the
                % GE element signal numbers at the connector will be 48 to 207; element signals 0:47
                % and 208:255 will not be used.
                
                % Note also that since Verasonics counts elements from 1 rather than zero, we will
                % be using elements signals 49 to 208 with the Verasonics numbering.
                
                % First create a connector array using the 1D array element
                % numbering GE uses for other probes, by selecting the
                % middle 160 elements as explained above.
                
               Connector1D = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                
                Connector1D = Connector1D(49:208); % 160 element probe uses center 160 channels, skipping 48 at each end
                
                % For the 1.5 D array of the M5ScD GE has mapped their
                % element numbers 1:16 to the first 16 center elements of
                % the 3-row array, and elements 17:32 to the first 16 outer
                % row element pairs.  This pattern repeats for each of the
                % five groups of 32 elements going across the array.  For
                % convenience in programming the system, we want to number
                % the center row as elements 1:80 and the outer rows as
                % elements 81:160.  Create Trans.Connector to reflect that
                % mapping:
                %Trans.Connector = zeros(160, 1); % create the 160 element column vector
                %Trans.Connector(1:80) = Connector1D([1:16, 33:48, 65:72, 81:88, 97:112, 129:144]); % center row of elements
                %Trans.Connector(81:160) = Connector1D([17:32, 49:64, 73:80, 89:96, 113:128, 145:160]); % outer rows of elements
                Trans.Connector = zeros(64, 1); % create the 64 element column vector
                tmp =  Connector1D([1:16, 33:48, 65:72, 81:88, 97:112, 129:144]); % center row of elements
                Trans.Connector(1:64) = tmp(9:72); % central 64
                                                   % elements of
                                                   % center row
                 Trans.Connector(1:64) = 1:64;
         case 'GEM5ScD' 
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.841; end % nominal frequency is 2.8 MHz from GE data sheet
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.7 4.2]; end % from frequency response plot on GE data sheet
                Trans.type = 0;  % specify 1.5D array as linear array with single row of elements.
                % Note this 1.5 D array cannot be steered in the Y direction, since the outer rows
                % of elements are connected together, allowing only the elevation focus to be adjusted.  
                % For the Verasonics implementation, there is no dynamic focusing in the elevation
                % dimension - instead a dynamic aperture approach is used, with the near field using
                % acquisitions with the central row only, and the far field using all three rows.
                Trans.connType = 7; % GE 408 pin connector on UTA (GE "D" family of probes)
                Trans.numelements = 160;
                % Note the actual transducer element array is 80 X 3 for a total of 240 active
                % elements, but the pairs of side elements are connected together within the probe
                % so only 160 element signals are present at the probe connector.
                Trans.spacingMm = 0.270;   % element pitch in mm; total active aperture is 21.6 mm (80 * 0.270)  
                Trans.elementWidth = 0.230; % width in mm from GE data sheet; kerf is 0.04 mm
                Trans.elevationApertureMm = 13; % active elevation aperture in mm from GE data sheet
                elevOffsetMm = 4.875; % distance in mm from center of the array to center of the side elements
                % total elevation aperture is 13 mm with a 6.5 mm center element and 3.25 mm side
                % element on each side, connected in parallel (and thus same total area; approx.
                % same impedance).  Center of a side element is therefore 4.875 mm from center of
                % the elevation aperture
                Trans.elevationFocusMm = 77; % nominal elevation focus depth from lens on face of 
                % transducer, from GE spec sheet. Note that by adjusting the delay and weighting of
                % the side element rows relative to the center row, the elevation focus can be
                % adjusted over some range.
                Trans.ElementPos = zeros(Trans.numelements,5);
                % If treated as a full 2D array, Trans.type would be set to 2 and the ElementPos array
                % would require five entries for each element (X, Y, Z position coordinates plus angular
                % orientation in two dimensions).  When treated as a type 0 linear, all elements are
                % at z = 0 and facing straight ahead, so both angular directions are zero. In addition,
                % elements 81:160, which represent the outer row pairs, are specified to be located
                % at y = 0 (similar to the inner row), since in the far field the 3 element rows are
                % treated as if they were a single element.
                xpos = Trans.spacingMm*(-((Trans.numelements-2)/4):((Trans.numelements-2)/4));
                Trans.ElementPos(:,1) = [xpos, xpos]; % same xpos for both rows 
                % If desiring to treat the GEMScD as a full 2D array (type 2), uncomment the lines below
                % to offset the outer row of elements in the y direction.  To provide zero delay focus
                % at elevationFocusMM,the outer rows can be offset in the z direction as well. This
                % allows computeTXDelays to place the elevation focus at the range focal point.
                % Trans.ElementPos(81:160,2) = elevOffsetMm; % add the y-offset for outer row
                % Trans.ElementPos(81:160,3) = Trans.elevationFocusMm - sqrt(elevOffsetMm^2 + Trans.elevationFocusMm^2);
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1; % in mm units  (wild guess; not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end
                
                % Notes on element numbering and element to channel mapping:  GE has defined element
                % signal numbers at the probe connector for 256 elements, with the convention that
                % probes using less than 256 elements will use the elements in the center of the
                % signal assignments for a 256 element array.  Thus for this 160 element probe the
                % GE element signal numbers at the connector will be 48 to 207; element signals 0:47
                % and 208:255 will not be used.
                
                % Note also that since Verasonics counts elements from 1 rather than zero, we will
                % be using elements signals 49 to 208 with the Verasonics numbering.
                
                % First create a connector array using the 1D array element
                % numbering GE uses for other probes, by selecting the
                % middle 160 elements as explained above.
                
               Connector1D = [ 97    98    99   100   101   102   103   104   105   106   111   108   109   110   107   112, ...
                    171   180   170   179   169   178   168   177   164   176   163   175   162   174   161   173, ...
                    120   121   116   122   115   123   114   124   113   125   119   126   118   127   117   128, ...
                    181   192   182   191   183   190   184   189   165   188   166   187   167   186   172   185, ...
                    33    44    34    43    35    42    36    41    37    48    38    47    39    46    40    45, ...
                    241   245   242   246   243   247   244   248   228   229   227   230   226   231   225   232, ...
                    49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64, ...
                    236   240   235   239   234   238   233   237   252   256   251   255   250   254   249   253, ...
                    8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13, ...
                    205   197   206   198   193   199   194   200   195   209   196   210   207   211   208   212, ...
                    17    28    18    27    19    26    20    25    21    32    22    31    23    30    24    29, ...
                    216   224   215   223   214   222   213   221   204   217   203   218   202   219   201   220, ...
                    83    84    82    85    81    86    80    87    68    69    67    70    66    71    65    72, ...
                    145   152   146   151   147   150   148   149   129   136   130   135   131   134   132   133, ...
                    92    96    91    95    90    94    89    93    88    79    75    78    74    77    73    76, ...
                    153   157   154   158   155   159   156   160   137   141   138   142   139   143   140   144 ]';
                
                Connector1D = Connector1D(49:208); % 160 element probe uses center 160 channels, skipping 48 at each end
                
                % For the 1.5 D array of the M5ScD GE has mapped their
                % element numbers 1:16 to the first 16 center elements of
                % the 3-row array, and elements 17:32 to the first 16 outer
                % row element pairs.  This pattern repeats for each of the
                % five groups of 32 elements going across the array.  For
                % convenience in programming the system, we want to number
                % the center row as elements 1:80 and the outer rows as
                % elements 81:160.  Create Trans.Connector to reflect that
                % mapping:
                Trans.Connector = zeros(160, 1); % create the 160 element column vector
                Trans.Connector(1:80) = Connector1D([1:16, 33:48, 65:72, 81:88, 97:112, 129:144]); % center row of elements
                Trans.Connector(81:160) = Connector1D([17:32, 49:64, 73:80, 89:96, 113:128, 145:160]); % outer rows of elements

         case 'GE6SD'
                if ~isfield(Trans,'frequency'), Trans.frequency = 4.6; end % nominal frequency in 4.6 MHz from GE data sheet
                if ~isfield(Trans,'Bandwidth'), BW = 0.85; Trans.Bandwidth = 4.6*[1-BW/2, 1+BW/2]; end   % 85 % relative bandwidth from GE data sheet
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 7; % GE 408 pin connector on UTA (GE "D" family of probes)
                Trans.numelements = 96; % was 256
                Trans.spacingMm = 0.16;   % element pitch in mm from GE data sheet
                Trans.elementWidth = 0.9 * Trans.spacingMm; % width in mm (wild guess of 10% kerf; not from GE data sheet)
                Trans.elevationApertureMm = 9; % active elevation aperture in mm (spec)
                Trans.elevationFocusMm = 50; % nominal elevation
                                             % focus depth from
                                             % lens on face of
                                             % transducer (spec),
                                             % was 22 for GEL3-12D

                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1; % in mm units  
                                          % (wild guess not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

                % JSM copied and modified from GE L3-12-D                               
                % Add HVMux control structure; Voltage rail settings are based on notes
                % from Marc for the GE L3-12-D.  The clockInvert, polarity, and
                % latchInvert signals are not used by the GE L3-12-D so the values
                % assigned here are meaningless placeholders that will be ignored by
                % the UTA baseboard
                numAperture = 1; % was 129
                Trans.HVMux = struct('highVoltageRails', 100, ...
                                     'logicRail', 3.8, ...
                                     'clock', 10, ...
                                     'clockInvert', 0, ...
                                     'polarity', 0, ...
                                     'latchInvert', 0, ...
                                     'Aperture', zeros(Trans.numelements, numAperture), ...
                                     'VDASAperture', zeros(3, numAperture,'uint8'));                                    

                % Define Trans.HVMux.Aperture for each of the 129 HVMux apertures,
                % to select the desired set of elements assuming a 1:1 mapping to
                % channels for the first 128 elements.  Then redefine the
                % HVMux.Aperture array to include the mapping from elements at the
                % connector to HW system channels done by the GE UTA adapter
                % module, to map element positions at the connector to the
                % connector channel of the HW system that is wired to that
                % connector pin.  This mapping is equivalent to the Trans.Connector
                % element to channel mapping used for non-HVMux probes, as given
                % here, for the GE UTA module while using 128 active connector
                % channels

                % connected RF channels are 1:32, 49:72, 77:84, 97:128
                
                EL2CH = [113   117   114   118   115   119   116   120   100   101    99   102    98   103    97   104 ...
                         49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64 ...
                         85:96 ...
                         124   128   123   127   122   126   121   125 ...
                         73:76 ...
                         108   112   107   111   106   110   105  109 ...
                           8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13 ...
                                                  33:48 ...
                          77    69    78    70    65    71    66    72    67    81    68    82    79    83    80    84 ...
                         17    28    18    27    19    26    20 ...
                        25   21    32    22    31    23    30    24 ...
                         29 ].';
                
                Trans.usedRFChannels = [1:32 49:72 77:84 97:128];
                i=0; 
                Trans.HVMux.Aperture(:,i+1) = EL2CH(Trans.usedRFChannels);
                %for i = 0:128
                    % first define the Aperture column as if there was 1:1 mapping, and
                    % then remap the Element entries using the EL2CH array
%                    Trans.HVMux.Aperture(:,i+1) = [zeros(1,i),EL2CH(mod((i:i+127),128)+1)',zeros(1,128-i)]';
 %               end

                % Defining the VDASAperture array:  The GE probe requires a 16 bit
                % control word to be sent to a CPLD in the probe, to tell it which
                % aperture to select through the HVMux chips.  Thus each
                % VDASAperture column will contain that 16 bit control word, with
                % the LS 8 bits in byte 3 and MS 8 bits in byte 2.  byte 1 is
                % set to two, the length of the table in bytes.  From the GE
                % documentation, the bit assignments are as follows in the control
                % word:
                % bits 15:12 Command, set to 0000 for legacy non-pipeline mode (9
                % usec load time)
                % bit ll Echo, set to 0 to disable echo back of the command received
                % bits 10:8 Row select (elevation aperture control) zero disables
                % all rows, 001 enables only the center row so we will assume that
                % is the value to use for L3-12 since it has only one elevation
                % row.
                % bits 7:0 Aperture select:  must be one less than the aperture
                % value we use in the script, since GE counts elements from zero
                % while we start at one.
                pldCmd = 0; % command value for bits 15:12
                pldEcho = 0; % value for bit 11
                pldRow = 1; % command valiue for bits 10:8
                % Set first byte for all apertures to two, the length of the table with two data bytes per aperture
                Trans.HVMux.VDASAperture(1, :) = uint8(2*ones(1, numAperture));
                % Set second byte to the 8 MSbits of the control word, which are
                % the same for all apertures since the only thing that changes is
                % the aperture select in bits 0:7
                MSbyte = 16*pldCmd + 8*pldEcho + pldRow;
                Trans.HVMux.VDASAperture(2, :) = uint8(MSbyte*ones(1, numAperture));
                % Set third byte to the aperture select value, one less than the
                % aperture index value we use which is in the range 1:129
                Trans.HVMux.VDASAperture(3, :) = uint8(0:numAperture-1);


         case 'GE6SD_64' % for 64 channel system
          if ~isfield(Trans,'frequency'), Trans.frequency = 4.6; end % nominal frequency in 4.6 MHz from GE data sheet
                if ~isfield(Trans,'Bandwidth'), BW = 0.85; Trans.Bandwidth = 4.6*[1-BW/2, 1+BW/2]; end   % 85 % relative bandwidth from GE data sheet
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI
                Trans.numelements = 64; % was 256
                Trans.spacingMm = 0.16;   % element pitch in mm from GE data sheet
                Trans.elementWidth = 0.9 * Trans.spacingMm; % width in mm (wild guess of 10% kerf; not from GE data sheet)
                Trans.elevationApertureMm = 9; % active elevation aperture in mm (spec)
                Trans.elevationFocusMm = 50; % nominal elevation
                                             % focus depth from
                                             % lens on face of
                                             % transducer (spec),
                                             % was 22 for GEL3-12D

                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1; % in mm units  
                                          % (wild guess not from GE data sheet)
                Trans.impedance = 50; % totally artificial made-up value not from GE data sheet
                if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

                EL2CH = [113   117   114   118   115   119   116   120   100   101    99   102    98   103    97   104 ...
                         49    57    51    58    53    59    55    61    50    63    52    60    54    62    56    64 ...
                         85:96 ...
                         124   128   123   127   122   126   121   125 ...
                         73:76 ...
                         108   112   107   111   106   110   105  109 ...
                           8    12     7    11     6    10     5     9     4    16     3    15     2    14     1    13 ...
                                                  33:48 ...
                          77    69    78    70    65    71    66    72    67    81    68    82    79    83    80    84 ...
                         17    28    18    27    19    26    20 ...
                        25   21    32    22    31    23    30    24 ...
                         29 ].';
                
                Trans.usedRFChannels = [1:32 49:72 77:84 97:128];
                
                %Trans.Connector = EL2CH(Trans.usedRFChannels(1:Trans.numelements));
                
            case 'P4-1'
                % The P4-1 is a 96 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.5, 3.5]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 96;
                Trans.Connector = [023 022 021 041 024 042 046 043 045 044 047 018 017 048 013 020 ...
                                   019 014 015 016 049 050 054 051 053 052 009 055 056 011 012 005 ...
                                   006 007 008 010 004 003 002 001 040 039 038 037 033 034 035 036 ...
                                   093 094 095 096 092 091 090 089 128 127 126 125 119 121 122 123 ...
                                   124 117 118 073 074 120 077 076 078 075 079 080 113 114 115 110 ...
                                   109 116 081 112 111 082 085 084 086 083 087 105 088 108 107 106]';
                kerf = .050;   % guess (in mm)
                Trans.spacingMm = 0.2950;   % Spacing between elements in mm.
                Trans.elementWidth = (Trans.spacingMm - kerf); % width in mm
                    Trans.elevationApertureMm = 16; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 80; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:96,1) = Trans.spacingMm*(-((96-1)/2):((96-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.464; % in mm units; was 4 wavelengths;
                Trans.impedance = [ 1 33.9-390i;  1.25 43.1-286i;  1.5 63.9-217i;  1.75 80.1-185i;  2 74.1-159i;...
                    2.25 73.9-126i;  2.5 84-99.2i;  2.75 98-93.3i;  3 89.5-99.4i;  3.25 63.6-89.5i;  3.5 44.2-63.6i;...
                    3.75 32.7-34.5i;  4 25.8-3.09i;  4.25 24.1+29.5i;  4.5 26.3+62.6i;  4.75 33.2+96.9i;  5 45.2+130i];   

        
            case 'P4-1'
                % The P4-1 is a 96 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [1.5, 3.5]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 96;
                Trans.Connector = [023 022 021 041 024 042 046 043 045 044 047 018 017 048 013 020 ...
                                   019 014 015 016 049 050 054 051 053 052 009 055 056 011 012 005 ...
                                   006 007 008 010 004 003 002 001 040 039 038 037 033 034 035 036 ...
                                   093 094 095 096 092 091 090 089 128 127 126 125 119 121 122 123 ...
                                   124 117 118 073 074 120 077 076 078 075 079 080 113 114 115 110 ...
                                   109 116 081 112 111 082 085 084 086 083 087 105 088 108 107 106]';
                kerf = .050;   % guess (in mm)
                Trans.spacingMm = 0.2950;   % Spacing between elements in mm.
                Trans.elementWidth = (Trans.spacingMm - kerf); % width in mm
                    Trans.elevationApertureMm = 16; % active elevation aperture in mm (estimate)
                    Trans.elevationFocusMm = 80; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:96,1) = Trans.spacingMm*(-((96-1)/2):((96-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.464; % in mm units; was 4 wavelengths;
                Trans.impedance = [ 1 33.9-390i;  1.25 43.1-286i;  1.5 63.9-217i;  1.75 80.1-185i;  2 74.1-159i;...
                    2.25 73.9-126i;  2.5 84-99.2i;  2.75 98-93.3i;  3 89.5-99.4i;  3.25 63.6-89.5i;  3.5 44.2-63.6i;...
                    3.75 32.7-34.5i;  4 25.8-3.09i;  4.25 24.1+29.5i;  4.5 26.3+62.6i;  4.75 33.2+96.9i;  5 45.2+130i];   

									    
            case 'P4-2'
                % The P4-2 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.  Elements 32-63 are
                %    wired to connector inputs 97-128.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.5; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.Connector = [1:32,97:128]';
                Trans.elementWidth = 0.2950; % width in mm
                Trans.spacingMm = 0.3200;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 3.08; % in mm units; was 5 wavelengths;
                Trans.impedance = [1.7, 98.0+91.0i; 2.0, 99.0; 2.5, 87.0+18.0i; 3.0, 77.0+16.0i; 3.5, 28.0+4.0i;...
                    4.0, 28.0+35.0i; 4.5, 42.0+83.0i; 5.0, 75.0+125.0i; 5.5, 128.0+159.0i; 6.0, 201.0+178.0i];

            case 'P4-2v'
                % The P4-2 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.  Elements 1-64 are
                %    wired to connector inputs 33-96.
                if ~isfield(Trans,'frequency'), Trans.frequency = 2.976; end % nominal frequency in MHz
                % Vantage:  2.976 is closest supported frequency to 3.0 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [2, 4]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.Connector = (33:96)';
                Trans.elementWidth = 0.250; % width in mm
                Trans.spacingMm = 0.300;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 0; % active elevation aperture in mm (unknown)
                    Trans.elevationFocusMm = 60; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 2.587; % in mm units; was 5 wavelengths;
                Trans.impedance = [ 1 33.8-366i;  1.25 49.4-252i;  1.5 64.1-196i;  1.75 75.5-140i;  2 99.5-118i;...
                    2.25 84.2-104i;  2.5 73.9-74.6i;  2.75 73.5-36i;  3 94.5-5.42i;  3.25 115-9.99i;  3.5 111-24.8i;...
                    3.75 74.9-14.8i;  4 66.9+17.5i;  4.25 70.2+22.4i;  4.5 37.7+37.8i;  4.75 26+75.3i;  5 24.8+107i;...
                    5.25 26.1+135i;  5.5 28.5+162i;  5.75 31.7+188i;  6 36.3+215i];   


            case 'P6-3'
                if ~isfield(Trans,'frequency'), Trans.frequency = 4.464; end % nominal frequency in MHz
                % Vantage:  4.464 is closest supported frequency to 4.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [3, 6]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = (.218-.025); % width in mm
                Trans.spacingMm = 0.218;   % Spacing between elements in mm
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:128,1) = Trans.spacingMm*(-((128-1)/2):((128-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.380; % in mm units; was 4 wavelengths;
                Trans.impedance = 48; % Z @ 3.4 Mhz

         case 'P6-3_64'
                if ~isfield(Trans,'frequency'), Trans.frequency = 4.464; end % nominal frequency in MHz
                % Vantage:  4.464 is closest supported frequency to 4.5 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [3, 6]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.elementWidth = (.218-.025); % width in mm
                Trans.spacingMm = 0.218;   % Spacing between elements in mm
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:128,1) = Trans.spacingMm*(- ...
                                                             ((128- ...
                                                               1)/2):((128-1)/2));
                keyboard
                Trans.ElementPos = Trans.ElementPos(1:64,:); 
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 1.380; % in mm units; was 4 wavelengths;
                Trans.impedance = 48; % Z @ 3.4 Mhz
                
            case 'P7-4'
                % The P7-4 is a 64 element array, so to use it with a 128 connector I/O system
                %    we have to define the connectivity in Trans.Connector.  Elements 32-63 are
                %    wired to connector inputs 97-128.
                if ~isfield(Trans,'frequency'), Trans.frequency = 5.208; end % nominal frequency in MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [4, 7]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 64;
                Trans.Connector = [1:32,97:128]';
                Trans.elementWidth = 0.1369; % width in mm
                Trans.spacingMm = 0.1711;   % Spacing between elements in mm.
                    % Trans.elevationApertureMm = 'unknown; % active elevation aperture in mm (estimate)
                    % Trans.elevationFocusMm = 'unknown; % nominal elevation focus depth from lens on face of transducer (estimate)
                Trans.ElementPos = zeros(Trans.numelements,4);
                %   Set default element x positions (units in mm).
                Trans.ElementPos(1:64,1) = Trans.spacingMm*(-((64-1)/2):((64-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 0.7; % in mm units; was 5 wavelengths;
                Trans.impedance = 50;%[1.7, 98.0+91.0i; 2.0, 99.0; 2.5, 87.0+18.0i; 3.0, 77.0+16.0i; 3.5, 28.0+4.0i;...
%                     4.0, 28.0+35.0i; 4.5, 42.0+83.0i; 5.0, 75.0+125.0i; 5.5, 128.0+159.0i; 6.0, 201.0+178.0i];

            case 'H235'
                if ~isfield(Trans,'frequency'), Trans.frequency = 0.6; end % nominal frequency in MHz
                % Vantage:  6.25 is closest supported frequency to 6.4286 MHz
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [.45 .75]; end 
                Trans.type = 0;     % Array geometry is linear (x values only).
                Trans.connType = 1; % HDI connector
                Trans.numelements = 128;
                Trans.elementWidth = 2.0; % width in mm
                Trans.spacingMm = 2.0;   % Spacing between elements in mm.
                    Trans.elevationApertureMm = 4; % active elevation aperture in mm (estimae)
                    Trans.elevationFocusMm = 35; % nominal elevation focus depth from lens on face of transducer (spec)
                Trans.ElementPos = zeros(Trans.numelements,4);
                Trans.ElementPos(:,1) = Trans.spacingMm*(-((Trans.numelements-1)/2):((Trans.numelements-1)/2));
                if ~isfield(Trans,'ElementSens')
                    % Set element sensitivity function (101 weighting values from -pi/2 to pi/2).
                    Theta = (-pi/2:pi/100:pi/2);
                    Theta(51) = 0.0000001; % set to almost zero to avoid divide by zero.
                    eleWidthWl = Trans.elementWidth * Trans.frequency/speedOfSound;
                    Trans.ElementSens = abs(cos(Theta).*(sin(eleWidthWl*pi*sin(Theta))./(eleWidthWl*pi*sin(Theta))));
                end
                Trans.lensCorrection = 10; % in mm units; was 5 wavelengths
                Trans.impedance =  60 ;   
%                 Trans.Connector = [ 21  22  23  24  25  26  27  28  29  30  31  32  33  34  50  51  ...
%                                     2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  ...
%                                     78  79  80  81  59  60  61  62  63  40  41  42  43  44  45  127 ...
%                                     109 91  92  73  74  75  55  56  57  58  37  38  39  19  20  1   ...
%                                     106 105 104 103 102 101 100 99  98  97  96  95  94  93  77  76  ...
%                                     125 124 123 122 121 120 119 118 117 116 115 114 113 112 111 110 ...
%                                     49  48  47  46  68  67  66  65  64  87  86  85  84  83  82  128 ...
%                                     18  36  35  54  53  52  72  71  70  69  90  89  88  108 107 126]'; el2ch
                 Trans.Connector = [64  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31 ...
                                    32 113  62  63   1   2   3   4   5   6   7   8   9  10  11  12 ...
                                    13  14 115 114  59  60  61  42  43  44  45  46  47 100  99  98 ...
                                    97  15  16 118 117 116  55  56  57  58  37  38  39  40  41 105 ...
                                   104 103 102 101 122 121 120 119  52  53  54  80  79  33  34  35 ...
                                    36 111 110 109 108 107 106 125 124 123  50  51  78  77  76  75 ...
                                    74  73  72  71  70  69  68  67  66  65 127 126  49  96  95  94 ...
                                    93  92  91  90  89  88  87  86  85  84  83  82  81 128  48 112]'; % ch2el 
          case 'MUSIC_Biobox64'
                %%
                if ~isfield(Trans,'frequency'), Trans.frequency = 13; end 
                if ~isfield(Trans,'Bandwidth'), Trans.Bandwidth = [5 18]; end
                Resource = evalin('base','Resource');
                Trans.type = 0; % linear array
                Trans.connType = 1;
                Trans.numelements = 64;
                Trans.spacingMm = 0.135; 
                Trans.ElementPos = zeros(Trans.numelements, 5);
                elemInArr = 64;
                Trans.ElementPos(1:elemInArr,1) = Trans.spacingMm*(-((elemInArr-1)/2):((elemInArr-1)/2));
                
                img1_32 = [36, 29, 35, 30, 34, 31, 33, 32, 40, 25, 39, 26, 38, 27, 37, 28, 44, 21, 43, 22, 42, 23, 41, 24, 48, 17, 47, 18, 46, 19, 45, 20];

                % IMG 33 - 64
                img33_64 = [52, 13, 51, 14, 50, 15, 49, 16, 56, 9, 55, 10, 54, 11, 53, 12, 60, 5, 59, 6, 58, 7, 57, 8, 64, 1, 63, 2, 62, 3, 61, 4];
                arr1 = [img1_32, img33_64];
                connectorArray = arr1;

                Trans.ConnectorES = (connectorArray)'; 
                Trans.elementWidth = 0.2; % width in mm 
                Trans.lensCorrection = 0; % in mm units;
                Trans.impedance = 65+49i; %Need to actually set this to match what's read from the Network Analyzer.
                Trans.elevationApertureMm = 2; % active elevation aperture in mm 
                Trans.elevationFocusMm = 7; % nominal elevation focus depth from lens on face of transducer
                Trans.maxHighVoltage = 96;
                Resource = evalin('base','Resource');
                Trans.elemInArr = 64;
                Trans.wvl = Resource.Parameters.speedOfSound/Trans.frequency/1000; % mm

            otherwise
                Trans.type = [];
                Trans.frequency = [];
                Trans.spacingMm = [];
                Trans.elementWidth = [];
                Trans.lensCorrection = [];
                Trans.ElementPos = [ 0 0 0 0 ];
                if verbose > 2
                    disp(' ');
                    disp(['computeTrans Status: Data not available for ', Trans.name]);
                    disp('Trans structure must be provided in user script.');
                    disp(' ');
                end
%                 fprintf(2, ['computeTrans: Data not available for ', Trans.name, ';\n']);
%                 fprintf(2, 'Trans structure must be provided in user script.\n');
%                 error(' ');

        end

        % Set a conservative value for maxHighVoltage, if not already defined
        if ~isfield(Trans,'maxHighVoltage'), Trans.maxHighVoltage = 50; end

        % Now convert all units as required, based on Trans.units
        scaleToWvl = Trans.frequency/speedOfSound; % conversion factor from mm to wavelengths

        % regardless of units, always provide spacing and radius in wavelengths
        Trans.spacing = Trans.spacingMm * scaleToWvl;   % Spacing between elements in wavelengths.
        if Trans.type == 1
            Trans.radius = Trans.radiusMm * scaleToWvl; % convert radiusMm to wavelengths
        end
        if strcmp(Trans.units, 'wavelengths')
            % convert all mm unit variables to wavelengths
            Trans.elementWidth = Trans.elementWidth * scaleToWvl;
            Trans.ElementPos(:,1) = Trans.ElementPos(:,1) * scaleToWvl;
            Trans.ElementPos(:,2) = Trans.ElementPos(:,2) * scaleToWvl;
            Trans.ElementPos(:,3) = Trans.ElementPos(:,3) * scaleToWvl;
            Trans.lensCorrection = Trans.lensCorrection * scaleToWvl;
        end
        
    case 2
        % Two inputs provided - 1st input is Trans.name, 2nd is parameter desired (currently only
        % 'maxHighVoltage' allowed).
        nameStr = varargin{1};
        if ischar(nameStr)
            probenum = find(strcmpi(nameStr, KnownTransducers(:, 1)), 1);
            % if 1st input is not a recognized transducer name, set probenum to
            % zero to flag as unrecognized transducer
            if isempty(probenum)
                probenum = 0; % special case of zero will trigger assignment of default value
            end
            Param = varargin{2};
            switch Param
                case 'maxHighVoltage'
                    if probenum == 0
                        % unrecognized transducer name, so assign default
                        % maxHighVoltage limit at hw max value
                        Trans = 96;
                    else
                        % return maxHighVoltage from KnownTransducers array
                        Trans = KnownTransducers{probenum, 3};
                    end
                case 'HVMux'
                    if probenum == 0
                        % unrecognized transducer name, so assign default
                        % of no HVMux
                        Trans = 0;
                    else
                        % return HVMux status from KnownTransducers array
                        Trans = KnownTransducers{probenum, 4};
                    end
                otherwise
                    error('computeTrans: Unrecognized parameter in 2nd argument.');
            end
        else
            error('computeTrans: When called with 2 inputs, 1st input must be transducer name.');
        end
        
    otherwise
        error('computeTrans: computeTrans accepts one or two input arguments.');

end
return
        

