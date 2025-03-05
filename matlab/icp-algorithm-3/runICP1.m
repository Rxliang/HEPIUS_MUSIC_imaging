function varargout = runICP1(varargin)
% RUNICP1 MATLAB code for runICP1.fig
%      RUNICP1, by itself, creates a new RUNICP1 or raises the existing
%      singleton*.
%
%      H = RUNICP1 returns the handle to a new RUNICP1 or the handle to
%      the existing singleton*.
%
%      RUNICP1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUNICP1.M with the given input arguments.
%
%      RUNICP1('Property','Value',...) creates a new RUNICP1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before runICP1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to runICP1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help runICP1

% Last Modified by GUIDE v2.5 22-Dec-2016 16:56:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @runICP1_OpeningFcn, ...
                   'gui_OutputFcn',  @runICP1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before runICP1 is made visible.
function runICP1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to runICP1 (see VARARGIN)
global figState

% initialize
set(handles.radioFig1, 'Value', 1);
figState.fig1.segment = 1;
figState.fig1.spectrum = 2;
figState.fig1.ARorder = 15;
figState.fig1.process = 1;
figState.fig1.greyThresh = 0.0001;
figState.fig1.SNSIThresh = 15;
figState.fig1.TVRThresh = 2.5;
figState.fig1.greyThreshInitial = [];

figState.fig2.segment = 2;
figState.fig2.spectrum = 2;
figState.fig2.ARorder = 15;
figState.fig2.process = 1;
figState.fig2.greyThresh = 0.0001;
figState.fig2.SNSIThresh = 15;
figState.fig2.TVRThresh = 2.5;
figState.fig2.greyThreshInitial = [];

set(handles.sliderSNSIThreshold, 'Min', 5);
set(handles.sliderSNSIThreshold, 'Max' ,35);
set(handles.sliderSNSIThreshold, 'SliderStep', [1 5]/30);
set(handles.sliderSNSIThreshold, 'Value', 15);

set(handles.sliderGreyScaleThreshold, 'Min', -7);
set(handles.sliderGreyScaleThreshold, 'Max' ,-2);
set(handles.sliderGreyScaleThreshold, 'SliderStep', [.01 .1]);
set(handles.sliderGreyScaleThreshold, 'Value', -4);

set(handles.sliderTVRFilter, 'Min', .5);
set(handles.sliderTVRFilter, 'Max' ,4.5);
set(handles.sliderTVRFilter, 'SliderStep', [.1 .5]/4);
set(handles.sliderTVRFilter, 'Value', 2.5);

set(handles.radiobuttonNormalize, 'Value', 1);

% Choose default command line output for runICP1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes runICP1 wait for user response (see UIRESUME)
% uiwait(handles.figureGetEnvelope);


% --- Outputs from this function are returned to the command line.
function varargout = runICP1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuSegment.
function popupmenuSegment_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSegment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSegment
global figState

if get(handles.radioFig1, 'value') == 1
    figState.fig1.segment = get(hObject, 'Value');
else
    figState.fig2.segment = get(hObject, 'Value');
end


% --- Executes during object creation, after setting all properties.
function popupmenuSegment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenuSelectSpectrum.
function popupmenuSelectSpectrum_Callback(hObject, eventdata, handles)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuSelectSpectrum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuSelectSpectrum
global figState

if get(handles.radioFig1, 'value') == 1
    figState.fig1.spectrum = get(hObject, 'Value');
else
    figState.fig2.spectrum = get(hObject, 'Value');
end


% --- Executes during object creation, after setting all properties.
function popupmenuSelectSpectrum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuSelectSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in radioFig1.
function radioFig1_Callback(hObject, eventdata, handles)

global figState

set(handles.radioFig1, 'value', 1);
set(handles.radioFig2, 'value', 0);

set(handles.popupmenuSegment, 'Value', figState.fig1.segment);
set(handles.popupmenuSelectSpectrum, 'Value', figState.fig1.spectrum);
set(handles.editARorder, 'String', num2str(figState.fig1.ARorder));
set(handles.popupmenuProcessingFn, 'Value', figState.fig1.process);
set(handles.sliderGreyScaleThreshold, 'Value', log10(figState.fig1.greyThresh));
set(handles.sliderSNSIThreshold, 'Value', figState.fig1.SNSIThresh);
set(handles.sliderTVRFilter, 'Value', figState.fig1.TVRThresh);


% --- Executes on button press in radioFig2.
function radioFig2_Callback(hObject, eventdata, handles)

global figState

set(handles.radioFig2, 'value', 1);
set(handles.radioFig1, 'value', 0);

set(handles.popupmenuSegment, 'Value', figState.fig2.segment);
set(handles.popupmenuSelectSpectrum, 'Value', figState.fig2.spectrum);
set(handles.editARorder, 'String', num2str(figState.fig2.ARorder));
set(handles.popupmenuProcessingFn, 'Value', figState.fig2.process);
set(handles.sliderGreyScaleThreshold, 'Value', log10(figState.fig2.greyThresh));
set(handles.sliderSNSIThreshold, 'Value', figState.fig2.SNSIThresh);
set(handles.sliderTVRFilter, 'Value', figState.fig2.TVRThresh);


% --- Executes on selection change in popupmenuProcessingFn.
function popupmenuProcessingFn_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuProcessingFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuProcessingFn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuProcessingFn
global figState

if get(handles.radioFig1, 'value') == 1
    figState.fig1.process = get(hObject, 'Value');
else
    figState.fig2.process = get(hObject, 'Value');
end


% --- Executes during object creation, after setting all properties.
function popupmenuProcessingFn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuProcessingFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderGreyScaleThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to sliderGreyScaleThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global figState

if get(handles.radioFig1, 'value') == 1
    figState.fig1.greyThresh = 10^(get(hObject, 'Value'));
else
    figState.fig2.greyThresh = 10^(get(hObject, 'Value'));
end

% --- Executes during object creation, after setting all properties.
function sliderGreyScaleThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderGreyScaleThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderSNSIThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to sliderSNSIThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global figState

if get(handles.radioFig1, 'value') == 1
    figState.fig1.SNSIThresh = get(hObject, 'Value');
else
    figState.fig2.SNSIThresh = get(hObject, 'Value');
end


% --- Executes during object creation, after setting all properties.
function sliderSNSIThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderSNSIThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderTVRFilter_Callback(hObject, eventdata, handles)
% hObject    handle to sliderTVRFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global figState

if get(handles.radioFig1, 'value') == 1
    figState.fig1.TVRThresh = get(hObject, 'Value');
else
    figState.fig2.TVRThresh = get(hObject, 'Value');
end


% --- Executes during object creation, after setting all properties.
function sliderTVRFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderTVRFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function editARorder_Callback(hObject, eventdata, handles)
% hObject    handle to editARorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editARorder as text
%        str2double(get(hObject,'String')) returns contents of editARorder as a double
global sonoSession sonoImg figState

arOrder = str2double(get(handles.editARorder, 'String'));
if get(handles.radioFig1, 'value') == 1
    figState.fig1.ARorder = str2double(get(hObject, 'String'));
    sonoSession.gensonofwdrevfn(1, 'burg', arOrder);
else
    figState.fig2.ARorder = str2double(get(hObject, 'String'));
    sonoSession.gensonofwdrevfn(2, 'burg', arOrder);
end

sonoImg.int.AR = sonoSession.sonoWav.int;
sonoImg.ext.AR = sonoSession.sonoWav.ext;



% --- Executes during object creation, after setting all properties.
function editARorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editARorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uipushtoolNewSonogram_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolNewSonogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sonoSession sonoImg figState

sonoImg = [];
figState.fig1.greyThreshInitial = [];
figState.fig2.greyThreshInitial = [];

sonoSession = sono;
sonoSession.msWinLen = 50;
winOffsetFrac = 0.05;
sonoSession.loadAdcFile;
sonoSession.gensonofwdrevfn(1, 'fft');
sonoSession.gensonofwdrevfn(2, 'fft');
sonoImg.int.fft = sonoSession.sonoWav.int;
sonoImg.ext.fft = sonoSession.sonoWav.ext;
sonoSession.gensonofwdrevfn(1, 'burg');
sonoSession.gensonofwdrevfn(2, 'burg');
sonoImg.int.AR = sonoSession.sonoWav.int;
sonoImg.ext.AR = sonoSession.sonoWav.ext;
sonoImg.newLoad = true; % use for just loaded but not processed. On 1st update call, set to false
sonoImg = trimSono(sonoImg, 0, 'fixed');

% auto trim on load
% seed dialog with guess of how many rows to cut
meanSono = mean(sonoImg.int.AR');
[m, I] = max(meanSono);
%sonoImg = trimSono(sonoImg, I-3, 'fixed');


set(handles.textFileName, 'String', ['Pressure :: ' num2str(sonoSession.pressure) ' mmHg          ' sonoSession.fileName]);


% --------------------------------------------------------------------
function uipushtoolSave_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sonoSession sonoImg figState
varsToSave = {'sonoSession', 'sonoImg'};
uisave(varsToSave, 'sono1');


% --------------------------------------------------------------------
function uipushtoolOpenFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtoolOpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sonoSession sonoImg

cwd = pwd;
cd  ~/Downloads/Aarau
uiopen('load')
cd(cwd)
set(handles.textFileName, 'String', sonoSession.fileName);

% --- Executes on button press in pushbuttonUpdate.
function pushbuttonUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global sonoSession sonoImg figState

h(1) = subplot(3,1,1, 'Parent', handles.uipanelPlots);
figState.fig1.ax = h(1);
if get(handles.radiobuttonFixFig1, 'Value')
    fig1String = 'figState.sonoImg';
else
    fig1String = 'sonoImg';
end

if figState.fig1.segment == 1
    fig1String = [fig1String '.' 'int'];
else
    fig1String = [fig1String '.' 'ext'];
end
if figState.fig1.spectrum == 1
    fig1String = [fig1String '.' 'fft'];
else
    fig1String = [fig1String '.' 'AR'];
end

keyboard

if ~sonoImg.newLoad
    ax = gca;
    sonoImg.sessionColorMap = colormap(ax); 
else 
    sonoImg.sessionColorMap = 'jet'; % default
end

imagesc(sonoSession.timeAx, sonoSession.fAxSigned, log10(eval(fig1String)));
xlabel('t (s)')
ylabel('f (Hz)')
%title(sonoArray(a).title);
axis xy
colormap(sonoImg.sessionColorMap);

% envelope
switch figState.fig1.process
    case 1 % binarize
        if isempty(figState.fig1.greyThreshInitial)
            [imBin, binThresh] = sonoSession.binarizeSono(eval(fig1String));
            figState.fig1.greyThreshInitial = binThresh;
            figState.fig1.greyThresh = binThresh;
            set(handles.sliderGreyScaleThreshold, 'Value', log10(binThresh));
            % TODO set(handles.sliderGreyScaleThreshold, 'Min',  );
        else
            imBin = sonoSession.binarizeSono(eval(fig1String), figState.fig1.greyThresh);
        end
        env1 = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        % optimize
        if sonoImg.newLoad
            env1 = optimizeThreshold(handles, eval(fig1String), env1, 1);
        end
    case 2  % modified SNSI
        env1 = sonoSession.getSimpleEnv(eval(fig1String), figState.fig1.SNSIThresh);
end

env1 = TVL1denoise(env1, figState.fig1.TVRThresh, 100);

if ~isempty(env1)
    hold on
    plot(sonoSession.timeAx, env1, 'lineWidth', 2);
    hold off
end

if figState.fig1.segment == 1
    sonoImg.int.env = env1;
else
    sonoImg.ext.env = env1;
end


% fig 2
fig2String = 'sonoImg';

if figState.fig2.segment == 1
    fig2String = [fig2String '.' 'int'];
else
    fig2String = [fig2String '.' 'ext'];
end
if figState.fig2.spectrum == 1
    fig2String = [fig2String '.' 'fft'];
else
    fig2String = [fig2String '.' 'AR'];
end


h(2) = subplot(3,1,2, 'Parent', handles.uipanelPlots);
figState.fig2.ax = h(2);
imagesc(sonoSession.timeAx, sonoSession.fAxSigned, log10(eval(fig2String)));
xlabel('t (s)')
ylabel('f (Hz)')
%title(sonoArray(a).title);
axis xy
colormap(sonoImg.sessionColorMap);

% envelope
switch figState.fig2.process
    case 1 % binarize
        if isempty(figState.fig2.greyThreshInitial)
            [imBin, binThresh] = sonoSession.binarizeSono(eval(fig2String));
            figState.fig2.greyThreshInitial = binThresh;
            figState.fig2.greyThresh = binThresh;
            set(handles.sliderGreyScaleThreshold, 'Value', log10(binThresh));
            % TODO set(handles.sliderGreyScaleThreshold, 'Min',  );
        else
            imBin = sonoSession.binarizeSono(eval(fig2String), figState.fig2.greyThresh);
        end
        env2 = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        % optimize
        if sonoImg.newLoad
            env2 = optimizeThreshold(handles, eval(fig2String), env2, 2);
        end
    case 2  % modified SNSI
        env2 = sonoSession.getSimpleEnv(eval(fig2String), figState.fig2.SNSIThresh);
end

env2 = TVL1denoise(env2, figState.fig2.TVRThresh, 100);
if ~isempty(env2)
    hold on
    plot(sonoSession.timeAx, env2, 'r', 'lineWidth', 2);
    hold off
end

if figState.fig2.segment == 1
    sonoImg.int.env = env2;
else
    sonoImg.ext.env = env2;
end

% matching
h(3) = subplot(3,1,3, 'Parent', handles.uipanelPlots);

matchingMethod = get(handles.popupmenuMatching, 'Value');
switch matchingMethod
    case 1
        matchingEnv1 = env1;
        matchingEnv2 = env2;
    case 2
        % get moment order
        momentOrder = get(handles.popupmenuMomentOrder, 'Value');
        %envSono1 = sonoSession.createEnvSono(eval(fig1String), env1);
        %envSono2 = sonoSession.createEnvSono(eval(fig2String), env1);
        %fig1String = ['log10(' fig1String ')'];
        %fig2String = ['log10(' fig2String ')'];
        matchingEnv1 = sonoSession.spectralMoment(eval(fig1String), env1, momentOrder);
        matchingEnv2 = sonoSession.spectralMoment(eval(fig2String), env2, momentOrder);
end

isNormalize = get(handles.radiobuttonNormalize, 'Value');
if isNormalize
    plot(sonoSession.timeAx, matchingEnv1./max(matchingEnv1));
    hold on;
    plot(sonoSession.timeAx, matchingEnv2./max(matchingEnv2), 'r');
    hold off;
    matchStats = sonoSession.matchWaveforms(matchingEnv1, matchingEnv2, env1, env2, 1); 
else
    plot(sonoSession.timeAx, matchingEnv1);
    hold on;
    plot(sonoSession.timeAx, matchingEnv2, 'r');
    hold off;
    matchStats = sonoSession.matchWaveforms(matchingEnv1, matchingEnv2, env1, env2, 0);
end

linkaxes(h, 'x');

% update matching summary stats in table

tableRow = [matchStats.SSE, matchStats.RMSSD, matchStats.numGoodPeriods, matchStats.PI1, matchStats.PI2];
set(handles.uitableMatchMetrics, 'Data', tableRow);

sonoImg.newLoad = false;

guidata(hObject, handles);

function env = optimizeThreshold(handles, img, env, fig)

global figState sonoSession 

% optimize
numPeaksOld = peakDetect(env,.3);
% Move slider to the left in big increments
%shifts = get(handles.sliderGreyScaleThreshold, 'SliderStep');
shifts = [0.1 0.5];
numPeaksNew = numPeaksOld;

if fig == 1
    while numPeaksNew == numPeaksOld
        figState.fig1.greyThresh = 10^(log10(figState.fig1.greyThresh) - shifts(2));
        imBin = sonoSession.binarizeSono(img, figState.fig1.greyThresh);
        env = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        %  env = TVL1denoise(env, figState.fig1.TVRThresh, 100);
        numPeaksNew = peakDetect(env,.3);
    end
    % move slider to the right in small increments
    numPeaksOld = numPeaksNew+1;
    numPeaksOld1 = numPeaksOld;
    while numPeaksOld > numPeaksNew || numPeaksOld1 > numPeaksOld || numPeaksNew > 30
        figState.fig1.greyThresh = 10^(log10(figState.fig1.greyThresh) + shifts(1));
        imBin = sonoSession.binarizeSono(img, figState.fig1.greyThresh);
        env = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        %    env = TVL1denoise(env, figState.fig1.TVRThresh, 100);
        numPeaksOld1 = numPeaksOld;
        numPeaksOld = numPeaksNew;
        numPeaksNew = peakDetect(env,.3);
    end
    set(handles.sliderGreyScaleThreshold, 'Value', log10(figState.fig1.greyThresh))
else
    while numPeaksNew == numPeaksOld
        figState.fig2.greyThresh = 10^(log10(figState.fig2.greyThresh) - shifts(2));
        imBin = sonoSession.binarizeSono(img, figState.fig2.greyThresh);
        env = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        %  env = TVL1denoise(env, figState.fig1.TVRThresh, 100);
        numPeaksNew = peakDetect(env,.3);
    end
    % move slider to the right in small increments
    % check that its stable for two moves
    numPeaksOld = numPeaksNew+1;
    numPeaksOld1 = numPeaksOld;
    while numPeaksOld > numPeaksNew || numPeaksOld1 > numPeaksOld || numPeaksNew > 30
        figState.fig2.greyThresh = 10^(log10(figState.fig2.greyThresh) + shifts(1));
        imBin = sonoSession.binarizeSono(img, figState.fig2.greyThresh);
        env = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        %    env = TVL1denoise(env, figState.fig1.TVRThresh, 100);
        numPeaksOld1 = numPeaksOld;
        numPeaksOld = numPeaksNew;
        numPeaksNew = peakDetect(env,.3);
    end
    set(handles.sliderGreyScaleThreshold, 'Value', log10(figState.fig2.greyThresh))
end

% --- Executes on button press in radiobuttonNormalize.
function radiobuttonNormalize_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonNormalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonNormalize


% --- Executes on selection change in popupmenuMatching.
function popupmenuMatching_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMatching (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuMatching contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMatching


% --- Executes during object creation, after setting all properties.
function popupmenuMatching_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMatching (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in popupmenuMomentOrder.
function popupmenuMomentOrder_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMomentOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuMomentOrder contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMomentOrder


% --- Executes during object creation, after setting all properties.
function popupmenuMomentOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMomentOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobuttonFixFig1.
function radiobuttonFixFig1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonFixFig1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonFixFig1

global figState sonoImg
if get(hObject, 'Value')
    set(handles.radioFig1, 'Enable', 'off');
    set(handles.popupmenuSegment, 'Value', figState.fig1.segment);
    set(handles.popupmenuSelectSpectrum, 'Value', figState.fig1.spectrum);
    set(handles.editARorder, 'String', num2str(figState.fig1.ARorder));
    set(handles.popupmenuProcessingFn, 'Value', figState.fig1.process);
    set(handles.sliderGreyScaleThreshold, 'Value', figState.fig1.greyThresh);
    set(handles.sliderSNSIThreshold, 'Value', figState.fig1.SNSIThresh);
    set(handles.sliderTVRFilter, 'Value', figState.fig1.TVRThresh);
    figState.sonoImg = sonoImg;
else
    set(handles.radioFig1, 'Enable', 'on');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function submenuSetGreyscaleRange_Callback(hObject, eventdata, handles)

prompt = {'Min greyscale:','Max greyscale:', 'Greyscale Step', 'Greyscale Value'};
dlg_title = 'Greyscale settings';
num_lines = 1;
currMin = num2str(get(handles.sliderGreyScaleThreshold, 'Min'));
currMax = num2str(get(handles.sliderGreyScaleThreshold, 'Max'));
currStep = num2str(get(handles.sliderGreyScaleThreshold, 'SliderStep'));
currValue = num2str(get(handles.sliderGreyScaleThreshold, 'Value'));

defaultans = {currMin, currMax, currStep, currValue};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

newMin = str2double(answer{1});
newMax = str2double(answer{2});
newStep = str2num(answer{3});
newValue = str2double(answer{4});

set(handles.sliderGreyScaleThreshold, 'Min', newMin);
set(handles.sliderGreyScaleThreshold, 'Max', newMax);
set(handles.sliderGreyScaleThreshold, 'SliderStep', newStep);
set(handles.sliderGreyScaleThreshold, 'Value', newValue);

if currValue >= newMin & currValue <= newMax
    % do nothing and leave
else
    newValue = mean([newMin newMax]);
    set(handles.sliderGreyScaleThreshold, 'Value', newValue);
end

guidata(hObject, handles);

% --------------------------------------------------------------------
function submenuSetSNSIRange_Callback(hObject, eventdata, handles)

prompt = {'Min SNSI:','Max SNSI:', 'SNSI Step'};
dlg_title = 'SNSI settings';
num_lines = 1;
currMin = num2str(get(handles.sliderSNSIThreshold, 'Min'));
currMax = num2str(get(handles.sliderSNSIThreshold, 'Max'));
currStep = num2str(get(handles.sliderSNSIThreshold, 'SliderStep'));
currValue = num2str(get(handles.sliderSNSIThreshold, 'Value'));

defaultans = {currMin, currMax, currStep};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

newMin = str2double(answer{1});
newMax = str2double(answer{2});
newStep = str2num(answer{3});

set(handles.sliderSNSIThreshold, 'Min', newMin);
set(handles.sliderSNSIThreshold, 'Max', newMax);
set(handles.sliderSNSIThreshold, 'SliderStep', newStep);

if currValue >= newMin & currValue <= newMax
    % do nothing and leave
else
    newValue = mean([newMin newMax]);
    set(handles.sliderSNSIThreshold, 'Value', newValue);
end

% --------------------------------------------------------------------
function submenuSetTVRRange_Callback(hObject, eventdata, handles)

prompt = {'Min TVR:','Max TVR:', 'TVR STep'};
dlg_title = 'TVR settings';
num_lines = 1;
currMin = num2str(get(handles.sliderTVRFilter, 'Min'));
currMax = num2str(get(handles.sliderTVRFilter, 'Max'));
currStep = num2str(get(handles.sliderTVRFilter, 'SliderStep'));
currValue = num2str(get(handles.sliderTVRFilter, 'Value'));

defaultans = {currMin, currMax, currStep};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

newMin = str2double(answer{1});
newMax = str2double(answer{2});
newStep = str2num(answer{3});

set(handles.sliderTVRFilter, 'Min', newMin);
set(handles.sliderTVRFilter, 'Max', newMax);
set(handles.sliderTVRFilter, 'SliderStep', newStep);

if currValue >= newMin & currValue <= newMax
    % do nothing and leave
else
    newValue = mean([newMin newMax]);
    set(handles.sliderTVRFilter, 'Value', newValue);
end

% --------------------------------------------------------------------
function subContextMenuFixed_Callback(hObject, eventdata, handles)

global sonoImg

prompt = {'Set number of rows to trim'};
dlg_title = 'Trim Sonogram';
num_lines = 1;

% seed dialog with guess of how many rows to cut
meanSono = mean(sonoImg.int.AROriginal');
[m, I] = max(meanSono);

defaultans = {num2str(I - 4)};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

sonoImg = trimSono(sonoImg, str2double(answer), 'fixed');


% --------------------------------------------------------------------
function subContextMenuPrctg_Callback(hObject, eventdata, handles)

global sonoImg

prompt = {'Set percentage of rows to trim'};
dlg_title = 'Trim Sonogram';
num_lines = 1;
defaultans = {'0'};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

sonoImg = trimSono(sonoImg, str2double(answer), 'prctg');


% --------------------------------------------------------------------
function TrimSonogram_Callback(hObject, eventdata, handles)
% hObject    handle to TrimSonogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonAverageSono.
function pushbuttonAverageSono_Callback(hObject, eventdata, handles)

global sonoSession sonoImg figState

fig1Img = 'sonoImg';
fig1Env = 'sonoImg';

if figState.fig1.segment == 1
    fig1Img = [fig1Img '.' 'int'];
    fig1Env = [fig1Env '.' 'int.env'];
else
    fig1Img = [fig1Img '.' 'ext'];
    fig2Env = [fig2Env '.' 'ext.env'];
end
if figState.fig1.spectrum == 1
    fig1Img = [fig1Img '.' 'fft'];
else
    fig1Img = [fig1Img '.' 'AR'];
end
fig1Img = ['flipud(' fig1Img ')'];

bDisplayLines = strcmp(get(handles.uitoggletoolShowLines, 'State'), 'on');
[averagedSono1, averagedEnv1] = sonoSession.averageSono(eval(fig1Img), eval(fig1Env), bDisplayLines);
f = figure;
set(f,'name','Averaged Sonograms','numbertitle','off')
h(1) = subplot(3,1,1);
imagesc(log10(averagedSono1));
colormap(sonoImg.sessionColorMap);

% envelope detection on averaged sono

imBin = sonoSession.binarizeSono(averagedSono1, figState.fig1.greyThresh);
env1 = sonoSession.getEnvelopeFromBinaryImg(1-imBin, 1);
env1 = TVL1denoise(env1, figState.fig1.TVRThresh, 100);

hold on;
plot(env1, 'lineWidth', 2);
plot(averagedEnv1, ':', 'lineWidth', 2);
colormapeditor;
hold off;

% save to sonoImg
if figState.fig1.segment == 1
    sonoImg.int.averagedSono = averagedSono1;
    sonoImg.int.averagedEnv = averagedEnv1;
else
    sonoImg.ext.averagedSono = averagedSono1;
    sonoImg.ext.averagedEnv = averagedEnv1;
end

% second figure
fig2Img = 'sonoImg';
fig2Env = 'sonoImg';
if figState.fig2.segment == 1
    fig2Img = [fig2Img '.' 'int'];
    fig2Env = [fig2Env '.' 'int.env'];
else
    fig2Img = [fig2Img '.' 'ext'];
    fig2Env = [fig2Env '.' 'ext.env'];
end
if figState.fig2.spectrum == 1
    fig2Img = [fig2Img '.' 'fft'];
else
    fig2Img = [fig2Img '.' 'AR'];
end
fig2Img = ['flipud(' fig2Img ')'];

[averagedSono2, averagedEnv2] = sonoSession.averageSono(eval(fig2Img), eval(fig2Env), bDisplayLines);
figure(f);
h(2) = subplot(3,1,2);
imagesc(log10(averagedSono2));
colormap(sonoImg.sessionColorMap);

% envelope detection on averaged sono

imBin = sonoSession.binarizeSono(averagedSono2, figState.fig2.greyThresh);
env2 = sonoSession.getEnvelopeFromBinaryImg(1-imBin, 1);
env2 = TVL1denoise(env2, figState.fig2.TVRThresh, 100);

hold on;
plot(env2, 'r', 'lineWidth', 2);
plot(averagedEnv2, 'r:', 'lineWidth', 2);
hold off;

% save to sonoImg
if figState.fig2.segment == 1
    sonoImg.int.averagedSono = averagedSono2;
    sonoImg.int.averagedEnv = averagedEnv2;
else
    sonoImg.ext.averagedSono = averagedSono2;
    sonoImg.ext.averagedEnv = averagedEnv2;
end

figure(f);
h(3) = subplot(3,1,3);

%env1 = size(averagedSono1,1) - env1;
%env2 = size(averagedSono2,1) - env2;
averagedEnv1 = size(averagedSono1,1) - averagedEnv1;
averagedEnv2 = size(averagedSono2,1) - averagedEnv2;

matchingMethod = get(handles.popupmenuMatching, 'Value');
switch matchingMethod
    case 1
        matchingEnv1 = averagedEnv1;
        matchingEnv2 = averagedEnv2;
    case 2
        % get moment order
        momentOrder = get(handles.popupmenuMomentOrder, 'Value');
        %envSono1 = sonoSession.createEnvSono(eval(fig1String), env1);
        %envSono2 = sonoSession.createEnvSono(eval(fig2String), env1);
        %fig1String = ['log10(' fig1String ')'];
        %fig2String = ['log10(' fig2String ')'];
        matchingEnv1 = sonoSession.spectralMoment(flipud(averagedSono1), averagedEnv1, momentOrder);
        matchingEnv2 = sonoSession.spectralMoment(flipud(averagedSono2), averagedEnv2, momentOrder);
end

isNormalize = get(handles.radiobuttonNormalize, 'Value');
if isNormalize
    plot(matchingEnv1./max(matchingEnv1));
    hold on;
    plot(matchingEnv2./max(matchingEnv2), 'r');
    hold off;

    matchStats = sonoSession.matchWaveforms(matchingEnv1./max(matchingEnv1), matchingEnv2./max(matchingEnv2), averagedEnv1, averagedEnv2);
else
    plot(matchingEnv1);
    hold on;
    plot(matchingEnv2, 'r');
    hold off;
    %matchStats = sonoSession.matchWaveforms(matchingEnv1, matchingEnv2);
end

linkaxes(h, 'x');

% --------------------------------------------------------------------
function subMenuColormap_Callback(hObject, eventdata, handles)

colormapeditor


% --------------------------------------------------------------------
function uipushtoolAppendSave_ClickedCallback(hObject, eventdata, handles)

global sonoSession sonoImg figState

% get folder name
splits = strsplit(sonoSession.fileName, '\');
folderName = splits(9);
folderName = folderName{1};

% get file name
[p,n,e] = fileparts(sonoSession.fileName);
fileNum = str2double(n(5:7));

if exist('studyStore.mat', 'file')
    load('studyStore.mat');
    % find number associated with name
    numFolders = length(studyStore);
    bExists = 0;
    for a = 1:numFolders
        if strmatch(folderName, studyStore{a}.name)
           % studyStore{a}.name = folderName;
            studyStore{a}.data(fileNum) = sonoImg;    
            bExists = 1;
        end
    end
    if ~bExists
        studyStore{numFolders+1}.name = folderName;
        studyStore{numFolders+1}.data(fileNum) = sonoImg;    
    end
else
    % 1st one
    studyStore{1}.name = folderName;
    studyStore{1}.data(fileNum) = sonoImg;
end

save('studyStore.mat', 'studyStore');
msgbox('Study Store Updated', 'Update Store');


% --------------------------------------------------------------------
function uitoggletoolShowLines_ClickedCallback(hObject, eventdata, handles)


% --- Executes on button press in pushbuttonMarkGoodEnv.
function pushbuttonMarkGoodEnv_Callback(hObject, eventdata, handles)

global sonoSession sonoImg figState

balIntObj = balanceFactors(sonoImg.int.AR, sonoImg.int.env);
balIntObj.calculatePeriod();

if ~isempty(balIntObj.cycles.PeriodStart)
    numCycles = length(balIntObj.cycles.PeriodStart);
    for a = 1:numCycles
        startIdx = balIntObj.cycles.PeriodStart(a);
        endIdx = balIntObj.cycles.PeriodEnd(a);
        axes(figState.fig1.ax);
        hold on;
        plot(sonoSession.timeAx(startIdx:endIdx), sonoImg.int.env(startIdx:endIdx), 'k', 'lineWidth', 2);
    end
    hold off;
end

balExtObj = balanceFactors(sonoImg.ext.AR, sonoImg.ext.env);
balExtObj.calculatePeriod();

if ~isempty(balExtObj.cycles.PeriodStart)
    numCycles = length(balExtObj.cycles.PeriodStart);
    for a = 1:numCycles
        startIdx = balExtObj.cycles.PeriodStart(a);
        endIdx = balExtObj.cycles.PeriodEnd(a);
        axes(figState.fig2.ax);
        hold on;
        plot(sonoSession.timeAx(startIdx:endIdx), sonoImg.ext.env(startIdx:endIdx), 'k', 'lineWidth', 2);
    end
    hold off;
end