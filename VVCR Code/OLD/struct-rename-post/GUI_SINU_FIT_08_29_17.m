function varargout = GUI_SINU_FIT_08_29_17 (varargin)
% GUI_SINU_FIT_08_29_17 MATLAB code for GUI_SINU_FIT_08_29_17.fig
%      GUI_SINU_FIT_08_29_17, by itself, creates a new GUI_SINU_FIT_08_29_17 or raises the existing
%      singleton*.
%
%      H = GUI_SINU_FIT_08_29_17 returns the handle to a new GUI_SINU_FIT_08_29_17 or the handle to
%      the existing singleton*.
%
%      GUI_SINU_FIT_08_29_17('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SINU_FIT_08_29_17.M with the given input arguments.
%
%      GUI_SINU_FIT_08_29_17('Property','Value',...) creates a new GUI_SINU_FIT_08_29_17 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_SINU_FIT_08_29_17_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_SINU_FIT_08_29_17_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_SINU_FIT_08_29_17

% Last Modified by GUIDE v2.5 03-Oct-2016 00:36:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SINU_FIT_08_29_17_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SINU_FIT_08_29_17_OutputFcn, ...
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
end

% --- Executes just before GUI_SINU_FIT_08_29_17 is made visible.
function GUI_SINU_FIT_08_29_17_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_SINU_FIT_08_29_17 (see VARARGIN)

% Choose default command line output for GUI_SINU_FIT_08_29_17
handles.output = hObject;

% set the input variable in the global handles environment
% passed as PeakStruct from VVCR_* script
handles.InVar = cell2mat(varargin);

% Extract variables from structure for a more clear workflow
% time and pressure vectors
time = handles.InVar.Time_D; % recall this is the 2x sampled time vector
Pres = handles.InVar.Pres_D; % recall this is the 2x sampled pressure vector

% extract isovolmic points (times). These are structures
isovoltime = handles.InVar.ivTime_D;

% Extract wave fits
waveFit = handles.InVar.fit.BadCyc;

% iso volumic points in array (for plotting)
totIsoTimePoints = handles.InVar.ivPlotTime;
totIsoPresPoints = handles.InVar.ivPlotPres;

% intialize these variables for the undo button
handles.UNDOivPlotTime = [];
handles.UNDOivPlotPres = [];

% EDP - end diastolic pressure
EDP = handles.InVar.Iso1StVal;
P_max2 = handles.InVar.fit.PIsoMax;

% regression constants
c_tot2 = handles.InVar.fit.RCoef;

% initial conditions
IC = handles.InVar.fit.InitIC;

% intialize old variable - used for the undo button
handles.UNDOivPlotTime = [];

% set editable text boxes with ICs
set(handles.Mean_txt, 'String',num2str(IC(1)));
set(handles.Amp_txt, 'String',num2str(IC(2)));
set(handles.Freq_txt, 'String',num2str(IC(3)));
set(handles.Phase_txt, 'String',num2str(IC(4)));

% plot pressure, sinusoid fits
[handles] = gui_sinu_plot (time, Pres, EDP, isovoltime, P_max2, c_tot2, ...
    totIsoTimePoints, totIsoPresPoints, hObject, eventdata, handles); 

% store the wavefit, so output tracks which wave forms did not have a good
% fit
handles.OutVar(1).output = waveFit;

% Store the pmax values
handles.OutVar(2).output = P_max2;

% store regression constants
handles.OutVar(3).output = c_tot2;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_SINU_FIT_08_29_17 wait for user response (see UIRESUME)
uiwait(handles.figure1);
end

% function that executes when user clicks on graph
function GraphCallBack(hObject, eventdata, handles)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% store the current isovolumic points as OldIsoT and OldIsoP
% intialize these variables for the undo button
handles.UNDOivPlotTime = handles.InVar.ivPlotTime;
handles.UNDOivPlotPres = handles.InVar.ivPlotPres;
handles.UNDOivTime_D = handles.InVar.ivTime_D;
handles.UNDOIso1StIdx_D = handles.InVar.Iso1StIdx_D;
handles.UNDOIso2StIdx_D = handles.InVar.Iso2StIdx_D;
handles.UNDOivPres_D = handles.InVar.ivPres_D;
handles.UNDOIso1StVal = handles.InVar.Iso1StVal;
handles.UNDOdPmaxIdx = handles.InVar.dPmaxIdx;
handles.UNDOdPminIdx = handles.InVar.dPminIdx;

% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
disp(['Time: ',num2str(cp(1))]);
disp(['Pressure: ',num2str(cp(2))]);

% obtain variables from InVar Struct for a clear workflow
% time and pressure vectors
time = handles.InVar.Time_D;
Pres = handles.InVar.Pres_D; % recall this is the 2x sampled pressure vector
Oldtime = handles.InVar.Time; % old time vector. 1/2 the points

EDP = handles.InVar.Iso1StVal;

% EDP - end diastolic pressure
EDP_T = handles.InVar.Iso1StIdx_D;
EDP_NT = handles.InVar.Iso2StIdx_D;

% pass the time indexes of minima and maxima
pksT = handles.InVar.dPmaxIdx;
MinIdx = handles.InVar.dPminIdx;

% find which waveform the interval was within. Note the click must be
% between EDP and Negative EDP. the following two lines find (1) all EDP
% times that are smaller than the time point of click (2) all negative EDP
% times that are greater than the time point of click.
WaveNumPosRm = find(time(EDP_T)<cp(1));
WaveNumNegRm = find(time(EDP_NT)>cp(1));

if ~isempty(WaveNumPosRm) && ~isempty(WaveNumNegRm)
    
    % find the common number. the last EDP that is smaller then 
    WaveRm = find(WaveNumPosRm==WaveNumNegRm(1));

    if ~isempty(WaveRm)
        disp(['Wave: ', num2str(WaveRm), ' is being removed']);
        
        % erase wave from isovoltime structure. make new structure with one
        % less row.
        EDP_T(WaveRm) = [];
        EDP_NT(WaveRm) = [];
        EDP(WaveRm) = [];
        pksT(WaveRm) = [];
        MinIdx(WaveRm) = [];

        [DatStr] = data_isovol (EDP_T, EDP_NT, Oldtime, time, Pres, pksT, ...
                                MinIdx, true);
        isovoltime = DatStr.T;
        isovol     = DatStr.P;

        % obtain current ICs
        Mea = str2double(get(handles.Mean_txt,'String'));
        Amp = str2double(get(handles.Amp_txt,'String'));
        Fre = str2double(get(handles.Freq_txt,'String'));
        Pha = str2double(get(handles.Phase_txt,'String'));

        ICS = [Mea Amp Fre Pha];
        [RetVal] = isovol_fit ( isovol, isovoltime, time, Pres, ICS );
        
        % update global handles - some from RetVal, others from above. If
        % the Vanderpool method isn't tripped, then nothing really has changed
        % from the call, so this is a just-in-case...
        handles.InVar.ivTime_D = RetVal.ivTime_D;
        handles.InVar.ivPres_D  = RetVal.ivPres_D;
        handles.InVar.ivPlotTime = RetVal.ivPlotTime;
        handles.InVar.ivPlotPres = RetVal.ivPlotPres;

        handles.InVar.Iso1StIdx_D = EDP_T;
        handles.InVar.Iso2StIdx_D = EDP_NT;
        handles.InVar.Iso1StVal = EDP;
        handles.InVar.dPmaxIdx = pksT;
        handles.InVar.dPminIdx = MinIdx;

        % Plot the results
        [handles] = gui_sinu_plot (time, Pres, EDP, isovoltime, ...
            RetVal.fit.PIsoMax, RetVal.fit.RCoef, RetVal.ivPlotTime, ...
            RetVal.ivPlotPres, hObject, eventdata, handles); 

        % store the wavefit, so output tracks which wave forms did not have a good
        % fit
        handles.OutVar(1).output = RetVal.fit.BadCyc;   % waveFit

        % Store the pmax values
        handles.OutVar(2).output = RetVal.fit.PIsoMax; % P_max2

        % store regression constants
        handles.OutVar(3).output = RetVal.fit.RCoef;   % c_tot2
    end
end

% update global handles
guidata(hObject,handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_SINU_FIT_08_29_17_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.OutVar;

% Destroy the GUI
delete(handles.figure1);
end

% --- Executes on button press in Next.
function Next_Callback(~, ~, handles)
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% call on uiresume so output function executes
uiresume(handles.figure1);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(~, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(handles.figure1); % output function does not execute after this!

% maybe flip some flag to stop execution of further code in VVCR code?
end

% --- Executes during object creation, after setting all properties.
function Mean_txt_CreateFcn(hObject, ~, ~)
% hObject    handle to Mean_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function Mean_txt_Callback(hObject, eventdata, handles)

% when user chnages the phase value and presses enter, evoke calculate
% function
calculate_Callback(hObject, eventdata, handles);

end

% --- Executes during object creation, after setting all properties.
function Amp_txt_CreateFcn(hObject, ~, ~)
% hObject    handle to Amp_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function Amp_txt_Callback(hObject, eventdata, handles)

% when user chnages the phase value and presses enter, evoke calculate
% function
calculate_Callback(hObject, eventdata, handles);

end

% --- Executes during object creation, after setting all properties.
function Freq_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function Freq_txt_Callback(hObject, eventdata, handles)

% when user chnages the phase value and presses enter, evoke calculate
% function
calculate_Callback(hObject, eventdata, handles);

end

function Phase_txt_Callback(hObject, eventdata, handles)

% when user chnages the phase value and presses enter, evoke calculate
% function
calculate_Callback(hObject, eventdata, handles);

end
% --- Executes during object creation, after setting all properties.
function Phase_txt_CreateFcn(hObject, ~, ~)
% hObject    handle to Phase_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% calculate sinusoids based on new ICs!
%%% obtain new ICs
Mea = str2double(get(handles.Mean_txt,'String'));
Amp = str2double(get(handles.Amp_txt,'String'));
Fre = str2double(get(handles.Freq_txt,'String'));
Pha = str2double(get(handles.Phase_txt,'String'));

%%% obtain variables from InVar Struct for a clear workflow
% time and pressure vectors
time = handles.InVar.Time_D;
Pres = handles.InVar.Pres_D;
Oldtime = handles.InVar.Time;

% extract isovolmic points (times). These are structures
isovoltime = handles.InVar.ivTime_D;

% extract isovolmic points (pressures). These are structures
isovol = handles.InVar.ivPres_D;

% iso volumic points in array (for plotting)
totIsoTimePoints = handles.InVar.ivPlotTime;
totIsoPresPoints = handles.InVar.ivPlotPres;

% EDP - end diastolic pressure
EDP = handles.InVar.Iso1StVal;

pksT = handles.InVar.dPmaxIdx;
MinIdx = handles.InVar.dPminIdx;

ICS = [Mea Amp Fre Pha];

[RetVal] = isovol_fit ( isovol, isovoltime, time, Pres, ICS );

% Update isovolumic points and global plotting vectors after return
% from isovol_fit.

handles.InVar.ivTime_D = RetVal.ivTime_D;
handles.InVar.ivPres_D  = RetVal.ivPres_D;
handles.InVar.ivPlotTime = RetVal.ivPlotTime;
handles.InVar.ivPlotPres = RetVal.ivPlotPres;

c_tot2 = RetVal.fit.RCoef;
P_max2 = RetVal.fit.PIsoMax;

[handles] = gui_sinu_plot (time, Pres, EDP, isovoltime, P_max2, c_tot2, ...
     totIsoTimePoints, totIsoPresPoints, hObject, eventdata, handles); 

% store the wavefit, so output tracks which wave forms did not have a good
% fit, the pmax values, and the regression constants. Note the latter two have
% been obtained from the RetVal structure above already.
handles.OutVar(1).output = RetVal.fit.BadCyc;
handles.OutVar(2).output = P_max2;
handles.OutVar(3).output = c_tot2;

% Update handles structure
guidata(hObject, handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end

% --- Executes on button press in Exit.
function Exit_Callback(hObject, ~, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% keep in mind when the exit button is pressed, the current
% patient, i, will not be evaluated
                
% set outputs to false
handles.OutVar(1).output = false;
handles.OutVar(2).output = false;

% update handles globally
guidata(hObject, handles);

% call on uiresume so output function executes
uiresume(handles.figure1);
end

% --- Executes on button press in Discard.
function Discard_Callback(hObject, ~, handles)
% hObject    handle to Discard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set outputs to true, indicating Discard button
handles.OutVar(1).output = true;
handles.OutVar(2).output = true;

% update handles globally
guidata(hObject, handles)

% call on uiresume so output function executes
uiresume(handles.figure1);
end

% --- Executes on button press in Undo.
function Undo_Callback(hObject, eventdata, handles)
% hObject    handle to Undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% if the handles.Old variables have been created (user has clicked on the
% plot and removed pressure waveform(s).
if ~isempty(handles.UNDOivPlotTime)
    % restore the old isovolumic points 
    handles.InVar.ivPlotTime = handles.UNDOivPlotTime;
    handles.InVar.ivPlotPres = handles.UNDOivPlotPres;
    handles.InVar.ivTime_D = handles.UNDOivTime_D;
    handles.InVar.Iso1StIdx_D = handles.UNDOIso1StIdx_D;
    handles.InVar.Iso2StIdx_D = handles.UNDOIso2StIdx_D;
    handles.InVar.ivPres_D = handles.UNDOivPres_D;
    handles.InVar.Iso1StVal = handles.UNDOIso1StVal;
    handles.InVar.dPmaxIdx = handles.UNDOdPmaxIdx;
    handles.InVar.dPminIdx = handles.UNDOdPminIdx;
    
    % calculate sinusoids based on current ICs!
    Mea = str2double(get(handles.Mean_txt,'String'));
    Amp = str2double(get(handles.Amp_txt,'String'));
    Fre = str2double(get(handles.Freq_txt,'String'));
    Pha = str2double(get(handles.Phase_txt,'String'));

    % obtain variables from InVar Struct for a clear workflow
    % this pulls data out of PeakStruct as passed to GUI_SINU_FIT
    time = handles.InVar.Time_D;
    Pres = handles.InVar.Pres_D;
    isovoltime = handles.InVar.ivTime_D;
    isovol = handles.InVar.ivPres_D;
    % iso volumic points in array (for plotting)
    totIsoTimePoints = handles.InVar.ivPlotTime;
    totIsoPresPoints = handles.InVar.ivPlotPres;
    % EDP - end diastolic pressure
    EDP = handles.InVar.Iso1StVal;

    % pre - allocate 
    ICS = [Mea Amp Fre Pha];
    [RetVal] = isovol_fit ( isovol, isovoltime, time, Pres, ICS );

    [handles] = gui_sinu_plot (time, Pres, EDP, isovoltime, ...
            RetVal.fit.PIsoMax, RetVal.fit.RCoef, RetVal.ivPlotTime, ...
            RetVal.ivPlotPres, hObject, eventdata, handles);

    % store the wavefit, so output tracks which wave forms did not have a good
    % fit
    handles.OutVar(1).output = RetVal.fit.BadCyc;   % waveFit

    % Store the pmax values
    handles.OutVar(2).output = RetVal.fit.PIsoMax; % P_max2

    % store regression constants
    handles.OutVar(3).output = RetVal.fit.RCoef;   % c_tot2
end

% Update handles structure
guidata(hObject, handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end

% --- Function that updates the main plot
function [handles] = gui_sinu_plot (time, Pres, EDP, isovoltime, P_max2, ...
                         c_tot2, totIsoTimePoints, totIsoPresPoints, ...
                         hObject, eventdata, handles);

axes(handles.pressure_axes);

h = plot(time,Pres,'b', ...
         totIsoTimePoints,totIsoPresPoints,'ro');
set(h, 'HitTest', 'off');

set(handles.pressure_axes,'ButtonDownFcn', ...
    @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',12);

title('Sinusoidal Fitting','FontSize',20);
xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);

hold on;

PmaxT = zeros(length(EDP),1);

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:length(EDP)

    % obtain the range of time of each peak
    interval =  ...
      time(isovoltime(i).PosIso(1,1)):0.002:time(isovoltime(i).NegIso(end,1));

    % plug into Naeiji equation that was just solved for
    FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + ...
      c_tot2(i,4));

    % find time point corresponding to Pmax
    [~, Idx] = min(abs(FitSinePres-P_max2(i)));

    PmaxT(i) = interval(Idx);

    plot(interval, FitSinePres, 'k--', PmaxT(i), P_max2(i), 'go');
    hold on;
end

% check the range of pressure values of Pmax. if the max p_max value is
% over 450, rescale y axis to (0, 300), so individual waveforms can be seen
if max(P_max2) > 450
    ylim([0, 300]);
end

legend('Pressure', 'Isovolumic Points',  'Sinusoid Fit','Pmax', ...
    'Location','southoutside', 'Orientation', 'horizontal');

box on;
grid on;
hold off;
end
