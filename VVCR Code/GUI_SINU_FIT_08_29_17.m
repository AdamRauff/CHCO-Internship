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

% extract isovolumic points (times). These are structures
isovolTime = handles.InVar.ivTime_D;

% Extract wave fits
waveFit = handles.InVar.fit.BadCyc;

% iso volumic points in array (for plotting)
totIsoTimePoints = handles.InVar.ivPlotTime;
totIsoPresPoints = handles.InVar.ivPlotPres;

% intialize these variables for the undo button
handles.UNDOivPlotTime = [];
handles.UNDOivPlotPres = [];

% EDP - end diastolic pressure
Iso1StVal = handles.InVar.Iso1StVal;

% regression constants
myfit = handles.InVar.fit;  

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
[handles] = gui_sinu_plot (time, Pres, Iso1StVal, isovolTime, myfit, ...
    totIsoTimePoints, totIsoPresPoints, hObject, eventdata, handles); 

% store first fit output into output structure; the output structure is
% essentially the entire fit struct created by isovol_fit.
handles.OutVar = myfit;

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
disp('GUI_SINU_FIT>GraphCallBack:');
disp(['    Time: ',num2str(cp(1))]);
disp(['    Pressure: ',num2str(cp(2))]);

% obtain variables from InVar Struct for a clear workflow
% time and pressure vectors
time = handles.InVar.Time_D;
Pres = handles.InVar.Pres_D; % recall this is the 2x sampled pressure vector
Oldtime = handles.InVar.Time; % old time vector. 1/2 the points

Iso1StVal = handles.InVar.Iso1StVal;

% EDP - end diastolic pressure
Iso1StIdx = handles.InVar.Iso1StIdx_D;
Iso2StIdx = handles.InVar.Iso2StIdx_D;

% pass the time indexes of minima and maxima
dPmaxIdx = handles.InVar.dPmaxIdx;
dPminIdx = handles.InVar.dPminIdx;

% find which waveform the interval was within. Note the click must be
% between EDP and Negative EDP. the following two lines find (1) all EDP
% times that are smaller than the time point of click (2) all negative EDP
% times that are greater than the time point of click.
WaveNumPosRm = find(time(Iso1StIdx)<cp(1));
WaveNumNegRm = find(time(Iso2StIdx)>cp(1));

if ~isempty(WaveNumPosRm) && ~isempty(WaveNumNegRm)
    
    % find the common number. the last EDP that is smaller then 
    WaveRm = find(WaveNumPosRm==WaveNumNegRm(1));

    if ~isempty(WaveRm)
        disp(['    Wave: ', num2str(WaveRm), ' is being removed']);
        
        % erase wave from isovolTime structure. make new structure with one
        % less row.
        Iso1StIdx(WaveRm) = [];
        Iso2StIdx(WaveRm) = [];
        Iso1StVal(WaveRm) = [];
        dPmaxIdx(WaveRm) = [];
        dPminIdx(WaveRm) = [];

        [DatStr] = data_isovol (Iso1StIdx, Iso2StIdx, Oldtime, time, Pres, ...
                                dPmaxIdx, dPminIdx, true);
                               
        isovolTime = DatStr.T;
        isovolPres = DatStr.P;

        % Store changes
        handles.InVar.Iso1StIdx_D = Iso1StIdx;
        handles.InVar.Iso2StIdx_D = Iso2StIdx;
        handles.InVar.Iso1StVal = Iso1StVal;
        handles.InVar.dPmaxIdx = dPmaxIdx;
        handles.InVar.dPminIdx = dPminIdx;

        % obtain current ICs
        Mea = str2double(get(handles.Mean_txt,'String'));
        Amp = str2double(get(handles.Amp_txt,'String'));
        Fre = str2double(get(handles.Freq_txt,'String'));
        Pha = str2double(get(handles.Phase_txt,'String'));

        ICS = [Mea Amp Fre Pha];
        [FitStr] = isovol_fit (isovolPres, isovolTime, time, Pres, ICS);
        
        % update global handles - some from FitStr, others from above. If
        % the Vanderpool method isn't tripped, then nothing really has changed
        % from the call, so this is a just-in-case...
        handles.InVar.ivTime_D = FitStr.ivTime_D;
        handles.InVar.ivPres_D = FitStr.ivPres_D;
        handles.InVar.ivPlotTime = FitStr.ivPlotTime;
        handles.InVar.ivPlotPres = FitStr.ivPlotPres;

        % store fit output into output structure
        handles.OutVar = FitStr.fit;

        % Plot the results
        [handles] = gui_sinu_plot (time, Pres, Iso1StVal, isovolTime, ...
            FitStr.fit, FitStr.ivPlotTime, FitStr.ivPlotPres, hObject, ...
            eventdata, handles); 
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

% extract isovolumic points (times). These are structures
isovolTime = handles.InVar.ivTime_D;

% extract isovolumic points (pressures). These are structures
isovolPres = handles.InVar.ivPres_D;

% iso volumic points in array (for plotting)
totIsoTimePoints = handles.InVar.ivPlotTime;
totIsoPresPoints = handles.InVar.ivPlotPres;

% EDP - end diastolic pressure
Iso1StVal = handles.InVar.Iso1StVal;

dPmaxIdx = handles.InVar.dPmaxIdx;
dPminIdx = handles.InVar.dPminIdx;

ICS = [Mea Amp Fre Pha];

[FitStr] = isovol_fit ( isovolPres, isovolTime, time, Pres, ICS );

% Update isovolumic points and global plotting vectors after return
% from isovol_fit.

handles.InVar.ivTime_D = FitStr.ivTime_D;
handles.InVar.ivPres_D = FitStr.ivPres_D;
handles.InVar.ivPlotTime = FitStr.ivPlotTime;
handles.InVar.ivPlotPres = FitStr.ivPlotPres;

[handles] = gui_sinu_plot (time, Pres, Iso1StVal, isovolTime, FitStr.fit, ...
     totIsoTimePoints, totIsoPresPoints, hObject, eventdata, handles); 

% store fit output into output structure
handles.OutVar = FitStr.fit;

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
                
% set output to false
handles.OutVar = false;

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
handles.OutVar = true;

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
    isovolTime = handles.InVar.ivTime_D;
    isovolPres = handles.InVar.ivPres_D;
    % iso volumic points in array (for plotting)
    totIsoTimePoints = handles.InVar.ivPlotTime;
    totIsoPresPoints = handles.InVar.ivPlotPres;
    % Iso1StVal - end diastolic pressure
    Iso1StVal = handles.InVar.Iso1StVal;

    % pre - allocate 
    ICS = [Mea Amp Fre Pha];
    [FitStr] = isovol_fit ( isovolPres, isovolTime, time, Pres, ICS );

    [handles] = gui_sinu_plot (time, Pres, Iso1StVal, isovolTime, ...
            FitStr.fit, FitStr.ivPlotTime, FitStr.ivPlotPres, hObject, ...
            eventdata, handles);

    % store fit output into output structure
    handles.OutVar = FitStr.fit;
end

% Update handles structure
guidata(hObject, handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end

% --- Function that updates the main plot
function [handles] = gui_sinu_plot (time, Pres, Iso1StVal, isovolTime, fit, ...
                         totIsoTimePoints, totIsoPresPoints, hObject, ...
                         eventdata, handles);

axes(handles.pressure_axes);

c_tot2 = fit.RCoef;
P_max2 = fit.PIsoMax;

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

PmaxT = zeros(length(Iso1StVal),1);

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:length(Iso1StVal)

    % obtain the range of time of each peak
    interval =  ...
      time(isovolTime(i).PosIso(1,1)):0.002:time(isovolTime(i).NegIso(end,1));

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
