function varargout = GUI_SINU_FIT_08_29_17 (varargin)
% GUI_SINU_FIT_08_29_17 MATLAB code for GUI_SINU_FIT_08_29_17.fig
%      GUI_SINU_FIT_08_29_17, by itself, creates a new GUI_SINU_FIT_08_29_17 
%      or raises the existing singleton*.
%
%      H = GUI_SINU_FIT_08_29_17 returns the handle to a new GUI_SINU_FIT_-
%      08_29_17 or the handle to the existing singleton*.
%
%      GUI_SINU_FIT_08_29_17('CALLBACK',hObject,eventData,handles,...) calls
%      the local function named CALLBACK in GUI_SINU_FIT_08_29_17.M with the 
%      given input arguments.
%
%      GUI_SINU_FIT_08_29_17('Property','Value',...) creates a new GUI_SINU_-
%      FIT_08_29_17 or raises the existing singleton*.  Starting from the left,
%      property value pairs are applied to the GUI before GUI_SINU_FIT_08_29_-
%      17_OpeningFcn gets called.  An unrecognized property name or invalid 
%      value makes property application stop.  All inputs are passed to GUI_-
%      SINU_FIT_08_29_17_OpeningFcn via varargin.
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

% Extract Data, Indices/Values, and Fit Segments from passed structures.
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivVal = handles.InVar.ivVal;
ivSeg = handles.InVar.ivSeg;
myfit = handles.InVar.Fit;  

% Initialize UNDO structure.
handles.UNDO.Fit = [];

% store first fit output into output structure.
handles.OutVar = myfit;

% set editable text boxes with ICs
IC = myfit.InitIC;
set(handles.Mean_txt, 'String',num2str(IC(1)));
set(handles.Amp_txt, 'String',num2str(IC(2)));
set(handles.Freq_txt, 'String',num2str(IC(3)));
set(handles.Phase_txt, 'String',num2str(IC(4)));

% plot pressure, sinusoid fits
[handles] = gui_sinu_plot (Data, ivVal, ivSeg, myfit, Plot, handles);

% Update handles.
guidata(hObject, handles);

% UIWAIT makes GUI_SINU_FIT_08_29_17 wait for user response (see UIRESUME)
uiwait(handles.figure1);
end

% function that executes when user clicks on graph
function GraphCallBack(hObject, eventdata, handles)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% obtain variables from InVar Struct for a clear workflow
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivIdx = handles.InVar.ivIdx;
ivVal = handles.InVar.ivVal;
ivSeg = handles.InVar.ivSeg;

% store the current structures in UNDO structure for the undo button.
handles.UNDO.Fit  = handles.OutVar;
handles.UNDO.Plot = Plot;
handles.UNDO.ivIdx = ivIdx;
handles.UNDO.ivVal = ivVal;
handles.UNDO.ivSeg = ivSeg;
       
% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
disp('GUI_SINU_FIT>GraphCallBack:');
disp(['    Time:     ',num2str(cp(1))]);
disp(['    Pressure: ',num2str(cp(2))]);

% find which waveform the interval was within. Note the click must be
% between EDP and Negative EDP. the following two lines find (1) all EDP
% times that are smaller than the time point of click (2) all negative EDP
% times that are greater than the time point of click.
WaveNumPosRm = find(Data.Time_D(ivIdx.Ps1_D)<cp(1));
WaveNumNegRm = find(Data.Time_D(ivIdx.Ne1_D)>cp(1));

if ~isempty(WaveNumPosRm) && ~isempty(WaveNumNegRm)
    
    % find the common number. the last EDP that is smaller then 
    WaveRm = find(WaveNumPosRm==WaveNumNegRm(1));

    if ~isempty(WaveRm)
        disp(['    Wave: ', num2str(WaveRm), ' is being removed']);
        
        % erase wave from iv1Time structure. make new structure with one
        % less row. THIS KILLS EVERYTHING IN ivIdx & ivVal, NOT JUST THE
        % ORIGINAL FIVE. UPDATE FOR KIND METHOD LATER.
        ivIdx.Ps1(WaveRm)   = [];
        ivIdx.Ne1(WaveRm)   = [];
        ivIdx.Ps1_D(WaveRm) = [];
        ivIdx.Ne1_D(WaveRm) = [];
        ivIdx.dPmax(WaveRm) = [];
        ivIdx.dPmin(WaveRm) = [];

        ivVal.Ps1(WaveRm)   = [];
        ivVal.Ne1(WaveRm)   = [];
        ivVal.dPmax(WaveRm) = [];
        ivVal.dPmin(WaveRm) = [];

        % Store changes
        handles.InVar.ivIdx = ivIdx;
        handles.InVar.ivVal = ivVal;

        % obtain current ICs
        Mea = str2double(get(handles.Mean_txt,'String'));
        Amp = str2double(get(handles.Amp_txt,'String'));
        Fre = str2double(get(handles.Freq_txt,'String'));
        Pha = str2double(get(handles.Phase_txt,'String'));

        % Recompute the segments for this new set of IV indicies
        [ivSeg] = data_isoseg (true, Data, ivIdx);

        ICS = [Mea Amp Fre Pha];
        [Fit, ivSeg, Plot] = isovol_fit (ivSeg, Data, ICS);
        
        % update global handles from isovol_returned values. If the Vanderpool
        % method isn't tripped, then ivSeg and Plot haven't changed, so this
        % is a just-in-case...
        handles.OutVar = Fit;

        handles.InVar.ivSeg = ivSeg;
        handles.IvVar.Plot  = Plot;

        % Plot the results
        [handles] = gui_sinu_plot (Data, ivVal, ivSeg, Fit, Plot, handles);

    end
end

% update global handles & set cursor back to normal
guidata(hObject,handles);
set(handles.figure1, 'pointer', 'arrow');

end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_SINU_FIT_08_29_17_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.OutVar;

% Destroy the GUI
delete(hObject);

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
function figure1_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, call UIRESUME
    uiresume(hObject);
 
    % If you close the figure, we understand that as stopping the analysis.
    handles.OutVar = false;
    guidata(hObject, handles);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

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
disp('GUI_SINU_FIT>calculate_Callback:');

% calculate sinusoids based on new ICs!
%%% obtain new ICs
Mea = str2double(get(handles.Mean_txt,'String'));
Amp = str2double(get(handles.Amp_txt,'String'));
Fre = str2double(get(handles.Freq_txt,'String'));
Pha = str2double(get(handles.Phase_txt,'String'));

% Extract Data, Indices/Values, and Fit Segments from passed structures.
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivIdx = handles.InVar.ivIdx;
ivVal = handles.InVar.ivVal;
ivSeg = handles.InVar.ivSeg;

ICS = [Mea Amp Fre Pha];

[Fit, ivSeg, Plot] = isovol_fit (ivSeg, Data, ICS );

% update global handles from isovol_returned values. If the Vanderpool
% method isn't tripped, then ivSeg and Plot haven't changed, so this
% is a just-in-case...
handles.OutVar = Fit;

handles.InVar.ivSeg = ivSeg;
handles.IvVar.Plot  = Plot;

[handles] = gui_sinu_plot (Data, ivVal, ivSeg, Fit, Plot, handles);

% update global handles & set cursor back to normal
guidata(hObject,handles);
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
if ~isempty(handles.UNDO.Fit)
    disp('GUI_SINU_FIT>Undo_Callback: Restoring Previous Fit & Plot');
    handles.OutVar = handles.UNDO.Fit;

    handles.InVar.Plot  = handles.UNDO.Plot;
    handles.InVar.ivIdx = handles.UNDO.ivIdx;
    handles.InVar.ivVal = handles.UNDO.ivVal;
    handles.InVar.ivSeg = handles.UNDO.ivSeg; 
    
    % Extract Data, Indices/Values, and Fit Segments from passed structures.
    Fit  = handles.OutVar;
    Data = handles.InVar.Data;
    Plot = handles.InVar.Plot;
    ivIdx = handles.InVar.ivIdx;
    ivVal = handles.InVar.ivVal;
    ivSeg = handles.InVar.ivSeg;

    [handles] = gui_sinu_plot (Data, ivVal, ivSeg, Fit, Plot, handles);

else

    disp('GUI_SINU_FIT>Undo_Callback: Nothing to Undo!');

end

% update global handles & set cursor back to normal
guidata(hObject,handles);
set(handles.figure1, 'pointer', 'arrow');

end

% --- Function that updates the main plot
function [handles] = gui_sinu_plot (Data, ivVal, ivSeg, Fit, Plot, handles);

axes(handles.pressure_axes);

h = plot(Data.Time_D,Data.Pres_D,'b', ...
         Plot.ivPlotTime,Plot.ivPlotPres,'ro');
set(h, 'HitTest', 'off');

set(handles.pressure_axes,'ButtonDownFcn', ...
    @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',12);

title('Sinusoidal Fitting','FontSize',20);
xlabel('Time [s]','FontSize',18);
ylabel('Data.Pres_Dsue [mmHg]','FontSize',18);

hold on;

mysz = length(ivVal.Ps1);
PmaxT = zeros(mysz,1);

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:mysz

    % obtain the range of time of each peak
    interval = Data.Time_D(ivSeg.iv1Time(i).PosIso(1,1)):0.002: ...
        Data.Time_D(ivSeg.iv1Time(i).NegIso(end,1));

    % plug into Naeiji equation that was just solved for
    FitSinePres = Fit.RCoef(i,1) + Fit.RCoef(i,2)*sin(Fit.RCoef(i,3)* ...
      interval + Fit.RCoef(i,4));

    % find time point corresponding to Pmax
    [~, Idx] = min(abs(FitSinePres-Fit.PIsoMax(i)));

    PmaxT(i) = interval(Idx);

    plot(interval, FitSinePres, 'k--', PmaxT(i), Fit.PIsoMax(i), 'go');
    hold on;
end

% check the range of pressure values of Pmax. if the max p_max value is
% over 450, rescale y axis to (0, 300), so individual waveforms can be seen
if max(Fit.PIsoMax) > 450
    ylim([0, 300]);
end

legend('Pressure', 'Isovolumic Points',  'Sinusoid Fit','Pmax', ...
    'Location','southoutside', 'Orientation', 'horizontal');

box on;
grid on;
hold off;

end
