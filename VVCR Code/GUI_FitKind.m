function varargout = GUI_FitKind (varargin)
% GUI_FitKind MATLAB code for GUI_FitKind.fig
%      GUI_FitKind, by itself, creates a new GUI_FitKind or raises the
%      existing singleton*.
%
%      H = GUI_FitKind returns the handle to a new GUI_FitKind or the
%      handle to the existing singleton*.
%
%      GUI_FitKind('CALLBACK',hObject,eventData,handles,...) calls
%      the local function named CALLBACK in GUI_FitKind.M with the 
%      given input arguments.
%
%      GUI_FitKind('Property','Value',...) creates a new GUI_FitKind or
%      raises the existing singleton*.  Starting from the left, property
%      value pairs are applied to the GUI before GUI_FitKind_OpeningFcn
%      gets called.  An unrecognized property name or invalid value 
%      makes property application stop.  All inputs are passed to GUI_-
%      FitKind_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only
%      one instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FitKind

% Last Modified by GUIDE v2.5 20-Sep-2017 13:31:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FitKind_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FitKind_OutputFcn, ...
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

% --- Executes just before GUI_FitKind is made visible.
function GUI_FitKind_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FitKind (see VARARGIN)

% Choose default command line output for GUI_FitKind
handles.output = hObject;

% set the input variable in the global handles environment
% passed as PeakStruct from VVCR_* script
handles.InVar = cell2mat(varargin(1));
GUIDat = cell2mat(varargin(2));
handles.InVar.Data  = GUIDat.Data;
handles.InVar.ivIdx = GUIDat.ivIdx;
handles.InVar.ivVal = GUIDat.ivVal;
handles.InVar.ivSeg = GUIDat.ivSeg;

handles.InVar.MeanTP = handles.InVar.FitK.MeanTP;

% Extract Data, Indices/Values, and Fit Segments from passed structures.
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivIdx = handles.InVar.ivIdx;
ivSeg = handles.InVar.ivSeg;
FitK = handles.InVar.FitK;

% Initialize UNDO structure.
handles.UNDO.FitK = [];

% store first fit output into output structure.
Res.FitK = handles.InVar.FitK;
handles.OutVar = Res;

% set editable text boxes with ICs
%IC = FitK.InitIC;
%set(handles.Mean_txt, 'String',num2str(IC(1)));
%set(handles.Amp_txt, 'String',num2str(IC(2)));
%set(handles.Freq_txt, 'String',num2str(IC(3)));
%set(handles.Phase_txt, 'String',num2str(IC(4)));

% plot pressure, sinusoid fits
[handles] = gui_kind_plot (Data, ivIdx, ivSeg, FitK, Plot, handles);

% Update handles.
guidata(hObject, handles);

% UIWAIT makes GUI_FitKind wait for user response (see UIRESUME)
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

MeanTP = handles.InVar.MeanTP;

% store the current structures in UNDO structure for the undo button.
handles.UNDO.Res  = handles.OutVar;
handles.UNDO.Plot = Plot;
handles.UNDO.ivIdx = ivIdx;
handles.UNDO.ivVal = ivVal;
handles.UNDO.ivSeg = ivSeg;
       
% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
disp('GUI_FitKind>GraphCallBack:');
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
        
        % Erase wave from (2 - Kind) ivIdx, ivVal structures.
        ivIdx.Ps2(WaveRm)   = [];
        ivIdx.Pe2(WaveRm)   = [];
        ivIdx.Ns2(WaveRm)   = [];
        ivIdx.Ne2(WaveRm)   = [];
        ivIdx.Ps2_D(WaveRm) = [];
        ivIdx.Pe2_D(WaveRm) = [];
        ivIdx.Ns2_D(WaveRm) = [];
        ivIdx.Ne2_D(WaveRm) = [];
        ivVal.Ps2(WaveRm)   = [];
        ivVal.Pe2(WaveRm)   = [];
        ivVal.Ns2(WaveRm)   = [];
        ivVal.Ne2(WaveRm)   = [];

        ivIdx.dPmax2(WaveRm)   = [];
        ivIdx.dPmin2(WaveRm)   = [];
        ivIdx.dPmin2_D(WaveRm) = [];
        ivVal.dPmax2(WaveRm)   = [];
        ivVal.dPmin2(WaveRm)   = [];

        % Store changes
        handles.InVar.ivIdx = ivIdx;
        handles.InVar.ivVal = ivVal;

        % obtain current ICs
        %Mea = str2double(get(handles.Mean_txt,'String'));
        %Amp = str2double(get(handles.Amp_txt,'String'));
        %Fre = str2double(get(handles.Freq_txt,'String'));
        %Pha = str2double(get(handles.Phase_txt,'String'));

        % Recompute the segments for this new set of IV indicies
        [ivSeg] = data_isoseg (true, Data, ivIdx);

        [FitK, PlotK] = fit_kind (ivSeg, ivIdx, Data, MeanTP);
        
        Res.FitK = FitK;
        handles.OutVar = Res;

        handles.IvVar.ivSeg = ivSeg;
        handles.IvVar.Plot  = PlotK;

        % Plot the results
        [handles] = gui_kind_plot (Data, ivIdx, ivSeg, FitK, Plot, handles);

    end
end

% update global handles & set cursor back to normal
guidata(hObject,handles);
set(handles.figure1, 'pointer', 'arrow');

end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_FitKind_OutputFcn(hObject, ~, handles) 
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
function calculate_Callback(hObject, ~, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;
disp('GUI_FitKind>calculate_Callback:');

% calculate sinusoids based on new ICs!
%Mea = str2double(get(handles.Mean_txt,'String'));
%Amp = str2double(get(handles.Amp_txt,'String'));
%Fre = str2double(get(handles.Freq_txt,'String'));
%Pha = str2double(get(handles.Phase_txt,'String'));
%ICS = [Mea Amp Fre Pha];

% Extract Data, Values, and Fit Segments from passed structures; get fits.
Data = handles.InVar.Data;
ivIdx = handles.InVar.ivIdx;
ivSeg = handles.InVar.ivSeg;

MeanTP = handles.InVar.MeanTP;

[FitK, PlotK] = fit_kind (ivSeg, ivIdx, Data, MeanTP);

Res.FitK = FitK;
handles.OutVar = Res;
handles.IvVar.Plot = PlotK;

[handles] = gui_kind_plot (Data, ivIdx, ivSeg, FitK, Plot, handles);

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
function Undo_Callback(hObject, ~, handles)
% hObject    handle to Undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% if the handles.Old variables have been created (user has clicked on the
% plot and removed pressure waveform(s).
if ~isempty(handles.UNDO.Res)
    disp('GUI_FitKind>Undo_Callback: Restoring Previous Fit & Plot');
    handles.OutVar = handles.UNDO.Res;

    handles.InVar.Plot  = handles.UNDO.Plot;
    handles.InVar.ivIdx = handles.UNDO.ivIdx;
    handles.InVar.ivVal = handles.UNDO.ivVal;
    handles.InVar.ivSeg = handles.UNDO.ivSeg; 
    
    % Extract Data, Values, Fit Segments, Plots, & Segments from handles.
    FitK = handles.OutVar.FitK;
    Data = handles.InVar.Data;
    Plot = handles.InVar.Plot;
    ivIdx = handles.InVar.ivIdx;
    ivSeg = handles.InVar.ivSeg;

    [handles] = gui_kind_plot (Data, ivIdx, ivSeg, FitK, Plot, handles);

else

    disp('GUI_FitKind>Undo_Callback: Nothing to Undo!');

end

% update global handles & set cursor back to normal
guidata(hObject,handles);
set(handles.figure1, 'pointer', 'arrow');

end

% --- Function that updates the main plot
function [handles] = gui_kind_plot (Data, ivIdx, ivSeg, Fit, Plot, handles);

axes(handles.pressure_axes);

h = plot(Data.Time_D,Data.Pres_D,'b', ...
         Plot.iv2PlotTime,Plot.iv2PlotPres,'ro');
set(h, 'HitTest', 'off');

set(handles.pressure_axes,'ButtonDownFcn', ...
    @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',12);

title('Sinusoidal Fitting','FontSize',20);
xlabel('Time [s]','FontSize',18);
ylabel('Data.Pres_Dsue [mmHg]','FontSize',18);

hold on;

mystp = Data.time_step/2;
mysz = length(ivSeg.iv2Time);
PmaxT = zeros(mysz,1);

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:mysz

    % obtain the range of time of each peak, then normalize to zero
    FitSineTime = Data.Time_D(ivSeg.iv2Time(i).PosIso(1,1)):mystp: ...
        Data.Time_D(ivSeg.iv2Time(i).NegIso(end,1))+Plot.iv2TShift(i);

    % plug into Kind equation
    dPtimes = [Data.Time(ivIdx.dPmax2(i)) Data.Time(ivIdx.dPmin2(i)) ...
        Data.time_per];
    FitSinePres = data_kind (Fit.RCoef(i,:), FitSineTime, dPtimes);

    % find time point corresponding to Pmax
    [~, Idx] = min(abs(FitSinePres-Fit.RCoef(i,1)));

    PmaxT(i) = FitSineTime(Idx);

    plot(FitSineTime, FitSinePres, 'k--', PmaxT(i), Fit.RCoef(i,1), 'go');
    hold on;
end

% check the range of pressure values of Pmax. if the max p_max value is
% over 450, rescale y axis to (0, 300), so individual waveforms can be seen
ylim([0, Inf]);
if max(Fit.RCoef(:,1)) > 450
    ylim([0, 300]);
end

legend('Pressure', 'Isovolumic Points',  'Sinusoid Fit','Pmax', ...
    'Location','southoutside', 'Orientation', 'horizontal');

box on;
grid on;
hold off;

end
