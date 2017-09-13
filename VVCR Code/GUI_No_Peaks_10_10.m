function varargout = GUI_No_Peaks_10_10(varargin)
% GUI_NO_PEAKS_10_10 MATLAB code for GUI_No_Peaks_10_10.fig
%      GUI_NO_PEAKS_10_10, by itself, creates a new GUI_NO_PEAKS_10_10 or 
%      raises the existing singleton*.
%
%      H = GUI_NO_PEAKS_10_10 returns the handle to a new GUI_NO_PEAKS_10_10 
%      or the handle to the existing singleton*.
%
%      GUI_NO_PEAKS_10_10('CALLBACK',hObject,eventData,handles,...) calls the
%      local function named CALLBACK in GUI_NO_PEAKS_10_10.M with the given 
%      input arguments.
%
%      GUI_NO_PEAKS_10_10('Property','Value',...) creates a new GUI_NO_PEAKS-
%      _10_10 or raises the existing singleton*.  Starting from the left, 
%      property value pairs are applied to the GUI before GUI_No_Peaks_10_10-
%      _OpeningFcn gets called.  An unrecognized property name or invalid 
%      value makes property application stop.  All inputs are passed to GUI-
%      _No_Peaks_10_10_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_No_Peaks_10_10

% Last Modified by GUIDE v2.5 10-Oct-2016 18:41:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_No_Peaks_10_10_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_No_Peaks_10_10_OutputFcn, ...
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

% --- Executes just before GUI_No_Peaks_10_10 is made visible.
function GUI_No_Peaks_10_10_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_No_Peaks_10_10 (see VARARGIN)

% Choose default command line output for GUI_No_Peaks_10_10
handles.output = hObject;

% set the input variable in the global handles environment
handles.InVar = cell2mat(varargin);
handles.OutVar = false;

% Extract variables from structure for a more clear workflow
Data = handles.InVar.Data;
Extr = handles.InVar.Extr;

% intialize these variables - used for undo button
handles.UNDO.Extr = [];

% set the total number of complete waveforms (editable text)
set(handles.NumWaves, 'String', num2str(length(Extr.dPmaxVal)));

% update number of maxima/minima
set(handles.Min_num, 'String', num2str(length(Extr.dPminVal)));
set(handles.Max_num, 'String', num2str(length(Extr.dPmaxVal)));

[handles] = gui_nopeaks_plot (Data, Extr, handles); 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_No_Peaks_10_10 wait for user response (see UIRESUME)
uiwait(handles.figure1);
end

function GraphCallBack(hObject, eventdata, handles)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% store current mins and maxs as old min,maxs. This is done so the undo
% button can function
handles.UNDO.Extr = handles.InVar.Extr;

% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
% disp('GUI_No_Peaks:');
% disp(['    Time: ',num2str(cp(1))]);
% disp(['    Pressure: ',num2str(cp(2))]);

% Extract variables from structure for a more clear workflow
Data = handles.InVar.Data;
Extr = handles.InVar.Extr;

% find indices of critical points within +- 0.5 seconds 
TMxIdx = find(Data.Time(Extr.dPmaxIdx)>cp(1)-0.5 & ...
    Data.Time(Extr.dPmaxIdx)<cp(1)+0.5);
TMnIdx = find(Data.Time(Extr.dPminIdx)>cp(1)-0.5 & ...
    Data.Time(Extr.dPminIdx)<cp(1)+0.5);

% find distances of critical point within the neighborhood
TMxDist = abs(Data.Time(Extr.dPmaxIdx(TMxIdx))-cp(1));
TMnDist = abs(Data.Time(Extr.dPminIdx(TMnIdx))-cp(1));

if isempty(TMnDist) 
    % remove the closest maximum (of Data.dPdt)
    [~, Idx] = min(TMxDist); % get index of minimum number
    MxIdx = TMxIdx(Idx); % get index of vector returned by find()
    
    % remove peak
    Extr.dPmaxIdx(MxIdx) = [];
    Extr.dPmaxVal(MxIdx) = [];
    
elseif min(TMxDist) < min(TMnDist)
    % remove the closest maximum (of Data.dPdt)
    [~, Idx] = min(TMxDist); % get index of minimum number
    MxIdx = TMxIdx(Idx); % get index of vector returned by find()
    
    % remove peak
    Extr.dPmaxIdx(MxIdx) = [];
    Extr.dPmaxVal(MxIdx) = [];
    
elseif isempty(TMxDist) 
    % remove closest minimum (of Data.dPdt)
    [~, Idx] = min(TMnDist); % get index of minimum number
    MnIdx = TMnIdx(Idx); % get index of vector returned by find()
    
    % remove min
    Extr.dPminIdx(MnIdx) = [];
    Extr.dPminVal(MnIdx) = [];
    
elseif min(TMxDist) > min(TMnDist)
    % remove closest minimum (of Data.dPdt)
    [~, Idx] = min(TMnDist); % get index of minimum number
    MnIdx = TMnIdx(Idx); % get index of vector returned by find()
    
    % remove min
    Extr.dPminIdx(MnIdx) = [];
    Extr.dPminVal(MnIdx) = [];
    
end
    
% update handles (global variable)
handles.InVar.Extr = Extr;

% update number of maxima/minima
set(handles.Max_num, 'String', num2str(length(Extr.dPmaxVal)));
set(handles.Min_num, 'String', num2str(length(Extr.dPminVal)));

[handles] = gui_nopeaks_plot (Data, Extr, handles); 

% update global handles
guidata(hObject,handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');

end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_No_Peaks_10_10_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.InVar;

% Destroy the GUI
delete(hObject);

end

% --- Executes on button press in Undo.
function Undo_Callback(hObject, ~, handles)
% hObject    handle to Undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.UNDO.Extr) 

    %Make the cusor a spinning wheel so user is aware program is busy
    set(handles.figure1, 'pointer', 'watch');
    drawnow;

    % retrieve the old critical points
    handles.InVar.Extr = handles.UNDO.Extr;

    % Extract variables from structure for a more clear workflow
    Data = handles.InVar.Data;
    Extr = handles.InVar.Extr;

    % update number of maxima/minima
    set(handles.Max_num, 'String', num2str(length(Extr.dPmaxVal)));
    set(handles.Min_num, 'String', num2str(length(Extr.dPminVal)));

    [handles] = gui_nopeaks_plot (Data, Extr, handles); 

    % update global handles
    guidata(hObject,handles);

    % Set cursor back to normal
    set(handles.figure1, 'pointer', 'arrow');
else

    msgbox('No point has been removed yet');

end

end

% --- Executes on button press in Next.
function Next_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% grab the total number of complete waveforms from the editable text
numWaves = uint8(str2double(get(handles.NumWaves, 'String')));
Extr = handles.InVar.Extr; 

handles.InVar.TotNumWaves = numWaves;

% update global handles
guidata(hObject,handles);
    
% check the status: check if # Minima == # Maxima
if length(Extr.dPminIdx) == length(Extr.dPmaxIdx)
    
    % call on uiresume so output function executes
    uiresume(handles.figure1);
else
    
    warndlg(sprintf(' Program cannot proceed unless \n number of minima and maxima match!'));
end

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
    handles.InVar = false;
    guidata(hObject, handles);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

end

% --- Executes on button press in Exit.
function Exit_Callback(hObject, ~, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set 1 output to false
handles.InVar.Ext.dPmaxIdx = false;

% update handles globally
guidata(hObject, handles)

% call on uiresume so output function executes
uiresume(handles.figure1);
end

% --- Executes on button press in Discard.
function Discard_Callback(hObject, ~, handles)
% hObject    handle to Discard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.InVar.TotNumWaves = true;

% update handles globally
guidata(hObject, handles)

% call on uiresume so output function executes
uiresume(handles.figure1);
end

function NumWaves_Callback(hObject, ~, handles)
% hObject    handle to NumWaves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumWaves as text
%        str2double(get(hObject,'String')) returns contents of NumWaves as a double
end

% --- Executes during object creation, after setting all properties.
function NumWaves_CreateFcn(hObject, ~, ~)
% hObject    handle to NumWaves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function [handles] = gui_nopeaks_plot (Data, Extr, handles);

% update axes of status
if length(Extr.dPmaxVal) == length(Extr.dPminVal)
    % display green check
    axes(handles.status_axes);
    imshow(handles.InVar.Green_Check); axis image; axis off
else
    % display red x
    axes(handles.status_axes);
    imshow(handles.InVar.Red_X); axis image; axis off
end

% plot pressure, Dp/dt, and minima and maxima on appropriate axes
axes(handles.pressure_axes);
h = plot(Data.Time, Data.Pres,'b', ...
        Data.Time(Extr.dPmaxIdx), Data.Pres(Extr.dPmaxIdx), 'ro', ...
        Data.Time(Extr.dPminIdx), Data.Pres(Extr.dPminIdx), 'ko');

set(h, 'HitTest', 'off');
set(handles.pressure_axes,'ButtonDownFcn', ...
    @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',11);

xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);
legend('Pressure', 'dP/dt Max', 'dP/dt Min', 'Location', 'northoutside', ...
    'Orientation', 'horizontal');

box on
grid on

axes(handles.dpdt_axes);
h2 = plot(Data.Time, Data.OrigdPdt, 'g', ...
        Data.Time, Data.dPdt, 'b', ...
        Data.Time(Extr.dPmaxIdx), Extr.dPmaxVal, 'ro', ...
        Data.Time(Extr.dPminIdx), Extr.dPminVal, 'ko');

set(h2, 'HitTest','off');
set(handles.dpdt_axes, 'ButtonDownFcn', ...
    @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.dpdt_axes,'fontsize',11);

ylabel('dP/dt [mmHg/s]','FontSize',18);
legend('dP/dt', 'Maxima', 'Minima', 'Location', 'northoutside', ...
    'Orientation', 'horizontal');

box on
grid on

end
