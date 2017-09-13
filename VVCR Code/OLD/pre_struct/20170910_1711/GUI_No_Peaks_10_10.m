function varargout = GUI_No_Peaks_10_10(varargin)
% GUI_NO_PEAKS_10_10 MATLAB code for GUI_No_Peaks_10_10.fig
%      GUI_NO_PEAKS_10_10, by itself, creates a new GUI_NO_PEAKS_10_10 or raises the existing
%      singleton*.
%
%      H = GUI_NO_PEAKS_10_10 returns the handle to a new GUI_NO_PEAKS_10_10 or the handle to
%      the existing singleton*.
%
%      GUI_NO_PEAKS_10_10('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_NO_PEAKS_10_10.M with the given input arguments.
%
%      GUI_NO_PEAKS_10_10('Property','Value',...) creates a new GUI_NO_PEAKS_10_10 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_No_Peaks_10_10_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_No_Peaks_10_10_OpeningFcn via varargin.
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
time = handles.InVar.Data.Time;
Pres = handles.InVar.Data.Pres;
dPdt = handles.InVar.Data.dPdt;

dPminIdx = handles.InVar.Ext.dPminIdx;
dPminVal = handles.InVar.Ext.dPminVal;
% update number of Minima
set(handles.Min_num, 'String',num2str(length(dPminVal)));

dPmaxIdx = handles.InVar.Ext.dPmaxIdx;
dPmaxVal = handles.InVar.Ext.dPmaxVal;
% update number of maxima
set(handles.Max_num, 'String', num2str(length(dPmaxVal)));

% update axes of status
if length(dPmaxVal) == length(dPminVal)
    % display green check
    axes(handles.status_axes);
    imshow(handles.InVar.Green_Check); axis image; axis off
else
    % display red x
    axes(handles.status_axes);
    imshow(handles.InVar.Red_X); axis image; axis off
end

% set the total number of complete waveforms (editable text)
set(handles.NumWaves, 'String',num2str(length(dPmaxVal)));

% intialize these variables - used for undo button
handles.UNDOdPminIdx = [];
handles.UNDOdPminVal = [];
handles.dPmaxIdx = [];
handles.dPmaxVal = [];

% plot pressure, Dp/dt, and minima and maxima on appropriate axes
axes(handles.pressure_axes);
h = plot(time,Pres,'b',time(dPmaxIdx), Pres(dPmaxIdx), 'ro', time(dPminIdx), Pres(dPminIdx), 'ko');
set(h, 'HitTest', 'off');
set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',11);
xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);
legend('Pressure', 'dP/dt Max', 'dP/dt Min', 'Location','northoutside', 'Orientation', 'horizontal');
box on
grid on

axes(handles.dpdt_axes);
h2 = plot(time,dPdt, 'b', time(dPmaxIdx), dPmaxVal, 'ro', time(dPminIdx), dPminVal, 'ko');
set(h2, 'HitTest','off');
set(handles.dpdt_axes, 'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.dpdt_axes,'fontsize',11);
ylabel('dP/dt [mmHg/s]','FontSize',18);
legend('dP/dt','Maxima', 'Minima', 'Location', 'northoutside', 'Orientation', 'horizontal');
box on
grid on

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
handles.UNDOdPminIdx = handles.InVar.Ext.dPminIdx;
handles.UNDOdPminVal = handles.InVar.Ext.dPminVal;

handles.dPmaxIdx = handles.InVar.Ext.dPmaxIdx;
handles.dPmaxVal = handles.InVar.Ext.dPmaxVal;

% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
% disp('GUI_No_Peaks:');
% disp(['    Time: ',num2str(cp(1))]);
% disp(['    Pressure: ',num2str(cp(2))]);

% Extract variables from structure for a more clear workflow
time = handles.InVar.Data.Time;
Pres = handles.InVar.Data.Pres;
dPdt = handles.InVar.Data.dPdt;
dPminIdx = handles.InVar.Ext.dPminIdx;
dPminVal = handles.InVar.Ext.dPminVal;
dPmaxIdx = handles.InVar.Ext.dPmaxIdx;
dPmaxVal = handles.InVar.Ext.dPmaxVal;

% find indices of critical points within +- 0.5 seconds 
TPksInds = find(time(dPmaxIdx)>cp(1)-0.5 & time(dPmaxIdx)<cp(1)+0.5);
TMnsInds = find(time(dPminIdx)>cp(1)-0.5 & time(dPminIdx)<cp(1)+0.5);

% find distances of critical point within the neighborhood
TPksdst = abs(time(dPmaxIdx(TPksInds))-cp(1));
TMnsdst = abs(time(dPminIdx(TMnsInds))-cp(1));

if isempty(TMnsdst) 
    % remove the closest maximum (of dPdt)
    [~, Idx] = min(TPksdst); % get index of minimum number
    dPmaxValInd = TPksInds(Idx); % get index of vector (or scalar) returned by find()
    
    % remove peak
    dPmaxIdx(dPmaxValInd) = [];
    dPmaxVal(dPmaxValInd) = [];
    
    % update handles (global variable)
    handles.InVar.Ext.dPmaxIdx = dPmaxIdx;
    handles.InVar.Ext.dPmaxVal = dPmaxVal;
    
    % update number of maxima
    set(handles.Max_num, 'String', num2str(length(dPmaxVal)));
    
elseif min(TPksdst) < min(TMnsdst)
    % remove the closest maximum (of dPdt)
    [~, Idx] = min(TPksdst); % get index of minimum number
    dPmaxValInd = TPksInds(Idx); % get index of vector (or scalar) returned by find()
    
    % remove peak
    dPmaxIdx(dPmaxValInd) = [];
    dPmaxVal(dPmaxValInd) = [];
    
    % update handles (global variable)
    handles.InVar.Ext.dPmaxIdx = dPmaxIdx;
    handles.InVar.Ext.dPmaxVal = dPmaxVal;
    
    % update number of maxima
    set(handles.Max_num, 'String', num2str(length(dPmaxVal)));
    
elseif isempty(TPksdst) 
    % remove closest minimum (of dPdt)
    [~, Idx] = min(TMnsdst); % get index of minimum number
    MnsInd = TMnsInds(Idx); % get index of vector returned by find()
    
    % remove min
    dPminIdx(MnsInd) = [];
    dPminVal(MnsInd) = [];
    
    % update handles (global variable)
    handles.InVar.Ext.dPminIdx = dPminIdx;
    handles.InVar.Ext.dPminVal = dPminVal;
    
    % update number of minima
    set(handles.Min_num, 'String',num2str(length(dPminVal)));
    
elseif min(TPksdst) > min(TMnsdst)
    % remove closest minimum (of dPdt)
    [~, Idx] = min(TMnsdst); % get index of minimum number
    MnsInd = TMnsInds(Idx); % get index of vector returned by find()
    
    % remove min
    dPminIdx(MnsInd) = [];
    dPminVal(MnsInd) = [];
    
    % update handles (global variable)
    handles.InVar.Ext.dPminIdx = dPminIdx;
    handles.InVar.Ext.dPminVal = dPminVal;
    
    % update number of minima
    set(handles.Min_num, 'String',num2str(length(dPminVal)));
end

% update axes of status
if length(dPmaxVal) == length(dPminVal)
    % display green check
    axes(handles.status_axes);
    imshow(handles.InVar.Green_Check); axis image; axis off
else
    % display red x
    axes(handles.status_axes);
    imshow(handles.InVar.Red_X); axis image; axis off
end

% update global handles
guidata(hObject,handles);

% re-plot pressure graph
axes(handles.pressure_axes);
h = plot(time,Pres,'b',time(dPmaxIdx), Pres(dPmaxIdx), 'ro', time(dPminIdx), Pres(dPminIdx), 'ko');
set(h, 'HitTest', 'off');
set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',11);
xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);
legend('Pressure', 'dP/dt Max', 'dP/dt Min', 'Location','northoutside', 'Orientation', 'horizontal');
box on
grid on

% re-plot dP/dt graph
axes(handles.dpdt_axes);
h2 = plot(time,dPdt, 'b', time(dPmaxIdx), dPmaxVal, 'ro', time(dPminIdx), dPminVal, 'ko');
set(h2, 'HitTest','off');
set(handles.dpdt_axes, 'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.dpdt_axes,'fontsize',11);
ylabel('dP/dt [mmHg/s]','FontSize',18);
legend('dP/dt','Maxima', 'Minima', 'Location', 'northoutside', 'Orientation', 'horizontal');
box on
grid on

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

if ~isempty(handles.UNDOdPminIdx) 

    %Make the cusor a spinning wheel so user is aware program is busy
    set(handles.figure1, 'pointer', 'watch');
    drawnow;

    % retrieve the old critical points
    handles.InVar.Ext.dPminIdx = handles.UNDOdPminIdx;
    handles.InVar.Ext.dPminVal = handles.UNDOdPminVal;
    handles.InVar.Ext.dPmaxIdx = handles.dPmaxIdx;
    handles.InVar.Ext.dPmaxVal = handles.dPmaxVal;

    % Extract variables from structure for a more clear workflow
    time = handles.InVar.Time_D;
    Pres = handles.InVar.Pres_D;
    dPdt = handles.InVar.dPdt_D;
    dPminIdx = handles.InVar.Ext.dPminIdx;
    dPminVal = handles.InVar.Ext.dPminVal;
    dPmaxIdx = handles.InVar.Ext.dPmaxIdx;
    dPmaxVal = handles.InVar.Ext.dPmaxVal;

    % update number of minima
    set(handles.Min_num, 'String',num2str(length(dPminVal)));

    % update number of maxima
    set(handles.Max_num, 'String', num2str(length(dPmaxVal)));

    % update axes of status
    if length(dPmaxVal) == length(dPminVal)
        % display green check
        axes(handles.status_axes);
        imshow(handles.InVar.Green_Check); axis image; axis off
    else
        % display red x
        axes(handles.status_axes);
        imshow(handles.InVar.Red_X); axis image; axis off
    end

    % update global handles
    guidata(hObject,handles);

    % re-plot pressure graph
    axes(handles.pressure_axes);
    h = plot(time,Pres,'b',time(dPmaxIdx), Pres(dPmaxIdx), 'ro', time(dPminIdx), Pres(dPminIdx), 'ko');
    set(h, 'HitTest', 'off');
    set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
    set(handles.pressure_axes,'fontsize',11);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressue [mmHg]','FontSize',18);
    legend('Pressure', 'dP/dt Max', 'dP/dt Min', 'Location','northoutside', 'Orientation', 'horizontal');
    box on;
    grid on;

    % re-plot dP/dt graph
    axes(handles.dpdt_axes);
    h2 = plot(time,dPdt, 'b', time(dPmaxIdx), dPmaxVal, 'ro', time(dPminIdx), dPminVal, 'ko');
    set(h2, 'HitTest','off');
    set(handles.dpdt_axes, 'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
    set(handles.dpdt_axes,'fontsize',11);
    ylabel('dP/dt [mmHg/s]','FontSize',18);
    legend('dP/dt','Maxima', 'dPminVal', 'Location', 'northoutside', 'Orientation', 'horizontal');
    box on;
    grid on;

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

handles.InVar.TotNumWaves = numWaves;

% update global handles
guidata(hObject,handles);
    
% check the status: check if # Minima == # Maxima
if length(handles.InVar.Ext.dPminIdx) == length(handles.InVar.Ext.dPmaxIdx)
    
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
