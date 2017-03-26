function varargout = GUI_No_Peaks(varargin)
% GUI_NO_PEAKS MATLAB code for GUI_No_Peaks.fig
%      GUI_NO_PEAKS, by itself, creates a new GUI_NO_PEAKS or raises the existing
%      singleton*.
%
%      H = GUI_NO_PEAKS returns the handle to a new GUI_NO_PEAKS or the handle to
%      the existing singleton*.
%
%      GUI_NO_PEAKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_NO_PEAKS.M with the given input arguments.
%
%      GUI_NO_PEAKS('Property','Value',...) creates a new GUI_NO_PEAKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_No_Peaks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_No_Peaks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_No_Peaks

% Last Modified by GUIDE v2.5 25-Sep-2016 13:20:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_No_Peaks_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_No_Peaks_OutputFcn, ...
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

% --- Executes just before GUI_No_Peaks is made visible.
function GUI_No_Peaks_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_No_Peaks (see VARARGIN)

% Choose default command line output for GUI_No_Peaks
handles.output = hObject;

% set the input variable in the global handles environment
handles.InVar = cell2mat(varargin);

% Extract variables from structure for a more clear workflow
time = handles.InVar(1).Data;
Pres = handles.InVar(2).Data;
dPdt = handles.InVar(3).Data;

MinIdx = handles.InVar(1).Min;
Minima = handles.InVar(2).Min;
% update number of Minima
set(handles.Min_num, 'String',num2str(length(Minima)));

pksT = handles.InVar(1).Max;
pks = handles.InVar(2).Max;
% update number of maxima
set(handles.Max_num, 'String', num2str(length(pks)));

% update axes of status
if length(pks) == length(Minima)
    % display green check
    axes(handles.status_axes);
    imshow(handles.InVar(1).IM); axis image; axis off
else
    % display red x
    axes(handles.status_axes);
    imshow(handles.InVar(2).IM); axis image; axis off
end

% intialize these variables - used for undo button
handles.OldMinIdx = [];
handles.OldMinima = [];
handles.OldpksT = [];
handles.Oldpks = [];

% plot pressure, Dp/dt, and minima and maxima on appropriate axes
axes(handles.pressure_axes);
h = plot(time,Pres,'b',time(pksT), Pres(pksT), 'ro', time(MinIdx), Pres(MinIdx), 'ko');
set(h, 'HitTest', 'off');
set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',11);
xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);
legend('Pressure', 'dP/dt Max', 'dP/dt Min', 'Location','northoutside', 'Orientation', 'horizontal');
box on
grid on

axes(handles.dpdt_axes);
h2 = plot(time,dPdt, 'b', time(pksT), pks, 'ro', time(MinIdx), Minima, 'ko');
set(h2, 'HitTest','off');
set(handles.dpdt_axes, 'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.dpdt_axes,'fontsize',11);
ylabel('dP/dt [mmHg/s]','FontSize',18);
legend('dP/dt','Maxima', 'Minima', 'Location', 'northoutside', 'Orientation', 'horizontal');
box on
grid on

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_No_Peaks wait for user response (see UIRESUME)
uiwait(handles.figure1);
end

function GraphCallBack(hObject, eventdata, handles)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% store current mins and maxs as old min,maxs. This is done so the undo
% button can function
handles.OldMinIdx = handles.InVar(1).Min;
handles.OldMinima = handles.InVar(2).Min;

handles.OldpksT = handles.InVar(1).Max;
handles.Oldpks = handles.InVar(2).Max;

% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
% disp(['Time: ',num2str(cp(1))]);
% disp(['Pressure: ',num2str(cp(2))]);

% Extract variables from structure for a more clear workflow
time = handles.InVar(1).Data;
Pres = handles.InVar(2).Data;
dPdt = handles.InVar(3).Data;
MinIdx = handles.InVar(1).Min;
Minima = handles.InVar(2).Min;
pksT = handles.InVar(1).Max;
pks = handles.InVar(2).Max;

% find indices of critical points within +- 0.5 seconds 
TPksInds = find(time(pksT)>cp(1)-0.5 & time(pksT)<cp(1)+0.5);
TMnsInds = find(time(MinIdx)>cp(1)-0.5 & time(MinIdx)<cp(1)+0.5);

% find distances of critical point within the neighborhood
TPksdst = abs(time(pksT(TPksInds))-cp(1));
TMnsdst = abs(time(MinIdx(TMnsInds))-cp(1));

if isempty(TMnsdst) 
    % remove the closest maximum (of dPdt)
    [~, Idx] = min(TPksdst); % get index of minimum number
    pksInd = TPksInds(Idx); % get index of vector (or scalar) returned by find()
    
    % remove peak
    pksT(pksInd) = [];
    pks(pksInd) = [];
    
    % update handles (global variable)
    handles.InVar(1).Max = pksT;
    handles.InVar(2).Max = pks;
    
    % update number of maxima
    set(handles.Max_num, 'String', num2str(length(pks)));
    
elseif min(TPksdst) < min(TMnsdst)
    % remove the closest maximum (of dPdt)
    [~, Idx] = min(TPksdst); % get index of minimum number
    pksInd = TPksInds(Idx); % get index of vector (or scalar) returned by find()
    
    % remove peak
    pksT(pksInd) = [];
    pks(pksInd) = [];
    
    % update handles (global variable)
    handles.InVar(1).Max = pksT;
    handles.InVar(2).Max = pks;
    
    % update number of maxima
    set(handles.Max_num, 'String', num2str(length(pks)));
    
elseif isempty(TPksdst) 
    % remove closest minimum (of dPdt)
    [~, Idx] = min(TMnsdst); % get index of minimum number
    MnsInd = TMnsInds(Idx); % get index of vector returned by find()
    
    % remove min
    MinIdx(MnsInd) = [];
    Minima(MnsInd) = [];
    
    % update handles (global variable)
    handles.InVar(1).Min = MinIdx;
    handles.InVar(2).Min = Minima;
    
    % update number of minima
    set(handles.Min_num, 'String',num2str(length(Minima)));
    
elseif min(TPksdst) > min(TMnsdst)
    % remove closest minimum (of dPdt)
    [~, Idx] = min(TMnsdst); % get index of minimum number
    MnsInd = TMnsInds(Idx); % get index of vector returned by find()
    
    % remove min
    MinIdx(MnsInd) = [];
    Minima(MnsInd) = [];
    
    % update handles (global variable)
    handles.InVar(1).Min = MinIdx;
    handles.InVar(2).Min = Minima;
    
    % update number of minima
    set(handles.Min_num, 'String',num2str(length(Minima)));
end

% update axes of status
if length(pks) == length(Minima)
    % display green check
    axes(handles.status_axes);
    imshow(handles.InVar(1).IM); axis image; axis off
else
    % display red x
    axes(handles.status_axes);
    imshow(handles.InVar(2).IM); axis image; axis off
end

% update global handles
guidata(hObject,handles);

% re-plot pressure graph
axes(handles.pressure_axes);
h = plot(time,Pres,'b',time(pksT), Pres(pksT), 'ro', time(MinIdx), Pres(MinIdx), 'ko');
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
h2 = plot(time,dPdt, 'b', time(pksT), pks, 'ro', time(MinIdx), Minima, 'ko');
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
function varargout = GUI_No_Peaks_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.InVar;

% Destroy the GUI
delete(handles.figure1);
end

% --- Executes on button press in Undo.
function Undo_Callback(hObject, ~, handles)
% hObject    handle to Undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.OldMinIdx) 

    %Make the cusor a spinning wheel so user is aware program is busy
    set(handles.figure1, 'pointer', 'watch');
    drawnow;

    % retrieve the old critical points
    handles.InVar(1).Min = handles.OldMinIdx;
    handles.InVar(2).Min = handles.OldMinima;
    handles.InVar(1).Max = handles.OldpksT;
    handles.InVar(2).Max = handles.Oldpks;

    % Extract variables from structure for a more clear workflow
    time = handles.InVar(1).Data;
    Pres = handles.InVar(2).Data;
    dPdt = handles.InVar(3).Data;
    MinIdx = handles.InVar(1).Min;
    Minima = handles.InVar(2).Min;
    pksT = handles.InVar(1).Max;
    pks = handles.InVar(2).Max;

    % update number of minima
    set(handles.Min_num, 'String',num2str(length(Minima)));

    % update number of maxima
    set(handles.Max_num, 'String', num2str(length(pks)));

    % update axes of status
    if length(pks) == length(Minima)
        % display green check
        axes(handles.status_axes);
        imshow(handles.InVar(1).IM); axis image; axis off
    else
        % display red x
        axes(handles.status_axes);
        imshow(handles.InVar(2).IM); axis image; axis off
    end

    % update global handles
    guidata(hObject,handles);

    % re-plot pressure graph
    axes(handles.pressure_axes);
    h = plot(time,Pres,'b',time(pksT), Pres(pksT), 'ro', time(MinIdx), Pres(MinIdx), 'ko');
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
    h2 = plot(time,dPdt, 'b', time(pksT), pks, 'ro', time(MinIdx), Minima, 'ko');
    set(h2, 'HitTest','off');
    set(handles.dpdt_axes, 'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
    set(handles.dpdt_axes,'fontsize',11);
    ylabel('dP/dt [mmHg/s]','FontSize',18);
    legend('dP/dt','Maxima', 'Minima', 'Location', 'northoutside', 'Orientation', 'horizontal');
    box on;
    grid on;

    % Set cursor back to normal
    set(handles.figure1, 'pointer', 'arrow');
else
    msgbox('No point has been removed yet');
end

end

% --- Executes on button press in Next.
function Next_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check the status: check if # Minima == # Maxima
if length(handles.InVar(1).Min) == length(handles.InVar(1).Max)
    
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

% Hint: delete(hObject) closes the figure
delete(handles.figure1); % output function does not execute after this!

% maybe flip some flag to stop execution of further code in VVCR code?
end


% --- Executes on button press in Exit.
function Exit_Callback(hObject, ~, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set 1 output to false
handles.InVar(1).Max = false;

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

handles.InVar(3).Max = true;

% update handles globally
guidata(hObject, handles)

% call on uiresume so output function executes
uiresume(handles.figure1);
end
