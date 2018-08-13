function varargout = GUI_FitVanderpool (varargin)
% GUI_FitVanderpool MATLAB code for GUI_FitVanderpool.fig
%      GUI_FitVanderpool, by itself, creates a new GUI_FitVanderpool 
%      or raises the existing singleton*.
%
%      H = GUI_FitVanderpool returns the handle to a new GUI_FitVanderpool
%      or the handle to the existing singleton*.
%
%      GUI_FitVanderpool('CALLBACK',hObject,eventData,handles,...) calls
%      the local function named CALLBACK in GUI_FitVanderpool.M with the 
%      given input arguments.
%
%      GUI_FitVanderpool('Property','Value',...) creates a new GUI_Fit-
%      Takeuchi or raises the existing singleton*.  Starting from the left,
%      property value pairs are applied to the GUI before GUI_FitVanderpool-
%      _OpeningFcn gets called.  An unrecognized property name or invalid 
%      value makes property application stop.  All inputs are passed to 
%      GUI_FitVanderpool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only
%      one instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_FitVanderpool

% Last Modified by GUIDE v2.5 23-Sep-2017 09:40:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_FitVanderpool_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_FitVanderpool_OutputFcn, ...
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

% --- Executes just before GUI_FitVanderpool is made visible.
function GUI_FitVanderpool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_FitVanderpool (see VARARGIN)

% Choose default command line output for GUI_FitVanderpool
handles.output = hObject;

% set the input variable in the global handles environment
% passed as PeakStruct from VVCR_* script
handles.InVar = cell2mat(varargin(1));
handles.InVar.Data  = cell2mat(varargin(2));
handles.InVar.ivIdx = cell2mat(varargin(3));
handles.InVar.ivVal = cell2mat(varargin(4));
handles.InVar.ivSeg = cell2mat(varargin(5));

handles.Cycle = 1;
handles.CycMx = length(handles.InVar.ivIdx.Ps1);
set(handles.CycleMinus, 'Enable', 'off');

Rsq = handles.InVar.FitV.Rsq;
Cyc = handles.Cycle;
set(handles.CycleInd, 'String', ['Cycle #' num2str(Cyc, '%02i')]);
set(handles.RsqTxt,   'String', ['Rsq = ' num2str(Rsq(Cyc),'%6.4f')]);

% Extract Data, Indices/Values, and Fit Segments from passed structures.
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivVal = handles.InVar.ivVal;
ivSeg = handles.InVar.ivSeg;
FitV = handles.InVar.FitV;

% store first fit output into output structure.
handles.OutVar.FitV = handles.InVar.FitV;
handles.OutVar.Exit = 'Good';

% plot pressure, sinusoid fits
[handles] = takeuchi_plot_single (Data, ivSeg, FitV, Plot, handles);
[handles] = open_all_plots (Data, ivSeg, FitV, Plot, handles);

% Update handles.
guidata(hObject, handles);

% UIWAIT makes GUI_FitVanderpool wait for user response (see UIRESUME)
uiwait(handles.figure1);
end

% function that executes when user clicks on graph
function MainGraphCallback(hObject, eventdata, handles)

% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
disp('GUI_FitVanderpool>MainGraphCallback:');
disp(['    Time:     ',num2str(cp(1))]);
disp(['    Pressure: ',num2str(cp(2))]);

end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_FitVanderpool_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.OutVar;

% Destroy the GUI
delete(hObject);

end

% --- Executes on button press in CyclePlus.
function CyclePlus_Callback(hObject, eventdata, handles)
% hObject    handle to CyclePlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Rsq = handles.OutVar.FitV.Rsq;
Cyc = handles.Cycle;

Cyc = Cyc + 1;
if Cyc > 1
    set(handles.CycleMinus, 'Enable', 'on');
end
if Cyc == handles.CycMx
    set(handles.CyclePlus, 'Enable', 'off');
end
set(handles.CycleInd, 'String', ['Cycle #' num2str(Cyc,'%02i')]);
set(handles.RsqTxt,   'String', ['Rsq = ' num2str(Rsq(Cyc),'%6.4f')]);

handles.Cycle = Cyc;

% Extract Data, Indices/Values, and Fit Segments from passed structures.
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivSeg = handles.InVar.ivSeg;
FitV = handles.OutVar.FitV;

% plot pressure, sinusoid fits, update indicator
[handles] = takeuchi_plot_single (Data, ivSeg, FitV, Plot, handles);
if ishandle(handles.figure2)
    [handles] = takeuchi_plot_all (Data, ivSeg, FitV, Plot, handles);
end

% Update handles.
guidata(hObject, handles);

end


% --- Executes on button press in CycleMinus.
function CycleMinus_Callback(hObject, eventdata, handles)
% hObject    handle to CycleMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Rsq = handles.OutVar.FitV.Rsq;
Cyc = handles.Cycle;

Cyc = Cyc - 1;
if Cyc == 1
    set(handles.CycleMinus, 'Enable', 'off');
end
if Cyc < handles.CycMx
    set(handles.CyclePlus, 'Enable', 'on');
end
set(handles.CycleInd, 'String', ['Cycle #' num2str(Cyc,'%02i')]);
set(handles.RsqTxt,   'String', ['Rsq = ' num2str(Rsq(Cyc),'%6.4f')]);

handles.Cycle = Cyc;

% Extract Data, Indices/Values, and Fit Segments from passed structures.
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivSeg = handles.InVar.ivSeg;
FitV = handles.OutVar.FitV;

% plot pressure, sinusoid fits, update indicator
[handles] = takeuchi_plot_single (Data, ivSeg, FitV, Plot, handles);
if ishandle(handles.figure2)
    [handles] = takeuchi_plot_all (Data, ivSeg, FitV, Plot, handles);
end

% Update handles.
guidata(hObject, handles);

end

% --- Executes on button press in Remove.
function Remove_Callback(hObject, eventdata, handles)
% hObject    handle to Remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

WaveRm = handles.Cycle;
disp(['GUI_FitVanderpool>Remove: wave ' num2str(WaveRm, '%02i') ...
    ' marked removed']);

handles.OutVar.FitV.BadCyc(WaveRm) = 1;

% Plot the results
ivSeg = handles.InVar.ivSeg;
FitV = handles.OutVar.FitV;
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
[handles] = takeuchi_plot_single (Data, ivSeg, FitV, Plot, handles);
if ishandle(handles.figure2)
    [handles] = takeuchi_plot_all (Data, ivSeg, FitV, Plot, handles);
end

set(handles.Include, 'Enable', 'on');
set(handles.Remove,  'Enable', 'off');

% update global handles & set cursor back to normal
guidata(hObject,handles);
set(handles.figure1, 'pointer', 'arrow');

end

% --- Executes on button press in Done.
function Done_Callback(~, ~, handles)
% hObject    handle to Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ishandle(handles.figure2)
    close(handles.figure2);
end

% call on uiresume so output function executes
uiresume(handles.figure1);
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ishandle(handles.figure2)
    close(handles.figure2);
end

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

% --- Executes on button press in Exit.
function Exit_Callback(hObject, ~, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% keep in mind when the exit button is pressed, the current
% patient, i, will not be evaluated
                
if ishandle(handles.figure2)
    close(handles.figure2);
end

% set output to false
handles.OutVar.Exit = false;

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

if ishandle(handles.figure2)
    close(handles.figure2);
end

% set outputs to true, indicating Discard button
handles.OutVar.Exit = true;

% update handles globally
guidata(hObject, handles)

% call on uiresume so output function executes
uiresume(handles.figure1);
end

% --- Executes on button press in Include.
function Include_Callback(hObject, ~, handles)
% hObject    handle to Include (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

disp(['GUI_FitVanderpool>Include: wave ' num2str(handles.Cycle, ...
   '%02i') ' included in final analysis']);

handles.OutVar.FitV.BadCyc(handles.Cycle) = 0;

% Extract Data, Values, Fit Segments, Plots, & Segments from handles.
FitV = handles.OutVar.FitV;
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivSeg = handles.InVar.ivSeg;

[handles] = takeuchi_plot_single (Data, ivSeg, FitV, Plot, handles);
if ishandle(handles.figure2)
    [handles] = takeuchi_plot_all (Data, ivSeg, FitV, Plot, handles);
end

set(handles.Include, 'Enable', 'off');
set(handles.Remove,  'Enable', 'on');

% update global handles
guidata(hObject,handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');

end

% --- Function that updates the main plot
function [handles] = takeuchi_plot_single (Data, ivSeg, Fit, Plot, handles);

cycid = handles.Cycle;

if Fit.BadCyc(cycid)
    set(handles.Include, 'Enable', 'on');
    set(handles.Remove,  'Enable', 'off');
else
    set(handles.Include, 'Enable', 'off');
    set(handles.Remove,  'Enable', 'on');
end

axes(handles.pressure_axes);

h = plot(Data.Time_D,Data.Pres_D,'b', ...
         Plot.iv1PlotTime,Plot.iv1PlotPres,'ro');
set(h, 'HitTest', 'off');

set(handles.pressure_axes,'ButtonDownFcn', ...
    @(hObject, eventdata)MainGraphCallback(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',12);

title('Takeuchi Sinusoidal Fitting','FontSize',16);
xlabel('Time [s]','FontSize',14);
ylabel('Pressure [mmHg]','FontSize',14);

hold on;

mystp = Data.time_step/2;

% Attain the sinusoid fit for this cycle (so Pmax can be visualized)
% obtain the range of time of each peak, then normalize to zero
FitSineTime = Data.Time_D(ivSeg.iv1Time(cycid).PosIso(1,1)):mystp: ...
    Data.Time_D(ivSeg.iv1Time(cycid).NegIso(end,1));

% plug into Naeiji equation that was just solved for; normalize range
% to start at one (as was done in fitting).
FitSinePres = Fit.RCoef(cycid,1) + Fit.RCoef(cycid,2)* ...
        sin(Fit.RCoef(cycid,3)*(FitSineTime-FitSineTime(1)) + ...
        Fit.RCoef(cycid,4));
      
% find time point corresponding to Pmax
[~, Idx] = min(abs(FitSinePres-Fit.PIsoMax(cycid)));

PmaxT = FitSineTime(Idx);

if Fit.BadCyc(cycid)
    plot(FitSineTime, FitSinePres, 'r--', PmaxT, Fit.PIsoMax(cycid), 'rx');
else
    plot(FitSineTime, FitSinePres, 'k--', PmaxT, Fit.PIsoMax(cycid), 'go');
end

% Set reasonable plot limits.
xmn = FitSineTime(1)-0.1;
xmx = FitSineTime(end)+0.1;
ymx = max(Fit.PIsoMax)+5;

xlim([xmn xmx]);
if ymx > 300
    ylim([0, 300]);
else
    ylim([0, ymx]);
end

legend('Pressure', 'Isovolumic Points', 'Sinusoid Fit', 'Pmax', ...
    'Location', 'southoutside', 'Orientation', 'horizontal');

box on;
grid on;
hold off;

end

% --- Executes on button press in AllPlotsGraph.
function AllPlotsGraph_Callback(hObject, ~, handles)
% hObject    handle to AllPlotsGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Extract Data, Indices/Values, and Fit Segments from passed structures.
Data = handles.InVar.Data;
Plot = handles.InVar.Plot;
ivSeg = handles.InVar.ivSeg;
FitV = handles.OutVar.FitV;

% Create all pressures figure
handles = open_all_plots (Data, ivSeg, FitV, Plot, handles);

% Update handles.
guidata(hObject, handles);

end

% --- Executes when the plot is pressed in figure2
function SubGraphCallback(hObject, eventdata, handles)

cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];

Rsq = handles.OutVar.FitV.Rsq;
Data = handles.InVar.Data;
ivIdx = handles.InVar.ivIdx;

WaveNumPosRm = find(Data.Time_D(ivIdx.Ps1_D)<cp(1));
WaveNumNegRm = find(Data.Time_D(ivIdx.Ne1_D)>cp(1));

if ~isempty(WaveNumPosRm) && ~isempty(WaveNumNegRm)

    % find the common number. the last EDP that is smaller then
    WavePk = find(WaveNumPosRm==WaveNumNegRm(1));

    if ~isempty(WavePk)

        handles.Cycle = WavePk;

        if WavePk == 1
            set(handles.CyclePlus,  'Enable', 'on');
            set(handles.CycleMinus, 'Enable', 'off');
        elseif WavePk == handles.CycMx
            set(handles.CyclePlus,  'Enable', 'off');
            set(handles.CycleMinus, 'Enable', 'on');
        else
            set(handles.CyclePlus,  'Enable', 'on');
            set(handles.CycleMinus, 'Enable', 'on');
        end

        Cyc = handles.Cycle;
        set(handles.CycleInd, 'String', ['Cycle #' num2str(Cyc, '%02i')]);
        set(handles.RsqTxt,   'String', ['Rsq = ' num2str(Rsq(Cyc),'%6.4f')]);

        Data = handles.InVar.Data;
        Plot = handles.InVar.Plot;
        ivSeg = handles.InVar.ivSeg;
        FitV = handles.OutVar.FitV;

        [handles] = takeuchi_plot_single (Data, ivSeg, FitV, Plot, handles);
        [handles] = takeuchi_plot_all (Data, ivSeg, FitV, Plot, handles);

        % update figure1 handles
        guidata(handles.figure1, handles);

    end
end

end

% --- Creates or brings up figure2
function [handles] = open_all_plots (Data, ivSeg, Fit, Plot, handles);

f2h = findall(0, 'tag', 'FAllPlots');

if isempty(f2h)

    handles.figure2 = figure ('Name', 'Takeuchi All Pressure Waveforms',...
        'Units', 'characters',...
        'Position', [30 35 140 30],...
        'NumberTitle', 'off', ...
        'MenuBar', 'none',...
        'Color', [0.94 0.94 0.94],...
        'Tag', 'FAllPlots');

    h = axes ('Position', [0.1 0.12 0.85 0.80]);

    handles = takeuchi_plot_all (Data, ivSeg, Fit, Plot, handles);

else

    figure(f2h);

end

end

% --- Pushes data to axes in figure2
function [handles] = takeuchi_plot_all (Data, ivSeg, Fit, Plot, handles);

axes(handles.figure2.CurrentAxes);

h = plot(Data.Time_D,Data.Pres_D,'b', ...
         Plot.iv1PlotTime,Plot.iv1PlotPres,'ro');
set(h, 'HitTest', 'off');

set(handles.figure2.CurrentAxes,'ButtonDownFcn', ...
    @(hObject, eventdata)SubGraphCallback(hObject, eventdata, handles));
set(handles.figure2.CurrentAxes,'fontsize',12);

xlabel('Time [s]','FontSize',12);
ylabel('Pressure [mmHg]','FontSize',12);

hold on;

mystp = Data.time_step/2;
mysz = length(ivSeg.iv1Time);
PmaxT = zeros(mysz,1);

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:mysz

    % obtain the range of time of each peak, then normalize to zero
    FitSineTime = Data.Time_D(ivSeg.iv1Time(i).PosIso(1,1)):mystp: ...
        Data.Time_D(ivSeg.iv1Time(i).NegIso(end,1));

    % plug into Naeiji equation that was just solved for; normalize range
    % to start at one (as was done in fitting).
    FitSinePres = Fit.RCoef(i,1) + Fit.RCoef(i,2)*sin(Fit.RCoef(i,3)* ...
      (FitSineTime-FitSineTime(1)) + Fit.RCoef(i,4));

    % find time point corresponding to Pmax
    [~, Idx] = min(abs(FitSinePres-Fit.PIsoMax(i)));

    PmaxT(i) = FitSineTime(Idx);

    if Fit.BadCyc(i)
        plot(FitSineTime, FitSinePres, 'r--', PmaxT(i), Fit.PIsoMax(i), 'rx');
    else
        plot(FitSineTime, FitSinePres, 'k--', PmaxT(i), Fit.PIsoMax(i), 'go');
    end
end

ymx = max(Fit.PIsoMax)+5;

% Bound the current cycle.
cycid = handles.Cycle;
xmn = Data.Time_D(ivSeg.iv1Time(cycid).PosIso(1,1))-0.05;
xmx = Data.Time_D(ivSeg.iv1Time(cycid).NegIso(end,1))+0.05;
plot([xmn xmn], [0, ymx], 'r--');
plot([xmx xmx], [0, ymx], 'r--');

% Set reasonable plot limits.
xlim([0 Data.Time_D(end)]);
if ymx > 300
    ylim([0, 300]);
else
    ylim([0, ymx]);
end

box on;
grid on;
hold off;

end
