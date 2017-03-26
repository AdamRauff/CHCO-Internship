function varargout = GUI_SINU_FIT_10_03(varargin)
% GUI_SINU_FIT_10_03 MATLAB code for GUI_SINU_FIT_10_03.fig
%      GUI_SINU_FIT_10_03, by itself, creates a new GUI_SINU_FIT_10_03 or raises the existing
%      singleton*.
%
%      H = GUI_SINU_FIT_10_03 returns the handle to a new GUI_SINU_FIT_10_03 or the handle to
%      the existing singleton*.
%
%      GUI_SINU_FIT_10_03('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SINU_FIT_10_03.M with the given input arguments.
%
%      GUI_SINU_FIT_10_03('Property','Value',...) creates a new GUI_SINU_FIT_10_03 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_SINU_FIT_10_03_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_SINU_FIT_10_03_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_SINU_FIT_10_03

% Last Modified by GUIDE v2.5 03-Oct-2016 00:36:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_SINU_FIT_10_03_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_SINU_FIT_10_03_OutputFcn, ...
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

% --- Executes just before GUI_SINU_FIT_10_03 is made visible.
function GUI_SINU_FIT_10_03_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_SINU_FIT_10_03 (see VARARGIN)

% Choose default command line output for GUI_SINU_FIT_10_03
handles.output = hObject;

% set the input variable in the global handles environment
handles.InVar = cell2mat(varargin);

% Extract variables from structure for a more clear workflow
% time and pressure vectors
time = handles.InVar(1).Data; % recall this is the 2x sampled time vector
Pres = handles.InVar(2).Data; % recall this is the 2x sampled pressure vector

% extract isovolmic points (times). These are structures
isovoltime = handles.InVar(1).ivt;

% Extract wave fits
waveFit = handles.InVar(2).iv;

% iso volumic points in array (for plotting)
totIsoTimePoints = handles.InVar(1).isoPts;
totIsoPresPoints = handles.InVar(2).isoPts;

% intialize these variables for the undo button
handles.OldIsoT = [];
handles.OldIsoP = [];

% EDP - end diastolic pressure
EDP = handles.InVar(1).Misc;
P_max2 = handles.InVar(2).Misc;

% regression constants
c_tot2 = handles.InVar(2).Cs;

% initial conditions
IC = handles.InVar(1).Cs;

% intialize old variable - used for the undo button
handles.OldIsoT = [];

% set editable text boxes with ICs
set(handles.Mean_txt, 'String',num2str(IC(1)));
set(handles.Amp_txt, 'String',num2str(IC(2)));
set(handles.Freq_txt, 'String',num2str(IC(3)));
set(handles.Phase_txt, 'String',num2str(IC(4)));

% plot pressure, sinusoid fits
axes(handles.pressure_axes);
h = plot(time,Pres,'b',totIsoTimePoints,totIsoPresPoints,'ro');
set(h, 'HitTest', 'off');
set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',12);
title('Sinusoidal Fitting','FontSize',20);
xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);
hold on;

PmaxT = zeros(length(EDP),1);

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:length(EDP)

    % obtain the range of time of each peak
    interval = time(isovoltime(i).PosIso(1,1)):0.002:time(isovoltime(i).NegIso(end,1));

    % plug into Naeiji equation that was just solved for 
    FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + c_tot2(i,4));
    
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
legend('Pressure', 'Isovolumic Points', 'Sinusoid Fit','Pmax', 'Location','southoutside', 'Orientation', 'horizontal');
box on;
grid on;
hold off;

% store the wavefit, so output tracks which wave forms did not have a good
% fit
handles.OutVar(1).output = waveFit;

% Store the pmax values
handles.OutVar(2).output = P_max2;

% store regression constants
handles.OutVar(3).output = c_tot2;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_SINU_FIT_10_03 wait for user response (see UIRESUME)
uiwait(handles.figure1);
end

% function that executes when user clicks on graph
function GraphCallBack(hObject, eventdata, handles)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% store the current isovolumic points as OldIsoT and OldIsoP
% intialize these variables for the undo button
handles.OldIsoT = handles.InVar(1).isoPts;
handles.OldIsoP = handles.InVar(2).isoPts;
handles.OldIsoVolT = handles.InVar(1).ivt;
handles.OldEDPT = handles.InVar(1).EDPs;
handles.OldEDPNT = handles.InVar(2).EDPs;
handles.OldIsoVol = handles.InVar(1).iv;
handles.OldEDP = handles.InVar(1).Misc;
handles.OldPksT = handles.InVar(1).Crit;
handles.OldMinIdx = handles.InVar(2).Crit;

% get the current point
cp(1,:) = [eventdata.IntersectionPoint(1), eventdata.IntersectionPoint(2)];
disp(['Time: ',num2str(cp(1))]);
disp(['Pressure: ',num2str(cp(2))]);

% obtain variables from InVar Struct for a clear workflow
% time and pressure vectors
time = handles.InVar(1).Data;
Pres = handles.InVar(2).Data; % recall this is the 2x sampled pressure vector
Oldtime = handles.InVar(2).ivt; % old time vector. 1/2 the points

EDP = handles.InVar(1).Misc;

% EDP - end diastolic pressure
EDP_T = handles.InVar(1).EDPs;
EDP_NT = handles.InVar(2).EDPs;

% pass the time indexes of minima and maxima
pksT = handles.InVar(1).Crit;
MinIdx = handles.InVar(2).Crit;

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
        
        % Each row holds the data for a single pressure wave
        isovol = struct('PosIso',cell(length(EDP),1),'NegIso',cell(length(EDP),1));
        isovoltime = struct('PosIso',cell(length(EDP),1),'NegIso',cell(length(EDP),1));

        for i = 1: length(EDP)
            % Positive
            % convert index to index of vector containg 2x data points 
            P2 = find(round(time,3)==round(Oldtime(pksT(i)),3));
            isovoltime(i).PosIso(:,1) = (EDP_T(i):1:P2)'; % keep in mind these are indices of the time vector, not real time points
            for j = 1:length(isovoltime(i).PosIso)
                isovol(i).PosIso(j,1) = Pres(isovoltime(i).PosIso(j,1)); % reall pressure points [mmHg]
            end
    
            %Negative
            % convert index to index of vector containg 2x data points
            P1 = find(round(time,3)==round(Oldtime(MinIdx(i)),3));
            isovoltime(i).NegIso(:,1) = (P1:1:EDP_NT(i))';
            for j = 1:length(isovoltime(i).NegIso)
                isovol(i).NegIso(j,1) = Pres(isovoltime(i).NegIso(j,1));
            end
        end
        
        % obtain current ICs
        Mea = str2double(get(handles.Mean_txt,'String'));
        Amp = str2double(get(handles.Amp_txt,'String'));
        Fre = str2double(get(handles.Freq_txt,'String'));
        Pha = str2double(get(handles.Phase_txt,'String'));
        
        % pre - allocate 
        waveFit = zeros(length(EDP),1);
        PmaxT = zeros(length(EDP),1);
        P_max2 = zeros(length(EDP),1);
        r_square2 = zeros(length(EDP),1);
        c_tot2 = zeros(length(EDP),4);
        
        totIsoTimePoints = [];
        totIsoPresPoints = [];
        
        for i = 1:length(EDP)
            WaveTs = [time(isovoltime(i).PosIso)'; time(isovoltime(i).NegIso)'];
            WavePs = [isovol(i).PosIso; isovol(i).NegIso];

            % from Naeiji et al, single beat method of VVC
            sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 

            % ICs
            c2=[Mea, Amp, Fre, Pha];

            [c,resnorm,~]=lsqnonlin(sin_fun2,c2); %least squares fitting

            Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));

            % r^2 value
            r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;

            % if the fit of the wave was bad, mark that wave
            if r_square2(i) <0.90
               waveFit(i) = 1; 
            end

            c_tot2(i,:)=c; %getting all the c values in a matrix

            P_max2(i)=c(1)+abs(c(2)); %first equation pmax, A+B
            %or the other definition of pmax

            % -------------------------------------------------------------
            % NOTE the absolute value of the amplitude is taken!!!!!!!
            % refer to patient HA002019, Wave 11 for an example of why

            % sometime amplitude of given equation solves for negative ( with a
            % significant phase shift, which makes a good fit (r^2 > 0.99).
            % ---------------------------------------------------------------

            % store the time points and pressure points in one array for easy
            % plotting 
            totIsoTimePoints = [totIsoTimePoints; WaveTs];
            totIsoPresPoints = [totIsoPresPoints; WavePs];
        end
        
        % update global handles
        handles.InVar(1).ivt = isovoltime;
        handles.InVar(1).EDPs = EDP_T;
        handles.InVar(2).EDPs = EDP_NT;
        handles.InVar(1).iv = isovol;
        handles.InVar(1).isoPts = totIsoTimePoints;
        handles.InVar(2).isoPts = totIsoPresPoints;
        handles.InVar(1).Misc = EDP;
        handles.InVar(1).Crit = pksT;
        handles.InVar(2).Crit = MinIdx;

        % print to command line the waves that were not fit correctly. This is used
        % as a debugger to check that the "bad" waves, the ones that don't have a
        % good fit, are not utilized in the VVCR calculation.
        indX = find(waveFit==1); % find indices of the bad waves
        if ~isempty(indX)
            disp('The following waves did NOT have a good fit (will not be included)');
            disp(['Wave(s): ', num2str(indX')]);
        else
            disp('All waves seemed to fit well!');
        end
        
        % plot pressure, sinusoid fits
        axes(handles.pressure_axes);
        h = plot(time,Pres,'b',totIsoTimePoints,totIsoPresPoints,'ro');
        set(h, 'HitTest', 'off');
        set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
        set(handles.pressure_axes,'fontsize',12);
        title('Sinusoidal Fitting','FontSize',20);
        xlabel('Time [s]','FontSize',18);
        ylabel('Pressue [mmHg]','FontSize',18);
        hold on;

        % Attain the sinusoid fit for all points (so Pmax can be visualized
        for i = 1:length(EDP)

            % obtain the range of time of each peak
            interval = time(isovoltime(i).PosIso(1,1)):0.002:time(isovoltime(i).NegIso(end,1));

            % plug into Naeiji equation that was just solved for 
            FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + c_tot2(i,4));

            % find time point corresponding to Pmax
            [~, Idx] = min(abs(FitSinePres-P_max2(i)));

            PmaxT(i) = interval(Idx);

            plot(interval, FitSinePres, 'k--', PmaxT(i), P_max2(i), 'go');
            hold on;
        end

        % check the range of pressure values of Pmax. if the max p_max value is
        % over 450, rescale y axis to (0, 300), so individual waveforms can be seen
        maxP = max(P_max2);
        if maxP > 450
            ylim([0, 300]);
        else
            ylim([0, abs(maxP)+5]);
        end
        legend('Pressure', 'Isovolumic Points', 'Sinusoid Fit', 'Pmax', 'Location','southoutside', 'Orientation', 'horizontal');
        box on;
        grid on;
        hold off;

        % store the wavefit, so output tracks which wave forms did not have a good
        % fit
        handles.OutVar(1).output = waveFit;

        % Store the pmax values
        handles.OutVar(2).output = P_max2;

        % store regression constants
        handles.OutVar(3).output = c_tot2;
    end
end

% update global handles
guidata(hObject,handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_SINU_FIT_10_03_OutputFcn(~, ~, handles) 
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
function calculate_Callback(hObject, ~, handles)
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
time = handles.InVar(1).Data;
Pres = handles.InVar(2).Data;

% extract isovolmic points (times). These are structures
isovoltime = handles.InVar(1).ivt;

% extract isovolmic points (pressures). These are structures
isovol = handles.InVar(1).iv;

% iso volumic points in array (for plotting)
totIsoTimePoints = handles.InVar(1).isoPts;
totIsoPresPoints = handles.InVar(2).isoPts;

% EDP - end diastolic pressure
EDP = handles.InVar(1).Misc;

% pre - allocate 
waveFit = zeros(length(EDP),1);
PmaxT = zeros(length(EDP),1);
P_max2 = zeros(length(EDP),1);
r_square2 = zeros(length(EDP),1);
c_tot2 = zeros(length(EDP),4);

%%% calculate regression per pressure wave
for i = 1:length(EDP)
    
    WaveTs = [time(isovoltime(i).PosIso)'; time(isovoltime(i).NegIso)'];
    WavePs = [isovol(i).PosIso; isovol(i).NegIso];

    % from Naeiji et al, single beat method of VVC
    sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 

    % ICs
    c2=[Mea, Amp, Fre, Pha];

    [c,resnorm,~]=lsqnonlin(sin_fun2,c2); %least squares fitting

    Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));

    % r^2 value
    r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;

    % if the fit of the wave was bad, mark that wave
    if r_square2(i) <0.90
       waveFit(i) = 1; 
    end

    c_tot2(i,:)=c; %getting all the c values in a matrix

    P_max2(i)=c(1)+abs(c(2)); %first equation pmax, A+B
    %or the other definition of pmax

    % -------------------------------------------------------------
    % NOTE the absolute value of the amplitude is taken!!!!!!!
    % refer to patient HA002019, Wave 11 for an example of why

    % sometime amplitude of given equation solves for negative ( with a
    % significant phase shift, which makes a good fit (r^2 > 0.99).
    % ---------------------------------------------------------------

%     % store the time points and pressure points in one array for easy
%     % plotting 
%     totIsoTimePoints = [totIsoTimePoints; WaveTs];
%     totIsoPresPoints = [totIsoPresPoints; WavePs];
end

% print to command line the waves that were not fit correctly. This is used
% as a debugger to check that the "bad" waves, the ones that don't have a
% good fit, are not utilized in the VVCR calculation.
indX = find(waveFit==1); % find indices of the bad waves
if ~isempty(indX)
    disp('The following waves did NOT have a good fit (will not be included)');
    disp(['Wave(s): ', num2str(indX')]);
else
    disp('All waves seemed to fit well!');
end

% plot pressure, sinusoid fits
axes(handles.pressure_axes);
h = plot(time,Pres,'b',totIsoTimePoints,totIsoPresPoints,'ro');
set(h, 'HitTest', 'off');
set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',12);
title('Sinusoidal Fitting','FontSize',20);
xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);
hold on;

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:length(EDP)
    % obtain the range of time of each peak
    interval = time(isovoltime(i).PosIso(1,1)):0.002:time(isovoltime(i).NegIso(end,1));

    % plug into Naeiji equation that was just solved for 
    FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + c_tot2(i,4));

    % find time point corresponding to Pmax
    [~, Idx] = min(abs(FitSinePres-P_max2(i)));

    PmaxT(i) = interval(Idx);

    plot(interval, FitSinePres, 'k--', PmaxT(i), P_max2(i), 'go');
    hold on;
end

% check the range of pressure values of Pmax. if the max p_max value is
% over 450, rescale y axis to (0, 300), so individual waveforms can be seen
maxP = max(P_max2);
if maxP > 450
    ylim([0, 300]);
else
    ylim([0, abs(maxP)+5]);
end
legend('Pressure', 'Isovolumic Points', 'Sinusoid Fit', 'Pmax', 'Location','southoutside', 'Orientation', 'horizontal');
box on;
grid on;
hold off;

% store the wavefit, so output tracks which wave forms did not have a good
% fit
handles.OutVar(1).output = waveFit;

% Store the pmax values
handles.OutVar(2).output = P_max2;

% store regression constants
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
function Undo_Callback(hObject, ~, handles)
% hObject    handle to Undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;

% if the handles.Old variables have been created (user has clicked on the
% plot and removed pressure waveform(s).
if ~isempty(handles.OldIsoT)
    % restore the old isovolumic points 
    handles.InVar(1).isoPts = handles.OldIsoT;
    handles.InVar(2).isoPts = handles.OldIsoP;
    handles.InVar(1).ivt = handles.OldIsoVolT;
    handles.InVar(1).EDPs = handles.OldEDPT;
    handles.InVar(2).EDPs = handles.OldEDPNT;
    handles.InVar(1).iv = handles.OldIsoVol;
    handles.InVar(1).Misc = handles.OldEDP;
    handles.InVar(1).Crit = handles.OldPksT;
    handles.InVar(2).Crit = handles.OldMinIdx;
    
    % calculate sinusoids based on current ICs!
    Mea = str2double(get(handles.Mean_txt,'String'));
    Amp = str2double(get(handles.Amp_txt,'String'));
    Fre = str2double(get(handles.Freq_txt,'String'));
    Pha = str2double(get(handles.Phase_txt,'String'));

    % obtain variables from InVar Struct for a clear workflow
    time = handles.InVar(1).Data;
    Pres = handles.InVar(2).Data;
    isovoltime = handles.InVar(1).ivt;
    isovol = handles.InVar(1).iv;
    % iso volumic points in array (for plotting)
    totIsoTimePoints = handles.InVar(1).isoPts;
    totIsoPresPoints = handles.InVar(2).isoPts;
    % EDP - end diastolic pressure
    EDP = handles.InVar(1).Misc;

    % pre - allocate 
    waveFit = zeros(length(EDP),1);
    PmaxT = zeros(length(EDP),1);
    P_max2 = zeros(length(EDP),1);
    r_square2 = zeros(length(EDP),1);
    c_tot2 = zeros(length(EDP),4);

    % calculate regression per pressure wave
    for i = 1:length(EDP)

        WaveTs = [time(isovoltime(i).PosIso)'; time(isovoltime(i).NegIso)'];
        WavePs = [isovol(i).PosIso; isovol(i).NegIso];

        % from Naeiji et al, single beat method of VVC
        sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 

        % ICs
        c2=[Mea, Amp, Fre, Pha];

        [c,resnorm,~]=lsqnonlin(sin_fun2,c2); %least squares fitting

        Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));

        % r^2 value
        r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;

        % if the fit of the wave was bad, mark that wave
        if r_square2(i) <0.90
           waveFit(i) = 1; 
        end

        c_tot2(i,:)=c; %getting all the c values in a matrix

        P_max2(i)=c(1)+abs(c(2)); %first equation pmax, A+B
        %or the other definition of pmax

        % -------------------------------------------------------------
        % NOTE the absolute value of the amplitude is taken!!!!!!!
        % refer to patient HA002019, Wave 11 for an example of why

        % sometime amplitude of given equation solves for negative ( with a
        % significant phase shift, which makes a good fit (r^2 > 0.99).
        % ---------------------------------------------------------------

    %     % store the time points and pressure points in one array for easy
    %     % plotting 
    %     totIsoTimePoints = [totIsoTimePoints; WaveTs];
    %     totIsoPresPoints = [totIsoPresPoints; WavePs];
    end

    % print to command line the waves that were not fit correctly. This is used
    % as a debugger to check that the "bad" waves, the ones that don't have a
    % good fit, are not utilized in the VVCR calculation.
    indX = find(waveFit==1); % find indices of the bad waves
    if ~isempty(indX)
        disp('The following waves did NOT have a good fit (will not be included)');
        disp(['Wave(s): ', num2str(indX')]);
    else
        disp('All waves seemed to fit well!');
    end

    % plot pressure, sinusoid fits
    axes(handles.pressure_axes);
    h = plot(time,Pres,'b',totIsoTimePoints,totIsoPresPoints,'ro');
    set(h, 'HitTest', 'off');
    set(handles.pressure_axes,'ButtonDownFcn', @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
    set(handles.pressure_axes,'fontsize',12);
    title('Sinusoidal Fitting','FontSize',20);
    xlabel('Time [s]','FontSize',18);
    ylabel('Pressue [mmHg]','FontSize',18);
    hold on;

    % Attain the sinusoid fit for all points (so Pmax can be visualized
    for i = 1:length(EDP)
        % obtain the range of time of each peak
        interval = time(isovoltime(i).PosIso(1,1)):0.002:time(isovoltime(i).NegIso(end,1));

        % plug into Naeiji equation that was just solved for 
        FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + c_tot2(i,4));

        % find time point corresponding to Pmax
        [~, Idx] = min(abs(FitSinePres-P_max2(i)));

        PmaxT(i) = interval(Idx);

        plot(interval, FitSinePres, 'k--', PmaxT(i), P_max2(i), 'go');
        hold on;
    end

    % check the range of pressure values of Pmax. if the max p_max value is
    % over 450, rescale y axis to (0, 300), so individual waveforms can be seen
    maxP = max(P_max2);
    if maxP > 450
        ylim([0, 300]);
    else
        ylim([0, abs(maxP)+5]);
    end
    legend('Pressure', 'Isovolumic Points', 'Sinusoid Fit', 'Pmax', 'Location','southoutside', 'Orientation', 'horizontal');
    box on;
    grid on;
    hold off;

    % store the wavefit, so output tracks which wave forms did not have a good
    % fit
    handles.OutVar(1).output = waveFit;

    % Store the pmax values
    handles.OutVar(2).output = P_max2;

    % store regression constants
    handles.OutVar(3).output = c_tot2;
end

% Update handles structure
guidata(hObject, handles);

% Set cursor back to normal
set(handles.figure1, 'pointer', 'arrow');
end