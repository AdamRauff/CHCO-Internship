function calculate_Callback(hObject, ~, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% NOTE THAT THE UNDO CALLBACK ALSO HAS FITTING FUNCTIONALITY!!! UGH

%Make the cusor a spinning wheel so user is aware program is busy
set(handles.figure1, 'pointer', 'watch');
drawnow;


%%% obtain variables from InVar Struct for a clear workflow
% time and pressure vectors
time = handles.InVar(1).Data;
Pres = handles.InVar(2).Data;
Oldtime = handles.InVar(2).ivt;

% extract isovolmic points (times). These are structures
isovoltime = handles.InVar(1).ivt;

% extract isovolmic points (pressures). These are structures
isovol = handles.InVar(1).iv;

% iso volumic points in array (for plotting)
totIsoTimePoints = handles.InVar(1).isoPts;
totIsoPresPoints = handles.InVar(2).isoPts;

% EDP - end diastolic pressure
EDP = handles.InVar(1).Misc;

pksT = handles.InVar(1).Crit;
MinIdx = handles.InVar(2).Crit;

%%%%%% PASS ICS TO isovol_fit either as null or as these.
% calculate sinusoids based on new ICs!
%%% obtain new ICs
Mea = str2double(get(handles.Mean_txt,'String'));
Amp = str2double(get(handles.Amp_txt,'String'));
Fre = str2double(get(handles.Freq_txt,'String'));
Pha = str2double(get(handles.Phase_txt,'String'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre - allocate 
c_tot2 = zeros(length(EDP),4);
P_max2 = zeros(length(EDP),1);
waveFit = zeros(length(EDP),1);
r_square2 = zeros(length(EDP),1);

PmaxT = zeros(length(EDP),1); % UNIQUE
WHILE_LOOP_FLAG = zeros(length(EDP),1); % UNIQUE
ADD_TPoints = []; % UNIQUE
ADD_PPoints = []; % UNIQUE

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

    % AR 6/5/17
    % adding points succesively to beginning of systole to make better fit
    % of sick patients with wide curves
    % obtain maximum pressure point on actual curve
%     P2 = find(round(time,3)==round(Oldtime(MinIdx(i),3)));
%     P1 = find(round(time,3)==round(Oldtime(pksT(i)),3));
    PresMax = max(Pres(isovoltime(i).PosIso(1,1):1:isovoltime(i).NegIso(end,1)));
    if r_square2(i) > 0.80 && P_max2(i) < PresMax
        
        % keep count of how many points added to systole side
        count = 0;
        
        temp_ADD_TPoints = [];
        temp_ADD_PPoints = [];
        while P_max2(i) < PresMax
            
            % add point to isovoltime(i).PosIso and corresponding isovol(i).PosIso
            isovoltime(i).PosIso = [(isovoltime(i).PosIso(1,1))-1; isovoltime(i).PosIso];
            isovol(i).PosIso = [Pres(isovoltime(i).PosIso(1,1)); isovol(i).PosIso];
            
            temp_ADD_TPoints = [ADD_Points, (isovoltime(i).PosIso(1,1))-1]; 
            temp_ADD_PPoints = [ADD_PPoints, PresD(isovoltime(i).PosIso(1,1))];
            
            % update Wave(x)s variables
            WaveTs = [time(isovoltime(i).PosIso)'; time(isovoltime(i).NegIso)'];
            WavePs = [isovol(i).PosIso; isovol(i).NegIso];
            
            % mark flag that points are added
            WHILE_LOOP_FLAG(i) = true;
            
            % re-fit sinusiod
            % equation from Naeiji et al, single beat method of VVC
            sin_fun2=@(P)(P(1)+P(2)*sin(P(3)*WaveTs+P(4)))-WavePs; 

            %least squares fitting
            [c,resnorm,~]=lsqnonlin(sin_fun2,c2); 

            Psine_RV2=(c(1)+c(2)*sin(c(3)*WaveTs+c(4)));

            % r^2 value
            r_square2(i)=1-resnorm/norm(Psine_RV2-mean(Psine_RV2))^2;

            % if the fit of the wave was bad, mark that wave
            if r_square2(i) <0.90
               waveFit(i) = 1;
               WHILE_LOOP_FLAG(i) = false;
            else
               waveFit(i) = 0;
               WHILE_LOOP_FLAG(i) = true;
            end
            
            %getting all the c values in a matrix
            c_tot2(i,:)=c; 
    
            %first equation pmax, A+B
            P_max2(i)=c(1)+abs(c(2));
            
            % increment count to keep track of added points
            count = count +1;
            
            % Do not let program add more than 10 points
            if count >= 10 && (P_max2(i) < PresMax || waveFit(i) == 1)
                waveFit(i) = 1;
                disp('Added nine points on systolic side of curve, and curve fit remains unsatisfying');
                disp(['Wave: ',num2str(i), 'is excluded']);
                WHILE_LOOP_FLAG(i) = false;
                break
            end
        end
        
        if WHILE_LOOP_FLAG(i) == true
            ADD_TPoints = [ADD_TPoints, temp_ADD_TPoints];
            ADD_PPoints = [ADD_PPoints, temp_ADD_PPoints];
        end
    end
    
    % -------------------------------------------------------------
    % NOTE the absolute value of the amplitude is taken!!!!!!!
    % refer to patient HA002019, Wave 11 for an example of why

    % sometime amplitude of given equation solves for negative ( with a
    % significant phase shift, which makes a good fit (r^2 > 0.99).
    % ---------------------------------------------------------------
end

% if iso points have been added, re-compose the totIsoPnts variables
if any(WHILE_LOOP_FLAG)

    % update handles global variable
    % extract isovolmic points (times). These are structures
    handles.InVar(1).ivt = isovoltime;

    % extract isovolmic points (pressures). These are structures
    handles.InVar(1).iv = isovol;

    % recompose totIsoTimePoints and totIsoPresPoints and 
    totIsoTimePoints = [totIsoTimePoints, ADD_TPoints];
    totIsoPresPoints = [totIsoPresPoints, ADD_PPoints];
    
    % update global structure
    handles.InVar(1).isoPts = totIsoTimePoints;
    handles.InVar(2).isoPts = totIsoPresPoints;
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
