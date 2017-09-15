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
    % this pulls data out of PeakStruct as passed to GUI_SINU_FIT
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
        
        % AR 6/5/17
        % adding points succesively to beginning of systole to make better fit
        % of sick patients with wide curves
        % PresMax = max(double(Pres(pksT(i):MinIdx(i))));
        % if r_square > 0.90 && P_max2 < PresMax
            % add point to isovoltime(i).PosIso and corresponding isovol(i).PosIso
        % end
    
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
