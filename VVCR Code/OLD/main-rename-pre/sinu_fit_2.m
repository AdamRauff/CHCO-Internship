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
        
%%%%% CAN THIS BE MADE INTO FUNCTION? IT'S ALSO IN VVCR_MULTIH
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
        
%%%%%% PASS ICS TO isovol_fit either as null or as these.
        % obtain current ICs
        Mea = str2double(get(handles.Mean_txt,'String'));
        Amp = str2double(get(handles.Amp_txt,'String'));
        Fre = str2double(get(handles.Freq_txt,'String'));
        Pha = str2double(get(handles.Phase_txt,'String'));
   
        ICS = [Mea, Amp, Fre, Pha];
        function [ PeakStruct ] = isovolPres_fit ( Iso1StVal, Iso1StVal_Doub, Iso2StIdx_Doub, TotNumWaves, time_end, isovolPres, isovolTime, time, timeDoub, dPmaxIdx, dPminIdx, Pres, PresDoub );        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pre - allocate 
        c_tot2 = zeros(length(EDP),4);
        P_max2 = zeros(length(EDP),1);
        waveFit = zeros(length(EDP),1);
        r_square2 = zeros(length(EDP),1);
        PmaxT = zeros(length(EDP),1); % UNIQUE can be moved outside
        
        totIsoTimePoints = [];
        totIsoPresPoints = [];
        
        for i = 1:length(EDP)
%%%%% HERE, time IS ALREADY = time_Doub
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
%%%%% THIS VERSION DOESN'T HAVE VANDERPOOL STUFF - WHICH IS FINE
            totIsoTimePoints = [totIsoTimePoints; WaveTs]; %UNIQUE
            totIsoPresPoints = [totIsoPresPoints; WavePs]; %UNIQUE
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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

            %%% UNIQUE %%%
            PmaxT(i) = interval(Idx);

            plot(interval, FitSinePres, 'k--', PmaxT(i), P_max2(i), 'go');
            hold on;
            %%% UNIQUE %%%
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
