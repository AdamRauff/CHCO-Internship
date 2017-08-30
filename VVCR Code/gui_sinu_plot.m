function [PmaxT] = gui_sinu_plot (time, Pres, EDP, isovoltime, P_max2, ...
                       c_tot2, totIsoTimePoints, totIsoPresPoints, handles)

% plot pressure, sinusoid fits
axes(handles.pressure_axes);

h = plot(time,Pres,'b', ...
         totIsoTimePoints,totIsoPresPoints,'ro');

set(h, 'HitTest', 'off');
set(handles.pressure_axes,'ButtonDownFcn', ...
    @(hObject, eventdata)GraphCallBack(hObject, eventdata, handles));
set(handles.pressure_axes,'fontsize',12);

title('Sinusoidal Fitting','FontSize',20);
xlabel('Time [s]','FontSize',18);
ylabel('Pressue [mmHg]','FontSize',18);

hold on;

PmaxT = zeros(length(EDP),1);

% Attain the sinusoid fit for all points (so Pmax can be visualized
for i = 1:length(EDP)

    % obtain the range of time of each peak
    interval =  ...
      time(isovoltime(i).PosIso(1,1)):0.002:time(isovoltime(i).NegIso(end,1));

    % plug into Naeiji equation that was just solved for
    FitSinePres = c_tot2(i,1) + c_tot2(i,2)*sin(c_tot2(i,3)*interval + ...
      c_tot2(i,4));

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

legend('Pressure', 'Isovolumic Points',  'Sinusoid Fit','Pmax', ...
    'Location','southoutside', 'Orientation', 'horizontal');

box on;
grid on;
hold off;
