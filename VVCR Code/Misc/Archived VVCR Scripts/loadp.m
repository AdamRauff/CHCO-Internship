function [Pres, dPdt, Rvals, nam, mrn, file, npath, marks] = ...
    loadp(npath, file, sm)

debug = 0;
marks = 0;

if debug == 1
    disp(['loadp pass ' num2str(sm)]);
end

if isa(file(1),'char')
    filename = strcat(npath, file);   % Get the full filename
    redname = basename(filename, 3);  % Get its basename for reporting
    fd0 = fopen(filename);            % open the file descriptor

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    tmp = textscan(fd0, '%s %s %s %s %s %s',1);
    nam = tmp{3};
    tmp = textscan(fd0, '%s %s %s',1);
    mrn = tmp{3};
    tmp = textscan(fd0, '%s %s %s',1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    lvppos = 0;
    ecgpos = 0;
    while ~strcmp (tmp{2},'Begin')
        [tmp, ~] = textscan(fd0, '%s %s %s %s',1);
    end
    tmp = textscan(fd0, '%s %s %s %s %s %s %s %s',1);
    for i = 1 : 2 : 7
        if strcmp(tmp{i},'RV')
            lvppos = i;  
        end
        if strcmp(tmp{i},'I')
            if ~strcmp(tmp{i},'III')
                ecgpos = -i;
            else
                ecgpos = i;
            end
        end
    end
 
    DataALL = textscan(fd0, '%s %s %s %s %s %s %s %s');  % Read all lines
    
    % check to see "RV" data is available. some of the text files don't
    % have an "RV" array
    if lvppos == 0 || isempty(lvppos)
        
        % set all outputs to 0
        Pres = 0;
        dPdt = 0; 
        Rvals = 0;
        nam = 0;
        mrn = 0;
        file = 0;
        npath = 0;
        marks = 0; 
        
        % exit the function
        return
    end
    
    tmp1 = DataALL{lvppos};
    tmp2 = DataALL{lvppos+1};
    Pres = zeros(length(tmp1),1);
    dPdt = zeros(length(tmp1),1);

    Pre_trig = [];
       
    for i = 1 : 1 : length(tmp1)
        temp1 = textscan(char(tmp1(i)), '%f(%c)');
        temp2 = textscan(char(tmp2(i)), '%f');
        Pres(i) = temp1{1};
        if ~isempty(temp2{1});
            dPdt(i) = temp2{1};
        end    
        if ~isempty(temp1{2})
            if strcmp(temp1{2},'S')
                Pre_trig = [Pre_trig i];
                marks = 1;
                if debug, disp(['got an S ' num2str(i)]), end;
            elseif strcmp(temp1{2},'D')
                Pre_trig = [Pre_trig -i];
                marks = 1;
                if debug, disp(['got an D ' num2str(i)]), end;
            end
        end
    end

    if sm == 1
        disp('    Using (S), (D) marks to gate on load');
        tmp = find(Pre_trig<0);
        Rvals = abs(Pre_trig(tmp));
        fclose(fd0);
        disp(['Witt Scientific Pressure File: [...]' redname ' loaded.']);
        if debug
            plot(Pres); hold on;
            cycMN = abs(Pre_trig(tmp));
            cycMX = Pre_trig(find(Pre_trig>0));
            plot(cycMX, Pres(cycMX),'bo','MarkerSize',5,'MarkerFaceColor','b');
            plot(cycMN, Pres(cycMN),'go','MarkerSize',5,'MarkerFaceColor','g');
            pause;
        end
        return;
    elsif sm == 2
        disp('    Pressure smoothing on load');
        time = [1:1:length(Pres)]*0.004;
        Pres = csaps(time,Pres,0.999995,time);
        dPdt = diff(Pres)/0.004;
        dPdt(length(Pres)) = 0;
    end
        
    if ecgpos > 0

        tmp3 = DataALL{ecgpos};
        tmp4 = DataALL{ecgpos+1};

        ECG  = zeros(length(tmp1),1);
        ECGp = zeros(length(tmp1),1);
        
        ECG_trig = [];
        
        % Loop over the 2,3 strings, looking for triggers
        for i = 1 : 1 : length(tmp1)
            temp3 = textscan(char(tmp3(i)), '%f');
            temp4 = textscan(char(tmp4(i)), '%f(%c)');
            ECG(i)  = temp3{1};
            ECGp(i) = temp4{1};
            if ~isempty(temp4{2})
                ECG_trig = [ECG_trig i];
                %disp(['got an R ' num2str(i)]);
            end
        end
    else
        peaks = []; percent = 0.8;
        while length(peaks) < 10
            mxdP = max(dPdt);
            peaks = find(dPdt>percent*mxdP);
            mndP = min(dPdt);
            valls = find(dPdt<percent*mndP);
            percent = percent-0.1;
        end
        
        cycMN(1) = valls(1);
        
        tmpv = valls;
        ccfind = 1;
        while ccfind > 0
            tmp = find(peaks>tmpv(1));         
            if ~isempty(tmp)
                cycMX(ccfind) = peaks(tmp(1));
                ccfind = ccfind+1;
                tmp = find(tmpv>peaks(tmp(1)));
                if isempty(tmp)
                    ccfind = -ccfind;
                end
                tmpv = tmpv(tmp);
                if ccfind == 2
                    valls = tmpv;
                end
            else
                ccfind = -ccfind;
            end
        end
        
        tmpp = peaks;
        ccfind = 2;
        while ccfind > 0
            tmp = find(valls>tmpp(1));
            if ~isempty(tmp)
                cycMN(ccfind) = valls(tmp(1));
                ccfind = ccfind+1;
                tmp = find(tmpp>valls(tmp(1)));
                if isempty(tmp)
                    ccfind = -ccfind;
                end
                tmpp = tmpp(tmp);
            else
                ccfind = -ccfind;
            end
        end
        
        % pre - allocate
        ECG_trig = zeros(length(cycMX),1);
        
        for i = 1: 1: length(cycMX)
            ECG_trig(i) = cycMN(i)+round(0.65*(cycMX(i)-cycMN(i)));
        end
        
        if debug == 1
            plot(Pres); hold on;
            plot(cycMX, Pres(cycMX),'bo','MarkerSize',5,'MarkerFaceColor','b');
            plot(cycMN, Pres(cycMN),'go','MarkerSize',5,'MarkerFaceColor','g');
            plot(ECG_trig, Pres(ECG_trig),'ro','MarkerSize',5,'MarkerFaceColor','r');
            pause;
        end
    end
    
    Rvals = ECG_trig;
    
    fclose(fd0);
    disp(['Witt Scientific Pressure File: [...]' redname ' loaded.']);
end
