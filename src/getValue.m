function [number,showMessage] = getValue(low, up, default, hObject, errormsg, show, opts)
% this function delivers the number of the edit text input --> number
% there are 2 types of message boxes: error and warning, both will be shown only once
% input parameters:
%     low:      lower input limit
%     up:       upper input limit
%     default:  edit text will be set to default if input is not correct
%     hObject:  leave this unchanged
%     errormsg: contains errormessages in case the lower or the upper limit is exceeded
%     show:     if show == 1 --> message box will be shown
%     opts:     in case the errormessage contains an index, the interpreter should be changed to latex
% -------------------------- default inputs - begin -----------------------------------
% default error/warning options
    if nargin == 4
        errormsg = ['Choose a value between ' num2str(low) ' and ' num2str(up) '.'];
        show.error = 1; show.warning = 1;
        opts = struct('WindowStyle','non-modal','Interpreter','none');
    elseif nargin == 5
        show.error = 1; show.warning = 1;
        opts = struct('WindowStyle','non-modal','Interpreter','none');
    elseif nargin == 6
        opts = struct('WindowStyle','non-modal','Interpreter','none');
    end
    if isempty(errormsg)
        errormsg = ['Choose a value between ' num2str(low) ' and ' num2str(up) '.'];
    end

% -------------------------- default inputs - end -------------------------------------

str = get(hObject,'string'); [number, status]=str2num(str);
% status = 1 --> string contains a number
% status = 0 --> string is not a number
    if status == 1 && isempty(number) == 0
        if length(number) == 2
            if contains(str,',')
                % --- ',' will be changed by '.' --- %
                idx = strfind(str,',');
                str(idx) = '.';
                number = str2double(str);
            else
                if show.error == 1
                    waitfor(errordlg('Please enter a number.','Error'))
                    showMessage.error = 0;
                end
                number = default;
            end
        elseif length(number) > 2
            if show.error == 1
                waitfor(errordlg('Please enter a number.','Error'))
                showMessage.error = 0;
            end
            number = default;
        end
        if number < low
            if show.warning == 1
                waitfor(warndlg(errormsg,'Attention',opts))
                showMessage.warning = 0;
            end
            number = low;
        elseif number > up
            if show.warning == 1
                waitfor(warndlg(errormsg,'Attention',opts))
                showMessage.warning = 0;
            end
            number = up;
        end
    else
        if show.error == 1
            waitfor(errordlg('Please enter a number.','Error'))
            showMessage.error = 0;
        end
        number = default;
    end
    if exist('showMessage','var')
        if isfield(showMessage,'error') == 0
            showMessage.error = [];
        elseif isfield(showMessage,'warning') == 0
            showMessage.warning = [];
        end
    else
        showMessage.error = [];
        showMessage.warning = [];
    end
    set(hObject,'string',num2str(number,12));
end