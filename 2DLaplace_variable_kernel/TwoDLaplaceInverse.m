function varargout = TwoDLaplaceInverse(varargin)
%% main program
% TwoDLaplaceInverse V1.1
% Developers: Sophie Godefroy and Brett Ryland

% TWODLAPLACEINVERSE Application M-file for TwoDLaplaceInverse.fig
%    FIG = TWODLAPLACEINVERSE launch TwoDLaplaceInverse GUI.
%    TWODLAPLACEINVERSE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 27-Nov-2008 11:29:06

% Last modified 03-Dec-2008 by Marcel Gratz (Matlab 2008a)
% - fixed label representations in new figures
% - added support for multiple kernels
% - fixed expressions | and & to || and &&
% - optimized / adapted use of several variables
% - commented source code a bit
% - added cells (starting with %%) to code for easier navigation in Matlab
%   to use it: Menubar -> 'Cell' -> 'Enable Cell Mode'
% - added some waitbars, for more 'control' for the users

% contains state of buttons such as nnls-smoothing
% becomes 1 when clicking 'Draw' Button and everything is fine
global unlocked

if nargin == 0  % LAUNCH GUI
	fig = openfig(mfilename,'new');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
	unlocked = 0;

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


% Check user-defined kernel for correct syntax and usage
% if fine, then set the corresponding parameters
function success = CheckSyntaxAndAssign(h,alignment)
%% Check Syntax of kernel
% will be set globally, to indicate type of experiment
% type is determined with the help of 'T' and 'D', respectively, in kernel
% NOTE: axis scaling also depends on these variables
global hd hr vd vr kernel1 kernel2

success = 0;
% get all used variables

% horizontal kernel check
if strcmp(alignment,'h')

    kernel1=get(h.edit37,'string');
    variabs = argnames(inline(kernel1));
    % is 'D' used as variable in kernel string? Use result (0 or 1) as input for hd
    hd = not(isempty(cell2mat(strfind(variabs,'D'))));
    % is 'T' used as variable in kernel string? Use result (0 or 1) as input for hr
    hr = not(isempty(cell2mat(strfind(variabs,'T'))));

    % either both relaxation and diffusion or nothing
    if ((hd == 0) && (hr == 0)) || ((hd == 1) && (hr == 1)) 
        warndlg('use EITHER diffusion OR relaxation for horizontal process');
        return
    end
    % the 'axisdata' variable is not used in kernel string
    if isempty(cell2mat(strfind(variabs,'h')))
       warndlg('you have to use ''h'' as representative for your horizontal axis');
       return
    end
    % we use exactly 2 variables in all kernels
    % either (h and T) or (h and D)
    if not(length(cell2mat(variabs)) == 2)
        warndlg('degree of freedom has to be 2 for horizontal process');
        return
    end
    % check for the substring 'exp'
    % Since we have inverse Laplace trafo, we should use exponential forms, shouldn't we??? :)
    if isempty(findstr(kernel1,'exp'))
        warndlg('you need exponential behavior for horizontal kernel');
        return
    end
    
    try
        T=1; D=1; h=1; % use sample data set to evaluate the formula
        evalc(kernel1);
    catch
        warndlg(lasterr);
        return
    end
else  % same stuff as above for vertical kernel

    kernel2=get(h.edit38,'string');
    variabs = argnames(inline(kernel2));
    vd = not(isempty(cell2mat(strfind(variabs,'D'))));
    vr = not(isempty(cell2mat(strfind(variabs,'T'))));
    % either both relaxation and diffusion or nothing
    if ((vd == 0) && (vr == 0)) || ((vd == 1) && (vr == 1)) 
        warndlg('use EITHER diffusion OR relaxation for vertical process');
        return
    end
    if isempty(cell2mat(strfind(variabs,'v')))
       warndlg('you have to use ''v'' as representative for your vertical axis');
       return
    end
    if not(length(cell2mat(variabs)) == 2)
        warndlg('degree of freedom has to be 2 for vertical process');
        return
    end    
    if isempty(findstr(kernel2,'exp'))
        warndlg('you need exponential behavior for vertical kernel');
        return
    end

    try
        T=1; D=1; v=1; % use sample data set to evaluate the formula
        evalc(kernel2);
    catch
        warndlg(lasterr);
        return
    end
    
end

% we only get here, if everything went well
success = 1;


% --------------------------------------------------------------------
function varargout = radiobutton7_Callback(h, eventdata, handles, varargin)
%% choice for a 1D experiment
% Stub for Callback of the uicontrol handles.radiobutton7.

global d1

if CheckSyntaxAndAssign(handles,'h') 
    d1 = 1;
    set(handles.radiobutton7,'value',1);
    set(handles.radiobutton8,'value',0);
    activate(handles);
end


% --------------------------------------------------------------------
function varargout = radiobutton8_Callback(h, eventdata, handles, varargin)
%% choice for a 2D experiment
% Stub for Callback of the uicontrol handles.radiobutton8.

global d2

if CheckSyntaxAndAssign(handles,'v') && CheckSyntaxAndAssign(handles,'h')
    d2 = 1;
    set(handles.radiobutton8,'value',1);
    set(handles.radiobutton7,'value',0);
    activate(handles);
end


% --------------------------------------------------------------------
% enable/disable various buttons depending on what is selected
function activate(handles)
%% activate handles

global unlocked

% ----------------
% current assignment of htype and vtype
% 1 - Diffusion
% 2 - T1-Relaxation
% 3 - T1-Saturation
% 4 - T2-Relaxation
% 5 - User-Defined
% ----------------

htype = get(handles.popupmenu1,'value');
vtype = get(handles.popupmenu2,'value');

if get(handles.radiobutton8,'value') == 1 % 2D experiment
	set(handles.popupmenu2,'enable','on');
    % Compare string instead of value for better handling when adding new
    % standard functions
    str = get(handles.popupmenu2,'String');
    if (strcmp(str(vtype),'User-defined...') == 1) 
    	set(handles.edit38,'enable','on');
    else
    	set(handles.edit38,'enable','off');
    end

    set(handles.edit3,'enable','on');
	set(handles.edit10,'enable','on');
	if unlocked == 1
		set(handles.edit16,'enable','on');
		set(handles.edit17,'enable','on');
		set(handles.edit18,'enable','on');
		set(handles.edit21,'enable','on');
		set(handles.checkbox2,'enable','on');
		set(handles.checkbox3,'enable','on');
		set(handles.checkbox9,'enable','on');
		set(handles.pushbutton12,'enable','on'); 
	end
	set(handles.edit36,'enable','on');
	set(handles.radiobutton2,'enable','on');
	set(handles.pushbutton3,'enable','on');
end
if get(handles.radiobutton7,'value') == 1 % 1D experiment
	set(handles.popupmenu2,'enable','off');
	set(handles.edit38,'enable','off');
	set(handles.edit3,'enable','off');
	set(handles.edit10,'enable','off');
	set(handles.edit16,'enable','off');
	set(handles.edit17,'enable','off');
	set(handles.edit18,'enable','off');
	set(handles.edit21,'enable','off');
	set(handles.edit36,'enable','off');
	set(handles.radiobutton2,'enable','off');
	set(handles.pushbutton3,'enable','off');
	set(handles.checkbox2,'enable','off');
	set(handles.checkbox3,'enable','off');
	set(handles.pushbutton12,'enable','off');
	set(handles.checkbox9,'enable','off');
end
set(handles.edit2,'enable','on');
set(handles.radiobutton1,'enable','on');
set(handles.edit9,'enable','on');
if (htype == 1) || ((vtype == 1) && (get(handles.radiobutton8,'value') == 1)) % diffusion experiment
	set(handles.checkbox4,'enable','on');
	set(handles.checkbox5,'enable','on');
	if (get(handles.checkbox4,'value') == 1) || (get(handles.checkbox5,'value') == 1) % generate q^2
		if(htype == 1)
			set(handles.edit2,'enable','off'); % time 1
			set(handles.radiobutton1,'enable','off');
			set(handles.edit9,'enable','off');
		end
		if((vtype == 1) && (get(handles.radiobutton8,'value') == 1))
			set(handles.edit3,'enable','off'); % time 2
			set(handles.radiobutton2,'enable','off');
			set(handles.edit10,'enable','off');
		end
		set(handles.edit4,'enable','on');
		set(handles.edit6,'enable','on');
		if get(handles.checkbox5,'value') == 1 % gradvar
			set(handles.edit34,'enable','on');
			set(handles.edit33,'enable','on');
			set(handles.edit5,'enable','off');
		else                                                      % delvar
			set(handles.edit5,'enable','on');
			set(handles.edit34,'enable','off');
			set(handles.edit33,'enable','off');
		end
	else
		set(handles.edit2,'enable','on'); % time 1
		set(handles.radiobutton1,'enable','on');
		set(handles.edit9,'enable','on');
		if(get(handles.radiobutton8,'value') == 1)
			set(handles.edit3,'enable','on'); % time 2
			set(handles.radiobutton2,'enable','on');
			set(handles.edit10,'enable','on');
		end
		set(handles.edit4,'enable','off');
		set(handles.edit5,'enable','off');
		set(handles.edit6,'enable','off');
		set(handles.edit34,'enable','off');
		set(handles.edit33,'enable','off');
	end
else
	set(handles.edit2,'enable','on'); % time 1
	if(get(handles.radiobutton8,'value') == 1)
		set(handles.edit3,'enable','on'); % time 2
	end
	set(handles.checkbox4,'enable','off');
	set(handles.checkbox5,'enable','off');
	set(handles.edit4,'enable','off');
	set(handles.edit5,'enable','off');
	set(handles.edit6,'enable','off');
	set(handles.edit34,'enable','off');
	set(handles.edit33,'enable','off');
end

% --------------------------------------------------------------------
% unlock all other buttons
function unlock(handles)
%% Unlock Handles

global unlocked

unlocked = 1;
set(handles.pushbutton5,'enable','on');
set(handles.pushbutton6,'enable','on');
set(handles.pushbutton8,'enable','on');
set(handles.pushbutton11,'enable','on');
set(handles.pushbutton12,'enable','on');
set(handles.pushbutton13,'enable','on');
set(handles.pushbutton14,'enable','on');
set(handles.pushbutton15,'enable','on');
set(handles.checkbox1,'enable','on');
set(handles.checkbox2,'enable','on');
set(handles.checkbox3,'enable','on');
set(handles.checkbox7,'enable','on');
set(handles.checkbox8,'enable','on');
set(handles.checkbox9,'enable','on');
set(handles.edit11,'enable','on');
set(handles.edit12,'enable','on');
set(handles.edit13,'enable','on');
set(handles.edit14,'enable','on');
set(handles.edit15,'enable','on');
set(handles.edit16,'enable','on');
set(handles.edit17,'enable','on');
set(handles.edit18,'enable','on');
set(handles.edit21,'enable','on');
set(handles.edit24,'enable','on');
set(handles.edit25,'enable','on');
set(handles.edit26,'enable','on');

% --------------------------------------------------------------------
% button to load data file
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
%% Load data file
% Stub for Callback of the uicontrol handles.pushbutton1.

global data filedata chemin nfig fichier OKh OKv videos videof dis trans

if not(CheckSyntaxAndAssign(handles,'h')), return, end
if not(CheckSyntaxAndAssign(handles,'v')), return, end

[fichier,chemin] = uigetfile('*.*');
set(handles.edit1,'string',fichier);
pwd2 = pwd;
cd(chemin);
data = load(fichier);
filedata=data;
cd(pwd2);

% initialization
OKh = 0;
OKv = 0;
nfig = 0;
set(handles.checkbox1,'value',0);
videos = 0;
videof = 0;
dis = 0;
trans = 0;
set(handles.checkbox8,'value',0);

% --------------------------------------------------------------------
% button to load "time" 1 (horizontal)
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)
%% load horizontal axis
% Stub for Callback of the uicontrol handles.pushbutton2.

global chemin filetimea timea fi1


pwd2 = pwd;
cd(chemin);
[fi1,chemint1] = uigetfile('*.*');
cd(chemint1);
timea = load(fi1);
filetimea = timea;
cd(pwd2);
set(handles.edit2,'string',fi1);
    

% --------------------------------------------------------------------
% button to load "time" 2 (vertical)
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)
%% load vertical axis
% Stub for Callback of the uicontrol handles.pushbutton3.

global chemin filetimeb timeb fi2

pwd2 = pwd;
cd(chemin);
[fi2,chemint2] = uigetfile('*.*');
cd(chemint2);
timeb = load(fi2);
filetimeb = timeb;
cd(pwd2);
set(handles.edit3,'string',fi2);


% --------------------------------------------------------------------
% edition of the horizontal size of the data 
function varargout = edit7_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit7.
% --------------------------------------------------------------------
% edition of the vertical size of the data
function varargout = edit8_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit8.
% --------------------------------------------------------------------
% choice button for manual "time" horizontal
function varargout = radiobutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton1.

global cmth

cmth = get(h,'value');


% --------------------------------------------------------------------
% choice button for manual "time" vertical
function varargout = radiobutton2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton2.

global cmtv

cmtv = get(h,'value');


% --------------------------------------------------------------------
% edition of the data file
function varargout = edit1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit1.
% --------------------------------------------------------------------
% edition of the "time" file 1
function varargout = edit2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit2.
% --------------------------------------------------------------------
% edition of the "time" file 2
function varargout = edit3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit3.

% --------------------------------------------------------------------
% T2 display checkbox for the T1 experiments
% function varargout = checkbox6_Callback(h, eventdata, handles, varargin)
% % Stub for Callback of the uicontrol handles.checkbox6.
% 
% global T1T2dis
% 
% T1T2dis = get(h,'value');


% --------------------------------------------------------------------
% draw data button
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)
%% Draw data
% Stub for Callback of the uicontrol handles.pushbutton4.

global cmth cmtv d1 d2 hr hd vr vd timea timeb data filetimea filetimeb filedata 
global Gmax delta DELTA nfig gradvar delvar
global Store_Thmin Store_Thmax Store_hsteps Store_Tvmin Store_Tvmax Store_vsteps
%global vsize hsize OKh OKv T1T2dis dis % currently not used


d1 = get(handles.radiobutton7,'value');
d2 = get(handles.radiobutton8,'value');
% T1T2dis = get(handles.checkbox6,'value');
cmth = get(handles.radiobutton1, 'value');

gradvar = get(handles.checkbox4,'value');
delvar = get(handles.checkbox5,'value');

Store_Thmin = 0.001;
Store_Thmax = 100;
Store_hsteps = 15;
Store_Tvmin = 0.001;
Store_Tvmax = 100;
Store_vsteps= 15;
colormap([[zeros(1,22) 0:.937/8:.937 .937:.02/10:.957 .482:-.398/30:.084]' [zeros(1,10) 0:.968/20:.968 .968:.012/10:.98 1:-.61/30:.39]' [.5:.5/10:1 1:-.2/20:.8 .4:.2/10:.6 .176:-.176/29:0]']);

unlock(handles);

timea = filetimea;
timeb = filetimeb;
data = filedata;

if not(CheckSyntaxAndAssign(handles,'h')), return, end

if d2 == 1
    if not(CheckSyntaxAndAssign(handles,'v')), return, end
    [sdb,sda] = size(data);
end
if d1 == 1
    [sda,sdb] = size(data);
end
set(handles.edit7,'string',sda);
set(handles.edit8,'string',sdb);

if (hr == 1) || ((hd == 1) && (gradvar == 0) && (delvar == 0))
    [sta,nimportequoi] = size(timea);
end
if (d2 == 1) && ((vr == 1) || ((vd == 1) && (gradvar == 0) && (delvar == 0)))
    [stb,nimportequoi] = size(timeb);
end

% if (d1 == 1) && (T1T2dis == 1) && (dis == 0)
%     dis = 1;
%     data = (max(abs(data(sda,1)),abs(data(1,1)))-data)/2;
% end

if cmth == 1
    tauvh = str2double(get(handles.edit9,'string'));
    timea = filetimea .* tauvh;
end
cmtv = get(handles.radiobutton2, 'value');
if cmtv == 1
    tauvv = str2double(get(handles.edit10,'string'));
    timeb = filetimeb .* tauvv;
end


% Check the validity of the data
if (hr == 1) || ((gradvar == 0) && (delvar == 0))
	if sta ~= sda
		warndlg('error: chose the good data and time set', '!!Sorry!!');
		return
	end
end
if (d2 == 1) && ((vr == 1) || ((gradvar == 0) && (delvar == 0)))
    if stb ~= sdb
        warndlg('error: chose the good data and time set', '!!Sorry!!');
		return
    end
end

if (hd == 1) || (vd == 1)
    Gmax = str2double(get(handles.edit4,'string'));
    delta = str2double(get(handles.edit5,'string'));
    DELTA = str2double(get(handles.edit6,'string'));
end

gradvar = get(handles.checkbox4,'value');
delvar = get(handles.checkbox5,'value');
firstdel = str2double(get(handles.edit33,'string'));
stepdel = str2double(get(handles.edit34,'string'));

if hd == 1
    if gradvar == 1
        timea = ((1:sda)/sda*Gmax*delta*10^-3*(2.67e8)).^2.*(DELTA-delta/3)*10^-3*10^-9;
        timea = timea';
    end
    if delvar == 1
        timea = (Gmax*((firstdel+(0:sda-1)*stepdel)*(10^-3))*(2.67e8)).^2.*(DELTA*10^-3-((firstdel+(0:sda-1)*stepdel)*(10^-3))/3)*10^-9;
        timea = timea';
    end
end
if vd == 1 && d2 == 1
    if gradvar == 1
        timeb = ((1:sdb)/sdb*Gmax*delta*10^-3*(2.67e8)).^2.*(DELTA-delta/3)*10^-3*10^-9;
        timeb = timeb';
    end
    if delvar == 1
        timeb = (Gmax*((firstdel+(0:sdb-1)*stepdel)*(10^-3))*(2.67e8)).^2.*(DELTA*10^-3-((firstdel+(0:sdb-1)*stepdel)*(10^-3))/3)*10^-9;
        timeb = timeb';
    end
end

if nfig == 1
    figure;
else
    axes(handles.axes1);
end
if d2 == 1
    surf(timea, timeb, data)
    va = timea; 
    vb = timeb;
    axis([va(1),va(sda),vb(1),vb(sdb)]);
	 shading interp;
    if nfig == 1
       colorbar;            
    end
end
if d1 == 1
    plot(timea,data)
    axis([timea(1),timea(sda),0,max(data)]); %modified
end

% --------------------------------------------------------------------
% tau value for the time 1
function varargout = edit9_Callback(h, eventdata, handles, varargin)
%% horizontal time multiplier
% Stub for Callback of the uicontrol handles.edit9.

global tauvh

tauvh = str2double(get(h,'string'));


% --------------------------------------------------------------------
% tau value for the time 2
function varargout = edit10_Callback(h, eventdata, handles, varargin)
%% vertical time multiplier
% Stub for Callback of the uicontrol handles.edit10.

global tauvv

tauvv = str2double(get(h,'string'));

% --------------------------------------------------------------------
% gradient experiment by the variation of the gradients
function varargout = checkbox4_Callback(h, eventdata, handles, varargin)
%% clicked Gradient calculation checkbox
% Stub for Callback of the uicontrol handles.checkbox4.

global gradvar delvar
%global hd vd % currently not used

gradvar = get(h,'value');
delvar = 0;
set(handles.checkbox5,'value',0);
activate(handles);

% --------------------------------------------------------------------
% gradient experiment by varying the duration of delta
function varargout = checkbox5_Callback(h, eventdata, handles, varargin)
%% clicked delta calculation checkbox
% Stub for Callback of the uicontrol handles.checkbox5.

global delvar gradvar 
%global hd vd % currently not used

delvar = get(h,'value');
gradvar = 0;
set(handles.checkbox4,'value',0);
activate(handles);

% --------------------------------------------------------------------
% value of G max
function varargout = edit4_Callback(h, eventdata, handles, varargin)
%% enter Gmax
% Stub for Callback of the uicontrol handles.edit4.

global Gmax

Gmax = str2double(get(h,'string'));


% --------------------------------------------------------------------
% value of delta (ms)
function varargout = edit5_Callback(h, eventdata, handles, varargin)
%% Enter delta
% Stub for Callback of the uicontrol handles.edit5.

global delta

delta = str2double(get(h,'string'));

% --------------------------------------------------------------------
% value of DELTA (ms)
function varargout = edit6_Callback(h, eventdata, handles, varargin)
%% Enter DELTA
% Stub for Callback of the uicontrol handles.edit6.

global DELTA

DELTA = str2double(get(h,'string'));

% --------------------------------------------------------------------
% first value of delta (ms)
function varargout = edit33_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit33.

% --------------------------------------------------------------------
%step of delta (ms)
function varargout = edit34_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit34.

% --------------------------------------------------------------------
% new figure checkbox
function varargout = checkbox1_Callback(h, eventdata, handles, varargin)
%% clicked new figure checkbox
% Stub for Callback of the uicontrol handles.checkbox1.

global nfig

nfig = get(h,'value');


% --------------------------------------------------------------------
% Function to remove rows or columns from the 2D data set
function varargout = pushbutton9_Callback(h, eventdata, handles, varargin)
%% remove columns and rows
% Stub for Callback of the uicontrol handles.pushbutton9.

global timea timeb data nfig d1 d2

if d2 == 1
    nbtop = str2double(get(handles.edit29,'string'));
    nbbottom = str2double(get(handles.edit30,'string'));
end
nbright = str2double(get(handles.edit31,'string'));
nbleft = str2double(get(handles.edit32,'string'));

if d2 == 1
    [sizeb,sizea] = size(data);
end
if d1 == 1
    [sizea,sizeb] = size(data);
end

if d2 == 1
    if nbbottom ~= 0
        data(1:nbbottom,:)=[];
        timeb(1:nbbottom)=[];
    end
    if nbtop ~= 0
        data(sizeb+1-nbtop:sizeb,:)=[];
        timeb(sizeb+1-nbtop:sizeb)=[];
    end
    if nbleft ~= 0
        data(:,1:nbleft)=[];
        timea(1:nbleft)=[];
    end
    if nbright ~= 0
        data(:,sizea+1-nbright:sizea)=[];
        timea(sizea+1-nbright:sizea)=[];
    end
end

if d1 == 1
    if nbleft ~= 0
        data(1:nbleft)=[];
        timea(1:nbleft)=[];
    end
    [sizea,sizeb] = size(data);
    if nbright ~= 0
        data(sizea+1-nbright:sizea)=[];
        timea(sizea+1-nbright:sizea)=[];
    end
end
        
if d2 == 1
    [sizeb,sizea] = size(data);
end
if d1 == 1
    [sizea,sizeb] = size(data);
end

set(handles.edit7,'string',sizea);
set(handles.edit8,'string',sizeb);
sda = sizea;
sdb = sizeb;
if nfig == 1
    figure;
else
    axes(handles.axes1);
end
if d2 == 1
    surf(timea, timeb, data)
    va = timea;
    vb = timeb;
    axis([va(1),va(sda),vb(1),vb(sdb)]);
    shading interp;
    if nfig == 1
       colorbar;            
    end
end
if d1 == 1
    plot(timea,data)
    xlim([timea(1),timea(sda)]);
end


% --------------------------------------------------------------------
% function to undo the removal of rows or columns from the data set
function varargout = pushbutton10_Callback(h, eventdata, handles, varargin)
%% undo removal
% Stub for Callback of the uicontrol handles.pushbutton10.

global timea timeb data filetimea filetimeb filedata nfig d1 d2 hr hd vr vd
 
if not(CheckSyntaxAndAssign(handles,'h')), return, end
if not(CheckSyntaxAndAssign(handles,'v')), return, end   %% correct this above and continue below!!!!

gradvar = get(handles.checkbox4,'value');
delvar = get(handles.checkbox5,'value');
Gmax = str2double(get(handles.edit4,'string'));
delta = str2double(get(handles.edit5,'string'));
DELTA = str2double(get(handles.edit6,'string'));

data = filedata;
if d2 == 1
	[sizeb,sizea] = size(data);
end
if d1 == 1
	[sizea,sizeb] = size(data);
end
set(handles.edit7,'string',sizea);
set(handles.edit8,'string',sizeb);
sda = sizea;
sdb = sizeb;

if (hr == 1) || ((hd == 1) && (gradvar == 0) && (delvar == 0))
	timea = filetimea;
else
    if gradvar == 1
        timea = ((1:sda)/sda*Gmax*delta*10^-3*(2.67e8)).^2.*(DELTA-delta/3)*10^-3*10^-9;
        timea = timea';
    end
    if delvar == 1
        timea = (Gmax*((firstdel+(0:sda-1)*stepdel)*(10^-3))*(2.67e8)).^2.*(DELTA*10^-3-((firstdel+(0:sda-1)*stepdel)*(10^-3))/3)*10^-9;
        timea = timea';
    end
end

if d2 == 1
	if (vr == 1) || ((vd == 1) && (gradvar == 0) && (delvar == 0))
		timeb = filetimeb;
	else
		if gradvar == 1
			timeb = ((1:sdb)/sdb*Gmax*delta*10^-3*(2.67e8)).^2.*(DELTA-delta/3)*10^-3*10^-9;
			timeb = timeb';
		end
		if delvar == 1
			timeb = (Gmax*((firstdel+(0:sdb-1)*stepdel)*(10^-3))*(2.67e8)).^2.*(DELTA*10^-3-((firstdel+(0:sdb-1)*stepdel)*(10^-3))/3)*10^-9;
			timeb = timeb';
		end
	end
end

if nfig == 1
    figure;
else
    axes(handles.axes1);
end
if d2 == 1
    surf(timea, timeb, data)
    va = timea;
    vb = timeb;
    axis([va(1),va(sda),vb(1),vb(sdb)]);
    shading interp;
    if nfig == 1
       colorbar;            
    end
end
if d1 == 1
    plot(timea,data)
    xlim([timea(1),timea(sda)]);
end



% --------------------------------------------------------------------
% number of top rows to remove
function varargout = edit29_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit29.
% --------------------------------------------------------------------
% number of bottom rows to remove
function varargout = edit30_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit30.
% --------------------------------------------------------------------
% number of right columns to remove
function varargout = edit31_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit31.
% --------------------------------------------------------------------
% number of left columns to remove
function varargout = edit32_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit32.


% --------------------------------------------------------------------
% --------------------------------------------------------------------
% ----------------- nnls - smoothing analysis ------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------


% --------------------------------------------------------------------
% button to perform and draw to distributions
function varargout = pushbutton5_Callback(h, eventdata, handles, varargin)
%% start NNLS smoothing
% Stub for Callback of the uicontrol handles.pushbutton5.

global data timea timeb chemin fichier fi1 fi2 spectrum tauv tauh 
global hr hd vr vd taulv taulh sta stb d1 d2 kernel1 kernel2
%global nfig spectrumi XI YI T1st T2st % unused currently

% "a" = "h" and "b" = "v"
% acquisition of the parameters of the nnls smoothing process
alpha = str2double(get(handles.edit11,'string'));
beta = str2double(get(handles.edit12,'string'));
Thmin = str2double(get(handles.edit13,'string'));
Thmax = str2double(get(handles.edit14,'string'));
Thmm = [Thmin,Thmax];
stepsh = str2double(get(handles.edit15,'string'));

if not(CheckSyntaxAndAssign(handles,'h')), return, end

if d2 == 1
    Tvmin = str2double(get(handles.edit16,'string'));
    Tvmax = str2double(get(handles.edit17,'string'));
    Tvmm = [Tvmin, Tvmax];
    stepsv = str2double(get(handles.edit18,'string'));
    if not(CheckSyntaxAndAssign(handles,'v')), return, end
end
orient = get(handles.edit21,'string');

if hd == 1
    a = Thmin;
    Thmin = 1/Thmax;
    Thmax = 1/a;
    Thmm = [Thmin,Thmax];
end
if vd == 1 && d2 == 1
    a = Tvmin;
    Tvmin = 1/Tvmax;
    Tvmax = 1/a;
    Tvmm = [Tvmin,Tvmax];
end

tic;

if d2 == 1
    [spectrum,tauh,tauv,chisq,compte]=upnnlsmooth3Dsvdfin(data,timea',timeb',Thmm,stepsh,Tvmm,stepsv,alpha,beta,orient,kernel1,kernel2);
end

if d1 == 1
    [spectrum,tauh,chisq,compte]=upnnlsmooth1D(data,timea,Thmin,Thmax,alpha,beta,stepsh,kernel1);
end
if hd == 1 
    spectrum = flipdim(spectrum,2);
    tauh = 1./tauh;
    tauh = flipdim(tauh,2);
end
if vd == 1 && d2 == 1
    spectrum = flipdim(spectrum,1);
    tauv = 1./tauv;
    tauv = flipdim(tauv,2);
end

if d2 == 1
    taulv = log10(tauv);
    stb = size(taulv);
end
taulh = log10(tauh);
sta = size(taulh);
% XI = taulh;
% if d2 == 1
%     YI = taulv;
% end
% T1st = sta;
% if d2 == 1
%     T2st = stb;
% end
savefile = [fichier '.out'];
pwd1 = pwd;
cd(chemin);
if d1 == 1
    spectrum = spectrum';
end
save(savefile,'spectrum','-ASCII');
if d1 == 1
    spectrum = spectrum';
end

gradvar = get(handles.checkbox4,'value');
delvar = get(handles.checkbox5,'value');

if (hr == 1) || ((hd == 1) && (gradvar == 0) && (delvar == 0))
    savefile = [fi1 '.out'];
    tauh = tauh';
    save(savefile,'tauh','-ASCII');
elseif hd == 1
    savefile = [fichier 'grad' '.out'];
    tauh = tauh';
    save (savefile,'tauh','-ASCII');
end
if d2 == 1
	if (vr == 1) || ((hd == 1) && (gradvar == 0) && (delvar == 0))
        savefile = [fi2 '.out'];
        tauv = tauv';
        save(savefile,'tauv','-ASCII');
    elseif vd == 1
        savefile = [fichier 'grad' '.out'];
        tauv = tauv';
        save(savefile,'tauv','-ASCII');
    end
end

cd(pwd1);

axes(handles.axes2);
if d2 == 1
    surf(taulh,taulv,spectrum)
	 if get(handles.checkbox9,'value')==0
		 caxis('auto');
	 else
		 caxis([0 1]);
	 end
    shading interp;
    set(gcf,'Renderer','zbuffer');
    axis([taulh(1),taulh(sta(2)),taulv(1),taulv(stb(2))]);
elseif d1 == 1
    plot(taulh,spectrum)
    xlim([taulh(1),taulh(sta(2))]);
end
t=toc
set(handles.edit22,'string', chisq);
set(handles.edit23,'string',t);
set(handles.edit20,'string',compte);


% --------------------------------------------------------------------
% alpha value
function varargout = edit11_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit11.
% --------------------------------------------------------------------
% beta value
function varargout = edit12_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit12.
% --------------------------------------------------------------------
% Tmin horizontal
function varargout = edit13_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit13.
% --------------------------------------------------------------------
% Tmax horizontal
function varargout = edit14_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit14.
% --------------------------------------------------------------------
% no. steps horizontal
function varargout = edit15_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit15.
% --------------------------------------------------------------------
% Tmin vertical
function varargout = edit16_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit16.
% --------------------------------------------------------------------
% Tmax vertical
function varargout = edit17_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit17.
% --------------------------------------------------------------------
% no. steps vertical
function varargout = edit18_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit18.
% --------------------------------------------------------------------
% number of loops
function varargout = edit20_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit20.
% --------------------------------------------------------------------
% orientation
function varargout = edit21_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit21.
% --------------------------------------------------------------------
% checkbox to get the transverse of the data
function varargout = checkbox8_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.checkbox8.



% --------------------------------------------------------------------
% redraw the distribution (for example) on another figure
function varargout = pushbutton6_Callback(h, eventdata, handles, varargin)
%% redraw
% Stub for Callback of the uicontrol handles.pushbutton6.

global spectrum taulv taulh nfig sta stb vd vr hd hr fichier videos videof chemin d1 d2

vd1 = vd;
vr1 = vr;
hd1 = hd;
hr1 = hr;
spectrum1 = spectrum;
trans = get(handles.checkbox8,'value');
taulv1 = taulv;
taulh1 = taulh;

if trans == 1
    if vd1 == hd1
    elseif vr1 == hr1
    else
        if vd1 == 1
            hd = 1;
            vd = 0;
        elseif hd1 == 1
            hd = 0;
            vd = 1;
        end
        if vr1 == 1
            vr = 0;
            hr = 1;
        elseif hr1 == 1
            hr = 0;
            vr = 1;        
        end
    end

    spectrum = spectrum'
    quelquechose = taulv;
    taulv = taulh;
    taulh = quelquechose;
end


if nfig == 1 && d2 == 1
    fig=figure;
    axes('FontSize',12);
    set(gcf,'Renderer','zbuffer');
    set(fig,'DoubleBuffer','on');
    set(gca,'NextPlot','replace','Visible','off')
    surf(taulh,taulv',spectrum);
	 if get(handles.checkbox9,'value')==0
		 caxis('auto');
	 else
		 caxis([0 1]);
	 end
    shading interp;
    axis([taulh(1),taulh(sta(2)),taulv(1),taulv(stb(2))]);
    colorbar;
    if hr == 1
        xlabel('Log_{10}(T) [s]','FontSize',14);
    elseif hd == 1
        xlabel('Log_{10}(D) [10^{-9} m^2/s]','FontSize',14);
    end
    if vr == 1 
        ylabel('Log_{10}(T) [s]','FontSize',14);
    elseif vd == 1
        ylabel('Log_{10}(D) [10^{-9} m^2/s]','FontSize',14);
    end
    if hd == 1 && vd == 1
        title('D-D correlation','FontSize',16);
    elseif hr == 1 && vr == 1
        title('T-T correlation','FontSize',16);
    else
        title('D-T correlation','FontSize',16);
    end
end

if nfig == 1 && d1 == 1
    fig = figure;
    axes('FontSize',20);
    plot(taulh,spectrum,'k')
    xlim([taulh(1),taulh(sta(2))]);
    if hr == 1
        xlabel('Log_{10}(T) [s]','FontSize',14);
    elseif hd == 1
        xlabel('Log_{10}(D) [10^{-9} m^2/s]','FontSize',14);
    end
end    

pwd1 = pwd;
cd(chemin);
savefile = [fichier 'm.fig'];
saveas(gcf,savefile);
savefile = [fichier 'm.jpg'];
saveas(gcf,savefile);


if nfig == 1 && videos == 1
    surf(taulh,taulv',spectrum)
	 if get(handles.checkbox9,'value')==0
		 caxis('auto');
	 else
		 caxis([0 1]);
	 end
	 shading interp;
    set(gcf,'Renderer','zbuffer');
    axis([taulh(1),taulh(sta(2)),taulv(1),taulv(stb(2)),min(min(spectrum)),max(max(spectrum))]);
    view([0 90]);
    set(gca,'visible','off');
    savefile = [fichier 's.avi'];
    mov = avifile(savefile, 'quality',100);
    F = getframe(gca);
    mov = addframe(mov,F);
    for i=1:2:90
        camorbit(2.19,-1.01,'data',[0.5 0.5 1.2]);
        drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);
    end
    mov = close(mov);

end

if nfig == 1 && videof == 1
    surf(taulh,taulv',spectrum)
    set(gcf,'Renderer','zbuffer');
	 if get(handles.checkbox9,'value')==0
		 caxis('auto');
	 else
		 caxis([0 1]);
	 end
    shading interp;
    axis([taulh(1),taulh(sta(2)),taulv(1),taulv(stb(2)),min(min(spectrum)),max(max(spectrum))]);
    view([0 90]);
    set(gca,'visible','off');
    savefile = [fichier 'f.avi'];
    mov = avifile(savefile, 'quality',100);
    F = getframe(gca);
    mov = addframe(mov,F);
    for i=1:2:90
        camorbit(0.0,-2.0,'data',[0 0.5 1]);
        drawnow;
        F = getframe(gca);
        mov = addframe(mov,F);
    end
    mov = close(mov);
    savefile = [fichier 'm.avi'];
end
cd(pwd1);
vd = vd1;
vr = vr1;
hd = hd1;
hr = hr1;
spectrum = spectrum1;
taulv = taulv1;
taulh = taulh1;


% --------------------------------------------------------------------
% edition of chisq
function varargout = edit22_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit22.
% --------------------------------------------------------------------
% elapse time of the nnls experiment
function varargout = edit23_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit23.



% --------------------------------------------------------------------
% ------------------------ loop on alpha -----------------------------
% --------------------------------------------------------------------

% --------------------------------------------------------------------
% Loop on alpha
function varargout = pushbutton8_Callback(h, eventdata, handles, varargin)
%% Alpha Loop
% Stub for Callback of the uicontrol handles.pushbutton8.

global data timea timeb spectrum tauv tauh hd vd d1 d2 kernel1 kernel2
%global nfig chemin fichier fi1 fi2 spectrumi XI YI T1st T2st  % currently not used

alphaf = str2double(get(handles.edit24,'string'));
alphal = str2double(get(handles.edit25,'string'));
alphanb = str2double(get(handles.edit26,'string'));

beta = str2double(get(handles.edit12,'string'));
Thmin = str2double(get(handles.edit13,'string'));
Thmax = str2double(get(handles.edit14,'string'));
Thmm = [Thmin,Thmax];
stepsh = str2double(get(handles.edit15,'string'));

if not(CheckSyntaxAndAssign(handles,'h')), return, end

if d2 == 1
    Tvmin = str2double(get(handles.edit16,'string'));
    Tvmax = str2double(get(handles.edit17,'string'));
    Tvmm = [Tvmin, Tvmax];
    stepsv = str2double(get(handles.edit18,'string'));
    if not(CheckSyntaxAndAssign(handles,'v')), return, end
end
orient = get(handles.edit21,'string');

if hd == 1
    a = Thmin;
    Thmin = 1/Thmax;
    Thmax = 1/a;
    Thmm = [Thmin,Thmax];
end
if vd == 1 && d2 == 1
    a = Tvmin;
    Tvmin = 1/Tvmax;
    Tvmax = 1/a;
    Tvmm = [Tvmin,Tvmax];
end

alphapas = (log10(alphal)-log10(alphaf))/(alphanb-1);
alphav = exp((log10(alphaf) : alphapas : log10(alphal))*log(10));

wbar=waitbar(0);
for na = 1 : alphanb
    waitbar(na./alphanb,wbar,strcat(['Performing step ',num2str(na),' of ',num2str(alphanb),'...']));
    set(handles.edit27,'string',na);
    tic;
    if d2 == 1
        [spectrum,tauh,tauv,chisq,compte]=upnnlsmooth3Dsvdfin(data,timea',timeb',Thmm,stepsh,Tvmm,stepsv,alphav(na),beta,orient,kernel1,kernel2);
    elseif d1 == 1
        [spectrum,tauh,chisq,compte]=upnnlsmooth1D(data,timea,Thmin,Thmax,alphav(na),beta,stepsh,kernel1);
    end
    chiv(na) = chisq;
    if hd == 1
        spectrum = flipdim(spectrum,2);
        tauh = 1./tauh;
        tauh = flipdim(tauh,2);
    end
    if vd == 1 && d2 == 1
        spectrum = flipdim(spectrum,1);
        tauv = 1./tauv;
        tauv = flipdim(tauv,2);
    end
    if d2 == 1
        taulv = log10(tauv);
        stb = size(taulv);
    end
    taulh = log10(tauh);
    sta = size(taulh);
    dra = get(handles.checkbox7,'value');
    if dra == 1
        figure
        if d2 == 1
            surf(taulh,taulv,spectrum);
				if get(handles.checkbox9,'value')==0
					caxis('auto');
				else
					caxis([0 1]);
				end
            shading interp;
            axis([taulh(1),taulh(sta(2)),taulv(1),taulv(stb(2))]);
            set(gcf,'Renderer','zbuffer');
            colorbar;
        end
        if d1 == 1
            plot(taulh,spectrum);
        end
    end
    t=toc
    set(handles.edit22,'string', chisq);
    set(handles.edit23,'string',t);
    set(handles.edit20,'string',compte);
end
close(wbar);
figure;
semilogx(alphav,chiv,'o')

% --------------------------------------------------------------------
% checkbox to draw all the distributions as a function of alpha, for the alpha loop
function varargout = checkbox7_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.checkbox7.
% --------------------------------------------------------------------
% first alpha 
function varargout = edit24_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit24.
% --------------------------------------------------------------------
% last alpha
function varargout = edit25_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit25.
% --------------------------------------------------------------------
% nb of alpha
function varargout = edit26_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit26.
% --------------------------------------------------------------------
% edition of the nbof the alpha
function varargout = edit27_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit27.



% --------------------------------------------------------------------
% Video side
function varargout = checkbox2_Callback(h, eventdata, handles, varargin)
%% VideoS Checkbox
% Stub for Callback of the uicontrol handles.checkbox2.

global videos

videos = get(h,'value');

% --------------------------------------------------------------------
% Video front
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)
%% VideoF Checkbox
% Stub for Callback of the uicontrol handles.checkbox3.

global videof

videof = get(h,'value');

% --------------------------------------------------------------------
% function to get the value of the D / T coordinates
function varargout = pushbutton11_Callback(h, eventdata, handles, varargin)
%% get coordinates
% Stub for Callback of the uicontrol handles.pushbutton11.

global d2 vd vr hd hr
%global d1 % currently not used
[H,V]=ginput(1);
if hr == 1
    H = 10^(H);
    set(handles.edit35,'string',H);
end
if hd == 1
    H = (10^(H))*1e-9;
    set(handles.edit35,'string',H);
end
if d2 == 1
    if vr == 1
        V = 10^(V);
        set(handles.edit36,'string',V);
    end
    if vd == 1
        V = (10^(V));
        V = V * 1e-9;
        set(handles.edit36,'string',V);
    end
end
        

% --------------------------------------------------------------------
% display of the horizontal value
function varargout = edit35_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit35.

% --------------------------------------------------------------------
% display of the vertical value
function varargout = edit36_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit36.

% --------------------------------------------------------------------
% Use 'ginput' to get new nnls values.
function varargout = pushbutton12_Callback(h, eventdata, handles, varargin)
%% get new borders for nnls
% Stub for Callback of the uicontrol handles.pushbutton12.
global d2

[Th,Tv]=ginput(2);
set(handles.edit13,'string',num2str(min(10.^Th)));
set(handles.edit14,'string',num2str(max(10.^Th)));
if d2 == 1
    set(handles.edit16,'string',num2str(min(10.^Tv)));
    set(handles.edit17,'string',num2str(max(10.^Tv)));
end



% --------------------------------------------------------------------
% Reset the nnls parameters to the last values stored.
function varargout = pushbutton13_Callback(h, eventdata, handles, varargin)
%% reset values for nnls
% Stub for Callback of the uicontrol handles.pushbutton13.
global d2 Store_Thmin Store_Thmax Store_hsteps Store_Tvmin Store_Tvmax Store_vsteps

set(handles.edit13,'string',num2str(Store_Thmin));
set(handles.edit14,'string',num2str(Store_Thmax));
set(handles.edit15,'string',num2str(Store_hsteps));

if d2 == 1
	set(handles.edit16,'string',num2str(Store_Tvmin));
	set(handles.edit17,'string',num2str(Store_Tvmax));
	set(handles.edit18,'string',num2str(Store_vsteps));
end



% --------------------------------------------------------------------
% Store the nnls values.
function varargout = pushbutton14_Callback(h, eventdata, handles, varargin)
%% store the values for nnls
% Stub for Callback of the uicontrol handles.pushbutton14.
global d2 Store_Thmin Store_Thmax Store_hsteps Store_Tvmin Store_Tvmax Store_vsteps

Store_Thmin = str2double(get(handles.edit13,'string'));
Store_Thmax = str2double(get(handles.edit14,'string'));
Store_hsteps = str2double(get(handles.edit15,'string'));

if d2 == 1
	Store_Tvmin = str2double(get(handles.edit16,'string'));
	Store_Tvmax = str2double(get(handles.edit17,'string'));
	Store_vsteps = str2double(get(handles.edit18,'string'));
end



% --------------------------------------------------------------------
% Choice of 'auto' or 'boolean' colour axis
function varargout = checkbox9_Callback(h, eventdata, handles, varargin)
%% change color axis
% Stub for Callback of the uicontrol handles.checkbox9.
axes(handles.axes2);
if get(handles.checkbox9,'value')==0
	caxis('auto');
else
	caxis([0 1]);
end



% --------------------------------------------------------------------
% Auto crop axes
function varargout = pushbutton15_Callback(h, eventdata, handles, varargin)
%% crop axis for nnls
% Stub for Callback of the uicontrol handles.pushbutton15.
global spectrum d2

axes(handles.axes2);
[m,n]=find(spectrum');
ax=axis;
Thmin = ax(1);
Thmax = ax(2);
% Thmin = log10(str2num(get(handles.edit13,'string')));
% Thmax = log10(str2num(get(handles.edit14,'string')));
hsteps = str2double(get(handles.edit15,'string'));

if d2 == 1
	Tvmin = ax(3);
	Tvmax = ax(4);
% 	Tvmin = log10(str2num(get(handles.edit16,'string')));
% 	Tvmax = log10(str2num(get(handles.edit17,'string')));
	vsteps = str2double(get(handles.edit18,'string'));
end

set(handles.edit13,'string',num2str(10^(Thmin+min(m-2)*(Thmax-Thmin)/(hsteps-1))));
set(handles.edit14,'string',num2str(10^(Thmin+max(m)*(Thmax-Thmin)/(hsteps-1))));

if d2 == 1
	set(handles.edit16,'string',num2str(10^(Tvmin+min(n-2)*(Tvmax-Tvmin)/(vsteps-1))));
	set(handles.edit17,'string',num2str(10^(Tvmin+max(n)*(Tvmax-Tvmin)/(vsteps-1))));
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
%% change the horizontal kernel
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

nr = get(hObject,'Value');
switch nr
    case 1 % Diffusion
        set(handles.edit37,'enable','off');
        set(handles.edit37,'string','exp(-h*D)');
    case 2 % T1-Relax
        set(handles.edit37,'enable','off');
        set(handles.edit37,'string','1-2*exp(-h/T)');
    case 3 % T1-Sat
        set(handles.edit37,'enable','off');
        set(handles.edit37,'string','1-exp(-h/T)');
    case 4 % T2-Relax
        set(handles.edit37,'enable','off');
        set(handles.edit37,'string','exp(-h/T)');
    case 5 % user-def
        set(handles.edit37,'enable','on');
end

CheckSyntaxAndAssign(handles,'h'); % will always be true, since formulae are set here
activate(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
%% change the vertical kernel
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

nr = get(hObject,'Value');
switch nr
    case 1 % Diffusion
        set(handles.edit38,'enable','off');
        set(handles.edit38,'string','exp(-v*D)');
    case 2 % T1-Relax
        set(handles.edit38,'enable','off');
        set(handles.edit38,'string','1-2*exp(-v/T)');
    case 3 % T1-Sat
        set(handles.edit38,'enable','off');
        set(handles.edit38,'string','1-exp(-v/T)');
    case 4 % T2-Relax
        set(handles.edit38,'enable','off');
        set(handles.edit38,'string','exp(-v/T)');
    case 5 % user-def
        set(handles.edit38,'enable','on');
end
CheckSyntaxAndAssign(handles,'v'); % will always be true, since formulae are set here
activate(handles);


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit37_Callback(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit37 as text
%        str2double(get(hObject,'String')) returns contents of edit37 as a double
CheckSyntaxAndAssign(handles,'h');


% --- Executes during object creation, after setting all properties.
function edit37_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit37 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit38_Callback(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit38 as text
%        str2double(get(hObject,'String')) returns contents of edit38 as a double
CheckSyntaxAndAssign(handles,'v');


% --- Executes during object creation, after setting all properties.
function edit38_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



