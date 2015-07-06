function varargout = twodlaplaceinverse(varargin)
% TwoDLaplaceInverse
% Last modified 16-May-2003
% Developers: Sophie Godefroy and Brett Ryland


% TWODLAPLACEINVERSE Application M-file for twodlaplaceinverse.fig
%    FIG = TWODLAPLACEINVERSE launch twodlaplaceinverse GUI.
%    TWODLAPLACEINVERSE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 04-Sep-2002 17:38:079

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


% --------------------------------------------------------------------
% choice for a 1D experiment
function varargout = radiobutton7_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton7.

global d1 d2

d1 = 1;
set(handles.radiobutton7,'value',1);
set(handles.radiobutton8,'value',0);
activate(handles);

% --------------------------------------------------------------------
% choice for a 2D experiment
function varargout = radiobutton8_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton8.

global d1 d2

d2 = 1;
set(handles.radiobutton8,'value',1);
set(handles.radiobutton7,'value',0);
activate(handles);

% --------------------------------------------------------------------
% enable/disable various buttons depending on what is selected
function activate(handles)
%

global unlocked

if get(handles.radiobutton8,'value') == 1 % 2D experiment
	set(handles.radiobutton5,'enable','on');
	set(handles.radiobutton6,'enable','on');
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
	set(handles.radiobutton5,'enable','off');
	set(handles.radiobutton6,'enable','off');
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
if (get(handles.radiobutton4,'value') == 1) | ((get(handles.radiobutton5,'value') == 1) & (get(handles.radiobutton8,'value') == 1)) % diffusion expt
	set(handles.checkbox4,'enable','on');
	set(handles.checkbox5,'enable','on');
	if (get(handles.checkbox4,'value') == 1) | (get(handles.checkbox5,'value') == 1) % generate q^2
		if(get(handles.radiobutton4,'value') == 1)
			set(handles.edit2,'enable','off'); % time 1
			set(handles.radiobutton1,'enable','off');
			set(handles.edit9,'enable','off');
		end
		if((get(handles.radiobutton5,'value') == 1) & (get(handles.radiobutton8,'value') == 1))
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
%

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
% choice for horizontal relaxation experiment
function varargout = radiobutton3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton3.

global hr hd vd

hr = 1;
hd = 0;
set(handles.radiobutton3, 'value', 1);
set(handles.radiobutton4, 'value', 0);
activate(handles);

% --------------------------------------------------------------------
% choice for horizontal diffusion experiment
function varargout = radiobutton4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton4.

global hr hd

hd = 1;
hr = 0;
set(handles.radiobutton4, 'value', 1);
set(handles.radiobutton3, 'value', 0);
activate(handles);

% --------------------------------------------------------------------
% choice for vertical relaxation experiment
function varargout = radiobutton6_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton6.

global vr vd hd

vr = 1;
vd = 0;
set(handles.radiobutton6, 'value', 1);
set(handles.radiobutton5, 'value', 0);
activate(handles);

% --------------------------------------------------------------------
% choice for vertical diffusion experiment
function varargout = radiobutton5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton5.

global vr vd

vd = 1;
vr = 0;
set(handles.radiobutton5, 'value', 1);
set(handles.radiobutton6, 'value', 0);
activate(handles);

% --------------------------------------------------------------------
% button to load data file
function varargout = pushbutton1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton1.

global data filedata chemin vr hr nfig fichier OKh OKv videos videof dis trans

[fichier,chemin] = uigetfile('*.*');
set(handles.edit1,'string',fichier);
pwd2 = pwd;
cd(chemin);
data = load(fichier);
filedata=data;
cd(pwd2);

vr = get(handles.radiobutton6,'value');
hr = get(handles.radiobutton3,'value');

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
function varargout = checkbox6_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.checkbox6.

global T1T2dis

T1T2dis = get(h,'value');


% --------------------------------------------------------------------
% draw data button
function varargout = pushbutton4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton4.

global T1T2dis dis cmth cmtv d1 d2 hr hd vr vd timea timeb data filetimea filetimeb filedata Gmax delta DELTA vsize hsize cmth cmtv nfig OKh OKv gradvar delvar Store_Thmin Store_Thmax Store_hsteps Store_Tvmin Store_Tvmax Store_vsteps


hr = get(handles.radiobutton3,'value');
vr = get(handles.radiobutton6,'value');
d1 = get(handles.radiobutton7,'value');
d2 = get(handles.radiobutton8,'value');
T1T2dis = get(handles.checkbox6,'value');
cmth = get(handles.radiobutton1, 'value');
vd = get(handles.radiobutton5,'value');
hd = get(handles.radiobutton4,'value');
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

if d2 == 1
    [sdb,sda] = size(data);
end
if d1 == 1
    [sda,sdb] = size(data);
end
set(handles.edit7,'string',sda);
set(handles.edit8,'string',sdb);

if (hr == 1) | ((hd == 1) & (gradvar == 0) & (delvar == 0))
    [sta,nimportequoi] = size(timea);
end
if (d2 == 1) & ((vr == 1) | ((vd == 1) & (gradvar == 0) & (delvar == 0)))
    [stb,nimportequoi] = size(timeb);
end
if (d1 == 1) & (T1T2dis == 1) & (dis == 0)
    dis = 1;
    data = (max(abs(data(sda,1)),abs(data(1,1)))-data)/2;
end

if cmth == 1
    tauvh = get(handles.edit9,'string');
    tauvh = str2num(tauvh);
    timea = filetimea .* tauvh;
end
cmtv = get(handles.radiobutton2, 'value');
if cmtv == 1
    tauvv = get(handles.edit10,'string');
    tauvv = str2num(tauvv);
    timeb = filetimeb .* tauvv;
end


% Check the validity of the data
if (hr == 1) | ((gradvar == 0) & (delvar == 0))
	if sta ~= sda
		warndlg('error: chose the good data and time set', '!!Sorry!!');
		return
	end
end
if (d2 == 1) & ((vr == 1) | ((gradvar == 0) & (delvar == 0)))
    if stb ~= sdb
        warndlg('error: chose the good data and time set', '!!Sorry!!');
		return
    end
end


if (hd == 1) | (vd == 1)
    Gmax = get(handles.edit4,'string');
    Gmax = str2num(Gmax);
    delta = get(handles.edit5,'string');
    delta = str2num(delta);
    DELTA = get(handles.edit6,'string');
    DELTA = str2num(DELTA);
end

gradvar = get(handles.checkbox4,'value');
delvar = get(handles.checkbox5,'value');
firstdel = get(handles.edit33,'string');
firstdel = str2num(firstdel);
stepdel = get(handles.edit34,'string');
stepdel = str2num(stepdel);

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
if vd == 1 & d2 == 1
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
% Stub for Callback of the uicontrol handles.edit9.

global tauvh

tauvh = get(h,'string');
tauvh = str2num(tauvh);


% --------------------------------------------------------------------
% tau value for the time 2
function varargout = edit10_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit10.

global tauvv

tauvv = get(h,'string');
tauvv = str2num(tauvv);

% --------------------------------------------------------------------
% gradient experiment by the variation of the gradients
function varargout = checkbox4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.checkbox4.

global gradvar delvar hd vd

gradvar = get(h,'value');
delvar = 0;
set(handles.checkbox5,'value',0);
activate(handles);

% --------------------------------------------------------------------
% gradient experiment by varying the duration of delta
function varargout = checkbox5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.checkbox5.

global delvar gradvar hd vd

delvar = get(h,'value');
gradvar = 0;
set(handles.checkbox4,'value',0);
activate(handles);

% --------------------------------------------------------------------
% value of G max
function varargout = edit4_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit4.

global Gmax

Gmax = get(h,'string');
Gmax = str2num(Gmax);


% --------------------------------------------------------------------
% value of delta (ms)
function varargout = edit5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit5.

global delta

delta = get(h,'string');
delta = str2num(delta);

% --------------------------------------------------------------------
% value of DELTA (ms)
function varargout = edit6_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit6.

global DELTA

DELTA = get(h,'string');
DELTA = str2num(DELTA);

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
% Stub for Callback of the uicontrol handles.checkbox1.

global nfig

nfig = get(h,'value');


% --------------------------------------------------------------------
% Function to remove rows or columns from the 2D data set
function varargout = pushbutton9_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton9.

global timea timeb data nfig d1 d2

if d2 == 1
    nbtop = get(handles.edit29,'string');
    nbtop = str2num(nbtop);
    nbbottom = get(handles.edit30,'string');
    nbbottom = str2num(nbbottom);
end
nbrigth = get(handles.edit31,'string');
nbrigth = str2num(nbrigth);
nbleft = get(handles.edit32,'string');
nbleft = str2num(nbleft);

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
    if nbrigth ~= 0
        data(:,sizea+1-nbrigth:sizea)=[];
        timea(sizea+1-nbrigth:sizea)=[];
    end
end

if d1 == 1
    if nbleft ~= 0
        data(1:nbleft)=[];
        timea(1:nbleft)=[];
    end
    [sizea,sizeb] = size(data);
    if nbrigth ~= 0
        data(sizea+1-nbrigth:sizea)=[];
        timea(sizea+1-nbrigth:sizea)=[];
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
% Stub for Callback of the uicontrol handles.pushbutton10.

global timea timeb data filetimea filetimeb filedata nfig d1 d2

hd = get(handles.radiobutton4,'value');
hr = get(handles.radiobutton3,'value');
vd = get(handles.radiobutton5,'value');
vr = get(handles.radiobutton6,'value');
gradvar = get(handles.checkbox4,'value');
delvar = get(handles.checkbox5,'value');
Gmax = str2num(get(handles.edit4,'string'));
delta = str2num(get(handles.edit5,'string'));
DELTA = str2num(get(handles.edit6,'string'));

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

if (hr == 1) | ((hd == 1) & (gradvar == 0) & (delvar == 0))
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
	if (vr == 1) | ((vd == 1) & (gradvar == 0) & (delvar == 0))
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
% Stub for Callback of the uicontrol handles.pushbutton5.

global data timea timeb nfig chemin fichier fi1 fi2 spectrum spectrumi tauv tauh XI YI T1st T2st hr hd vr vd taulv taulh sta stb d1 d2

% "a" = "h" and "b" = "v"
% acquisition of the parameters of the nnls smoothing process
alpha = get(handles.edit11,'string');
alpha = str2num(alpha);
beta = get(handles.edit12,'string');
beta = str2num(beta);
Thmin = get(handles.edit13,'string');
Thmin = str2num(Thmin);
Thmax = get(handles.edit14,'string');
Thmax = str2num(Thmax);
Thmm = [Thmin,Thmax];
stepsh = get(handles.edit15,'string');
stepsh = str2num(stepsh);
if d2 == 1
    Tvmin = get(handles.edit16,'string');
    Tvmin = str2num(Tvmin);
    Tvmax = get(handles.edit17,'string');
    Tvmax = str2num(Tvmax);
    Tvmm = [Tvmin, Tvmax];
    stepsv = get(handles.edit18,'string');
    stepsv = str2num(stepsv);
end
orient = get(handles.edit21,'string');

if hd == 1
    a = Thmin;
    Thmin = 1/Thmax;
    Thmax = 1/a;
    Thmm = [Thmin,Thmax];
end
if vd == 1 & d2 == 1
    a = Tvmin;
    Tvmin = 1/Tvmax;
    Tvmax = 1/a;
    Tvmm = [Tvmin,Tvmax];
end

tic;

if d2 == 1
    [spectrum,tauh,tauv,chisq,compte]=upnnlsmooth3Dsvdfin(data,timea',timeb',Thmm,stepsh,Tvmm,stepsv,alpha,beta,orient);
end

if d1 == 1
    [spectrum,tauh,chisq,compte]=upnnlsmooth1D(data,timea,Thmin,Thmax,alpha,beta,stepsh);
end
if hd == 1 
    spectrum = flipdim(spectrum,2);
    tauh = 1./tauh;
    tauh = flipdim(tauh,2);
end
if vd == 1 & d2 == 1
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

if (hr == 1) | ((hd == 1) & (gradvar == 0) & (delvar == 0))
    savefile = [fi1 '.out'];
    tauh = tauh';
    save(savefile,'tauh','-ASCII');
elseif hd == 1
    savefile = [fichier 'grad' '.out'];
    tauh = tauh';
    save (savefile,'tauh','-ASCII');
end
if d2 == 1
	if (vr == 1) | ((hd == 1) & (gradvar == 0) & (delvar == 0))
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
% Stub for Callback of the uicontrol handles.pushbutton6.

global spectrum taulv taulh nfig sta stb vd vr hd hr fichier videos videof chemin d1 d2

vd1 = vd;
vr1 = vr;
hd1 = hd;
hr1 = hr;
spectrum1 = spectrum;
trans = get(handles.checkbox8,'value')
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


if nfig == 1 & d2 == 1
    fig=figure;
    axes('FontSize',20);
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
        xlabel('Log10(T) (s)','FontSize',20);
    elseif hd == 1
        xlabel('Log10(D) (10^-9 m^2/s)','FontSize',20);
    end
    if vr == 1 
        ylabel('Log10(T) (s)','FontSize',20);
    elseif vd == 1
        ylabel('Log10(D) (10^-9 m^2/s)','FontSize',20);
    end
    if hd == 1 & vd == 1
        title('D-D correlation','FontSize',20);
    elseif hr == 1 & vr == 1
        title('T-T correlation','FontSize',20);
    else
        title('D-T correlation','FontSize',20);
    end
end

if nfig == 1 & d1 == 1
    fig = figure;
    axes('FontSize',20);
    plot(taulh,spectrum,'k')
    xlim([taulh(1),taulh(sta(2))]);
    if hr == 1
        xlabel('Log10(T) (s)','FontSize',20);
    elseif hd == 1
        xlabel('Log10(D) (10^-9 m^2/s)','FontSize',20);
    end
end    

pwd1 = pwd;
cd(chemin);
savefile = [fichier 'm.fig'];
saveas(gcf,savefile);
savefile = [fichier 'm.jpg'];
saveas(gcf,savefile);


if nfig == 1 & videos == 1
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

if nfig == 1 & videof == 1
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
% Stub for Callback of the uicontrol handles.pushbutton8.

global data timea timeb nfig chemin fichier fi1 fi2 spectrum spectrumi tauv tauh XI YI T1st T2st hd vd d1 d2

alphaf = get(handles.edit24,'string');
alphaf = str2num(alphaf);
alphal = get(handles.edit25,'string');
alphal = str2num(alphal);
alphanb = get(handles.edit26,'string');
alphanb = str2num(alphanb);

beta = get(handles.edit12,'string');
beta = str2num(beta);
Thmin = get(handles.edit13,'string');
Thmin = str2num(Thmin);
Thmax = get(handles.edit14,'string');
Thmax = str2num(Thmax);
Thmm = [Thmin,Thmax];
stepsh = get(handles.edit15,'string');
stepsh = str2num(stepsh);
if d2 == 1
    Tvmin = get(handles.edit16,'string');
    Tvmin = str2num(Tvmin);
    Tvmax = get(handles.edit17,'string');
    Tvmax = str2num(Tvmax);
    Tvmm = [Tvmin, Tvmax];
    stepsv = get(handles.edit18,'string');
    stepsv = str2num(stepsv);
end
orient = get(handles.edit21,'string');

if hd == 1
    a = Thmin;
    Thmin = 1/Thmax;
    Thmax = 1/a;
    Thmm = [Thmin,Thmax];
end
if vd == 1 & d2 == 1
    a = Tvmin;
    Tvmin = 1/Tvmax;
    Tvmax = 1/a;
    Tvmm = [Tvmin,Tvmax];
end


alphapas = (log10(alphal)-log10(alphaf))/(alphanb-1);
alphav = exp((log10(alphaf) : alphapas : log10(alphal))*log(10));

for na = 1 : alphanb
    set(handles.edit27,'string',na);
    tic;
    if d2 == 1
        [spectrum,tauh,tauv,chisq,compte]=upnnlsmooth3Dsvdfin(data,timea',timeb',Thmm,stepsh,Tvmm,stepsv,alphav(na),beta,orient);
    elseif d1 == 1
        [spectrum,tauh,chisq,compte]=upnnlsmooth1D(data,timea,Thmin,Thmax,alphav(na),beta,stepsh);
    end
    chiv(na) = chisq;
    if hd == 1
        spectrum = flipdim(spectrum,2);
        tauh = 1./tauh;
        tauh = flipdim(tauh,2);
    end
    if vd == 1 & d2 == 1
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
% Stub for Callback of the uicontrol handles.checkbox2.

global videos

videos = get(h,'value');

% --------------------------------------------------------------------
% Video front
function varargout = checkbox3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.checkbox3.

global videof

videof = get(h,'value');

% --------------------------------------------------------------------
% function to get the value of the D / T coordinates
function varargout = pushbutton11_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton11.

global d1 d2 vd vr hd hr

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
% Stub for Callback of the uicontrol handles.pushbutton14.
global d2 Store_Thmin Store_Thmax Store_hsteps Store_Tvmin Store_Tvmax Store_vsteps

Store_Thmin = str2num(get(handles.edit13,'string'));
Store_Thmax = str2num(get(handles.edit14,'string'));
Store_hsteps = str2num(get(handles.edit15,'string'));

if d2 == 1
	Store_Tvmin = str2num(get(handles.edit16,'string'));
	Store_Tvmax = str2num(get(handles.edit17,'string'));
	Store_vsteps = str2num(get(handles.edit18,'string'));
end



% --------------------------------------------------------------------
% Choice of 'auto' or 'boolean' colour axis
function varargout = checkbox9_Callback(h, eventdata, handles, varargin)
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
% Stub for Callback of the uicontrol handles.pushbutton15.
global spectrum d2

axes(handles.axes2);
[m,n]=find(spectrum');
ax=axis;
Thmin = ax(1);
Thmax = ax(2);
% Thmin = log10(str2num(get(handles.edit13,'string')));
% Thmax = log10(str2num(get(handles.edit14,'string')));
hsteps = str2num(get(handles.edit15,'string'));

if d2 == 1
	Tvmin = ax(3);
	Tvmax = ax(4);
% 	Tvmin = log10(str2num(get(handles.edit16,'string')));
% 	Tvmax = log10(str2num(get(handles.edit17,'string')));
	vsteps = str2num(get(handles.edit18,'string'));
end

set(handles.edit13,'string',num2str(10^(Thmin+min(m-2)*(Thmax-Thmin)/(hsteps-1))));
set(handles.edit14,'string',num2str(10^(Thmin+max(m)*(Thmax-Thmin)/(hsteps-1))));

if d2 == 1
	set(handles.edit16,'string',num2str(10^(Tvmin+min(n-2)*(Tvmax-Tvmin)/(vsteps-1))));
	set(handles.edit17,'string',num2str(10^(Tvmin+max(n)*(Tvmax-Tvmin)/(vsteps-1))));
end
