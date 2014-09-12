function varargout = PCBCTmain(varargin)
%% PCBTmain  main function
% PCBCTMAIN MATLAB code for PCBCTmain.fig
%      PCBCTMAIN, by itself, creates a new PCBCTMAIN or raises the existing
%      singleton*.
%
%      H = PCBCTMAIN returns the handle to a new PCBCTMAIN or the handle to
%      the existing singleton*.
%
%      PCBCTMAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCBCTMAIN.M with the given input arguments.
%
%      PCBCTMAIN('Property','Value',...) creates a new PCBCTMAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCBCTmain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCBCTmain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCBCTmain

% Last Modified by GUIDE v2.5 08-Aug-2014 13:56:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PCBCTmain_OpeningFcn, ...
                   'gui_OutputFcn',  @PCBCTmain_OutputFcn, ...
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


%% --- Executes just before PCBCTmain is made visible.
function PCBCTmain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCBCTmain (see VARARGIN)

% Choose default command line output for PCBCTmain
handles.output = hObject;

%Display logo
[pic,map] = imread('NSRL_Logo.gif');
SoftwareLogo = ind2rgb(pic,map);
axes(handles.axes1);
image(SoftwareLogo);
axis off
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PCBCTmain wait for user response (see UIRESUME)
% uiwait(handles.mainfig);


% --- Outputs from this function are returned to the command line.
function varargout = PCBCTmain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in reconstruction.
function reconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to reconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
recon;
    


% --- Executes on button press in preprocess.
function preprocess_Callback(hObject, eventdata, handles)
% hObject    handle to preprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
projpreproc;


% --- Executes on button press in quit.
function quit_Callback(hObject, eventdata, handles)
% hObject    handle to quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% % Construct a quest dialogue
% choice = questdlg('Quit?', 'Quit menu', ...
% 	'Yes','No','No');
% % Handle response
% switch choice
%     case 'Yes'
%         delete(gcbf);
%     case 'No'
%         
% end

% Get the current position of the GUI from the handles structure
% to pass to the modal dialog.
pos_size = get(handles.mainfig,'Position');

% Call modaldlg with the argument 'Position'.
user_response = modaldlg('Title','Main window');
switch user_response
case {'No'}
	% take no action
case 'Yes'
	% Prepare to close GUI application window
	%                  .
	%                  .
	%                  .
	delete(handles.mainfig)
end
