function varargout = projpreproc(varargin)

%% TODO Add these function
%%  Dark-flood correction
%%  Logarithmic Projection normalization
%%  Cosine Weighting
%%  Parker Weighing
%%  Ramp Filter

%% Ref: 
%%  1:Jeffrey, H. S., et al. (2013). Advances in 3D Image Reconstruction
%%      Image Processing in Radiation Therapy, CRC Press: 171-192.

% PROJPREPROC MATLAB code for projpreproc.fig
%      PROJPREPROC, by itself, creates a new PROJPREPROC or raises the existing
%      singleton*.
%
%      H = PROJPREPROC returns the handle to a new PROJPREPROC or the handle to
%      the existing singleton*.
%
%      PROJPREPROC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJPREPROC.M with the given input arguments.
%
%      PROJPREPROC('Property','Value',...) creates a new PROJPREPROC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before projpreproc_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to projpreproc_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help projpreproc

% Last Modified by GUIDE v2.5 12-Aug-2014 17:27:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @projpreproc_OpeningFcn, ...
                   'gui_OutputFcn',  @projpreproc_OutputFcn, ...
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


% --- Executes just before projpreproc is made visible.
function projpreproc_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to projpreproc (see VARARGIN)

% Choose default command line output for projpreproc
handles.output = hObject;

%Display logo
[pic,map] = imread('NSRL_Logo.gif');
SoftwareLogo = ind2rgb(pic,map);
axes(handles.logo);
image(SoftwareLogo);
axis off

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes projpreproc wait for user response (see UIRESUME)
% uiwait(handles.ProjGUI);


% --- Outputs from this function are returned to the command line.
function varargout = projpreproc_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loaddata.
function loaddata_Callback(hObject, eventdata, handles)
% hObject    handle to loaddata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[csv_filename, DefaultPathname] = uigetfile( {'*.csv','csv Files'}, 'Pick a csv file');
fid_Name = strcat(DefaultPathname,csv_filename);


fid = fopen(fid_Name,'r'); % Open CSV file for reading
C = textscan( fid,' %s %f %f %f %f %d','Delimiter', ',' );
    proj_param.proj_filenames = C{1};
    proj_param.theta = C{2};
    proj_param.u_off = C{3};
    proj_param.v_off = C{4};
    proj_param.I_0 = C{5};
    proj_param.weights = C{6};

    charName = fullfile(DefaultPathname,proj_param.proj_filenames{1});

    I1 = imread(charName);
    [nu,nv]=size(I1);
    N = size(proj_param.proj_filenames,1);
    proj_param.N_proj=N;
    proj_param.N_row=nu;
    proj_param.N_col=nv;
    
    mem = OSCaRProjGUIMemory(nu,nv,N);
    set(handles.mem_require,'String',mem);
    
    try
        proj_data.P = zeros(nu,nv,N);
    catch
        error('Out of Memory!');
    end

    for i=1:N
%         tmp = proj_param.proj_filenames(i);
%         FormStr=findstr(tmp{1},'.dcm');

        tmp = fullfile(DefaultPathname,proj_param.proj_filenames{i});
        FormStr=findstr(tmp,'.dcm');
        if isempty(FormStr)
            %This will handle all the file formats covered by imread.
            %I=imread(tmp{1});
            I = imread(tmp);
        else
            %The gray DICOM images will be stored in 2-D matrices.
            %I=dicomread(tmp{1});
            I = dicomread(tmp);
        end

        proj_data.P(:,:,i)=I;
    end

for k=1:proj_param.N_proj  %% normalization
     proj_data.P(:,:,k)=log(proj_param.I_0(k))-log(proj_data.P(:,:,k));
     A=proj_data.P(:,:,k);
     ind=find(isinf(A));
     A(ind)=log(proj_param.I_0(k));
     proj_data.P(:,:,k)=A;     
end

h=findobj('Tag','ProjGUI');
setappdata(h,'MyData',proj_data);
setappdata(h,'MyParam',proj_param);

global l1;
global l2;
global l3;

nu=proj_param.N_row;
nv=proj_param.N_col;
N = proj_param.N_proj; %size(proj_param.proj_filenames,1);
set(handles.projslider,'Min',1);
set(handles.projslider,'Max',N);

set(handles.projslider,'SliderStep',[5/N,5/N]);

val = get(handles.projslider,'Value');
if isnumeric(val) && length(val) == 1 && ...
        (val > N || val < 1)
    warndlg(sprintf('The value of the slider is out of range. It must be between 1 and %d.',N),'Out of Range');
else

    if isnumeric(val) && length(val) == 1 && ...
            val >= get(handles.projslider,'Min') && ...
            val <= get(handles.projslider,'Max')
        
        axes(handles.projImg);
        cla;
        try
            imagesc(squeeze(proj_data.P(:,:,val)));
        catch
            error('Pixel Type not right!');
         
        end

    
        colormap(gray(236));
        xlabel('u (pixels)');
        ylabel('v (pixels)');
        axis xy;
        axis equal;
        axis tight;
        colormap bone;
        colorbar;
        cross=line('XData',proj_param.u_off(val),'YData',proj_param.v_off(val),'Marker','x','Color','r', 'MarkerSize',18);
        
        
%         axes(handles.u_off);
%         cla;
%         plot(proj_param.u_off);
%         axis tight
    l1=line('XData',val,'YData',proj_param.u_off(val),'Marker','x','Color','r');
%         xlabel('N');
%         ylabel('u_{off} (pixels)');
%         
%         axes(handles.v_off);
%         cla;
%         plot(proj_param.v_off);
%         axis tight
    l2=line('XData',val,'YData',proj_param.v_off(val),'Marker','x','Color','r');
%         xlabel('N');
%         ylabel('v_{off} (pixels)');
%         
%         axes(handles.theta);
%         cla;
%         plot(proj_param.theta);
%         axis tight
    l3=line('XData',val,'YData',proj_param.theta(val),'Marker','x','Color','r');
%         xlabel('N');
%         ylabel('\theta (degrees)');
%         
%         set(handles.uoff,'String',num2str(proj_param.u_off(val)));
%         set(handles.voff,'String',num2str(proj_param.v_off(val)));
%         set(handles.projangle,'String',num2str(proj_param.theta(val)));
%         set(handles.weight,'String',num2str(proj_param.weights(val)));
%         set(handles.airnormal,'String',num2str(proj_param.I_0(val)));
    else
        set(handles.projno,'String',...
            ['You have entered an invalid entry ']);
    end
end

function memstr=OSCaRProjGUIMemory(nu,nv,N)
%Nargol Rezvani
%July 2008
%OSCaRMemory approximates the MEM_REQUIRE needed to OSCaRPreprocess.
%------------
%Notes:
% 1KB = 1024B.
% 1MB = 1048576B.
% 1GB = 1073741824B.
% 

memo = 8*(2*N*nu*nv);
if memo>1023 && memo< 1048576
        units='KB';
        memo=memo/1024;
else
    if memo> 1048575 && memo< 1073741824
        units=' MB';
        memo=memo/1048576;
    else
        if memo>1073741823
            units=' GB';
            memo=memo/1073741824;
        else
            units=' Bytes';
        end
    end        
end
memstr=strcat(num2str(2^nextpow2(memo)),units);


function edit_du_Callback(hObject, eventdata, handles)
% hObject    handle to edit_du (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_du as text
%        str2double(get(hObject,'String')) returns contents of edit_du as a double
h = findobj('Tag','ProjGUI');
proj_param=getappdata(h,'MyParam');

proj_param.du=str2num(get(hObject,'String'));
setappdata(h,'MyParam',proj_param);



% --- Executes during object creation, after setting all properties.
function edit_du_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_du (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dv_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dv as text
%        str2double(get(hObject,'String')) returns contents of edit_dv as a double

h = findobj('Tag','ProjGUI');
proj_param=getappdata(h,'MyParam');

proj_param.dv=str2num(get(hObject,'String'));
setappdata(h,'MyParam',proj_param);

% --- Executes during object creation, after setting all properties.
function edit_dv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SAD_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SAD as text
%        str2double(get(hObject,'String')) returns contents of edit_SAD as a double
h = findobj('Tag','ProjGUI');
proj_param=getappdata(h,'MyParam');

proj_param.SAD=str2num(get(hObject,'String'));
setappdata(h,'MyParam',proj_param);

% --- Executes during object creation, after setting all properties.
function edit_SAD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SAD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SDD_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SDD as text
%        str2double(get(hObject,'String')) returns contents of edit_SDD as a double
h = findobj('Tag','ProjGUI');
proj_param=getappdata(h,'MyParam');

proj_param.SDD=str2num(get(hObject,'String'));
setappdata(h,'MyParam',proj_param);


% --- Executes during object creation, after setting all properties.
function edit_SDD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SDD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on slider movement.
% function slider1_Callback(hObject, eventdata, handles)
% % hObject    handle to slider1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% 
% 
% % --- Executes during object creation, after setting all properties.
% function slider1_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to slider1 (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end


% --- Executes on slider movement.
function projslider_Callback(hObject, eventdata, handles)
% hObject    handle to projslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



h = findobj('Tag','ProjGUI');
proj_data = getappdata(h,'MyData');
proj_param = getappdata(h,'MyParam');

set(0,'RecursionLimit',proj_param.N_proj);

%global cross;
global l1;
global l2;
global l3;

val = round(get(hObject,'Value'));
set(handles.projno,'String',num2str(val));
set(hObject,'Value',val)
if ~exist('DefaultPathname','var')
    DefaultPathname = [];
end

% set(handles.uoff,'String',num2str(proj_param.u_off(val)));
% set(handles.voff,'String',num2str(proj_param.v_off(val)));
% set(handles.projangle,'String',num2str(proj_param.theta(val)));
% set(handles.weight,'String',num2str(proj_param.weights(val)));
% set(handles.airnormal,'String',num2str(proj_param.I_0(val)));

axes(handles.projImg);
imagesc(squeeze(proj_data.P(:,:,val)));
cross=line('XData',proj_param.u_off(val),'YData',proj_param.v_off(val),'Marker','x','Color','r', 'MarkerSize',18);

xlabel('u (pixels)');
ylabel('v (pixels)');
axis xy;
axis equal;
axis tight;
colormap bone;
colorbar;

set(l1,'XData',val,'YData',proj_param.u_off(val),'Marker','x','Color','r');
set(l2,'XData',val,'YData',proj_param.v_off(val),'Marker','x','Color','r');
set(l3,'XData',val,'YData',proj_param.theta(val),'Marker','x','Color','r');


% --- Executes during object creation, after setting all properties.
function projslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current position of the GUI from the handles structure
% to pass to the modal dialog.
pos_size = get(handles.ProjGUI,'Position');

% Call modaldlg with the argument 'Position'.
user_response = modaldlg('Title','Confirm Close');
switch user_response
case {'No'}
	% take no action
case 'Yes'
	% Prepare to close GUI application window
	%                  .
	%                  .
	%                  .
	delete(handles.ProjGUI)
end



% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h = findobj('Tag','ProjGUI');
proj_data=getappdata(h,'MyData');
proj_param=getappdata(h,'MyParam');



tag=1;
fnames=fieldnames(proj_param);
%the number of the fields is 15.

if size(fnames,1)==13
    
    for i=1:13
        if isempty(getfield(proj_param,fnames{i}))
            tag=0;
        end
    end
else
    tag=0;
end

if tag==0
     errordlg('Check the fields again!','Information Missing');
else
    proj_param.IAD=proj_param.SDD-proj_param.SAD;
    [mat_filename, DefaultPathname] = uiputfile( '*.mat', 'Save as');
 
    FileName = strcat(DefaultPathname,mat_filename);
    experiment.param=proj_param;
    experiment.data=proj_data;
    str = [ 'save ''' FileName ''' experiment' ];
    eval(str);
    h = findobj('Tag','ProjGUI');
    setappdata(h, 'MyData', proj_data);
    setappdata(h, 'MyParam', proj_param);
    msgbox('Now proceed to the reconstruct stage','Next Stage');
end

function mem_require_Callback(hObject, eventdata, handles)
% hObject    handle to mem_require (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mem_require as text
%        str2double(get(hObject,'String')) returns contents of mem_require as a double


% --- Executes during object creation, after setting all properties.
function mem_require_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mem_require (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function projno_Callback(hObject, eventdata, handles)
% hObject    handle to projno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of projno as text
%        str2double(get(hObject,'String')) returns contents of projno as a double

val = round(get(hObject,'Value'));
set(handles.projno,'String',num2str(val));
set(hObject,'Value',val)

h = findobj('Tag','OSCaRProjGUI');
proj_data=getappdata(h,'MyData');
proj_param=getappdata(h,'MyParam');

global l1;
global l2;
global l3;
%global l4;

val = str2double(get(hObject,'String'));

% Determine whether val is a number between the
% slider's Min and Max. If it is, set the slider Value.
if isnumeric(val) && length(val) == 1 && ...
        val >= get(handles.proj_slider,'Min') && ...
        val <= get(handles.proj_slider,'Max')
    set(handles.proj_slider,'Value',val);
%     set(handles.uoff,'String',num2str(proj_param.u_off(val)));
%     set(handles.voff,'String',num2str(proj_param.v_off(val)));
%     set(handles.projangle,'String',num2str(proj_param.theta(val)));
%     set(handles.weight,'String',num2str(proj_param.weights(val)));
%     set(handles.airnormal,'String',num2str(proj_param.I_0(val)));

    
    set(l1,'XData',val,'YData',proj_param.u_off(val),'Marker','x','Color','r');
    set(l2,'XData',val,'YData',proj_param.v_off(val),'Marker','x','Color','r');
    set(l3,'XData',val,'YData',proj_param.theta(val),'Marker','x','Color','r');
    if ~exist('DefaultPathname','var')
        DefaultPathname = [];
    end

    axes(handles.projimg);
    imagesc(squeeze(proj_data.P(:,:,val)));
    cross=line('XData',proj_param.u_off(val),'YData',proj_param.v_off(val),'Marker','x','Color','r', 'MarkerSize',18);
    xlabel('u (pixels)');
    ylabel('v (pixels)');
    axis xy;
    axis equal;
    axis tight;
    colormap bone;
    colorbar;
    
else
    % Increment the error count, and display it.
    set(hObject,'String',...
        ['You have entered an invalid entry ']);
    
end

% --- Executes during object creation, after setting all properties.
function projno_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projno (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function projangle_Callback(hObject, eventdata, handles)
% hObject    handle to projangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of projangle as text
%        str2double(get(hObject,'String')) returns contents of projangle as a double
h = findobj('Tag','ProjGUI');
proj_param=getappdata(h,'MyParam');
ind = get (handles.projslider,'Value');

proj_param.theta(1,ind) = str2num(get(hObject, 'String'));
setappdata(h,'MyParam',proj_param);

% --- Executes during object creation, after setting all properties.
function projangle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to projangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
