function varargout = recon(varargin)
% RECON MATLAB code for recon.fig
%      RECON, by itself, creates a new RECON or raises the existing
%      singleton*.
%
%      H = RECON returns the handle to a new RECON or the handle to
%      the existing singleton*.
%
%      RECON('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECON.M with the given input arguments.
%
%      RECON('Property','Value',...) creates a new RECON or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before recon_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to recon_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help recon

% Last Modified by GUIDE v2.5 14-Aug-2014 16:52:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @recon_OpeningFcn, ...
                   'gui_OutputFcn',  @recon_OutputFcn, ...
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


% --- Executes just before recon is made visible.
function recon_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to recon (see VARARGIN)

% Choose default command line output for recon
handles.output = hObject;

%Display logo
[pic,map] = imread('NSRL_Logo.gif');
SoftwareLogo = ind2rgb(pic,map);
axes(handles.logo);
image(SoftwareLogo);
axis off

% Update handles structure
guidata(hObject, handles);

if ~exist('DefaultPathname','var')
    DefaultPathname = [];
end

[mat_filename, DefaultPathname] = uigetfile( {'*.mat','MAT Files'}, 'Pick a MAT-file');
name = strcat(DefaultPathname,mat_filename);
load(name);

proj_param=experiment.param;
proj_data=experiment.data;

%Initialize the first figure.
axes(handles.axes2);
ind_theta = find(abs(diff(proj_param.theta))>0.1);
[t,ind2] = unique(proj_param.theta(ind_theta));
proj_param.k_RL = max(1,min(round(interp1(t,ind_theta(ind2),0,'nearest','extrap')),proj_param.N_proj));
proj_param.k_AP = max(1,min(round(interp1(t,ind_theta(ind2),90,'nearest','extrap')),proj_param.N_proj));

proj_data.log_P0=proj_data.P(:,:,proj_param.k_RL);
proj_data.log_P1=proj_data.P(:,:,proj_param.k_AP);

proj_param.u_off_ML = proj_param.u_off;
proj_param.v_off_ML = proj_param.v_off;
proj_param.MF = 1 + (proj_param.IAD)./(proj_param.SAD);

proj_param.u_aper = round(proj_param.u_off_ML(proj_param.k_RL))+[-20:20];
proj_param.v_aper = round(proj_param.v_off_ML(proj_param.k_RL))+[-20:20];
m = mean2(proj_data.log_P0(proj_param.v_aper,proj_param.u_aper));
s = std2(proj_data.log_P0(proj_param.v_aper,proj_param.u_aper));


imagesc(proj_data.log_P0,m+s*[-3,3])
axis xy
axis equal;
axis tight;

% Handles for the borders of reconstruction volume
proj_param.l1 = line([proj_param.u_off_ML(proj_param.k_RL)] ,[proj_param.v_off_ML(proj_param.k_RL)]);
proj_param.l2 = line([proj_param.u_off_ML(proj_param.k_RL)] ,[proj_param.v_off_ML(proj_param.k_RL)]);
proj_param.l3 = line([proj_param.u_off_ML(proj_param.k_RL)] ,[proj_param.v_off_ML(proj_param.k_RL)]);
proj_param.l4 = line([proj_param.u_off_ML(proj_param.k_RL)] ,[proj_param.v_off_ML(proj_param.k_RL)]);

colormap(gray(236))
xlabel('u [pixels]','Color',[0 0 1])
ylabel('v [pixels]','Color',[0 0 1])
title(['Projection #1 - ',num2str(proj_param.theta(1)),' degrees'])

%Handle for the cross on the image.
cross1 = line( [proj_param.u_off_ML(proj_param.k_RL)] ,[proj_param.v_off_ML(proj_param.k_RL)], 'Color','r', ...
      'MarkerSize',18, 'Marker', 'x') ;
set(gca,'XColor',[0 0 1],'YColor',[0 0 1])

%Initializeing the second figure.
axes(handles.axes3);
proj_param.u_aper = round(proj_param.u_off_ML(proj_param.k_AP))+[-20:20];
proj_param.v_aper = round(proj_param.v_off_ML(proj_param.k_AP))+[-20:20];

m = mean2(proj_data.log_P1(proj_param.v_aper,proj_param.u_aper));
s = std2(proj_data.log_P1(proj_param.v_aper,proj_param.u_aper));

imagesc(proj_data.log_P1,m+s*[-3,3])
axis xy
axis equal;
axis tight;

% Handles for the borders of reconstruction volume
proj_param.l5 = line([proj_param.u_off_ML(proj_param.k_AP)], [proj_param.v_off_ML(proj_param.k_AP)]);
proj_param.l6 = line([proj_param.u_off_ML(proj_param.k_AP)], [proj_param.v_off_ML(proj_param.k_AP)]);
proj_param.l7 = line([proj_param.u_off_ML(proj_param.k_AP)], [proj_param.v_off_ML(proj_param.k_AP)]);
proj_param.l8 = line([proj_param.u_off_ML(proj_param.k_AP)], [proj_param.v_off_ML(proj_param.k_AP)]);
xlabel('u [pixels]','Color',[0 0 1])
ylabel('v [pixels]','Color',[0 0 1])
title(['Projection #',num2str(proj_param.k_AP),' - ',num2str(proj_param.theta(proj_param.k_AP)),' degrees'])
cross2 = line( [proj_param.u_off_ML(proj_param.k_AP)], [proj_param.v_off_ML(proj_param.k_AP)], ...
         'Color', 'r','MarkerSize', 18, 'Marker', 'x' );
set(gca,'XColor',[0 0 1],'YColor',[0 0 1])

proj_param.dx = round(100*proj_param.du/mean(proj_param.MF))/100;
proj_param.dy = round(100*proj_param.du/mean(proj_param.MF))/100;
proj_param.dz = round(100*proj_param.dv/mean(proj_param.MF))/100;


set(handles.dx,'String',num2str(proj_param.dx));
set(handles.dy,'String',num2str(proj_param.dy));
set(handles.dz,'String',num2str(proj_param.dz));


%Set the Default Value of Filter,"No Filter"
%It means we are not filtering at all.

filterName='No Filter';
du=proj_param.du;
nu=proj_param.N_col;
nv=proj_param.N_row;

%d = str2num(get(handles.edit9,'String'));
d=1.0;%frequency
Phihat = OSCaRFilter( filterName, nu, du, d );
proj_param.Phihat=Phihat;
H=nv*proj_param.dv;
L=nu*proj_param.du;
%Limits of the reconstrcution domain.
proj_param.XMIN = -proj_param.SAD/proj_param.SDD*L/2;
proj_param.XMAX = proj_param.SAD/proj_param.SDD*L/2;
proj_param.YMIN = -proj_param.SAD/proj_param.SDD*L/2;
proj_param.YMAX = proj_param.SAD/proj_param.SDD*L/2;
proj_param.ZMIN = -proj_param.SAD/proj_param.SDD*H/2;
proj_param.ZMAX = proj_param.SAD/proj_param.SDD*H/2;

%Set the default values for xMin xMax, ...
set(handles.x_min,'String',num2str(round(proj_param.XMIN/2)));
set(handles.x_max,'String',num2str(round(proj_param.XMAX/2)));
set(handles.y_min,'String',num2str(round(proj_param.YMIN/2)));
set(handles.y_max,'String',num2str(round(proj_param.YMAX/2)));
set(handles.z_min,'String',num2str(0));
set(handles.z_max,'String',num2str(0));
setappdata(hObject,'MyParam',proj_param);
setappdata(hObject,'MyData',proj_data);

update_projGUI(handles);

function update_projGUI( handles )
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');
proj_param.dx = str2num(get(handles.dx,'String'));
proj_param.dy = str2num(get(handles.dy,'String'));
proj_param.dz = str2num(get(handles.dz,'String'));


%Get the grid

proj_param.xmin = str2num(get(handles.x_min, 'String'));
proj_param.x0 = proj_param.xmin/proj_param.dx+round(proj_param.u_off_ML(proj_param.k_AP));

proj_param.xmax = str2num(get(handles.x_max, 'String'));
proj_param.x1 = proj_param.xmax/proj_param.dx+round(proj_param.u_off_ML(proj_param.k_AP));

proj_param.ymin = str2num(get(handles.y_min, 'String'));
proj_param.y0 = proj_param.ymin/proj_param.dy+round(proj_param.u_off_ML(proj_param.k_RL));

proj_param.ymax = str2num(get(handles.y_max, 'String'));
proj_param.y1 = proj_param.ymax/proj_param.dy+round(proj_param.u_off_ML(proj_param.k_RL));

proj_param.zmin = str2num(get(handles.z_min, 'String'));
proj_param.z0 = proj_param.zmin/proj_param.dz+round(proj_param.v_off_ML(proj_param.k_RL));

proj_param.zmax = str2num(get(handles.z_max, 'String'));
proj_param.z1 = proj_param.zmax/proj_param.dz+round(proj_param.v_off_ML(proj_param.k_RL));

proj_param.x = -flipud(proj_param.dx*( [proj_param.x0:proj_param.x1]'-round(proj_param.u_off_ML(proj_param.k_AP))));         % Account for A-P view since +ve u is -ve x
proj_param.y = proj_param.dy*( [proj_param.y0:proj_param.y1]'-round(proj_param.u_off_ML(proj_param.k_RL)));
proj_param.z = proj_param.dz*( [proj_param.z0:proj_param.z1]'-round(proj_param.v_off_ML(proj_param.k_RL)));

if (str2num(get(handles.x_min,'String')) < proj_param.XMIN) ||...
        (str2num(get(handles.x_max,'String'))>proj_param.XMAX) ||...
        (str2num(get(handles.y_min,'String')) <proj_param.YMIN) ||...
        (str2num(get(handles.y_max,'String'))>proj_param.YMAX) || ...
        (str2num(get(handles.z_min,'String'))<proj_param.ZMIN) ||...
        (str2num(get(handles.z_max,'String'))>proj_param.ZMAX)
    errordlg('Values out of Range.','Wrong Range');

else if (str2num(get(handles.x_min,'String')) > str2num(get(handles.x_max,'String'))) || ...
            (str2num(get(handles.y_min,'String')) > str2num(get(handles.y_max,'String'))) || ...
            (str2num(get(handles.z_min,'String')) > str2num(get(handles.z_max,'String')))

        errordlg('Min value must be smaller than Max!','Wrong Range');
    else


        proj_param.x_size = length(proj_param.x);
        proj_param.y_size = length(proj_param.y);
        proj_param.z_size = length(proj_param.z);
        axes(handles.axes2);
        axis xy
        %This part plots the line on the right image.
        % Add two more handles to make a box...
        set( proj_param.l1, 'XData', [proj_param.y0;proj_param.y1], 'YData', [proj_param.z0;proj_param.z0], 'Color', 'c');
        set( proj_param.l2, 'XData', [proj_param.y0;proj_param.y1], 'YData', [proj_param.z1;proj_param.z1], 'Color', 'c');
        set( proj_param.l3, 'XData', [proj_param.y0;proj_param.y0], 'YData', [proj_param.z0;proj_param.z1], 'Color', 'c');
        set( proj_param.l4, 'XData', [proj_param.y1;proj_param.y1], 'YData', [proj_param.z0;proj_param.z1], 'Color', 'c');

        axes(handles.axes3);
        axis xy
        set( proj_param.l5, 'XData', [proj_param.x0;proj_param.x1], 'YData', [proj_param.z0;proj_param.z0],'Color','c');
        set( proj_param.l6, 'XData', [proj_param.x0;proj_param.x1], 'YData', [proj_param.z1;proj_param.z1],'Color','c');
        set( proj_param.l7, 'XData', [proj_param.x0;proj_param.x0], 'YData', [proj_param.z0;proj_param.z1],'Color','c');
        set( proj_param.l8, 'XData', [proj_param.x1;proj_param.x1], 'YData', [proj_param.z0;proj_param.z1],'Color','c') ;

        set(handles.nx,'String',num2str(proj_param.x_size));
        set(handles.ny,'String',num2str(proj_param.y_size));
        set(handles.nz,'String',num2str(proj_param.z_size));
    end
end

setappdata(h,'MyParam',proj_param);
mem = GUIMemory(length(proj_param.y),length(proj_param.x),length(proj_param.z),proj_param.N_row,proj_param.N_col);
set(handles.mem_required,'String',mem);
setappdata(h,'MyParam',proj_param);
set( handles.axes2, 'ButtonDownFcn', 'OSCaRReconstruct( ''axes2_ButtonDownFcn'', handles.axes2,[],handles )' )

% UIWAIT makes recon wait for user response (see UIRESUME)
% uiwait(handles.projGUI);

function memstr=GUIMemory(nx,ny,nz,nu,nv)
%Nargol Rezvani
%July 4th 2008
%OSCaRMemory approximates the RAM needed to OSCaRGUI.
%------------
%Notes:
% 1KB = 1024B.
% 1MB = 1048576B.
% 1GB = 1073741824B.
% 

memo = 8*(4*nx*ny*nz + 2*nu*nv);
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

function axes2_ButtonDownFcn( hObject, eventdata, handles )

h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');

currpt = get(hObject, 'CurrentPoint');
disp(currpt)
y = currpt(1,1);
z = currpt(1,2);
% coordinates in terms of pixels
ymin = proj_param.y0;
ymax = proj_param.y1;
zmin = proj_param.z0;
zmax = proj_param.z1;

disp( [y z ] )
disp( [ymin, ymax, zmin, zmax])
% Reset values of bounding box of FOV
if abs(y-ymin) <= abs(y-ymax)
    proj_param.y0 = y;
else
    proj_param.y1 = y;
end
if abs(z-zmin) <= abs(z-zmax)
    proj_param.z0 = z;
else
    proj_param.z1 = z;
end

% Update text boxes
yminstr=num2str(proj_param.dy*( proj_param.y0-round(proj_param.u_off_ML(proj_param.k_RL))));
set(handles.y_min, 'String', yminstr);

ymaxstr=num2str(proj_param.dy*( proj_param.y1-round(proj_param.u_off_ML(proj_param.k_RL))));
set(handles.y_max, 'String', ymaxstr);

zminstr=num2str(proj_param.dz*( proj_param.z0-round(proj_param.v_off_ML(proj_param.k_RL))));
set(handles.z_min, 'String', zminstr);

zmaxstr=num2str(proj_param.dz*( proj_param.z1-round(proj_param.v_off_ML(proj_param.k_RL))));
set(handles.z_max, 'String', zmaxstr);

% Update data and GUI appearance to reflect changes
setappdata(h,'MyParam',proj_param);
update_projGUI(handles);

% --- Outputs from this function are returned to the command line.
function varargout = recon_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filter
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');

set(hObject, 'String', {'No Filter','Ram-Lak','Shepp-Logan','Cosine','Hamming','Hann','Phase-contrast'});
popup_sel_index1 = get(hObject, 'Value');

switch popup_sel_index1
    case 1
        filterName = 'No Filter';
    case 2
        filterName = 'ram-lak';
%        set(handles.slider2,'Visible','on');
%         set(handles.edit9,'Visible','on');
%         set(handles.text12,'Visible','on');
    case 3
        filterName = 'shepp-logan';
%     %    set(handles.slider2,'Visible','on');
%         set(handles.edit9,'Visible','on');
%         set(handles.text12,'Visible','on');
    case 4
        filterName = 'cosine';
%    %     set(handles.slider2,'Visible','on');
%         set(handles.edit9,'Visible','on');
%         set(handles.text12,'Visible','on');
    case 5
        filterName = 'hamming';
    %    set(handles.slider2,'Visible','on');
%         set(handles.edit9,'Visible','on');
%         set(handles.text12,'Visible','on');
    case 6
        filterName = 'hann';
   %     set(handles.slider2,'Visible','on');
%         set(handles.edit9,'Visible','on');
%         set(handles.text12,'Visible','on');
    case 7
        filterName = 'Phase-contrast';
    %    set(handles.slider2,'Visible','on');
%         set(handles.edit9,'Visible','on');
%         set(handles.text12,'Visible','on');   
    case 8
%         filterName = 'New Filter';
   %     set(handles.slider2,'Value',1);
%          set(handles.edit9,'String','1');
%      %   set(handles.slider2,'Visible','off');
%         set(handles.edit9,'Visible','off');
%         set(handles.text12,'Visible','off');
end

du=proj_param.du;
nu=proj_param.N_col;
nv=proj_param.N_row;
%d = str2num(get(handles.edit9,'String'));
d=1; %fequency
if d<0||d>1
    errordlg('Value of d must be between 0 and 1.','Wrong Value');
else
    Phihat = OSCaRFilter( filterName, nu, du, d );
end
proj_param.Phihat=Phihat;

%Plotting the filter in use on axes6.

if ~isempty(Phihat)
    
    axes(handles.axes1);
    cla;
    plot(linspace(0,1,round(length(Phihat)/2)),Phihat(1:round(end/2)));
    axis tight;
    title([filterName]);
    ylabel('Filter Magnitude');
    xlabel('Frequency');
else
    
    axes(handles.axes1);
    cla;
    title('');
    
end

setappdata(h,'MyParam',proj_param);
% --- Executes during object creation, after setting all properties.

% --- Executes during object creation, after setting all properties.
function filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in execute_button.
function execute_button_Callback(hObject, eventdata, handles)
% hObject    handle to execute_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This is the "Execute" button. The reconstruction
%begins by clicking this button.

h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');
proj_data=getappdata(h,'MyData');


proj_param.x_size = str2num(get(handles.nx,'String'));
proj_param.y_size = str2num(get(handles.ny,'String'));
proj_param.z_size = str2num(get(handles.nz,'String'));

% Construct Projection Matrices
A = projectionMatrix(proj_param.theta,proj_param.du,proj_param.dv,proj_param.u_off,proj_param.v_off,proj_param.SDD,proj_param.SAD);

% Allocate space for Voxel Matrix and Projection matrix 

try
    R = zeros(length(proj_param.y), length(proj_param.x), length(proj_param.z));
    P = zeros(proj_param.N_row,proj_param.N_col);
    dR=R;
catch
    error('Out of Memory!');
end

% % Weight factors for non-equal angles
dtheta = [proj_param.theta(2)-proj_param.theta(1) ; ...
        ( proj_param.theta(3:end)-proj_param.theta(1:end-2) )/2 ; ...
        proj_param.theta(end)-proj_param.theta(end-1)]; 
dtheta_bar = mean(dtheta);
Wt = dtheta/dtheta_bar;
Wt = Wt./mean(Wt);
Wr = pi/proj_param.N_proj;


% d = str2num(get(handles.edit9,'String'));
d=1  %frequency
nu=proj_param.N_col;
nv=proj_param.N_row;

global button;    
divisor=proj_param.N_proj;

%"tag" is a global var that shows whether the computation must go on
% or we want to stop it. If we hit the "Stop" button, then the value of tag
%will be 0 and it exits the loop.
global tag;
tag=1;

axes(handles.axes4);
cla;
imagesc(proj_param.x,proj_param.y,R(:,:,round(proj_param.z_size/2)));
xlabel('x [cm]','Color',[0 0 1])
ylabel('y [cm]','Color',[0 0 1])
title('x vs. y [cm]');
axis xy;
axis equal;
axis tight;
colormap bone;
set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
hold on;

axes(handles.axes5);
cla;
R_cor = reshape( R(round(proj_param.y_size/2),:,:), proj_param.x_size,proj_param.z_size )';
imagesc( proj_param.x, proj_param.z, 1500 * ( R_cor-min(R_cor(:)) ) / ...
    ( max(R_cor(:))-min(R_cor(:)) ) );
xlabel('x [cm]','Color',[0 0 1])
ylabel('z [cm]','Color',[0 0 1])
title('x vs. z [cm]');
axis xy;
axis equal;
axis tight;
colormap bone;
set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
hold on;

axes(handles.axes6);
cla;
imagesc( proj_param.y, proj_param.z, reshape( ...
    R(:,round(proj_param.x_size/2),:),proj_param.y_size,proj_param.z_size)' );
xlabel('y [cm]','Color',[0 0 1])
ylabel('z [cm]','Color',[0 0 1])
title('y vs. z [cm]');
axis xy;
axis equal;
axis tight;
colormap bone;
set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
hold on;
set(handles.totalproj,'String',num2str(proj_param.N_proj));

%Begin reconstruction
for k=1:proj_param.N_proj
    
    if ( tag==1) % Checks whether the "Stop" button is hit.
        logP=proj_data.P(:,:,k);
        if Wr ~= 0 
            dR = phaseFDK( proj_param.x, proj_param.y, proj_param.z, proj_param.u_off(k), ...
                proj_param.v_off(k), proj_param.du, proj_param.dv, ...
                proj_param.theta(k), logP, A(:,:,k), ...
                proj_param.SDD, proj_param.SAD, proj_param.Phihat );
            R = R + Wr*dR;
            Cmin=min(min(min(R)));
            Cmax=max(max(max(R)));
        end
        
 %       if ~mod(k,5)
            axes(handles.axes4);
            cla;
            imagesc(proj_param.x,proj_param.y,R(:,:,round(proj_param.z_size/2)))
            caxis([Cmin Cmax]);

            axes(handles.axes5);
            cla;
%            R_cor = reshape( R(round(proj_param.y_size/2),:,:), proj_param.x_size,proj_param.z_size )';
%            imagesc( proj_param.x, proj_param.z, 1500 * ( R_cor-min(R_cor(:)) ) / ...
%                ( max(R_cor(:))-min(R_cor(:)) ) );
           imagesc(proj_param.x,proj_param.z,flipud(reshape(R(round(proj_param.y_size/2),:,:),...
               proj_param.x_size,proj_param.z_size))');
            caxis([Cmin,Cmax]);


            axes(handles.axes6);
            cla;
            imagesc( proj_param.y, proj_param.z, reshape( ...
                R(:,round(proj_param.x_size/2),:),proj_param.y_size,proj_param.z_size)' );
            caxis([Cmin Cmax]);

%             axes (handles.axes7);
%             cla;
%             plot(1,1);
%             axis off;
%             caxis([Cmin Cmax]);
%             colormap bone;
%             colorbar;

           
 %       end
        
%         if k==proj_param.N_proj
% 
%             set(handles.edit22,'String',num2str(k));
%             set(handles.edit25,'ForegroundColor',[1,0,0],'String','Complete');
%         else 
%            if ~mod(k,5)
%                 str=['Projection ',num2str(k),' of ', ...
%                     num2str(proj_param.N_proj),', \theta_G =', ...
%                     num2str(proj_param.theta(k),'%.1f'), ...
%                     '[deg]'];
%                
% 
%                 set(handles.status,'ForegroundColor',[0,0,0],'FontSize',10,'String',str);
                set(handles.projno,'String',num2str(k));
                set(handles.projdeg,'String',num2str(proj_param.theta(k)));
%            end
%        end

    end

end

Save_data.P=R;
Save_data.xmax=proj_param.xmax;
Save_data.xmin=proj_param.xmin;
Save_data.ymax=proj_param.ymax;
Save_data.ymin=proj_param.ymin;
Save_data.zmax=proj_param.zmax;
Save_data.zmin=proj_param.zmin;
setappdata(handles.export_button,'SaveData',Save_data);

% --- Executes on button press in cancel_button.
function cancel_button_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%This is the "Close" button.

%We should first stop the computation.
global tag;
tag=0; 
set(handles.stop_button,'Value',1);

% Get the current position of the GUI from the handles structure
% to pass to the modal dialog.
pos_size = get(handles.projGUI,'Position');

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
	delete(handles.projGUI)
end

% --- Executes on button press in export_button.
function export_button_Callback(hObject, eventdata, handles)
% hObject    handle to export_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');
% Save_data is a structure with different fields such as 
% P the reconstructed matrix and the reconstruction limits.
Save_data = getappdata(hObject,'SaveData');

if isempty(Save_data)
    errordlg('Reconstruction Not Complete!','Error');
else
    %Save button
    %By hitting this button the reconstructed data will be
    %saved in a .mat file.
    [mat_filename, DefaultPathname] = uiputfile( '*.mat', 'Save as','Untitled');
     FileName = strcat(DefaultPathname,mat_filename);
     try
        if  isequal(mat_filename,0)
            errordlg('Data not saved.','Mat File Not Saved');
        else
        str = [ 'save ''' FileName ''' Save_data' ];
        eval(str);
        msgbox('Data saved to a .mat file!','*.mat File');
        end
        catch
        errordlg('Data not saved. Try again!','Mat File Not Saved');
        end
end

function x_min_Callback(hObject, eventdata, handles)
% hObject    handle to x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_min as text
%        str2double(get(hObject,'String')) returns contents of x_min as a double
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');
proj_param.xmin = str2num(get(hObject, 'String'));
proj_param.x0 = proj_param.xmin/proj_param.dx+round(proj_param.u_off_ML(proj_param.k_AP));
setappdata(h,'MyParam',proj_param);
% call update function
update_projGUI(handles);

% --- Executes during object creation, after setting all properties.
function x_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_min_Callback(hObject, eventdata, handles)
% hObject    handle to y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_min as text
%        str2double(get(hObject,'String')) returns contents of y_min as a double
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');

proj_param.ymin = str2num(get(hObject, 'String'));
proj_param.y0 = proj_param.ymin/proj_param.dy+round(proj_param.u_off_ML(proj_param.k_RL));
setappdata(h,'MyParam',proj_param);
update_projGUI(handles);

% --- Executes during object creation, after setting all properties.
function y_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_min_Callback(hObject, eventdata, handles)
% hObject    handle to z_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_min as text
%        str2double(get(hObject,'String')) returns contents of z_min as a double
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');
proj_param.zmin = str2num(get(hObject, 'String'));
proj_param.z0 = proj_param.zmin/proj_param.dz+round(proj_param.v_off_ML(proj_param.k_RL));
setappdata(h,'MyParam',proj_param);
update_projGUI(handles);

% --- Executes during object creation, after setting all properties.
function z_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_max_Callback(hObject, eventdata, handles)
% hObject    handle to x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_max as text
%        str2double(get(hObject,'String')) returns contents of x_max as a double
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');

proj_param.xmax = str2num(get(hObject, 'String'));
proj_param.x1 = proj_param.xmax/proj_param.dx+round(proj_param.u_off_ML(proj_param.k_AP));
setappdata(h,'MyProj',proj_param);
update_projGUI(handles);

% --- Executes during object creation, after setting all properties.
function x_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_max_Callback(hObject, eventdata, handles)
% hObject    handle to y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_max as text
%        str2double(get(hObject,'String')) returns contents of y_max as a double
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');
proj_param.ymax = str2num(get(hObject, 'String'));
proj_param.y1 = proj_param.ymax/proj_param.dy+round(proj_param.u_off_ML(proj_param.k_RL));
setappdata(h,'MyParam',proj_param);
update_projGUI(handles);

% --- Executes during object creation, after setting all properties.
function y_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_max_Callback(hObject, eventdata, handles)
% hObject    handle to z_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_max as text
%        str2double(get(hObject,'String')) returns contents of z_max as a double
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');
proj_param.zmax = str2num(get(hObject, 'String'));
proj_param.z1 = proj_param.zmax/proj_param.dz+round(proj_param.v_off_ML(proj_param.k_RL));
setappdata(h,'MyParam',proj_param);
update_projGUI(handles);

% --- Executes during object creation, after setting all properties.
function z_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fov_select.
function fov_select_Callback(hObject, eventdata, handles)
% hObject    handle to fov_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%This is "Enter the borders graphically" button.
%Once you hit it, it gives you the mouse to enter the bordes
%graphically.

% set(hObject,'Visible','on');
% %=========================
h = findobj('Tag','projGUI');
proj_param=getappdata(h,'MyParam');


[U,V] = ginput(2);
proj_param.y0 = max(1,round(U(1)));
proj_param.y1 = min(max(proj_param.y0,round(U(2))),proj_param.N_col);

proj_param.z0 = max(1,round(V(1)));
proj_param.z1 = min(max(proj_param.z0,round(V(2))),proj_param.N_row);

yminstr=num2str(proj_param.dy*( proj_param.y0-round(proj_param.u_off_ML(proj_param.k_RL))));
set(handles.y_min, 'String', yminstr);

ymaxstr=num2str(proj_param.dy*( proj_param.y1-round(proj_param.u_off_ML(proj_param.k_RL))));
set(handles.y_max, 'String', ymaxstr);

zminstr=num2str(proj_param.dz*( proj_param.z0-round(proj_param.v_off_ML(proj_param.k_RL))));
set(handles.z_min, 'String', zminstr);

zmaxstr=num2str(proj_param.dz*( proj_param.z1-round(proj_param.v_off_ML(proj_param.k_RL))));
set(handles.z_max, 'String', zmaxstr);
set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
pause(0.5)


setappdata(h,'MyParam',proj_param);
update_projGUI(handles);

set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
% Select region of interest for rest of image set
[U,V] = ginput(2);
proj_param.x0 = max(1,round(U(1)));
proj_param.x1 = min(max(proj_param.x0,round(U(2))),proj_param.N_col);
set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
pause(0.5)

xminstr=num2str(proj_param.dx*(proj_param.x0-round(proj_param.u_off_ML(proj_param.k_AP))));
set(handles.x_min, 'String', xminstr);

xmaxstr=num2str(proj_param.dx*(proj_param.x1-round(proj_param.u_off_ML(proj_param.k_AP))));
set(handles.x_max, 'String', xmaxstr);

setappdata(h,'MyParam',proj_param);
update_projGUI(handles);


% --- Executes on button press in stop_button.
function stop_button_Callback(hObject, eventdata, handles)
% hObject    handle to stop_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global tag;
if get(hObject,'Value')==1
       tag=0;
end
