function recon = PhaseRecon(projfile,varain)
%% Independent Reconstruction mode. Non-GUI mode.
%Syntax
% projfile is created by  projpreproc
% varain
    % filter   (default: R-L)
    % FOV
         

%% Input arguments initialization
% if filterName not exist
        filterName='R-L';
  %else
  %    filterName= filterName;
graphicalTag=1; 
VoxelTag=0;
verboseTag=1; %Display Reconstruction process.
d=1;
%%  Load projection parameters and data. projfile is created by projpreproc();
load(projfile);
proj_param = experiment.param;
proj_data = experiment.data;
%%  Display sagital and cornal projection (also used for graphical FOV selection)

ind_theta = find(abs(diff(proj_param.theta))>0.1);
[t,ind2] = unique(proj_param.theta(ind_theta));
proj_param.k_RL = max(1,min(round(interp1(t,ind_theta(ind2),0,'nearest','extrap')),proj_param.N_proj))
                    %0度投影编号
proj_param.k_AP = max(1,min(round(interp1(t,ind_theta(ind2),90,'nearest','extrap')),proj_param.N_proj));
                    %90度投影编号

proj_data.log_P0=proj_data.P(:,:,proj_param.k_RL);%0度投影
proj_data.log_P1=proj_data.P(:,:,proj_param.k_AP);%90度投影

proj_param.u_off_ML = proj_param.u_off;%horizontal offset
proj_param.v_off_ML = proj_param.v_off;%vertical offset
proj_param.MF = 1 + (proj_param.IAD)./(proj_param.SAD);

proj_param.u_aper = round(proj_param.u_off_ML(proj_param.k_RL))+[-20:20];
proj_param.v_aper = round(proj_param.v_off_ML(proj_param.k_RL))+[-20:20];
m = mean2(proj_data.log_P0(proj_param.v_aper,proj_param.u_aper));
s = std2(proj_data.log_P0(proj_param.v_aper,proj_param.u_aper));

if verboseTag==1
    subplot(1,2,1)
    imagesc(proj_data.log_P0,m+s*[-3,3])%Show 0 degree projection
    axis xy
    axis equal;
    axis tight;
    colormap(gray(236))
    xlabel('u [pixels]','Color',[0 0 1])
    ylabel('v [pixels]','Color',[0 0 1])
    title(['Projection #1 - ',num2str(proj_param.theta(1)),' degrees'])

    hold on
    plot(proj_param.u_off_ML(proj_param.k_RL),proj_param.v_off_ML(proj_param.k_RL),'rx','MarkerSize',18)
                % mark on offset position

    subplot(1,2,2)

    proj_param.u_aper = round(proj_param.u_off_ML(proj_param.k_AP))+[-20:20];
    proj_param.v_aper = round(proj_param.v_off_ML(proj_param.k_AP))+[-20:20];

    m = mean2(proj_data.log_P1(proj_param.v_aper,proj_param.u_aper));
    s = std2(proj_data.log_P1(proj_param.v_aper,proj_param.u_aper));

    imagesc(proj_data.log_P1,m+s*[-3,3])%Show 90 degree projection
    axis xy
    axis equal;
    axis tight;

    xlabel('u [pixels]','Color',[0 0 1])
    ylabel('v [pixels]','Color',[0 0 1])
    title(['Projection #',num2str(proj_param.k_AP),' - ',num2str(proj_param.theta(proj_param.k_AP)),' degrees'])
    hold on
    plot(proj_param.u_off_ML(proj_param.k_AP),proj_param.v_off_ML(proj_param.k_AP),'rx','MarkerSize',18)
    hold off
    set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
end
%% Set Voxel size (autoset or user-defined)
if VoxelTag==0
    proj_param.dx = round(100*proj_param.du/mean(proj_param.MF))/100;
    proj_param.dy = round(100*proj_param.du/mean(proj_param.MF))/100;
    proj_param.dz = round(100*proj_param.dv/mean(proj_param.MF))/100;
else
    proj_param.dx=VoxelSize(1);
    proj_param.dy=VoxelSize(2);
    proj_param.dz=VoxelSize(3);
end
%% Set FOV (Graphical input or user-defined)
if graphicalTag==1 && verboseTag==1
    %Enter the borders graphically...
    [U,V] = ginput(2);
    y0 = max(1,round(U(1)));
    y1 = min(max(y0,round(U(2))),proj_param.N_col);

    z0 = max(1,round(V(1)));
    z1 = min(max(z0,round(V(2))),proj_param.N_row);


    set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
    % Select region of interest for rest of image set
    [U,V] = ginput(2);
    x0 = max(1,round(U(1)));
    x1 = min(max(x0,round(U(2))),proj_param.N_col);
    set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
elseif graphicalTag==1 && verboseTag==0
    error('Verbose is off. Therefore you need to enter the Reconstruction Volume.');
else
    %Use the entered values az reconstruction borders...
    xmin=RecVol(1);
    xmax=RecVol(2);
    ymin=RecVol(3);
    ymax=RecVol(4);
    zmin=RecVol(5);
    zmax=RecVol(6);
    x0=xmin/proj_param.dx+round(proj_param.u_off_ML(proj_param.k_AP));
    x1=xmax/proj_param.dx+round(proj_param.u_off_ML(proj_param.k_AP));
    y0=ymin/proj_param.dy+round(proj_param.u_off_ML(proj_param.k_RL));
    y1=ymax/proj_param.dy+round(proj_param.u_off_ML(proj_param.k_RL));
    z0=zmin/proj_param.dz+round(proj_param.v_off_ML(proj_param.k_RL));
    z1=zmax/proj_param.dz+round(proj_param.v_off_ML(proj_param.k_RL));
end

if verboseTag==1
    subplot(1,2,1)
    hold on
    plot([y0;y1],[z0;z0],'c-')
    plot([y0;y1],[z1;z1],'c-')
    plot([y0;y0],[z0;z1],'c-')
    plot([y1;y1],[z0;z1],'c-')
    hold off


    set(gca,'XColor',[0 0 1],'YColor',[0 0 1])


    subplot(1,2,2)
    hold on
    plot([x0;x1],[z0;z0],'c-')
    plot([x0;x1],[z1;z1],'c-')
    plot([x0;x0],[z0;z1],'c-')
    plot([x1;x1],[z0;z1],'c-')


    set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
    pause(0.5)
end

proj_param.x = -flipud(proj_param.dx*( [x0:x1]'-round(proj_param.u_off_ML(proj_param.k_AP))));         % Account for A-P view since +ve u is -ve x
proj_param.y = proj_param.dy*( [y0:y1]'-round(proj_param.u_off_ML(proj_param.k_RL)));
proj_param.z = proj_param.dz*( [z0:z1]'-round(proj_param.v_off_ML(proj_param.k_RL)));

proj_param.x_size = length(proj_param.x);
proj_param.y_size = length(proj_param.y);
proj_param.z_size = length(proj_param.z);
