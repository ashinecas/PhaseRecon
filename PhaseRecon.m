function R = PhaseRecon(projfile,varain)
%% Independent Reconstruction mode. Non-GUI mode.
%Syntax
% projfile is created by  projpreproc
% varain
    % filter   (default: R-L)
    % FOV
         

%% Input arguments initialization
% if filterName not exist
        %filterName='ram-lak';
        filterName='Phase-contrast';
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
proj_param.k_RL = max(1,min(round(interp1(t,ind_theta(ind2),0,'nearest','extrap')),proj_param.N_proj));
                    %0 degree projection No.
proj_param.k_AP = max(1,min(round(interp1(t,ind_theta(ind2),90,'nearest','extrap')),proj_param.N_proj));
                    %90degree projection No.

proj_data.log_P0=proj_data.P(:,:,proj_param.k_RL);%0 degree projection
proj_data.log_P1=proj_data.P(:,:,proj_param.k_AP);%90 degree projection

proj_param.u_off_ML = proj_param.u_off;%horizontal offset
proj_param.v_off_ML = proj_param.v_off;%vertical offset
proj_param.MF = 1 + (proj_param.IAD)./(proj_param.SAD);%Magnitude ratio

proj_param.u_aper = round(proj_param.u_off_ML(proj_param.k_RL))+[-20:20];
proj_param.v_aper = round(proj_param.v_off_ML(proj_param.k_RL))+[-20:20];
m = mean2(proj_data.log_P0(proj_param.v_aper,proj_param.u_aper));
s = std2(proj_data.log_P0(proj_param.v_aper,proj_param.u_aper));

if verboseTag==1 %Display reconstruction process
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

%% Set Filter
du=proj_param.du;
nu=proj_param.N_col;
nv=proj_param.N_row;

proj_param.Phihat = Filter( filterName, nu, du, d );
disp(['Using filter ', filterName,' with d=',num2str(d)]);

%% Create Projection Matrix
A = projectionMatrix(proj_param.theta,proj_param.du,proj_param.dv,proj_param.u_off,proj_param.v_off,proj_param.SDD,proj_param.SAD);

%% Allocate space for Reconstruction Matrix
R = zeros(length(proj_param.y), length(proj_param.x), length(proj_param.z));

%% Weighting for non-central slices.
dtheta = [proj_param.theta(2)-proj_param.theta(1) ; ( proj_param.theta(3:end)-proj_param.theta(1:end-2) )/2 ; proj_param.theta(end)-proj_param.theta(end-1)]; 
dtheta_bar = mean(dtheta);
Wt = dtheta/dtheta_bar;
Wt = Wt./mean(Wt);
Wr = pi/proj_param.N_proj;

%% Begin Reconstruction!
nu=proj_param.N_col;
nv=proj_param.N_row;

for k=1:proj_param.N_proj
        logP=proj_data.P(:,:,k);%load projection
        dR = phaseFDK( proj_param.x, proj_param.y, proj_param.z, proj_param.u_off(k), proj_param.v_off(k), proj_param.du, proj_param.dv, ...
            proj_param.theta(k), logP, A(:,:,k), ...
            proj_param.SDD, proj_param.SAD, proj_param.Phihat );
        R = R + Wr*dR;
     %display Reconstruction processes
        if verboseTag==1
            if ~mod(k,5)


                figure(2)
                subplot(2,2,1) %axial image
                imagesc(proj_param.x,proj_param.y,R(:,:,round(proj_param.z_size/2)))
                xlabel('x [cm]','Color',[0 0 1])
                ylabel('y [cm]','Color',[0 0 1])
                title('Axial','FontSize',12,'Color',[0 0 1]);

                axis xy;
                axis equal;
                axis tight;
                colormap bone;
                set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
                colorbar;


                hold on;
                subplot(2,2,3) % Cornal image
                R_cor = reshape( R(round(proj_param.y_size/2),:,:), proj_param.x_size,proj_param.z_size )';
                imagesc( proj_param.x, proj_param.z, 1500 * ( R_cor-min(R_cor(:)) ) / ...
                    ( max(R_cor(:))-min(R_cor(:)) ) );
                xlabel('x [cm]','Color',[0 0 1])
                ylabel('z [cm]','Color',[0 0 1])
                title('Coronal','FontSize',12,'Color',[0 0 1]);
                hold on;
                axis xy;
                axis equal;
                axis tight;
                colormap bone;
                set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
                colorbar;

                hold on;
                subplot(2,2,4) %Sagital Image
                imagesc( proj_param.y, proj_param.z, reshape( ...
                    R(:,round(proj_param.x_size/2),:),round(proj_param.y_size),round(proj_param.z_size))' );
                xlabel('y [cm]','Color',[0 0 1])
                ylabel('z [cm]','Color',[0 0 1])
                title('Sagital','FontSize',12,'Color',[0 0 1]);
                axis xy;
                axis equal;
                axis tight;
                colormap bone;
                set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
                colorbar;
                hold on;

                subplot(2,2,2) %Projection image
                proj_param.u_aper = round(proj_param.u_off_ML(k))+[-20:20];
                proj_param.v_aper = round(proj_param.v_off_ML(k))+[-20:20];
                m = mean2(logP(proj_param.v_aper,proj_param.u_aper));
                s = std2(logP(proj_param.v_aper,proj_param.u_aper));

                imagesc(logP,m+s*[-3,3])
                hold on
                plot(proj_param.u_off_ML(k),proj_param.v_off_ML(k),'rx','MarkerSize',14,'LineWidth',1.5)
                hold off
                xlabel('u [pixels]','Color',[0 0 1])
                ylabel('v [pixels]','Color',[0 0 1])
                title(['Projection ',num2str(k),' of ',num2str(proj_param.N_proj),', \theta_G=',num2str(proj_param.theta(k),'%.1f'),'[deg]'],'Color',[0 0 1],'FontSize',10)

                axis xy
                axis equal;
                axis tight;
                colormap bone;
                colormap(gray(236))
                set(gca,'XColor',[0 0 1],'YColor',[0 0 1])
                pause(0.2)
            end

        end
    
end

