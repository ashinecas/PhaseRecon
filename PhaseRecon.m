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

if verboseTag==1
    subplot(1,2,1)
    imagesc(proj_data.log_P0,m+s*[-3,3])
    axis xy
    axis equal;
    axis tight;
    colormap(gray(236))
    xlabel('u [pixels]','Color',[0 0 1])
    ylabel('v [pixels]','Color',[0 0 1])
    title(['Projection #1 - ',num2str(proj_param.theta(1)),' degrees'])

    hold on
    plot(proj_param.u_off_ML(proj_param.k_RL),proj_param.v_off_ML(proj_param.k_RL),'rx','MarkerSize',18)

    subplot(1,2,2)

    proj_param.u_aper = round(proj_param.u_off_ML(proj_param.k_AP))+[-20:20];
    proj_param.v_aper = round(proj_param.v_off_ML(proj_param.k_AP))+[-20:20];

    m = mean2(proj_data.log_P1(proj_param.v_aper,proj_param.u_aper));
    s = std2(proj_data.log_P1(proj_param.v_aper,proj_param.u_aper));

    imagesc(proj_data.log_P1,m+s*[-3,3])
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
