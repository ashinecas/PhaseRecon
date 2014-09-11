function dR = phaseFDK(x,y,z,u_off,v_off,du,dv,Beta,P,A,SDD,SAD,varargin)
%% dR = OSCaRFDK(x,y,z,u_off,v_off,du,dv,Beta,P,A,Xs,SDD,PhiHat);
%%*************************************************************************
%% Module name:      OSCaRFDK.m
%% Nargol Rezvani
%% May 2007
%% Revised July 2007, D. Aruliah
%%  Inputs:
%%      Beta:       the projection angle which is 
%%                  fixed for the whole function.
%%      P:          2-D Matrix of size nu*nv
%%      u_off, v_off: give us the farthest point from
%%                  the origin of the detector 
%%                  in horizontal and vertical
%%                  axes respectively.
%%                  {In other words, our detector is of size
%%                  [2*u_off]*[2*v_off]}(?)
%%      du, dv:     step size
%%      x,y,z:      The 3-D grid for reconstructing the image. 
%%      A:          Projection Matrix
%%      SDD:        distance from source to detector
%%      radius:     radius of the field-of-view
%%      Xs:         ???
%%      Phihat:       1-D Ramp Filter, Different options of filters 
%%                       -Will be added later.
%%                       -will be available (in Fourier space)
%%                       -It has to be of the same length as u
%%                       -length(ghat)=nu
%%
%%==================================================================
%%
%%  Outputs:
%%      After fixing the projection angle, the routine returns the
%%      value calculated by FDK equation for each voxel.
%%  Description:
%%      Reconstruct single Cone-Beam projection using 3D Feldkamp

% Dimensions of reconstruction grid
nx = length(x);
ny = length(y);
nz = length(z);

% Remark: There may be some savings in fft computation if P is transposed
[nv,nu] = size(P);
u = (-u_off + [0:nu-1]')*du;
v = (-v_off + [0:nv-1]')*dv;
% Remark: if P is transposed, 'meshgrid' should be changed to 'ndgrid'
[uu,vv] = meshgrid(u,v);
weight = SDD./sqrt(SDD^2 + vv.^2 + uu.^2);

% Remark: P is over-written to conserve memory in following computations.
%         The interpretation at each stage should be clear from context
P = P .* weight;

% Filtering:
Phihat = varargin{1};
if ~isempty(Phihat), % Skip to backprojection is filter is empty
    % Remark: This call to fft can likely be optimised using the transpose
    %         of P instead to capitalise on Matlab's column-oriented storage.
    P = fft( P, length(Phihat), 2 ); % Rows zero-padded automatically
    % Remark: This loop can likely be vectorised using repmat
    for j=1:nv
        P(j,:) = P(j,:) .* Phihat ;
    end;
    % Remark: This call to fft can likely be optimised using the transpose
    %         of P instead to capitalise on Matlab's column-oriented storage.
    P = real( ifft( P, [], 2 ) );
    P = P(:,1:nu); % Trim zero-padding on rows
end % of if-block

% Increment to add to reconstruction in backprojection stage
dR = zeros(ny,nx,nz);

% Vectorised computation of backprojection

[yy, xx, zz] = ndgrid( y, x, z);
% Use projection matrix to project reconstruction (x,y,z) grid
% into detector (u,v) grid
UV = A * [ xx(:)'; yy(:)'; zz(:)'; ones(1,nx*ny*nz) ];

% Nearest neighbour interpolation to find detector coordinates (u,v)
% that are closest to projections of voxels (x,y,z)
U = round( UV(1,:)./UV(3,:) ) + 1 ;
V = round( UV(2,:)./UV(3,:) ) + 1 ;

% Identify indices of voxels with projections strictly within detector grid
ind = find( (U>=1) & (V>=1) & (U<=nu) & (V<=nv) );
% Remark: This will have to be changed if P is transposed
P_ind = ( U(ind) - 1 )*nv + V(ind) ;

co = cos(Beta*pi/180);
si = sin(Beta*pi/180);
% Grid-dependent weighting factors (only computed for voxels needed)
W = ( SAD ./ ( SAD - co*xx(ind) - si*yy(ind) ) ).^2;
dR( ind ) = W .* P( P_ind );
return
