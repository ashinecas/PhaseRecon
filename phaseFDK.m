function dR = phaseFDK(x,y,z,u_off,v_off,du,dv,Beta,P,A,SDD,SAD,filter)
%% Description: Reconstruct single Cone-Beam projection using FDK algorithm
%       x,y,z:      The 3-D grid for reconstructing the image. 
%       u_off, v_off: Offset position in each Projection. 
%                       If there is no offset:  (u_off,v_off)  -> (0.5*len(u),0.5*len(v))
%       du, dv:     pixel sizes  (unit: mm)
%       Beta:   rotation angle
%       P:          single projection (nu,nv).
%       A:          Projection Matrix
%       SDD:         source to detector distance
%       SAD:        source to axis distance
%       filter:     1D filter
%  Ref: OSCARFDK

%%  Dimensions of reconstructed Volume
nx = length(x);  
ny = length(y);
nz = length(z);
%% Non-central slice weighting 
% TODO: We can do this weighting in preceding steps.
%       since the weight function in our system is the same for every projection

% Remark: There may be some savings in fft computation if P is transposed
[nv,nu] = size(P); 
u = (-u_off + [0:nu-1]')*du;% u coordinate of pixels  (unit: mm)  vector
v = (-v_off + [0:nv-1]')*dv;% v coordinate of pixels  (unit: mm)  vector
% Remark: if P is transposed, 'meshgrid' should be changed to 'ndgrid'
[uu,vv] = meshgrid(u,v); % the coordinate of each pixel is specified by [uu,vv].

weight = SDD./sqrt(SDD^2 + vv.^2 + uu.^2);%non-central slice modify ->cosine weighting
                                          % other vesion is replace SDD by  SAD 
                                                                                                                         

P = P .* weight;  % weighted projection

%%   Filtering:  TODO: we can move the part in preceding steps.
%Phihat = varargin{1};
%if ~isempty(Phihat), % Filtered
if ~isempty(filter)
    % Remark: This call to fft can likely be optimised using the transpose
    %         of P instead to capitalise on Matlab's column-oriented storage.
    P = fft( P, length(filter), 2 ); % Rows zero-padded automatically
    % Remark: This loop can likely be vectorised using repmat
    for j=1:nv
        P(j,:) = P(j,:) .* filter ;
    end;
    % Remark: This call to fft can likely be optimised using the transpose
    %         of P instead to capitalise on Matlab's column-oriented storage.
    P = real( ifft( P, [], 2 ) );
    P = P(:,1:nu); % Trim zero-padding on rows
end % of if-block
%% allocate space for reconstruction volume

% Increment to add to reconstruction in backprojection stage
dR = zeros(ny,nx,nz);


[yy, xx, zz] = ndgrid( y, x, z); % 3D coordinate of dexels.

% Use projection matrix to project reconstruction (x,y,z) grid
% into detector (u,v) grid
UV = A * [ xx(:)'; yy(:)'; zz(:)'; ones(1,nx*ny*nz) ];% Matrices multiplication vectorlization. 
% UV -> homogeneous coordinates

% Nearest neighbour interpolation to find detector coordinates (u,v)
% that are closest to projections of voxels (x,y,z)
U = round( UV(1,:)./UV(3,:) ) + 1 ;
V = round( UV(2,:)./UV(3,:) ) + 1 ;

% Identify indices of voxels with projections strictly within detector grid
ind = find( (U>=1) & (V>=1) & (U<=nu) & (V<=nv) );

P_ind = ( U(ind) - 1 )*nv + V(ind) ; %

co = cos(Beta*pi/180);
si = sin(Beta*pi/180);
% Grid-dependent weighting factors (only computed for voxels needed)
W = ( SAD ./ ( SAD - co*xx(ind) - si*yy(ind) ) ).^2;
dR( ind ) = W .* P( P_ind );
return
