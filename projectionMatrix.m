function ProjM = projectionMatrix(thetaP,du, dv, uoff, voff, SDD,SAD)

N= length(thetaP);

ProjM = zeros(3,4,N);

for i = 1:N
    
    %K camera intrinsic calibration matrix  3*3   u = K*u_i   P^2->P^2
%     K=[1/du 0   -uoff(i)
%         0   g   -voff(i)
%         0   0   1];

    K =[     1/du    0         uoff(i)  
             0       1/dv      voff(i) 
             0        0        1    ] ;
    f=SDD;
    %   u_i = P*X_c  P^3->P^2
    P= [-f 0 0 0
        0 -f 0 0
        0 0 1 0];
    
    
    %extrinsic camerac calibration matrix   X_c = M*X   P^3 -> P^3
    
    %R=rotationMatrix(thetaP(i),0,0); %% TODO  ?
    R=rotationMatrix(90,0,thetaP(i)+90);  %Rotation matrix 
    t=[SAD*cosd(thetaP(i)); SAD*sind(thetaP(i)); 0];  % Translation matrix
    
    M =[R -R*t
        zeros(1,3) 1];  
    
  
    ProjM(:,:,i)=K*P*M;
    
    ProjM(:,:,i) = ProjM(:,:,i)./ProjM(3,4,i);
    
    
end

end

function R=rotationMatrix(phi,theta,psi)

Theta = theta*pi/180;% x rotation angle
Phi = phi*pi/180;   %y rotation  angle
Psi = psi*pi/180;   % z rotation angle

R = [ 1           0         0
      0           cos(Phi)  sin(Phi)
      0          -sin(Phi)  cos(Phi)]  ...
    * ...
    [ cos(Theta)  0        -sin(Theta)
      0           1         0
      sin(Theta)  0         cos(Theta)] ...
    * ...
    [ cos(Psi)    sin(Psi)  0
     -sin(Psi)    cos(Psi)  0
      0           0         1         ];
end

