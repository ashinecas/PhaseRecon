function filt = Filter( filter, N, dx, d )
% filt = OSCaRFilter( filter, N, dx, d )
% Nargol Rezvani
% D. Aruliah, July 2007
% Returns the Fourier Transform of the filter which will be 
% used to filter the projections
%
% INPUT ARGS:   filter - string specifying the filter 
%               N      - length of the vector of projection data
%               dx     - spatial resolution of projection data
%               d      - fraction of frequencies below the nyquist
%                        which we want to pass (0<d<=1)
%
% OUTPUT ARGS:  filt   - filter to use on the projections
% Notes:
% (1) Scaling of time & frequency domains computed according to reciprocity
%     relations as in Briggs & Henson's book
%     "The DFT: An Owners' Manual for the Discrete Fourier Transform" (1995)
% (2) Definitions of windows from James S. Walker's book
%     "Fast Fourier Transforms" (1996, 2nd ed)

% No work to do if no filter required
if strcmp(filter,'No Filter')
   filt = [];
   return
end

Nfilt = max( 64, 2^nextpow2( 2*N ) );  % Length of desired filter
Omega = 1/dx;                          % Omega = 2 x Nyquist frequency
domega = Omega/Nfilt;                  % Effective frequency resolution
omega = domega*(0:Nfilt/2);            % Effective frequency domain
filt = omega;                          % Initialise ramp filter

% Generate windowed filter according to choice
switch filter
    case 'ram-lak'
         % Do nothing
    case 'shepp-logan'
          filt(2:end) = filt(2:end).*(sin(pi*omega(2:end)/(d*Omega))...
                        ./ (pi*omega(2:end)/(d*Omega)));
    case 'cosine'
          filt(2:end) = filt(2:end) .* cos(pi*omega(2:end)/(d*Omega));
    case 'hamming'  
          filt(2:end) = filt(2:end) ...
                        .* (.54 + .46*cos(2*pi*omega(2:end)/(d*Omega)));
    case 'hann'
          filt(2:end) = 0.5* filt(2:end) ...
                        .* (1 + cos(2*pi*omega(2:end)/(d*Omega)));
    case 'New Filter'
        prompt={'Enter the apodizing window as a function of variable x that will be multiplied by the ramp filter (x is a vector):'};
        name='Input for the apodizing window';
        numlines=1;
        defaultanswer={''};
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        
        answer=inputdlg(prompt,name,numlines,defaultanswer,options);
        str=char(answer);
        if ~isempty(answer)
            try
                apod=OSCaRApod(str,omega(2:end));
                filt(2:end) = filt(2:end).*apod;
            catch
                errordlg('The apodizing window must be a function of x.','Wrong Input');
                filt = [];
                return
            end
        else
            filt=[];
            return
        end
        
        
        % so basically here we want to convert this input to a function of x
        % and then replace x by 2*pi*omega(2:end)/(d*Omega) and use it as
        % the apodizing function and multiply it by filt(2:end)
    otherwise
          eid = sprintf('Images:%s:invalidFilter',mfilename);
          msg = 'Invalid filter selected.';
          error(eid,'%s',msg);
end
filt(omega>0.5*Omega*d) = 0;       % Crop frequencies at d*Nyquist freq
filt = [ filt, filt(end-1:-1:1) ]; % Frequencies ordered [0:N/2,-1:-N/2+1]
return

function apodV=OSCaRApod(str,x)
apodV=eval(str);