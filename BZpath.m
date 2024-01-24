function [qpath,tix,til]=BZpath(step,qs,qe,qs_str,qe_str,varargin)

%% GENERATE PATH IN BRILLOUIN ZONE (BZ)
% - in 2*pi/spacing units (where spacing = lattice spacing);
% - according to traditional solid state physics representation;
% - providing a fictitious linearly increasing coordinate
%   for the purpose of a later graphical representation 
%   of the dispersion relations of quasiparticles;
% - providing abscissa tick labels and their positions for
%   plotting.

%% INPUT
%  step = Step along path in BZ
%  qs(1:3,i) = 3XN matrix = Start point of segment i=1:N
%  qe(1:3,i) = 3XN matrix = End   point of segment i=1:N
%  qs_str{i} = label of start     point of segment i=1:N
%  qe_str{i} = label of end       point of segment i=1:N
 
%% OUTPUT
%  qpath(1:3,i) = Cartesian coordinates of wavevector i
%  qpath(4,i) = magnitude of wavevector i 
%  qpath(5,i) = fictitious linearly increasing coordinate 
%               for the purpose linear graphical representation 
%               of vector qpath(1:3,i)
%  qpath(6,i) = segment identifier of vector qpath(1:3,i)
%
%  til{1:2*N+1} = string array = tick labels along abscissa axis 
%  tix(1:2*N+1) = tick positions of abscissa labels
%
%  The arrays tix and til shall be used as follows in the script
%  drawing the quasiparticle dispersion relations after the plot
%  command:
%  set (gca,'xtick',tix);
%  set (gca,'xticklabel',til,'FontSize',10,'FontWeight','bold'); 

%% Recognized options in varargin 
% (uppercases for readability are optional): 
%
% 'PlotPath'       to plot the path segments in 3D plot
%
% if varargin{k} = 'File', then varargin{k+1} = Character string 
%                          = name of output file
%                                     
% if varargin{k} = 'Colors', then varargin{k+1 } = 7x3 matrix
%                            = Color order of successive graphical objects
% if varargin{k} = 'ShowLegend', 
%                   then varargin{k+1} = Boolean to show legend or not

%% DEFAUT VALUES OF OPTIONAL ARGUMENTS

plot_path=false; default_color=true; show_legend=true;
outputfile=[];

%% PARSE OPTIONAL ARGUMENT LIST

name_value_pair=false;
for k = 1:length(varargin);
    if (name_value_pair)
        name_value_pair=false;
    else
        switch lower(varargin{k}) % varargin is a "cell array"
          case {'file'}
            name_value_pair=true; outputfile=varargin{k+1}; 
          case {'plotpath'}
            plot_path=true;
          case {'colors'}
            color_scheme = varargin{k+1}; default_color=false; 
            name_value_pair=true;
          case {'showlegend'}
            show_legend=varargin{k+1}; name_value_pair=true;
          otherwise
            error('function BZpath: option %s not recognized.\n',...
                  varargin{k}); return;
        end
    end
end

%% CORE JOB

[~,s] = size(qs); [~,t] = size(qe); 
if ( s ~= t  )
    error('BZpath: numbers of start and end points do not match.\n');
else
    fprintf('\nFunction BZpath: %d segments\n',s);
end

%% Initialization

qu=qe-qs;     % Direction vectors
qu=sign(qu);  % Traditional solid state physics
	      % cristallography vectors along each segment

for is=1:s
    fprintf('Segment %d: direction [%d %d %d]''\n',is,qu(1:3,is))
end

q=zeros(5,1); % Initialize q vector  
p=zeros(5,1); % Initialize vector preceding q
iq=0;         % Initialize number of points for the path in BZ  

tix(1)=0; til{1}=qs_str{1}; l=1; % Label at the abscissa origin

%% Loop overs segments

for is=1:s    
  ip=0; cond=0;
  while (cond != 1)
    xsi=ip*step;
    buf=xsi*qu(:,is);
    % cond=1 if all vector components are non-zero 
    %        <=> if going beyond end point
    cond=all(abs(buf)>=abs((qe(:,is)-qs(:,is)))); 
    
    iq=iq+1;
    if(cond==1)           % Last point of segment
      q=qe(1:3,is); q(4)=sqrt(q'*q);
      gap=max(abs(qe(1:3,is)-p(1:3)));
      q(5)=p(5)+gap;
      % Manage labels when last point of segment is found
      l=l+1; tix(l)=tix(l-1)+(q(5)-tix(l-1))/2; 
      til{l}=sprintf('\n\n[%d,%d,%d]',qu(1,is),qu(2,is),qu(3,is));
      l=l+1; tix(l)=q(5); til{l}=qe_str{is}; 
    else                  % Standard point of segment
      q=buf(1:3)+qs(1:3,is); q(4)=sqrt(q'*q);
      if(iq == 1)         % First point of path in BZ
	q(5)=0; p(5)=0;
      else
	if(ip == 0) 
	  q(5)=p(5);      % First point of segment
	else
	  q(5)=p(5)+step; % Standard point of segment
	end
      end
    end
    p=q; ip=ip+1;
    qpath(1:5,iq)=q(1:5);
    qpath(6,iq)=is;       % Segment identifier of point iq
  end % End while
  
end % End of loop over segments

%% Formatted output

if ~isempty(outputfile)
    printtable(qpath','LineNumber',true,'Integer',[6],'ArrayName','q',...
               'file',outputfile)
end

%% OPTIONAL PLOT

if(plot_path)

    h = findall(0,'type','figure');
    if isempty(h)
       fig=figure('NumberTitle','off','name','Path in Brillouin Zone');
    end
    if (default_color) 
        color_scheme = get(gca,'colororder');
    end
    
    for i=1:s
        plot3(qs(1,i),qs(2,i),qs(3,i),...
              'o','color',color_scheme(i,1:3),...
              'DisplayName',qs_str{i}); hold on; 
    end
    
    for i=1:s
        pt(1:3,1)=qs(1:3,i); pt(1:3,2)=qe(1:3,i);
        plot3(pt(1,:),pt(2,:),pt(3,:),...
              '-','color',color_scheme(i,1:3),...
              'DisplayName',sprintf('segment %d',i)); hold on; 
    end
    
    axis equal; grid on; legend('Location','EastOutside');
    xlabel('q_1'); ylabel('q_2'); zlabel('q_3');
    if(~show_legend)
        legend(gca,'off');
    end
end

end  % End of function BZpath
