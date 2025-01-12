function [smooth_value,dd]=inverse_smooth(file_asc_surf,varargin)
%
% inverse_smooth	smoothing surface values 
%
% [smooth_value,dd]=inverse_smooth(file_asc_surf,varargin)
% file_asc_surf		surface file in ASCII 
% 	[option, option_value]
%	'wfile':	w-file name/path
%	'stcfile': 	stc-file name/path
%	'value':	2D value matrix to be smoothed
%	'step':		the number of smoothing steps (default: 5)
%
% fhlin@feb. 26, 2004
%

wfile='';
stcfile='';
value=[];
step=5;
vertex=[];
face=[];

for i=1:length(varargin)/2
	option=varargin{i*2-1};
	option_value=varargin{i*2};
	switch option
	case 'wfile'
		wfile=option_value;
	case 'stcfile'
		stcfile=option_value;
	case 'value'
		value=option_value;
	case 'step'
		step=option_value;
	case 'vertex'
		vertex=option_value;
	case 'face'
		face=option_value;
	otherwise
		fprintf('unknown option [%s]\n',option);
		return;
	end;
end;

if(~isempty(file_asc_surf))
	if(strcmp('asc',file_asc_surf(end-2:end)))
		[nv,nf,vertex,face]=inverse_read_surf_asc(file_asc_surf);
	else
		[vertex,face]=read_surf(file_asc_surf); %freesurfer dev toolbox 061706
		vertex=vertex';
		face=face';
		nv=size(vertex,2);
		nf=size(face,2);
	end;
else
	nv=size(vertex,2);
	nf=size(face,2);
end;

fprintf('constructing connection graph\n');
%assume triangle surface!!
connection=face(1:3,:)+1; %shift zero-based dipole indices to 1-based dipole indices

d1=[connection(1,:);[1:nf];ones(1,size(connection,2))]';
d2=[connection(2,:);[1:nf];ones(1,size(connection,2))]';
d3=[connection(3,:);[1:nf];ones(1,size(connection,2))]';
A=spones(spconvert([d1;d2;d3]));
xx=sum(A,2);
B=spones(A*A'+speye(size(A,1),size(A,1)));
yy=sum(B,2);

%yy1=sum(spones(A*A'+speye(size(A,1),size(A,1))),2);
%yy2=sum(spones(A*A'*A*A'+speye(size(A,1),size(A,1))),2);
%yy3=sum(spones(A*A'*A*A'*A*A'+speye(size(A,1),size(A,1))),2);
%yy4=sum(spones(A*A'*A*A'*A*A'*A*A'+speye(size(A,1),size(A,1))),2);

neighbor=find(B(:,52363));

if(~isempty(wfile))
	[dvalue,dec_dip] = inverse_read_wfile(wfile);
	value=zeros(nv,size(dvalue,2));
	value(dec_dip+1,:)=dvalue;
end;

if(~isempty(stcfile))
	[dvalue,dec_dip]=inverse_read_stc(stcfile);
	value=zeros(nv,size(dvalue,2));
	value(dec_dip+1,:)=dvalue;
end;

if(isempty(value))
	fprintf('nothing to be smoothed!\n');
	return;
end;

if(min(size(value))==1)
	value=reshape(value,[length(value),1]);
end;

smooth_value=zeros(size(value));
dd=zeros(size(value,2),step);
for tt=1:size(value,2)
	fprintf('smoothing [%d|%d]',tt,size(value,2));
	w=value(:,tt);
	w0=w;
	%pidx=find(w0>max(w0).*0.99);
	%nidx=find(w0<min(w0).*0.5);
	pidx=find(w0>eps);
	nidx=find(w0<-eps);
	
	non_zero=find(abs(w)>eps);
	fprintf(' gridding...');
	w=griddatan(vertex(1:3,non_zero)',w(non_zero),vertex(1:3,:)','nearest');
	fprintf('done!\n');
	
	for ss=1:step
		fprintf('.');
		
		%w=(A*(mean(w(connection),1))')./xx;
		%w=mean([w,w0],2);

		paint(:,ss)=w(neighbor);
		
%         size(w)
%         size(B)
%         size(yy)
		w=B*w./yy;
		
		%w(pidx)=max([w(pidx) w0(pidx)],[],2);
		%w(nidx)=min([w(nidx) w0(nidx)],[],2);
		w(pidx)=w0(pidx);
		w(nidx)=w0(nidx);
		
		%dd(tt,ss)=sum(abs(w-w0).^2);
		
	end;
	fprintf('\n');

	smooth_value(:,tt)=fmri_scale(w,max(w0),min(w0));
end;

return;


