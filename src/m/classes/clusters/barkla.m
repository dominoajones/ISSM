%barkla cluster class definition
%
%   Usage:
%      cluster=hpc();
%      cluster=hpc('np',3);
%      cluster=hpc('np',3,'login','username');
%
%   PARTITION | NODES | MEMORY | TIMELIMIT
%     nodes*     53     380000   3*24*60*60
%     long       8      380000   7*24*60*60
%     himem      2      1100000  6*24*60*60
%     phi        4      192000   3*24*60*60
%     gpu        5      380000   3*24*60*60

classdef barkla
	properties (SetAccess=public)
		% {{{
		name='barkla6.liv.ac.uk'
		login='domjon';
		numtasks       = 1;
		cpuspertask    = 8;
		cpuspernode=4; 
		port=0;
		queue='long';
		codepath='/mnt/data1/users/domjon/ISSM/bin/';
		executionpath='/mnt/data1/users/domjon/ISSM/execution';
		interactive=0;
        time= '7-00:00:00';
		memory= 2;
        email= 'domino.jones@liverpool.ac.uk';
		mailtype= 'ALL';
	end
	%}}}
	methods
		function cluster=barkla(varargin) % {{{

			%initialize cluster using default settings if provided
			if (exist('barkla_settings')==2), barkla_settings; end

			%use provided options to change fields
			cluster=AssignObjectFields(pairoptions(varargin{:}),cluster);
		end
		%}}}
				 function disp(cluster) % {{{
			 %  display the object
			 disp(sprintf('class ''%s'' object ''%s'' = ',class(cluster),inputname(1)));
			 disp(sprintf('    name: %s',cluster.name));
			 disp(sprintf('    login: %s',cluster.login));
			 disp(sprintf('    port: %i',cluster.port));
			 disp(sprintf('    numtasks: %i',cluster.numtasks));
			 disp(sprintf('    cpuspertask: %i',cluster.cpuspertask));
			 disp(sprintf('    codepath: %s',cluster.codepath));
			 disp(sprintf('    executionpath: %s',cluster.executionpath));
			 disp(sprintf('    time: %i',cluster.time));
			 disp(sprintf('    memory: %i',cluster.memory));
			 disp(sprintf('    email: %s', cluster.email));
			 disp(sprintf('    mailtype: %s', cluster.mailtype));
			 
		 end
		 %}}}
		 function numprocs=np(cluster) % {{{
			 %compute number of processors
			 numprocs=cluster.numtasks*cluster.cpuspertask;
		 end
		 %}}}
		 function md = checkconsistency(cluster,md,solution,analyses) % {{{
			 if isempty(cluster.name), md = checkmessage(md,'name empty'); end
			 if isempty(cluster.login), md = checkmessage(md,'login empty'); end
			 if ~(cluster.numtasks > 0), md = checkmessage(md,'numtasks must be > 0'); end
			 if ~(cluster.cpuspertask > 0), md = checkmessage(md,'cpuspertask must be > 0'); end
			 if ~(cluster.port >= 0), md = checkmessage(md,'port must be >= 0'); end
			 if isempty(cluster.email), md = checkmessage(md,'email empty'); end
			 if isempty(cluster.mailtype), md = checkmessage(md,'mailtype empty'); end
			 if isempty(cluster.codepath), md = checkmessage(md,'codepath empty'); end
			 if isempty(cluster.executionpath), md = checkmessage(md,'executionpath empty'); end
			 if ~(cluster.time > 0), md = checkmessage(md,'time must be > 0'); end
			 if ~(cluster.memory > 0), md = checkmessage(md,'memory must be > 0'); end
		 end
		 %}}}
		 function BuildKrigingQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{
			 error('not implemented yet');
		 end
		 %}}}
		 function BuildQueueScript(cluster,dirname,modelname,solution,io_gather,isvalgrind,isgprof,isdakota,isoceancoupling) % {{{

			 if(isvalgrind), disp('valgrind not supported by cluster, ignoring...'); end
			 if(isgprof),    disp('gprof not supported by cluster, ignoring...'); end

			 %write queuing script 
			 fid=fopen([modelname '.queue'],'w');
			 fprintf(fid,'#!/bin/bash\n');
			 fprintf(fid,'#SBATCH --job-name=%s\n',modelname);
			 fprintf(fid,'#SBATCH --ntasks=%i  \n',cluster.numtasks);
			 fprintf(fid,'#SBATCH --cpus-per-task=%i\n',cluster.cpuspertask);
			 fprintf(fid,'#SBATCH --time=3-00:00:00\n'); 
			 fprintf(fid,'#SBATCH --mem-per-cpu=%igb\n',cluster.memory); %memory in in gigabytes
			 fprintf(fid,'#SBATCH --mail-user=%s\n',cluster.email); %email
			 fprintf(fid,'#SBATCH --mail-type=%s\n',cluster.mailtype); 
			 fprintf(fid,'#SBATCH --output=%s.outlog \n',modelname);
			 fprintf(fid,'#SBATCH --error=%s.errlog \n\n',modelname);
			 fprintf(fid,'export ISSM_DIR="%s/"\n',cluster.codepath); %FIXME
			 fprintf(fid,'cd %s/%s\n\n',cluster.executionpath,dirname);
			 fprintf(fid,'srun %s/issm.exe %s %s %s\n',cluster.codepath,solution,[cluster.executionpath '/' dirname],modelname);
			 if ~io_gather, %concatenate the output files:
				 fprintf(fid,'cat %s.outbin.* > %s.outbin',modelname,modelname);
			 end
			 fclose(fid);
		 end %}}}
		 function UploadQueueJob(cluster,modelname,dirname,filelist)% {{{
             % filelist
			 %compress the files into one zip.
			 compressstring=['tar -zcf ' dirname '.tar.gz '];
			 for i=1:numel(filelist),
				 compressstring = [compressstring ' ' filelist{i}];
			 end
			 system(compressstring);

			 disp('uploading input file and queuing script');
			 issmscpout(cluster.name,cluster.executionpath,cluster.login,cluster.port,{[dirname '.tar.gz']});

		 end %}}}
		 function LaunchQueueJob(cluster,modelname,dirname,filelist,restart,batch)% {{{

			 disp('launching solution sequence on remote cluster');
			 if ~isempty(restart)
				 launchcommand=['cd ' cluster.executionpath ' && cd ' dirname ' && hostname && sbatch ' modelname '.queue '];
			 else
				 launchcommand=['cd ' cluster.executionpath ' && rm -rf ./' dirname ' && mkdir ' dirname ...
					 ' && cd ' dirname ' && mv ../' dirname '.tar.gz ./ && tar -zxf ' dirname '.tar.gz  && hostname && sbatch ' modelname '.queue '];
			 end
			 issmssh(cluster.name,cluster.login,cluster.port,launchcommand)
		 end %}}}
		 function Download(cluster,dirname,filelist)% {{{

             %copy files from cluster to current directory
			 directory=[cluster.executionpath '/' dirname '/'];
			 issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
             
            % 
            % directory = sprintf('%s/%s/', cluster.executionpath, dirname);
            % % fileliststr = ['{' strjoin(arrayfun(@(x) char(x), filelist, 'UniformOutput', false), ',') '}'];
            % fileliststr = 'inverted.outbin'
            % fullpathtofiles = fullfile(directory, fileliststr);
            % fullpathtofilesunix = fullpathtofiles;
            % fullpathtofilesunix(strfind(fullpathtofilesunix,'\'))='/';
            % downloadcommand = sprintf('scp %s@%s:%s %s/.', cluster.login, cluster.name, fullpathtofilesunix, pwd);
            % % downloadcommand = sprintf('scp %s@%s:%s %s/.', cluster.login, cluster.name, fullfile(directory, fileliststr), pwd)
            % status = system(downloadcommand);
            % 
			 % %copy files from cluster to current directory
			 % directory=[cluster.executionpath '/' dirname '/'];
			 % issmscpin(cluster.name,cluster.login,cluster.port,directory,filelist);
            % if numel(filelist)==1,
	        %     filestring=filelist{1};
            % else
	        %     filestring='\{';
	        % for i=1:numel(filelist)-1,
		    %     filestring=[filestring filelist{i} ','];
	        % end
	        %     filestring=[filestring filelist{end} '\}'];
            % end
            % 
            % % downloadstring=['scp -i ' cluster.idfile ' ' cluster.login '@' cluster.name ':' directory '/' filestring ' ./'];
            % downloadstring=['scp -i ' cluster.login '@' cluster.name ':' directory '/' filestring ' ./'];
            % [status,result]=system(downloadstring);
            % if status, 
	        %     error(['cluster.Download error message: ' status]);
            % end
            % 
            % %check scp worked
            % for i=1:numel(filelist),
	        %     if ~exist(['./' filelist{i}]),
		    %         warning(['cluster.Download error message: could not scp ' filelist{i}]);
	        %     end
            % end

		 end %}}}
	end
end
