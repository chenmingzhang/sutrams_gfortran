function [o,o2]=readNOD(varargin)
  % readNOD reads NOD file   
  %
  % INPUT
  %   filename     -- if file is named as 'abc.nod', a filename='abc'
  %                   is required
  %   outputnumber -- number of result extracted, this is useful when
  %                   output file is huge
  %   outputstart  -- (not implimented yet) the start of the result
  %
  % OUTPUT
  % o  -- a struct the same size as the number of output.
  % o2 -- a struct extracting headers with the extraction inf
  %
  % Example:
  % [noddata,nodhead]=readNOD('project','outputnumber',3);
  %    Purpose: parsing 'project.nod' (or 'project.NOD')
  %            only the first three result gets extracted

  % a string storing the caller functions
  caller = dbstack('-completenames'); caller = caller.name;

  o2.varargin       = varargin;
  [fname, varargin] = getNext(varargin,'char','');
  % an option to see whether use inp contents to guide the reading process
  %   a hard reading process will be conducted if left empty
  [output_no,  varargin]   = getProp(varargin,'outputnumber',0);
  [output_from,  varargin] = getProp(varargin,'outputfrom',0);
  [inpObj,  varargin]      = getProp(varargin,'inpObj',[]);
  o2.output_no             = output_no;
  fn                       = fopen([fname,'.NOD']);

  if fn==-1 
    fprintf(1,'%s : Trying to open %s .nod\n',caller,fname);
    fn=fopen([fname,'.nod']);
    if fn==-1
      fprintf('%s: file nod found!!\n',caller,fname);
      o=-1;o2=-1;
      return
    end
  end
  
  o2.title1 = getNextLine(fn,'criterion',...
                         'with','keyword','## ','operation','delete');
  o2.title2 = getNextLine(fn,'criterion',...
                         'with','keyword','## ','operation','delete');

  % ---------------- Parsing the line with node, element info-----------------
  o2.MeshInfo =getNextLine(fn,'criterion','equal','keyword','## ');
  tmp=regexprep(o2.MeshInfo,{'#','(',')','\,','*','='},{'','','','','',''});
  tmp2=textscan(tmp,'%s %s %s %f %f %f %*s');
  o2.mshtyp{1}                    = tmp2{1}{1};
  o2.mshtyp{2}                    = tmp2{2}{1};
  if strcmp(o2.mshtyp{1},'2-D') && strcmp(o2.mshtyp{2},'REGULAR')
    tmp=textscan(tmp,'%s %s %s %f %f %f %*s');
  % how to realize this by one-liner
  %  [o2.mshtyp{1} o2.mshtyp{2} ] = deal(tmp{1:2}{1});
    [o2.nn1,o2.nn2,o2.nn]    = deal(tmp{4:6});
  elseif strcmp(o2.mshtyp{1},'3-D') && strcmp(o2.mshtyp{2},'BLOCKWISE')
    tmp=textscan(tmp,'%s %s %s %f %f %f %f %*s %f %*s');
    [o2.nn1,o2.nn2,o2.nn3,o2.nn,o2.ne ]    = deal(tmp{4:8});
  end
  o.nn = o2.nn;
  % ---------------- parsing the number of results    ------------------------
  tmp = getNextLine(fn,'criterion','with','keyword',...
                 '## NODEWISERESULTS','operation','delete');
  tmp           = textscan(tmp,'%f ');
  o2.ktprn      = tmp{1};  % expected no. time steps
  if output_no ~= 0;
    output_no   = min(o2.ktprn,output_no);
  else
    output_no = o2.ktprn;
  end

  % ---------------- parsing expected results    ----------------------------
  % Refering to OUTNOD.......19900
  %tmp       = getNextLine(fn,'criterion','with','keyword','##   --');
  tmp_table = textscan(fn,'##  %f %f %s %f %s %f %s %f %s %f ',o2.ktprn);
  [o2.itt,o2.tt,o2.cphorp,o2.ishorp,o2.cptorc1,o2.istorc1,o2.cptorc2,o2.istorc2,o2.cpsatu,o2.issatu]=...
    deal(tmp_table{:});

  % ---------------- Parsing simulation results -----------------------------
  fprintf(1,'%s is parsing the %g of %g outputs\n', caller,output_no,o2.ktprn);
  for n=1:output_no
    fprintf('.');
    if rem(n,50)==0; fprintf('%d\n',n);   end
    tmp       = getNextLine(fn,'criterion','with','keyword','## TIME STEP');
    if tmp  ~= -1
      tmp = regexprep(tmp,{'## TIME STEP','Duration:','sec','Time:'}...
            ,{'','','',''});
      tmp   = textscan(tmp,'%f %f %f');
      [o(n).itout,o(n).durn,o(n).tout] = deal(tmp{:});
      tmp = getNextLine(fn,'criterion','with'...
            ,'keyword','##  ','operation','delete');
%       tmp        = textscan(tmp,'%s');
%       o(n).label = tmp{1}';
      o(n).label = getOutputLabelName( tmp );
      fmt        = repmat('%f ',1, length(o(n).label));
      o(n).terms = textscan(fn,fmt,o2.nn);
    else
      fprintf(1,['WARNING FROM %s: Simulation is not completed\n %g'...
             'out of %g outputs extracted\n'],caller,n,o2.ktprn);
      return
    end % if condition
  end  % n loops
  fprintf('%s: Parsed %g of %g outputs\n', caller,output_no,o2.ktprn);
