function segmenta()
% SEGMENTA - Read/Modifies or Create the mask for a given cubic image using manual
%            segmentations.
%
% Last update: Sep 11, 2020
% J. Mura - V. Canales - G. Vasquez - V. Silva

global MD4D MASK time dataDICOM % no encontre una solucion mejor para pasar datos a GUI :-/

DEBUG_ = true; % modo debug ?



%% Seleccion de datos de entrada: DICOM o MAT

[Files,Path] = uigetfile({'*.dcm;*.mat',...
  'Archivos admisibles (*.dcm,*.mat)';'*.dcm','DICOM (*.dcm)';'*.mat','Archivos MAT (*.mat)'},...
  'Seleccione todos los archivos DICOM o una segmentacion previa',...
  'Multiselect','on');

if isequal(Files,0)
  % --- no leido ---
  disp('Lectura cancelada. Adios!');
  return;
else

  if iscell(Files)
    numfiles = length(Files); % cuando son varios vienen en celdas
    file1 = Files{1};
  else
    numfiles = 1;
    file1 = Files;
  end
  
  if isequal(file1(end-2:end),'dcm');
    
    if DEBUG_
      fid = fopen('log.txt','w');
    end
    
    % --- DICOM ---
    fprintf('lectura de %d archivos en formato DICOM ...',numfiles);
    
    szIMGdyn = size(dicomread(fullfile(Path,Files{1})));
    
    dataDICOM = dicominfo(fullfile(Path,Files{1}));
    ntp = dataDICOM.NumberOfTemporalPositions;
    
    num_slices = numfiles/ntp; % cantidad de slices en eje Z
    
    if (num_slices ~= round(num_slices))
      errordlg(sprintf('Cantidad de archivos (%d) no es multiplo de %d',numfile,ntp));
    end
    
    all_img = zeros([szIMGdyn,num_slices]); % pre-alocacion

    Data_ = zeros(numfiles,2); % guarda todo para ordenar despues
    
    for i=1:numfiles
      all_img(:,:,i) = dicomread(fullfile(Path,Files{i})); % dicomread( ...,'frames','all')??
      
      bb = dicominfo(fullfile(Path,Files{i}));
      
      Data_(i,1) = bb.SliceLocation; % position
      Data_(i,2) = str2double(bb.AcquisitionTime); % time
      
      if DEBUG_
        fprintf(fid,'file %3d: %s | slice: %3.3f | time: %3.3f \n',...
          i,Files{i},Data_(i,1),Data_(i,2));
      end
      
    end

    % calculo enredado pero es correcto.
    for i=1:ntp % Se adquiere el vector con los tiempos generalizados
      aa = dicominfo(fullfile(Path,Files{i}));
      a = str2double(aa.AcquisitionTime);
      
      if DEBUG_
        fprintf(fid,'time[%2d] = (num) %7.3f  / (str) %s \n',...
          i,a,aa.AcquisitionTime);
      end
      
      horasec = floor(a/10000)*3600;
      minutosec = floor((a-floor(a/10000)*10000)/100)*60;
      sec = a-floor((a/100))*100;
      time(i) = horasec + minutosec + sec; % lo unico que hace esto es pasar todo a segundos
    end
    time = time - time(1); % lleva al instante inicial a 0
    
    
    %{  
    %  0,12.7,25,57,...,402.86,415.64 !!! demasiado largo?? explicacion??
    Time_ = unique(Data_(:,2));
    time = Time_ - min(Time_); % fija t=0
    if DEBUG_
      for ii=1:length(time)
        fprintf(fid,'time[%2d]: %3.2f = %3.2f\n',ii,Time_(ii),time(ii));
      end
    end
    %}
    
    Sl_ = unique(Data_(:,1));
    slices = Sl_ - min(Sl_);
    if DEBUG_
      for ii=1:length(slices)
        fprintf(fid,'pos[%2d]: %3.2f = %3.2f\n',ii,Sl_(ii),slices(ii));
      end
    end
    
    
    if DEBUG_
      fclose(fid);
    end
    
    nz = numfiles/ntp; % << OJO: siempre da entero
    
    MD4D = uint16( zeros(dataDICOM.Rows,dataDICOM.Columns,nz,ntp)); % quitamos 'gpuArray'

    for z = 1:nz
      for t = 1:ntp
        pos = ntp*(z-1) + t; % OJO! aqui habia un error antes (12 en lugar de ntp)?
        MD4D(:,:,z,t) = uint16( all_img(:,:,pos) );
      end
    end
    
    % no hay mascara
    MASK = []; % empty
    
    fprintf(' listo.\n');
    
  elseif numfiles == 1
    % --- MAT ---
    fprintf('lectura de una segmentacion en formato MAT ...');
    data = load(fullfile(Path,Files));
    
    MD4D = data.MD4D;
    time = data.time;
    MASK = uint8(data.MASK);
    dataDICOM = data.dataDICOM;
    
    fprintf(' listo.\n');
  else
    % --- demasiados MAT ---
    errordlg('No se puede leer mas de 1 archivo *.mat al mismo tiempo.');
    return;
  end
end
    

%% Procede con la segmentacion
segmentGUI();

%{
ojo: 
voxelsize = [ dataDICOM.PixelSpacing , dataDICOM.SliceThickness] 
en mm

dataDICOM.SpacingBetweenSlices ???
%}


end % function