function varargout = segmentGUI(varargin)
% SEGMENTGUI MATLAB code for segmentGUI.fig
%      SEGMENTGUI, by itself, creates a new SEGMENTGUI or raises the existing
%      singleton*.
%
%      H = SEGMENTGUI returns the handle to a new SEGMENTGUI or the handle to
%      the existing singleton*.
%
%      SEGMENTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENTGUI.M with the given input arguments.
%
%      SEGMENTGUI('Property','Value',...) creates a new SEGMENTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segmentGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segmentGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segmentGUI

% Last Modified by GUIDE v2.5 09-Sep-2020 13:27:03

%%

%global MD4D MASK

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segmentGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @segmentGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end





% --- Executes just before segmentGUI is made visible.
function segmentGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segmentGUI (see VARARGIN)

global slice cube nrows ncols nz nt MD4D MASK clim hCube 
global mouse_click
global act_write act_erase act_MARCAR xlimIMG ylimIMG
global selIZQ
global colorMat mcolors
global recienGuardado

% Choose default command line output for segmentGUI
handles.output = hObject;

% extras: INICIALIZACION ----------------
slice = 1;

[nrows,ncols,nz,nt] = size(MD4D);
xlimIMG = [0.5,ncols+0.5];
ylimIMG = [0.5,nrows+0.5];


mcolors = uint8([102,204,0; 0,0,204]); % colores marcas 1 y 2 resp (green - lightblue)
%colorMat = zeros(nrows,ncols,3); % RGB
colorMat = uint8( zeros(nrows,ncols,3) ); % RGB & UInt32 / quita gpuArray

% Crea mascaras si no existen previamente
if isempty(MASK)
  disp('- crea nueva mascara');
  MASK = uint8(zeros([nrows,ncols,nz])); %szMD4D(1:3)); % Mascara (estatica) quita gpuArray
else
  % renueva colores
  m1 = uint8(MASK(:,:,slice) == 1); % localiza indices
  m2 = uint8(MASK(:,:,slice) == 2);
    
  for k=1:3
    colorMat(:,:,k) = m1*mcolors(1,k);
    colorMat(:,:,k) = colorMat(:,:,k) + m2*mcolors(2,k); % asumo cada sub-mascara no se superpone a la otra :-)
  end
end


% Calcula imagen promediada en el tiempo una sola vez
cube = mean(single(MD4D),4);
  
% fija limites en intensidad de la imagen
c1 = min(cube(:));
c2 = max(cube(:));
clim = [c1 + 0.1*(c2-c1),c1 + 0.7*(c2-c1)];
  
% escribe 1/nz
handles.text_slice.String{1} = sprintf('1 / %d',nz);

axes(handles.ax_imagen); % gca
hCube = imshow(cube(:,:,1),clim); % muestra figura por primera vez

set(handles.mark1RButton,'Enable','off'); % por defecto no hay marcas
set(handles.mark2RButton,'Enable','off');

selIZQ = false; % parte con el boton-radio derecho
act_write = false;
act_erase = false;
act_MARCAR = false; % cuando apreta MARCAR
mouse_click = false;
recienGuardado = true;

set(0, 'units', 'normalized');
set(gcf, ...
        'windowbuttondownfcn', @mouse_down_call, ...
        'windowbuttonupfcn',   @mouse_up_call, ...
        'menubar', 'none', ...
        'pointershapecdata',   NaN(16,16), ...
        'units',   'normalized', ...
        'color',   [0 0 0], ...
        'numbertitle', 'off');

set (gcf, 'WindowButtonMotionFcn', {@mouseMove,hObject});

% Update handles structure
guidata(hObject, handles);

end



% --- Outputs from this function are returned to the command line.
function varargout = segmentGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output; % se ejecuta al inicio! << declaracion

end




%% Presion sobre boton corte +
function plusSLButton_Callback(hObject, eventdata, handles)
% hObject    handle to plusSLButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global slice nz cube hCube clim MASK colorMat
global xlimIMG ylimIMG mcolors

slice = slice + 1;

if slice > nz
  slice = nz;
end

handles.text_slice.String{1} = sprintf('%2d / %d',slice,nz);

% ---- actualiza colorMat
m1 = uint8(MASK(:,:,slice) == 1); % localiza indices
m2 = uint8(MASK(:,:,slice) == 2);
    
for k=1:3
  colorMat(:,:,k) = m1*mcolors(1,k);
  colorMat(:,:,k) = colorMat(:,:,k) + m2*mcolors(2,k); % asumo cada sub-mascara no se superpone a la otra :-)
end

xlimIMG = handles.ax_imagen.XLim; % actualizacion (fundamental)
ylimIMG = handles.ax_imagen.YLim;

% muestra imagen
axes(handles.ax_imagen); % gca

% muestra la figura y la mascara
hCube = imshow(cube(:,:,slice),clim); % <<< antes: [50 500]
hold on;
image(hCube.Parent,colorMat,'AlphaData',0.3*(MASK(:,:,slice)>0));
hold off;

set(gca,'xlim',xlimIMG); % actualiza ejes
set(gca,'ylim',ylimIMG);

% Update handles structure
guidata(hObject, handles);

end




%% Presion sobre boton corte -
function minusSLButton_Callback(hObject, eventdata, handles)
% hObject    handle to minusSLButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global slice nz cube hCube clim MASK colorMat
global xlimIMG ylimIMG mcolors

slice = slice - 1;

if slice < 1
  slice = 1;
end

handles.text_slice.String{1} = sprintf('%2d / %d',slice,nz);

% ---- actualiza colorMat
m1 = uint8(MASK(:,:,slice) == 1); % localiza indices
m2 = uint8(MASK(:,:,slice) == 2);
    
for k=1:3
  colorMat(:,:,k) = m1*mcolors(1,k);
  colorMat(:,:,k) = colorMat(:,:,k) + m2*mcolors(2,k); % asumo cada sub-mascara no se superpone a la otra :-)
end

% muestra la figura y la mascara
xlimIMG = handles.ax_imagen.XLim; % actualizacion (fundamental)
ylimIMG = handles.ax_imagen.YLim;

axes(handles.ax_imagen); % gca
hCube = imshow(cube(:,:,slice),clim); % <<< antes: [50 500]
hold on;
image(hCube.Parent,colorMat,'AlphaData',0.3*(MASK(:,:,slice)>0));
hold off;

set(gca,'xlim',xlimIMG); % actualiza ejes
set(gca,'ylim',ylimIMG);

% Update handles structure
guidata(hObject, handles);

end




%% Creacion de la imagen
% --- Executes during object creation, after setting all properties.
function ax_imagen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax_imagen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax_imagen
% >>> Nota: esta funcion se ejecuta antes de OpeningFcn
end





%% PANEL MARCAS
% --- Executes when selected object is changed in markGroup.
function markGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in markGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global selIZQ

switch handles.markGroup.SelectedObject.String
  case 'Izquierda'
    selIZQ = true;

  case 'Derecha'
    selIZQ = false;
    
end

end




%% PANEL ACCIONES
% --- Executes when selected object is changed in accionesGroup.
function accionesGroup_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in accionesGroup 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global act_MARCAR MASK

switch handles.accionesGroup.SelectedObject.String
  case 'Zoom'
    %disp('Zoom');
    act_MARCAR = false;
    
  case 'Pan'
    %disp('Pan');
    act_MARCAR = false;
        
  case 'Marcar'
    %disp('Mark');
    act_MARCAR = true;
    
  case 'Salir'
    disp('Fin de la segmentacion');
    varargout{2} = MASK; % a ver si resulta
    
    return;
    
end

end


% --- Executes on button press in zoomButton.
function zoomButton_Callback(hObject, eventdata, handles)
% hObject    handle to zoomButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomButton

% mientras estï¿½este boton activado (los otros desactivados)
% pan se fija
zoom on;
pan off;

axes(handles.ax_imagen); % gca

% desativa ventana de 'Marcas'
set(handles.mark1RButton,'Enable','off');
set(handles.mark2RButton,'Enable','off');

% Update handles structure
guidata(hObject, handles);

end


% --- Executes on button press in panButton.
function panButton_Callback(hObject, eventdata, handles)
% hObject    handle to panButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of panButton

% mientras estï¿½este boton activado (los otros desactivados)
% zoom se fija
zoom off;
pan on;

axes(handles.ax_imagen); % gca

% desativa ventana de 'Marcas'
set(handles.mark1RButton,'Enable','off');
set(handles.mark2RButton,'Enable','off');

%set(gca,'xlim',xlimIMG); % fundamental ?
%set(gca,'ylim',ylimIMG);

% Update handles structure
guidata(hObject, handles);

end


% --- Executes on button press in markButton.
function markButton_Callback(hObject, eventdata, handles)
% hObject    handle to markButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of markButton

global xlimIMG ylimIMG

% los efectos de movimiento se fijan
zoom off;
pan off;

%axes(handles.ax_imagen); % gca

%disp(' --- markButton');
xlimIMG = handles.ax_imagen.XLim;  % importante
ylimIMG = handles.ax_imagen.YLim; 

axes(handles.ax_imagen); % gca

% activa ventana de 'Marcas'
set(handles.mark1RButton,'Enable','on');
set(handles.mark2RButton,'Enable','on');

% Update handles structure
guidata(hObject, handles);

end


% --- Executes on button press in mark1RButton.
function mark1RButton_Callback(hObject, eventdata,handles)
% hObject    handle to mark1RButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mark1RButton

global esIZQ

esIZQ = true;

end



% --- Executes on button press in mark2RButton.
function mark2RButton_Callback(hObject, eventdata, handles)
% hObject    handle to mark2RButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mark2RButton

global esIZQ

esIZQ = false;

end



%% Clicks del mouse

% ----
% ---- ACTUA CUANDO SE APRETA UN BOTON DEL MOUSE
function mouse_down_call(hObject,eventdata)
global act_write act_erase act_MARCAR
global mouse_click recienGuardado


mouse_click = true;

button_show = get(gcf,'SelectionType');

switch button_show
  case 'normal'
    %disp('dibujar');
    act_write = true && act_MARCAR;
    act_erase = false;
    recienGuardado = false;
    
  case 'alt'
    %disp('borrar');
    act_write = false;
    act_erase = true && act_MARCAR;
    recienGuardado = false;
    
  case 'extend'
    % cambiamos apretando la rueda del mouse
    act_write = false;
    act_erase = false;
    
  case 'open'
    % repite la ultima accion
    %[act_write, act_erase]
    
end

end

% ------
% ---- ACTUA CUANDO SE SUELTA AL BOTON DEL MOUSE
%      funciona a consecuencia de mouse_down_call
function mouse_up_call(hObject,eventdata) 
global mouse_click act_write act_erase 

mouse_click = false;

% cuando se suelta el boton, deja de escribir o borrar
act_write = false;
act_erase = false;

end


% -
function mouseMove(hObject,eventdata,hFigure)
%function timerCallback(source,event,hFigure)

global selIZQ act_write act_erase MASK slice cube clim
global ncols nrows colorMat mcolors act_MARCAR xlimIMG ylimIMG
global mouse_click

if act_MARCAR && mouse_click
  handles = guidata(hFigure);
  
  if act_write
    val = 1.0*(selIZQ) + 2.0*(~selIZQ);
    
    C = get(handles.ax_imagen,'CurrentPoint');
    RX=ceil( C(1,1) );
    RY=ceil( C(1,2) );
    
    if RX<2
      RX = 2;
    end
    if RY<2
      RY = 2;
    end
    if RX>ncols-1
      RX = ncols-1;
    end
    if RY>nrows-1
      RY = nrows-1;
    end
    
    MASK(RY-1:RY+1,RX-1:RX+1,slice) = uint8(val);
    
    cc(1,1,1:3) = uint8(selIZQ)*mcolors(1,:) + uint8(~selIZQ)*mcolors(2,:);
    colorMat(RY-1:RY+1,RX-1:RX+1,:) = uint8( repmat(cc,3,3,1) );
    
  end
  
  if act_erase
    C = get(handles.ax_imagen,'CurrentPoint');
    RX=ceil( C(1,1) );
    RY=ceil( C(1,2) );
    
    if RX<1
      RX = 1;
    end
    if RY<1
      RY = 1;
    end
    if RX>ncols
      RX = ncols;
    end
    if RY>nrows
      RY = nrows;
    end
    
    MASK(RY-1:RY+1,RX-1:RX+1,slice) = uint8(0); % puede borrar seleccion 1 o 2 (izq o der)
    colorMat(RY-1:RY+1,RX-1:RX+1,:) = uint8( zeros(3,3,3) );
  end
  
  
  % --- dibuja
  % Atencion: "handles" aqui es una copia de "handles" original
  % muestra la figura
  hCube = imshow(cube(:,:,slice),clim); 
  
  hold on;
  image(hCube.Parent,colorMat,'AlphaData',0.3*(MASK(:,:,slice)>0));
  hold off;
  
  set(gca,'xlim',xlimIMG);
  set(gca,'ylim',ylimIMG);
    
end % if marcar



end % function mouseMove


% ---

% --- Executes on button press in borrarMaskButton.
function borrarMaskButton_Callback(hObject, eventdata, handles)
% hObject    handle to borrarMaskButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global MASK colorMat recienGuardado nrows ncols nz cube slice clim

% reinicia mascaras
MASK     = uint8( zeros([nrows,ncols,nz]) ); % quita gpuArray
colorMat = uint8( zeros(nrows,ncols,3) ); % RGB & UInt8 / quita gpuArray

recienGuardado = false;

% muestra actualizacion de la figura
hCube = imshow(cube(:,:,slice),clim);
hold on;
image(hCube.Parent,colorMat,'AlphaData',0.3*(MASK(:,:,slice)>0));
hold off;
set(gca,'xlim',xlimIMG);
set(gca,'ylim',ylimIMG);

end



%% ----------------- CIERRE -----------------

%{
% --- Executes on button press in slicePMButton.
function slicePMButton_Callback(hObject, eventdata, handles)
% hObject    handle to slicePMButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slicePMButton

% Aqui cuando se apreta "salir"

disp('exit_button_callback');
closereq(); % esto es lo ultimo que hace, y borra el handle (deleted figure)
end
%}


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% aqui cuando se cierra la ventana desde "X"

disp('cierre forzoso (close request)');

% Hint: delete(hObject) closes the figure
delete(hObject);

end


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global MASK recienGuardado mcolors colorMat slice xlimIMG ylimIMG cube clim

[file,path,filteridx] = uigetfile({'*.mat','Archivos MAT'},'Seleccione una segmentación que sobreescribe a la actual.');

if ~isequal(file,0) || ~isequal(path,0)
  myfile = fullfile(path,file);
  if isfile(myfile)
    aa = load(myfile,'MASK');
    %MASK = gpuArray(aa.MASK); % Ojo con GPU
    MASK = aa.MASK; % Ojo con GPU
    
    recienGuardado = true;
    
    % ---- actualiza colorMat
    m1 = uint8(MASK(:,:,slice) == 1); % localiza indices
    m2 = uint8(MASK(:,:,slice) == 2);
    
    for k=1:3
      colorMat(:,:,k) = m1*mcolors(1,k);
      colorMat(:,:,k) = colorMat(:,:,k) + m2*mcolors(2,k); % asumo cada sub-mascara no se superpone a la otra :-)
    end
        
    % --- muestra la figura y la mascara
    xlimIMG = handles.ax_imagen.XLim; % actualizacion (fundamental)
    ylimIMG = handles.ax_imagen.YLim;
    
    axes(handles.ax_imagen); % gca
    hCube = imshow(cube(:,:,slice),clim); % <<< antes: [50 500]
    hold on;
    image(hCube.Parent,colorMat,'AlphaData',0.3*(MASK(:,:,slice)>0));
    hold off;

    set(gca,'xlim',xlimIMG); % actualiza ejes
    set(gca,'ylim',ylimIMG);
    
    % listo.
    
  else
    opts = struct('WindowStyle','modal', ...
    'Interpreter','tex');
  warndlg(['No se encuentra el archivo "',myfile,'"'],...
    'segmentGUI', opts);
  end
end
%else
%  disp('cancela carga de segmentacion');
end




% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

grabarDatos();

end

% --- Executes on button press in exitButton.
function exitButton_Callback(hObject, eventdata, handles)
% hObject    handle to exitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global recienGuardado

if ~recienGuardado
  
  answer = questdlg('¿Desea guardar a la segmentación y a la imagen antes de salir?', ...
    'Salir de segmentGUI' , 'Si' , 'No' , 'No');
  
  if strcmp(answer,'Si')
    %disp(' grabando ');
    grabarDatos();
  end
  
end

closereq(); % esto es lo ultimo que hace, y borra el handle (deleted figure)

end


% --- GRABAR ---
function grabarDatos(v)

global MASK MD4D time dataDICOM

[file,ppath] = uiputfile({'*.mat','Archivos MAT'},'Elija donde guardar la segmentación actual');

if ~isequal(file,0) || ~isequal(ppath,0)
  
  %MASK = gather(MASK); si hubiera gpu en todas partes
  %MD4D = gather(MD4D);
    
  save(fullfile(ppath,file),'MASK','MD4D','time','dataDICOM');
  
  %MASK = gpuArray(MASK); % se devuelve a GPU
  %MD4D = gpuArray(MD4D);
  
  
else
  opts = struct('WindowStyle','modal', ...
    'Interpreter','tex');
  warndlg('\color{blue} No se guardaron datos',...
    'segmentGUI', opts);
end
end



%% ========================================================================
%  APENDICE: Funciones o partes de codigo interesante pero quedaron sin uso


%{
% ---  ESTATICO --- sirve mientras el mouse no se mueve, solo si pincha
******
% --- en OpenFnc ...

% 'closerequestfcn',     @close_req_call,  'pointer', 'custom', ...

% - recommended by
% https://www.mathworks.com/matlabcentral/answers/275722-windowbuttondownfcn-cant-detect-handles-data-from-other-callbacks
%set(hObject, ...
%    'WindowButtonDownFcn', @mouse_down_call, ...
%    'WindowButtonUpFcn', @mouse_up_call);


%{
con boton apretado, pero no al mover
% --- mientras el boton del mouse esta apretado, llama a timer:
handles.timer = timer; 

set(handles.timer, 'ExecutionMode', 'fixedrate');
set(handles.timer, 'Period', .1);
set(handles.timer, 'TasksToExecute', inf);
set(handles.timer, 'TimerFcn', {@timerCallback,hObject});
%}

%{
myTimer = timer('Name','MyMouseButtonTimer',    ...
  'Period',0.3,                   ...
  'StartDelay',0.001,             ...
  'TasksToExecute',inf,           ...
  'ExecutionMode','fixedSpacing', ...
  'TimerFcn',@timerCallback); % pasa 3er arg a timerCallback
%'TimerFcn',{@timerCallback,handles}); % pasa 3er arg a timerCallback
%}



% ---
en

function timerCallback(hObject,eventdata,hFigure)
%function timerCallback(source,event,hFigure)

global selIZQ act_write act_erase MASK slice cube clim
global ncols nrows

disp('  timecallback ---- ');

%{
% coordenadas normalizadas en la figura (complicado)
mouse_pos  = get(0, 'pointerlocation');
figure_pos = get(gcf, 'position');
pos_x = (mouse_pos(1) - figure_pos(1))/(figure_pos(3));        
pos_y = (mouse_pos(2) - figure_pos(2))/(figure_pos(4)); 
%}

handles = guidata(hFigure);

%[xx,yy] = ginput(1); % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%fprintf('GINPUT x=%g  y=%g\n',xx,yy); % <<<<<<<< TEST

if act_write
  val = 1.0*(selIZQ) + 2.0*(~selIZQ);

  %C = get(handles.ax_imagen,'CurrentPoint');
  C = get(gca,'CurrentPoint');
  RX=ceil( C(1,1) );
  RY=ceil( C(1,2) );
  
  if RX<1 
    RX = 1;
  end
  if RY<1
    RY = 1;
  end
  if RX>ncols 
    RX = ncols;
  end
  if RY>nrows
    RY = nrows;
  end
  
  fprintf('pinta x=%g  y=%g\n',RX,RY); % <<<<<<<< TEST

  MASK(RY,RX,slice) = val; % puede borrar seleccion 1 o 2 (izq o der)
end

if act_erase
  C = get(handles.ax_imagen,'CurrentPoint');
  RX=ceil( C(1,1) );
  RY=ceil( C(1,2) );
  
  if RX<1
    RX = 1;
  end
  if RY<1
    RY = 1;
  end
  if RX>ncols
    RX = ncols;
  end
  if RY>nrows
    RY = nrows;
  end
  
  fprintf('borra x=%g  y=%g\n',RX,RY); % <<<<<<<< TEST

  MASK(RY,RX,slice) = 0; % puede borrar seleccion 1 o 2 (izq o der)
end

%{
fprintf('      mouse button is still down!\n'); %% <<<<<<
if selIZQ disp('[IZQ]'); else disp('[DER]'); end  %% <<<<

if act_write disp('dibujar'); end
if act_erase disp('borrar'); end
%}

% --- dibuja
%imshow(cube(:,:,slice),clim); % <<< antes: [50 500]
  
%hold on;
%image(MASK(:,:,slice),'AlphaData',0.3*(MASK(:,:,slice)~=0));
%hold off;


%guidata(hObject, handles);
guidata(hFigure, handles);

end
%}
