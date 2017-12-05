function varargout = SNR(varargin)
% SNR MATLAB code for SNR.fig
%      SNR, by itself, creates a new SNR or raises the existing
%      singleton*.
%
%      H = SNR returns the handle to a new SNR or the handle to
%      the existing singleton*.
%
%      SNR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SNR.M with the given input arguments.
%
%      SNR('Property','Value',...) creates a new SNR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SNR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SNR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SNR

% Last Modified by GUIDE v2.5 13-Feb-2017 13:25:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SNR_OpeningFcn, ...
                   'gui_OutputFcn',  @SNR_OutputFcn, ...
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


% --- Executes just before SNR is made visible.
function SNR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SNR (see VARARGIN)

% Choose default command line output for SNR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SNR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SNR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDA
global dwi mask_noise mask_signal


[filename, pathname] = uigetfile({'*.nhdr;*.nrrd','NRRD file'},...
          'Select NRRD file for Quality Control');
    
dwi= loadNrrdStructure([pathname,filename]); % this part is to load the DWI image selected by the user


tmp=([pathname, '/tmp']);
system(['mkdir ' tmp]);
outputMaskName=([tmp '/dwi_Mask_Noise.nii.gz'])


%%%%%%%%%%%%%%%%%%% THIS PART IS TO CREATE MASK FOR MONKEY DATA

%%%% this part is to create a mask for noise component
% 
%  system(['/rfanfs/pnl-zorro/projects/Lyall_R03/Slicer-build2/Slicer-build/Slicer  --launch DiffusionWeightedVolumeMasking ' pathname '/' filename ' ' tmp '/dwi_b0.nrrd ' tmp '/dwi_mask_noise.nrrd --baselineBValueThreshold 1000 --removeislands']);
%  system(['ConvertBetweenFileFormats ' tmp '/dwi_b0.nrrd ' tmp '/dwi_bo.nii.gz']);
%  mask_dwi=strrep('dwi_mask_noise.nrrd','nrrd','nii.gz');
%  system(['ConvertBetweenFileFormats ' tmp '/dwi_mask_noise.nrrd ' tmp '/' mask_dwi]);
%  system(['fslmaths ' tmp '/dwi_mask_noise.nii.gz -dilM ' tmp '/dwi_mask_noise.nii.gz']);
%  system(['fslmaths '  tmp '/dwi_mask_noise.nii.gz -fillh  '  outputMaskName ]);
%  FinalMaskNoise=strrep(outputMaskName,'.nii.gz', '.nhdr'); %% this is the final name of the file that gets uploaded for the noise component
%  system(['ConvertBetweenFileFormats '  outputMaskName ' '  FinalMaskNoise]);
%  
%%% I have changed this part to fit for Stiffy

% system(['/rfanfs/pnl-zorro/software/Slicer-4.5.0-1-linux-amd64/Slicer  --launch DiffusionWeightedVolumeMasking ' pathname '/' filename ' ' tmp '/dwi_b0.nrrd ' tmp '/dwi_mask.nrrd --otsuomegathreshold 0']);
% system(['ConvertBetweenFileFormats ' tmp '/dwi_b0.nrrd ' tmp '/dwi_bo.nii.gz']);
% mask_dwi=strrep('dwi_mask.nrrd','nrrd','nii.gz');
% system(['ConvertBetweenFileFormats ' tmp '/dwi_mask.nrrd ' tmp '/' mask_dwi]);
% system(['bet ' tmp '/dwi_bo.nii.gz ' tmp '/dwi-bet -t -s -m -A ']);
% system(['fslmaths ' tmp '/dwi-bet_skull_mask.nii.gz -dilM ' tmp '/dwi-bet_skull_mask_1.nii.gz']);
% system(['fslmaths ' tmp '/dwi-bet_skull_mask_1.nii.gz -dilM ' tmp '/dwi-bet_skull_mask_2.nii.gz']);
% system(['fslmaths ' tmp '/dwi-bet_skull_mask_2.nii.gz -dilM ' tmp '/dwi-bet_skull_mask_3.nii.gz']);
% system(['fslmaths ' tmp '/dwi-bet_skull_mask_3.nii.gz -dilM ' tmp '/dwi-bet_skull_mask_4.nii.gz']);
% system(['fslmaths ' tmp '/dwi-bet_skull_mask_4.nii.gz -add  ' tmp '/' mask_dwi ' ' outputMaskName]);
% system(['fslmaths '  outputMaskName ' -bin  -fillh '   outputMaskName ]);
% FinalMaskNoise=strrep(outputMaskName,'.nii.gz', '.nhdr');
% system(['ConvertBetweenFileFormats '  outputMaskName ' '  FinalMaskNoise]);

%%%% THIS PART IS TO CREATE SIGNAL MASK FOR THE MONKEY DATA


%%%%%%%%%%%%%%%%%% the part above is for the monkey data set DO not edit

% %%%%%%%%%%%%%%%% THIS PART IS FOR THE HUMAN DATA

%%%% this part is to create a mask for noise component

system(['/rfanfs/pnl-zorro/software/Slicer-4.5.0-1-linux-amd64/Slicer  --launch DiffusionWeightedVolumeMasking ' pathname '/' filename ' ' tmp '/dwi_b0.nrrd ' tmp '/dwi_mask.nrrd --otsuomegathreshold 0']);
%system(['/rfanfs/pnl-zorro/projects/Lyall_R03/Slicer-build2/Slicer-build/Slicer  --launch DiffusionWeightedVolumeMasking ' pathname '/' filename ' ' tmp '/dwi_b0.nrrd ' tmp '/dwi_mask.nrrd --baselineBValueThreshold 0 --removeislands']);
system(['ConvertBetweenFileFormats ' tmp '/dwi_b0.nrrd ' tmp '/volume.nii']);
mask_dwi=strrep('dwi_mask.nrrd','nrrd','nii.gz');
system(['ConvertBetweenFileFormats ' tmp '/dwi_mask.nrrd ' tmp '/' mask_dwi]);
system(['bet ' tmp '/volume.nii ' tmp '/dwi-bet -t -s -m -A ']);
system(['fast -v -t 2 ' tmp '/seg '  tmp '/volume.nii']);
system(['fslmaths ' tmp '/volume_pve_0.nii.gz -mas ' tmp '/volume_pve_0.nii.gz']);
system(['fslmaths ' tmp '/volume_pve_1.nii.gz -mas ' tmp '/volume_pve_1.nii.gz']);
system(['fslmaths ' tmp '/volume_pve_0.nii.gz -add ' tmp '/volume_pve_1.nii.gz ' tmp '/added_noise.nii.gz']);
system(['fslmaths ' tmp '/dwi-bet_outskin_mask.nii.gz -dilM ' tmp '/dwi-bet_outskin_mask.nii.gz']);
system(['fslmaths ' tmp '/dwi-bet_outskin_mask.nii.gz -dilM ' tmp '/dwi-bet_outskin_mask.nii.gz']);
system(['fslmaths ' tmp '/dwi-bet_outskin_mask.nii.gz -dilM ' tmp '/dwi-bet_outskin_mask.nii.gz']);
system(['fslmaths ' tmp '/added_noise.nii.gz -add  ' tmp '/dwi-bet_outskin_mask.nii.gz  ' tmp '/dwi-Nosie_mask.nii.gz']);
system(['fslmaths ' tmp '/dwi-Nosie_mask.nii.gz  -add  ' tmp '/' mask_dwi ' ' outputMaskName]);
system(['fslmaths '  outputMaskName ' -bin  -fillh '   outputMaskName ]);
FinalMaskNoise=strrep(outputMaskName,'.nii.gz', '.nhdr');
system(['ConvertBetweenFileFormats '  outputMaskName ' '  FinalMaskNoise]);

%%% this part is to create mask for the signal component 

system(['bet ' tmp '/volume.nii ' tmp '/dwi-signal-bet -f 0.3 -m -R']);
system(['fslmaths ' tmp '/dwi-signal-bet_mask.nii.gz -ero ' tmp '/dwi-signal-bet_mask.nii.gz']);
FinalMaskSignal=strrep('dwi-signal-bet_mask.nii.gz', 'nii.gz','nhdr');
system(['ConvertBetweenFileFormats ' tmp '/dwi-signal-bet_mask.nii.gz ' tmp '/' FinalMaskSignal]);


% %%%%%%%%%%%%%%%%%%%% DO NOT EDIT THE PART ABOVE 

noise_mask=(FinalMaskNoise);
signal_mask=([tmp '/' FinalMaskSignal]);

mask_noise=loadNrrdStructure(noise_mask); %% noise mask
mask_signal=loadNrrdStructure(signal_mask); %% signal mask

%system(['rm -r ' tmp ]);

% this part loads two different masks for the signal and noise components


%axes(handles,axes2);
%for iii =1:nz
 %   imagesc(dwi.data(:,:,iii,1)) ;pause(0.1)
%end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future verbval = zeros(size(dwi.gradients, 1),1);
global tmp dwi mask_noise mask_signal signal_to_noise bval 


bval = zeros(size(dwi.gradients, 1),1);

[nx ny nz nt]=size(dwi.data);

N=length(bval);

for j =1:N
    bval(j) = dwi.bvalue * norm(dwi.gradients(j,:), 2)^2;
end

bval=round(bval);

std_signal=zeros(nt,1);
std_noise=zeros(nt, 1);
signal_to_noise=zeros(nt,1);

%mask_data=mask.data;
%se = strel('disk',6);
%mask_data = imdilate(mask_data,se);

for i =1:N
    
    signal=dwi.data(:,:,:,i);
    signal=signal(mask_signal.data==1);
    signal=double(signal);
    std_signal(i,:)=std(signal(:));


    noise=dwi.data(:,:,:,i);
    noise=noise(mask_noise.data~=1);
    noise=double(noise);
    std_noise(i,:)=std(noise(:));

    signal_to_noise(i,:)=std(signal(:))/std(noise(:));
end

axes(handles.axes1)
hold on;
grid on
plot(std_signal, '-', 'color', [1 0 1],'linewidth', 0.5);
plot(std_noise, '--', 'color', [122 16 228]/255,'linewidth', 0.5);
axes(handles.axes1)
title([' Standard Deviation of noise and signal over the number of gradients :' num2str(nt)], ...
'fontsize',15,'fontweight','bold');
legend({'Standard Deviation Signal (Brain)', 'Standard Deviation Noise (Background)'}, ...
'fontsize', 10, 'location', 'northeast');
axes(handles.axes1)
ylabel('Standard Deviation','fontsize',10);
axes(handles.axes1)
xlabel('Diffusion Weighted Image: gradient number','fontsize',10);
box off

axes(handles.axes2)
hold on 
grid on 
plot(signal_to_noise, '-','color', [122 19 223]/255, 'linewidth', 1);
axes(handles.axes2)
title(['SIGNAL TO NOISE'],.....
    'fontsize',15,'fontweight','bold');
axes(handles.axes2)
xlabel('Gradients of diffusion image','fontsize', 12, 'fontweight','bold');
axes(handles.axes2)
ylabel('Signal to Noise', 'fontsize',12,'fontweight','bold');
legend({'Signal to noise ratio'},'fontsize', 12,'location','northeast');

box off

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  signal_to_noise bval

bval=round(bval);
num_bval=length(unique(bval));
bvalues=unique(bval);
for k=1:num_bval  
    c{k}=find(bval==bvalues(k,:));
end

C=numel(c);
meanSignal=zeros(C,1);

for d=1:C
    AllSignal=signal_to_noise(c{d});
    meanSignal(d,:)=mean(AllSignal);
end
meanSignal_to_noise=[meanSignal,bvalues];

axes(handles.axes3)
hold on
grid on 
scatter(bvalues,meanSignal);
axes(handles.axes3)
title(['Mean Signal to Noise Ratio for unique B-Value: ' num2str(num_bval)],....
    'fontsize', 15, 'fontweight', 'bold')
axes(handles.axes3)
xlabel('Unique B-Value','fontsize', 12, 'fontweight','bold');
axes(handles.axes3)
ylabel('Mean Signal to Noise Ratio', 'fontsize',12,'fontweight','bold');
legend({'Mean Signal to noise ratio'},'fontsize', 12,'location','northeast');

    
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)

argout.user_data_path = handles.user_data_path;
argout.user_volume_name = handles.user_volume_name;

% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.

function figure1_WindowButtonDownFcn(hObject, eventdata, handles)

% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
