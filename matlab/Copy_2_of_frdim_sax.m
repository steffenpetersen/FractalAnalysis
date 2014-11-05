function [raw_dcm_data]=Copy_2_of_frdim_sax();


tic 

data_folder=uigetdir('*,*','select folder containing dcm data');
files=dir(strcat(data_folder,'/*.dcm'));
sl_pos=0;
sl_loc=-10000;
phs=1;

wh=waitbar(0,'DICOMS Loading - Please Wait');

% Modality,'MR'
% SeriesDescription,'SA CINE'
for fl=1:length(files)
    flnom=files(fl).name;
    if flnom=='.'
    else
       im_pn=fullfile(data_folder,flnom);
       dcm_head=dicominfo(im_pn);
       if strcmp(dcm_head.Modality,'MR')==1
           if dcm_head.SliceLocation==sl_loc
               sl_pos=sl_pos;
               phs=dcm_head.InstanceNumber;
           else
               sl_pos=sl_pos+1;
               phs=dcm_head.InstanceNumber;
               sl_loc=dcm_head.SliceLocation;
           end
           raw_dcm_data(:,:,sl_pos,phs)=dicomread(dcm_head);
       end
    end
    waitbar(fl/length(files));
end

Puid=dcm_head.PatientID;
pix_space=dcm_head.PixelSpacing;
sc=round(pix_space(1)/0.5);
%%% DEBUG %%%%
% sc=2;
%%% 
figure
colormap(gray)
imagesc(sum(sum(raw_dcm_data,3),4));
title('Create Rectangle Around LV');
rect=getrect;


for slc=1:size(raw_dcm_data,3)
    for frm=1:size(raw_dcm_data,4)
        img=raw_dcm_data(:,:,slc,frm);
        Im_cropped=imcrop(img,rect);
       Im_int=double(imresize(Im_cropped,sc,'bicubic'));
       % Im_int=double(Im_cropped);
        dcm_data(:,:,slc,frm)=Im_int;
    end
end


SL=1;
PH=1;
x=cell(size(raw_dcm_data,3),size(raw_dcm_data,4));
y=cell(size(raw_dcm_data,3),size(raw_dcm_data,4));
data={dcm_data,SL,PH,x,y,Puid};

fh=figure('Position',[100 100 500 500]);

colormap(gray)
set(fh,'userdata',data)
imagesc(dcm_data(:,:,1,1));
title(strcat('SL=',num2str(SL),'  PH=',num2str(PH)),'FontSize', 16);

fgui = uicontrol;
set(fgui,'style','frame');
set(fgui,'position',[8 3 105 18],'foregroundcolor',[0 0 0]);

fgui = uicontrol;
set(fgui,'style','frame');
set(fgui,'position',[148 3 105 18],'foregroundcolor',[0 0 0]);

sh = uicontrol(fh,'Style','slider','callback',@slider_callback,'Max',size(dcm_data,3),'Min',1,'Value',1,'SliderStep',[1/size(dcm_data,3) 1/size(dcm_data,3)],'Position',[10 12 100 30]);
sh2 = uicontrol(fh,'Style','slider','callback',@slider_callback2,'Max',size(dcm_data,4),'Min',1,'Value',1,'SliderStep',[1/size(dcm_data,4) 1/size(dcm_data,4)],'Position',[150 12 100 30]);
shtxt= uicontrol(fh,'Style','text','string','SLICE','Position',[10 5 100 14]);
shtxt2= uicontrol(fh,'Style','text','string','PHASE','Position',[150 5 100 14]);

sh3 = uicontrol(fh,'Style','pushbutton','callback',@button_callback,'string','SEGMENT','Position',[260 10 80 30]);
sh4 = uicontrol(fh,'Style','pushbutton','callback',@button2_callback,'string','FRACTAL ANALYSIS','Position',[350 10 120 30]);

function slider_callback(source,events)

SL=round(get(source,'VALUE'));
fh=get(source,'parent');
data=get(fh,'userdata');
dcm_data=data{1};
PH=data{3};
imagesc(dcm_data(:,:,SL,PH))
title(strcat('SL=',num2str(SL),'  PH=',num2str(PH)),'FontSize', 16);
x=data{4};
y=data{5};
hold on
plot(x{SL,PH},y{SL,PH},'-r')
hold off
data{2}=SL;
set(fh,'userdata',data);


function slider_callback2(source,events)

PH=round(get(source,'VALUE'));
fh=get(source,'parent');
data=get(fh,'userdata');
dcm_data=data{1};
SL=data{2};
imagesc(dcm_data(:,:,SL,PH))
title(strcat('SL=',num2str(SL),'  PH=',num2str(PH)),'FontSize', 16);
x=data{4};
y=data{5};
hold on
plot(x{SL,PH},y{SL,PH},'-r');
hold off
data{3}=PH;
set(fh,'userdata',data);

function button_callback(source,eventdata)

fh=get(source,'parent');
data=get(fh,'userdata');

dcm_data=data{1};
SL=data{2};
PH=data{3};
x=data{4};
y=data{5};

Image=dcm_data(:,:,SL,PH);
[xd,yd] = MESA_GUI_Segmentation_V1(Image);
x{SL,PH}=xd;
y{SL,PH}=yd;
data{4}=x;
data{5}=y;
set(fh,'userdata',data);

function button2_callback(source,eventdata)

fh=get(source,'parent');
data=get(fh,'userdata');

x=data{4};
y=data{5};
ID=data{6};

sz_cell=size(x);
index_cell=sz_cell(1)*sz_cell(2);

pt=1;

for fr=1:index_cell
    if isempty(x{fr})==0
        [FD]=FD2D(x{fr},y{fr});
        [SL, PH]=ind2sub(sz_cell,fr);
        tab_data(pt,:)=[SL,PH,FD]
        pt=pt+1;
    else
    end
end
%strcat - folder to collect all the mesa data strcat('dir/',ID,'txt')
fileID = fopen(strcat('/Users/catactg/code/osirix/myosirixplugins/fractal_test_images/',ID,'.txt'),'w');
fprintf(fileID,'%6s %6s %12s\n','SL','PH','FD');
fprintf(fileID,'%6.0f %6.0f %12.5f\n',tab_data');
fclose(fileID);


save (strcat('/Users/catactg/code/osirix/myosirixplugins/fractal_test_images/',ID),'data');
 toc
