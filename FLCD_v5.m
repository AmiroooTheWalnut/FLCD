function varargout = FLCD_v5(varargin)
% Written by Amir Mohammad Esmaieeli Sikaroudi
% This code is provided for researchers, although it's not suitable for
% mass running, but the user can cut the parts from this integrated code.
% 
% WARNING:
% This code is intended for academic use only, please cite the following
% paper:
% "Facility layout by collision detection and force exertion heuristics"
% https://doi.org/10.1016/j.jmsy.2016.07.001
% IN ORDER TO ACHIEVE A RELIABLE BENCHMARK FOR THE ALGORITHM, DO NOT USE 
% THE CODE MIXED WITH GUI. PLEASE USE CUT THE PARTS OF THIS CODE.
%
% FLCD_v5 MATLAB code for FLCD_v5.fig
%      FLCD_v5, by itself, creates a new FLCD_v5 or raises the existing
%      singleton*.
%
%      H = FLCD_v5 returns the handle to a new FLCD_v5 or the handle to
%      the existing singleton*.
%
%      FLCD_v5('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLCD_v5.M with the given input arguments.
%
%      FLCD_v5('Property','Value',...) creates a new FLCD_v5 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FLCD_v5_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FLCD_v5_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FLCD_v5

% Last Modified by GUIDE v2.5 28-Nov-2014 15:26:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FLCD_v5_OpeningFcn, ...
                   'gui_OutputFcn',  @FLCD_v5_OutputFcn, ...
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


% --- Executes just before FLCD_v5 is made visible.
function FLCD_v5_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FLCD_v5 (see VARARGIN)

% Choose default command line output for FLCD_v5
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


global method_selection
global nonmesh_method
global mesh_method
global fixed_size
global fixed_ratio
global flexible
global factory_limits
global fixed_size_x
global fixed_size_y
global fixed_ratio_ratio
global flexible_compress
global flexible_release
global collision_penalty
global outer_walls_penalty
global load_relation_matrix
global load_size_matrix
global relation_matrix
global size_matrix
global point_axis
global surface_axis
global exit
global cancel
global presentation
global run
global about
global save_results


method_selection=findobj(gcf,'tag','method_panel');
nonmesh_method=findobj(gcf,'tag','nonmesh_method');
mesh_method=findobj(gcf,'tag','mesh_method');
fixed_size=findobj(gcf,'tag','fixed_lw');
fixed_ratio=findobj(gcf,'tag','ratio_lw');
flexible=findobj(gcf,'tag','flexible_lw');
factory_limits=findobj(gcf,'tag','factory_limits');
fixed_size_x=findobj(gcf,'tag','fixed_size_x');
fixed_size_y=findobj(gcf,'tag','fixed_size_y');
fixed_ratio_ratio=findobj(gcf,'tag','fixed_ratio_ratio');
flexible_compress=findobj(gcf,'tag','flexible_compress');
flexible_release=findobj(gcf,'tag','flexible_release');
collision_penalty=findobj(gcf,'tag','collision_penalty');
outer_walls_penalty=findobj(gcf,'tag','outer_walls_penalty');
load_relation_matrix=findobj(gcf,'tag','load_relation_matrix');
load_size_matrix=findobj(gcf,'tag','load_size_matrix');
relation_matrix=findobj(gcf,'tag','relation_matrix');
size_matrix=findobj(gcf,'tag','size_matrix');
point_axis=findobj(gcf,'tag','point_axis');
surface_axis=findobj(gcf,'tag','surface_axis');
exit=findobj(gcf,'tag','exit');
cancel=findobj(gcf,'tag','cancel');
presentation=findobj(gcf,'tag','presentation');
run=findobj(gcf,'tag','run');
about=findobj(gcf,'tag','about');
save_results=findobj(gcf,'tag','save_results');

global cancel_val
global outer_walls_penalty_val
global collision_penalty_val
global compress_speed_val
global release_speed_val
global fixed_size_radio
global fixed_ratio_radio
global flexible_radio
global nonmesh_radio
global mesh_radio

fixed_size_radio=0;
fixed_ratio_radio=0;
flexible_radio=1;
nonmesh_radio=1;
mesh_radio=0;

cancel_val=0;
collision_penalty_val=0.9;
outer_walls_penalty_val=0.8;
compress_speed_val=7;
release_speed_val=6;

set(factory_limits,'SelectionChangeFcn',@limit_selectionfcn);
set(method_selection,'SelectionChangeFcn',@method_selectionfcn);

set(collision_penalty,'value',collision_penalty_val-0.5)
set(outer_walls_penalty,'value',outer_walls_penalty_val-0.5)
set(flexible_compress,'value',(compress_speed_val-5)/15)
set(flexible_release,'value',(release_speed_val-5)/15)


%***********************initializing
set(save_results,'enable','off')
set(cancel,'enable','off')
set(run,'enable','on')
set(relation_matrix,'data',{  blanks(0) blanks(0); blanks(0) blanks(0); blanks(0) blanks(0); blanks(0) blanks(0) });
set(size_matrix,'data',{  blanks(0) blanks(0); blanks(0) blanks(0); blanks(0) blanks(0); blanks(0) blanks(0) });
cla(point_axis)
cla(surface_axis)
radio_flexible=findobj(gcf,'tag','flexible');
radio_fixed_ratio=findobj(gcf,'tag','fixed_ratio');
radio_fixed_size=findobj(gcf,'tag','fixed_size');
set(radio_flexible,'value',flexible_radio)
set(radio_fixed_ratio,'value',fixed_ratio_radio)
set(radio_fixed_size,'value',fixed_size_radio)
set(nonmesh_method,'value',nonmesh_radio)
set(mesh_method,'value',mesh_radio)
set(fixed_size_x,'string','')
set(fixed_size_y,'string','')
set(fixed_ratio_ratio,'string','')
set(fixed_size,'visible','off')
set(fixed_ratio,'visible','off')
set(flexible,'visible','on')
%***********************initializing

%set(nonmesh_method,'visible','off')



% UIWAIT makes FLCD_v5 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FLCD_v5_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fixed_size_radio
global fixed_ratio_radio
global flexible_radio
global numpoints
global relation_matrix_val
global points
global point_axis
global surface_axis
global collision_penalty_val
global outer_walls_penalty_val
global recdatax
global recdatay
global factoryx
global factoryy
global tempxy
global compress_speed_val
global release_speed_val
global rectangeles
global z
global total_score
global nonmesh_radio
global mesh_radio
global cancel
global run
global save_results
global cancel_val
global size_table
global wall
global vertex_data
global finished_method
global sizes
global facility_projection_x
global facility_projection_y
set(run,'enable','off')
set(cancel,'enable','on')
cancel_val=0;


if mesh_radio==1

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed size^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if fixed_size_radio==1
numpoints=size(relation_matrix_val,1);
initpoints=rand(numpoints,2);
points=initpoints;
for i=1:numpoints
    for j=1:numpoints
        relation_matrix_val(j,i)=relation_matrix_val(i,j);
    end
end
step=0.001;
totalx=0;
totaly=0;
iteration_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************MOVE POINTS*******************************
while iteration_p <200
for i=1:numpoints
    for j=1:numpoints
        if i~=j
            totalx=totalx+((points(j,1)-points(i,1)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
            totaly=totaly+((points(j,2)-points(i,2)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
        end
    end
    points(i,1)=points(i,1)+totalx*step;
    points(i,2)=points(i,2)+totaly*step;
    totalx=0;
    totaly=0;
end
pause(0.0001)
scatter(point_axis,points(:,1),points(:,2),20,[iteration_p/200,1-(iteration_p/200),0])
iteration_p=iteration_p+1;
if cancel_val==1
    break;
end
end
%********************************MOVE POINTS*******************************
%********************************NORMALIZE*********************************
move=min(points);
tempxy=factoryx/factoryy;
for i=1:numpoints
            points(i,1)=points(i,1)-move(1,1);
            points(i,2)=points(i,2)-move(1,2);
end
if tempxy>=1
    normalize_val=max(points)*tempxy;
end
if tempxy<1
    normalize_val=max(points)/tempxy;
end
for i=1:numpoints
        points(i,1)=points(i,1)/normalize_val(1,1);
        points(i,2)=points(i,2)/normalize_val(1,2);
end
%********************************NORMALIZE*********************************
%%manual snap to grid
gridsize=0.01;
%%manual snap to grid

scatter(point_axis,points(:,1),points(:,2),24,[0,0,1])
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',point_axis)
end
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iteration_m=0;
max_iteration_m=50;
%iteration_s=0;
%max_iteration_s=20;
[sx,sy]=meshgrid(-factoryx:gridsize:2*factoryx,-factoryy:gridsize:2*factoryy);
gridacu(1)=size(sx,1)/3;
gridacu(2)=size(sy,2)/3;
%**************************VERTEX ADJUSTMENT*******************************
vertex_data=cell(numpoints,2);
z(4,numpoints)=0;
if size(size_table,1)==2
    sizes(numpoints)=0;
    for i=1:numpoints
        vertex_data{i,1}(1,1)=points(i,1)-(size_table(1,i)/2);
        vertex_data{i,2}(1,1)=points(i,2)-(size_table(2,i)/2);
        vertex_data{i,1}(2,1)=points(i,1)+(size_table(1,i)/2);
        vertex_data{i,2}(2,1)=points(i,2)-(size_table(2,i)/2);
        vertex_data{i,1}(3,1)=points(i,1)+(size_table(1,i)/2);
        vertex_data{i,2}(3,1)=points(i,2)+(size_table(2,i)/2);
        vertex_data{i,1}(4,1)=points(i,1)-(size_table(1,i)/2);
        vertex_data{i,2}(4,1)=points(i,2)+(size_table(2,i)/2);
        z(i)=i;
        sizes(i)=4;
    end
elseif size(size_table,1)>2
    sizes(numpoints)=0;
    for j=1:2:size(size_table,2)
        for i=1:size(size_table,1)
            if isnan(size_table(i,j))==1
                sizes(round(j/2))=i-1;
                break
            elseif i==size(size_table,1)
                sizes(round(j/2))=i;
                break
            end
        
        end
    end
    for i=1:size(sizes,2)
        for j=1:sizes(1,i)
            vertex_data{i,2}(j,1)=size_table(j,2*i);
            vertex_data{i,1}(j,1)=size_table(j,2*i-1);
        end
        z(i)=i;
    end
end
%**************************VERTEX ADJUSTMENT*******************************
%freespace=0;
%while iteration_s<max_iteration_s
while iteration_m<max_iteration_m
%************************PLOTS*********************************************
if cancel_val==1
    break;
end
pause(0.0001)
grid on;
%{
centeredx(4,numpoints)=0;
centeredy(4,numpoints)=0;
for i=1:numpoints
centeredx(1,i)=points(i,1)-(recdatax(i)/2);
centeredy(1,i)=points(i,2)-(recdatay(i)/2);
centeredx(2,i)=points(i,1)+(recdatax(i)/2);
centeredy(2,i)=points(i,2)-(recdatay(i)/2);
centeredx(3,i)=points(i,1)+(recdatax(i)/2);
centeredy(3,i)=points(i,2)+(recdatay(i)/2);
centeredx(4,i)=points(i,1)-(recdatax(i)/2);
centeredy(4,i)=points(i,2)+(recdatay(i)/2);
for j=1:4
z(j,i)=i;
end
end
%}
for i=1:numpoints
    center_x_move=points(i,1)-mean(vertex_data{i,1});
    vertex_data{i,1}=vertex_data{i,1}+center_x_move;
    
    center_y_move=points(i,2)-mean(vertex_data{i,2});
    vertex_data{i,2}=vertex_data{i,2}+center_y_move;
end
axis([-0.1 factoryx+0.1 -0.1 factoryy+0.1])
cla(surface_axis)
for i=1:numpoints
    patch(vertex_data{i,1},vertex_data{i,2},z(i),'parent',surface_axis)
end
%patch(centeredx,centeredy,z,'parent',surface_axis)
rectangle('Position',[0,0,factoryx,factoryy],'EdgeColor',[1,0,0],'parent',surface_axis)
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',surface_axis)
end
%************************PLOTS*********************************************
%************************SURFACE*******************************************
IN(size(sx,1),size(sy,2),numpoints)=0;
for i=1:numpoints
    IN(:,:,i)=double(inpolygon(sx,sy,vertex_data{i,1},vertex_data{i,2}));
end
IN_total(size(sx,1),size(sy,2))=0;
for i=1:numpoints
    IN_total=IN_total+double(inpolygon(sx,sy,vertex_data{i,1},vertex_data{i,2}));
end
%************************SURFACE*******************************************
%************************FACILITY SURFACE**********************************
if cancel_val==1
    break;
end
IN_weighted(size(IN_total,1),size(IN_total,2),numpoints)=0;
for k=1:numpoints
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if IN(i,j,k)>0
                if IN_total(i,j)>1
                    IN_weighted(i,j,k)=1-(collision_penalty_val*IN_total(i,j));
                end
                if IN_total(i,j)==1
                    IN_weighted(i,j,k)=1;
                end
                if i<gridacu(1)+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(gridacu(1)+1-i)*0.5;
                end
                if j<gridacu(2)+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(gridacu(2)+1-j)*0.5;
                end
                if i>(2*gridacu(1))+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(i-(2*gridacu(1))+1)*0.5;
                end
                if j>(2*gridacu(2))+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(j-(2*gridacu(2))+1)*0.5;
                end
            end
            if IN_weighted(i,j,k)==0
                IN_weighted(i,j,k)=NaN;
            end
        end
    end
end
%************************FACILITY SURFACE**********************************
%************************COLLISION AMOUNT**********************************
%count=0;
%sum_w=0;
%collision_amount(numpoints)=0;
%S(numpoints)=0;
%for k=1:numpoints
%    for i=1:size(sx)
%        for j=1:size(sy)
%            if IN_weighted(i,j,k)>0
%                count=count+1;
%                sum_w=sum_w+IN_weighted(i,j,k);
%            end
%        end
%    end
%    collision_amount(k)=sum_w/count;
%    S(k)=count;%surface of each workstation;
%    count=0;
%    sum_w=0;
%end
%************************COLLISION AMOUNT**********************************
%**************************NEW POINT***************************************
for I=1:numpoints
newx=0;
m=1;
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if isnan(IN_weighted(i,j,I))==0
                temp(m)=IN_weighted(i,j,I);
                tempsx(m)=sx(i,j);
                m=m+1;
            end
        end
    end
    tempsx=tempsx-points(I,1);
    newx=mean(tempsx.*temp);
    points(I,1)=points(I,1)+newx;

temp=[];
tempsx=[];
newy=0;
m=1;
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if isnan(IN_weighted(i,j,I))==0
                temp(m)=IN_weighted(i,j,I);
                tempsy(m)=sy(i,j);
                m=m+1;
            end
        end
    end
    tempsy=tempsy-points(I,2);
    newy=mean(tempsy.*temp);
    points(I,2)=points(I,2)+newy;
    temp=[];
    tempsy=[];
end
%**************************NEW POINT***************************************

%**************************FREE SPACE AMOUNT*******************************
%freespace=0;
%    for i=round(gridacu(1)):round((2*gridacu(1)))
%        for j=round(gridacu(2)):round((2*gridacu(2)))
%            if IN_total(i,j)==0
%                freespace=freespace+1;
%            end
%        end
%    end
%freespace=freespace/(size(sx,1)*size(sy,2));
%**************************FREE SPACE AMOUNT*******************************

if cancel_val==1
    break;
end

iteration_m=iteration_m+1;
if iteration_m<max_iteration_m
IN=[];
IN_total=[];
IN_weighted=[];
end
end
%**************************RESIZE FACTORY**********************************
%iteration_s=iteration_s+1;
%if freespace>0.038
%    factoryx=factoryx*0.95;
%    factoryy=factoryy*0.95;
%elseif freespace<0.005
%    factoryx=factoryx/0.95;
%    factoryy=factoryy/0.95;
%else
%    break
%end
%**************************RESIZE FACTORY**********************************
%iteration_m=0;
%[sx,sy]=meshgrid(-factoryx:gridsize:2*factoryx,-factoryy:gridsize:2*factoryy);
%gridacu(1)=size(sx,1)/3;
%gridacu(2)=size(sy,2)/3;
%IN=[];
%IN_total=[];
%IN_weighted=[];
%end
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed size^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%#########################################################################
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed ratio^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if fixed_ratio_radio==1
numpoints=size(relation_matrix_val,1);
initpoints=rand(numpoints,2);
points=initpoints;
for i=1:numpoints
    for j=1:numpoints
        relation_matrix_val(j,i)=relation_matrix_val(i,j);
    end
end
step=0.001;
totalx=0;
totaly=0;
iteration_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************MOVE POINTS*******************************
while iteration_p <200
for i=1:numpoints
    for j=1:numpoints
        if i~=j
            totalx=totalx+((points(j,1)-points(i,1)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
            totaly=totaly+((points(j,2)-points(i,2)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
        end
    end
    points(i,1)=points(i,1)+totalx*step;
    points(i,2)=points(i,2)+totaly*step;
    totalx=0;
    totaly=0;
end
pause(0.0001)
scatter(point_axis,points(:,1),points(:,2),20,[iteration_p/200,1-(iteration_p/200),0])
iteration_p=iteration_p+1;
if cancel_val==1
    break;
end
end
%********************************MOVE POINTS*******************************
%********************************NORMALIZE*********************************
move=min(points);
for i=1:numpoints
            points(i,1)=points(i,1)-move(1,1);
            points(i,2)=points(i,2)-move(1,2);
end
if tempxy>=1
    normalize_val=max(points)*tempxy;
end
if tempxy<1
    normalize_val=max(points)/tempxy;
end
for i=1:numpoints
        points(i,1)=points(i,1)/normalize_val(1,1);
        points(i,2)=points(i,2)/normalize_val(1,2);
end
%********************************NORMALIZE*********************************
%%manual snap to grid
gridsize=0.01;
%points=points*(100/4);
%points=round(points);
%points=points/(100/4);
%%manual snap to grid

scatter(point_axis,points(:,1),points(:,2),24,[0,0,1])
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',point_axis)
end
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************VERTEX ADJUSTMENT*******************************
vertex_data=cell(numpoints,2);
z(4,numpoints)=0;
if size(size_table,1)==2
    sizes(numpoints)=0;
    for i=1:numpoints
        vertex_data{i,1}(1,1)=points(i,1)-(size_table(1,i)/2);
        vertex_data{i,2}(1,1)=points(i,2)-(size_table(2,i)/2);
        vertex_data{i,1}(2,1)=points(i,1)+(size_table(1,i)/2);
        vertex_data{i,2}(2,1)=points(i,2)-(size_table(2,i)/2);
        vertex_data{i,1}(3,1)=points(i,1)+(size_table(1,i)/2);
        vertex_data{i,2}(3,1)=points(i,2)+(size_table(2,i)/2);
        vertex_data{i,1}(4,1)=points(i,1)-(size_table(1,i)/2);
        vertex_data{i,2}(4,1)=points(i,2)+(size_table(2,i)/2);
        z(i)=i;
        sizes(i)=4;
    end
elseif size(size_table,1)>2
    sizes(numpoints)=0;
    for j=1:2:size(size_table,2)
        for i=1:size(size_table,1)
            if isnan(size_table(i,j))==1
                sizes(round(j/2))=i-1;
                break
            elseif i==size(size_table,1)
                sizes(round(j/2))=i;
                break
            end
        
        end
    end
    for i=1:size(sizes,2)
        for j=1:sizes(1,i)
            vertex_data{i,2}(j,1)=size_table(j,2*i);
            vertex_data{i,1}(j,1)=size_table(j,2*i-1);
        end
        z(i)=i;
    end
end
%**************************VERTEX ADJUSTMENT*******************************
facilities_surfaces(numpoints)=0;
for i=1:numpoints
    facilities_surfaces(i)=polyarea(vertex_data{i,1},vertex_data{i,2});
end

if sum(facilities_surfaces)>1 && tempxy>=1
    for i=1:numpoints
        vertex_data{i,1}=((vertex_data{i,1})./sqrt(sum(facilities_surfaces)))./tempxy;
        vertex_data{i,2}=((vertex_data{i,2})./sqrt(sum(facilities_surfaces)))./tempxy;
    end
end
if sum(facilities_surfaces)>1 && tempxy<1
    for i=1:numpoints
        vertex_data{i,1}=((vertex_data{i,1})./sqrt(sum(facilities_surfaces))).*tempxy;
        vertex_data{i,2}=((vertex_data{i,2})./sqrt(sum(facilities_surfaces))).*tempxy;
    end
end

for i=1:numpoints
    facilities_surfaces(i)=polyarea(vertex_data{i,1},vertex_data{i,2});
end
if tempxy>=1
    factoryx=1;
    factoryy=1/tempxy;
end
if tempxy<1
    factoryx=tempxy;
    factoryy=1;
end
iteration_m=0;
max_iteration_m=40;
iteration_s=0;
max_iteration_s=60;
[sx,sy]=meshgrid(-factoryx:gridsize:2*factoryx,-factoryy:gridsize:2*factoryy);
gridacu(1)=size(sx,1)/3;
gridacu(2)=size(sy,2)/3;
freespace=0;
move_x=0;%amount of x forces of compressed elements, to decide whether to move factory walls or not
move_y=0;%amount of y forces of compressed elements, to decide whether to move factory walls or not
while iteration_s<max_iteration_s
if cancel_val==1
    break;
end
while iteration_m<max_iteration_m
if cancel_val==1
    break;
end
%************************PLOTS*********************************************
pause(0.0001)
grid on;
%{
centeredx(4,numpoints)=0;
centeredy(4,numpoints)=0;
z(4,numpoints)=0;
for i=1:numpoints
centeredx(1,i)=points(i,1)-(recdatax(i)/2);
centeredy(1,i)=points(i,2)-(recdatay(i)/2);
centeredx(2,i)=points(i,1)+(recdatax(i)/2);
centeredy(2,i)=points(i,2)-(recdatay(i)/2);
centeredx(3,i)=points(i,1)+(recdatax(i)/2);
centeredy(3,i)=points(i,2)+(recdatay(i)/2);
centeredx(4,i)=points(i,1)-(recdatax(i)/2);
centeredy(4,i)=points(i,2)+(recdatay(i)/2);
for j=1:4
z(j,i)=i;
end
end
%}
for i=1:numpoints
    center_x_move=points(i,1)-mean(vertex_data{i,1});
    vertex_data{i,1}=vertex_data{i,1}+center_x_move;
    
    center_y_move=points(i,2)-mean(vertex_data{i,2});
    vertex_data{i,2}=vertex_data{i,2}+center_y_move;
end
axis([-0.1 factoryx+0.1 -0.1 factoryy+0.1])
cla(surface_axis)
for i=1:numpoints
    patch(vertex_data{i,1},vertex_data{i,2},z(i),'parent',surface_axis)
end
rectangle('Position',[0,0,factoryx,factoryy],'EdgeColor',[1,0,0],'parent',surface_axis)
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',surface_axis)
end
%************************PLOTS*********************************************
%************************SURFACE*******************************************
IN(size(sx,1),size(sy,2),numpoints)=0;
for i=1:numpoints
    IN(:,:,i)=double(inpolygon(sx,sy,vertex_data{i,1},vertex_data{i,2}));
end
IN_total(size(sx,1),size(sy,2))=0;
for i=1:numpoints
    IN_total=IN_total+double(inpolygon(sx,sy,vertex_data{i,1},vertex_data{i,2}));
end
if cancel_val==1
    break;
end
%************************SURFACE*******************************************
%************************FACILITY SURFACE**********************************
if cancel_val==1
    break;
end
IN_weighted(size(IN_total,1),size(IN_total,2),numpoints)=0;
for k=1:numpoints
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if IN(i,j,k)>0
                if IN_total(i,j)>1
                    IN_weighted(i,j,k)=1-(collision_penalty_val*IN_total(i,j));
                end
                if IN_total(i,j)==1
                    IN_weighted(i,j,k)=1;
                end
                if i<gridacu(1)+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(gridacu(1)+1-i)*0.5;
                    move_y=move_y+(gridacu(1)+1-i)*0.5;
                end
                if j<gridacu(2)+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(gridacu(2)+1-j)*0.5;
                    move_x=move_x+(gridacu(2)+1-j)*0.5;
                end
                if i>(2*gridacu(1))+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(i-(2*gridacu(1))+1)*0.5;
                    move_y=move_y+(i-(2*gridacu(1))+1)*0.5;
                end
                if j>(2*gridacu(2))+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(j-(2*gridacu(2))+1)*0.5;
                    move_x=move_x+(j-(2*gridacu(2))+1)*0.5;
                end
            end
            if IN_weighted(i,j,k)==0
                IN_weighted(i,j,k)=NaN;
            end
        end
    end
end
%************************FACILITY SURFACE**********************************
%************************COLLISION AMOUNT**********************************
%count=0;
%sum_w=0;
%collision_amount(numpoints)=0;
%S(numpoints)=0;
%for k=1:numpoints
%    for i=1:size(sx)
%        for j=1:size(sy)
%            if IN_weighted(i,j,k)>0
%                count=count+1;
%                sum_w=sum_w+IN_weighted(i,j,k);
%            end
%        end
%    end
%    collision_amount(k)=sum_w/count;
%    S(k)=count;%surface of each workstation;
%    count=0;
%    sum_w=0;
%end
%************************COLLISION AMOUNT**********************************
%**************************NEW POINT***************************************
for I=1:numpoints
newx=0;
m=1;
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if isnan(IN_weighted(i,j,I))==0
                temp(m)=IN_weighted(i,j,I);
                tempsx(m)=sx(i,j);
                m=m+1;
            end
        end
    end
    tempsx=tempsx-points(I,1);
    newx=mean(tempsx.*temp);
    points(I,1)=points(I,1)+newx;

temp=[];
tempsx=[];
newy=0;
m=1;
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if isnan(IN_weighted(i,j,I))==0
                temp(m)=IN_weighted(i,j,I);
                tempsy(m)=sy(i,j);
                m=m+1;
            end
        end
    end
    tempsy=tempsy-points(I,2);
    newy=mean(tempsy.*temp);
    points(I,2)=points(I,2)+newy;
    temp=[];
    tempsy=[];
end
%**************************NEW POINT***************************************

%**************************FREE SPACE AMOUNT*******************************
freespace=0;
    for i=round(gridacu(1)):round((2*gridacu(1)))
        for j=round(gridacu(2)):round((2*gridacu(2)))
            if IN_total(i,j)==0
                freespace=freespace+1;
            end
        end
    end
freespace=9*freespace/(size(sx,1)*size(sy,2));
%**************************FREE SPACE AMOUNT*******************************
if cancel_val==1
    break;
end
iteration_m=iteration_m+1;
if iteration_m<max_iteration_m
IN=[];
IN_total=[];
IN_weighted=[];
end
end
%**************************RESIZE FACTORY**********************************
iteration_s=iteration_s+1;
if move_x+move_y==0 || freespace>0.04
    factoryx=factoryx*0.95;
    factoryy=factoryy*0.95;
elseif move_x+move_y>200
    factoryx=factoryx/0.96;
    factoryy=factoryy/0.96;
else
    break
end
move_x=0;
move_y=0;
%**************************RESIZE FACTORY**********************************
iteration_m=0;
[sx,sy]=meshgrid(-factoryx:gridsize:2*factoryx,-factoryy:gridsize:2*factoryy);
gridacu(1)=size(sx,1)/3;
gridacu(2)=size(sy,2)/3;
IN=[];
IN_total=[];
IN_weighted=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed ratio^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%#########################################################################
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Flexible^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if flexible_radio==1
numpoints=size(relation_matrix_val,1);
initpoints=rand(numpoints,2);
points=initpoints;
for i=1:numpoints
    for j=1:numpoints
        relation_matrix_val(j,i)=relation_matrix_val(i,j);
    end
end
step=0.001;
totalx=0;
totaly=0;
iteration_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************MOVE POINTS*******************************
while iteration_p <200
for i=1:numpoints
    for j=1:numpoints
        if i~=j
            totalx=totalx+((points(j,1)-points(i,1)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
            totaly=totaly+((points(j,2)-points(i,2)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
        end
    end
    points(i,1)=points(i,1)+totalx*step;
    points(i,2)=points(i,2)+totaly*step;
    totalx=0;
    totaly=0;
end
refresh
pause(0.0001)
scatter(point_axis,points(:,1),points(:,2),20,[iteration_p/200,1-(iteration_p/200),0])
hold on;
iteration_p=iteration_p+1;
if cancel_val==1
    break;
end
end
%********************************MOVE POINTS*******************************
%********************************NORMALIZE*********************************
move=min(points);
for i=1:numpoints

            points(i,1)=points(i,1)-move(1,1);
            points(i,2)=points(i,2)-move(1,2);

end
normalize_val=max(points);
for i=1:numpoints
        points(i,1)=points(i,1)/normalize_val(1,1);
        points(i,2)=points(i,2)/normalize_val(1,2);
end
%********************************NORMALIZE*********************************
%%manual snap to grid
gridsize=0.01;
%%manual snap to grid

scatter(point_axis,points(:,1),points(:,2),24,[0,0,1])
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',point_axis)
end
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************VERTEX ADJUSTMENT*******************************
vertex_data=cell(numpoints,2);
z(4,numpoints)=0;
if size(size_table,1)==2
    sizes(numpoints)=0;
    for i=1:numpoints
        vertex_data{i,1}(1,1)=points(i,1)-(size_table(1,i)/2);
        vertex_data{i,2}(1,1)=points(i,2)-(size_table(2,i)/2);
        vertex_data{i,1}(2,1)=points(i,1)+(size_table(1,i)/2);
        vertex_data{i,2}(2,1)=points(i,2)-(size_table(2,i)/2);
        vertex_data{i,1}(3,1)=points(i,1)+(size_table(1,i)/2);
        vertex_data{i,2}(3,1)=points(i,2)+(size_table(2,i)/2);
        vertex_data{i,1}(4,1)=points(i,1)-(size_table(1,i)/2);
        vertex_data{i,2}(4,1)=points(i,2)+(size_table(2,i)/2);
        z(i)=i;
        sizes(i)=4;
    end
elseif size(size_table,1)>2
    sizes(numpoints)=0;
    for j=1:2:size(size_table,2)
        for i=1:size(size_table,1)
            if isnan(size_table(i,j))==1
                sizes(round(j/2))=i-1;
                break
            elseif i==size(size_table,1)
                sizes(round(j/2))=i;
                break
            end
        
        end
    end
    for i=1:size(sizes,2)
        for j=1:sizes(1,i)
            vertex_data{i,2}(j,1)=size_table(j,2*i);
            vertex_data{i,1}(j,1)=size_table(j,2*i-1);
        end
        z(i)=i;
    end
end
%**************************VERTEX ADJUSTMENT*******************************
tempx=0;
tempy=0;
for i=1:numpoints
    tempx=tempx+mean(vertex_data{i,1});
    tempy=tempx+mean(vertex_data{i,2});
end
tempx=tempx/numpoints;
tempy=tempy/numpoints;
tempxy=tempy/tempx;
facilities_surfaces(numpoints)=0;
for i=1:numpoints
    facilities_surfaces(i)=polyarea(vertex_data{i,1},vertex_data{i,2});
end

if sum(facilities_surfaces)>1 && tempxy>=1
    for i=1:numpoints
        vertex_data{i,1}=((vertex_data{i,1})./sqrt(sum(facilities_surfaces)))./tempxy;
        vertex_data{i,2}=((vertex_data{i,2})./sqrt(sum(facilities_surfaces)))./tempxy;
    end
end
if sum(facilities_surfaces)>1 && tempxy<1
    for i=1:numpoints
        vertex_data{i,1}=((vertex_data{i,1})./sqrt(sum(facilities_surfaces))).*tempxy;
        vertex_data{i,2}=((vertex_data{i,2})./sqrt(sum(facilities_surfaces))).*tempxy;
    end
end

for i=1:numpoints
    facilities_surfaces(i)=polyarea(vertex_data{i,1},vertex_data{i,2});
end
factoryx=1;
factoryy=tempxy;
iteration_m=0;
max_iteration_m=compress_speed_val;
iteration_s=0;
max_iteration_s=500;
[sx,sy]=meshgrid(-factoryx:gridsize:2*factoryx,-factoryy:gridsize:2*factoryy);
gridacu(1)=size(sx,1)/3;
gridacu(2)=size(sy,2)/3;
freespace=0;
free_comarioson=0;
compress_release=0;
move_x=0;%amount of x forces of compressed elements, to decide whether move factory walls in x direction or not
move_y=0;%amount of y forces of compressed elements, to decide whether move factory walls in y direction or not
while iteration_s<max_iteration_s
if cancel_val==1
    break;
end
while iteration_m<max_iteration_m
%************************PLOTS*********************************************
pause(0.0001)
grid on;
%{
centeredx(4,numpoints)=0;
centeredy(4,numpoints)=0;
z(4,numpoints)=0;
for i=1:numpoints
centeredx(1,i)=points(i,1)-(recdatax(i)/2);
centeredy(1,i)=points(i,2)-(recdatay(i)/2);
centeredx(2,i)=points(i,1)+(recdatax(i)/2);
centeredy(2,i)=points(i,2)-(recdatay(i)/2);
centeredx(3,i)=points(i,1)+(recdatax(i)/2);
centeredy(3,i)=points(i,2)+(recdatay(i)/2);
centeredx(4,i)=points(i,1)-(recdatax(i)/2);
centeredy(4,i)=points(i,2)+(recdatay(i)/2);
for j=1:4
z(j,i)=i;
end
end
%}
for i=1:numpoints
    center_x_move=points(i,1)-mean(vertex_data{i,1});
    vertex_data{i,1}=vertex_data{i,1}+center_x_move;
    
    center_y_move=points(i,2)-mean(vertex_data{i,2});
    vertex_data{i,2}=vertex_data{i,2}+center_y_move;
end
axis([-0.1 factoryx+0.1 -0.1 factoryy+0.1])
cla(surface_axis)
for i=1:numpoints
    patch(vertex_data{i,1},vertex_data{i,2},z(i),'parent',surface_axis)
end
rectangle('Position',[0,0,factoryx,factoryy],'EdgeColor',[1,0,0],'parent',surface_axis)
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',surface_axis)
end
%************************PLOTS*********************************************
%************************SURFACE*******************************************
IN(size(sx,1),size(sy,2),numpoints)=0;
for i=1:numpoints
    IN(:,:,i)=double(inpolygon(sx,sy,vertex_data{i,1},vertex_data{i,2}));
end
IN_total(size(sx,1),size(sy,2))=0;
for i=1:numpoints
    IN_total=IN_total+double(inpolygon(sx,sy,vertex_data{i,1},vertex_data{i,2}));
end
if cancel_val==1
    break;
end
%************************SURFACE*******************************************
%************************FACILITY SURFACE**********************************
IN_weighted(size(IN_total,1),size(IN_total,2),numpoints)=0;
for k=1:numpoints
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if IN(i,j,k)>0
                if IN_total(i,j)>1
                    IN_weighted(i,j,k)=1-(collision_penalty_val*IN_total(i,j));
                end
                if IN_total(i,j)==1
                    IN_weighted(i,j,k)=1;
                end
                if i<gridacu(1)+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(gridacu(1)+1-i)*0.5;
                    move_y=move_y+(gridacu(1)+1-i)*0.5;
                end
                if j<gridacu(2)+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(gridacu(2)+1-j)*0.5;
                    move_x=move_x+(gridacu(2)+1-j)*0.5;
                end
                if i>(2*gridacu(1))+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(i-(2*gridacu(1))+1)*0.5;
                    move_y=move_y+(i-(2*gridacu(1))+1)*0.5;
                end
                if j>(2*gridacu(2))+1
                    IN_weighted(i,j,k)=-outer_walls_penalty_val-(j-(2*gridacu(2))+1)*0.5;
                    move_x=move_x+(j-(2*gridacu(2))+1)*0.5;
                end
            end
            if IN_weighted(i,j,k)==0
                IN_weighted(i,j,k)=NaN;
            end
        end
    end
end
%************************FACILITY SURFACE**********************************
%************************COLLISION AMOUNT**********************************
%count=0;
%sum_w=0;
%collision_amount(numpoints)=0;
%S(numpoints)=0;
%for k=1:numpoints
%    for i=1:size(sx)
%        for j=1:size(sy)
%            if IN_weighted(i,j,k)>0
%                count=count+1;
%                sum_w=sum_w+IN_weighted(i,j,k);
%            end
%        end
%    end
%    collision_amount(k)=sum_w/count;
%    S(k)=count;%surface of each workstation;
%    count=0;
%    sum_w=0;
%end
%************************COLLISION AMOUNT**********************************
%**************************NEW POINT***************************************
for I=1:numpoints
newx=0;
m=1;
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if isnan(IN_weighted(i,j,I))==0
                temp(m)=IN_weighted(i,j,I);
                tempsx(m)=sx(i,j);
                m=m+1;
            end
        end
    end
    tempsx=tempsx-points(I,1);
    newx=mean(tempsx.*temp);
    points(I,1)=points(I,1)+newx;

temp=[];
tempsx=[];
newy=0;
m=1;
    for i=1:size(sx,1)
        for j=1:size(sy,2)
            if isnan(IN_weighted(i,j,I))==0
                temp(m)=IN_weighted(i,j,I);
                tempsy(m)=sy(i,j);
                m=m+1;
            end
        end
    end
    tempsy=tempsy-points(I,2);
    newy=mean(tempsy.*temp);
    points(I,2)=points(I,2)+newy;
    temp=[];
    tempsy=[];
end
%**************************NEW POINT***************************************

%**************************FREE SPACE AMOUNT*******************************
if cancel_val==1
    break;
end
freespace=0;
    for i=round(gridacu(1)):round((2*gridacu(1)))
        for j=round(gridacu(2)):round((2*gridacu(2)))
            if IN_total(i,j)==0
                freespace=freespace+1;
            end
        end
    end
freespace=9*freespace/(size(sx,1)*size(sy,2));
if abs(freespace-free_comarioson)<0.0001
    break
else
    free_comarioson=freespace;
end
if freespace<0.005
    compress_release=1;%release factory
    max_iteration_m=release_speed_val;
end
%if freespace>0.005
    %compress_release=0;%compress factory
%end
%**************************FREE SPACE AMOUNT*******************************

iteration_m=iteration_m+1;
if iteration_m<max_iteration_m
IN=[];
IN_total=[];
IN_weighted=[];
end
end
%**************************RESIZE FACTORY**********************************
iteration_s=iteration_s+1;
if compress_release==1 && move_x<=6 && move_y<=6
    break
end
if compress_release==0 && (move_x>0 || move_y>0)
    temp=max(move_x,move_y);
    move_x=move_x/temp;
    move_y=move_y/temp;
        factoryx=factoryx*(1-(move_y)*0.02);
        factoryy=factoryy*(1-(move_x)*0.02);
end
if compress_release==1 && (move_x>0 || move_y>0)
    temp=max(move_x,move_y);
    move_x=move_x/temp;
    move_y=move_y/temp;
        factoryx=factoryx/(1-(move_x)*0.01);
        factoryy=factoryy/(1-(move_y)*0.01);
end
if compress_release==0 && move_x==0 && move_y==0
        factoryx=factoryx*(1-0.02);
        factoryy=factoryy*(1-0.02);
end
if compress_release==1 && move_x==0 && move_y==0
        factoryx=factoryx/(1-0.01);
        factoryy=factoryy/(1-0.01);
end
move_x=0;
move_y=0;
%**************************RESIZE FACTORY**********************************
iteration_m=0;
[sx,sy]=meshgrid(-factoryx:gridsize:2*factoryx,-factoryy:gridsize:2*factoryy);
gridacu(1)=size(sx,1)/3;
gridacu(2)=size(sy,2)/3;
IN=[];
IN_total=[];
IN_weighted=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Flexible^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%Score layout%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***********************Departments distance score*************************
score_dist=0;
for i=1:numpoints
    for j=i+1:numpoints
        score_dist=score_dist+relation_matrix_val(i,j)*(1/(sqrt(((points(i,1)-points(j,1))^2)+((points(i,2)-points(j,2))^2))));
    end
end
%***********************Departments distance score*************************
%************************Surface usage score*******************************
score_surf=1/(factoryx*factoryy);
%************************Surface usage score*******************************
total_score=score_dist+score_surf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%Score layout%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

if nonmesh_radio==1
if size(size_table,1)>2
    errordlg('Facilities are not rectangulare or not in row format!You can select mesh format to support non-rectangles','Bad Input','modal')
    set(run,'enable','on')
    set(cancel,'enable','off')
    return
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed size^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if fixed_size_radio==1
numpoints=size(relation_matrix_val,1);
initpoints=rand(numpoints,2);
points=initpoints;
for i=1:numpoints
    for j=1:numpoints
        relation_matrix_val(j,i)=relation_matrix_val(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%new
recdatax=size_table(1,:);
recdatay=size_table(2,:);
rectangeles_init_data=([recdatax;recdatay]);%first row is lengths(x) of rectangels and second is widths(y)
rectangles_surfaces(numpoints)=0;
for i=1:numpoints
    rectangles_surfaces(i)=rectangeles_init_data(1,i)*rectangeles_init_data(2,i);
end

%COMPASS
%    1        4
%     --------
%    |        |
%    |        |
%    |        |
%     --------
%    2        3
%%%%%%%%%%%%%%%%%%%%%new

score_dist=0;
%score_surf=0;
%total_score=0;
%for test=1:test_number
%%%%%%%%%%%%%%%%%%%%%new
rectangeles(4,numpoints,2)=0;%index 1 is compass index 2 is point index 3 is x or y
z(4,numpoints)=0;
rectangles_move(numpoints,2)=0;
wall_init=sum(rectangeles_init_data,2)/0.1;
wall=[0,0;factoryx,factoryy];%START_X-START_Y-WIDTH-HEIGHT
%wall_force(2,2)=0;%FIRST ROW IS IN HORIZENTAL DIRECTION(X) FIRST COLUMN IN
%FIRST ROW IS LEFT WALL AND SECOND COLUMN IS RIGHT WALL
%SECOND ROW IS IN VERTICAL DIRECTION(Y) FIRST COLUMN IS LOWER WALL AND
%SECOND COLUMN IS UPPER WALL
%coverted_wall_force(2,2)=0;
facility_projection_x(numpoints,2)=0;
facility_projection_y(numpoints,2)=0;
target_facility_projection_x(1,2)=0;
target_facility_projection_y(1,2)=0;
status=0;%ZERO MEANS LARGING-ONE MEANS SHRINKING
max_num_adjustment=1;
num_adjustment=0;
%%%%%%%%%%%%%%%%%%%%%new

step=0.001;
totalx=0;
totaly=0;
iteration_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************MOVE POINTS*******************************
while iteration_p <200
for i=1:numpoints
    for j=1:numpoints
        if i~=j
            totalx=totalx+((points(j,1)-points(i,1)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
            totaly=totaly+((points(j,2)-points(i,2)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
        end
    end
    points(i,1)=points(i,1)+totalx*step;
    points(i,2)=points(i,2)+totaly*step;
    totalx=0;
    totaly=0;
end
refresh
pause(0.0001)
scatter(point_axis,points(:,1),points(:,2),20,[iteration_p/200,1-(iteration_p/200),0])
hold on;
iteration_p=iteration_p+1;
if cancel_val==1
    break;
end
end
%********************************MOVE POINTS*******************************
%********************************NORMALIZE*********************************
move=min(points);
for i=1:numpoints

            points(i,1)=points(i,1)-move(1,1);
            points(i,2)=points(i,2)-move(1,2);

end
normalize_val=max(points);
for i=1:numpoints
        points(i,1)=points(i,1)/normalize_val(1,1);
        points(i,2)=points(i,2)/normalize_val(1,2);
end
%********************************NORMALIZE*********************************

scatter(point_axis,points(:,1),points(:,2),24,[0,0,1])
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',point_axis)
end
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:10000
if cancel_val==1
    break;
end
for i=1:numpoints
    rectangeles(1,i,1)=points(i,1)-rectangeles_init_data(1,i)/2;
    rectangeles(2,i,1)=points(i,1)-rectangeles_init_data(1,i)/2;
    rectangeles(3,i,1)=points(i,1)+rectangeles_init_data(1,i)/2;
    rectangeles(4,i,1)=points(i,1)+rectangeles_init_data(1,i)/2;
    
    rectangeles(1,i,2)=points(i,2)+rectangeles_init_data(2,i)/2;
    rectangeles(2,i,2)=points(i,2)-rectangeles_init_data(2,i)/2;
    rectangeles(3,i,2)=points(i,2)-rectangeles_init_data(2,i)/2;
    rectangeles(4,i,2)=points(i,2)+rectangeles_init_data(2,i)/2;
    for j=1:4
        z(j,i)=i;
    end
end
cla
patch(rectangeles(:,:,1),rectangeles(:,:,2),z,'parent',surface_axis)
rectangle('Position',[wall(1,1),wall(1,2),wall(2,1),wall(2,2)],'EdgeColor',[1,0,0],'parent',surface_axis)
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',surface_axis)
end
pause(0.000001);
%*****************************PROJECTIONS**********************************
for i=1:numpoints
    facility_projection_x(i,:)=[rectangeles(1,i,1),rectangeles(3,i,1)];
    facility_projection_y(i,:)=[rectangeles(2,i,2),rectangeles(4,i,2)];
end
%*****************************PROJECTIONS**********************************
%***************************DETECT COLLISION*******************************
%wall_force(2,2)=0;
for i=1:numpoints
    if cancel_val==1
        break;
    end
    collision_data=cell(1,numpoints+1);%LAST COLLISION IS FOR FACTORY WALLS
    collision_counter=1;
    %************************PEOJECTIONS PEOCESSING************************
    target_facility_projection_x=[rectangeles(1,i,1),rectangeles(3,i,1)];
    target_facility_projection_y=[rectangeles(2,i,2),rectangeles(4,i,2)];
    %************************PEOJECTIONS PEOCESSING************************
    temp_inside_edges=non_mesh_inrectangle(facility_projection_x,facility_projection_y,target_facility_projection_x,target_facility_projection_y);
    if cancel_val==1
        break;
    end
    for j=1:numpoints
        temp_inside_posibility=temp_inside_edges(:,j);
        if sum(temp_inside_posibility)==1
            if temp_inside_posibility(1)==1
                x=(rectangeles(3,i,1)+rectangeles(1,j,1))/2;
                y=(rectangeles(3,i,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(3,i,1)+rectangeles(1,j,1))/2,(rectangeles(3,i,2)+rectangeles(1,j,2))/2);
                collision_surface=abs(rectangeles(3,i,1)-rectangeles(1,j,1))*abs(rectangeles(3,i,2)-rectangeles(1,j,2));
                %collision_data:FOURTH INDEX MEANS WHAT KIND OF COLLISION
                %IS OCCURRED:1 MEANS ONE POINT IS INSIDE SO IT'S FORCE MUST
                %BE DIVIDED BY TWO-2 MEANS THAT TWO POINTS ARE INSIDE SO
                %FORCE MUST BE EXERTED TO ANOTHER ELEMENT TOO-4 IS SIMILAR
                %TO 2 BECAUSE IT CAN'T BE IDENTIFIED BY ANOTHER ELEMENT.
                collision_data{1,collision_counter}=[x,y,collision_surface,1,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(2)==1
                x=(rectangeles(4,i,1)+rectangeles(2,j,1))/2;
                y=(rectangeles(4,i,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(4,i,1)+rectangeles(2,j,1))/2,(rectangeles(4,i,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(4,i,1)-rectangeles(2,j,1))*abs(rectangeles(4,i,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,2,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(3)==1
                x=(rectangeles(1,i,1)+rectangeles(3,j,1))/2;
                y=(rectangeles(1,i,2)+rectangeles(3,j,2))/2;
                %scatter((rectangeles(1,i,1)+rectangeles(3,j,1))/2,(rectangeles(1,i,2)+rectangeles(3,j,2))/2);
                collision_surface=abs(rectangeles(1,i,1)-rectangeles(3,j,1))*abs(rectangeles(1,i,2)-rectangeles(3,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,3,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(4)==1
                x=(rectangeles(2,i,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(4,j,2))/2;
                %scatter((rectangeles(2,i,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(4,j,2))/2);
                collision_surface=abs(rectangeles(2,i,1)-rectangeles(4,j,1))*abs(rectangeles(2,i,2)-rectangeles(4,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,4,j];
                collision_counter=collision_counter+1;
            end
        elseif sum(temp_inside_posibility)==2
            if temp_inside_posibility(1)==1 && temp_inside_posibility(2)==1
                x=(rectangeles(4,i,1)+rectangeles(1,j,1))/2;
                y=(rectangeles(1,j,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(4,i,1)+rectangeles(1,j,1))/2,(rectangeles(1,j,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(4,i,1)-rectangeles(1,j,1))*abs(rectangeles(1,j,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,5,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(2)==1 && temp_inside_posibility(3)==1
                x=(rectangeles(2,j,1)+rectangeles(3,j,1))/2;
                y=(rectangeles(1,i,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(2,j,1)+rectangeles(3,j,1))/2,(rectangeles(1,i,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(2,j,1)-rectangeles(3,j,1))*abs(rectangeles(1,i,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,6,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(3)==1 && temp_inside_posibility(4)==1
                x=(rectangeles(1,i,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(3,j,2)+rectangeles(4,j,2))/2;
                %scatter((rectangeles(1,i,1)+rectangeles(4,j,1))/2,(rectangeles(3,j,2)+rectangeles(4,j,2))/2);
                collision_surface=abs(rectangeles(1,i,1)-rectangeles(4,j,1))*abs(rectangeles(3,j,2)-rectangeles(4,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,7,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(4)==1 && temp_inside_posibility(1)==1
                x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=abs(rectangeles(1,j,1)-rectangeles(4,j,1))*abs(rectangeles(2,i,2)-rectangeles(1,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,8,j];
                collision_counter=collision_counter+1;
            end
        elseif sum(temp_inside_posibility)==4 && j~=i
            x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
            y=(rectangeles(1,j,2)+rectangeles(2,j,2))/2;
            %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(1,j,2)+rectangeles(2,j,2))/2);
            collision_surface=abs(rectangeles(1,j,1)-rectangeles(4,j,1))*abs(rectangeles(1,j,2)-rectangeles(2,j,2));
            collision_data{1,collision_counter}=[x,y,collision_surface,9,j];
            collision_counter=collision_counter+1;
        elseif sum(temp_inside_posibility)==0 && j~=i
            if facility_projection_x(j,1)>target_facility_projection_x(1,1) && facility_projection_x(j,2)<target_facility_projection_x(1,2) && facility_projection_y(j,1)<target_facility_projection_y(1,1) && facility_projection_y(j,2)>target_facility_projection_y(1,2)
                x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(1,i,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=(abs(rectangeles(1,j,1)+rectangeles(4,j,1))*abs(rectangeles(2,i,2)+rectangeles(1,i,2)))/5;
                collision_data{1,collision_counter}=[x,y,collision_surface,10,j];
                collision_counter=collision_counter+1;
            elseif facility_projection_y(j,1)>target_facility_projection_y(1,1) && facility_projection_y(j,2)<target_facility_projection_y(1,2) && facility_projection_x(j,1)<target_facility_projection_x(1,1) && facility_projection_x(j,2)>target_facility_projection_x(1,2)
                x=(rectangeles(1,i,1)+rectangeles(4,i,1))/2;
                y=(rectangeles(2,j,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=(abs(rectangeles(1,i,1)+rectangeles(4,i,1))*abs(rectangeles(2,j,2)+rectangeles(1,j,2)))/5;
                collision_data{1,collision_counter}=[x,y,collision_surface,11,j];
                collision_counter=collision_counter+1;
            end
        end
    end
   %***************************FACTORY WALLS*******************************
    temp_inside_edges=inpolygon(rectangeles(:,i,1),rectangeles(:,i,2),[wall(1,1),wall(1,1),wall(1,1)+wall(2,1),wall(1,1)+wall(2,1)],[wall(1,2)+wall(2,2),wall(1,2),wall(1,2),wall(1,2)+wall(2,2)]);
        if sum(temp_inside_edges)~=4
            if sum(temp_inside_edges)==1
                if temp_inside_edges(1,1)==1
                    x=(rectangeles(1,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(1,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(1,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(1,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,1,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(2,1)==1
                    x=(rectangeles(2,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(2,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(2,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(2,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,2,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(3,1)==1
                    x=(rectangeles(3,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(3,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(3,i,1)-(wall(1,1)))*abs(rectangeles(3,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,3,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(4,1)==1
                    x=(rectangeles(4,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(4,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(4,i,1)-(wall(1,1)))*abs(rectangeles(4,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,4,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                end
            elseif sum(temp_inside_edges)==2
                if temp_inside_edges(1,1)==1 && temp_inside_edges(2,1)==1
                    x=(rectangeles(1,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(1,i,2)+(rectangeles(2,i,2)))/2;
                    collision_surface=abs(rectangeles(1,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(1,i,2)-(rectangeles(2,i,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,5,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(2,1)==1 && temp_inside_edges(3,1)==1
                    x=(rectangeles(2,i,1)+rectangeles(3,i,1))/2;
                    y=(rectangeles(2,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(2,i,1)-rectangeles(3,i,1))*abs(rectangeles(2,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,6,0];
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(3,1)==1 && temp_inside_edges(4,1)==1
                    x=(rectangeles(3,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(3,i,2)+(rectangeles(4,i,2)))/2;
                    collision_surface=abs(rectangeles(3,i,1)-wall(1,1))*abs(rectangeles(3,i,2)-rectangeles(4,i,2));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,7,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(4,1)==1 && temp_inside_edges(1,1)==1
                    x=(rectangeles(4,i,1)+rectangeles(1,i,1))/2;
                    y=(rectangeles(4,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(4,i,1)-rectangeles(1,i,1))*abs(rectangeles(4,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,8,0];
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                end
            elseif sum(temp_inside_edges)==0
                x=(wall(1,1)+wall(2,1))/2;
                y=(wall(1,2)+wall(2,2))/2;
                collision_surface=rectangles_surfaces(i)/10;
                %scatter(x,y);
                collision_data{1,collision_counter}=[x,y,collision_surface,9,0];
            end
        end
   %***************************FACTORY WALLS*******************************
   %****************************MOVE AMOUNT********************************
    if cancel_val==1
        break;
    end
    for j=1:numpoints
        if isempty(collision_data{j})==0
            if collision_data{j}(4)<5 && collision_data{j}(5)~=0
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            elseif collision_data{j}(4)==5 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=0;
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==6 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=0;
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==7 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=0;
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==8 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=0;
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==9 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(5)==0
                total_x=(collision_data{j}(1)-points(i,1))*collision_data{j}(3);
                total_y=(collision_data{j}(2)-points(i,2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            elseif collision_data{j}(4)==10 || collision_data{j}(4)==11
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            end
        else
            break;
        end
    end
    clear collision_data
    %collision_counter=1;
   %****************************MOVE AMOUNT********************************
end
%***************************DETECT COLLISION*******************************
%******************************MOVE POINT**********************************
for i=1:numpoints
    if abs(rectangles_move(i,1))<0.0001
        rectangles_move(i,1)=0.001*sign(rectangles_move(i,1));
    end
    if abs(rectangles_move(i,2))<0.0001
        rectangles_move(i,2)=0.001*sign(rectangles_move(i,2));
    end
end
if max(max(rectangles_move(:,:)))==0 && status==3
    break;
end
for i=1:numpoints
    points(i,:)=points(i,:)+4*rectangles_move(i,:);
    %max(max(rectangles_move));
    rectangles_move(i,:)=[0,0];
end
if cancel_val==1
    break;
end
%******************************MOVE POINT**********************************
%******************************MOVE WALLS**********************************
%{
if status==0
    if max(max(wall_force))<0.001 && max(max(wall_force))~=0
        for o=1:2
            for p=1:2
                coverted_wall_force(o,p)=0.001*wall_force(o,p)/max(max(wall_force));
            end
        end
    else
        for o=1:2
            for p=1:2
                coverted_wall_force(o,p)=wall_force(o,p);
            end
        end
    end
    wall(1,1)=wall(1,1)-0.7*coverted_wall_force(1,1);
    wall(2,1)=wall(2,1)+0.7*coverted_wall_force(1,1);

    wall(1,2)=wall(1,2)-0.7*coverted_wall_force(2,1);
    wall(2,2)=wall(2,2)+0.7*coverted_wall_force(2,1);

    wall(2,1)=wall(2,1)+0.7*coverted_wall_force(1,2);

    wall(2,2)=wall(2,2)+0.7*coverted_wall_force(2,2);
    if max(wall_force)==0
        if num_adjustment>=max_num_adjustment
            break;
        end
        status=1;
    end
elseif status==1
    if wall_force(1,1)==min(min(wall_force))
        wall(1,1)=wall(1,1)+0.005*wall(2,1);
        wall(2,1)=wall(2,1)-0.005*wall(2,1);
    elseif wall_force(1,2)==min(min(wall_force))
        wall(2,1)=wall(2,1)-0.005*wall(2,1);
    elseif wall_force(2,1)==min(min(wall_force))
        wall(1,2)=wall(1,2)+0.005*wall(2,2);
        wall(2,2)=wall(2,2)-0.005*wall(2,2);
    elseif wall_force(2,2)==min(min(wall_force))
        wall(2,2)=wall(2,2)-0.005*wall(2,2);
    end
    free_space=(wall(2,1))*(wall(2,2))-sum(rectangles_surfaces);
    if free_space<0.005*sum(rectangles_surfaces) && num_adjustment<max_num_adjustment
        status=0;
        num_adjustment=num_adjustment+1;
    end
end
clear wall_force
%}
%******************************MOVE WALLS**********************************
end
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed size^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed ratio^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if fixed_ratio_radio==1
numpoints=size(relation_matrix_val,1);
initpoints=rand(numpoints,2);
points=initpoints;
for i=1:numpoints
    for j=1:numpoints
        relation_matrix_val(j,i)=relation_matrix_val(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%new
recdatax=size_table(1,:);
recdatay=size_table(2,:);
rectangeles_init_data=([recdatax;recdatay]);%first row is lengths(x) of rectangels and second is widths(y)
rectangles_surfaces(numpoints)=0;
for i=1:numpoints
    rectangles_surfaces(i)=rectangeles_init_data(1,i)*rectangeles_init_data(2,i);
end

if sum(rectangles_surfaces)>1
    for i=1:numpoints
        rectangeles_init_data(1,i)=rectangeles_init_data(1,i)/sqrt(sum(rectangles_surfaces));
        rectangeles_init_data(2,i)=rectangeles_init_data(2,i)/sqrt(sum(rectangles_surfaces));
    end
end

for i=1:numpoints
    rectangles_surfaces(i)=rectangeles_init_data(1,i)*rectangeles_init_data(2,i);
end

%COMPASS
%    1        4
%     --------
%    |        |
%    |        |
%    |        |
%     --------
%    2        3
%%%%%%%%%%%%%%%%%%%%%new

score_dist=0;
%score_surf=0;
%total_score=0;
%for test=1:test_number
%%%%%%%%%%%%%%%%%%%%%new
rectangeles(4,numpoints,2)=0;%index 1 is compass index 2 is point index 3 is x or y
z(4,numpoints)=0;
rectangles_move(numpoints,2)=0;
wall_init=sum(rectangeles_init_data,2)/0.1;
wall=[0,0;0,0];%START_X-START_Y-WIDTH-HEIGHT
if tempxy>=1
    wall(2,1)=1;
    wall(2,2)=1/tempxy;
end
if tempxy<1
    wall(2,1)=tempxy;
    wall(2,2)=1;
end
wall_force(2,2)=0;%FIRST ROW IS IN HORIZENTAL DIRECTION(X) FIRST COLUMN IN
%FIRST ROW IS LEFT WALL AND SECOND COLUMN IS RIGHT WALL
%SECOND ROW IS IN VERTICAL DIRECTION(Y) FIRST COLUMN IS LOWER WALL AND
%SECOND COLUMN IS UPPER WALL
coverted_wall_force(2,2)=0;
facility_projection_x(numpoints,2)=0;
facility_projection_y(numpoints,2)=0;
target_facility_projection_x(1,2)=0;
target_facility_projection_y(1,2)=0;
status=0;%ZERO MEANS LARGING-ONE MEANS SHRINKING
max_num_adjustment=1;
num_adjustment=0;
%%%%%%%%%%%%%%%%%%%%%new

step=0.001;
totalx=0;
totaly=0;
iteration_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************MOVE POINTS*******************************
while iteration_p <200
for i=1:numpoints
    for j=1:numpoints
        if i~=j
            totalx=totalx+((points(j,1)-points(i,1)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
            totaly=totaly+((points(j,2)-points(i,2)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
        end
    end
    points(i,1)=points(i,1)+totalx*step;
    points(i,2)=points(i,2)+totaly*step;
    totalx=0;
    totaly=0;
end
refresh
pause(0.0001)
scatter(point_axis,points(:,1),points(:,2),20,[iteration_p/200,1-(iteration_p/200),0])
hold on;
iteration_p=iteration_p+1;
if cancel_val==1
    break;
end
end
%********************************MOVE POINTS*******************************
%********************************NORMALIZE*********************************
move=min(points);
for i=1:numpoints

            points(i,1)=points(i,1)-move(1,1);
            points(i,2)=points(i,2)-move(1,2);

end
normalize_val=max(points);
for i=1:numpoints
        points(i,1)=points(i,1)/normalize_val(1,1);
        points(i,2)=points(i,2)/normalize_val(1,2);
end
%********************************NORMALIZE*********************************

scatter(point_axis,points(:,1),points(:,2),24,[0,0,1])
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',point_axis)
end
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:10000
if cancel_val==1
    break;
end
for i=1:numpoints
    rectangeles(1,i,1)=points(i,1)-rectangeles_init_data(1,i)/2;
    rectangeles(2,i,1)=points(i,1)-rectangeles_init_data(1,i)/2;
    rectangeles(3,i,1)=points(i,1)+rectangeles_init_data(1,i)/2;
    rectangeles(4,i,1)=points(i,1)+rectangeles_init_data(1,i)/2;
    
    rectangeles(1,i,2)=points(i,2)+rectangeles_init_data(2,i)/2;
    rectangeles(2,i,2)=points(i,2)-rectangeles_init_data(2,i)/2;
    rectangeles(3,i,2)=points(i,2)-rectangeles_init_data(2,i)/2;
    rectangeles(4,i,2)=points(i,2)+rectangeles_init_data(2,i)/2;
    for j=1:4
        z(j,i)=i;
    end
end
cla
patch(rectangeles(:,:,1),rectangeles(:,:,2),z,'parent',surface_axis)
rectangle('Position',[wall(1,1),wall(1,2),wall(2,1),wall(2,2)],'EdgeColor',[1,0,0],'parent',surface_axis)
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',surface_axis)
end
pause(0.000001);
%*****************************PROJECTIONS**********************************
for i=1:numpoints
    facility_projection_x(i,:)=[rectangeles(1,i,1),rectangeles(3,i,1)];
    facility_projection_y(i,:)=[rectangeles(2,i,2),rectangeles(4,i,2)];
end
%*****************************PROJECTIONS**********************************
%***************************DETECT COLLISION*******************************
wall_force(2,2)=0;
for i=1:numpoints
    if cancel_val==1
        break;
    end
    collision_data=cell(1,numpoints+1);%LAST COLLISION IS FOR FACTORY WALLS
    collision_counter=1;
    %************************PEOJECTIONS PEOCESSING************************
    target_facility_projection_x=[rectangeles(1,i,1),rectangeles(3,i,1)];
    target_facility_projection_y=[rectangeles(2,i,2),rectangeles(4,i,2)];
    %************************PEOJECTIONS PEOCESSING************************
    temp_inside_edges=non_mesh_inrectangle(facility_projection_x,facility_projection_y,target_facility_projection_x,target_facility_projection_y);
    if cancel_val==1
        break;
    end
    for j=1:numpoints
        temp_inside_posibility=temp_inside_edges(:,j);
        if sum(temp_inside_posibility)==1
            if temp_inside_posibility(1)==1
                x=(rectangeles(3,i,1)+rectangeles(1,j,1))/2;
                y=(rectangeles(3,i,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(3,i,1)+rectangeles(1,j,1))/2,(rectangeles(3,i,2)+rectangeles(1,j,2))/2);
                collision_surface=abs(rectangeles(3,i,1)-rectangeles(1,j,1))*abs(rectangeles(3,i,2)-rectangeles(1,j,2));
                %collision_data:FOURTH INDEX MEANS WHAT KIND OF COLLISION
                %IS OCCURRED:1 MEANS ONE POINT IS INSIDE SO IT'S FORCE MUST
                %BE DIVIDED BY TWO-2 MEANS THAT TWO POINTS ARE INSIDE SO
                %FORCE MUST BE EXERTED TO ANOTHER ELEMENT TOO-4 IS SIMILAR
                %TO 2 BECAUSE IT CAN'T BE IDENTIFIED BY ANOTHER ELEMENT.
                collision_data{1,collision_counter}=[x,y,collision_surface,1,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(2)==1
                x=(rectangeles(4,i,1)+rectangeles(2,j,1))/2;
                y=(rectangeles(4,i,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(4,i,1)+rectangeles(2,j,1))/2,(rectangeles(4,i,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(4,i,1)-rectangeles(2,j,1))*abs(rectangeles(4,i,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,2,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(3)==1
                x=(rectangeles(1,i,1)+rectangeles(3,j,1))/2;
                y=(rectangeles(1,i,2)+rectangeles(3,j,2))/2;
                %scatter((rectangeles(1,i,1)+rectangeles(3,j,1))/2,(rectangeles(1,i,2)+rectangeles(3,j,2))/2);
                collision_surface=abs(rectangeles(1,i,1)-rectangeles(3,j,1))*abs(rectangeles(1,i,2)-rectangeles(3,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,3,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(4)==1
                x=(rectangeles(2,i,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(4,j,2))/2;
                %scatter((rectangeles(2,i,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(4,j,2))/2);
                collision_surface=abs(rectangeles(2,i,1)-rectangeles(4,j,1))*abs(rectangeles(2,i,2)-rectangeles(4,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,4,j];
                collision_counter=collision_counter+1;
            end
        elseif sum(temp_inside_posibility)==2
            if temp_inside_posibility(1)==1 && temp_inside_posibility(2)==1
                x=(rectangeles(4,i,1)+rectangeles(1,j,1))/2;
                y=(rectangeles(1,j,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(4,i,1)+rectangeles(1,j,1))/2,(rectangeles(1,j,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(4,i,1)-rectangeles(1,j,1))*abs(rectangeles(1,j,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,5,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(2)==1 && temp_inside_posibility(3)==1
                x=(rectangeles(2,j,1)+rectangeles(3,j,1))/2;
                y=(rectangeles(1,i,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(2,j,1)+rectangeles(3,j,1))/2,(rectangeles(1,i,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(2,j,1)-rectangeles(3,j,1))*abs(rectangeles(1,i,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,6,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(3)==1 && temp_inside_posibility(4)==1
                x=(rectangeles(1,i,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(3,j,2)+rectangeles(4,j,2))/2;
                %scatter((rectangeles(1,i,1)+rectangeles(4,j,1))/2,(rectangeles(3,j,2)+rectangeles(4,j,2))/2);
                collision_surface=abs(rectangeles(1,i,1)-rectangeles(4,j,1))*abs(rectangeles(3,j,2)-rectangeles(4,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,7,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(4)==1 && temp_inside_posibility(1)==1
                x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=abs(rectangeles(1,j,1)-rectangeles(4,j,1))*abs(rectangeles(2,i,2)-rectangeles(1,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,8,j];
                collision_counter=collision_counter+1;
            end
        elseif sum(temp_inside_posibility)==4 && j~=i
            x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
            y=(rectangeles(1,j,2)+rectangeles(2,j,2))/2;
            %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(1,j,2)+rectangeles(2,j,2))/2);
            collision_surface=abs(rectangeles(1,j,1)-rectangeles(4,j,1))*abs(rectangeles(1,j,2)-rectangeles(2,j,2));
            collision_data{1,collision_counter}=[x,y,collision_surface,9,j];
            collision_counter=collision_counter+1;
        elseif sum(temp_inside_posibility)==0 && j~=i
            if facility_projection_x(j,1)>target_facility_projection_x(1,1) && facility_projection_x(j,2)<target_facility_projection_x(1,2) && facility_projection_y(j,1)<target_facility_projection_y(1,1) && facility_projection_y(j,2)>target_facility_projection_y(1,2)
                x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(1,i,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=(abs(rectangeles(1,j,1)+rectangeles(4,j,1))*abs(rectangeles(2,i,2)+rectangeles(1,i,2)))/5;
                collision_data{1,collision_counter}=[x,y,collision_surface,10,j];
                collision_counter=collision_counter+1;
            elseif facility_projection_y(j,1)>target_facility_projection_y(1,1) && facility_projection_y(j,2)<target_facility_projection_y(1,2) && facility_projection_x(j,1)<target_facility_projection_x(1,1) && facility_projection_x(j,2)>target_facility_projection_x(1,2)
                x=(rectangeles(1,i,1)+rectangeles(4,i,1))/2;
                y=(rectangeles(2,j,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=(abs(rectangeles(1,i,1)+rectangeles(4,i,1))*abs(rectangeles(2,j,2)+rectangeles(1,j,2)))/5;
                collision_data{1,collision_counter}=[x,y,collision_surface,11,j];
                collision_counter=collision_counter+1;
            end
        end
    end
   %***************************FACTORY WALLS*******************************
    temp_inside_edges=inpolygon(rectangeles(:,i,1),rectangeles(:,i,2),[wall(1,1),wall(1,1),wall(1,1)+wall(2,1),wall(1,1)+wall(2,1)],[wall(1,2)+wall(2,2),wall(1,2),wall(1,2),wall(1,2)+wall(2,2)]);
        if sum(temp_inside_edges)~=4
            if sum(temp_inside_edges)==1
                if temp_inside_edges(1,1)==1
                    x=(rectangeles(1,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(1,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(1,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(1,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,1,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(2,1)==1
                    x=(rectangeles(2,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(2,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(2,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(2,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,2,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(3,1)==1
                    x=(rectangeles(3,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(3,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(3,i,1)-(wall(1,1)))*abs(rectangeles(3,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,3,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(4,1)==1
                    x=(rectangeles(4,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(4,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(4,i,1)-(wall(1,1)))*abs(rectangeles(4,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,4,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                end
            elseif sum(temp_inside_edges)==2
                if temp_inside_edges(1,1)==1 && temp_inside_edges(2,1)==1
                    x=(rectangeles(1,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(1,i,2)+(rectangeles(2,i,2)))/2;
                    collision_surface=abs(rectangeles(1,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(1,i,2)-(rectangeles(2,i,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,5,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(2,1)==1 && temp_inside_edges(3,1)==1
                    x=(rectangeles(2,i,1)+rectangeles(3,i,1))/2;
                    y=(rectangeles(2,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(2,i,1)-rectangeles(3,i,1))*abs(rectangeles(2,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,6,0];
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(3,1)==1 && temp_inside_edges(4,1)==1
                    x=(rectangeles(3,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(3,i,2)+(rectangeles(4,i,2)))/2;
                    collision_surface=abs(rectangeles(3,i,1)-wall(1,1))*abs(rectangeles(3,i,2)-rectangeles(4,i,2));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,7,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(4,1)==1 && temp_inside_edges(1,1)==1
                    x=(rectangeles(4,i,1)+rectangeles(1,i,1))/2;
                    y=(rectangeles(4,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(4,i,1)-rectangeles(1,i,1))*abs(rectangeles(4,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,8,0];
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                end
            elseif sum(temp_inside_edges)==0
                x=(wall(1,1)+wall(2,1))/2;
                y=(wall(1,2)+wall(2,2))/2;
                collision_surface=rectangles_surfaces(i)/10;
                %scatter(x,y);
                collision_data{1,collision_counter}=[x,y,collision_surface,9,0];
            end
        end
   %***************************FACTORY WALLS*******************************
   %****************************MOVE AMOUNT********************************
    if cancel_val==1
        break;
    end
    for j=1:numpoints
        if isempty(collision_data{j})==0
            if collision_data{j}(4)<5 && collision_data{j}(5)~=0
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            elseif collision_data{j}(4)==5 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=0;
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==6 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=0;
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==7 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=0;
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==8 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=0;
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==9 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(5)==0
                total_x=(collision_data{j}(1)-points(i,1))*collision_data{j}(3);
                total_y=(collision_data{j}(2)-points(i,2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            elseif collision_data{j}(4)==10 || collision_data{j}(4)==11
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            end
        else
            break;
        end
    end
    clear collision_data
    %collision_counter=1;
   %****************************MOVE AMOUNT********************************
end
%***************************DETECT COLLISION*******************************
%******************************MOVE POINT**********************************
for i=1:numpoints
    if abs(rectangles_move(i,1))<0.0001
        rectangles_move(i,1)=0.001*sign(rectangles_move(i,1));
    end
    if abs(rectangles_move(i,2))<0.0001
        rectangles_move(i,2)=0.001*sign(rectangles_move(i,2));
    end
end
if max(max(rectangles_move(:,:)))==0 && status==3
    break;
end
for i=1:numpoints
    points(i,:)=points(i,:)+4*rectangles_move(i,:);
    %max(max(rectangles_move));
    rectangles_move(i,:)=[0,0];
end
%******************************MOVE POINT**********************************
%******************************MOVE WALLS**********************************
if status==0
    if max(max(wall_force))<0.001 && max(max(wall_force))~=0
        for o=1:2
            for p=1:2
                coverted_wall_force(o,p)=0.001*wall_force(o,p)/max(max(wall_force));
            end
        end
    else
        for o=1:2
            for p=1:2
                coverted_wall_force(o,p)=wall_force(o,p);
            end
        end
    end
    wall(1,1)=wall(1,1)-0.5*sum(sum(coverted_wall_force));
    wall(2,1)=wall(2,1)+0.5*sum(sum(coverted_wall_force));

    wall(1,2)=wall(1,2)-0.5*sum(sum(coverted_wall_force));
    wall(2,2)=wall(2,2)+0.5*sum(sum(coverted_wall_force));

    wall(2,1)=wall(2,1)+0.5*sum(sum(coverted_wall_force));

    wall(2,2)=wall(2,2)+0.5*sum(sum(coverted_wall_force));
    if max(wall_force)==0
        if num_adjustment>=max_num_adjustment
            break;
        end
        status=1;
    end
elseif status==1
    %if wall_force(1,1)==min(min(wall_force))
        wall(1,1)=wall(1,1)+0.005*wall(2,1);
        wall(2,1)=wall(2,1)-0.005*wall(2,1);
    %elseif wall_force(1,2)==min(min(wall_force))
        wall(2,1)=wall(2,1)-0.005*wall(2,1);
    %elseif wall_force(2,1)==min(min(wall_force))
        wall(1,2)=wall(1,2)+0.005*wall(2,2);
        wall(2,2)=wall(2,2)-0.005*wall(2,2);
    %elseif wall_force(2,2)==min(min(wall_force))
        wall(2,2)=wall(2,2)-0.005*wall(2,2);
    %end
    free_space=(wall(2,1))*(wall(2,2))-sum(rectangles_surfaces);
    if free_space<0.005*sum(rectangles_surfaces) && num_adjustment<max_num_adjustment
        status=0;
        num_adjustment=num_adjustment+1;
    end
end
clear wall_force
%******************************MOVE WALLS**********************************
if cancel_val==1
    break;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^Fixed ratio^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Flexible^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if flexible_radio==1
numpoints=size(relation_matrix_val,1);
initpoints=rand(numpoints,2);
points=initpoints;
for i=1:numpoints
    for j=1:numpoints
        relation_matrix_val(j,i)=relation_matrix_val(i,j);
    end
end
%%%%%%%%%%%%%%%%%%%%%new
recdatax=size_table(1,:);
recdatay=size_table(2,:);
rectangeles_init_data=([recdatax;recdatay]);%first row is lengths(x) of rectangels and second is widths(y)
rectangles_surfaces(numpoints)=0;
for i=1:numpoints
    rectangles_surfaces(i)=rectangeles_init_data(1,i)*rectangeles_init_data(2,i);
end

if sum(rectangles_surfaces)>1
    for i=1:numpoints
        rectangeles_init_data(1,i)=rectangeles_init_data(1,i)/sqrt(sum(rectangles_surfaces));
        rectangeles_init_data(2,i)=rectangeles_init_data(2,i)/sqrt(sum(rectangles_surfaces));
    end
end

for i=1:numpoints
    rectangles_surfaces(i)=rectangeles_init_data(1,i)*rectangeles_init_data(2,i);
end

%COMPASS
%    1        4
%     --------
%    |        |
%    |        |
%    |        |
%     --------
%    2        3
%%%%%%%%%%%%%%%%%%%%%new

score_dist=0;
%score_surf=0;
%total_score=0;
%for test=1:test_number
%%%%%%%%%%%%%%%%%%%%%new
rectangeles(4,numpoints,2)=0;%index 1 is compass index 2 is point index 3 is x or y
z(4,numpoints)=0;
rectangles_move(numpoints,2)=0;
wall_init=sum(rectangeles_init_data,2)/0.1;
wall=[-1,-1;3,3];%START_X-START_Y-WIDTH-HEIGHT
wall_force(2,2)=0;%FIRST ROW IS IN HORIZENTAL DIRECTION(X) FIRST COLUMN IN
%FIRST ROW IS LEFT WALL AND SECOND COLUMN IS RIGHT WALL
%SECOND ROW IS IN VERTICAL DIRECTION(Y) FIRST COLUMN IS LOWER WALL AND
%SECOND COLUMN IS UPPER WALL
coverted_wall_force(2,2)=0;
facility_projection_x(numpoints,2)=0;
facility_projection_y(numpoints,2)=0;
target_facility_projection_x(1,2)=0;
target_facility_projection_y(1,2)=0;
status=0;%ZERO MEANS LARGING-ONE MEANS SHRINKING
max_num_adjustment=1;
num_adjustment=0;
%%%%%%%%%%%%%%%%%%%%%new

step=0.001;
totalx=0;
totaly=0;
iteration_p=1;
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%********************************MOVE POINTS*******************************
while iteration_p <200
for i=1:numpoints
    for j=1:numpoints
        if i~=j
            totalx=totalx+((points(j,1)-points(i,1)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
            totaly=totaly+((points(j,2)-points(i,2)))*relation_matrix_val(j,i)*(-1*(pdist([points(j,:);points(i,:)])-5)/(1*(pdist([points(j,:);points(i,:)]))^2+1));
        end
    end
    points(i,1)=points(i,1)+totalx*step;
    points(i,2)=points(i,2)+totaly*step;
    totalx=0;
    totaly=0;
end
refresh
pause(0.0001)
scatter(point_axis,points(:,1),points(:,2),20,[iteration_p/200,1-(iteration_p/200),0])
hold on;
iteration_p=iteration_p+1;
if cancel_val==1
    break;
end
end
%********************************MOVE POINTS*******************************
%********************************NORMALIZE*********************************
move=min(points);
for i=1:numpoints

            points(i,1)=points(i,1)-move(1,1);
            points(i,2)=points(i,2)-move(1,2);

end
normalize_val=max(points);
for i=1:numpoints
        points(i,1)=points(i,1)/normalize_val(1,1);
        points(i,2)=points(i,2)/normalize_val(1,2);
end
%********************************NORMALIZE*********************************

scatter(point_axis,points(:,1),points(:,2),24,[0,0,1])
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',point_axis)
end
%%%%%%%%%%%%%%%%%%%%%%%%%POINTS ADJUSTED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:10000
if cancel_val==1
    break;
end
for i=1:numpoints
    rectangeles(1,i,1)=points(i,1)-rectangeles_init_data(1,i)/2;
    rectangeles(2,i,1)=points(i,1)-rectangeles_init_data(1,i)/2;
    rectangeles(3,i,1)=points(i,1)+rectangeles_init_data(1,i)/2;
    rectangeles(4,i,1)=points(i,1)+rectangeles_init_data(1,i)/2;
    
    rectangeles(1,i,2)=points(i,2)+rectangeles_init_data(2,i)/2;
    rectangeles(2,i,2)=points(i,2)-rectangeles_init_data(2,i)/2;
    rectangeles(3,i,2)=points(i,2)-rectangeles_init_data(2,i)/2;
    rectangeles(4,i,2)=points(i,2)+rectangeles_init_data(2,i)/2;
    for j=1:4
        z(j,i)=i;
    end
end
cla
patch(rectangeles(:,:,1),rectangeles(:,:,2),z,'parent',surface_axis)
rectangle('Position',[wall(1,1),wall(1,2),wall(2,1),wall(2,2)],'EdgeColor',[1,0,0],'parent',surface_axis)
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i),'parent',surface_axis)
end
pause(0.000001);
%*****************************PROJECTIONS**********************************
for i=1:numpoints
    facility_projection_x(i,:)=[rectangeles(1,i,1),rectangeles(3,i,1)];
    facility_projection_y(i,:)=[rectangeles(2,i,2),rectangeles(4,i,2)];
end
%*****************************PROJECTIONS**********************************
%***************************DETECT COLLISION*******************************
wall_force(2,2)=0;
for i=1:numpoints
    collision_data=cell(1,numpoints+1);%LAST COLLISION IS FOR FACTORY WALLS
    collision_counter=1;
    %************************PEOJECTIONS PEOCESSING************************
    target_facility_projection_x=[rectangeles(1,i,1),rectangeles(3,i,1)];
    target_facility_projection_y=[rectangeles(2,i,2),rectangeles(4,i,2)];
    %************************PEOJECTIONS PEOCESSING************************
    temp_inside_edges=non_mesh_inrectangle(facility_projection_x,facility_projection_y,target_facility_projection_x,target_facility_projection_y);
    if cancel_val==1
        break;
    end
    for j=1:numpoints
        temp_inside_posibility=temp_inside_edges(:,j);
        if sum(temp_inside_posibility)==1
            if temp_inside_posibility(1)==1
                x=(rectangeles(3,i,1)+rectangeles(1,j,1))/2;
                y=(rectangeles(3,i,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(3,i,1)+rectangeles(1,j,1))/2,(rectangeles(3,i,2)+rectangeles(1,j,2))/2);
                collision_surface=abs(rectangeles(3,i,1)-rectangeles(1,j,1))*abs(rectangeles(3,i,2)-rectangeles(1,j,2));
                %collision_data:FOURTH INDEX MEANS WHAT KIND OF COLLISION
                %IS OCCURRED:1 MEANS ONE POINT IS INSIDE SO IT'S FORCE MUST
                %BE DIVIDED BY TWO-2 MEANS THAT TWO POINTS ARE INSIDE SO
                %FORCE MUST BE EXERTED TO ANOTHER ELEMENT TOO-4 IS SIMILAR
                %TO 2 BECAUSE IT CAN'T BE IDENTIFIED BY ANOTHER ELEMENT.
                collision_data{1,collision_counter}=[x,y,collision_surface,1,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(2)==1
                x=(rectangeles(4,i,1)+rectangeles(2,j,1))/2;
                y=(rectangeles(4,i,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(4,i,1)+rectangeles(2,j,1))/2,(rectangeles(4,i,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(4,i,1)-rectangeles(2,j,1))*abs(rectangeles(4,i,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,2,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(3)==1
                x=(rectangeles(1,i,1)+rectangeles(3,j,1))/2;
                y=(rectangeles(1,i,2)+rectangeles(3,j,2))/2;
                %scatter((rectangeles(1,i,1)+rectangeles(3,j,1))/2,(rectangeles(1,i,2)+rectangeles(3,j,2))/2);
                collision_surface=abs(rectangeles(1,i,1)-rectangeles(3,j,1))*abs(rectangeles(1,i,2)-rectangeles(3,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,3,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(4)==1
                x=(rectangeles(2,i,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(4,j,2))/2;
                %scatter((rectangeles(2,i,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(4,j,2))/2);
                collision_surface=abs(rectangeles(2,i,1)-rectangeles(4,j,1))*abs(rectangeles(2,i,2)-rectangeles(4,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,4,j];
                collision_counter=collision_counter+1;
            end
        elseif sum(temp_inside_posibility)==2
            if temp_inside_posibility(1)==1 && temp_inside_posibility(2)==1
                x=(rectangeles(4,i,1)+rectangeles(1,j,1))/2;
                y=(rectangeles(1,j,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(4,i,1)+rectangeles(1,j,1))/2,(rectangeles(1,j,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(4,i,1)-rectangeles(1,j,1))*abs(rectangeles(1,j,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,5,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(2)==1 && temp_inside_posibility(3)==1
                x=(rectangeles(2,j,1)+rectangeles(3,j,1))/2;
                y=(rectangeles(1,i,2)+rectangeles(2,j,2))/2;
                %scatter((rectangeles(2,j,1)+rectangeles(3,j,1))/2,(rectangeles(1,i,2)+rectangeles(2,j,2))/2);
                collision_surface=abs(rectangeles(2,j,1)-rectangeles(3,j,1))*abs(rectangeles(1,i,2)-rectangeles(2,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,6,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(3)==1 && temp_inside_posibility(4)==1
                x=(rectangeles(1,i,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(3,j,2)+rectangeles(4,j,2))/2;
                %scatter((rectangeles(1,i,1)+rectangeles(4,j,1))/2,(rectangeles(3,j,2)+rectangeles(4,j,2))/2);
                collision_surface=abs(rectangeles(1,i,1)-rectangeles(4,j,1))*abs(rectangeles(3,j,2)-rectangeles(4,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,7,j];
                collision_counter=collision_counter+1;
            elseif temp_inside_posibility(4)==1 && temp_inside_posibility(1)==1
                x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=abs(rectangeles(1,j,1)-rectangeles(4,j,1))*abs(rectangeles(2,i,2)-rectangeles(1,j,2));
                collision_data{1,collision_counter}=[x,y,collision_surface,8,j];
                collision_counter=collision_counter+1;
            end
        elseif sum(temp_inside_posibility)==4 && j~=i
            x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
            y=(rectangeles(1,j,2)+rectangeles(2,j,2))/2;
            %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(1,j,2)+rectangeles(2,j,2))/2);
            collision_surface=abs(rectangeles(1,j,1)-rectangeles(4,j,1))*abs(rectangeles(1,j,2)-rectangeles(2,j,2));
            collision_data{1,collision_counter}=[x,y,collision_surface,9,j];
            collision_counter=collision_counter+1;
        elseif sum(temp_inside_posibility)==0 && j~=i
            if facility_projection_x(j,1)>target_facility_projection_x(1,1) && facility_projection_x(j,2)<target_facility_projection_x(1,2) && facility_projection_y(j,1)<target_facility_projection_y(1,1) && facility_projection_y(j,2)>target_facility_projection_y(1,2)
                x=(rectangeles(1,j,1)+rectangeles(4,j,1))/2;
                y=(rectangeles(2,i,2)+rectangeles(1,i,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=(abs(rectangeles(1,j,1)+rectangeles(4,j,1))*abs(rectangeles(2,i,2)+rectangeles(1,i,2)))/5;
                collision_data{1,collision_counter}=[x,y,collision_surface,10,j];
                collision_counter=collision_counter+1;
            elseif facility_projection_y(j,1)>target_facility_projection_y(1,1) && facility_projection_y(j,2)<target_facility_projection_y(1,2) && facility_projection_x(j,1)<target_facility_projection_x(1,1) && facility_projection_x(j,2)>target_facility_projection_x(1,2)
                x=(rectangeles(1,i,1)+rectangeles(4,i,1))/2;
                y=(rectangeles(2,j,2)+rectangeles(1,j,2))/2;
                %scatter((rectangeles(1,j,1)+rectangeles(4,j,1))/2,(rectangeles(2,i,2)+rectangeles(1,j,2))/2);
                collision_surface=(abs(rectangeles(1,i,1)+rectangeles(4,i,1))*abs(rectangeles(2,j,2)+rectangeles(1,j,2)))/5;
                collision_data{1,collision_counter}=[x,y,collision_surface,11,j];
                collision_counter=collision_counter+1;
            end
        end
    end
   %***************************FACTORY WALLS*******************************
    temp_inside_edges=inpolygon(rectangeles(:,i,1),rectangeles(:,i,2),[wall(1,1),wall(1,1),wall(1,1)+wall(2,1),wall(1,1)+wall(2,1)],[wall(1,2)+wall(2,2),wall(1,2),wall(1,2),wall(1,2)+wall(2,2)]);
        if sum(temp_inside_edges)~=4
            if sum(temp_inside_edges)==1
                if temp_inside_edges(1,1)==1
                    x=(rectangeles(1,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(1,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(1,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(1,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,1,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(2,1)==1
                    x=(rectangeles(2,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(2,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(2,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(2,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,2,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(3,1)==1
                    x=(rectangeles(3,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(3,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(3,i,1)-(wall(1,1)))*abs(rectangeles(3,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,3,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(4,1)==1
                    x=(rectangeles(4,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(4,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(4,i,1)-(wall(1,1)))*abs(rectangeles(4,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,4,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                end
            elseif sum(temp_inside_edges)==2
                if temp_inside_edges(1,1)==1 && temp_inside_edges(2,1)==1
                    x=(rectangeles(1,i,1)+(wall(1,1)+wall(2,1)))/2;
                    y=(rectangeles(1,i,2)+(rectangeles(2,i,2)))/2;
                    collision_surface=abs(rectangeles(1,i,1)-(wall(1,1)+wall(2,1)))*abs(rectangeles(1,i,2)-(rectangeles(2,i,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,5,0];
                    wall_force(1,2)=wall_force(1,2)+abs(abs(x-(wall(1,1)+wall(2,1)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(2,1)==1 && temp_inside_edges(3,1)==1
                    x=(rectangeles(2,i,1)+rectangeles(3,i,1))/2;
                    y=(rectangeles(2,i,2)+(wall(1,2)+wall(2,2)))/2;
                    collision_surface=abs(rectangeles(2,i,1)-rectangeles(3,i,1))*abs(rectangeles(2,i,2)-(wall(1,2)+wall(2,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,6,0];
                    wall_force(2,2)=wall_force(2,2)+abs(abs(y-(wall(1,2)+wall(2,2)))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(3,1)==1 && temp_inside_edges(4,1)==1
                    x=(rectangeles(3,i,1)+(wall(1,1)))/2;
                    y=(rectangeles(3,i,2)+(rectangeles(4,i,2)))/2;
                    collision_surface=abs(rectangeles(3,i,1)-wall(1,1))*abs(rectangeles(3,i,2)-rectangeles(4,i,2));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,7,0];
                    wall_force(1,1)=wall_force(1,1)+abs(abs(x-wall(1,1))*(rectangles_surfaces(i)-collision_surface));
                elseif temp_inside_edges(4,1)==1 && temp_inside_edges(1,1)==1
                    x=(rectangeles(4,i,1)+rectangeles(1,i,1))/2;
                    y=(rectangeles(4,i,2)+(wall(1,2)))/2;
                    collision_surface=abs(rectangeles(4,i,1)-rectangeles(1,i,1))*abs(rectangeles(4,i,2)-(wall(1,2)));
                    %scatter(x,y);
                    collision_data{1,collision_counter}=[x,y,collision_surface,8,0];
                    wall_force(2,1)=wall_force(2,1)+abs(abs(y-wall(1,2))*(rectangles_surfaces(i)-collision_surface));
                end
            elseif sum(temp_inside_edges)==0
                x=(wall(1,1)+wall(2,1))/2;
                y=(wall(1,2)+wall(2,2))/2;
                collision_surface=rectangles_surfaces(i)/10;
                %scatter(x,y);
                collision_data{1,collision_counter}=[x,y,collision_surface,9,0];
            end
        end
   %***************************FACTORY WALLS*******************************
   %****************************MOVE AMOUNT********************************
    if cancel_val==1
        break;
    end
    for j=1:numpoints
        if isempty(collision_data{j})==0
            if collision_data{j}(4)<5 && collision_data{j}(5)~=0
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            elseif collision_data{j}(4)==5 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=0;
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==6 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=0;
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==7 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=0;
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==8 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=0;
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(4)==9 && collision_data{j}(5)~=0
                another=collision_data{j}(5);
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                total_x_another=(points(another,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y_another=(points(another,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                if abs(total_x_another) > mean(rectangeles_init_data(:,another)) || abs(total_y_another) > mean(rectangeles_init_data(:,another))
                    total_x_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    total_y_another=max(abs(total_x_another),abs(total_y_another))*mean(rectangeles_init_data(:,another))*sign(total_x_another);
                    %disp('max_reached_another')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y;
                rectangles_move(another,1)=rectangles_move(another,1)+total_x_another;
                rectangles_move(another,2)=rectangles_move(another,2)+total_y_another;
            elseif collision_data{j}(5)==0
                total_x=(collision_data{j}(1)-points(i,1))*collision_data{j}(3);
                total_y=(collision_data{j}(2)-points(i,2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            elseif collision_data{j}(4)==10 || collision_data{j}(4)==11
                total_x=(points(i,1)-collision_data{j}(1))*collision_data{j}(3);
                total_y=(points(i,2)-collision_data{j}(2))*collision_data{j}(3);
                if abs(total_x) > mean(rectangeles_init_data(:,i)) || abs(total_y) > mean(rectangeles_init_data(:,i))
                    total_x=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_x);
                    total_y=max(abs(total_x),abs(total_y))*mean(rectangeles_init_data(:,i))*sign(total_y);
                    %disp('max_reached')
                end
                rectangles_move(i,1)=rectangles_move(i,1)+total_x/2;
                rectangles_move(i,2)=rectangles_move(i,2)+total_y/2;
            end
        else
            break;
        end
    end
    clear collision_data
    %collision_counter=1;
   %****************************MOVE AMOUNT********************************
end
%***************************DETECT COLLISION*******************************
%******************************MOVE POINT**********************************
for i=1:numpoints
    if abs(rectangles_move(i,1))<0.0001
        rectangles_move(i,1)=0.001*sign(rectangles_move(i,1));
    end
    if abs(rectangles_move(i,2))<0.0001
        rectangles_move(i,2)=0.001*sign(rectangles_move(i,2));
    end
end
if max(max(rectangles_move(:,:)))==0 && status==3
    break;
end
for i=1:numpoints
    points(i,:)=points(i,:)+4*rectangles_move(i,:);
    %max(max(rectangles_move));
    rectangles_move(i,:)=[0,0];
end
%******************************MOVE POINT**********************************
%******************************MOVE WALLS**********************************
if status==0
    if max(max(wall_force))<0.001 && max(max(wall_force))~=0
        for o=1:2
            for p=1:2
                coverted_wall_force(o,p)=0.001*wall_force(o,p)/max(max(wall_force));
            end
        end
    else
        for o=1:2
            for p=1:2
                coverted_wall_force(o,p)=wall_force(o,p);
            end
        end
    end
    wall(1,1)=wall(1,1)-release_speed_val*0.03*coverted_wall_force(1,1);
    wall(2,1)=wall(2,1)+release_speed_val*0.03*coverted_wall_force(1,1);

    wall(1,2)=wall(1,2)-release_speed_val*0.03*coverted_wall_force(2,1);
    wall(2,2)=wall(2,2)+release_speed_val*0.03*coverted_wall_force(2,1);

    wall(2,1)=wall(2,1)+release_speed_val*0.03*coverted_wall_force(1,2);

    wall(2,2)=wall(2,2)+release_speed_val*0.03*coverted_wall_force(2,2);
    if max(wall_force)==0
        if num_adjustment>=max_num_adjustment
            break;
        end
        status=1;
    end
elseif status==1
    if wall_force(1,1)==min(min(wall_force))
        wall(1,1)=wall(1,1)+compress_speed_val*0.0004*wall(2,1);
        wall(2,1)=wall(2,1)-compress_speed_val*0.0004*wall(2,1);
    elseif wall_force(1,2)==min(min(wall_force))
        wall(2,1)=wall(2,1)-compress_speed_val*0.0004*wall(2,1);
    elseif wall_force(2,1)==min(min(wall_force))
        wall(1,2)=wall(1,2)+compress_speed_val*0.0004*wall(2,2);
        wall(2,2)=wall(2,2)-compress_speed_val*0.0004*wall(2,2);
    elseif wall_force(2,2)==min(min(wall_force))
        wall(2,2)=wall(2,2)-compress_speed_val*0.0004*wall(2,2);
    end
    free_space=(wall(2,1))*(wall(2,2))-sum(rectangles_surfaces);
    if free_space<0.005*sum(rectangles_surfaces) && num_adjustment<max_num_adjustment
        status=0;
        num_adjustment=num_adjustment+1;
    end
end
clear wall_force
%******************************MOVE WALLS**********************************
if cancel_val==1
    break;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%SURFACE ADJUSTMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Flexible^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%%%%%%%%%%%%%%%%%%%%%%%%%%%Score layout%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***********************Departments distance score*************************
score_dist=0;
for i=1:numpoints
    for j=i+1:numpoints
        score_dist=score_dist+relation_matrix_val(i,j)*(1/(sqrt(((points(i,1)-points(j,1))^2)+((points(i,2)-points(j,2))^2))));
    end
end
%***********************Departments distance score*************************
%************************Surface usage score*******************************
score_surf=1/(wall(2,1)*wall(2,2));
%************************Surface usage score*******************************
total_score=score_dist+score_surf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%Score layout%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
set(cancel,'enable','off')
set(run,'enable','on')
if cancel_val==0
    set(save_results,'enable','on')
    if nonmesh_radio==1
        finished_method='nonmesh';
    elseif mesh_radio==1
        finished_method='mesh';
    end
else
    set(save_results,'enable','off')
end



function output=non_mesh_inrectangle(facility_projection_x,facility_projection_y,target_facility_projection_x,target_facility_projection_y)
numobject=size(facility_projection_x,1);
output(4,numobject)=0;
for i=1:numobject
    if facility_projection_x(i,1)<target_facility_projection_x(1,2) && facility_projection_x(i,1)>target_facility_projection_x(1,1) && facility_projection_y(i,2)<target_facility_projection_y(1,2) && facility_projection_y(i,2)>target_facility_projection_y(1,1)
        output(:,i)=output(:,i)+[1;0;0;0];
    end
    if facility_projection_x(i,1)<target_facility_projection_x(1,2) && facility_projection_x(i,1)>target_facility_projection_x(1,1) && facility_projection_y(i,1)<target_facility_projection_y(1,2) && facility_projection_y(i,1)>target_facility_projection_y(1,1)
        output(:,i)=output(:,i)+[0;1;0;0];
    end
    if facility_projection_x(i,2)<target_facility_projection_x(1,2) && facility_projection_x(i,2)>target_facility_projection_x(1,1) && facility_projection_y(i,1)<target_facility_projection_y(1,2) && facility_projection_y(i,1)>target_facility_projection_y(1,1)
        output(:,i)=output(:,i)+[0;0;1;0];
    end
    if facility_projection_x(i,2)<target_facility_projection_x(1,2) && facility_projection_x(i,2)>target_facility_projection_x(1,1) && facility_projection_y(i,2)<target_facility_projection_y(1,2) && facility_projection_y(i,2)>target_facility_projection_y(1,1)
        output(:,i)=output(:,i)+[0;0;0;1];
    end
end



% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear
close


% --- Executes on button press in load_relation_matrix.
function load_relation_matrix_Callback(hObject, eventdata, handles)
% hObject    handle to load_relation_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global relation_matrix
global relation_matrix_val
[FileName,PathName,FilterIndex] = uigetfile('*.xls');
if FilterIndex~=0
oldFolder = cd(PathName);
[relation_matrix_val,~,raw]=xlsread(FileName);
cd(oldFolder)
set(relation_matrix,'data',relation_matrix_val)
end



% --- Executes on button press in load_size_matrix.
function load_size_matrix_Callback(hObject, eventdata, handles)
% hObject    handle to load_size_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global size_matrix
global size_table
[FileName,PathName,FilterIndex] = uigetfile('*.xls');
if FilterIndex~=0
oldFolder = cd(PathName);
[size_table,~,~]=xlsread(FileName);
cd(oldFolder)
set(size_matrix,'data',size_table)
end

% --- Executes on button press in presentation.
function presentation_Callback(hObject, eventdata, handles)
% hObject    handle to presentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
winopen('v43.WRL')


% --- Executes on button press in save_results.
function save_results_Callback(hObject, eventdata, handles)
% hObject    handle to save_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numpoints
global points
global factoryx
global factoryy
global centeredx
global centeredy
global z
%global total_score
global vertex_data
global finished_method
global sizes
global rectangeles
global wall
[FileName,PathName,FilterIndex] = uiputfile({'*.xls';'*.txt'},'Save results');
if FilterIndex==2 && FilterIndex~=0
oldFolder = cd(PathName);
temp=figure;
switch finished_method
    case 'mesh'
  axis([-0.1 factoryx+0.1 -0.1 factoryy+0.1])
  for i=1:numpoints
      patch(vertex_data{i,1},vertex_data{i,2},z(i))
  end
  rectangle('Position',[0,0,factoryx,factoryy],'EdgeColor',[1,0,0])
  for i=1:numpoints
      text(points(i,1),points(i,2),num2str(i))
  end
  FileName_pic = strrep(FileName, '.txt', '');
  saveas(temp,FileName_pic,'jpg')
  close(temp)
  fid = fopen(strcat(FileName), 'w');
  for i=1:numpoints
      fprintf(fid, '%d:\t\n', i);
      for l=1:2
          if l==1
              fprintf(fid, 'x:\t');
          elseif l==2
              fprintf(fid, 'y:\t');
          end
          for j=1:sizes(i)
              fprintf(fid, '%f\t', vertex_data{i,l}(j,1));
          end
          fprintf(fid, '\n');
      end
  end
  fprintf(fid, 'Factory width:%f\t', factoryx);
  fprintf(fid, 'Factory height:%f\t', factoryy);
  %fprintf(fid, 'Score:%f', total_score);
  fclose(fid);
  cd(oldFolder)
    case 'nonmesh'
  patch(rectangeles(:,:,1),rectangeles(:,:,2),z)
  rectangle('Position',[wall(1,1),wall(1,2),wall(2,1),wall(2,2)],'EdgeColor',[1,0,0])
  for i=1:numpoints
      text(points(i,1),points(i,2),num2str(i))
  end
  FileName_pic = strrep(FileName, '.txt', '');
  saveas(temp,FileName_pic,'jpg')
  close(temp)
  fid = fopen(strcat(FileName), 'w');
  for i=1:numpoints
      fprintf(fid, '%d:\t\n', i);
      for l=1:2
          if l==1
              fprintf(fid, 'x:\t');
          elseif l==2
              fprintf(fid, 'y:\t');
          end
          for j=1:4
              fprintf(fid, '%f\t', rectangeles(j,i,l));
          end
          fprintf(fid, '\n');
      end
  end
  fprintf(fid, 'Factory width:%f\t', wall(2,1));
  fprintf(fid, 'Factory height:%f\t', wall(2,2));
  %fprintf(fid, 'Score:%f', total_score);
  fclose(fid);
  cd(oldFolder)
 end
end
if FilterIndex==1 && FilterIndex~=0
    oldFolder = cd(PathName);
    temp=figure;
    switch finished_method
        case 'mesh'
            axis([-0.1 factoryx+0.1 -0.1 factoryy+0.1])
            for i=1:numpoints
                patch(vertex_data{i,1},vertex_data{i,2},z(i))
            end
            rectangle('Position',[0,0,factoryx,factoryy],'EdgeColor',[1,0,0])
            for i=1:numpoints
                text(points(i,1),points(i,2),num2str(i))
            end
            FileName_pic = strrep(FileName, '.xls', '');
            saveas(temp,FileName_pic,'jpg')
            close(temp)
            
            data_mesh=cell(2*numpoints+1,max(sizes)+2);
            for i=1:numpoints
                data_mesh{(i*2)-1,1}=i;
                data_mesh{(i*2)-1,2}='x';
                data_mesh{(i*2),2}='y';
                for j=1:sizes(i)
                    if isnan(vertex_data{i,1}(j,1))==0
                        data_mesh{(i*2)-1,j+2}=vertex_data{i,1}(j,1);
                        data_mesh{(i*2),j+2}=vertex_data{i,2}(j,1);
                    end
                end
            end
            xlswrite(FileName,data_mesh);
            
            cd(oldFolder)
            
        case 'nonmesh'
            patch(rectangeles(:,:,1),rectangeles(:,:,2),z)
            rectangle('Position',[wall(1,1),wall(1,2),wall(2,1),wall(2,2)],'EdgeColor',[1,0,0])
            for i=1:numpoints
                text(points(i,1),points(i,2),num2str(i))
            end
            FileName_pic = strrep(FileName, '.xls', '');
            saveas(temp,FileName_pic,'jpg')
            close(temp)
            
            data_mesh=cell(2*numpoints+1,4+2);
            for i=1:numpoints
                data_mesh{(i*2)-1,1}=i;
                data_mesh{(i*2)-1,2}='x';
                data_mesh{(i*2),2}='y';
                for j=1:4
                    data_mesh{(i*2)-1,j+2}=rectangeles(j,i,1);
                    data_mesh{(i*2),j+2}=rectangeles(j,i,2);
                end
            end
            xlswrite(FileName,data_mesh);
            
            cd(oldFolder)
            
    end
end
%{
if size(size_table,1)==2
axis([-0.1 factoryx+0.1 -0.1 factoryy+0.1])
patch(centeredx,centeredy,z)
rectangle('Position',[0,0,factoryx,factoryy],'EdgeColor',[1,0,0])
for i=1:numpoints
    text(points(i,1),points(i,2),num2str(i))
end
saveas(temp,FileName,'jpg')
close(temp)
fid = fopen(strcat(FileName,'.txt'), 'w');
for i=1:numpoints
fprintf(fid, '%d:\t', i);
fprintf(fid, '%f\t', centeredx(1,i));
fprintf(fid, '%f\t', centeredy(1,i));
fprintf(fid, '%f\t', centeredx(2,i));
fprintf(fid, '%f\t', centeredy(2,i));
fprintf(fid, '%f\t', centeredx(3,i));
fprintf(fid, '%f\t', centeredy(3,i));
fprintf(fid, '%f\t', centeredx(4,i));
fprintf(fid, '%f\t', centeredy(4,i));
end
fprintf(fid, 'Factory x:%f\t', factoryx);
fprintf(fid, 'Factory y:%f\t', factoryy);
%fprintf(fid, 'Score:%f', total_score);
elseif size(size_table,1)>2

end
fclose(fid);
cd(oldFolder)
%}


% --- Executes when selected object is changed in uipanel8.
function limit_selectionfcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel8 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global fixed_size
global fixed_ratio
global flexible
global fixed_size_radio
global fixed_ratio_radio
global flexible_radio
switch get(eventdata.NewValue,'Tag')
    case 'fixed_size'
        fixed_size_radio=1;
        fixed_ratio_radio=0;
        flexible_radio=0;
        set(fixed_size,'visible','on')
        set(fixed_ratio,'visible','off')
        set(flexible,'visible','off')
    case 'fixed_ratio'
        fixed_size_radio=0;
        fixed_ratio_radio=1;
        flexible_radio=0;
        set(fixed_size,'visible','off')
        set(fixed_ratio,'visible','on')
        set(flexible,'visible','off')
    case 'flexible'
        fixed_size_radio=0;
        fixed_ratio_radio=0;
        flexible_radio=1;
        set(fixed_size,'visible','off')
        set(fixed_ratio,'visible','off')
        set(flexible,'visible','on')
end


% --- Executes when selected object is changed in uipanel8.
function method_selectionfcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel8 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global nonmesh_radio
global mesh_radio
switch get(eventdata.NewValue,'Tag')
    case 'nonmesh_method'
        nonmesh_radio=1;
        mesh_radio=0;
    case 'mesh_method'
        nonmesh_radio=0;
        mesh_radio=1;
end



% --- Executes on button press in fixed_size.
function fixed_size_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixed_size


% --- Executes on button press in fixed_ratio.
function fixed_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fixed_ratio


% --- Executes on button press in flexible.
function flexible_Callback(hObject, eventdata, handles)
% hObject    handle to flexible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of flexible



function fixed_size_x_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_size_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global factoryx
factoryx=str2double(get(hObject,'string'));
if isnan(factoryx)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
end
if factoryx>1 || factoryx<0.1
  errordlg('Factory length can be 0.1 to 1','Bad Input','modal')
  uicontrol(hObject)
end
% Hints: get(hObject,'String') returns contents of fixed_size_x as text
%        str2double(get(hObject,'String')) returns contents of fixed_size_x as a double


% --- Executes during object creation, after setting all properties.
function fixed_size_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixed_size_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fixed_size_y_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_size_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global factoryy
factoryy=str2double(get(hObject,'string'));
if isnan(factoryy)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
end
if factoryy>1 || factoryy<0.1
  errordlg('Factory width can be 0.1 to 1','Bad Input','modal')
  uicontrol(hObject)
end
% Hints: get(hObject,'String') returns contents of fixed_size_y as text
%        str2double(get(hObject,'String')) returns contents of fixed_size_y as a double


% --- Executes during object creation, after setting all properties.
function fixed_size_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixed_size_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fixed_ratio_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to fixed_ratio_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tempxy
tempxy=str2double(get(hObject,'string'));
if isnan(tempxy)
  errordlg('You must enter a numeric value','Bad Input','modal')
  uicontrol(hObject)
end
if tempxy>10 || tempxy<0.1
  errordlg('Length to width ratio can be 0.1 to 10','Bad Input','modal')
  uicontrol(hObject)
end
% Hints: get(hObject,'String') returns contents of fixed_ratio_ratio as text
%        str2double(get(hObject,'String')) returns contents of fixed_ratio_ratio as a double


% --- Executes during object creation, after setting all properties.
function fixed_ratio_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fixed_ratio_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function flexible_release_Callback(hObject, eventdata, handles)
% hObject    handle to flexible_release (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flexible_release
global release_speed_val
release_speed_val=round((get(flexible_release,'value'))*15+5);
set(flexible_release,'value',(release_speed_val-5)/15)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function flexible_release_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flexible_release (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function flexible_compress_Callback(hObject, eventdata, handles)
% hObject    handle to flexible_compress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flexible_compress
global compress_speed_val
compress_speed_val=round((get(flexible_compress,'value'))*15+5);
set(flexible_compress,'value',(compress_speed_val-5)/15)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function flexible_compress_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flexible_compress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cancel_val
cancel_val=1;

% --- Executes on slider movement.
function outer_walls_penalty_Callback(hObject, eventdata, handles)
% hObject    handle to outer_walls_penalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global outer_walls_penalty
global outer_walls_penalty_val
outer_walls_penalty_val=(get(outer_walls_penalty,'value'))+0.5;
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function outer_walls_penalty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outer_walls_penalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function collision_penalty_Callback(hObject, eventdata, handles)
% hObject    handle to collision_penalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global collision_penalty
global collision_penalty_val
collision_penalty_val=(get(collision_penalty,'value'))+0.5;
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function collision_penalty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to collision_penalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in about.
function about_Callback(hObject, eventdata, handles)
% hObject    handle to about (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global about_fig
global page_number
page_number=0
about_fig=figure('menubar','none','toolbar','none','name','About','NumberTitle','off');
annotation('textbox',[0.15 0.15 0.7 0.7],'Color',[0,0,0],'EdgeColor',[1,0.2,0.2],'BackgroundColor',[0.94,0.94,1],'FontSize',19,'string','This GUI is designed for running facility layout collision detection algorithm, loading input data, tunning parameters and saving results. This GUI is designed by Amir Mohammad Esmaieeli Sikaroudi. Email: a.a.easyplot@gmail.com')
h1 = uicontrol(about_fig,'Style','PushButton','Units','normalized',...
               'String','Hints','Position',[.4 .03 .2 .1],'callback', @pushbutton_cb);

function pushbutton_cb(hcbo, eventStruct)
global about_fig
clf(about_fig)




% --- Executes when entered data in editable cell(s) in relation_matrix.
function relation_matrix_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to relation_matrix (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
global relation_matrix
global relation_matrix_val
relation_matrix_val=get(relation_matrix,'data');


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all
FLCD_v5
