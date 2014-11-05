% Function [FD]=FD2D(xd,yd)
% calculates the fractal dimension of a set in the 2D plane. 
% using the box counting method.
%
% INPUT:        x,y           coordinates of the set points
%         
% DEFINITIONS:  maxbox        maximum box size (set at 45% of ROI)
%               (minbox       minimum box size is set at 1 pixels)
%                   
% OUTPUT:       FD            the estimated fractal dimension

% AUTHOR: Vivek Muthurangu, James C Moon, Gaby Captur - UCL
% VERSION: 3.0, 25.10.2014
% EMAIL: gabriella.captur.11@ucl.ac.uk

function [FD]=FD2D(xd,yd)

tic 
 
% convert the cells to a numeric array using the cell2mat function.


coordmat= [xd,yd];


% creating the bounding box.
glob_llx = min(coordmat(:,1)); % lower left corner of the initial box, x-coordinate (ll= lower limits).
glob_lly = min(coordmat(:,2)); % lower left corner of the initial box, y-coordinate (ll= lower limits).
glob_urx = max(coordmat(:,1)); % upper right corner of the initial box, x-coordinate (ul= upper limits).
glob_ury = max(coordmat(:,2)); % upper right corner of the initial box, y-coordinate (ul= upper limits).
glob_width = glob_urx - glob_llx; % width of the bounding box.
glob_height = glob_ury - glob_lly; % height of the bounding box.
    if (glob_width < 49) && (glob_height < 49)
        msgbox('Slice too apical - Reimport and do not segment this slice.','Fractal Analysis Warning','warn')
    end
sz_maxx=glob_width*1.5; % expands the size of the bounding box, by factor 1.5 to ensure complete coverage of the ROI.
sz_maxy=glob_height*1.5;% expands the size of the bounding box, by factor 1.5 to ensure complete coverage of the ROI.
sz_minx=glob_llx-(sz_maxx-glob_urx);% expands the size of the bounding box, by factor 1.5 to ensure complete coverage of the ROI.
sz_miny=glob_lly-(sz_maxy-glob_ury);% expands the size of the bounding box, by factor 1.5 to ensure complete coverage of the ROI.
maxsz = max(glob_width, glob_height); % finds which out of height/width is the largest value for the bounding box.
maxbox = ceil(maxsz*0.45); % the maximum box size will be 45% of the maximum size of the original ROI.

%% ADDED by CTG 
% four fixed points along a diagonal
init_minx  =   sz_minx; 
init_miny  =   sz_miny;
step = ceil(maxbox / 3);

for start_pos=1:4
     sz_minx(start_pos) = floor(init_minx + (start_pos-1) * step);
     sz_miny(start_pos) = floor(init_miny + (start_pos-1) * step);
    for bx_sz=1:maxbox % establishes the number of grid calibres in this grid (iterations).
    
    new_im=[]; % pre-allocates an empty array which will later be populated by px_x/y data elements. 
    [mx,my]=meshgrid(sz_minx(start_pos):bx_sz:sz_maxx,sz_miny(start_pos):bx_sz:sz_maxy); % dictates the max and min bounds of the meshgrid.    
        for pt=1:length(xd) % x/y are of same original starting lenght (eg 730 elements, so this function iterates 730 times)
            z=sqrt(((mx-xd(pt)).^2)+((my-yd(pt)).^2));
            % every value of x is subtracted from every value of mx and the result squared.
            % every value of y is subtracted from every value of my and the result squared.
            % the pixel distance of xy from each and every mx/my point is calculated.
            [gp_x,gp_y]=find(z==min(min(z))); % the minimum distance, identifies the nearest mx/my to the xy point of interest.
            % this mx/my point identified by searching for min z, is labelled gridpoint x, gridpoint y.
            new_im(gp_x, gp_y)=1;
            % all gp_x/y become 1, all other mx/my points are ignored (0).
            % the numbers of 1s, is equivalent to the number of boxes, at that grid calibre, which contain a pixel data point.
        end
    boxcount(bx_sz,start_pos)=nnz(new_im);% nnz(new_im) returns the number of nonzero elements in matrix new_im.
    % that is the boxcount at every iteration of this gridset (bx_sz).
    end
end

bc=mean(boxcount,2);
eps=[1:maxbox]/maxsz;
epsilon=eps.'; % transposes row vector of box-size into a row vector.
Lnbc=log(bc(4:end)); % log(x) is an alias for the natural logarithm ln(x).
Lnepsilon=log(epsilon(4:end));
%Lnbc=log(bc(1:end)); % log(x) is an alias for the natural logarithm ln(x).
%Lnepsilon=log(epsilon(1:end));


frac_coeff=polyfit(Lnepsilon,Lnbc,1);
FD=-(frac_coeff(1));

toc

end



