

function im1 = myTextInsert(im,xy_pxl_pos,txt,fnt_size)

%% Create the text mask 

% Make an image the same size and put text in it 

hf = figure('color','white','units','normalized','position',[.1 .1 .8 .8]); 

image(ones(size(im))); 

set(gca,'units','pixels','position',[5 5 size(im,2)-1 size(im,1)-1],'visible','off');

% Text at arbitrary position 

text('units','pixels','position',xy_pxl_pos,'fontsize',fnt_size,'string',txt,'fontweight','bold'); 

% Capture the text image 

% Note that the size will have changed by about 1 pixel 

tim = getframe(gca); 

close(hf);

% Extract the cdata

tim2 = tim.cdata;

% Make a mask with the negative of the text 

tmask = tim2==0; 

% Place white text 

% Replace mask pixels with UINT8 max 

im1 = im;
im1(tmask) = uint8(255); 

image(im1);

axis off
