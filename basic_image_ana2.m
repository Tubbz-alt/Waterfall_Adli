function IM_ANA = basic_image_ana2(im_struct,use_bg,roi,header, idx)

%ntotal = im_struct.N_IMGS;
ntotal = size(idx,1);

if use_bg
    bg = load([header im_struct.background_dat{1}]);
    %bg = fliplr(bg.img);
end
    

nv = roi.bottom - roi.top + 1;
nh = roi.right - roi.left + 1;

x_axis = im_struct.RESOLUTION(1)*((1:nh) - nh/2);
y_axis = im_struct.RESOLUTION(1)*((1:nv) - nv/2);



IM_ANA.x_profs = zeros(nh,ntotal);
IM_ANA.y_profs = zeros(nv,ntotal);
IM_ANA.x_max = zeros(1,ntotal);
IM_ANA.y_max = zeros(1,ntotal);
IM_ANA.x_cent = zeros(1,ntotal);
IM_ANA.y_cent = zeros(1,ntotal);
IM_ANA.x_rms = zeros(1,ntotal);
IM_ANA.y_rms = zeros(1,ntotal);
IM_ANA.sum   = zeros(1,ntotal);

disp(['Analyzing ' num2str(ntotal) ' images.']);

nc = 0;
for i = idx',
  nc=nc+1;
    if(mod(nc,100) == 0)
      disp(['Done analysing image # ' num2str(nc)]);
    end% if
    image = imread([header im_struct.dat{i}]);
    
    if roi.rot
        image = rot90(image,roi.rot);
    end
    
    if use_bg
        image = image - bg.img;
    end
    
    % substract static background, by avg box in 4 corners
    n_box = 50; % pixels
    stc_back = mean ([mean2(image(1:n_box, 1:n_box))  mean2(image(end-n_box:end, end-n_box:end))    mean2(image(end-n_box:end, 1:n_box))   mean2(image(1:n_box, end-n_box:end))  ]);
    image = image - round(stc_back);
%    imagesc(image); colorbar;
%    pause;
    
    image = image(roi.top:roi.bottom,roi.left:roi.right);
    x_prof = mean(image);
    y_prof = mean(image,2);
    
    [mx,ix] = max(x_prof);
    [my,iy] = max(y_prof);
    
    xc = wm(x_axis,x_prof,1);
    xrms = wm(x_axis,x_prof,2);
    
    yc = wm(y_axis,y_prof,1);
    yrms = wm(y_axis,y_prof,2);
    
    IM_ANA.x_max(nc) = x_axis(ix);
    IM_ANA.y_max(nc) = y_axis(iy);
    
    IM_ANA.x_cent(nc) = xc;
    IM_ANA.x_rms(nc) = xrms;
    
    IM_ANA.y_cent(nc) = yc;
    IM_ANA.y_rms(nc) = yrms;
    
    IM_ANA.x_profs(:,nc) = x_prof;
    IM_ANA.y_profs(:,nc) = y_prof;
    
    IM_ANA.sum(nc) = sum(image(:));
end

IM_ANA.x_axis = x_axis;
IM_ANA.y_axis = y_axis;
IM_ANA.roi = roi;
