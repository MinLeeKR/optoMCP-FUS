classdef MinStack
    
    properties
        TIF;
        tags;
        out_pkfnd;
        out_cntrd;
        out_track;
        out_track_plt;
    end
    
    % instance methods
    methods
        % constructor
        function obj = MinStack(inputArg1)

            if nargin==1
                obj.TIF = inputArg1;
            else
                obj.TIF = uint16([]);
            end                        
        end
        

        
            
        
        
%% Size Function
        function outputArg = size(obj)
            %METHOD1 
            outputArg = size(obj.TIF);
        end
        
        
%% Save Stack Function
        
        
        function outputArg = save(obj,filename)
            
            obj.TIF = uint16(obj.TIF);
            fiji_descr = ['ImageJ=1.52p' newline ...
            'images=' num2str(size(obj.TIF,3)*...
                              size(obj.TIF,4)*...
                              size(obj.TIF,5)) newline... 
            'channels=' num2str(size(obj.TIF,3)) newline...
            'slices=' num2str(size(obj.TIF,4)) newline...
            'frames=' num2str(size(obj.TIF,5)) newline... 
            'hyperstack=true' newline...
            'mode=RGB' newline...  
            'loop=false' newline...  
            'min=0.0' newline...      
            'max=65535.0'];  % change this to 256 if you use an 8bit image
        
            
            filename = convertCharsToStrings(filename);
            t = Tiff(filename+".tif",'w')
            tagstruct.ImageLength = size(obj.TIF,1);
            tagstruct.ImageWidth = size(obj.TIF,2);
            tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
            tagstruct.BitsPerSample = 16;
            tagstruct.SamplesPerPixel = 1;
            tagstruct.Compression = Tiff.Compression.LZW;
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
            tagstruct.ImageDescription = fiji_descr;
            for frame = 1:size(obj.TIF,5)
                for slice = 1:size(obj.TIF,4)
                    for channel = 1:size(obj.TIF,3)
                        t.setTag(tagstruct)
                        t.write(im2uint16(obj.TIF(:,:,channel,slice,frame)));
                        t.writeDirectory(); % saves a new page in the tiff file
                    end
                end
            end
            t.close()
        end
        
        
%% Save data Function

        function outputArg = saveData(obj,filename)
        
        obj.TIF = [];
        save(filename,'obj');
            
        end

%% Bpass
    
        function obj2 = bpass(obj,channel,lnoise,lobject,threshold)

            % 
            % NAME:
            %               bpass
            % PURPOSE:
            %               Implements a real-space bandpass filter that suppresses 
            %               pixel noise and long-wavelength image variations while 
            %               retaining information of a characteristic size.
            % 
            % CATEGORY:
            %               Image Processing
            % CALLING SEQUENCE:
            %               res = bpass( image_array, lnoise, lobject )
            % INPUTS:
            %               image:  The two-dimensional array to be filtered.
            %               lnoise: Characteristic lengthscale of noise in pixels.
            %                       Additive noise averaged over this length should
            %                       vanish. May assume any positive floating value.
            %                       May be set to 0 or false, in which case only the
            %                       highpass "background subtraction" operation is 
            %                       performed.
            %               lobject: (optional) Integer length in pixels somewhat 
            %                       larger than a typical object. Can also be set to 
            %                       0 or false, in which case only the lowpass 
            %                       "blurring" operation defined by lnoise is done,
            %                       without the background subtraction defined by
            %                       lobject.  Defaults to false.
            %               threshold: (optional) By default, after the convolution,
            %                       any negative pixels are reset to 0.  Threshold
            %                       changes the threshhold for setting pixels to
            %                       0.  Positive values may be useful for removing
            %                       stray noise or small particles.  Alternatively, can
            %                       be set to -Inf so that no threshholding is
            %                       performed at all.
            %
            % OUTPUTS:
            %               res:    filtered image.
            % PROCEDURE:
            %               simple convolution yields spatial bandpass filtering.
            % NOTES:
            % Performs a bandpass by convolving with an appropriate kernel.  You can
            % think of this as a two part process.  First, a lowpassed image is
            % produced by convolving the original with a gaussian.  Next, a second
            % lowpassed image is produced by convolving the original with a boxcar
            % function. By subtracting the boxcar version from the gaussian version, we
            % are using the boxcar version to perform a highpass.
            % 
            % original - lowpassed version of original => highpassed version of the
            % original
            % 
            % Performing a lowpass and a highpass results in a bandpassed image.
            % 
            % Converts input to double.  Be advised that commands like 'image' display 
            % double precision arrays differently from UINT8 arrays.

            % MODIFICATION HISTORY:
            %               Written by David G. Grier, The University of Chicago, 2/93.
            %
            %               Greatly revised version DGG 5/95.
            %
            %               Added /field keyword JCC 12/95.
            % 
            %               Memory optimizations and fixed normalization, DGG 8/99.
            %               Converted to Matlab by D.Blair 4/2004-ish
            %
            %               Fixed some bugs with conv2 to make sure the edges are
            %               removed D.B. 6/05
            %
            %               Removed inadvertent image shift ERD 6/05
            % 
            %               Added threshold to output.  Now sets all pixels with
            %               negative values equal to zero.  Gets rid of ringing which
            %               was destroying sub-pixel accuracy, unless window size in
            %               cntrd was picked perfectly.  Now centrd gets sub-pixel
            %               accuracy much more robustly ERD 8/24/05
            %
            %               Refactored for clarity and converted all convolutions to
            %               use column vector kernels for speed.  Running on my 
            %               macbook, the old version took ~1.3 seconds to do
            %               bpass(image_array,1,19) on a 1024 x 1024 image; this
            %               version takes roughly half that. JWM 6/07
            %
            %       This code 'bpass.pro' is copyright 1997, John C. Crocker and 
            %       David G. Grier.  It should be considered 'freeware'- and may be
            %       distributed freely in its original form when properly attributed.  
            if nargin < 4, lobject = false; end
            if nargin < 5, threshold = 0; end
            obj2 = obj;
%             obj2.TIF = double(obj2.TIF);
            
            if max(size(lnoise))==1
                lnoise = lnoise*ones(size(channel));
            end
            
            if lobject
                if max(size(lobject))==1
                    lobject = lobject*ones(size(channel));
                end            
            end
            
            for cc = 1 : max(size(channel))
                for z = 1 : size(obj.TIF,4)
                    for t = 1 : size(obj.TIF,5)

                        image_array = double(obj.TIF(:,:,channel(cc),z,t));


                        normalize = @(x) x/sum(x);

                        image_array = double(image_array);

                        if lnoise(cc) == 0
                          gaussian_kernel = 1;
                        else      
                          gaussian_kernel = normalize(...
                            exp(-((-ceil(5*lnoise(cc)):ceil(5*lnoise(cc)))/(2*lnoise(cc))).^2));
                        end

                        if lobject(cc)  
                          boxcar_kernel = normalize(...
                              ones(1,length(-round(lobject(cc)):round(lobject(cc)))));
                        end

                        % JWM: Do a 2D convolution with the kernels in two steps each.  It is
                        % possible to do the convolution in only one step per kernel with 
                        %
                          % gconv = conv2(gaussian_kernel',gaussian_kernel,image_array,'same');
                          % bconv = conv2(boxcar_kernel', boxcar_kernel,image_array,'same');
                        % 
                        % but for some reason, this is slow.  The whole operation could be reduced
                        % to a single step using the associative and distributive properties of
                        % convolution:
                        %
                          % filtered = conv2(image_array,...
                          %   gaussian_kernel'*gaussian_kernel - boxcar_kernel'*boxcar_kernel,...
                          %   'same');
                        %
                        % But this is also comparatively slow (though inexplicably faster than the
                        % above).  It turns out that convolving with a column vector is faster than
                        % convolving with a row vector, so instead of transposing the kernel, the
                        % image is transposed twice.

                        gconv = conv2(image_array',gaussian_kernel','same');
                        gconv = conv2(gconv',gaussian_kernel','same');


                        if lobject(cc)
                          bconv = conv2(image_array',boxcar_kernel','same');
                          bconv = conv2(bconv',boxcar_kernel','same');

                          filtered = gconv - bconv;
                        else
                          filtered = gconv;
                        end

                        % Zero out the values on the edges to signal that they're not useful.     
                        lzero = max(lobject(cc),ceil(5*lnoise(cc)));

                        filtered(1:(round(lzero)),:) = 0;
                        filtered((end - lzero + 1):end,:) = 0;
                        filtered(:,1:(round(lzero))) = 0;
                        filtered(:,(end - lzero + 1):end) = 0;

                        % JWM: I question the value of zeroing out negative pixels.  It's a
                        % nonlinear operation which could potentially mess up our expectations
                        % about statistics.  Is there data on 'Now centroid gets subpixel accuracy
                        % much more robustly'?  To choose which approach to take, uncomment one of
                        % the following two lines.
                        % ERD: The negative values shift the peak if the center of the cntrd mask
                        % is not centered on the particle.

                        % res = filtered;
                        filtered(filtered < threshold) = 0;
                        obj2.TIF(:,:,channel(cc),z,t) = uint16(filtered);
    %                     obj2 = MinStack(filtered);
                    end
                end
            end
        end

    
%% PKFND

        function obj2=pkfnd(obj,channel,th,sz)
        % finds local maxima in an image to pixel level accuracy.   
        %  this provides a rough guess of particle
        %  centers to be used by cntrd.m.  Inspired by the lmx subroutine of Grier
        %  and Crocker's feature.pro
        % INPUTS:
        % im: image to process, particle should be bright spots on dark background with little noise
        %   ofen an bandpass filtered brightfield image (fbps.m, fflt.m or bpass.m) or a nice
        %   fluorescent image
        % th: the minimum brightness of a pixel that might be local maxima. 
        %   (NOTE: Make it big and the code runs faster
        %   but you might miss some particles.  Make it small and you'll get
        %   everything and it'll be slow.)
        % sz:  if your data's noisy, (e.g. a single particle has multiple local
        % maxima), then set this optional keyword to a value slightly larger than the diameter of your blob.  if
        % multiple peaks are found withing a radius of sz/2 then the code will keep
        % only the brightest.  Also gets rid of all peaks within sz of boundary
        %OUTPUT:  a N x 2 array containing, [row,column] coordinates of local maxima
        %           out(:,1) are the x-coordinates of the maxima
        %           out(:,2) are the y-coordinates of the maxima
        %CREATED: Eric R. Dufresne, Yale University, Feb 4 2005
        %MODIFIED: ERD, 5/2005, got rid of ind2rc.m to reduce overhead on tip by
        %  Dan Blair;  added sz keyword 
        % ERD, 6/2005: modified to work with one and zero peaks, removed automatic
        %  normalization of image
        % ERD, 6/2005: due to popular demand, altered output to give x and y
        %  instead of row and column
        % ERD, 8/24/2005: pkfnd now exits politely if there's nothing above
        %  threshold instead of crashing rudely
        % ERD, 6/14/2006: now exits politely if no maxima found
        % ERD, 10/5/2006:  fixed bug that threw away particles with maxima
        %  consisting of more than two adjacent points

        obj2 = obj;
        
        if max(size(th))==1
            th = th*ones(size(channel));
        end

        if sz
            if max(size(sz))==1
                sz = sz*ones(size(channel));
            end            
        end    
        
        
        for cc = 1 : max(size(channel))                
            for z = 1 : size(obj.TIF,4)
                for t = 1 : size(obj.TIF,5)

                    im = double(obj.TIF(:,:,channel(cc),z,t));

                    %find all the pixels above th(cc)reshold
                    %im=im./max(max(im)); 
                    ind=find(im > th(cc));
                    [nr,nc]=size(im);
                    tst=zeros(nr,nc);
                    n=length(ind);
                    if n==0
                        out=[];
                        obj2.out_pkfnd{channel(cc),z,t} = out;
                        display('nothing above threshold');
                        continue;
                    end
                    mx=[];
                    %convert index from find to row and column
                    rc=[mod(ind,nr),floor(ind/nr)+1];
                    for i=1:n
                        r=rc(i,1);c=rc(i,2);
                        %check each pixel above threshold to see if it's brighter than it's neighbors
                        %  THERE'S GOT TO BE A FASTER WAY OF DOING THIS.  I'M CHECKING SOME MULTIPLE TIMES,
                        %  BUT th(cc)IS DOESN'T SEEM THAT SLOW COMPARED TO THE OTHER ROUTINES, ANYWAY.
                        if r>1 & r<nr & c>1 & c<nc
                            if im(r,c)>=im(r-1,c-1) & im(r,c)>=im(r,c-1) & im(r,c)>=im(r+1,c-1) & ...
                             im(r,c)>=im(r-1,c)  & im(r,c)>=im(r+1,c) &   ...
                             im(r,c)>=im(r-1,c+1) & im(r,c)>=im(r,c+1) & im(r,c)>=im(r+1,c+1)
                            mx=[mx,[r,c]']; 
                            %tst(ind(i))=im(ind(i));
                            end
                        end
                    end
                    %out=tst;
                    mx=mx';

                    [npks,crap]=size(mx);

                    %if size is specified, th(cc)en get ride of pks with(cc)in size of boundary
                    if nargin==3 & npks>0
                       %th(cc)row out all pks with(cc)in sz(cc) of boundary;
                        ind=find(mx(:,1)>sz(cc) & mx(:,1)<(nr-sz(cc)) & mx(:,2)>sz(cc) & mx(:,2)<(nc-sz(cc)));
                        mx=mx(ind,:); 
                    end

                    %prevent from finding peaks with(cc)in size of each oth(cc)er
                    [npks,crap]=size(mx);
                    if npks > 1 
                        %CREATE AN IMAGE WIth(cc) ONLY PEAKS
                        nmx=npks;
                        tmp=0.*im;
                        for i=1:nmx
                            tmp(mx(i,1),mx(i,2))=im(mx(i,1),mx(i,2));
                        end
                        %LOOK IN NEIGHBORHOOD AROUND EACH PEAK, PICK th(cc)E BRIGHTEST
                        for i=1:nmx
                            roi=tmp( (mx(i,1)-floor(sz(cc)/2)):(mx(i,1)+(floor(sz(cc)/2)+1)),(mx(i,2)-floor(sz(cc)/2)):(mx(i,2)+(floor(sz(cc)/2)+1))) ;
                            [mv,indi]=max(roi);
                            [mv,indj]=max(mv);
                            tmp( (mx(i,1)-floor(sz(cc)/2)):(mx(i,1)+(floor(sz(cc)/2)+1)),(mx(i,2)-floor(sz(cc)/2)):(mx(i,2)+(floor(sz(cc)/2)+1)))=0;
                            tmp(mx(i,1)-floor(sz(cc)/2)+indi(indj)-1,mx(i,2)-floor(sz(cc)/2)+indj-1)=mv;
                        end
                        ind=find(tmp>0);
                        mx=[mod(ind,nr),floor(ind/nr)+1];
                    end

                    if size(mx)==[0,0]
                        out=[];
                    else
                        out(:,2)=mx(:,1);
                        out(:,1)=mx(:,2);
                    end
                    obj2.out_pkfnd{channel(cc),z,t} = out;
                    clear out;

                end
            end
        end
        end
        



%% CNTRD
        function obj2 = cntrd(obj,channel,lobject,interactive)
            
            % out=cntrd(im,threshold,lobject,interactive)
            % 
            % PURPOSE:  calculates the centroid of bright spots to sub-pixel accuracy.
            %  Inspired by Grier & Crocker's feature for IDL, but greatly simplified and optimized
            %  for matlab
            % 
            % INPUT:
            % im: image to process, particle should be bright spots on dark background with little noise
            %   ofen an bandpass filtered brightfield image or a nice fluorescent image
            %
            % obj.out_cntrd: locations of local maxima to pixel-level accuracy from pkfnd.m
            %
            % lobject: diamter of the window over which to average to calculate the centroid.  
            %     should be big enough
            %     to capture the whole particle but not so big that it captures others.  
            %     if initial guess of center (from pkfnd) is far from the centroid, the
            %     window will need to be larger than the particle size.  RECCOMMENDED
            %     size is the long lengthscale used in bpass plus 2.
            %     
            %
            % interactive:  OPTIONAL INPUT set this variable to one and it will show you the image used to calculate  
            %    each centroid, the pixel-level peak and the centroid
            %
            % NOTE:
            %  - if pkfnd, and cntrd return more then one location per particle then
            %  you should try to filter your input more carefully.  If you still get
            %  more than one peak for particle, use the optional lobject parameter in pkfnd
            %  - If you want sub-pixel accuracy, you need to have a lot of pixels in your window (lobject>>1). 
            %    To check for pixel bias, plot a histogram of the fractional parts of the resulting locations
            %  - It is HIGHLY recommended to run in interactive mode to adjust the parameters before you
            %    analyze a bunch of images.
            %
            % OUTPUT:  a N x 4 array containing, x, y and brightness for each feature
            %           out(:,1) is the x-coordinates
            %           out(:,2) is the y-coordinates
            %           out(:,3) is the brightnesses
            %           out(:,4) is the sqare of the radius of gyration
            %
            % CREATED: Eric R. Dufresne, Yale University, Feb 4 2005
            %  5/2005 inputs diamter instead of radius
            %  Modifications:
            %  D.B. (6/05) Added code from imdist/dist to make this stand alone.
            %  ERD (6/05) Increased frame of reject locations around edge to 1.5*lobject
            %  ERD 6/2005  By popular demand, 1. altered input to be formatted in x,y
            %  space instead of row, column space  2. added forth column of output,
            %  rg^2
            %  ERD 8/05  Outputs had been shifted by [0.5,0.5] pixels.  No more!
            %  ERD 8/24/05  Woops!  That last one was a red herring.  The real problem
            %  is the "ringing" from the output of bpass.  I fixed bpass (see note),
            %  and no longer need this kludge.  Also, made it quite nice if threshold=[];
            %  ERD 6/06  Added size and brightness output ot interactive mode.  Also 
            %   fixed bug in calculation of rg^2
            %  JWM 6/07  Small corrections to documentation 
            
            obj2 = obj;
            
            if nargin==3
               interactive=0; 
            end

            if lobject/2 == floor(lobject/2)
            warning('lobject must be odd, like bpass');
            end

            if isempty(obj.out_pkfnd)
                warning('there were no positions inputted into cntrd. check your pkfnd theshold')
                out=[];
                return;
            end
            
            if max(size(lobject))==1
                lobject = lobject*ones(size(channel));
            end



            for cc = 1 : max(size(channel))                
                for z = 1 : size(obj.TIF,4)
                    for t = 1 : size(obj.TIF,5)
                        
                        mx = obj.out_pkfnd{channel(cc),z,t};
                        if max(size(mx))==0
                            continue;
                        end
                        
                        im = double(obj.TIF(:,:,channel(cc),z,t));
                        r=(lobject(cc)+1)/2;
                        %create mask - window around trial location over which to calculate the centroid
                        m = 2*r;
                        x = 0:(m-1) ;
                        cent = (m-1)/2;
                        x2 = (x-cent).^2;
                        dst=zeros(m,m);
                        for i=1:m
                            dst(i,:)=sqrt((i-1-cent)^2+x2);
                        end


                        ind=find(dst < r);

                        msk=zeros([2*r,2*r]);
                        msk(ind)=1.0;
                        %msk=circshift(msk,[-r,-r]);

                        dst2=msk.*(dst.^2);
                        ndst2=sum(sum(dst2));

                        [nr,nc]=size(im);
                        %remove all potential locations within distance lobject(cc) from edges of image
                        ind=find(mx(:,2) > 1.5*lobject(cc) & mx(:,2) < nr-1.5*lobject(cc));
                        mx=mx(ind,:);
                        ind=find(mx(:,1) > 1.5*lobject(cc) & mx(:,1) < nc-1.5*lobject(cc));
                        mx=mx(ind,:);

                        [nmx,crap] = size(mx);

                        %inside of the window, assign an x and y coordinate for each pixel
                        xl=zeros(2*r,2*r);
                        for i=1:2*r
                            xl(i,:)=(1:2*r);
                        end
                        yl=xl';

                        pts=[];
                        %loop through all of the candidate positions
                        for i=1:nmx
                            %create a small working array around each candidate location, and apply the window function
                            tmp=msk.*im((mx(i,2)-r+1:mx(i,2)+r),(mx(i,1)-r+1:mx(i,1)+r));
                            %calculate the total brightness
                            norm=sum(sum(tmp));
                            %calculate the weigthed average x location
                            xavg=sum(sum(tmp.*xl))./norm;
                            %calculate the weighted average y location
                            yavg=sum(sum(tmp.*yl))./norm;
                            %calculate the radius of gyration^2
                            %rg=(sum(sum(tmp.*dst2))/ndst2);
                            rg=(sum(sum(tmp.*dst2))/norm);

                            %concatenate it up
                            pts=[pts,[mx(i,1)+xavg-r,mx(i,2)+yavg-r,norm,rg]'];

                            %OPTIONAL plot things up if you're in interactive mode
                            if interactive==1
                             imagesc(tmp)


                             hold on;
                             plot(xavg,yavg,'x')
                             plot(xavg,yavg,'o')
                             plot(r,r,'.')
                             hold off
                             title(['brightness ',num2str(norm),' size ',num2str(sqrt(rg))])
                             pause
                            end
                        end
                        obj2.out_cntrd{channel(cc),z,t}=pts';
                    end
                end
            end
            


        end
        
        
        
        
%% Track
        function obj2 = track(obj,channel,maxdisp,param)

        %;
        % ; see http://glinda.lrsm.upenn.edu/~weeks/idl
        % ;   for more information
        % ;
        % ;+
        % ; NAME:
        % ; track
        % ; PURPOSE:
        % ; Constructs n-dimensional trajectories from a scrambled list of
        % ; particle coordinates determined at discrete times (e.g. in
        % ; consecutive video frames).
        % ; CATEGORY:
        % ; Image Processing
        % ; CALLING SEQUENCE:
        % ; result = track( positionlist, maxdisp, param )
        % ;  set all keywords in the space below
        % ; INPUTS:
        % ; positionlist: an array listing the scrambled coordinates and data 
        % ;     of the different particles at different times, such that:
        % ;  positionlist(0:d-1,*): contains the d coordinates and
        % ;     data for all the particles, at the different times. must be positve
        % ;  positionlist(d,*): contains the time t that the position 
        % ;     was determined, must be integers (e.g. frame number.  These values must 
        % ;               be monotonically increasing and uniformly gridded in time.
        % ; maxdisp: an estimate of the maximum distance that a particle 
        % ;     would move in a single time interval.(see Restrictions)
        %  OPTIONAL INPUT:
        %   param:  a structure containing a few tracking parameters that are
        %       needed for many applications.  If param is not included in the
        %       function call, then default values are used.  If you set one value
        %       make sure you set them all:
        % ;         param.mem: this is the number of time steps that a particle can be
        % ;             'lost' and then recovered again.  If the particle reappears
        % ;             after this number of frames has elapsed, it will be
        % ;             tracked as a new particle. The default setting is zero.
        % ;             this is useful if particles occasionally 'drop out' of
        % ;             the data.
        % ;         param.dim: if the user would like to unscramble non-coordinate data
        % ;             for the particles (e.g. apparent radius of gyration for
        % ;             the particle images), then positionlist should
        % ;             contain the position data in positionlist(0:param.dim-1,*)
        % ;             and the extra data in positionlist(param.dim:d-1,*). It is then
        % ;             necessary to set dim equal to the dimensionality of the
        % ;             coordinate data to so that the track knows to ignore the
        % ;             non-coordinate data in the construction of the 
        % ;             trajectories. The default value is two.
        % ;         param.good: set this keyword to eliminate all trajectories with
        % ;             fewer than param.good valid positions.  This is useful
        % ;             for eliminating very short, mostly 'lost' trajectories
        % ;             due to blinking 'noise' particles in the data stream.
        %;          param.quiet: set this keyword to 1 if you don't want any text
        % ; OUTPUTS:
        % ; result:  a list containing the original data rows sorted 
        % ;     into a series of trajectories.  To the original input 
        % ;     data structure there is appended an additional column 
        % ;     containing a unique 'id number' for each identified 
        % ;     particle trajectory.  The result array is sorted so 
        % ;     rows with corresponding id numbers are in contiguous 
        % ;     blocks, with the time variable a monotonically
        % ;     increasing function inside each block.  For example:
        % ;     
        % ;     For the input data structure (positionlist):
        % ;         (x)      (y)      (t)
        % ;     pos = 3.60000      5.00000      0.00000
        % ;           15.1000      22.6000      0.00000
        % ;           4.10000      5.50000      1.00000 
        % ;           15.9000      20.7000      2.00000
        % ;           6.20000      4.30000      2.00000
        % ;
        % ;     IDL> res = track(pos,5,mem=2)
        % ;
        % ;     track will return the result 'res'
        % ;         (x)      (y)      (t)          (id)
        % ;     res = 3.60000      5.00000      0.00000      0.00000
        % ;           4.10000      5.50000      1.00000      0.00000
        % ;           6.20000      4.30000      2.00000      0.00000
        % ;           15.1000      22.6000      0.00000      1.00000
        % ;           15.9000      20.7000      2.00000      1.00000
        % ;
        % ;     NB: for t=1 in the example above, one particle temporarily
        % ;     vanished.  As a result, the trajectory id=1 has one time
        % ;     missing, i.e. particle loss can cause time gaps to occur 
        % ;     in the corresponding trajectory list. In contrast:
        % ;
        % ;     IDL> res = track(pos,5)
        % ;
        % ;     track will return the result 'res'
        % ;         (x)      (y)      (t)          (id)
        % ;     res = 15.1000      22.6000      0.00000      0.00000
        % ;                   3.60000      5.00000      0.00000      1.00000
        % ;               4.10000      5.50000      1.00000      1.00000
        % ;               6.20000      4.30000      2.00000      1.00000
        % ;               15.9000      20.7000      2.00000      2.00000
        % ; 
        % ;     where the reappeared 'particle' will be labelled as new
        % ;     rather than as a continuation of an old particle since
        % ;     mem=0.  It is up to the user to decide what setting of 
        % ;     'mem' will yeild the highest fidelity .
        % ; 
        % ; SIDE EFFECTS:
        % ; Produces informational messages.  Can be memory intensive for
        % ; extremely large data sets.
        % ; RESTRICTIONS:
        % ; maxdisp should be set to a value somewhat less than the mean 
        % ; spacing between the particles. As maxdisp approaches the mean
        % ; spacing the runtime will increase significantly. The function 
        % ; will produce an error message: "Excessive Combinatorics!" if
        % ; the run time would be too long, and the user should respond 
        % ; by re-executing the function with a smaller value of maxdisp.
        % ; Obviously, if the particles being tracked are frequently moving
        % ; as much as their mean separation in a single time step, this
        % ; function will not return acceptable trajectories.
        % ; PROCEDURE:
        % ; Given the positions for n particles at time t(i), and m possible
        % ; new positions at time t(i+1), this function considers all possible 
        % ; identifications of the n old positions with the m new positions,
        % ; and chooses that identification which results in the minimal total
        % ; squared displacement. Those identifications which don't associate
        % ; a new position within maxdisp of an old position ( particle loss )
        % ; penalize the total squared displacement by maxdisp^2. For non-
        % ; interacting Brownian particles with the same diffusivity, this
        % ; algorithm will produce the most probable set of identifications 
        % ; ( provided maxdisp >> RMS displacement between frames ).
        % ; In practice it works reasonably well for systems with oscillatory,
        % ; ballistic, correlated and random hopping motion, so long as single 
        % ; time step displacements are reasonably small.  NB: multidimensional
        % ; functionality is intended to facilitate tracking when additional
        % ; information regarding target identity is available (e.g. size or 
        % ; color).  At present, this information should be rescaled by the
        % ; user to have a comparable or smaller (measurement) variance than 
        % ; the spatial displacements.
        % ;
        % ; MODIFICATION HISTORY:
        % ;  2/93 Written by John C. Crocker, University of Chicago (JFI).
        % ;  7/93 JCC fixed bug causing particle loss and improved performance
        % ;     for large numbers of (>100) particles.
        % ; 11/93 JCC improved speed and memory performance for large
        % ;     numbers of (>1000) particles (added subnetwork code).
        % ;  3/94 JCC optimized run time for trivial bonds and d<7. (Added
        % ;     d-dimensional raster metric code.)
        % ;  8/94 JCC added functionality to unscramble non-position data
        % ;     along with position data.
        % ;  9/94 JCC rewrote subnetwork code and wrote new, more efficient 
        % ;     permutation code.
        % ;  5/95 JCC debugged subnetwork and excessive combinatorics code.
        % ; 12/95 JCC added memory keyword, and enabled the tracking of
        % ;     newly appeared particles.
        % ;  3/96 JCC made inipos a keyword, and disabled the adding of 'new'
        % ;     particles when inipos was set.
        % ;  3/97 JCC added 'add' keyword, since Chicago users didn't like 
        % ;     having particle addition be the default. 
        % ;  9/97 JCC added 'goodenough' keyword to improve memory efficiency
        % ;     when using the 'add' keyword and to filter out bad tracks.
        % ;       10/97 JCC streamlined data structure to speed runtime for >200 
        % ;               timesteps.  Changed 'quiet' keyword to 'verbose'. Made
        % ;               time labelling more flexible (uniform and sorted is ok).
        % ;  9/98 JCC switched trajectory data structure to a 'list' form,
        % ;     resolving memory issue for large, noisy datasets.
        % ;  2/99 JCC added Eric Weeks's 'uberize' code to post-facto 
        % ;     rationalize the particle id numbers, removed 'add' keyword.
        % ;  1/05 Transmuted to MATLAB by D. Blair
        % ;  5/05  ERD Added the param structure to simplify calling.
        %    6/05  ERD Added quiet to param structure
        %    7/05  DLB Fixed slight bug in trivial bond code
        %    3/07  DLB Fixed bug with max disp pointed out by Helene Delanoe-Ayari
        %
        % ; This code 'track.pro' is copyright 1999, by John C. Crocker. 
        % ; It should be considered 'freeware'- and may be distributed freely 
        % ; (outside of the military-industrial complex) in its original form 
        % ; when properly attributed.
        % ;
        % ;-

        obj2 = obj;

        for cc = 1 : max(size(channel))
            xyzs = [];
            for t = 1 : size(obj.TIF,5)
                if max(size(obj.out_cntrd{channel(cc),1,t}))==0
                    continue;
                end
            xyzs = [xyzs;[obj.out_cntrd{channel(cc),1,t}(:,[1 2]), t*ones(size(obj.out_cntrd{channel(cc),1,t},1),1)]];
            end
            
            dd = length(xyzs(1,:));

            %use default parameters if none given
            if nargin==3
                %default values
                memory_b=0; % if mem is not needed set to zero
                goodenough = 0;  % if goodenough is not wanted set to zero
                dim = dd - 1;
                quiet=0;
            else
                memory_b    =   param.mem;
                goodenough  =   param.good;
                dim         =   param.dim;
                quiet       =   param.quiet;
            end


            % checking the input time vector
            t = xyzs(:,dd);
            st = circshift(t,1);
            st = t(2:end) - st(2:end);
            if  sum(st(find(st < 0))) ~= 0
                disp('The time vectors is not in order')
                return
            end
            info = 1;

            w = find(st > 0);
            z = length(w);
            z = z +1;
            if isempty(w)
                disp('All positions are at the same time... go back!')
                return
            end

            % partitioning the data with unique times

            %res = unq(t);
            % implanting unq directly
                indices = find(t ~= circshift(t,-1));
                    count = length(indices);
                    if count > 0
                        res = indices;
                    else  
                        res = length(t) -1;
                    end
             %%%%%%%%%%%%%%%%%%%%%%%       

            res = [1,res',length(t)];
            ngood = res(2) - res(1) + 1;
            eyes = 1:ngood;
            pos = xyzs(eyes,1:dim);
            istart = 2;
            n = ngood;

            zspan = 50;
            if n > 200 
                zspan = 20;
            end
            if n > 500 
                zspan = 10;
            end
            resx = zeros(zspan,n) - 1;

            bigresx = zeros(z,n) - 1;
            mem = zeros(n,1);
            %  whos resx
            %  whos bigresx
            uniqid = 1:n;
            maxid = n;
            olist = [0.,0.];

            if goodenough > 0 
                dumphash = zeros(n,1);
                nvalid = ones(n,1);
            end

            %  whos eyes;
            resx(1,:) = eyes;
            % setting up constants
            maxdisq = maxdisp^2;

            % John calls this the setup for "fancy code" ???
            notnsqrd = (sqrt(n*ngood) > 200) & (dim < 7);
            notnsqrd = notnsqrd(1);

            if notnsqrd
                %;   construct the vertices of a 3x3x3... d-dimensional hypercube

                cube = zeros(3^dim,dim);


                for d=0:dim-1,
                    numb = 0;
                    for j=0:(3^d):(3^dim)-1,
                        cube(j+1:j+(3^(d)),d+1) = numb;
                        numb = mod(numb+1,3);
                    end
                end    

                %   calculate a blocksize which may be greater than maxdisp, but which
                %   keeps nblocks reasonably small.  

                volume = 1;
                for d = 0:dim-1
                    minn = min(xyzs(w,d+1));
                    maxx = max(xyzs(w,d+1));
                    volume = volume * (maxx-minn);
                end
                volume;
                blocksize = max( [maxdisp,((volume)/(20*ngood))^(1.0/dim)] );
            end
            %   Start the main loop over the frames.
            for i=istart:z
                ispan = mod(i-1,zspan)+1;
                %disp(ispan)
                % get new particle positions
                m = res(i+1) - res(i);
                res(i);
                eyes = 1:m;
                eyes = eyes + res(i);

                if m > 0

                    xyi = xyzs(eyes,1:dim);
                    found = zeros(m,1);

                    % THE TRIVIAL BOND CODE BEGINS   

                    if notnsqrd
                        %Use the raster metric code to do trivial bonds

                        % construct "s", a one dimensional parameterization of the space 
                        % which consists of the d-dimensional raster scan of the volume.)

                        abi = fix(xyi./blocksize);
                        abpos = fix(pos./blocksize);
                        si = zeros(m,1);
                        spos = zeros(n,1);
                        dimm = zeros(dim,1);
                        coff = 1.;

                        for j=1:dim
                            minn = min([abi(:,j);abpos(:,j)]);
                            maxx = max([abi(:,j);abpos(:,j)]);
                            abi(:,j) = abi(:,j) - minn;
                            abpos(:,j) = abpos(:,j) - minn;
                            dimm(j) = maxx-minn + 1;
                            si = si + abi(:,j).*coff;
                            spos = spos + abpos(:,j).*coff;
                            coff = dimm(j).*coff;
                        end
                        nblocks = coff;
                        % trim down (intersect) the hypercube if its too big to fit in the
                        % particle volume. (i.e. if dimm(j) lt 3)

                        cub = cube;
                        deg = find( dimm < 3);
                        if ~isempty(deg)
                            for j = 0:length(deg)-1
                                cub = cub(find(cub(:,deg(j+1)) < dimm(deg(j+1))),:);
                            end
                        end 

                        % calculate the "s" coordinates of hypercube (with a corner @ the origin)
                        scube = zeros(length(cub(:,1)),1);
                        coff = 1;
                        for j=1:dim
                            scube = scube + cub(:,j).*coff;
                            coff = coff*dimm(j);      
                        end

                        % shift the hypercube "s" coordinates to be centered around the origin

                        coff = 1;
                        for j=1:dim
                            if dimm(j) > 3
                                scube = scube - coff;
                            end
                            coff = dimm(j).* coff;
                        end
                        scube = mod((scube + nblocks),nblocks);
                        % get the sorting for the particles by their "s" positions.
                        [ed,isort] = sort(si);

                        % make a hash table which will allow us to know which new particles
                        % are at a given si.
                        strt = zeros(nblocks,1) -1;
                        fnsh = zeros(nblocks,1);
                        h = find(si == 0);
                        lh = length(h);
                        if lh > 0

                        si(h) = 1;  
                        end

                        for j=1:m
                            if strt(si(isort(j))) == -1
                                strt(si(isort(j))) = j;
                                fnsh(si(isort(j))) = j;
                            else
                                fnsh(si(isort(j))) = j;
                            end
                        end
                        if lh > 0
                        si(h) = 0;   
                        end
                        coltot = zeros(m,1);
                        rowtot = zeros(n,1);
                        which1 = zeros(n,1);
                        for j=1:n


                            map = fix(-1);

                            scub_spos = scube + spos(j);
                            s = mod(scub_spos,nblocks);
                            whzero = find(s == 0 );
                            if ~isempty(whzero)
                                nfk = find(s ~=0);
                                s = s(nfk);
                            end

                            w = find(strt(s) ~= -1);

                            ngood = length(w);
                            ltmax=0;
                            if ngood ~= 0

                                s = s(w);
                                for k=1:ngood
                                    map = [map;isort( strt(s(k)):fnsh(s(k)))];
                                end
                                map = map(2:end);
            %                     if length(map) == 2
            %                         if (map(1) - map(2)) == 0
            %                             map = unique(map);
            %                          end
            %                     end
                                %   map = map(umap);
                                %end
                                % find those trival bonds
                                distq = zeros(length(map),1);
                                for d=1:dim     
                                    distq = distq + (xyi(map,d) - pos(j,d)).^2;
                                end
                                ltmax = distq < maxdisq;

                                rowtot(j) = sum(ltmax);

                                if rowtot(j) >= 1 
                                    w = find(ltmax == 1);
                                    coltot( map(w) ) = coltot( map(w)) +1;
                                    which1(j) = map( w(1) );
                                end
                            end

                        end


                        ntrk = fix(n - sum(rowtot == 0));

                        w = find( rowtot == 1);
                        ngood = length(w);


                        if ngood ~= 0 
                            ww = find(coltot( which1(w) ) == 1);
                            ngood = length(ww);
                            if ngood ~= 0 
                                 %disp(size(w(ww)))
                                resx(ispan,w(ww)) = eyes( which1(w(ww)));
                                found( which1( w(ww))) = 1;
                                rowtot( w(ww)) = 0;
                                coltot( which1(w(ww))) = 0;
                            end
                        end

                        labely = find( rowtot > 0);
                        ngood = length(labely);
                        if ngood ~= 0 
                            labelx = find( coltot > 0);

                            nontrivial = 1;
                        else
                            nontrivial = 0;
                        end

                    else 

                        %   or: Use simple N^2 time routine to calculate trivial bonds      

                        % let's try a nice, loopless way!
                        % don't bother tracking perm. lost guys.
                        wh = find( pos(:,1) >= 0);
                        ntrack = length(wh);
                        if ntrack == 0 
                            'There are no valid particles to track idiot!'
                            break
                        end
                        xmat = zeros(ntrack,m);
                        count = 0;
                        for kk=1:ntrack
                            for ll=1:m
                                xmat(kk,ll) = count;
                                count = count+1;
                            end
                        end
                        count = 0;
                        for kk=1:m
                            for ll=1:ntrack
                                ymat(kk,ll) = count;
                                count = count+1;
                            end
                        end

                        xmat = (mod(xmat,m) + 1);
                        ymat = (mod(ymat,ntrack) +1)';
                        [lenxn,lenxm] = size(xmat);
            %            whos ymat
            %            whos xmat
            %            disp(m)

                        for d=1:dim
                            x = xyi(:,d);
                            y = pos(wh,d);
                            xm = x(xmat);
                            ym = y(ymat(1:lenxn,1:lenxm));
                            if size(xm) ~= size(ym)
                                xm = xm';
                            end

                            if d == 1
                                dq = (xm -ym).^2;
                                %dq = (x(xmat)-y(ymat(1:lenxn,1:lenxm))).^2;
                            else
                                dq = dq + (xm-ym).^2;
                                %dq = dq + (x(xmat)-y(ymat(1:lenxn,1:lenxm)) ).^2;
                            end
                        end

                        ltmax = dq < maxdisq;

                        % figure out which trivial bonds go with which

                        rowtot = zeros(n,1);
                        rowtot(wh) = sum(ltmax,2);


                        if ntrack > 1 
                            coltot = sum(ltmax,1);
                        else
                            coltot = ltmax;
                        end
                        which1 = zeros(n,1);
                        for j=1:ntrack 
                            [mx, w] = max(ltmax(j,:));
                            which1(wh(j)) = w;
                        end

                        ntrk = fix( n - sum(rowtot == 0));
                        w= find( rowtot == 1) ;
                        ngood = length(w);
                        if ngood ~= 0
                            ww = find(coltot(which1(w)) == 1);
                            ngood = length(ww);
                            if ngood ~= 0 
                                resx( ispan, w(ww) ) = eyes( which1( w(ww)));
                                found(which1( w(ww))) = 1;
                                rowtot(w(ww)) = 0;
                                coltot(which1(w(ww))) = 0;
                            end
                        end

                        labely = find( rowtot > 0);
                        ngood = length(labely);

                        if ngood ~= 0
                            labelx = find( coltot > 0);
                            nontrivial = 1;
                        else
                            nontrivial = 0;
                        end
                    end

                    %THE TRIVIAL BOND CODE ENDS

                    if nontrivial

                        xdim = length(labelx);
                        ydim = length(labely);

                        %  make a list of the non-trivial bonds            

                        bonds = zeros(1,2);
                        bondlen = 0;

                        for j=1:ydim
                            distq = zeros(xdim,1);

                            for d=1:dim
                                %distq
                                distq = distq + (xyi(labelx,d) - pos(labely(j),d)).^2;
                                %distq    
                            end

                            w= find(distq <  maxdisq)' - 1;
                            ngood = length(w);
                            newb = [w;(zeros(1,ngood)+j)];


                            bonds = [bonds;newb'];

                            bondlen = [ bondlen;distq( w + 1) ];

                        end
                        bonds = bonds(2:end,:);

                        bondlen = bondlen(2:end);
                        numbonds = length(bonds(:,1));
                        mbonds = bonds;
                        max([xdim,ydim]);


                        if max([xdim,ydim]) < 4
                            nclust = 1;
                            maxsz = 0;
                            mxsz = xdim;
                            mysz = ydim;
                            bmap = zeros(length(bonds(:,1)+1),1) - 1;

                        else


                            %   THE SUBNETWORK CODE BEGINS            
                            lista = zeros(numbonds,1);
                            listb = zeros(numbonds,1);
                            nclust = 0;
                            maxsz = 0;
                            thru = xdim;

                            while thru ~= 0
                                %  the following code extracts connected 
                                %   sub-networks of the non-trivial 
                                %   bonds.  NB: lista/b can have redundant entries due to 
                                %   multiple-connected subnetworks      


                                w = find(bonds(:,2) >= 0);
               %                 size(w)

                                lista(1) = bonds(w(1),2);
                                listb(1) = bonds(w(1),1);
                                bonds(w(1),:) = -(nclust+1);
                                bonds;
                                adda = 1;
                                addb = 1;
                                donea = 0;
                                doneb = 0;
                                if (donea ~= adda) | (doneb ~= addb)
                                    true = 0;
                                else
                                true = 1;   
                                end

                                while ~true

                                    if (donea ~= adda)
                                        w = find(bonds(:,2) == lista(donea+1));
                                        ngood = length(w);
                                        if ngood ~= 0 
                                            listb(addb+1:addb+ngood,1) = bonds(w,1);
                                            bonds(w,:) = -(nclust+1);
                                            addb = addb+ngood;
                                        end
                                        donea = donea+1;
                                    end
                                    if (doneb ~= addb) 
                                        w = find(bonds(:,1) == listb(doneb+1));
                                        ngood = length(w);
                                        if ngood ~= 0
                                            lista(adda+1:adda+ngood,1) = bonds(w,2);
                                            bonds(w,:) = -(nclust+1);
                                            adda = adda+ngood;
                                        end
                                        doneb = doneb+1;
                                    end
                                  if (donea ~= adda) | (doneb ~= addb) 
                                      true = 0;
                                  else  
                                      true = 1;
                                  end
                                end

                                [pp,pqx] = sort(listb(1:doneb));
                                %unx =  unq(listb(1:doneb),pqx);
                                %implanting unq directly
                                    arr = listb(1:doneb);
                                    q = arr(pqx);
                                    indices = find(q ~= circshift(q,-1));
                                    count = length(indices);
                                    if count > 0
                                        unx = pqx(indices);
                                    else
                                        unx = length(q) -1;
                                    end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                xsz = length(unx);
                                [pp,pqy] = sort(lista(1:donea));
                                %uny =  unq(lista(1:donea),pqy);
                                %implanting unq directly
                                    arr = lista(1:donea);
                                    q = arr(pqy);
                                    indices = find(q ~= circshift(q,-1));
                                    count = length(indices);
                                    if count > 0
                                        uny = pqy(indices);
                                    else
                                        uny = length(q) -1;
                                    end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   




                                ysz = length(uny);
                                if xsz*ysz > maxsz
                                    maxsz = xsz*ysz;
                                    mxsz = xsz;
                                    mysz = ysz; 
                                end


                                thru = thru -xsz;
                                nclust = nclust + 1;
                            end
                            bmap = bonds(:,2);                    
                        end
                        % THE SUBNETWORK CODE ENDS
                        % put verbose in for Jaci

                        %   THE PERMUTATION CODE BEGINS

                        for nc =1:nclust
                            w = find( bmap == -1*(nc));

                            nbonds = length(w);
                            bonds = mbonds(w,:);
                            lensq = bondlen(w);
                            [pq,st] = sort( bonds(:,1));
                            %un = unq(bonds(:,1),st);
                               %implanting unq directly     
                                    arr = bonds(:,1);
                                    q = arr(st);
                                    indices = find(q ~= circshift(q,-1));
                                    count = length(indices);
                                    if count > 0
                                        un = st(indices);
                                    else
                                        un = length(q) -1;
                                    end
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                            uold = bonds(un,1);

                            nold = length(uold);

                            %un = unq(bonds(:,2));

                            %implanting unq directly  
                            indices = find(bonds(:,2) ~= circshift(bonds(:,2),-1));
                            count = length(indices);
                                if count > 0
                                    un = indices;
                                else  
                                    un = length(bonds(:,2)) -1;
                                end
                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

                            unew = bonds(un,2);
                            nnew = length(unew);

                            if nnew > 5
                                rnsteps = 1;
                                for ii =1:nnew
                                    rnsteps = rnsteps * length( find(bonds(:,2) == ...
                                        unew(ii)));
                                    if rnsteps > 5.e+4
                                        disp('Warning: difficult combinatorics encountered.')
                                    end
                                    if rnsteps > 2.e+5
                                        disp(['Excessive Combinitorics you FOOL LOOK WHAT YOU HAVE' ...
                                                ' DONE TO ME!!!'])
                                        return
                                    end
                                end
                            end
                            st = zeros(nnew,1);
                            fi = zeros(nnew,1);
                            h = zeros(nbonds,1);
                            ok = ones(nold,1);
                            nlost = (nnew - nold) > 0;


                            for ii=1:nold 
                                h(find(bonds(:,1) == uold(ii))) = ii;
                            end
                            st(1) = 1 ;
                            fi(nnew) = nbonds; % check this later
                            if nnew > 1 
                                sb = bonds(:,2);
                                sbr = circshift(sb,1);
                                sbl = circshift(sb,-1);
                                st(2:end) = find( sb(2:end) ~= sbr(2:end)) + 1;
                                fi(1:nnew-1) = find( sb(1:nbonds-1) ~= sbl(1:nbonds-1));
                            end
            %                if i-1 == 13
            %                    hi
            %                end
                            checkflag = 0;
                            while checkflag ~= 2

                                pt = st -1;
                                lost = zeros(nnew,1);
                                who = 0;
                                losttot = 0;
                                mndisq = nnew*maxdisq;


                                while who ~= -1

                                    if pt(who+1) ~= fi(who+1)


                                        w = find( ok( h( pt( who+1 )+1:fi( who+1 ) ) ) ); % check this -1
                                        ngood = length(w);
                                        if ngood > 0
                                            if pt(who+1) ~= st(who+1)-1
                                                ok(h(pt(who+1))) = 1;
                                            end
                                            pt(who+1) = pt(who+1) + w(1);
                                            ok(h(pt(who+1))) = 0;
                                            if who == nnew -1
                                                ww = find( lost == 0);
                                                dsq = sum(lensq(pt(ww))) + losttot*maxdisq;

                                                if dsq < mndisq 
                                                    minbonds = pt(ww);
                                                    mndisq = dsq;
                                                end
                                            else
                                                who = who+1;
                                            end
                                        else
                                            if ~lost(who+1) & (losttot ~= nlost)
                                                lost(who+1) = 1;
                                                losttot = losttot + 1;
                                                if pt(who+1) ~= st(who+1) -1;
                                                    ok(h(pt(who+1))) = 1;
                                                end
                                                if who == nnew-1
                                                    ww = find( lost == 0);
                                                    dsq = sum(lensq(pt(ww))) + losttot*maxdisq;
                                                    if dsq < mndisq
                                                        minbonds = pt(ww);
                                                        mndisq = dsq;
                                                    end
                                                else    
                                                   who = who + 1;
                                                end

                                            else
                                                if pt(who+1) ~= (st(who+1) -1) 
                                                    ok(h(pt(who+1))) = 1;
                                                end
                                                pt(who+1) = st(who+1) -1;
                                                if lost(who+1) 
                                                    lost(who+1) = 0;
                                                    losttot = losttot -1;
                                                end
                                                who = who -1;
                                            end
                                        end
                                    else  
                                        if ~lost(who+1) & (losttot ~= nlost)
                                            lost(who+1) = 1;
                                            losttot = losttot + 1;
                                            if pt(who+1) ~= st(who+1)-1
                                                ok(h(pt(who+1))) = 1;
                                            end
                                            if who == nnew -1
                                                ww = find( lost == 0);
                                                dsq = sum(lensq(pt(ww))) + losttot*maxdisq;

                                                if dsq < mndisq
                                                    minbonds = pt(ww);
                                                    mndisq = dsq;
                                                end
                                            else
                                                who = who + 1;
                                            end
                                        else
                                            if pt(who+1) ~= st(who+1) -1
                                                ok(h(pt(who+1))) = 1;
                                            end
                                            pt(who+1) = st(who+1) -1;
                                            if lost(who+1) 
                                                lost(who+1) = 0;
                                                losttot = losttot -1;
                                            end
                                            who = who -1;
                                        end
                                    end
                                end

                                checkflag = checkflag + 1;
                                if checkflag == 1
                                    plost = min([fix(mndisq/maxdisq) , (nnew -1)]);
                                    if plost > nlost 
                                        nlost = plost; 
                                    else
                                        checkflag = 2;
                                    end
                                end

                            end  
                            %   update resx using the minimum bond configuration               

                            resx(ispan,labely(bonds(minbonds,2))) = eyes(labelx(bonds(minbonds,1)+1));
                            found(labelx(bonds(minbonds,1)+1)) = 1;

                        end

                        %   THE PERMUTATION CODE ENDS
                    end

                    w = find(resx(ispan,:) >= 0);
                    nww = length(w);

                    if nww > 0 
                        pos(w,:) = xyzs( resx(ispan,w) , 1:dim);
                        if goodenough > 0 
                            nvalid(w) = nvalid(w) + 1;
                        end
                    end  %go back and add goodenough keyword thing   
                    newguys = find(found == 0);

                   nnew = length(newguys);

                    if (nnew > 0) % & another keyword to workout inipos
                        newarr = zeros(zspan,nnew) -1;
                        resx = [resx,newarr];

                        resx(ispan,n+1:end) = eyes(newguys);
                        pos = [[pos];[xyzs(eyes(newguys),1:dim)]];
                        nmem = zeros(nnew,1);
                        mem = [mem;nmem];
                        nun = 1:nnew;
                        uniqid = [uniqid,((nun) + maxid)];
                        maxid = maxid + nnew;
                        if goodenough > 0 
                            dumphash = [dumphash;zeros(1,nnew)'];
                            nvalid = [nvalid;zeros(1,nnew)'+1];
                        end
                        % put in goodenough 
                        n = n + nnew;

                    end

                else
                    ' Warning- No positions found for t='
                end
                w = find( resx(ispan,:) ~= -1);
                nok = length(w);
                if nok ~= 0
                    mem(w) =0;
                end

                mem = mem + (resx(ispan,:)' == -1);
                wlost = find(mem == memory_b+1);
                nlost =length(wlost);

                if nlost > 0 
                    pos(wlost,:) = -maxdisp;
                    if goodenough > 0
                        wdump = find(nvalid(wlost) < goodenough);
                        ndump = length(wdump);
                        if ndump > 0
                            dumphash(wlost(wdump)) = 1;
                        end
                    end
                    % put in goodenough keyword stuff if 
                end
                if (ispan == zspan) | (i == z)
                    nold = length(bigresx(1,:));
                    nnew = n-nold;
                    if nnew > 0
                        newarr = zeros(z,nnew) -1;
                        bigresx = [bigresx,newarr];
                    end
                    if goodenough > 0  
                        if (sum(dumphash)) > 0
                            wkeep = find(dumphash == 0);
                            nkeep = length(wkeep);
                            resx = resx(:,wkeep);
                            bigresx = bigresx(:,wkeep);
                            pos = pos(wkeep,:);
                            mem = mem(wkeep);
                            uniqid = uniqid(wkeep);
                            nvalid = nvalid(wkeep);
                            n = nkeep;
                            dumphash = zeros(nkeep,1);
                        end
                    end

                    % again goodenough keyword
                    if quiet~=1
                        disp(strcat(num2str(i), ' of ' ,num2str(z), ' done.  Tracking  ',num2str(ntrk),' particles  ', num2str(n),' tracks total'));
                    end
                    bigresx(i-(ispan)+1:i,:) = resx(1:ispan,:);
                    resx = zeros(zspan,n) - 1;


                    wpull = find(pos(:,1) == -maxdisp);
                    npull = length(wpull);

                    if npull > 0
                        lillist = zeros(1,2);
                        for ipull=1:npull
                            wpull2 = find(bigresx(:,wpull(ipull)) ~= -1);
                            npull2 = length(wpull2);
                            thing = [bigresx(wpull2,wpull(ipull)),zeros(npull2,1)+uniqid(wpull(ipull))];
                            lillist = [lillist;thing];

                        end
                        olist = [[olist];[lillist(2:end,:)]];

                    end



                    wkeep = find(pos(:,1) >= 0);
                    nkeep = length(wkeep);
                    if nkeep == 0 
                            'Were going to crash now, no particles....'
                    end
                    resx = resx(:,wkeep);
                    bigresx = bigresx(:,wkeep);
                    pos = pos(wkeep,:);
                    mem = mem(wkeep);
                    uniqid = uniqid(wkeep);
                    n = nkeep;
                    dumphash = zeros(nkeep,1);
                    if goodenough > 0
                        nvalid = nvalid(wkeep);
                    end
                end

            end

            if goodenough > 0 
                nvalid = sum(bigresx >= 0 ,1);
                wkeep = find(nvalid >= goodenough);
                nkeep = length(wkeep);
                if nkeep == 0
                    for i=1:10
                    disp('You are not going any further, check your params and data')
                    end
                    disp('the code broke at line 1045')
                    return
                end
                if nkeep < n
                    bigresx = bigresx(:,wkeep);
                    n = nkeep;
                    uniqid = uniqid(wkeep);
                    pos = pos(wkeep,:);
                end
            end


            wpull = find( pos(:,1) ~= -2*maxdisp);
            npull = length(wpull);
            if npull > 0
                lillist = zeros(1,2);
                for ipull=1:npull
                    wpull2 = find(bigresx(:,wpull(ipull)) ~= -1);
                    npull2 = length(wpull2);   
                    thing = [bigresx(wpull2,wpull(ipull)),zeros(npull2,1)+uniqid(wpull(ipull))];
                    lillist = [lillist;thing];
                end
                olist = [olist;lillist(2:end,:)];
            end

            olist = olist(2:end,:);
            %bigresx = 0;
            %resx = 0;

            nolist = length(olist(:,1));
            res = zeros(nolist,dd+1);
            for j=1:dd
                res(:,j) = xyzs(olist(:,1),j);
            end
            res(:,dd+1) = olist(:,2);

            % this is uberize included for simplicity of a single monolithic code

            ndat=length(res(1,:));
            newtracks=res;


            %u=unq(newtracks(:,ndat));

            % inserting unq
            indices = find(newtracks(:,ndat) ~= circshift(newtracks(:,ndat),-1));
                    count = length(indices);
                    if count > 0
                        u = indices;
                    else  
                        u = length(newtracks(:,ndat)) -1;
                    end


            ntracks=length(u);
            u=[0;u];
            for i=2:ntracks+1
                newtracks(u(i-1)+1:u(i),ndat) = i-1;
            end

            obj2.out_track{channel(cc)} = newtracks;
            
            for t = 1 : size(obj.TIF,5)
                obj2.out_track_plt{channel(cc),t} = newtracks(newtracks(:,3)==t,:);
            end
            


        end
        end



%% Transform

        function obj2 = RtoGtform(obj,FileName)
            
            obj2 = MinStack;
            obj2.TIF = im2uint16(obj.TIF);
            load(FileName);
            
            %TIFFStack : xytzc --> saveStacks : xyczt
            for t=1:size(obj.TIF,5)
                tempImage_mCh = obj.TIF(:,:,1,:,t);
                for z = 1 : size(obj.TIF,4)
                    image(:,:,z) = tempImage_mCh(:,:,1,z);
                    transImageBuffrd = uint16(zeros(size(tempImage_mCh,1),size(tempImage_mCh,2)));
                    transImage = imtransform(image(:,:,z), tform2, 'XData', [1 (size(image(:,:,z),2)+tform2.tdata.T(3,1))],'YData', [1 (size(image(:,:,z),1)+ tform2.tdata.T(3,2))]);
                    transImageBuffrd(1:size(transImage,1),1:size(transImage,2)) = transImage;
                    tempImage_mCh(:,:,z) = transImageBuffrd((1:size(image(:,:,z),1)), (1:size(image(:,:,z),2)));
                end
                obj2.TIF(:,:,1,1:size(obj.TIF,4),t) = tempImage_mCh;    
            end
            
        end
        

    
    
%% Draw Localized Points
        function outputArg = drawpoint(obj,channel,filename,radius)

        obj.TIF = uint16(obj.TIF);
        if nargin==3
           radius=7; 
        end
        
        
        
         
        theta = linspace(0, 2*pi, round(4 * pi * radius));
        
        oriChannel = size(obj.TIF,3);
        
        for cc = 1 : max(size(channel))                
            for z = 1 : size(obj.TIF,4)
                for t = 1 : size(obj.TIF,5)
                    points = obj.out_cntrd{channel(cc),z,t};
                    if isempty(points)
                        break
                    end
                    myImage1 = zeros(size(obj.TIF,1,2),'uint16');
                    for j = 1:size(points,1)
                            xCenter = points(j,1);
                            yCenter = points(j,2);
                            x = radius * cos(theta) + xCenter;
                            y = radius * sin(theta) + yCenter;

                            for k = 1 : length(x)
                                row = round(y(k));
                                col = round(x(k));
                                myImage1(row, col) = 2^15-1;
                            end                        
                            
                    end
                    
                    obj.TIF(:,:,oriChannel+channel(cc),z,t) = myImage1;
                    
                end
            end
        end
        
        obj.save(filename);
        
        end
        
        
        
%% Draw Tracked Points
        function outputArg = drawpoint_track(obj,channel,channel_track,filename,radius)

        obj.TIF = uint16(obj.TIF);
        if nargin==4
           radius=7; 
        end
        
        
      
        
        theta = linspace(0, 2*pi, round(4 * pi * radius));
        
        oriChannel = size(obj.TIF,3);
        
        for cc = 1 : max(size(channel))                
            for z = 1 : size(obj.TIF,4)
                for t = 1 : size(obj.TIF,5)
                    points = obj.out_cntrd{channel(cc),z,t};
                    if isempty(points)
                        break
                    end
                      myImage1 = zeros(size(obj.TIF,1,2),'uint16');
                    for j = 1:size(points,1)
                            xCenter = points(j,1);
                            yCenter = points(j,2);
                            x = radius * cos(theta) + xCenter;
                            y = radius * sin(theta) + yCenter;

                            for k = 1 : length(x)
                                row = round(y(k));
                                col = round(x(k));
                                myImage1(row, col) = 2^15-1;
                            end                        
                            
                    end
                    
                    obj.TIF(:,:,oriChannel+channel(cc),z,t) = myImage1;
                    
                end
            end
        end
        
        for cc = 1 : max(size(channel_track))                
            for z = 1 : size(obj.TIF,4)
                for t = 1 : size(obj.TIF,5)
                    
                    if isempty(obj.out_track_plt)
                        myImage2 = zeros(size(obj.TIF,1,2),'uint16');
                        obj.TIF(:,:,oriChannel+max(size(channel))+cc,z,t) = myImage2;
                        continue;
                    end
                    points = obj.out_track_plt{channel_track(cc),t};
                    myImage2 = zeros(size(obj.TIF,1,2),'uint16');
                    for j = 1:size(points,1)
                            xCenter = points(j,1);
                            yCenter = points(j,2);
                            x = radius * cos(theta) + xCenter;
                            y = radius * sin(theta) + yCenter;

                            for k = 1 : length(x)
                                row = round(y(k));
                                col = round(x(k));
                                myImage2(row, col) = 2^15-1;
                            end                        
                            
                    end
                    
                    obj.TIF(:,:,oriChannel+max(size(channel))+cc,z,t) = myImage2;
                    
                end
            end
        end        
        obj.save(filename);
        
        end
        

        
        
        
%% Find Tracked Molecules

        function Track = findmol(obj,xyt,channel)
        
        xy = xyt([1 2]);
        cand = obj.out_track_plt{channel,xyt(3)}(:,:); 
        [M I] = min(vecnorm((cand(:,[1 2])-xy)'));
        I = cand(I,4);
        Track = obj.out_track{channel}(:,4)==I;
        Track = obj.out_track{channel}(Track,:);

        end



%%  End of methods        
    end
    
    
%%   Static methods
    methods(Static, Access = public)
        
        
            %% Load Stack
            function obj = loadStack(filename)
            
            tsStack = TIFFStack(filename);
            %TIFFStack : xytzc --> saveStacks : xyczt
            obj = MinStack;
%             obj.TIF = zeros(size(permute(tsStack(:,:,:,:,:),[1 2 5 4 3])));
            obj.TIF(:,:,:,:,:) = permute(tsStack(:,:,:,:,:),[1 2 5 4 3]);                    
                        
            
            end
            
            %% Draw Tube Plot
            
            function s = tubeplot(x,y,res,R,Rres,color,dt,tstart)
            
                
                
                
                
                % x = rand(1,30);
            % y = rand(1,30);
            % res = 25; R = 5; Rres = 25;    


            %%

            if nargin < 3, res = 25; R = 5; Rres = 25;  color = 'r'; dt = false; end
            if nargin < 4, R = 5; Rres = 25;  color = 'r'; dt = false; end
            if nargin < 5, Rres = 25;  color = 'r'; dt = false; end
            if nargin < 6, color = 'r'; dt = false; end
            if nargin < 7, dt = false; end
            if nargin < 8, tstart = false; end


            [X,trash] = meshgrid(x,1:Rres);
            [Y,trash] = meshgrid(y,1:Rres);
%             surface.u = linspace(0,100,max(size(x)));
            surface.u = 1:max(size(x));
            if tstart, surface.u = tstart:tstart+max(size(x))-1; end

            if dt, surface.u = surface.u *dt; end
            surface.v = linspace(0, 2*pi, Rres);     % Tube surface
            [surface.u, surface.v] = meshgrid(surface.u, surface.v);
            surface.x = X + R*cos(surface.v);
            surface.z = Y + R*sin(surface.v);
            surface.y = surface.u;


            %% Plot
%             figure()
            grid on
            hold on
            s = surf(surface.x, surface.y, surface.z)
            s.FaceColor = color;
            s.FaceAlpha = 0.2;
            s.EdgeColor = color;
            s.MeshStyle = 'column';
            view(70, 20);  
            axis equal
            xlabel('X');
            zlabel('Y');
            ylabel('Time (sec)');
            
            end

%% Fill missing frames in tracked molecules

            function Track = fillframe(Track)

            for t = Track(1,3) : Track(end,3)
                if ~max(Track(:,3)==t)
                   t
                   t = t - Track(1,3)+1;
                   Track = [Track(1:t-1,:); mean([Track(t-1,:);Track(t,:)]) ; Track(t:end,:)];
                end   

            end

            end
        
    end
    
    
end

