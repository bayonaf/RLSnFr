function exportToTIF(imgdata, selpath, folder, fileName)
%EXPORTTOTIFF Summary of this function goes here
%   Detailed explanation goes here
    
    mkdir(selpath + folder)

    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    
    if iscell(imgdata)
        disp("Exporting " + folder + "... ")
        tic
        f = waitbar(0, ' Exporting " + folder + " ... 0%');
        for s = 1: length(imgdata)
            tagstruct.ImageLength = size(imgdata{s},1);
            tagstruct.ImageWidth = size(imgdata{s},2);
            t = Tiff(selpath + folder + fileName + (s-1) + ".tif",'w');
            setTag(t,tagstruct)
            write(t,single(imgdata{s}(:,:,1)));
            close(t);

            for i = 2 : size(imgdata{s},3)
                t = Tiff(selpath + folder + fileName + (s-1) + ".tif",'a');
                setTag(t,tagstruct)
                write(t,single(imgdata{s}(:,:,i)));
            end
            close(t);
            waitbar(s/length(imgdata), f, "Exporting " + folder...
                +" ...(" +string(fix(s*100/length(imgdata)))+"%)"); % Progress Bar
        end
        close(f)
        elapsedTime = toc;
        disp("All " + folder + " exported. Elapsed time> "+ string(elapsedTime)+ "Seconds")
    else
        disp("Exporting " + folder + "... ")
        tic
        tagstruct.ImageLength = size(imgdata,1);
        tagstruct.ImageWidth = size(imgdata,2);
        
        t = Tiff(selpath + folder + fileName + ".tif",'w');
        setTag(t,tagstruct)
        write(t,single(imgdata(:,:,1)));
        close(t);

        for i = 2 : size(imgdata,3)
            t = Tiff(selpath + folder + fileName + ".tif",'a');
            setTag(t,tagstruct)
            write(t,single(imgdata(:,:,i)));
        end
        close(t);
        elapsedTime = toc;
        disp("All " + folder + " exported. Elapsed time> "+ string(elapsedTime)+ "Seconds")
    end
    
end

