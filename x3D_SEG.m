%Single color (red) object detection and stack-organization in Matlab

tic
LayerDistance = 3; %distance of layers from each other, used in the volumetric model (at the end)

seg_input = 'G:\test_segmentation_out\SEGTLvsAB_2\'; %3D segmentation input, use the parent folder containing all of the 2D segmentation time point folders
cd(seg_input);

files = dir;
N = length(files);
set(0,'DefaultFigureRenderer','opengl'); %utilise gpu
parfor a = 3:N
    BASEPATH = seg_input;
    files = dir;
    tic
    currD = files(a).name;
    disp(currD);
    BASEPATH = strcat(BASEPATH,currD,'\');
    FILES = dir([BASEPATH '*.tif']); %haseeb set as tif

    %Initializing image dimensions and matrix
    FILENAME = FILES(1).name;
    A = (imread([BASEPATH FILENAME])); 
    Rows = size(A,1); 
    Cols = size(A,2);
    Layers = size(FILES,1);
    disp([Rows,Cols])
    M = zeros(Rows,Cols,Layers); %the Matrix of virtual stack
    disp(['Dimensions: ' num2str(Rows) ' x ' num2str(Cols) ' x ' num2str(Layers)]);

    for i = 1:size(FILES,1)

            FILENAME = FILES(i).name;
            disp(['Processing image #' num2str(i) ' of ' num2str(size(FILES,1)) '...']);

            A = imread([BASEPATH FILENAME]);

    
            bw4_perim = bwperim(A);
    %        overlay1 = imoverlay(I_eq, bw4_perim, [1 0 0]); %[.3 1 .3]

            %imwrite(overlay1,[BASEPATH FILENAME(1:end-4) '_det.tif']); %haseeb
            %comment out as this is breaking the output
            
            temp_name = tempname
            M(:,:,i) = bw4_perim;
            %writematrix(M(:,:,i), temp_name);

    %        figure,imshow(overlay1)
                %set(gcf,'Position',[5 1 1000 1000]);      %Enlarge figure and axes as much as possible
                %set(gca,'Position',[0.01 0.01 0.98 0.98]);


    %         stats = regionprops(bw4_perim, 'BoundingBox','Centroid','Area','MajorAxisLength','MinorAxisLength');
    % 
    % %        figure,imshow(A);
    %         hold on
    % 
    %         set(gcf,'Position',[5 1 1000 1000]);      %Enlarge figure and axes as much as possible
    %         set(gca,'Position',[0.01 0.01 0.98 0.98]);
    % 
    % %        disp([num2str(size(stats,1)) ' objects found']);
    % 
    %         % Stats below
    % 
    % 
    %         for object = 1:length(stats)
    %             bb = stats(object).BoundingBox;
    %             bc = stats(object).Centroid;
    %             diameters = mean([stats(object).MajorAxisLength stats(object).MinorAxisLength],2);
    %             radii = diameters/2;
    %             viscircles(bc,radii,'LineWidth',1,'LineStyle','-');
    %             %rectangle('Position',bb,'EdgeColor','r','LineWidth',1)
    %             plot(bc(1),bc(2), '-m+')
    %             %a=text(bc(1)+15,bc(2), strcat('X: ', num2str(round(bc(1))), '    Y: ', num2str(round(bc(2)))));
    %             %set(a, 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 12, 'Color', 'yellow');
    %             %pause(0.01); drawnow;
    %         end

    end
    FILENAME = [BASEPATH 'stack.mat'];
    %save(FILENAME,'M','-v7.3'); %haseeb v7.3 for files >2GB save param
    %stack  as table
    %disp(['Processing took ' num2str(ceil(toc/60)) ' minutes.']);
    disp(['improc done'])
    disp('Creating 3D Virtual Stack...');
    %tic;
    %--3D OUTPUT VISUALISATION BLOCK--
    %{
    fig=figure; hold on;
    sort(Layers); %multithreaded parallelisation disorders processed images so reorder before drawing them. Takes about the same time to do this.
    %no need to parallelise next steps as this is not rate limiting.
    for l=1:Layers

        disp(['Processing Layer #' num2str(l) ' of ' num2str(Layers) '...']);

        for r=1:Rows %%%%%%%%%%%%%%%%%%%%%%%%%%
            idx = find(M(r,:,l)>0);

            if ~isempty(idx)
                for i=1:length(idx)
                    plot3(idx(i),r,l,'.'); %drawnow;
                end
            end


        end


    end
    set(gcf, 'renderer', 'opengl') %haseeb perhaps use gpu for figure rendering??
    view(30,30);
    grid on;
    disp('Processing... Please wait...');
    FIGFILENAME = [BASEPATH 'stack.fig'];
    %saveas(fig, FIGFILENAME, 'fig'); %%%%%%%%%%%%%%%%%%%%%

    figure;
    %}
    %--END 3D OUTPUT VISUALISATION BLOCK--
    MCC = bwconncomp(M);
    stats3D = regionprops3(MCC,'all');
    %relevant_props3 = table(stats3D.Volume,stats3D.Centroid(:,1),stats3D.Centroid(:,2),stats3D.Centroid(:,3), stats3D.BoundingBox(:,1),stats3D.BoundingBox(:,2),stats3D.BoundingBox(:,3),stats3D.Solidity, stats3D.MaxIntensity, stats3D.MeanIntensity, stats3D.MinIntensity, stats3D.WeightedCentroid, stats3D.Extent, stats3D.EquivDiameter)
    %relevant_props3 = table(["Volume","X_Centroid","Y_Centroid","Z_Centroid","Bounding_Box_X","Bounding_Box_Y","Bounding_Box_Z","Solidity","Extent","EquivDiameter"], [stats3D.Volume,stats3D.Centroid(:,1),stats3D.Centroid(:,2),stats3D.Centroid(:,3), stats3D.BoundingBox(:,1),stats3D.BoundingBox(:,2),stats3D.BoundingBox(:,3),stats3D.Solidity, stats3D.Extent, stats3D.EquivDiameter])
    relevant_props3 = table(stats3D.Volume,stats3D.Centroid(:,1),stats3D.Centroid(:,2),stats3D.Centroid(:,3), stats3D.BoundingBox(:,1),stats3D.BoundingBox(:,2),stats3D.BoundingBox(:,3),stats3D.BoundingBox(:,4),stats3D.BoundingBox(:,5),stats3D.BoundingBox(:,6),stats3D.Solidity, stats3D.Extent, stats3D.EquivDiameter)
    csvname = strcat("stats3D",currD,".csv");
    writetable(relevant_props3,csvname);
    %stats2D = regionprops(MCC);
    %scatter3sph(stats3D.Centroid(:,1),stats3D.Centroid(:,2),stats3D.Centroid(:,3),'size',stats3D.EquivDiameter * LayerDistance,'trans',0.3);
    %grid on;
    title([num2str(size(stats3D,1)) ' objects identified']);
    %disp(['All done. Building model took ' num2str(ceil(toc/60)) ' minutes. 3D Figure saved.']);
    disp(['all done']);
    toc
    %clear
    clc
    
    tic

    LayerDistance = 3; %distance of layers from each other, used in the volumetric model (at the end)

    
    cd(seg_input);
    BASEPATH = seg_input;
    files = dir;
    N = length(files);
end
toc

