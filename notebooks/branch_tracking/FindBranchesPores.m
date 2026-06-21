%This is (hopfeully) the final version of the code to find the branches.
%This code finds all branches and saves them. Calculating the diameter of
%the branch is done in another program.

function [network_edges,network_vertices,branch,imbranch] = FindBranchesPores(basepath, skel,BW)

    im=skel;
    
    %some branches might touch the border of the image - to avoid problems, we set the border of the image to 0.
    im(1,:)=0;
    im(end,:)=0;
    im(:,1)=0;
    im(:,end)=0;

    im=(im>0);
    
    %there are sometimes isolated pixels, remove them
    cc=bwconncomp(im);
    if (cc.NumObjects>1)
        fprintf(1, '%d objects, keeping only the biggest\n', cc.NumObjects);
        nbpix=cellfun(@numel, cc.PixelIdxList);
        [~, idx]=sort(nbpix, 'descend');
        mask=zeros(size(im));
        mask(cc.PixelIdxList{idx(1)})=1;
        im=mask>0;
    end
    
    im2=im;
    siz_im=size(im);
    [~,~]=find(im2);

    %find endpoints
    bwend=bwmorph(im, 'endpoints');
    endptslist=find(bwend);
    %find branchpoints
    bwbranch=bwmorph(im, 'branchpoints');
    branchptslist=find(bwbranch);
    
    %if there are no branchpoints (single vein), put the end points in the
    %branch points and empty the endpoints
    if (isempty(branchptslist))
        branchptslist=endptslist;
        endptslist=[];
    end
    
    %branchptslist = [branchptslist;endptslist];
    
%     look for cluster of touching branchpoints, and keep only the centroid
%     of each branchpoint
%     [rb, cb]=ind2sub(siz_im, branchptslist);
%     imbp=zeros(siz_im);
%     imbp(branchptslist)=1;
%     imbp=imbp>0;
%     cc=bwconncomp(imbp);
%     bplist2=[];
%     for ii=1:cc.NumObjects
%         if (length(cc.PixelIdxList{ii})==1)
%             bplist2=[bplist2; cc.PixelIdxList{ii}];
%         else
%             ct= ceil(cc.PixelIdxList{ii}/siz_im(1));
%             rt= cc.PixelIdxList{ii}-(ct-1)*siz_im(1);
%             meanc=mean(ct);
%             meanr=mean(rt);
%             [~, miind]=min(sqrt((rb-meanr).^2+(cb-meanc).^2));
%             rt=rb(miind);
%             ct=cb(miind);
%             bplist2=[bplist2; (ct-1)*siz_im(1) + rt]; 
%         end
%     end
%         
%     branchptslist=bplist2;

    
      
   
    
    % make a matrix with the location of all branchpoints,
    % where the value of the branchpoint is the number of branches that
    % leave this branchpoint
    matbpts=zeros(size(im));
    matbpts(branchptslist)=1;
    torem=[];
    for imbp=1:length(branchptslist)
        [rb, cb]=ind2sub(size(im), branchptslist(imbp));
        %neighborhood of the branchpoint
        tmpim=im(rb-1:rb+1, cb-1:cb+1);
        %number of branches that leave this branchpoint
        matbpts(branchptslist(imbp))=sum(sum(tmpim))-1;
        
        %if the branchpoint is an isolated pixel between 2 branches, remove
        %it
        if ( sum(sum(tmpim))==3)
            torem=[torem; imbp];
        end
        
    end
    savematbpts = matbpts;
    branchptslist(torem)=[];
    
    [rb, cb]=ind2sub(size(im2), branchptslist);

    %initialize matrices with vertices
    vert=zeros(0,4); %will have the coordinates of the vertices for each branch
    vertind=zeros(0,2); %will have the indices of the vertices for each branch

    %initalize cell with branches
    branch=cell(0);
    ibranch=0;
    
    %scan all branchpoints 
    for ibpl=1:1:length(branchptslist)
        [rib, cib]=ind2sub(size(im2), branchptslist(ibpl));
        im2(rib, cib)=0;
        
        ri=rib;
        ci=cib;
        
        
        %while there are branches going out of the vertex, continue looking
        %for them
        while(matbpts(branchptslist(ibpl))>0)
            
            nohit=1;
            ibranch=ibranch+1;
            branch{ibranch}=[];
            branch{ibranch}=[branch{ibranch}; ri ci branchptslist(ibpl)];

            %save beginning of branch
            vert(ibranch,1:2)=[ri, ci];
            vertind(ibranch,1)=sub2ind(size(im2), ri, ci);
            wrong=0;
            weird=0;
            
            
            
            % while we're on the same branch, look for neighboring point
            while(nohit)
%                 if (ibranch==274)
%                 figure(1)
%             imagesc(im2(ri-10:ri+10, ci-10:ci+10)+im(ri-10:ri+10, ci-10:ci+10))
%             hold on
%             plot(11, 11, 'ok')
%             hold off
%             pause
%                 end
                
                %walk along one of the veins
                tmpim=im2(ri-1:ri+1, ci-1:ci+1); %neighborhood of the point

                
                %if there are no points in the neighborhood,
                %it means that we're close to a branchpoint that has been
                %previously visited (and therefore erased from im2)
                if (sum(tmpim(:))==0)
                    %look at the distance between the current point and
                    %all branchpoints
                    d=sqrt((ri-rb).^2+(ci-cb).^2);
                    %if the current point is a branchpoint, don't consider
                    %it
                    indz=find(d==0);
                    if (~isempty(indz))
                        d(indz)=NaN; 
                    end
                    %sort the distances
                    [d,ixd]=sort(d);
                    if (d(1)<2)
                        %if there is a branchpoint close enough (1 pixel
                        %away in horizontal or vertical or diagonal), then go
                        %to it.
                        indi=branchptslist(ixd(1));
                        [ri,ci]=ind2sub(size(im2), indi);
                    else 
                        if (sum(im2(:))~=0)
                            if (size(branch{ibranch},1)>1)
                        fprintf(1, 'something is wrong\n');
                        fprintf(1, 'ibranch=%d\t ri=%d\t ci=%d\n', ibranch, ri, ci);
                        fprintf(1, 'distance to closest branchpoint:%f\n', d(1));
                        fprintf(1, 'current size of the branch:%d\n', size(branch{ibranch},1));
                        figure;imagesc(im(ri-10:ri+10,ci-10:ci+10)+im2(ri-10:ri+10, ci-10:ci+10)+...
                            matbpts(ri-10:ri+10, ci-10:ci+10));
                        %pause
                            end
                        nohit=0;
                        wrong=1;
                        else
                            fprintf(1, 'imnb=%d\t , something is weird but not too much, ibranch=%d\n', imnb, ibranch);
                            weird=1;
                            if (sum(im2(:))==0)
                                nohit=0;
                            end
                        end
%                        nohit=0; %we have hit a branchpoint in diagonal
%                         if (length( branch{ibranch}(:,1)==1))
%                             branch{ibranch}=[];
%                             ibranch=ibranch-1;
%                         else
%                             'what happens????'
%                             pause
%                         end
                    end
                    
                    
                else
                    
                    %otherwise we find the points in the neighborhood
                    [ri_n,ci_n]=find(tmpim);
                    
                    if (length(ri_n)>1)
                        %calculate the distance between current point and
                        %neighboring points
                        d=sqrt((ri_n-2).^2+(ci_n-2).^2);
                        %sort by ascending order
                        %[ds, ids]=sort(d);
                        
                        %find all the closest pixels
                        ids=find(d==min(d));
                        
                        %if there are more than 1 possibility, avoid
                        %cutting off another branch from a branchpoint 
                        if (length(ids)>1)
                            %possible pixels
                            riposs=zeros(length(ids),1);
                            ciposs=zeros(length(ids),1);
                            %indiposs=zeros(length(ids),1);
                            for jjj=1:length(ids)
                                riposs(jjj) = ri+(ri_n(ids(jjj))-2);
                                ciposs(jjj) = ci+(ci_n(ids(jjj))-2);
                            end
                            indiposs= siz_im(1)*(ciposs-1)+riposs; %convert to index 

                            %if there is a possibility of ending the branch
                            %at a branch point, choose it
                            Lia=ismember(indiposs, branchptslist);
                            % but if already visited previously then go
                            % ahead and use another one
                            if (sum(Lia)) && sum(ismember(vert,[rib,cib,riposs(Lia),ciposs(Lia)],'rows')) == 0
                                india=find(Lia);
                                %in case there are several branchpoints, we just choose one
                                india=india(1); 
                                ri_n=ri_n(ids(india));
                                ci_n=ci_n(ids(india));
                           
                            else
                            
                                dposs=zeros(length(ids),1);
                                for jjj=1:length(ids)
                                    dposs(jjj)=min(sqrt( (riposs(jjj)-rb).^2 + (ciposs(jjj)-cb).^2));
                                end
                               
                                %sort
                                [dd, idd]=sort(dposs);
                                %choose the 1st pixel that is not
                                %4-connected, if it exists
                                indd=find(dd>1);
                                if (~isempty(indd))
                                    okd=idd(indd(1));
                                    ri_n=ri_n(okd);
                                    ci_n=ci_n(okd);
                                else
                                    %just pick one arbitrarily
                                    ri_n=ri_n(1);
                                    ci_n=ci_n(1);
                                end
                            end
                      
                        
                        else %if there's only one closest pixel
                            [~, id]=min(d);
                            ri_n=ri_n(id);
                            ci_n=ci_n(id);
                        end
                        
                        
                    end
                        
                    ri=ri+(ri_n-2);
                    ci=ci+(ci_n-2);
                    indi=siz_im(1)*(ci-1)+ri; %convert to index - faster than sub2ind
                   
                end
                
                if (~wrong && ~weird)
                %add this point to the current branch
                branch{ibranch}=[branch{ibranch}; ri ci indi];
                end
                
                %if we hit a branchpoint, increment branch nb 
                if (sum(branchptslist==indi))
                    nohit=0;
                    %save that the branch has finished there
                    vert(ibranch,3:4)=[ri, ci];
                    vertind(ibranch,2)=indi; %sub2ind(size(im2), ri, ci);
                   
                    %subtract one branch from the other end
                    matbpts(indi)=matbpts(indi)-1;
                    %reinitialize to original branchpoint
%                     [ri,ci]=ind2sub(size(im2), branchptslist(ibpl));
                    ri=rib;
                    ci=cib;

                else
                    %remove this point from the image
                    im2(ri,ci)=0;
                end
                
                %if we hit an endpoint, increment branch nb
                if (sum(endptslist==indi))
                    nohit=0;
                    %save that the branch has finished here
                    vert(ibranch,3:4)=[ri, ci];
                    vertind(ibranch,2)=sub2ind(size(im2), ri, ci);
                    %reinitialize to original branchpoint
                    %[ri,ci]=ind2sub(size(im2), branchptslist(ibpl));
                    ri=rib;
                    ci=cib;
                end
                
                if (wrong)
                    if (size(branch{ibranch},1)>1)
               %find the closest branchpoint to the last pixel of the branch 
                rici=branch{ibranch}(end,:);
                ri=rici(1);
                ci=rici(2);
                indi=siz_im(1)*(ci-1)+ri;
                dist=(rb-ri).^2+(cb-ci).^2;
                [~, mind]=min(dist);
                branch{ibranch}=[branch{ibranch}; rb(mind) cb(mind)  siz_im(1)*(cb(mind)-1)+ rb(mind)];
                vert(ibranch,3:4)=[rb(mind) cb(mind)];
                vertind(ibranch,2)=siz_im(1)*(cb(mind)-1)+ rb(mind);
                    end
                end
                
                %do not take into account branches with 1 point...
                if (nohit==0)
                    u=unique(branch{ibranch}(:,3));
                    if (length(u)==1)
                        %reinitialize
                         branch{ibranch}=[];
                        vert(ibranch,:)=[];
                        vertind(ibranch,:)=[];
                        ibranch=ibranch-1;
                        
                    end
                end
                
               
                
            end
            
            matbpts(branchptslist(ibpl))=matbpts(branchptslist(ibpl))-1;
            
        end
        
        %when all branches from the branchpoint have been found, put the
        %pixel to 0
        im2(branchptslist(ibpl))=0;
        
    end
    

    %we have found all branches, now save their pixel positions, and also
    %save the skeleton as a color image, where color codes for the branch
    %number
    totnbbr=ibranch;
    imtest=zeros(size(im));
    imbranch = zeros(size(im));
    im2 = 5*totnbbr*double(BW);
    
    abr=[];
    
    for i=1:totnbbr
        mybranch=branch{i};
        %add branch to image
        im2(mybranch(:,3))= i;
        imbranch(mybranch(:,3)) = i;
        imtest(mybranch(:,3)) = 1;
        %save branch
        if (~isempty(abr))
            s1=size(abr,1);
            s2=size(mybranch,1);
            if (s1<s2)
                abr=[abr; zeros(s2-s1,size(abr,2))];
            else
                mybranch=[mybranch; zeros(s1-s2,3)];
            end
        end
        abr=[abr mybranch];
    end          
    
    %save all branches in a big file at each time step
%     name2save=strcat(basepath, 'allbranches.txt');
%     sib=size(abr,2);
%     fid=fopen(name2save, 'w');
%     strf=repmat('%d\t', [1 sib]);
%     strf=sprintf('%s\n', strf);
%     fprintf(fid, strf, abr');
%     fclose(fid);
    h = figure;imagesc(im2);colormap(gray);axis equal
    name2save=strcat(basepath, '_skeleton.png');
    saveas(h,name2save)
    close
    
    %save matrix
%     im2=uint8(im2);
%     cm=jet(8);
%     cm=repmat(cm,[32 1]);
%     cm(1,:)=[1 1 1];
%     axis equal;
%     name2save=strcat(basepath, '_skeleton.png');
%     imwrite(im2, cm, name2save, 'png');
    
%     figure;imagesc(savematbpts);axis equal
    % make sure that all branches have been traced
%     figure;imagesc(im);axis equal
figure;imagesc(im-imtest);axis equal
    assert(sum(sum(im-imtest))==0)
    
    

    %at each time step, save information about the connectivity
    vertices=unique(vertind);
    %make it a column vector if it's not already
    if (size(vertices,2)>size(vertices,1))
        vertices=vertices';
    end
    network_edges=zeros(size(vertind,1), 2);

    for i=1:length(vertices)
        [x,y]=find(vertind==vertices(i));
        ind=sub2ind([size(vertind,1), 2], x, y);
        network_edges(ind)=i;
    end

    network_vertices=zeros(length(vertices), 2);
    for i=1:length(vertices)
        [x,y]=find(vertind==vertices(i));
        %we don't need the multiple indices
        x=x(1);
        y=y(1);
        %get the vertex position - remember that the data is stored like
        %this:
        %**** vert(ibranch,3:4)=[ri, ci];
        %**** vertind(ibranch,2)=sub2ind(size(im2), ri, ci);
        
        vpos=vert(x,2*(y-1)+1:2*y)  ;      
        network_vertices(i,:)=[vpos(1) vpos(2)];
    end
    
    % delete edges between branch points occuring twice within the network
    i = 0;
    while(i < size(network_edges,1))
        i = i + 1;
        if size(branch{i},1) <3
            % detect loops
            [~,seso] = find(ismember(network_edges(i+1:end,:),network_edges(i,:),'rows')>0);
            seso = seso + i;
            [~,sero] = find(ismember(network_edges(1:end,:),fliplr(network_edges(i,:)),'rows')>0);
            % if loop detected delete original edge
            if ~isempty(seso) || ~isempty(sero)
                network_edges(i,:) = [];
                branch(i) = [];
                i = i - 1;
            end
        end
    end

end
