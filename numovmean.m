function output = numovmean(X,Y,xx,w,dim,shape)
% Y - input y-values
% X - timestamps for Y-data
% xx- timepoints to resample at (default X)
% w - window size (does not need to be integer)
% dim-dimension to run moving average (default 1, max = 3)
% shape-window shape for weights, can be set to rect or tri

narginchk(2,6)
if nargin<3 || isempty(xx)
    xx = X;
end
if nargin<4 || isempty(w)
    w = 2*mean(diff(xx));
end
if nargin<5 || isempty(dim)
    dim = 1;
end
if nargin<6 || isempty(shape)
    shape = 'rect';
end

if length(X)~=size(Y,dim)
    error('X and Y vectors must be the same length');
end

iter = numel(xx);
shape = lower(shape);
switch shape
    case 'rect'
        switch dim
            case 1
                output = NaN([iter size(Y,2:ndims(Y))]);
                if ~ismatrix(Y)
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        output(i,:,:) = mean(Y(idx,:,:),dim,'omitnan');
                    end
                else
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        output(i,:) = mean(Y(idx,:),dim,'omitnan');
                    end
                end
            case 2
                if ~ismatrix(Y)
                    output = NaN([size(Y,1) iter size(Y,3:ndims(Y))]);
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        output(:,i,:) = mean(Y(:,idx,:),dim,'omitnan');
                    end
                else
                    output = NaN([size(Y,1) iter]);
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        output(i,:) = mean(Y(:,idx),dim,'omitnan');
                    end
                end
            case 3
                output = NaN([size(Y,1,2) iter]);
                for i=1:iter
                    idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                    output(:,:,i) = mean(Y(:,:,idx),dim,'omitnan'); %#ok<*FNDSB>
                end
        end
    case 'tri'
        switch dim
            case 1
                output = NaN([iter size(Y,2:ndims(Y))]);
                if ~ismatrix(Y)
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        weight = w/2-abs(X(idx)-xx(i));
                        tot = sum(weight);
                        output(i,:,:) = sum(Y(idx,:,:).*weight,dim,'omitnan')/tot;
                    end
                else
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        weight = w/2-abs(X(idx)-xx(i));
                        tot = sum(weight);
                        output(i,:) = sum(Y(idx,:).*weight,dim,'omitnan')/tot;
                    end
                end
            case 2
                if ~ismatrix(Y)
                    output = NaN([size(Y,1) iter size(Y,3:ndims(Y))]);
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        weight = w/2-abs(X(idx)-xx(i));
                        tot = sum(weight);
                        output(:,i,:) = sum(Y(:,idx,:).*weight,dim,'omitnan')/tot;
                    end
                else
                    output = NaN([size(Y,1) iter]);
                    for i=1:iter
                        idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                        weight = w/2-abs(X(idx)-xx(i));
                        tot = sum(weight);
                        output(i,:) = sum(Y(:,idx).*weight,dim,'omitnan')/tot;
                    end
                end
            case 3
                output = NaN([size(Y,1,2) iter]);
                for i=1:iter
                    idx = find(X>=(xx(i)-w/2) & X<=(xx(i)+w/2));
                    weight = (w/2-abs(X(idx)-xx(i)));
                    tot = sum(weight);
                    output(:,:,i) = sum(Y(:,:,idx).*reshape(weight,1,1,numel(idx)),dim,'omitnan')/tot;
                end
        end
end