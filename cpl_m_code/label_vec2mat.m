function Z = label_vec2mat(z,varargin)

if nargin > 1
    K = varargin{1};
else
    K = max(z);
end

temp = eye(K);
Z = temp(z(:),:);

