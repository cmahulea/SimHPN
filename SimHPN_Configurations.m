%SimHPN_Configurations - List of configurations of a given TCPN
%
%   This function returns a list of the possible configurations of a TCPN
%   in form of a vector. The number that appears in the entry conf(i)
%   indicates that the place p_(conf(i)) constrains the flow of transition
%   t_i
%
%   conf = Configurations(Pre)
function conf = SimHPN_Configurations(Pre)

[num_p, num_t]=size(Pre);

for i=1:num_t                          %% Obtaining the places before each
    pre_t{i}=[];                       %% transition (pre_t).
    for j=1:num_p
        if Pre(j,i)>0
            pre_t{i}=[pre_t{i} j];
        end
    end
end

conf = cartesian(pre_t{:});