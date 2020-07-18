function R = get_R()
    num_states = 14;
    num_reactions = 34;
      R = sparse(num_states, num_reactions);
    % returns a matrix of net changes in species levels caused by firing of
    % the reactions 
    % (1st dimension: states, 2nd dimension: reaction)
    
   % e.g. R(2,20) = -1 represents that the second state is changed by -1
   % units when the second reaction is fired.
    R(1,1)=1;
    R(2,2)=1;
    R(3,3)=1;
    R(4,4)=1;
    R(5,5)=1;
    R(6,6)=1;
    R(7,7)=1;
    R(8,8)=1;
    R(1,9)=-1;
    R(2,10)=-1;
    R(2,11)=-2;
    R(13,11)=1;
    R(13,12)=-1;
    R(2,12)=2;
    R(2,13)=-1;
    R(3,13)=-1;
    R(10,13)=1;
    R(10,14)=-1;
    R(2,14)=1;
    R(3,14)=1;
    R(1,15)=-2;
    R(9,15)=1;
    R(3,16)=-1;
    R(3,17)=-2;
    R(14,17)=1;
    R(14,18)=-1;
    R(3,18)=2;
    R(9,19)=-1;
    R(11,20)=-1;
    R(12,21)=-1;
    R(13,22)=-1;
    R(10,23)=-1;
    R(14,24)=-1;
    R(9,25)=-1;
    R(1,25)=2;
    R(4,26)=-1;
    R(1,27)=-1;
    R(2,27)=-1;
    R(11,27)=1;
    R(5,28)=-1;
    R(11,29)=-1;
    R(1,29)=1;
    R(2,29)=1;
    R(6,30)=-1;
    R(1,31)=-1;
    R(3,31)=-1;
    R(12,31)=1;
    R(7,32)=-1;
    R(12,33)=-1;
    R(1,33)=1;
    R(3,33)=1;
    R(8,34)=-1;
    
end