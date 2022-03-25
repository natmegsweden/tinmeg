


a = [1 2 1 1 2 1];

for i = 1:numel(a)
    
    if i == 1
        disp(a(i))
    elseif i > 1
        
        if a(i-1) == 1
            disp(a(i)+1);        
        elseif a(i) == 1
            display(a(i))
        end
        
    end
    
end