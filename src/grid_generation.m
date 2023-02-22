for n=1:length(blocks)
    eval(['[' blocks{n} ']=block_grid(' blocks{n} ', mesh, flowopt);'])
    
end