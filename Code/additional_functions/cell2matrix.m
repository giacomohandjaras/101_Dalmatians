%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% script by Giacomo Handjaras, Francesca Setti %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function matrix=cell2matrix(cell_data)

dim_matrix=numel(cell_data{1});
dim_cell_data=size(cell_data,1);

matrix=nan(dim_cell_data,dim_matrix);

for i=1:dim_cell_data
    temp_matrix=cell_data{i};
    if isempty(temp_matrix)
        temp_matrix=nan(1,dim_matrix);
    end
    
    matrix(i,:)=temp_matrix(:);
end


end
