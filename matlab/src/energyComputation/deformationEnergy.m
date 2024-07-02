function [E] = deformationEnergy(shape1, shape2, combsAll)
    % energy is computed with feature energy function
    curr_featureEnergy = featureEnergy(shape1, shape2, ...
        combsAll(:, 1:3), ...
        combsAll(:, 4:6));
    E = curr_featureEnergy;
end
