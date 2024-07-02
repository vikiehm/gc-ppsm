function e = featureEnergy(shape1, shape2, FaCombo, FbCombo)
    Feata = cell2mat(shape1('features'));
    Featb = cell2mat(shape2('features'));
    e = get_energy_vector(Feata, Featb, FaCombo, FbCombo);
end

function e = get_energy_vector(Feata, Featb, FaCombo, FbCombo)
    e = zeros(size(FaCombo, 1), 1);
    for i = 1:size(FaCombo, 1)
        currComboA = FaCombo(i, :);
        currComboB = FbCombo(i, :);
        currFeatA = Feata(currComboA, :);
        currFeatB = Featb(currComboB, :);
        e(i) = mean(vecnorm((currFeatA - currFeatB), 2, 2));
    end
end
