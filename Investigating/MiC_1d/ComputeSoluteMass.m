function mass = ComputeSoluteMass(zMark, cMark, ModelDim)
    mass = sum(ComputeNodalValues(zMark, cMark, ModelDim) .* abs(ModelDim.dzin));
end