function compiledSchnitzCells = initializeParticleFields(compiledSchnitzCells,nucleusCounter,...
                                                          n_entries,has3DSpotInfo,hasAPInfo)

    compiledSchnitzCells(nucleusCounter).particleID = NaN;
    compiledSchnitzCells(nucleusCounter).xPosParticle = NaN(1,n_entries);
    compiledSchnitzCells(nucleusCounter).yPosParticle = NaN(1,n_entries);
    compiledSchnitzCells(nucleusCounter).zPosParticle = NaN(1,n_entries);      
    
    if has3DSpotInfo
      compiledSchnitzCells(nucleusCounter).xPosParticle3D = NaN(1,n_entries);
      compiledSchnitzCells(nucleusCounter).yPosParticle3D = NaN(1,n_entries);
      compiledSchnitzCells(nucleusCounter).zPosParticle3D = NaN(1,n_entries); 
      compiledSchnitzCells(nucleusCounter).fluo3D = NaN(1,n_entries);
    end
    if hasAPInfo 
      compiledSchnitzCells(nucleusCounter).APPosParticle = NaN(1,n_entries);
    end

    compiledSchnitzCells(nucleusCounter).zPosParticleDetrended = NaN(1,n_entries);

    compiledSchnitzCells(nucleusCounter).fluo = NaN(1,n_entries); 
    compiledSchnitzCells(nucleusCounter).fluoOffset = NaN(1,n_entries); 
    compiledSchnitzCells(nucleusCounter).spotFrames = NaN(1,n_entries); 