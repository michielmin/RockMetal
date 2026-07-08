The RockMetal subroutine and jupyter notebook provide a simple and fast internal structure model for cool, rocky exoplanets.
The physical model is a simple three-layered structure with an iron core, a silicate mantle and an h2o top-layer.
The values of f_core and f_h2o can be set to compute the radius corresponding to a given mass.
References for the parameters of the EOS are given in the source code.

f_core: this is the ratio between the iron core and the total mass in iron and silicate
f_h2o: this is the ratio between the h2o layer and the total mass of the planet

For Earth the values should be roughly:
f_core=0.325
f_h2o=0.005

This gives R=0.99 for M=1.0.

The model is a crude approximation to internal structure computations and especially useful for first estimates or retrieval studies.
When more accuracte values are needed, we recommend using a more complete internal structure model.
