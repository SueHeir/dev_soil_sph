# dev_soil_sph Example Validation

See [`../validation/README.md`](../validation/README.md) for the full validation
set. Result figures committed under individual examples are embedded there and in
the example READMEs.

## Footpad Bearing/Sinkage

[`footpad`](footpad/README.md) compares the seated force-sinkage branch of a
driven SPH footpad against the Bekker/Wong pressure-sinkage form as independently
validated by DIRT's DEM plate-sinkage benchmark. The same oracle rejects the
zero-gravity control.

![footpad pressure-sinkage validation](footpad/plots/footpad_bekker_validation.png)
