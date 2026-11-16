# Provenance and evidence limits

Run 0: single loudspeaker, absolute lateral position not independently recorded.
Run 1: sources (-0.20,0) and (0.00,0) m. Run 2: (-0.20,0) and (+0.20,0) m.
Run 3: two-bladed propeller, absolute lateral position not independently recorded.
Run identities began with user recall, supported by diagnostics; actual dual-source
geometry was confirmed by user and teammate. It is not inferred from the map
then reused as its own localization truth. F7B a,d are qualitative schematics.

Metadata: 56 channels; 500000 samples; fs≈50000.0000101864 Hz (floating timestamp
representation), nominal 50000 Hz; about 10 s. Sound Pressure_0..55 map in numeric
stored order to coordinate rows 1..56. Independent physical wiring and viewing
direction are unverified. The nominal coordinate system is x horizontal,
y vertical and z toward the scan plane. Recorded units are Pascals and NI
reference 20 µPa, but sensitivity calibration is not independently certified.
These limitations affect SPL accuracy and nominal localization interpretation.

The original audit found prefixed `._*` files to be AppleDouble sidecars, not
alternate code. The supplied CBF code includes edits of unknown earlier origin;
it is the actual local audited reference, not a claim of a pristine TA version.
The legacy fft_p uses peak amplitude; fft_p_blockaverage has a bounds/coherent-
averaging defect. Neither replaces the CBF block processing in this package.
Original TDMS/MAT spotchecks and baseline evidence are summarized in QA.

Only project-owned MATLAB/Python helpers and derived scientific arrays are
published. Original recordings, auxiliary FT/Motor files, teacher code, handout,
slides and original coordinate text remain local and read-only. Input hashes
identify evidence without redistributing those materials. No report text or
Section 8 figures are generated or included.
