-- SQLite statements for setting up a "New TIMS data format" database.
--

-- List of global metadata describing this analysis.
-- 
-- There is a core set of "essential global metadata" that is always set in a valid
-- analysis.
--
CREATE TABLE GlobalMetadata (
    Key TEXT PRIMARY KEY,
    Value TEXT
    );
-- This table defines a set of "generic properties" that are used in this analysis.
--
-- The acquisition engine can define a generic property to attach arbitrary metadata to
-- frames (or to frame groups such as "segments") that are not an integral part of the
-- officially documented file format schema.
--
CREATE TABLE PropertyDefinitions (
    -- Unique ID for this property. Property IDs may change between versions of
    -- otofControl; always use the 'PermanentName' of a property as a handle.
    Id INTEGER PRIMARY KEY, 

    -- Name of the generic property.
    PermanentName TEXT NOT NULL, 

    -- Indicates the type of the data stored for this property:
    -- 
    --    0 = int32             (SQLite type: INTEGER)
    --    1 = double            (SQLite type: REAL)
    --    2 = string            (SQLite type: TEXT)
    --   10 = array of int32    (SQLite type: BLOB; containing array of little-endian int32)
    --   11 = array of double   (SQLite type: BLOB; containing array of little-endian
    --                                              standard x86 floating-point double /
    --                                              IEEE-754 binary64)
    --   12 = array of string   (SQLite type: BLOB; containing zero-terminated UTF-8 strings)
    --
    -- When a property has a value of NULL, it means it hasn't been set at all. This is
    -- different from a TEXT or BLOB of length zero, which represent empty strings or
    -- vectors, respectively.
    --
    -- Notes: Since SQLite doesn't have native array
    -- types, we serialize the array types into BLOBs to avoid the hassle associated with
    -- having extra "offload" tables to store the individual array entries.
    --
    Type INTEGER NOT NULL,

    -- For user convenience, DataAnalysis groups all properties according to their
    -- 'DisplayGroupName'.
    DisplayGroupName TEXT NOT NULL, 

    -- Display name [English]. Used in GUIs when displaying the property value.
    DisplayName TEXT NOT NULL,

    -- Maps integer values to descriptive strings for display purposes, e.g.,
    -- "0:Positive;1:Negative" for the property with permanent name 'Mode_IonPolarity'.
    DisplayValueText TEXT NOT NULL,

    -- A printf()-compatible format string that suggests the form in which to best display
    -- the value of this property in a GUI.
    DisplayFormat TEXT NOT NULL,

    -- The physical dimension of this quantity, e.g. "µs", "V", or "TOF cycles".
    DisplayDimension TEXT NOT NULL,

    -- Unused, yet
    Description TEXT NOT NULL
);

CREATE UNIQUE INDEX PropertyDefinitionsIndex ON PropertyDefinitions (PermanentName);

-- This table just contains ids for groups of properties. Such a table is required
-- because every frame references no or one group of properties
CREATE TABLE PropertyGroups (
    -- Unique ID for this property group.
    Id INTEGER PRIMARY KEY
) WITHOUT ROWID;

-- This table contains settings for "slowly-varying" properties that apply to an entire
-- group of spectra. Each of these properties may be overridden at any time using the
-- 'FrameProperties' table. The 'Properties' view facilitates data extraction for the
-- user.
CREATE TABLE GroupProperties (
    PropertyGroup INTEGER NOT NULL,
    Property INTEGER NOT NULL,
    Value NOT NULL,
    PRIMARY KEY (PropertyGroup, Property),
    FOREIGN KEY (PropertyGroup) REFERENCES PropertyGroups (Id),
    FOREIGN KEY (Property) REFERENCES PropertyDefinitions (Id)
) WITHOUT ROWID;

-- This table contains settings for "quickly-varying" properties that apply only to a
-- single frame and override a possible setting in 'GroupProperties'.
CREATE TABLE FrameProperties (
    Frame INTEGER NOT NULL,
    Property INTEGER NOT NULL,
    Value NOT NULL,
    PRIMARY KEY (Frame, Property),
    FOREIGN KEY (Frame) REFERENCES Frames (Id)
    FOREIGN KEY (Property) REFERENCES PropertyDefinitions (Id)
) WITHOUT ROWID;

-- Additional MS/MS meta-information.
CREATE TABLE FrameMsMsInfo (
    -- The frame to which this information applies. Should be an MS^2 frame, i.e., the
    -- corresponding Frames.MsMsType should be 2.
    Frame INTEGER PRIMARY KEY,

    -- Links to the corresponding "parent" MS^1 frame, in which this precursor was
    -- found. (Due to possible out-of-order / asynchronous scan-task execution in the
    -- engine, this is not necessarily the first MS^1 frame preceding this MS^2 frame in
    -- the analysis.)
    -- Parent is NULL for MRM scans.
    Parent INTEGER,

    -- The mass to which the quadrupole has been tuned for isolating this particular
    -- precursor. (in the m/z calibration state that was used during acquisition). May or
    -- may not coincide with one of the peaks in the parent frame.
    TriggerMass REAL NOT NULL,

    -- The total 3-dB width of the isolation window (in m/z units), the center of which is
    -- given by 'TriggerMass'.
    IsolationWidth REAL NOT NULL,

    -- The charge state of the precursor as estimated by the precursor selection code that
    -- controls the DDA acquisition. Can be NULL, which means that the charge state
    -- could not be determined, e.g., because only one isotope peak could be detected.
    PrecursorCharge INTEGER,

    -- Collision energy (in eV) using which this frame was produced.
    CollisionEnergy REAL NOT NULL,

    FOREIGN KEY (Frame) REFERENCES Frames (Id)
);

-- List of calibration information for this analysis.
-- 
-- There is a core set of "essential calibration information" that is always set in a valid TDF
-- analysis.
CREATE TABLE CalibrationInfo (
    KeyPolarity CHAR(1) CHECK (KeyPolarity IN ('+', '-')),
    KeyName TEXT,
    Value TEXT,
    PRIMARY KEY (KeyPolarity, KeyName)
    );

-- This table collects information about the various segments defined in
-- otofControl. Every frame is member of at most one segment.
CREATE TABLE Segments (
    -- Segment number. Normally numbered consecutively starting from one; however, "empty"
    -- segments for which no spectra have been acquired will not have a corresponding
    -- entry in this table.
    Id INTEGER PRIMARY KEY,

    -- Number of first frame included in this segment.
    FirstFrame INTEGER NOT NULL,

    -- Number of last frame included in this segment.
    LastFrame INTEGER NOT NULL,

    -- Has this been marked as a calibration segment by the user?
    IsCalibrationSegment BOOLEAN NOT NULL,

    FOREIGN KEY (FirstFrame) REFERENCES Frames (Id),
    FOREIGN KEY (LastFrame) REFERENCES Frames (Id)
);

--
-- This is the top-level view on the generic properties which most users would use.
--
-- Example:
--
--   SELECT Frame, Property, Value FROM Properties
--      WHERE Property IN (1,2,3) AND Frame IN (25145,1,23);
--
-- Possible output:
--
--   1|1|
--   1|2|872
--   1|3|
--   23|1|273647
--   23|2|1744
--   23|3|2
--   25145|1|273647
--   25145|2|1744
--   25145|3|2
--
-- Non-set variables are SQL "NULL", as in the first and third result rows above.
--
-- Don't rely on any implicit ordering of the result. Add an ORDER BY clause to the
-- SQL statement if required.
--
-- -------------------------------------------------------------------------------
--
-- Another example: a "property chromatogram", i.e., the dependence of a property's value
-- on time.
--
-- SELECT f.Id, f.Time, p.Value FROM Frames f
--  JOIN Properties p ON p.Frame=f.Id
--    AND p.Property=(SELECT Id FROM PropertyDefinitions WHERE PermanentName='TOF_DeviceTempCurrentValue1')
-- ORDER BY f.Time;
--
CREATE VIEW Properties AS
    SELECT s.Id Frame, pd.Id Property, COALESCE(fp.Value, gp.Value) Value
    FROM Frames s
    JOIN PropertyDefinitions pd
    LEFT JOIN GroupProperties gp ON gp.PropertyGroup=s.PropertyGroup AND gp.Property=pd.Id
    LEFT JOIN FrameProperties fp ON fp.Frame=s.Id AND fp.Property=pd.Id;

-- Logs warnings and errors occurring during acquisition.
CREATE TABLE ErrorLog (
    -- Id of frame to which this log entry belongs.
    Frame INTEGER NOT NULL,

    -- Number of scan inside the frame to which this log entry belongs.
    Scan INTEGER,

    -- Free-text message giving more information on what happened. Not intended for
    -- parsing by software.
    Message TEXT NOT NULL
);

-- Bruker-internal information on m/z calibration for individual frames/spectra. 
-- For digitizer delay / timebase information, refer to the global metadata.
CREATE TABLE MzCalibration (
-- Unique ID for this mz calibration.
    Id INTEGER PRIMARY KEY,
    ModelType INTEGER NOT NULL,
    DigitizerTimebase REAL NOT NULL,
    DigitizerDelay REAL NOT NULL,
    T1 REAL NOT NULL,
    T2 REAL NOT NULL,
    dC1 REAL NOT NULL,
    dC2 REAL NOT NULL,
    C0
-- C1,
-- etc. - as many columns as required to represent all calibrations according to
-- 'ModelType'.
);

INSERT INTO GlobalMetadata VALUES ('SchemaType', 'TDF');
INSERT INTO GlobalMetadata VALUES ('SchemaVersionMajor', '3');
INSERT INTO GlobalMetadata VALUES ('SchemaVersionMinor', '4');
-- An analysis consists of several frames (or "mobility scans"). A frame represents all
-- TOF scans acquired during a single TIMS-voltage ramp.
CREATE TABLE Frames (
    -- Unique ID for this frame. Numbered consecutively starting from one.
    Id INTEGER PRIMARY KEY,

    -- Time (in seconds), relative to the start time of the acquisition.
    Time REAL NOT NULL,

    -- Ionization mode used when acquiring this frame.
    Polarity CHAR(1) CHECK (Polarity IN ('+', '-')) NOT NULL,

    -- This enumaration type describes the mode in which MS and MS/MS frames where acquired:
    --
    --   0 = MS
    --   1 = AutoMSMS
    --   2 = MRM
    --   3 = in-source CID
    --   4 = broadband CID
    --   8 = PASEF
    --   9 = DIA
    --  10 = PRM
    --  20 = Maldi
    ScanMode INTEGER NOT NULL,

    -- Type of this frame:
    -- 
    --   0 = MS frame
    --   2 = MS/MS fragment frame
    --   8 = PASEF frame
    --   9 = DIA frame
    --  10 = PRM frame
    MsMsType INTEGER NOT NULL,

    -- ID for accessing mobility-resolved data for this frame. May be
    -- NULL in case this frame does not have any mobility-resolved data.
    TimsId INTEGER,

    -- Maximum intensity occurring in all data belonging to this frame (do not use to 
    -- generate a BPC directly from this!).
    MaxIntensity INTEGER NOT NULL,

    -- Sum of all intensities occurring in the data belonging to this frame (can quickly
    -- generate a TIC from this).
    SummedIntensities INTEGER NOT NULL,

    -- The number of TOF scans that contributed to this frame. If 0, this frame does not
    -- contain any data.
    NumScans INTEGER NOT NULL,

    -- The number of peaks stored in this frame (total over all scans).
    NumPeaks INTEGER NOT NULL,

    -- ID of the mz calibration for this frame. Every frame has exactly one mz calibration.h
    MzCalibration INTEGER NOT NULL,

    -- Measured Temperature1 of the device. Required to perform a temperature compensated mz calibration
    T1 REAL NOT NULL,

    -- Measured Temperature2 of the device. Required to perform a temperature compensated mz calibration
    T2 REAL NOT NULL,

    -- ID of the TIMS calibration for this frame. Every Frame has exactly one TIMS calibration.
    TimsCalibration INTEGER NOT NULL,

    -- The property group for this frame. May be overridden by properties for this
    -- specific frame in table 'FrameProperties'. Use the 'Properties' view for easy
    -- access to frame properties. This field may be NULL, in which case only the
    -- FrameProperties apply to this frame.
    PropertyGroup INTEGER,

    AccumulationTime REAL NOT NULL,

    RampTime REAL NOT NULL,

    FOREIGN KEY (MzCalibration) REFERENCES MzCalibration (Id),
    FOREIGN KEY (TimsCalibration) REFERENCES TimsCalibration (Id),
    FOREIGN KEY (PropertyGroup) REFERENCES PropertyGroups (Id)
);

CREATE UNIQUE INDEX FramesTimeIndex ON Frames (Time);

-- Bruker-internal information on TIMS calibration for individual frames.
CREATE TABLE TimsCalibration (
-- Unique ID for this TIMS calibration.
    Id INTEGER PRIMARY KEY,
    ModelType INTEGER NOT NULL,
    C0
-- C1,
-- etc. - as many columns as required to represent all calibrations according to
-- 'ModelType'.
);

-- Multiple DIA PASEF windows are grouped using an Id from this table.
-- A DIA PASEF window group represents a list of non-overlapping scan-number ranges.
-- For each of the scans inside a given range, the quadrupole has been configured
-- to be at the same specified isolation m/z and width.
CREATE TABLE DiaFrameMsMsWindowGroups (
"   Id INTEGER PRIMARY KEY"
);

-- This table describes DIA PASEF windows. A DIA PASEF window represents a scan-number range.
-- For each of the scans inside a given range, the quadrupole has been configured to be at
-- the same specified isolation m/z and width.
CREATE TABLE DiaFrameMsMsWindows (

    -- The DIA window group to which this information applies.
   WindowGroup INTEGER NOT NULL,

    -- Beginning of scan-number range where quadrupole was tuned to 'IsolationMz'. (R5)
   ScanNumBegin INTEGER NOT NULL,

    -- End (exclusive) of scan-number range where quadrupole was tuned to 'IsolationMz'.
   ScanNumEnd INTEGER NOT NULL,

    -- The isolation m/z (in the m/z calibration state that was used during
    -- acquisition). The quadrupole has been tuned to this mass during fragmentation
    -- between 'ScanNumBegin' and 'ScanNumEnd'.
   IsolationMz REAL NOT NULL,

    -- Specifies the total 3-dB width of the isolation window (in m/z units), the center
    -- of which is given by 'IsolationMz'.
   IsolationWidth REAL NOT NULL,

    -- Collision energy (in eV) set between 'ScanNumBegin' and 'ScanNumEnd'.
   CollisionEnergy REAL NOT NULL,

   PRIMARY KEY(WindowGroup, ScanNumBegin),

   FOREIGN KEY (WindowGroup) REFERENCES DiaFrameMsMsWindowGroups (Id)
) WITHOUT ROWID;

-- This table describes DIA PASEF frames. For every DIA PASEF frame a DIA PASEF window group is referenced,
-- that represents a list of non-overlapping scan-number ranges. For each of the scans inside a given range, the
-- quadrupole has been configured to be at the same specified isolation m/z and width.
--
-- NOTE: every frame for which there is an entry in this table, will have Frames.MsMsType = 9
--
-- NOTE: DIA PASEF only allows a single DIA PASEF window group per frame. The same window group will be used
-- repeatedly throughout a DIA PASEF measurment.
CREATE TABLE DiaFrameMsMsInfo (

    -- The DIA PASEF frame to which this information applies.
   Frame INTEGER PRIMARY KEY,

    -- The DIA window group to which this information applies.
   WindowGroup INTEGER NOT NULL,

   FOREIGN KEY (Frame) REFERENCES Frames (Id),
   FOREIGN KEY (WindowGroup) REFERENCES DiaFrameMsMsWindowGroups (Id)
);

-- Rationale behind this table is to include some metadata for every target that might be useful for
-- generic processing software (i.e., processing software which does not know anything about the
-- software that set up the method / target list in the first place).
CREATE TABLE PrmTargets (
 
    -- Number assigned by the acquisition software, uniquely identifies this target (only within
    -- this raw-data file).
    Id INTEGER PRIMARY KEY,
 
    -- Optionally, an ID assignable by external software setting up a PRM method. Not interpreted by
    -- Bruker acquisition software. (The UNIQUE constraint implicitly adds an index for this column,
    -- so efficient queries are possible.)
    ExternalId TEXT UNIQUE,

    -- The retention time (in seconds)
    Time REAL NOT NULL,
 
    -- Inverse reduced mobility (1/K0)
    OneOverK0 REAL NOT NULL,
 
    -- The monoisotopic m/z
    MonoisotopicMz REAL NOT NULL,
 
    -- The charge state
    Charge INTEGER NOT NULL,
 
    -- Optionally, a free-text description for display purposes (might contain, e.g., peptide
    -- sequence, sum formula, substance name, SMILES string, etc.)
    Description TEXT
);
 
-- This table describes PRM PASEF frames. For every PRM PASEF frame, it represents a list of
-- non-overlapping scan-number ranges. For each of the scans inside a given range, the
-- quadrupole has been configured to be at the same specified isolation m/z and width.
--
-- NOTE: every frame for which there is an entry in this table, will have Frames.MsMsType
-- = 10, and no entry in any other *FrameMsMsInfo table.
CREATE TABLE PrmFrameMsMsInfo (
    -- The PRM PASEF frame to which this information applies.
    Frame INTEGER NOT NULL,
 
    -- Beginning of scan-number range where quadrupole was tuned to 'IsolationMz'.
    ScanNumBegin INTEGER NOT NULL,
 
    -- End (exclusive) of scan-number range where quadrupole was tuned to 'IsolationMz'.
    ScanNumEnd INTEGER NOT NULL,
 
    -- The isolation m/z (in the m/z calibration state that was used during
    -- acquisition). The quadrupole has been tuned to this mass during fragmentation
    -- between 'ScanNumBegin' and 'ScanNumEnd'.
    IsolationMz REAL NOT NULL,
 
    -- Specifies the total 3-dB width of the isolation window (in m/z units), the center
    -- of which is given by 'IsolationMz'.
    IsolationWidth REAL NOT NULL,
 
    -- Collision energy (in eV) set between 'ScanNumBegin' and 'ScanNumEnd'.
    CollisionEnergy REAL NOT NULL,
 
    -- The target measured, reference to the PrmTargets table
    Target INTEGER NOT NULL,
 
    PRIMARY KEY (Frame, ScanNumBegin),
 
    FOREIGN KEY (Frame) REFERENCES Frames(Id),
    FOREIGN KEY (Target) REFERENCES PrmTargets(Id)
) WITHOUT ROWID;
 
 
-- this is to enable efficient queries like "give me all frames where the given target has been measured":
CREATE INDEX PrmFrameMsMsInfoTargetIndex ON PrmFrameMsMsInfo (Target);


-- Table containing all precursors / features detected by an on-line precursor-selection algorithm.
--
-- 'Precursors' table lists only the information derived from
-- MS^1 spectra. How a precursor is actually measured (= quadrupole settings
-- and scan-number begin/end) remains in the 'PasefFrameMsMsInfo' table.
-- 
CREATE TABLE Precursors (
    -- Number that uniquely identifies this precursor in this analysis.
   Id INTEGER PRIMARY KEY,

    -- m/z of the largest (most intensive) peak in this precursor's isotope pattern.
   LargestPeakMz REAL NOT NULL,

    -- Intensity-weighted average m/z of this precursor's isotope pattern. If only one
    -- peak was detected, this will be the m/z of that peak, and identical to
    -- 'LargestPeakMz'.
   AverageMz REAL NOT NULL,

    -- An estimate for the monoisotopic m/z derived from the isotope pattern of the
    -- precursor. May be NULL when detection failed.
   MonoisotopicMz REAL,

    -- The charge state of the precursor (a positive integer) as estimated from the
    -- isotope pattern of the precursor. May be NULL when detection failed.
   Charge INTEGER,

    -- Mobility (in scan-number units) of this precursor in the corresponding MS^1 frame.
   ScanNumber REAL NOT NULL,

    -- Intensity of this precursor in the corresponding MS^1 frame.
   Intensity REAL NOT NULL,

    -- The corresponding MS^1 frame in which this precursor was detected. In the case that
    -- MS^1 frames were repeatedly measured and averaged to improve SNR for precursor
    -- detection, the TDF stores those frames individually, and this field points to the
    -- last of that set of frames. (Field can be NULL, which means that the parent MS^1 is
    -- not included in the TDF; e.g., because recording started in the middle of the DDA
    -- cycle, or due to an error writing the MS^1 data.)
   Parent INTEGER,

   FOREIGN KEY(Parent) REFERENCES Frames(Id)
);
            
CREATE INDEX IF NOT EXISTS PrecursorsParentIndex ON Precursors (Parent);


-- This table describes PASEF frames. For every PASEF frame, it represents a list of
-- non-overlapping scan-number ranges. For each of the scans inside a given range, the
-- quadrupole has been configured to be at the same specified isolation m/z and width.
--
-- NOTE: every frame for which there is an entry in this table, will have Frames.MsMsType
-- = 8, and no entry in any other *FrameMsMsInfo table.
--
-- NOTE: PASEF acquisition allows, in principle, to change the quadrupole isolation window
-- on a per-scan basis. Therefore, the compound primary key (Frame, ScanNumBegin) uniquely
-- identifies each fragmentation region.
CREATE TABLE PasefFrameMsMsInfo (

    -- The PASEF frame to which this information applies.
   Frame INTEGER NOT NULL,

    -- Beginning of scan-number range where quadrupole was tuned to 'IsolationMz'. (R5)
   ScanNumBegin INTEGER NOT NULL,

    -- End (exclusive) of scan-number range where quadrupole was tuned to 'IsolationMz'.
   ScanNumEnd INTEGER NOT NULL,

    -- The isolation m/z (in the m/z calibration state that was used during
    -- acquisition). The quadrupole has been tuned to this mass during fragmentation
    -- between 'ScanNumBegin' and 'ScanNumEnd'.
   IsolationMz REAL NOT NULL,

    -- Specifies the total 3-dB width of the isolation window (in m/z units), the center
    -- of which is given by 'IsolationMz'.
   IsolationWidth REAL NOT NULL,

    -- Collision energy (in eV) set between 'ScanNumBegin' and 'ScanNumEnd'.
   CollisionEnergy REAL NOT NULL,

    -- Optionally, the ID of a precursor that was measured in this frame and scan-number
    -- range. May be NULL in case this measurement was not based on a precursor search,
    -- but manually programmed instead.
   Precursor INTEGER,

   PRIMARY KEY(Frame, ScanNumBegin),

   FOREIGN KEY(Frame) REFERENCES Frames(Id),
   FOREIGN KEY(Precursor) REFERENCES Precursors(Id)
) WITHOUT ROWID;
            
CREATE INDEX PasefFrameMsMsInfoPrecursorIndex ON PasefFrameMsMsInfo (Precursor);


