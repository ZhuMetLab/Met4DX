-- SQLite statements for setting up a "New TIMS spectrum format" database.
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

INSERT INTO GlobalMetadata VALUES ('SchemaType', 'TSF');
INSERT INTO GlobalMetadata VALUES ('SchemaVersionMajor', '3');
INSERT INTO GlobalMetadata VALUES ('SchemaVersionMinor', '2');
-- An analysis consists of several spectra (called frames in the db schema to be able to use the same queries
-- for tdf and tsf data.
CREATE TABLE Frames (
    -- Unique ID for this spectrum. Numbered consecutively starting from one.
    Id INTEGER PRIMARY KEY,

    -- Time (in seconds), relative to the start time of the acquisition.
    Time REAL NOT NULL,

    -- Ionization mode used when acquiring this spectrum.
    Polarity CHAR(1) CHECK (Polarity IN ('+', '-')) NOT NULL,

    -- This enumaration type describes the mode in which MS and MS/MS spectra where acquired:
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

    -- Type of this spectrum:
    -- 
    --   0 = MS frame
    --   2 = MS/MS fragment frame
    --   8 = PASEF frame
    --   9 = DIA frame
    --  10 = PRM frame
    MsMsType INTEGER NOT NULL,

    -- ID for accessing mz and intensity data for this spectrum. May be
    -- NULL in case this spectrum does not have any data.
    TimsId INTEGER,

    -- Maximum intensity occurring in all data belonging to this spectrum (can quickly
    -- generate a BPC from this).
    MaxIntensity INTEGER NOT NULL,

    -- Sum of all intensities occurring in the data belonging to this spectrum (can quickly
    -- generate a TIC from this).
    SummedIntensities INTEGER NOT NULL,

    -- The number of peaks stored in this line spectrum, NULL if no line spectrum is stored.
    NumPeaks INTEGER,

    -- ID of the mz calibration for this spectrum. Every spectrum has exactly one mz calibration.h
    MzCalibration INTEGER NOT NULL,

    -- Measured Temperature1 of the device. Required to perform a temperature compensated mz calibration
    T1 REAL NOT NULL,

    -- Measured Temperature2 of the device. Required to perform a temperature compensated mz calibration
    T2 REAL NOT NULL,

    -- The property group for this spectrum. May be overridden by properties for this
    -- specific spectrum in table 'FrameProperties'. Use the 'Properties' view for easy
    -- access to spectrum properties. This field may be NULL, in which case only the
    -- FrameProperties apply to this spectrum.
    PropertyGroup INTEGER,

    FOREIGN KEY (MzCalibration) REFERENCES MzCalibration (Id),
    FOREIGN KEY (PropertyGroup) REFERENCES PropertyGroups (Id)
);

CREATE UNIQUE INDEX FramesTimeIndex ON Frames (Time);

