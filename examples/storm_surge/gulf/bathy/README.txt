README for Gulf Storm Surge simulation

gulf_carribean.tt3 - Bathymetry covering larger region including Carribean and
    part of the Atlantic sea-board.
    Lower left corner = (99W,8N), upper right corner = (50W,17N)

gulf_coarse_bathy.tt3 - Bathymetry covering entire gulf at 2 minute resolution
    Lower left corner = (99W,17N), upper right corner = (80W,17N)

houston_ship_channel.xyz - Bathymetry of the Houston ship channel from:

    USACE surveys of dredged shipping channels from the Galveston District.
    Coordinate system and datum have been converted to Geographic NAD'83 and NAVD'88

    Data has been located by Hugh Roberts

old_NOAA_Galveston_Houston.tt3 - Bathymetry covering Galveston and Houston ship 
    channel area with 3 second resolution.
    
    Lower left corner  (95º 26' W, 29º 06' N)
    Upper right corner (94º 25' W, 29º 55' N)

NOAA_Galveston_Houston.tt3 - New larger bathymetry coveraged for Galveston and 
    Houston area with 3 second resolution. (from NOAA Design-a-grid)
                30º 12' N = 30.2
              ______________        
    95º 52' W |            | 93º 24' W = 93.4
   = 95.8666  |____________|
                28º 38' N = 28.63333

extract_bathy.py - Python script for extracting requested subsections of xyz
    bathymetry file and turn them into gridded topography type 3.

