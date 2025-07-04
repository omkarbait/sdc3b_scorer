from astropy.io import fits
import numpy as np

def get_volume_from_coords(fits_file, x, y, z, cube_widths=(0.1, 0.1, 0.1), team_name='Team_SKAO'):
    """
    Integrates pixel values within a cube centered at (x, y, z) with specified physical side widths per axis.

    Parameters:
    - fits_file: Path to the FITS cube
    - x, y, z: Coordinates in physical units (same units as CRVAL/CDELT)
    - cube_widths: Tuple of cube widths (width_x, width_y, width_z) in physical units
    - team_name: Team name for axis-order handling

    Returns:
    - Sum of pixel values within the cube region
    """
    width_x, width_y, width_z = cube_widths

    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header


        if team_name == 'ReionYuga':

            data = data/np.sum(data)
        else:
            None

        def to_index(coord, crval, crpix, cdelt):
            return int(round((coord - crval) / cdelt + crpix - 1))

        # Convert physical coords to voxel indices
        ix = to_index(x, header['CRVAL1'], header['CRPIX1'], header['CDELT1'])
        iy = to_index(y, header['CRVAL2'], header['CRPIX2'], header['CDELT2'])
        iz = to_index(z, header['CRVAL3'], header['CRPIX3'], header['CDELT3'])

        # Half-widths in voxel units
        hx = int(round(width_x / abs(header['CDELT1']) / 2))
        hy = int(round(width_y / abs(header['CDELT2']) / 2))
        hz = int(round(width_z / abs(header['CDELT3']) / 2))

        shape = data.shape

        if team_name in ['LoreliB', 'EoR-PIE-MC', 'EoR-PIE']:
            # Axis order: (z, y, x)
            x_min, x_max = np.clip(ix - hx, 0, shape[2]-1), np.clip(ix + hx + 1, 0, shape[2])
            y_min, y_max = np.clip(iy - hy, 0, shape[1]-1), np.clip(iy + hy + 1, 0, shape[1])
            z_min, z_max = np.clip(iz - hz, 0, shape[0]-1), np.clip(iz + hz + 1, 0, shape[0])
            cube = data[z_min:z_max, y_min:y_max, x_min:x_max]
        else:
            # Axis order: (x, y, z)
            x_min, x_max = np.clip(ix - hx, 0, shape[0]-1), np.clip(ix + hx + 1, 0, shape[0])
            y_min, y_max = np.clip(iy - hy, 0, shape[1]-1), np.clip(iy + hy + 1, 0, shape[1])
            z_min, z_max = np.clip(iz - hz, 0, shape[2]-1), np.clip(iz + hz + 1, 0, shape[2])
            cube = data[x_min:x_max, y_min:y_max, z_min:z_max]

        

        return np.sum(cube)

