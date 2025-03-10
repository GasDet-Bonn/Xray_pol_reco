import numpy as np
from tqdm import tqdm
import multiprocessing
import sys
import os
import h5py
import argparse

def copy_hdf5_file(source_file, destination_file, nodes_to_copy=None):
    if nodes_to_copy is None:
        nodes_to_copy = []

    # Open the source file
    with h5py.File(source_file, 'r') as src:
        # Create the destination file or overwrite it
        with h5py.File(destination_file, 'w') as dst:
            # Copy the data
            def recursive_copy(src_group, dst_group, first_level=False):
                for name, obj in src_group.items():
                    if first_level and name not in nodes_to_copy:
                        continue
                    if isinstance(obj, h5py.Group):
                        dst_group.create_group(name)
                        recursive_copy(obj, dst_group[name])
                    elif isinstance(obj, h5py.Dataset):
                        dst_group.create_dataset(name, data=obj[:])
                    # Copy the attributes
                    for attr_name, attr_value in obj.attrs.items():
                        dst_group[name].attrs[attr_name] = attr_value

            # Copy the attributes of the root group
            for attr_name, attr_value in src.attrs.items():
                dst.attrs[attr_name] = attr_value

            # Start copying from the root group with first_level flag
            recursive_copy(src, dst, first_level=True)

def reco(coords, charges, radius_min = None, radius_max = None, phi_1 = None):
    # Calculate the center of charge for 2D or 3D and shift the coordinates by the center of charge
    if coords.shape[0] == 2:
        x_pos, y_pos = coords
        center = np.average(coords, axis=1, weights=charges)
        X = np.vstack((x_pos - center[0], y_pos - center[1]))
    elif coords.shape[0] == 3:
        x_pos, y_pos, z_pos = coords
        center = np.average(coords, axis=1, weights=charges)
        X = np.vstack((x_pos - center[0], y_pos - center[1], z_pos - center[2]))
    else:
        raise ValueError("Coordinates must be 2D or 3D")

    # Create the covariance matrix for the charge distribution
    M = np.dot(X * charges, X.T)

    # Get the eigenvalues and eigenvectors for the covariance matrix
    eigenvalues, eigenvectors = np.linalg.eig(M)

    # The axis that maximises the second moment corresponds to the eigenvector with the biggest eigenvalue
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]

    # Project the new axis on the xy plane and calculate its angle to the x axis
    projection_xy = np.array([principal_axis[0], principal_axis[1]])
    angle = np.arctan2(projection_xy[1], projection_xy[0])

    # Project the data points on the projection of the axis in the xy plane
    projection_xy_fit = np.dot(X[:2].T, principal_axis[:2])

    # Calculate the skewness in the xy plane
    skewness_xy = np.sum(charges * projection_xy_fit**3) / np.sum(charges)

    # For the 3D case project the data points also on the 3D axis and calculate the skewness for this
    if coords.shape[0] == 3:
        projection_new_plane = np.dot(X.T, principal_axis)
        skewness_new_plane = np.sum(charges * projection_new_plane**3) / np.sum(charges)
    else:
        skewness_new_plane = 0
        angle_new_plane = 0

    # For the 3D case calculate the angle of the axis in the xz-plane (currently not used)
    if coords.shape[0] == 3:
        projection_xz = np.array([principal_axis[0], principal_axis[2]])
        angle_new_plane = np.arctan2(projection_xz[1], projection_xz[0])
    else:
        angle_new_plane = 0

    # For the 2D and the 3D case get list of indices in the "left" and the "right" part of the track
    left_indices = np.where(projection_xy_fit < 0)[0]
    right_indices = np.where(projection_xy_fit >= 0)[0]
    if coords.shape[0] == 3:
        upper_indices = np.where(projection_new_plane < 0)[0]
        lower_indices = np.where(projection_new_plane >= 0)[0]

    # Based on symmetry extend the angle to -pi to +pi
    if phi_1 == None:
        if skewness_xy > 0:
            if angle > 0:
                angle = -np.pi + angle
            else:
                angle = np.pi + angle
        else:
            angle = angle
    else:
        if (angle > (np.pi/2) and phi_1 < (-np.pi/2)) or (angle < (-np.pi/2) and phi_1 > (np.pi/2)):
            if (np.fabs(phi_1-angle)<(3*np.pi/2)):
                if angle > 0:
                    angle = angle - np.pi
                else:
                    angle = np.pi + angle
        elif (np.fabs(phi_1-angle)>(np.pi/2)):
            if angle > 0:
                angle = angle - np.pi
            else:
                angle = np.pi + angle

        return angle, angle_new_plane, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    # Decide based on the third moment which part to keep. skewness_xy is the 3rd moment in the xy plane, skewness_new_plane for the full 3D information
    if skewness_xy <= 0:
        if coords.shape[0] == 3:
            if skewness_new_plane <= 0:
                start_indices = np.intersect1d(left_indices, upper_indices)
                end_indices = np.setdiff1d(np.arange(len(x_pos)), start_indices)
            else:
                start_indices = np.intersect1d(left_indices, lower_indices)
                end_indices = np.setdiff1d(np.arange(len(x_pos)), start_indices)
        else:
            start_indices = left_indices
            end_indices = right_indices
    else:
        if coords.shape[0] == 3:
            if skewness_new_plane <= 0:
                start_indices = np.intersect1d(right_indices, upper_indices)
                end_indices = np.setdiff1d(np.arange(len(x_pos)), start_indices)
            else:
                start_indices = np.intersect1d(right_indices, lower_indices)
                end_indices = np.setdiff1d(np.arange(len(x_pos)), start_indices)
        else:
            start_indices = right_indices
            end_indices = left_indices

    if coords.shape[0] == 3:
        m2_max = np.sum(charges * np.power(projection_new_plane, 2)) / np.sum(charges)
        m2_min = np.sum(charges * np.power(np.dot(X.T, eigenvectors[:, np.argmin(eigenvalues)]), 2)) / np.sum(charges)
        quality = m2_max / m2_min
    else:
        m2_max = np.sum(charges * np.power(projection_xy_fit, 2)) / np.sum(charges)
        projection_min = eigenvectors[:, np.argmin(eigenvalues)]
        m2_min = np.sum(charges * np.power(np.dot(X[:2].T, projection_min[:2]), 2)) / np.sum(charges)
        quality = m2_max / m2_min

    # Additional selection of start pixel based on a radius around the center of charge.
    # The radius is normalised to the second moment to reduce the energy dependence of it.
    if radius_min is not None and radius_max is not None:
        if coords.shape[0] == 3:
            radii = np.power(x_pos - center[0], 2) + np.power(y_pos-center[1], 2) + np.power(z_pos-center[2], 2)
        else:
            radii = np.power(x_pos - center[0], 2) + np.power(y_pos-center[1], 2)
        circle_indices_inner = np.where(radii > np.power(radius_min * np.sqrt(m2_max), 2))[0]
        circle_indices_outer = np.where(radii < np.power(radius_max * np.sqrt(m2_max), 2))[0]
        circle_indices = np.intersect1d(circle_indices_inner, circle_indices_outer)
        absorption_indices = np.array(np.intersect1d(start_indices, circle_indices), dtype=int)
        absorption_point_x = np.average(x_pos[absorption_indices], weights=charges[absorption_indices])
        absorption_point_y = np.average(y_pos[absorption_indices], weights=charges[absorption_indices])
        if coords.shape[0] == 3:
            absorption_point_z = np.average(z_pos[absorption_indices], weights=charges[absorption_indices])
    else:
        absorption_point_x = np.nan
        absorption_point_y = np.nan
        absorption_point_z = np.nan

    if coords.shape[0] == 3:
        return angle, angle_new_plane, start_indices, end_indices, skewness_xy, center[0], center[1], center[2], skewness_new_plane, absorption_point_x, absorption_point_y, absorption_point_z, quality
    else:
        return angle, angle_new_plane, start_indices, end_indices, skewness_xy, center[0], center[1], 0, 0, absorption_point_x, absorption_point_y, 0, quality

# Two step reconstruction for Timepix data
def reco_tpx(x, y, charges_raw, matrix, radius_min, radius_max, weighting):
    try:
        if matrix:
            charges = charges_raw * matrix[np.array(y, dtype=int), np.array(x, dtype=int)]
        else:
            charges = charges_raw
        phi_1, _, start_indices, end_indices, m3, xc, yc, _, _, absorption_point_x, absorption_point_y, _, quality = reco(np.array([x, y]), charges, radius_min=radius_min, radius_max=radius_max)
        if weighting:
            # Perform the second step with all pixels but weight the charges by the distance to the absorption point
            distance = np.sqrt(np.power(x-absorption_point_x, 2) + np.power(y-absorption_point_y, 2))
            weighted_charges = charges * np.exp(-(distance/weighting))
            phi_2, _, _, _, _, xc2, yc2, _, _, _, _, _, _ = reco(np.array([x, y]), weighted_charges, phi_1=phi_1)
        else:
            # Perform the second step with the start pixels
            phi_2, _, _, _, _, xc2, yc2, _, _, _, _, _, _ = reco(np.array([x[start_indices], y[start_indices]]), charges[start_indices], phi_1=phi_1)
        if radius_min is None and radius_max is None:
            absorption_point_x = xc2
            absorption_point_y = yc2
        if phi_2 == 0 or phi_2 == np.pi or phi_2 == -np.pi or phi_2 == np.pi/2 or phi_2 == -np.pi/2 or phi_2 == np.pi/4 or phi_2 == -np.pi/4 or phi_2 == 3*np.pi/4 or phi_2 == -3*np.pi/4: 
            return np.nan, np.nan, np.empty(0), np.empty(0), np.nan, np.nan, np.nan
        else:
            return phi_1, phi_2, start_indices, end_indices, absorption_point_x, absorption_point_y, quality
        return np.nan, np.nan, np.empty(0), np.empty(0), np.nan, np.nan, np.nan

# Two step reconstruction for Timepix3 data
def reco_tpx3(x, y, toa, ftoa, charges_raw, full3d, velocity, matrix, radius_min, radius_max, weighting):
    try:
        if matrix is not None:
            charges = charges_raw * matrix[np.array(y, dtype=int), np.array(x, dtype=int)]
        else:
            charges = charges_raw
        z = ((toa * 25 + 1) - ftoa * 1.5625)*velocity/55.
        phi_1, theta_1, start_indices, end_indices, m3, xc, yc, zc, m3_z, absorption_point_x, absorption_point_y, absorption_point_z, quality = reco(np.array([x, y, z]), charges, radius_min=radius_min, radius_max=radius_max)
        if weighting:
            # Perform the second step with all pixels but weight the charges by the distance to the absorption poisnt
            if full3d:
                distance = np.sqrt(np.power(x-absorption_point_x, 2) + np.power(y-absorption_point_y, 2), + np.power(z-absorption_point_z, 2))
                weighted_charges = charges * np.exp(-(distance/weighting))
                phi_2, theta_2, _, _, _, xc2, yc2, _, _, _, _, _, _ = reco(np.array([x, y, z]), weighted_charges, phi_1=phi_1)
            else:
                distance = np.sqrt(np.power(x-absorption_point_x, 2) + np.power(y-absorption_point_y, 2))
                weighted_charges = charges * np.exp(-(distance/weighting))
                phi_2, theta_2, _, _, _, xc2, yc2, _, _, _, _, _, _ = reco(np.array([x, y]), weighted_charges, phi_1=phi_1)
        else:
            # Perform the second step with the start pixels
            if full3d:
                phi_2, theta_2, _, _, _, xc2, yc2, _, _, _, _, _, _ = reco(np.array([x[start_indices], y[start_indices], z[start_indices]]), charges[start_indices], phi_1=phi_1)
            else:
                phi_2, theta_2, _, _, _, xc2, yc2, _, _, _, _, _, _ = reco(np.array([x[start_indices], y[start_indices]]), charges[start_indices], phi_1=phi_1)
        if radius_min is None and radius_max is None:
            absorption_point_x = xc2
            absorption_point_y = yc2
        if phi_2 == 0 or phi_2 == np.pi or phi_2 == -np.pi or phi_2 == np.pi/2 or phi_2 == -np.pi/2 or phi_2 == np.pi/4 or phi_2 == -np.pi/4 or phi_2 == 3*np.pi/4 or phi_2 == -3*np.pi/4: 
            return np.nan, np.nan, np.nan, np.nan, np.empty(0), np.empty(0), np.nan, np.nan, np.nan
        else:
            return phi_1, phi_2, theta_1, theta_2, start_indices, end_indices, absorption_point_x, absorption_point_y, quality
        return np.nan, np.nan, np.nan, np.nan, np.empty(0), np.empty(0), np.nan, np.nan, np.nan

def tpx_wrapper(inputs):
    return reco_tpx(*inputs)

def tpx3_wrapper(inputs):
    return reco_tpx3(*inputs)

def write_dataset(h5file, path, data, overwrite=False, dtype=None):
    # Check if dataset already exists
    if path in h5file:
        if overwrite:
            # Delete dataset if overwrite option is true
            del h5file[path]
        else:
            # Raise an error is dataset exists but overwrite is false
            raise ValueError(f"Error: Dataset {path} exists and overwrite is set to False.")
    
    # Write data
    if dtype:
        h5file.create_dataset(path, (len(data),), dtype=dtype)
        h5file[path][...] = data
    else:
        h5file.create_dataset(path, data=data)

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Angular reconstruction of Timepix and Timepix3 polarimetry data')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file')
    parser.add_argument('--velocity', type=float, help='Drift velocity of electrons in the drift field in µm/ns. Default is 1 µm/ns. Only influences 3D reconstruction.', default=1)
    parser.add_argument('--rotation', type=float, help='Rotation of the input coordinates in the xy plane with respect to the x axis. Given in degrees.', default=0)
    parser.add_argument('--shiftcenter', type=float, nargs=2, help='Shift each hit by x/y pixels.', default=None)
    parser.add_argument('--full2d', action='store_true', help='Analyze Timepix3 data only in 2D')
    parser.add_argument('--full3d', action='store_true', help='Analyze Timepix3 data in the full 3D approach instead of 3D for the first step and 2D for the second step')
    parser.add_argument('--radius', type=float, nargs=2, help='Set the lower and upper normalised radius for the pixel selection around the center of charge in the first step.', default=None)
    parser.add_argument('--weighting', type=float, help='Instead of cutting pixels for the second step weight them to by the distance to the absorption point. Requires radius argument.')
    parser.add_argument('--overwrite', action='store_true', help='If the hdf5 file already contains reconstructed angular data, this option activates overwriting it.')
    parser.add_argument('--weights', type=str, help='Path to a txt file that contains a 256 by 256 matrix with weights for all pixels.')
    parser.add_argument('--iweights', type=str, help='Path to a txt file that contains a 256 by 256 matrix with weights for all pixels. The weights are inverted')
    parser.add_argument('--output', type=str, help='Instead of storing the data in the input datafile the data is stored in a given output file.')
    args = parser.parse_args()

    if args.full2d and args.full3d:
        parser.error('You can not use the options --full2d and --full3d at the same time.')

    # Process the arguments
    run = args.runpath
    if run.endswith('h5'):
        datafile = os.path.basename(run)
    else:
        print("Please choose a correct data file")
        return
    angle_offset = args.rotation
    tpx3_2d = args.full2d
    tpx3_full3d = args.full3d
    velocity = args.velocity

    # Open the corresponding datafile
    f = h5py.File(run, 'r+')
    filename = datafile.replace('.h5', '')

    if args.weights and args.iweights:
        print("You can not use 'weights' and 'iweights' at the same time")
        return
    if args.weights:
        path = args.weights
        matrix = np.loadtxt(path, delimiter='\t')
        matrix = matrix / np.nanmedian(matrix)
        matrix[np.where(matrix == np.inf)] = 0
        matrix = np.nan_to_num(matrix, nan=0.0)
        matrix = matrix
    elif args.iweights:
        path = args.iweights
        matrix = np.loadtxt(path, delimiter='\t')
        matrix = np.nanmean(matrix) / matrix
        matrix[np.where(matrix == np.inf)] = 0
        matrix = np.nan_to_num(matrix, nan=0.0)
        matrix = matrix
    else:
        matrix = None

    if args.weighting and not args.radius:
        print("When the argument 'weighting' is used, also 'radius' must be used.")
        return

    if args.radius:
        radius_min, radius_max = args.radius
    else:
        radius_min = None
        radius_max = None

    if args.shiftcenter:
        x_shift, y_shift = args.shiftcenter
    else:
        x_shift = None
        y_shift = None

    timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
    reconstruction = f['reconstruction']
    for name in reconstruction:
        # Load the x and y coordinates of the pixel with the option to rotate
        posx_raw = f.get('reconstruction/' + name + '/chip_0/x')[:]
        posy_raw = f.get('reconstruction/' + name + '/chip_0/y')[:]
        rot = np.radians(angle_offset)
        posx = (posx_raw * np.cos(rot) - posy_raw * np.sin(rot))
        posy = (posx_raw * np.sin(rot) + posy_raw * np.cos(rot))
        if args.shiftcenter:
            posx = posx - x_shift
            posy = posy - y_shift
            filtered_x_coords = []
            filtered_y_coords = []
            print("Shift coordinates")
            for x_arr, y_arr in tqdm(zip(posx, posy), total=len(posx)):
                mask = (x_arr >= 0) & (x_arr <= 255) & (y_arr >= 0) & (y_arr <= 255)
                filtered_x_coords.append(x_arr[mask])
                filtered_y_coords.append(y_arr[mask])

            posx = np.array(filtered_x_coords, dtype=object)
            posy = np.array(filtered_y_coords, dtype=object)

        # Fot Timepix3 use additionally the timestamp per pixel
        if timepix_version == 'Timepix3':
            toa = f.get('reconstruction/' + name + '/chip_0/ToA')[:]
            ftoa = f.get('reconstruction/' + name + '/chip_0/fToA')[:]

        # If there is calibrated charge data use it, otherwise use the ToT
        try:
            charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
        except:
            charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]

        # Get the number of cpu cores for the multiprocessing
        num_threads = multiprocessing.cpu_count() - 1

        print(timepix_version)
        if timepix_version == 'Timepix1' or tpx3_2d:
            inputs = list(zip(posx, posy, charge, [matrix]*len(posx), [radius_min]*len(posx), [radius_max]*len(posx), [args.weighting]*len(posx)))
            # Perform the reconstruction per event in multithreading
            with multiprocessing.Pool(processes=num_threads) as pool:
                results = list(tqdm(pool.imap(tpx_wrapper, inputs), total=len(inputs)))
            pool.close()
            pool.join()

            # Split up the results in separate arrays
            results = np.array(results, dtype=object)
            phi1 = np.array(results[:, 0], dtype=np.float64)
            phi2 = np.array(results[:, 1], dtype=np.float64)
            start = results[:, 2]
            end = results[:, 3]
            absorption_point_x = np.array(results[:, 4], dtype=np.float64)
            absorption_point_y = np.array(results[:, 5], dtype=np.float64)
            quality = np.array(results[:, 6], dtype=np.float64)

        elif timepix_version == 'Timepix3':
            inputs = list(zip(posx, posy, toa, ftoa, charge, [tpx3_full3d]*len(posx), [velocity]*len(posx), [matrix]*len(posx), [radius_min]*len(posx), [radius_max]*len(posx), [args.weighting]*len(posx)))
            # Perform the reconstruction per event in multithreading
            with multiprocessing.Pool(processes=num_threads) as pool:
                results = list(tqdm(pool.imap(tpx3_wrapper, inputs), total=len(inputs)))
            pool.close()
            pool.join()

            # Split up the results in separate arrays
            results = np.array(results, dtype=object)
            phi1 = np.array(results[:, 0], dtype=np.float64)
            phi2 = np.array(results[:, 1], dtype=np.float64)
            theta1 = np.array(results[:, 2], dtype=np.float64)
            theta2 = np.array(results[:, 3], dtype=np.float64)
            start = results[:, 4]
            end = results[:, 5]
            absorption_point_x = np.array(results[:, 6], dtype=np.float64)
            absorption_point_y = np.array(results[:, 7], dtype=np.float64)
            quality = np.array(results[:, 8], dtype=np.float64)

        dataset_path_firststage = f'reconstruction/{name}/chip_0/angle_firststage'
        dataset_path_secondstage = f'reconstruction/{name}/chip_0/angle_secondstage'
        dataset_path_indices_start = f'reconstruction/{name}/chip_0/start_indices'
        dataset_path_indices_end = f'reconstruction/{name}/chip_0/end_indices'
        dataset_path_absorption_point_x = f'reconstruction/{name}/chip_0/absorption_point_x'
        dataset_path_absorption_point_y = f'reconstruction/{name}/chip_0/absorption_point_y'
        dataset_path_quality= f'reconstruction/{name}/chip_0/quality'
        dt = h5py.special_dtype(vlen=np.dtype('float64'))

        if args.output:
            f.close()
            nodes_to_copy = ['reconstruction']
            copy_hdf5_file(args.runpath, args.output, nodes_to_copy)
            f = h5py.File(args.output, 'r+')

        # Save data to the hdf5 file
        try:
            write_dataset(f, dataset_path_firststage, phi1, overwrite=args.overwrite)
            write_dataset(f, dataset_path_secondstage, phi2, overwrite=args.overwrite)
            write_dataset(f, dataset_path_absorption_point_x, absorption_point_x, overwrite=args.overwrite)
            write_dataset(f, dataset_path_absorption_point_y, absorption_point_y, overwrite=args.overwrite)
            write_dataset(f, dataset_path_quality, quality, overwrite=args.overwrite)
            write_dataset(f, dataset_path_indices_start, start, overwrite=args.overwrite, dtype=dt)
            write_dataset(f, dataset_path_indices_end, end, overwrite=args.overwrite, dtype=dt)
        except ValueError as e:
            print(e)

        if timepix_version == 'Timepix3' and not tpx3_2d:
            dataset_path_firststage_theta = f'reconstruction/{name}/chip_0/theta_firststage'
            dataset_path_secondstage_theta = f'reconstruction/{name}/chip_0/theta_secondstage'
            try:
                write_dataset(f, dataset_path_firststage_theta, theta1, overwrite=args.overwrite)
                write_dataset(f, dataset_path_secondstage_theta, theta2, overwrite=args.overwrite)
            except ValueError as e:
                print(e)

if __name__ == "__main__":
    main()
