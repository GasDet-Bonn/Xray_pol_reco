import numpy as np
from tqdm import tqdm
import multiprocessing
import sys
import os
import h5py
import argparse

def reco(coords, charges, phi_1 = None):
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

    # The axis that minimises the second moment corresponds to the eigenvector with the biggest eigenvalue
    principal_axis = eigenvectors[:, np.argmax(eigenvalues)]

    # Project the new axis on the xy plane and calculate its angle to the x axis
    projection_xy = np.array([principal_axis[0], principal_axis[1]])
    angle = np.arctan2(projection_xy[1], projection_xy[0])

    # Project the data points on the projection of the axis in the xy plane
    projection_xy_fit = np.dot(np.vstack((x_pos, y_pos)).T - center[:2], principal_axis[:2])

    # Calculate the skewness in the xy plane
    skewness_xy = np.sum(charges * (projection_xy_fit - np.average(projection_xy_fit, weights=charges))**3) / np.sum(charges)

    # For the 3D case project the data points also on the 3D axis and calculate the skewness for this
    if coords.shape[0] == 3:
        projection_new_plane = np.dot(coords.T - center, principal_axis)
        skewness_new_plane = np.sum(charges * (projection_new_plane - np.average(projection_new_plane, weights=charges))**3) / np.sum(charges)
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

    if coords.shape[0] == 3:
        return angle, angle_new_plane, left_indices, right_indices, skewness_xy, center[0], center[1], center[2], skewness_new_plane, upper_indices, lower_indices
    else:
        return angle, angle_new_plane, left_indices, right_indices, skewness_xy, center[0], center[1], 0, 0, 0, 0

# Two step reconstruction for Timepix data
def reco_tpx(x, y, charges):
    phi_1, _, d_i_left_indices, d_i_right_indices, m3, xc, yc, _, _, _, _ = reco(np.array([x, y]), charges)
    try:
        # Decide based on the third moment which part to keep
        if m3 <= 0:
            phi_2, _, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, _, _, _, _ = reco(np.array([x[d_i_left_indices], y[d_i_left_indices]]), charges[d_i_left_indices], phi_1)
            start_indices = d_i_left_indices
            end_indices = d_i_right_indices
        else:
            phi_2, _, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, _, _, _, _ = reco(np.array([x[d_i_right_indices], y[d_i_right_indices]]), charges[d_i_right_indices], phi_1)
            start_indices = d_i_right_indices
            end_indices = d_i_left_indices
        if phi_2 == 0 or phi_2 == np.pi or phi_2 == -np.pi or phi_2 == np.pi/2 or phi_2 == -np.pi/2: 
            return np.nan, np.nan, np.empty(0), np.empty(0)
        else:
            return phi_1, phi_2, start_indices, end_indices
    except:
        return np.nan, np.nan, np.empty(0), np.empty(0)

# Two step reconstruction for Timepix3 data
def reco_tpx3(x, y, toa, ftoa, charges, full3d):
    try:
        z = ((toa * 25 + 1) - ftoa * 1.5625)*8.3/55.
        phi_1, theta_1, d_i_left_indices, d_i_right_indices, m3, xc, yc, zc, m3_z, d_i_upper_indices, d_i_lower_indices = reco(np.array([x, y, z]), charges)
        # Decide based on the third moment which part to keep. m3 is the 3rd moment in the xy plane, m3_z for the full 3D information
        if m3 <= 0:
            if m3_z <= 0:
                if full3d:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_left_indices, d_i_upper_indices)], y[np.intersect1d(d_i_left_indices, d_i_upper_indices)], z[np.intersect1d(d_i_left_indices, d_i_upper_indices)]]), charges[np.intersect1d(d_i_left_indices, d_i_upper_indices)], phi_1)
                else:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_left_indices, d_i_upper_indices)], y[np.intersect1d(d_i_left_indices, d_i_upper_indices)]]), charges[np.intersect1d(d_i_left_indices, d_i_upper_indices)], phi_1)
                start_indices = np.intersect1d(d_i_left_indices, d_i_upper_indices)
                end_indices = np.setdiff1d(np.arange(len(x)), start_indices)
            else:
                if full3d:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_left_indices, d_i_lower_indices)], y[np.intersect1d(d_i_left_indices, d_i_lower_indices)], z[np.intersect1d(d_i_left_indices, d_i_lower_indices)]]), charges[np.intersect1d(d_i_left_indices, d_i_lower_indices)], phi_1)
                else:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_left_indices, d_i_lower_indices)], y[np.intersect1d(d_i_left_indices, d_i_lower_indices)]]), charges[np.intersect1d(d_i_left_indices, d_i_lower_indices)], phi_1)
                start_indices = np.intersect1d(d_i_left_indices, d_i_lower_indices)
                end_indices = np.setdiff1d(np.arange(len(x)), start_indices)
        else:
            if m3_z <= 0:
                if full3d:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_right_indices, d_i_upper_indices)], y[np.intersect1d(d_i_right_indices, d_i_upper_indices)], z[np.intersect1d(d_i_right_indices, d_i_upper_indices)]]), charges[np.intersect1d(d_i_right_indices, d_i_upper_indices)], phi_1)
                else:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_right_indices, d_i_upper_indices)], y[np.intersect1d(d_i_right_indices, d_i_upper_indices)]]), charges[np.intersect1d(d_i_right_indices, d_i_upper_indices)], phi_1)
                start_indices = np.intersect1d(d_i_right_indices, d_i_upper_indices)
                end_indices = np.setdiff1d(np.arange(len(x)), start_indices)
            else:
                if full3d:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_right_indices, d_i_lower_indices)], y[np.intersect1d(d_i_right_indices, d_i_lower_indices)], z[np.intersect1d(d_i_right_indices, d_i_lower_indices)]]), charges[np.intersect1d(d_i_right_indices, d_i_lower_indices)], phi_1)
                else:
                    phi_2, theta_2, d_i_left_indices2, d_i_right_indices2, m32, xc2, yc2, zc2, m32_z, d_i_upper_indices2, d_i_lower_indices2 = reco(np.array([x[np.intersect1d(d_i_right_indices, d_i_lower_indices)], y[np.intersect1d(d_i_right_indices, d_i_lower_indices)]]), charges[np.intersect1d(d_i_right_indices, d_i_lower_indices)], phi_1)
                start_indices = np.intersect1d(d_i_right_indices, d_i_lower_indices)
                end_indices = np.setdiff1d(np.arange(len(x)), start_indices)
        if phi_2 == 0 or phi_2 == np.pi or phi_2 == -np.pi or phi_2 == np.pi/2 or phi_2 == -np.pi/2: 
            return np.nan, np.nan, np.empty(0), np.empty(0)
        else:
            return phi_1, phi_2, start_indices, end_indices
    except:
        return np.nan, np.nan, np.empty(0), np.empty(0)

def tpx_wrapper(inputs):
    return reco_tpx(*inputs)

def tpx3_wrapper(inputs):
    return reco_tpx3(*inputs)

def main():
    # Get the arguments
    parser = argparse.ArgumentParser(description='Angular reconstruction of Timepix and Timepix3 polarimetry data')
    parser.add_argument('runpath', type=str, help='Path to the hdf5 file')
    parser.add_argument('--rotation', type=float, help='Rotation of the input coordinates in the xy plane with respect to the x axis. Given in degrees.', default=0)
    parser.add_argument('--full2d', action='store_true', help='Analyze Timepix3 data only in 2D')
    parser.add_argument('--full3d', action='store_true', help='Analyze Timepix3 data in the full 3D approach instead of 3D for the first step and 2D for the second step')
    args = parser.parse_args()

    # Process the arguments
    run = args.runpath
    if run.endswith('h5'):
        datafile = os.path.basename(run)
    else:
        print("Please choose a correct data file")
    angle_offset = args.rotation
    tpx3_2d = args.full2d
    tpx3_full3d = args.full3d

    # Open the corresponding datafile
    f = h5py.File(run, 'r+')
    filename = datafile.replace('.h5', '')

    timepix_version = f['reconstruction'].attrs['TimepixVersion'][0].decode('utf-8')
    reconstruction = f['reconstruction']
    for name in reconstruction:
        # Load the x and y coordinates of the pixel with the option to rotate
        posx_raw = f.get('reconstruction/' + name + '/chip_0/x')[:]
        posy_raw = f.get('reconstruction/' + name + '/chip_0/y')[:]
        rot = np.radians(angle_offset)
        posx = posx_raw * np.cos(rot) - posy_raw * np.sin(rot)
        posy = posx_raw * np.sin(rot) + posy_raw * np.cos(rot)

        # Fot Timepix3 use additionally the timestamp per pixel
        if timepix_version == 'Timepix3':
            toa = f.get('reconstruction/' + name + '/chip_0/ToA')[:]
            ftoa = f.get('reconstruction/' + name + '/chip_0/fToA')[:]

        # If there is calibrated charge data use it, otherwise use the ToT
        try:
            charge = f.get('reconstruction/' + name + '/chip_0/charge')[:]
        except:
            charge = f.get('reconstruction/' + name + '/chip_0/ToT')[:]

        print(timepix_version)
        if timepix_version == 'Timepix1' or tpx3_2d:
            inputs = list(zip(posx, posy, charge))
            # Perform the reconstruction per event in multithreading
            with multiprocessing.Pool(processes=8) as pool:
                results = list(tqdm(pool.imap(tpx_wrapper, inputs), total=len(inputs)))
            pool.close()
            pool.join()

            # Split up the results in separate arrays
            results = np.array(results, dtype=object)
            phi1 = np.array(results[:, 0], dtype=np.float64)
            phi2 = np.array(results[:, 1], dtype=np.float64)
            start = results[:, 2]
            end = results[:, 3]

        elif timepix_version == 'Timepix3':
            inputs = list(zip(posx, posy, toa, ftoa, charge, [tpx3_full3d]*len(posx)))
            # Perform the reconstruction per event in multithreading
            with multiprocessing.Pool(processes=8) as pool:
                results = list(tqdm(pool.imap(tpx3_wrapper, inputs), total=len(inputs)))
            pool.close()
            pool.join()

            # Split up the results in separate arrays
            results = np.array(results, dtype=object)
            phi1 = np.array(results[:, 0], dtype=np.float64)
            phi2 = np.array(results[:, 1], dtype=np.float64)
            start = results[:, 2]
            end = results[:, 3]

        # Save data to the hdf5 file
        f.create_dataset('reconstruction/' + name + '/chip_0/angle_fiststage', data=phi1)
        f.create_dataset('reconstruction/' + name + '/chip_0/angle_secondstage', data=phi2)
        dt = h5py.special_dtype(vlen=np.dtype('float64'))
        f.create_dataset('reconstruction/' + name + '/chip_0/start_indices', (len(start),), dtype=dt)
        f['reconstruction/' + name + '/chip_0/start_indices'][...] = start
        f.create_dataset('reconstruction/' + name + '/chip_0/end_indices', (len(end),), dtype=dt)
        f['reconstruction/' + name + '/chip_0/end_indices'][...] = end

if __name__ == "__main__":
    main()
