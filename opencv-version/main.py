import copy
import cv2
import h5py
import numpy as np
import os
from skimage import io
import scipy.io as sio
import matplotlib.pyplot as plt
import tifffile as tiff
from utils import *

PWS_Path = "Dataset1/PWSimages"
DV_Path = "Dataset1/DVimages"
Result_Path = "Dataset1/Results"


### TO DO:
# 0. Check if any files is missing
# 1. Read PWS and Confocal Image and Mask
# 2. Use the PWS and Confocal mask to do the image matching
# 3. Image * Mask -> Cell Image
# 4. Use the bounding box of cell image to crop
# 5. Check the angle of cell in pws and confocal image, rotate them to the same angle
# 6. Generate same size cell image



for folder in os.listdir(PWS_Path):
    print("Current folder:", folder)

    ############################################### 0. CHECK ###############################################
    pws_filepath = os.path.join(PWS_Path, folder, "PWS/analyses/analysisResults_p0.h5")
    roi_filepath = os.path.join(PWS_Path, folder, "ROI_nuc.h5")

    dv_filepath = os.path.join(DV_Path, folder, "R3D_D3D-processed.tif")
    label_filepath = os.path.join(DV_Path, folder, "labels.mat")

    if not(os.path.exists(pws_filepath)) or not(os.path.exists(roi_filepath)) or not(os.path.exists(dv_filepath)) or not(os.path.exists(label_filepath)):
        print("Check if the files exist!")
        continue



    ############################################### 1. READ ###############################################
    ### PWS Image
    ## Raw PWS Image
    h5file = h5py.File(pws_filepath, 'r')
    pws = h5file['rms']
    pws = np.array(pws)

    ## ROI Mask
    h5file = h5py.File(roi_filepath, 'r')
    keys = [key for key in h5file.keys()]
    roi = h5file[keys[0]]['mask']
    roi = np.array(roi)

    ### DV Image
    dv = io.imread(dv_filepath)

    ## Mask
    label = sio.loadmat(label_filepath)
    label = label['labels']
    label = np.transpose(label, (2, 0, 1))

    label_index = []
    for i in range((label.shape)[0]):
        if label[i].max() != 0:
            label_index.append(i)



    ############################################### 2. Matching ###############################################
    pws_mask = copy.deepcopy(roi)
    dv_mask = copy.deepcopy(label)

    hu_distance = []
    for idx in label_index:
        d2 = cv2.matchShapes(pws_mask, dv_mask[idx], cv2.CONTOURS_MATCH_I2,0)
        hu_distance.append((idx, d2))
    hu_distance = np.asarray(hu_distance)



    ############################################### 3. Cell Image ###############################################
    ### PWS Cell Image
    pws_cell = np.multiply(pws, pws_mask.astype('float32'))

    ### DV Cell Image
    hu_distance = hu_distance[hu_distance[:, 1].argsort()]

    if len(hu_distance) == 1:
        idx_1 = int(hu_distance[0][0])
        dv_cell_1 = np.multiply(dv[idx_1], dv_mask[idx_1].astype('float32'))
        ave_cell = dv_cell_1
        ave_mask = dv_mask[idx_1]
    elif len(hu_distance) == 2:
        idx_1 = int(hu_distance[0][0])
        idx_2 = int(hu_distance[1][0])
        dv_cell_1 = np.multiply(dv[idx_1], dv_mask[idx_1].astype('float32'))
        dv_cell_2 = np.multiply(dv[idx_2], dv_mask[idx_2].astype('float32'))
        ave_cell = np.add(dv_cell_1, dv_cell_2) / 2.0
        ave_mask = np.add(dv_mask[idx_1], dv_mask[idx_2]) / 2.0
    else:
        idx_1 = int(hu_distance[0][0])
        idx_2 = int(hu_distance[1][0])
        idx_3 = int(hu_distance[2][0])
        dv_cell_1 = np.multiply(dv[idx_1], dv_mask[idx_1].astype('float32'))
        dv_cell_2 = np.multiply(dv[idx_2], dv_mask[idx_2].astype('float32'))
        dv_cell_3 = np.multiply(dv[idx_3], dv_mask[idx_3].astype('float32'))
        ave_cell = np.add(np.add(dv_cell_1, dv_cell_2), dv_cell_3) / 3.0
        ave_mask = np.add(np.add(dv_mask[idx_1], dv_mask[idx_2]), dv_mask[idx_3]) / 3.0

    ## make the mask to cover as many cell elements as possible
    ave_mask[ave_mask > 0] = 1



    ############################################### 4. Crop Image ###############################################
    ### Use the mask to find contours.
    ## PWS Mask
    pws_mask_contours, _ = cv2.findContours(pws_mask, cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
    pws_x, pws_y, pws_w, pws_h = cv2.boundingRect(pws_mask_contours[0])

    # Crop
    cropped_pws = pws_cell[pws_y:pws_y+pws_h, pws_x:pws_x+pws_w]


    ## Confocal Mask
    ave_mask = ave_mask.astype(np.uint8).copy()
    kernel = np.ones((5, 5), np.uint8)
    ave_mask = cv2.dilate(ave_mask, kernel, iterations=3)

    confocal_contours, _ = cv2.findContours(ave_mask, cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    max_contour = max(confocal_contours, key=cv2.contourArea)
    confocal_x, confocal_y, confocal_w, confocal_h = cv2.boundingRect(max_contour)

    # Crop
    cropped_confocal = ave_cell[confocal_y:confocal_y+confocal_h, confocal_x:confocal_x+confocal_w]

    ## Make sure pws and confocal in the same size
    resized_cropped_pws = cv2.resize(cropped_pws, [confocal_w, confocal_h])



    ###############################################5. Angle (BUG)#####################################################
    ## PWS Image
    ret, bw_pws = cv2.threshold(resized_cropped_pws, 0.00000001, 1, cv2.THRESH_BINARY)
    bw_pws = bw_pws.astype(np.uint8)
    pws_contours, _ = cv2.findContours(bw_pws, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    pws_max_contour = max(pws_contours, key=cv2.contourArea)
    pws_angle = getAngle(pws_max_contour, bw_pws)

    ## Confocal Image
    ret, bw_confocal = cv2.threshold(cropped_confocal, 0.00000001, 1, cv2.THRESH_BINARY)
    bw_confocal = bw_confocal.astype(np.uint8)
    confocal_contours, _ = cv2.findContours(bw_confocal, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    confocal_max_contour = max(confocal_contours, key=cv2.contourArea)
    confocal_angle = getAngle(confocal_max_contour, bw_confocal)

    ## Rotate PWS Image -> size will change
    # clockwise is negative, otherwise positive.
    angle_dif = confocal_angle - pws_angle
    angle_dif = angle_dif if angle_dif > 0 else -angle_dif
    rotated_pws = rotate_bound(resized_cropped_pws, angle_dif)

    ## Check PWS angle
    ret, bw_pws_check = cv2.threshold(rotated_pws, 0.00000001, 1, cv2.THRESH_BINARY)
    bw_pws_check = bw_pws_check.astype(np.uint8)
    pws_contours_check, _ = cv2.findContours(bw_pws_check, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    pws_max_contour_check = max(pws_contours_check, key=cv2.contourArea)
    pws_angle_check = getAngle(pws_max_contour_check, bw_pws_check)

    if abs(pws_angle_check - confocal_angle) >= 5:
        continue



    ###############################################6. Same Size ###############################################
    # ## Padding to (480,480)
    desired_size = [512, 512]

    resized_cropped_pws_h, resized_cropped_pws_w = rotated_pws.shape[0], rotated_pws.shape[1]
    cropped_confocal_h, cropped_confocal_w = cropped_confocal.shape[0], cropped_confocal.shape[1]

    pws_delta_h = desired_size[0] - resized_cropped_pws_h
    pws_delta_w = desired_size[1] - resized_cropped_pws_w
    pws_top, pws_bottom = pws_delta_h//2, pws_delta_h - (pws_delta_h//2)
    pws_left, pws_right = pws_delta_w//2, pws_delta_w - (pws_delta_w//2)

    confocal_delta_h = desired_size[0] - cropped_confocal_h
    confocal_delta_w = desired_size[1] - resized_cropped_pws_w
    confocal_top, confocal_bottom = confocal_delta_h//2, confocal_delta_h//2
    confocal_left, confocal_right = confocal_delta_w//2, confocal_delta_w//2

    padding_pws = cv2.copyMakeBorder(resized_cropped_pws, pws_top, pws_bottom, pws_left, pws_right, cv2.BORDER_CONSTANT, None, value = 0)
    padding_confocal = cv2.copyMakeBorder(cropped_confocal, confocal_top, confocal_bottom, confocal_left, confocal_right, cv2.BORDER_CONSTANT, None, value = 0)

    os.mkdir(os.path.join(Result_Path, folder))
    pwssavepath = os.path.join(Result_Path, folder, "pws.tif")
    cv2.imwrite(pwssavepath, padding_pws)

    confocalsavepath = os.path.join(Result_Path, folder, "confocal.tif")
    cv2.imwrite(confocalsavepath, padding_confocal)