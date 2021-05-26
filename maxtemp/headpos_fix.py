#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 11:30:57 2021
@author: mikkel
"""
#%% PATHS
raw_path = '/archive/20061_tinnitus/MEG'
out_path = '/home/mikkel/tinmeg/MEG'               #### ! CHANGE THIS !

# ADD SUBJECTS/DATES HERE
subj_dates = [
    #     'NatMEG_0756/210224',
        
    # 'NatMEG_0538'
    # 'NatMEG_0539
    # 'NatMEG_0697
    # 'NatMEG_0750
    # 'NatMEG_0756
    # 'NatMEG_0832
    # 'NatMEG_0835
    # 'NatMEG_0836
    # 'NatMEG_0838
    # 'NatMEG_0839
    
    'NatMEG_0839/210511'
    # NatMEG_0840/210512
    # NatMEG_0841/210512
    # 'NatMEG_0842/210513'
    # 'NatMEG_0844/210517'
    # 'NatMEG_0845/210517'
    # 'NatMEG_0847/210519'
    ]

# ADD FILENAME
condition = 'tinmeg'

#%% IMPORT STUFF
import os.path as op
import numpy as np
import re
import warnings
from os import listdir, mkdir
from mne.transforms import rotation_angles, rotation3d, write_trans, quat_to_rot
from mne.io import read_raw_fif

#%% Find files
def find_condition_files(folder, string):
    allfiles = listdir(folder)
    strfiles = [f for f in allfiles if string in f and f.find('-') == -1 and not 'sss' in f and not '_avg' in f]
    strfiles.sort()
    return strfiles

#%% RUN
for ii, subj_date in enumerate(subj_dates):
    subj, date = subj_date.split('/')
    print('Processing subj: '+subj+' ('+str(ii+1)+' of '+str(len(subj_dates))+').')
    # Check that the method works 
    method = 'median'
    
    # Get and set folders
    rawdir =  op.join(raw_path, subj_date)                                  # [!] Match up with bash script !
    
    print(rawdir)
    quatdir = op.join(rawdir,'quat_files')
    
    mean_trans_folder = op.join(out_path, subj, 'trans_files')
    if not op.exists(op.join(out_path, subj)):
        mkdir(op.join(out_path, subj))
    if not op.exists(op.join(out_path, subj_date)):
        mkdir(op.join(out_path, subj_date))     
    if not op.exists(mean_trans_folder):                      # Make sure output folder exists
        mkdir(mean_trans_folder)
        
    mean_trans_file = op.join(mean_trans_folder, condition+'-trans.fif')
    if op.isfile(mean_trans_file):
        warnings.warn('N"%s\" already exists is %s. Delete if you want to rerun' % (mean_trans_file, mean_trans_folder), RuntimeWarning)
        continue
    
    # Change to subject dir     
    files2combine = find_condition_files(quatdir, condition)     
    files2combine.sort()
    
    if not files2combine:
        raise RuntimeError('No files called \"%s\" found in %s' % (condition, quatdir))
        
    allfiles = []
    for ff in files2combine:
        fl = ff.split('_')[0]
        tmplist = [f for f in listdir(quatdir) if fl in f and '_quat' in f]
        
        #Fix order
        if len(tmplist) > 1:
            tmplist.sort()
            if any("-" in f for f in tmplist):
                firstfile = tmplist[-1]  # The file without a number will always be last!  
                tmpfs = sorted(tmplist[:-1], key=lambda a: int(re.split('-|.fif', a)[-2]) )  # Assuming consistent naming!!!
                tmplist[0] = firstfile
                tmplist[1:] = tmpfs
                allfiles = allfiles + tmplist
        
    
    if len(allfiles) > 1:
        print('Files used for average head pos:')    
        for ib in range(len(allfiles)):
            print('{:d}: {:s}'.format(ib + 1, allfiles[ib]))
    else:
        print('Will find average head pos in %s' % files2combine)    
    
    # LOAD DATA
    # raw = read_raw_fif(op.join(quatdir,firstfile), preload=True, allow_maxshield=True, verbose=False).pick_types(meg=False, chpi=True)
    # Use files2combine instead of allfiles as MNE will find split files automatically.
    for idx, ffs in enumerate(files2combine):
        if idx == 0:
            raw = read_raw_fif(op.join(quatdir,ffs), preload=False, allow_maxshield=True).pick_types(meg=False, chpi=True)
        else:
            raw.append(read_raw_fif(op.join(quatdir,ffs), preload=False, allow_maxshield=True).pick_types(meg=False, chpi=True))
        
    quat, times = raw.get_data(return_times=True)
    gof = quat[6,]                                              # Godness of fit channel
    # fs = raw.info['sfreq']
    
    # In case "record raw" started before "cHPI"
    if np.any(gof < 0.98):
        begsam = np.argmax(gof>0.98)
        raw.crop(tmin=raw.times[begsam])
        quat = quat[:,begsam:].copy()
        times = times[begsam:].copy()
        
    # Get continous transformation    
    print('Reading transformation. This will take a while...')
    H = np.empty([4,4,len(times)])                              # Initiate transforms
    init_rot_angles = np.empty([len(times),3])
        
    for i,t in enumerate(times):
        Hi = np.eye(4,4)
        Hi[0:3,3] = quat[3:6,i].copy()
        Hi[:3,:3] = quat_to_rot(quat[0:3,i])
        init_rot_angles[i,:] = rotation_angles(Hi[:3,:3])
        assert(np.sum(Hi[-1]) == 1.0)                           # sanity check result
        H[:,:,i] = Hi.copy()
    

    H_mean = np.median(H, axis=2)                 # stack, then average over new dim
    mean_rot_xfm = rotation3d(*tuple(np.median(init_rot_angles, axis=0)))  # stack, then average, then make new xfm        
        
    H_mean[:3,:3] = mean_rot_xfm
    assert(np.sum(H_mean[-1]) == 1.0)  # sanity check result
    
    # Create the mean structure and save as .fif    
    mean_trans = raw.info['dev_head_t']  # use the last info as a template
    mean_trans['trans'] = H_mean.copy()
    
    # Write file
    write_trans(mean_trans_file, mean_trans)
    print("Wrote "+mean_trans_file)

#END